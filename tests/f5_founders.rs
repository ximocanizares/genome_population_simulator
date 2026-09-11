//! F5.3.1. Modo por bloque×subpoblación, cuotas controladas y muestreo
//! probabilístico.
//!
//! Las pruebas recorren la cadena productiva completa —configuración,
//! admisión sellada, layout, marginales, búsqueda y fundadores— y contrastan
//! el reparto controlado con un oráculo racional escrito en enteros de 128
//! bits que no llama a ninguna función del productor. Los perfiles usan pesos
//! **diádicos** y `C` potencia de dos, de modo que `C*q_h`, sus pisos y sus
//! restos son exactos en binary64 y la comparación con el oráculo es una
//! igualdad, no una tolerancia.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_founders::{
    F5FounderAssignment, F5SamplingMode, F5SamplingSource, assign_f5_founders,
};
use genome_population_simulator::f5_layout::{F5LayoutRequest, build_f5_physical_layout};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_search::{F5HaplotypeSearch, search_f5_haplotypes};
use genome_population_simulator::identity::{ExecutionId, ScenarioId};
use genome_population_simulator::randomness::{
    F5RandomDecision, RootSeed, derive_f5_subseed, f5_uniform_index,
};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";
const CONTINUOUS_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n";

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-founders-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Subpoblación declarada del escenario.
#[derive(Clone)]
struct Subpopulation<'a> {
    id: &'a str,
    individual_count: u64,
    /// `None` deja el predeterminado del contrato, que es `probabilistic`.
    sampling: Option<&'a str>,
    inbreeding_f: Option<&'a str>,
    /// Declara `founder_keys` en lugar de `founder_ids`.
    keys: bool,
}

impl<'a> Subpopulation<'a> {
    const fn new(id: &'a str, individual_count: u64) -> Self {
        Self {
            id,
            individual_count,
            sampling: None,
            inbreeding_f: None,
            keys: false,
        }
    }

    const fn sampling(mut self, mode: &'a str) -> Self {
        self.sampling = Some(mode);
        self
    }

    const fn inbreeding(mut self, value: &'a str) -> Self {
        self.inbreeding_f = Some(value);
        self
    }

    const fn keys(mut self) -> Self {
        self.keys = true;
        self
    }
}

/// Fila `[[sampling_overrides]]` declarada.
#[derive(Clone)]
struct Substitution<'a> {
    population: &'a str,
    block_index: u64,
    mode: &'a str,
    /// Literal decimal exacto; se escribe tal cual en el TOML.
    tolerance: Option<&'a str>,
}

struct Fixture<'a> {
    subpopulations: Vec<Subpopulation<'a>>,
    overrides: Vec<Substitution<'a>>,
    block_count: u64,
    block_length: u64,
    weights: Vec<&'a str>,
    gap: u64,
    spectrum: &'a str,
    bins: Vec<(u64, u64, f64, f64)>,
    max_haplotype_attempts: u64,
    max_swap_proposals: u64,
    ploidy: u64,
    seed: &'a str,
    slack: u64,
}

impl Default for Fixture<'_> {
    fn default() -> Self {
        Self {
            subpopulations: vec![Subpopulation::new("P1", 1)],
            overrides: Vec::new(),
            block_count: 2,
            block_length: 3,
            // `H=2` mantiene `r²` en 0 o 1, de modo que el objetivo 1,0 es
            // alcanzable y la búsqueda termina en pocos intentos.
            weights: vec!["0.5", "0.5"],
            gap: 1,
            spectrum: CONTINUOUS_SPECTRUM,
            bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05)],
            max_haplotype_attempts: 32,
            max_swap_proposals: 256,
            ploidy: 2,
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
        let mut sorted = self
            .subpopulations
            .iter()
            .map(|subpopulation| subpopulation.id)
            .collect::<Vec<_>>();
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
                    self.gap, self.gap
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
                "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{}\"\nploidy = {}\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
                self.seed,
                self.ploidy,
                sha256(fasta.as_bytes())
            ),
        )
        .expect("general");

        let mut subpopulations = String::new();
        let mut founder_count = 0_u64;
        for subpopulation in &self.subpopulations {
            founder_count += subpopulation.individual_count;
            let founders = (0..subpopulation.individual_count)
                .map(|index| format!("\"{}-I{index}\"", subpopulation.id))
                .collect::<Vec<_>>()
                .join(", ");
            let field = if subpopulation.keys {
                "founder_keys"
            } else {
                "founder_ids"
            };
            subpopulations.push_str(&format!(
                "\n[[subpopulations]]\npopulation_id = \"{}\"\nindividual_count = {}\n{field} = [{founders}]\n",
                subpopulation.id, subpopulation.individual_count
            ));
            if let Some(mode) = subpopulation.sampling {
                subpopulations.push_str(&format!("sampling = \"{mode}\"\n"));
            }
            if let Some(value) = subpopulation.inbreeding_f {
                subpopulations.push_str(&format!("inbreeding_f = {value}\n"));
            }
        }
        let mut overrides = String::new();
        for substitution in &self.overrides {
            overrides.push_str(&format!(
                "\n[[sampling_overrides]]\npopulation_id = \"{}\"\nblock_index = {}\nmode = \"{}\"\n",
                substitution.population, substitution.block_index, substitution.mode
            ));
            if let Some(tolerance) = substitution.tolerance {
                overrides.push_str(&format!("controlled_max_abs_error = {tolerance}\n"));
            }
        }
        fs::write(
            scenario_directory.join("scenario.toml"),
            format!(
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {founder_count}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = {}\nmax_swap_proposals = {}\nmax_haplotype_cells = 65536\nblock_count = {}\n{subpopulations}{overrides}",
                sha256(metadata.as_bytes()),
                self.max_haplotype_attempts,
                self.max_swap_proposals,
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
    search: F5HaplotypeSearch,
    root: PathBuf,
}

impl Drop for Chain {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

/// Cadena completa hasta la búsqueda; exige que la búsqueda haya terminado,
/// porque sin matriz publicada no hay fundadores que asignar.
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
    let search = search_f5_haplotypes(&config, &sealed, &execution, &layout)
        .expect("búsqueda haplotípica terminada");
    Chain {
        config,
        sealed,
        search,
        root: root.to_path_buf(),
    }
}

fn assign(chain: &Chain) -> F5FounderAssignment {
    assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
        .expect("asignación fundadora completa")
}

/// Perfil **diádico** que termina: `H=4` con pesos `1/2, 1/4, 1/8, 1/8`.
///
/// Es el fixture de las cuotas: `C*q_h` y sus restos son exactos en binary64 y
/// el oráculo entero puede comprobarlos sin tolerancia.
fn dyadic_fixture<'a>() -> Fixture<'a> {
    Fixture {
        weights: vec!["0.5", "0.25", "0.125", "0.125"],
        block_length: 3,
        bins: vec![(1, 2, 0.25, 0.05)],
        max_haplotype_attempts: 512,
        max_swap_proposals: 2048,
        ..Fixture::default()
    }
}

/// Perfil de cinco haplotipos equiprobables, que también termina.
///
/// Con `C` pequeño ninguna cuota entera cabe dentro de `0,5/C`: es el fixture
/// de la cuota imposible.
fn five_haplotype_fixture<'a>() -> Fixture<'a> {
    Fixture {
        weights: vec!["0.2", "0.2", "0.2", "0.2", "0.2"],
        block_length: 4,
        bins: vec![(1, 3, 0.2, 0.05)],
        max_haplotype_attempts: 512,
        max_swap_proposals: 2048,
        ..Fixture::default()
    }
}

// ---------------------------------------------------------------------------
// Oráculo racional exacto del reparto por mayores restos
// ---------------------------------------------------------------------------

/// Pesos declarados como numeradores enteros sobre un denominador diádico.
///
/// Los perfiles de esta batería declaran pesos que ya suman exactamente 1 y
/// son múltiplos de una potencia de dos, de modo que la normalización DTC-049
/// los divide por `1.0` y los deja intactos: el peso efectivo es exactamente
/// `numerador/denominador` y el oráculo trabaja sin perder un bit. La
/// aserción de la suma lo comprueba en cada llamada.
fn integer_weights(declared: &[&str]) -> (Vec<u128>, u128) {
    let values: Vec<f64> = declared
        .iter()
        .map(|weight| weight.parse::<f64>().expect("peso decimal"))
        .collect();
    let mut denominator: u128 = 1;
    while values
        .iter()
        .any(|value| (value * denominator as f64).fract() != 0.0)
    {
        denominator *= 2;
        assert!(denominator <= 1 << 40, "el oráculo exige pesos diádicos");
    }
    let numerators: Vec<u128> = values
        .iter()
        .map(|value| (value * denominator as f64) as u128)
        .collect();
    assert_eq!(
        numerators.iter().sum::<u128>(),
        denominator,
        "los pesos declarados deben sumar exactamente 1"
    );
    (numerators, denominator)
}

/// Suma de Neumaier escrita desde su definición, para el oráculo del paseo de
/// prefijos.
#[derive(Default)]
struct NeumaierOracle {
    sum: f64,
    compensation: f64,
}

impl NeumaierOracle {
    fn add(&mut self, value: f64) -> f64 {
        let next = self.sum + value;
        self.compensation += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
        self.sum + self.compensation
    }
}

/// Comprueba en enteros exactos todo lo que el reparto por mayores restos
/// debe cumplir, sin reproducir ningún sorteo:
///
/// 1. cada `n_h` es el piso o el piso más uno de `C*q_h`;
/// 2. la suma es exactamente `C`;
/// 3. el número de incrementos es el residual `C - Σ pisos`;
/// 4. ningún rango incrementado tiene un resto menor que uno no incrementado;
/// 5. cada `n_h` respeta la cota `|n_h/C - q_h| <= t` como racional exacto.
fn assert_largest_remainder(
    counts: &[u64],
    declared: &[&str],
    copies: u64,
    tolerance: (u128, u128),
) {
    let (numerators, denominator) = integer_weights(declared);
    let (tolerance_numerator, tolerance_denominator) = tolerance;
    assert_eq!(counts.len(), numerators.len(), "un recuento por rango");
    let copies = u128::from(copies);
    let mut residual = 0_u128;
    let mut increments = 0_u128;
    let mut incremented_remainders: Vec<u128> = Vec::new();
    let mut untouched_remainders: Vec<u128> = Vec::new();
    let mut total = 0_u128;
    for (rank, numerator) in numerators.iter().enumerate() {
        let ideal = copies * numerator;
        let floor = ideal / denominator;
        let remainder = ideal % denominator;
        residual += remainder;
        let count = u128::from(counts[rank]);
        assert!(
            count == floor || count == floor + 1,
            "rango {rank}: {count} no es el piso {floor} ni su sucesor"
        );
        if count == floor + 1 {
            increments += 1;
            incremented_remainders.push(remainder);
        } else {
            untouched_remainders.push(remainder);
        }
        total += count;

        // `|n/C - w/D| <= tn/td`  ⟺  `|n*D - C*w| * td <= tn * C * D`.
        let left = (count * denominator) as i128 - (copies * numerator) as i128;
        assert!(
            left.unsigned_abs() * tolerance_denominator
                <= tolerance_numerator * copies * denominator,
            "rango {rank}: la cuota {count} incumple la cota exacta"
        );
    }
    assert_eq!(total, copies, "la cuota no suma C");
    assert_eq!(
        residual % denominator,
        0,
        "la suma de restos debe ser un entero de copias"
    );
    assert_eq!(
        increments,
        residual / denominator,
        "el número de incrementos debe ser el residual"
    );
    let minimum_incremented = incremented_remainders.iter().min().copied().unwrap_or(0);
    let maximum_untouched = untouched_remainders.iter().max().copied().unwrap_or(0);
    if !incremented_remainders.is_empty() && !untouched_remainders.is_empty() {
        assert!(
            minimum_incremented >= maximum_untouched,
            "un rango incrementado tiene menos resto que otro sin incrementar"
        );
    }
}

// ---------------------------------------------------------------------------
// Resolución de modo y tolerancia
// ---------------------------------------------------------------------------

/// Fila resuelta tal como la comparan las pruebas de modo: modo solicitado,
/// efectivo, origen y las dos tolerancias.
type ResolvedRow<'a> = (&'a str, &'a str, &'a str, Option<f64>, Option<f64>);

/// Dos bloques × dos subpoblaciones con modos distintos, con y sin
/// sustitución.
#[test]
fn two_blocks_and_two_subpopulations_resolve_their_own_mode() {
    let root = temporary("modes");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 1).sampling("controlled"),
            Subpopulation::new("P2", 1),
        ],
        overrides: vec![
            Substitution {
                population: "P1",
                block_index: 1,
                mode: "probabilistic",
                tolerance: None,
            },
            Substitution {
                population: "P2",
                block_index: 0,
                mode: "controlled",
                tolerance: Some("0.25"),
            },
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    assert_eq!(assignment.sampling().len(), 4, "una fila por unidad");

    let observed: BTreeMap<(String, u64), ResolvedRow<'_>> = assignment
        .sampling()
        .iter()
        .map(|row| {
            (
                (row.population_id().as_str().to_owned(), row.block_index()),
                (
                    row.requested_mode().as_str(),
                    row.effective_mode().as_str(),
                    row.source().as_str(),
                    row.controlled_max_abs_error_declared(),
                    row.controlled_max_abs_error_effective(),
                ),
            )
        })
        .collect();

    assert_eq!(
        observed[&("P1".to_owned(), 0)],
        (
            "controlled",
            "controlled",
            "subpopulation_default",
            None,
            Some(0.25)
        ),
        "el bloque sin sustitución hereda el modo de su subpoblación y la tolerancia 0,5/C"
    );
    assert_eq!(
        observed[&("P1".to_owned(), 1)],
        (
            "probabilistic",
            "probabilistic",
            "block_override",
            None,
            None
        ),
        "la sustitución cambia el modo del bloque y anula ambas tolerancias"
    );
    assert_eq!(
        observed[&("P2".to_owned(), 0)],
        (
            "controlled",
            "controlled",
            "block_override",
            Some(0.25),
            Some(0.25)
        ),
        "la sustitución declara su tolerancia y la efectiva la respeta"
    );
    assert_eq!(
        observed[&("P2".to_owned(), 1)],
        (
            "probabilistic",
            "probabilistic",
            "subpopulation_default",
            None,
            None
        ),
        "sin sustitución ni `sampling` el predeterminado es probabilístico"
    );

    for unit in assignment.units() {
        let expected = observed[&(unit.population_id().as_str().to_owned(), unit.block_index())];
        assert_eq!(unit.mode().as_str(), expected.1, "modo de la unidad");
        assert_eq!(
            unit.quota().is_some(),
            unit.mode() == F5SamplingMode::Controlled,
            "solo el controlado publica cuota"
        );
    }
}

/// Sin ninguna fila `[[sampling_overrides]]`, cada unidad conserva el
/// predeterminado de su subpoblación y lo declara así.
#[test]
fn an_absent_substitution_keeps_the_subpopulation_default() {
    let root = temporary("absent-substitution");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 1).sampling("controlled"),
            Subpopulation::new("P2", 1).sampling("probabilistic"),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for row in assignment.sampling() {
        assert_eq!(
            row.source(),
            F5SamplingSource::SubpopulationDefault,
            "ninguna fila puede venir de una sustitución inexistente"
        );
        let expected = if row.population_id().as_str() == "P1" {
            F5SamplingMode::Controlled
        } else {
            F5SamplingMode::Probabilistic
        };
        assert_eq!(row.requested_mode(), expected);
        assert_eq!(row.effective_mode(), expected);
    }
}

/// Una sustitución fuera de `0..B-1` se rechaza con `B` ya resuelto, antes de
/// que exista ningún layout y mucho antes de fundadores.
///
/// La comprobación vive en el preflight, que es quien resuelve `B`; la
/// resolución de modos vuelve a exigir que cada sustitución corresponda a una
/// unidad real y no puede consumir una que nadie reclame.
#[test]
fn a_substitution_out_of_range_is_rejected_before_any_layout() {
    let root = temporary("out-of-range");
    let fixture = Fixture {
        overrides: vec![Substitution {
            population: "P1",
            block_index: 2,
            mode: "controlled",
            tolerance: None,
        }],
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let failure = load_and_validate_f5_config(&general).expect_err("bloque inexistente");
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("sampling_overrides.block_index")
    );
    assert_eq!(
        failure.context().get("detail").map(String::as_str),
        Some("debe pertenecer a 0..1"),
        "el diagnóstico nombra el rango legal"
    );
    let _ = fs::remove_dir_all(&root);
}

/// El último índice legal sí se acepta: la frontera es `B-1`, no `B-2`.
#[test]
fn the_last_legal_block_index_is_accepted() {
    let root = temporary("last-index");
    let fixture = Fixture {
        overrides: vec![Substitution {
            population: "P1",
            block_index: 1,
            mode: "controlled",
            tolerance: None,
        }],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let row = assignment
        .sampling()
        .iter()
        .find(|row| row.block_index() == 1)
        .expect("bloque 1");
    assert_eq!(row.source(), F5SamplingSource::BlockOverride);
    assert_eq!(row.effective_mode(), F5SamplingMode::Controlled);
}

/// Dos sustituciones para la misma pareja ya las rechaza la admisión, antes de
/// que exista ninguna unidad que resolver.
#[test]
fn a_duplicated_substitution_is_rejected_by_admission() {
    let root = temporary("duplicated-substitution");
    let fixture = Fixture {
        overrides: vec![
            Substitution {
                population: "P1",
                block_index: 0,
                mode: "controlled",
                tolerance: None,
            },
            Substitution {
                population: "P1",
                block_index: 0,
                mode: "probabilistic",
                tolerance: None,
            },
        ],
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let failure = load_and_validate_f5_config(&general).expect_err("sustitución duplicada");
    assert!(
        failure
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("duplicada")),
        "la admisión nombra la duplicación: {:?}",
        failure.context()
    );
    let _ = fs::remove_dir_all(&root);
}

/// Una sustitución cuya población no existe la rechaza la admisión.
#[test]
fn a_substitution_for_an_unknown_population_is_rejected_by_admission() {
    let root = temporary("unknown-population");
    let fixture = Fixture {
        overrides: vec![Substitution {
            population: "P9",
            block_index: 0,
            mode: "controlled",
            tolerance: None,
        }],
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let failure = load_and_validate_f5_config(&general).expect_err("población desconocida");
    assert!(
        failure
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("desconocida")),
        "la admisión nombra la población ausente: {:?}",
        failure.context()
    );
    let _ = fs::remove_dir_all(&root);
}

/// En probabilístico las dos tolerancias son `null` y la unidad no publica
/// cuota.
#[test]
fn probabilistic_units_publish_null_tolerances_and_no_quota() {
    let root = temporary("null-tolerances");
    let chain = build_chain(&root, &Fixture::default());
    let assignment = assign(&chain);
    for row in assignment.sampling() {
        assert_eq!(row.effective_mode(), F5SamplingMode::Probabilistic);
        assert_eq!(
            row.controlled_max_abs_error_declared(),
            None,
            "el probabilístico prohíbe declarar tolerancia"
        );
        assert_eq!(
            row.controlled_max_abs_error_effective(),
            None,
            "el probabilístico no tiene tolerancia efectiva"
        );
    }
    for unit in assignment.units() {
        assert!(unit.quota().is_none(), "el probabilístico no reparte cuota");
    }
}

/// Las filas se publican ordenadas por `population_id` y `block_id`, que es el
/// orden literal de la lista `sampling` del efectivo.
///
/// No es el orden de `block_index`: `block_id` es un identificador derivado con
/// BLAKE3, de modo que ambos órdenes difieren, y este fixture lo demuestra.
#[test]
fn rows_are_ordered_by_population_and_block_id() {
    let root = temporary("row-order");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 2).sampling("controlled"),
            Subpopulation::new("P2", 5).sampling("controlled"),
        ],
        block_count: 3,
        slack: 2,
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let keys: Vec<(String, String)> = assignment
        .sampling()
        .iter()
        .map(|row| {
            (
                row.population_id().as_str().to_owned(),
                row.block_id().as_str().to_owned(),
            )
        })
        .collect();
    let mut sorted = keys.clone();
    sorted.sort();
    assert_eq!(keys, sorted, "las filas no siguen población+block_id");
    let indices: Vec<u64> = assignment
        .sampling()
        .iter()
        .take(3)
        .map(|row| row.block_index())
        .collect();
    assert_ne!(
        indices,
        vec![0, 1, 2],
        "este fixture debe separar el orden de block_id del de block_index"
    );
    let units: Vec<(String, u64)> = assignment
        .units()
        .iter()
        .map(|unit| (unit.population_id().as_str().to_owned(), unit.block_index()))
        .collect();
    let rows: Vec<(String, u64)> = assignment
        .sampling()
        .iter()
        .map(|row| (row.population_id().as_str().to_owned(), row.block_index()))
        .collect();
    assert_eq!(units, rows, "unidades y filas comparten orden");
}

// ---------------------------------------------------------------------------
// Modo probabilístico
// ---------------------------------------------------------------------------

/// Cada copia recibe exactamente una fila y su sorteo se reproduce desde el
/// ámbito publicado `founder_sample`, con un U53 por copia.
#[test]
fn every_probabilistic_copy_reproduces_its_own_scope() {
    // Con pesos diádicos y con pesos que no lo son: el segundo perfil es el
    // único en el que la suma compensada y una ingenua pueden divergir.
    for (label, profile) in [
        ("dyadic", dyadic_fixture()),
        ("five", five_haplotype_fixture()),
    ] {
        let root = temporary(&format!("probabilistic-scope-{label}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", 3)],
            ..profile
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
        let root_seed = RootSeed::parse("seed", SEED).expect("seed");
        let weights = chain.config.profile().haplotype_weights().to_vec();

        let mut checked = 0_u64;
        let mut observed = BTreeSet::new();
        for unit in assignment.units() {
            assert_eq!(unit.assignments().len(), 6, "tres individuos diploides");
            for copy in unit.assignments() {
                let mut rng = derive_f5_subseed(
                    root_seed,
                    F5RandomDecision::founder_sample(
                        &scenario_id,
                        unit.population_id(),
                        unit.block_id(),
                        copy.individual_id(),
                        copy.copy_index(),
                    ),
                )
                .expect("subsemilla")
                .local_rng();
                // U53 independiente: `bits53 / 2^53`.
                let threshold = (rng.next_u64() >> 11) as f64 / (1_u64 << 53) as f64;
                // Recorrido secuencial con la **suma compensada** que exige
                // DTC-089 y con el último rango positivo como residuo: una
                // suma ingenua o el último rango a secas no comprobarían ni la
                // compensación ni esa divergencia deliberada.
                let mut prefix = NeumaierOracle::default();
                let mut expected = weights
                    .iter()
                    .rposition(|weight| *weight > 0.0)
                    .expect("algún peso positivo");
                for (rank, weight) in weights.iter().enumerate() {
                    if threshold < prefix.add(*weight) {
                        expected = rank;
                        break;
                    }
                }
                assert_eq!(
                    copy.haplotype_rank(),
                    expected as u64,
                    "{label}: la copia {}#{} no reproduce su ámbito",
                    copy.individual_id(),
                    copy.copy_index()
                );
                observed.insert(copy.haplotype_rank());
                checked += 1;
            }
        }
        assert_eq!(checked, 12, "{label}: dos bloques × seis copias");
        assert!(
            observed.len() > 1,
            "{label}: el muestreo debe repartir copias entre varios haplotipos"
        );
    }
}

/// Un peso nulo no llega a fundadores: el perfil lo rechaza al admitirlo.
///
/// La robustez interna frente a un `q_h = 0` —prefijo que no crece y cuota
/// ideal cero— se prueba aparte, en las unitarias del módulo.
#[test]
fn the_profile_rejects_a_null_weight_before_founders() {
    let root = temporary("null-weight");
    let fixture = Fixture {
        weights: vec!["0.5", "0.5", "0"],
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let failure = load_and_validate_f5_config(&general).expect_err("peso nulo");
    assert_eq!(failure.code(), "GPS-FMT-101");
    assert_eq!(
        failure.context().get("detail").map(String::as_str),
        Some("peso inválido")
    );
    let _ = fs::remove_dir_all(&root);
}

// ---------------------------------------------------------------------------
// Modo controlado
// ---------------------------------------------------------------------------

/// El reparto suma `C`, lo realizan las copias y coincide con el oráculo
/// entero en cada unidad, con `C` pequeño y con `H` mayor que `C`.
#[test]
fn controlled_quotas_sum_c_and_match_the_integer_oracle() {
    for (label, individuals, ploidy) in [
        ("c1", 1_u64, 1_u64),
        ("c2", 1, 2),
        ("c4", 2, 2),
        ("c10", 5, 2),
    ] {
        let root = temporary(&format!("oracle-{label}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", individuals).sampling("controlled")],
            ploidy,
            ..dyadic_fixture()
        };
        let weights = fixture.weights.clone();
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let copies = individuals * ploidy;
        for unit in assignment.units() {
            let quota = unit.quota().expect("cuota controlada");
            assert_eq!(quota.copies(), copies);
            assert_eq!(
                quota.counts().iter().sum::<u64>(),
                copies,
                "{label}: la cuota debe sumar C"
            );
            assert_eq!(
                unit.assignments().len() as u64,
                copies,
                "{label}: una asignación por copia"
            );
            let mut realized = vec![0_u64; quota.counts().len()];
            for copy in unit.assignments() {
                realized[copy.haplotype_rank() as usize] += 1;
            }
            assert_eq!(
                realized,
                quota.counts(),
                "{label}: las copias no realizan la cuota"
            );
            // Cota efectiva predeterminada `0,5/C`, como racional `1/(2C)`.
            assert_largest_remainder(
                quota.counts(),
                &weights,
                copies,
                (1, u128::from(2 * copies)),
            );
        }
    }
}

/// La cota predeterminada `0,5/C` es **inclusiva**: una cuota cuyo error vale
/// exactamente `0,5/C` se acepta, y basta bajar la tolerancia declarada para
/// que la misma cuota se rechace.
#[test]
fn the_default_bound_is_inclusive_and_a_stricter_one_rejects_the_same_quota() {
    let inclusive = temporary("inclusive-bound");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2).sampling("controlled")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&inclusive, &fixture);
    let assignment = assign(&chain);
    let mut borderline = 0_u64;
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        assert_eq!(quota.max_abs_error_effective(), 0.125, "0,5/C con C=4");
        // Los dos rangos de peso 1/8 valen 0,5 copias: cualquiera de sus dos
        // estados queda exactamente sobre la cota.
        for rank in [2_usize, 3] {
            let error = (quota.counts()[rank] as f64 / 4.0 - 0.125).abs();
            assert_eq!(
                error,
                quota.max_abs_error_effective(),
                "la cuota admitida está justo en el borde de la cota"
            );
            borderline += 1;
        }
    }
    assert_eq!(borderline, 4, "dos rangos fronterizos en cada bloque");

    let stricter = temporary("stricter-bound");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2).sampling("controlled")],
        overrides: vec![Substitution {
            population: "P1",
            block_index: 0,
            mode: "controlled",
            tolerance: Some("0.0625"),
        }],
        ..dyadic_fixture()
    };
    let chain = build_chain(&stricter, &fixture);
    let failure = assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
        .expect_err("la tolerancia estricta no admite ninguna cuota entera");
    assert_eq!(failure.code(), "GPS-DOMAIN-103");
    assert_eq!(
        failure.context().get("population_id").map(String::as_str),
        Some("P1")
    );
    assert!(
        failure.context().get("detail").is_some_and(
            |detail| detail.contains("tolerance=0.0625") && detail.contains("0.25 copias")
        ),
        "el diagnóstico nombra la tolerancia efectiva por haplotipo y en copias: {:?}",
        failure.context()
    );
}

/// Una cuota imposible con la cota predeterminada se rechaza con
/// `GPS-DOMAIN-103`, sin corregirla ni relajarla.
///
/// Con cinco haplotipos equiprobables y dos copias, dar una copia a un rango
/// deja `|1/2-1/5| = 0,3` frente a `0,5/C = 0,25`, y no dársela deja
/// `|0-1/5| = 0,2`: la combinación no existe, y el grupo fronterizo lo
/// demuestra comprobando sus dos estados.
#[test]
fn an_impossible_quota_is_rejected_with_its_own_diagnostic() {
    let root = temporary("impossible-quota");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1).sampling("controlled")],
        ..five_haplotype_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let failure = assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
        .expect_err("cinco haplotipos no caben en dos copias dentro de 0,5/C");
    assert_eq!(failure.code(), "GPS-DOMAIN-103");
    assert_eq!(
        failure.context().get("population_id").map(String::as_str),
        Some("P1")
    );
    assert!(
        failure
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("block_id=")),
        "el detalle nombra el bloque, que el contexto del catálogo no lleva: {:?}",
        failure.context()
    );
}

/// Con `C` suficiente, los mismos cinco haplotipos sí reparten cuota exacta y
/// su grupo de restos iguales queda **incluido por completo**: no hay
/// frontera empatada y el desempate no llega a abrirse.
#[test]
fn a_group_included_in_full_never_opens_a_tie() {
    let root = temporary("no-tie");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 5).sampling("controlled")],
        ..five_haplotype_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        assert_eq!(quota.counts(), &[2, 2, 2, 2, 2], "diez copias entre cinco");
        assert!(
            quota.tie().is_none(),
            "un grupo incluido por completo no cruza ninguna frontera"
        );
    }
}

/// El desempate del grupo fronterizo se reproduce desde su ámbito publicado:
/// una única Fisher–Yates sobre los rangos empatados en orden ascendente y el
/// prefijo necesario.
#[test]
fn the_boundary_tie_is_reproducible_from_its_published_scope() {
    let root = temporary("tie-replay");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2).sampling("controlled")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");

    let mut ties = 0_u64;
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        let tie = quota.tie().expect("la frontera está empatada");
        assert_eq!(
            tie.candidate_ranks(),
            &[2, 3],
            "los dos rangos de peso 1/8 empatan"
        );
        assert_eq!(tie.select_count(), 1, "el residual solo alcanza a uno");
        assert_eq!(tie.remainder(), 0.5, "ambos restos valen 0,5");

        let mut candidates = tie.candidate_ranks().to_vec();
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::controlled_quota_tie(
                &scenario_id,
                unit.population_id(),
                unit.block_id(),
            ),
        )
        .expect("subsemilla")
        .local_rng();
        // Fisher–Yates literal: `i = len-1..1` con `j` en `0..=i`.
        for index in (1..candidates.len()).rev() {
            let other =
                usize::try_from(f5_uniform_index(&mut rng, (index + 1) as u64).expect("entero"))
                    .expect("índice");
            candidates.swap(index, other);
        }
        assert_eq!(
            tie.selected_ranks(),
            &candidates[..1],
            "el prefijo elegido no reproduce su sorteo"
        );
        let selected = tie.selected_ranks()[0] as usize;
        let excluded = if selected == 2 { 3 } else { 2 };
        assert_eq!(
            quota.counts()[selected],
            1,
            "el ganador recibe piso más uno"
        );
        assert_eq!(quota.counts()[excluded], 0, "el perdedor se queda en piso");
        ties += 1;
    }
    assert_eq!(ties, 2, "los dos bloques resuelven su propia frontera");
}

/// La permutación de copias se reproduce desde `founder_controlled_permutation`
/// sobre las `C` identidades ordenadas por `individual_id` y `copy_index`.
#[test]
fn the_copy_permutation_is_reproducible_and_assigns_the_sorted_multiset() {
    let root = temporary("permutation");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 4).sampling("controlled")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");

    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        let mut seats: Vec<(String, u32)> = unit
            .assignments()
            .iter()
            .map(|copy| (copy.individual_id().as_str().to_owned(), copy.copy_index()))
            .collect();
        seats.sort();
        let mut permuted = seats.clone();
        let mut rng = derive_f5_subseed(
            root_seed,
            F5RandomDecision::founder_controlled_permutation(
                &scenario_id,
                unit.population_id(),
                unit.block_id(),
            ),
        )
        .expect("subsemilla")
        .local_rng();
        for index in (1..permuted.len()).rev() {
            let other =
                usize::try_from(f5_uniform_index(&mut rng, (index + 1) as u64).expect("entero"))
                    .expect("índice");
            permuted.swap(index, other);
        }
        let mut multiset = Vec::new();
        for (rank, count) in quota.counts().iter().enumerate() {
            for _ in 0..*count {
                multiset.push(rank as u64);
            }
        }
        assert!(
            multiset.windows(2).all(|pair| pair[0] <= pair[1]),
            "el multiconjunto se ordena por haplotype_rank"
        );
        let expected: BTreeMap<(String, u32), u64> = permuted
            .into_iter()
            .zip(multiset)
            .collect::<BTreeMap<_, _>>();
        for copy in unit.assignments() {
            let key = (copy.individual_id().as_str().to_owned(), copy.copy_index());
            assert_eq!(
                copy.haplotype_rank(),
                expected[&key],
                "la copia {key:?} no recibió la posición que le tocaba"
            );
        }
    }
}

/// Regresión del bloqueante de la revisión F5.3.1: tres fundadores haploides y
/// dos haplotipos equiprobables.
///
/// El error exacto vale `|1/3 - 1/2| = 1/6` y la cota `0,5/3` vale también
/// `1/6`: la cota inclusiva se cumple y `[2,1]` es la única cuota razonable.
/// Evaluada como `|n/C - q|` en binary64 daba `0.16666666666666669` frente a
/// `0.16666666666666666` y se rechazaba por un ULP.
#[test]
fn the_minimal_haploid_case_meets_its_bound_exactly() {
    let root = temporary("haploid-minimal");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 3).sampling("controlled")],
        ploidy: 1,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        assert_eq!(quota.copies(), 3);
        assert_eq!(
            quota.counts().iter().sum::<u64>(),
            3,
            "la cuota debe sumar C"
        );
        let mut sorted = quota.counts().to_vec();
        sorted.sort_unstable();
        assert_eq!(sorted, vec![1, 2], "tres copias entre dos haplotipos");
        assert_largest_remainder(quota.counts(), &fixture.weights, 3, (1, 6));
    }
}

/// El mismo bloqueante por la vía de las cuotas diádicas: `C=6` reparte
/// `[3,1,1,1]` y ninguno de sus cuatro rangos se sale de la cota.
#[test]
fn six_copies_split_the_dyadic_profile_without_rejection() {
    let root = temporary("six-copies");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 3).sampling("controlled")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        assert_eq!(
            quota.counts(),
            &[3, 1, 1, 1],
            "el reparto exacto de 6 copias entre 1/2, 1/4, 1/8 y 1/8"
        );
        assert!(
            quota.tie().is_none(),
            "los restos 0, 0,5, 0,75 y 0,75 no dejan frontera partida"
        );
        assert_largest_remainder(quota.counts(), &fixture.weights, 6, (1, 12));
    }
}

/// Toda cardinalidad legal de una subpoblación diádica reparte cuota: la
/// aceptación no puede depender de que `C` sea potencia de dos.
#[test]
fn every_legal_copy_count_accepts_its_quota() {
    for individuals in 1_u64..=8 {
        let root = temporary(&format!("copy-count-{individuals}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", individuals).sampling("controlled")],
            ..dyadic_fixture()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
            .unwrap_or_else(|error| {
                panic!("C={} rechazada: {:?}", individuals * 2, error.context())
            });
        let copies = individuals * 2;
        for unit in assignment.units() {
            let quota = unit.quota().expect("cuota controlada");
            assert_eq!(quota.counts().iter().sum::<u64>(), copies, "C={copies}");
            assert_largest_remainder(
                quota.counts(),
                &fixture.weights,
                copies,
                (1, u128::from(2 * copies)),
            );
        }
    }
}

/// `F>0` es incompatible con el modo controlado: la combinación se rechaza
/// antes de admitir el escenario, tanto en el predeterminado como en la
/// sustitución.
#[test]
fn a_positive_f_is_incompatible_with_the_controlled_mode() {
    for (label, subpopulations, overrides) in [
        (
            "default",
            vec![
                Subpopulation::new("P1", 1)
                    .sampling("controlled")
                    .inbreeding("0.25"),
            ],
            Vec::new(),
        ),
        (
            "override",
            vec![Subpopulation::new("P1", 1).inbreeding("0.25")],
            vec![Substitution {
                population: "P1",
                block_index: 0,
                mode: "controlled",
                tolerance: None,
            }],
        ),
    ] {
        let root = temporary(&format!("incompatible-f-{label}"));
        let fixture = Fixture {
            subpopulations,
            overrides,
            ..Fixture::default()
        };
        let general = fixture.write(&root);
        let failure =
            load_and_validate_f5_config(&general).expect_err("controlado con F>0 es incompatible");
        assert!(
            failure
                .context()
                .get("detail")
                .is_some_and(|detail| detail.contains("F efectivo igual a cero")),
            "{label}: {:?}",
            failure.context()
        );
        let _ = fs::remove_dir_all(&root);
    }
}

/// `F=0` declarado explícitamente sí es compatible con el controlado.
#[test]
fn a_zero_f_is_compatible_with_the_controlled_mode() {
    let root = temporary("zero-f");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 1)
                .sampling("controlled")
                .inbreeding("0.0"),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    assert!(
        assignment
            .units()
            .iter()
            .all(|unit| unit.mode() == F5SamplingMode::Controlled),
        "F=0 no impide el modo controlado"
    );
}

// ---------------------------------------------------------------------------
// Identidad, orden y reproducibilidad
// ---------------------------------------------------------------------------

/// Las asignaciones se publican ordenadas por `individual_id` y `copy_index`,
/// con una fila exacta por copia y sin repetir ninguna.
#[test]
fn assignments_are_canonical_and_cover_every_copy_once() {
    let root = temporary("canonical-order");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 3).sampling("controlled"),
            Subpopulation::new("P2", 2),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for unit in assignment.units() {
        let keys: Vec<(String, u32)> = unit
            .assignments()
            .iter()
            .map(|copy| (copy.individual_id().as_str().to_owned(), copy.copy_index()))
            .collect();
        let mut sorted = keys.clone();
        sorted.sort();
        assert_eq!(keys, sorted, "las copias no están en orden canónico");
        assert_eq!(
            keys.iter().collect::<BTreeSet<_>>().len(),
            keys.len(),
            "una copia aparece dos veces"
        );
    }
    let rosters = assignment.rosters();
    assert_eq!(rosters.len(), 2);
    assert_eq!(rosters[0].population_id().as_str(), "P1");
    assert_eq!(rosters[0].copies(), 6);
    assert_eq!(rosters[1].copies(), 4);
    assert!(
        rosters[0]
            .individuals()
            .iter()
            .all(|individual| individual.founder_key().is_none()),
        "los fundadores aportados por ID no llevan clave"
    );
}

/// La asignación es reproducible bajo la misma semilla y cambia con otra.
#[test]
fn the_assignment_is_reproducible_under_its_seed() {
    let fixture = |seed: &'static str| Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 4).sampling("controlled"),
            Subpopulation::new("P2", 4),
        ],
        seed,
        ..dyadic_fixture()
    };
    let first = temporary("seed-a");
    let second = temporary("seed-a-again");
    let third = temporary("seed-b");
    let ranks = |chain: &Chain| {
        assign(chain)
            .units()
            .iter()
            .flat_map(|unit| {
                unit.assignments()
                    .iter()
                    .map(|copy| copy.haplotype_rank())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>()
    };
    let alternate = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef";
    let left = ranks(&build_chain(&first, &fixture(SEED)));
    let right = ranks(&build_chain(&second, &fixture(SEED)));
    let other = ranks(&build_chain(&third, &fixture(alternate)));
    assert_eq!(left, right, "la misma semilla debe producir lo mismo");
    assert_ne!(
        left, other,
        "otra semilla debe mover al menos una asignación"
    );
}

/// La asignación exige la admisión sellada de su propia configuración.
#[test]
fn a_foreign_sealed_admission_cannot_assign_founders() {
    let first = temporary("binding-a");
    let second = temporary("binding-b");
    let left = build_chain(&first, &Fixture::default());
    let right = build_chain(
        &second,
        &Fixture {
            block_count: 3,
            slack: 2,
            ..Fixture::default()
        },
    );
    let failure = assign_f5_founders(&left.config, &right.sealed, &left.search)
        .expect_err("la admisión sellada no corresponde");
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("preflight.binding")
    );
}

/// Los fundadores declarados por clave derivan su `IndividualId` con el
/// contrato población×clave, y el censo conserva ambas cosas.
#[test]
fn founder_keys_derive_their_individual_identity() {
    let root = temporary("founder-keys");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2).keys()],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let roster = &assignment.rosters()[0];
    assert_eq!(roster.individuals().len(), 2);
    for individual in roster.individuals() {
        assert!(
            individual.individual_id().as_str().starts_with("IND_"),
            "un ID derivado lleva su prefijo de clase: {}",
            individual.individual_id()
        );
        assert!(
            individual.founder_key().is_some(),
            "el censo conserva la clave que derivó la identidad"
        );
    }
    let identities: BTreeSet<&str> = roster
        .individuals()
        .iter()
        .map(|individual| individual.individual_id().as_str())
        .collect();
    assert_eq!(identities.len(), 2, "dos claves derivan dos identidades");
    for unit in assignment.units() {
        for copy in unit.assignments() {
            assert!(
                identities.contains(copy.individual_id().as_str()),
                "una copia nombra un individuo ajeno al censo"
            );
        }
    }
}

/// Una tolerancia declarada `0` admite exactamente la cuota sin error y
/// rechaza cualquier otra.
///
/// Con `C=8` los cuatro ideales del perfil diádico son enteros y el error es
/// nulo; con `C=4` los dos rangos de peso `1/8` valen media copia y ninguna
/// cuota entera puede tener error cero.
#[test]
fn a_declared_zero_tolerance_admits_only_an_exact_quota() {
    let exact = temporary("zero-tolerance-exact");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 4).sampling("controlled")],
        overrides: vec![Substitution {
            population: "P1",
            block_index: 0,
            mode: "controlled",
            tolerance: Some("0.0"),
        }],
        ..dyadic_fixture()
    };
    let chain = build_chain(&exact, &fixture);
    let assignment = assign(&chain);
    let row = assignment
        .sampling()
        .iter()
        .find(|row| row.block_index() == 0)
        .expect("bloque 0");
    assert_eq!(row.controlled_max_abs_error_declared(), Some(0.0));
    assert_eq!(row.controlled_max_abs_error_effective(), Some(0.0));
    let unit = assignment
        .unit(row.population_id(), 0)
        .expect("unidad del bloque 0");
    assert_eq!(
        unit.quota().expect("cuota").counts(),
        &[4, 2, 1, 1],
        "la única cuota con error nulo para C=8"
    );

    let impossible = temporary("zero-tolerance-impossible");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2).sampling("controlled")],
        overrides: vec![Substitution {
            population: "P1",
            block_index: 0,
            mode: "controlled",
            tolerance: Some("0.0"),
        }],
        ..dyadic_fixture()
    };
    let chain = build_chain(&impossible, &fixture);
    let failure = assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
        .expect_err("media copia no cabe en una tolerancia nula");
    assert_eq!(failure.code(), "GPS-DOMAIN-103");
    assert!(
        failure
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("tolerance=0;")),
        "el diagnóstico nombra la tolerancia declarada: {:?}",
        failure.context()
    );
}

/// Las claves fundadoras también derivan identidad en modo controlado, donde
/// además ordenan la permutación de copias.
#[test]
fn founder_keys_also_work_in_controlled_mode() {
    let root = temporary("keys-controlled");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 4).sampling("controlled").keys()],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let roster = &assignment.rosters()[0];
    assert!(
        roster
            .individuals()
            .iter()
            .all(|individual| individual.founder_key().is_some()),
        "el censo conserva las claves"
    );
    let identities: Vec<&str> = roster
        .individuals()
        .iter()
        .map(|individual| individual.individual_id().as_str())
        .collect();
    let mut sorted = identities.clone();
    sorted.sort_unstable();
    assert_eq!(identities, sorted, "el censo se publica ordenado por ID");
    for unit in assignment.units() {
        let quota = unit.quota().expect("cuota controlada");
        assert_eq!(quota.counts(), &[4, 2, 1, 1]);
        let keys: Vec<(String, u32)> = unit
            .assignments()
            .iter()
            .map(|copy| (copy.individual_id().as_str().to_owned(), copy.copy_index()))
            .collect();
        let mut expected = keys.clone();
        expected.sort();
        assert_eq!(keys, expected, "las copias salen en orden canónico");
    }
}

/// El coste contado crece **linealmente** con `C`: ni el reparto ni la
/// asignación ordenan las copias, de modo que la razón por copia no puede
/// crecer con `log C`.
///
/// La medida se contabiliza dentro de cada bucle del productor —rangos
/// visitados, comparaciones de la selección parcial, sorteos de las dos
/// permutaciones y copias escritas—, no junto a él.
#[test]
fn the_counted_cost_stays_linear_in_the_copy_count() {
    let mut measured: Vec<(u64, u64)> = Vec::new();
    for individuals in [1_u64, 4, 16, 64, 256] {
        for mode in ["controlled", "probabilistic"] {
            let root = temporary(&format!("cost-{mode}-{individuals}"));
            let fixture = Fixture {
                subpopulations: vec![Subpopulation::new("P1", individuals).sampling(mode)],
                ..dyadic_fixture()
            };
            let chain = build_chain(&root, &fixture);
            let assignment = assign(&chain);
            let copies = individuals * 2;
            for unit in assignment.units() {
                assert!(
                    unit.operations() <= 8 * (copies + 4),
                    "{mode} con C={copies}: {} visitas superan la cota lineal",
                    unit.operations()
                );
                if mode == "controlled" {
                    measured.push((copies, unit.operations()));
                }
            }
        }
    }
    // El coste marginal por copia es constante: si volviera una ordenación de
    // las `C` asignaciones, crecería con `log2(C)`.
    let smallest = measured.first().copied().expect("medida mínima");
    let largest = measured.last().copied().expect("medida máxima");
    let marginal = (largest.1 - smallest.1) as f64 / (largest.0 - smallest.0) as f64;
    assert!(
        marginal <= 4.0,
        "el coste marginal por copia ({marginal}) debe ser constante y pequeño"
    );
}

/// F5.3.1 no adelanta nada de F5.3.2: ninguna ruta del módulo nombra el ámbito
/// `founder_f`, de modo que no puede consumirlo.
#[test]
fn the_founder_f_scope_is_never_opened_in_this_milestone() {
    let source = include_str!("../src/f5_founders.rs");
    assert!(
        !source.contains("founder_f("),
        "el muestreo de F es de F5.3.2 y este módulo no puede abrir su ámbito"
    );
    assert!(
        source.contains("F5RandomDecision::founder_sample"),
        "la búsqueda negativa solo vale si el módulo sí abre los ámbitos suyos"
    );
}
