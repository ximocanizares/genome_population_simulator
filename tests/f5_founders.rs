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
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::f5_founders::{
    F5BlockAssignment, F5FounderAssignment, F5FounderCopy, F5FounderIdentities, F5FounderUnit,
    F5RealizedFStatus, F5SamplingMode, F5SamplingSource, F5UnitHaplotypes, assign_f5_founders,
    f5_haplotype_alleles, materialize_f5_founder_identities, resolve_f5_founder_roster,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_search::{F5HaplotypeSearch, search_f5_haplotypes};
use genome_population_simulator::identity::{
    ExecutionId, F5BlockId, F5HaplotypeId, FounderCopyId, IndividualId, PopulationId, ScenarioId,
};
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
    /// Identificadores fundadores literales; sin ellos se generan `<id>-I<n>`.
    founder_ids: Option<&'a [&'a str]>,
}

impl<'a> Subpopulation<'a> {
    const fn new(id: &'a str, individual_count: u64) -> Self {
        Self {
            id,
            individual_count,
            sampling: None,
            inbreeding_f: None,
            keys: false,
            founder_ids: None,
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

    /// Fija los identificadores fundadores literales de la subpoblación.
    const fn founders(mut self, founder_ids: &'a [&'a str]) -> Self {
        self.founder_ids = Some(founder_ids);
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
    /// Secuencias de la referencia sintética; `block_count` debe ser múltiplo
    /// suyo, de modo que cada secuencia reciba exactamente los mismos bloques.
    sequences: u64,
    weights: Vec<&'a str>,
    gap: u64,
    spectrum: &'a str,
    bins: Vec<(u64, u64, f64, f64)>,
    max_haplotype_attempts: u64,
    max_swap_proposals: u64,
    ploidy: u64,
    seed: &'a str,
    slack: u64,
    /// `[f5] inbreeding_f`: el F **de población** de DTC-089, que cada
    /// subpoblación puede sustituir con el suyo.
    scenario_inbreeding: Option<&'a str>,
}

impl Default for Fixture<'_> {
    fn default() -> Self {
        Self {
            subpopulations: vec![Subpopulation::new("P1", 1)],
            overrides: Vec::new(),
            block_count: 2,
            block_length: 3,
            sequences: 1,
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
            scenario_inbreeding: None,
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
        assert!(
            self.sequences >= 1 && self.block_count.is_multiple_of(self.sequences),
            "cada secuencia debe recibir el mismo número entero de bloques"
        );
        let span =
            usize::try_from(self.block_length * (self.block_count / self.sequences) * self.slack)
                .expect("span");
        let mut fasta = String::new();
        for ordinal in 0..self.sequences {
            // Cada secuencia empieza en una base distinta: dos secuencias con
            // el mismo contenido seguirían siendo dos secuencias, pero así el
            // fixture no puede confundirlas por accidente.
            let sequence = (0..span)
                .map(|index| char::from(bases[(index + ordinal as usize) % 4]))
                .collect::<String>();
            fasta.push_str(&format!(">chr{}\n{sequence}\n", ordinal + 1));
        }
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
            let founders = subpopulation.founder_ids.map_or_else(
                || {
                    (0..subpopulation.individual_count)
                        .map(|index| format!("\"{}-I{index}\"", subpopulation.id))
                        .collect::<Vec<_>>()
                        .join(", ")
                },
                |ids| {
                    ids.iter()
                        .map(|id| format!("\"{id}\""))
                        .collect::<Vec<_>>()
                        .join(", ")
                },
            );
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
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {founder_count}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = {}\nmax_swap_proposals = {}\nmax_haplotype_cells = 65536\nblock_count = {}\n{}{subpopulations}{overrides}",
                sha256(metadata.as_bytes()),
                self.max_haplotype_attempts,
                self.max_swap_proposals,
                self.block_count,
                self.scenario_inbreeding
                    .map_or(String::new(), |value| format!("inbreeding_f = {value}\n")),
            ),
        )
        .expect("scenario");
        general
    }
}

struct Chain {
    config: EffectiveF5Config,
    sealed: SealedF5Admission,
    layout: F5PhysicalLayout,
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
        layout,
        search,
        root: root.to_path_buf(),
    }
}

fn assign(chain: &Chain) -> F5FounderAssignment {
    assign_f5_founders(&chain.config, &chain.sealed, &chain.search)
        .expect("asignación fundadora completa")
}

/// Fila haplotípica tal como `f5-haplotypes.tsv` la publicaría.
///
/// El repositorio conserva el `HAP_` y no la cadena de alelos; quien escribe el
/// TSV la recompone desde la matriz ganadora, que es la misma de la que salió
/// el material del identificador. Esta batería hace lo mismo.
#[derive(Clone, Debug, Eq, PartialEq)]
struct PublishedHaplotype {
    population_id: PopulationId,
    block_id: F5BlockId,
    block_index: u64,
    haplotype_rank: u64,
    alleles: String,
    haplotype_id: F5HaplotypeId,
}

impl PublishedHaplotype {
    const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }
    const fn block_id(&self) -> &F5BlockId {
        &self.block_id
    }
    const fn block_index(&self) -> u64 {
        self.block_index
    }
    const fn haplotype_rank(&self) -> u64 {
        self.haplotype_rank
    }
    fn alleles(&self) -> &str {
        &self.alleles
    }
    const fn haplotype_id(&self) -> &F5HaplotypeId {
        &self.haplotype_id
    }
}

/// El repositorio más lo que un escritor habría recogido de sus dos recorridos.
///
/// Las pruebas sí pueden retener las `B*C_total` asignaciones: sus fixtures son
/// diminutos. El productor no, y por eso `for_each_assignment` las entrega una
/// a una. Este tipo concentra esa diferencia en un solo sitio.
struct Materialized {
    identities: F5FounderIdentities,
    haplotypes: Vec<PublishedHaplotype>,
    assignments: Vec<F5BlockAssignment>,
}

impl Materialized {
    fn copies(&self) -> &[F5FounderCopy] {
        self.identities.copies()
    }
    fn units(&self) -> &[F5UnitHaplotypes] {
        self.identities.units()
    }
    fn haplotypes(&self) -> &[PublishedHaplotype] {
        &self.haplotypes
    }
    fn assignments(&self) -> &[F5BlockAssignment] {
        &self.assignments
    }
    fn copy(&self, founder_copy_id: &FounderCopyId) -> Option<&F5FounderCopy> {
        self.identities.copy(founder_copy_id)
    }
}

/// Recompone las filas de `f5-haplotypes.tsv` desde el repositorio y la matriz.
fn published_haplotypes(
    chain: &Chain,
    identities: &F5FounderIdentities,
) -> Vec<PublishedHaplotype> {
    let mut rows = Vec::new();
    for unit in identities.units() {
        let matrix = chain
            .search
            .block(unit.population_id(), unit.block_index())
            .expect("bloque ganador")
            .matrix();
        for (rank, haplotype_id) in unit.haplotype_ids().iter().enumerate() {
            rows.push(PublishedHaplotype {
                population_id: unit.population_id().clone(),
                block_id: unit.block_id().clone(),
                block_index: unit.block_index(),
                haplotype_rank: rank as u64,
                alleles: f5_haplotype_alleles(matrix, rank).expect("alelos de la fila"),
                haplotype_id: haplotype_id.clone(),
            });
        }
    }
    rows
}

/// Identidades fundadoras de la cadena completa, con la identidad de ejecución
/// que la produjo y sus dos recorridos ya recogidos.
fn materialize(chain: &Chain) -> Materialized {
    let assignment = assign(chain);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let mut identities = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &execution,
        &chain.layout,
        &chain.search,
        &assignment,
    )
    .expect("identidades fundadoras materializadas");
    let haplotypes = published_haplotypes(chain, &identities);
    let mut assignments = Vec::new();
    identities
        .for_each_assignment(&assignment, |row| {
            assignments.push(row.clone());
            Ok(())
        })
        .expect("recorrido de asignaciones");
    Materialized {
        identities,
        haplotypes,
        assignments,
    }
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

/// Pesos efectivos de un perfil de la batería.
///
/// Los perfiles de estas pruebas suman exactamente 1 en binary64 con la suma
/// compensada, de modo que la normalización DTC-049 los divide por `1.0` y el
/// peso efectivo es el declarado. La aserción lo comprueba en cada llamada: un
/// perfil que no lo cumpliera exigiría replicar aquí la normalización, y el
/// oráculo dejaría de ser independiente.
fn effective_weights(declared: &[&str]) -> Vec<f64> {
    let values: Vec<f64> = declared
        .iter()
        .map(|weight| weight.parse::<f64>().expect("peso decimal"))
        .collect();
    let mut compensated = NeumaierOracle::default();
    let mut sum = 0.0_f64;
    for value in &values {
        sum = compensated.add(*value);
    }
    assert_eq!(
        sum, 1.0,
        "los pesos declarados deben sumar exactamente 1 en binary64"
    );
    values
}

/// Un `f64` finito y no negativo como racional **exacto**.
///
/// Todo binary64 es un racional diádico: mantisa por una potencia de dos. Es lo
/// que permite al oráculo trabajar sin perder un bit con pesos que **no** son
/// diádicos en decimal, como `0,3`, cuya representación sí lo es en binario.
fn exact_rational(value: f64) -> (u128, u128) {
    assert!(value.is_finite() && value >= 0.0, "{value} no es admisible");
    if value == 0.0 {
        return (0, 1);
    }
    let bits = value.to_bits();
    let raw_exponent = ((bits >> 52) & 0x7ff) as i32;
    let raw_mantissa = bits & ((1_u64 << 52) - 1);
    let (mantissa, exponent) = if raw_exponent == 0 {
        (u128::from(raw_mantissa), -1074_i32)
    } else {
        (
            u128::from(raw_mantissa | (1_u64 << 52)),
            raw_exponent - 1075,
        )
    };
    if exponent >= 0 {
        assert!(exponent < 64, "el oráculo acota el numerador");
        (mantissa << exponent, 1)
    } else {
        let shift = u32::try_from(-exponent).expect("exponente");
        assert!(shift < 100, "el oráculo acota el denominador a 2^100");
        (mantissa, 1_u128 << shift)
    }
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

/// Comprueba en racionales exactos todo lo que el reparto por mayores restos
/// debe cumplir, sin reproducir ningún sorteo:
///
/// 1. cada `n_h` es el piso o el piso más uno del ideal normativo `fl(C*q_h)`;
/// 2. la suma es exactamente `C`;
/// 3. el número de incrementos es el residual `C - Σ pisos`, contado desde los
///    pisos y no desde los restos: con ideales `fl(C*q_h)` la suma de los
///    ideales no vale exactamente `C` —`fl(7*0,7)+fl(7*0,3)` se queda un ulp
///    corto— y solo los pisos son enteros exactos;
/// 4. ningún rango incrementado tiene un resto menor que uno no incrementado;
/// 5. cada `n_h` respeta la cota **en copias** `|n_h - fl(C*q_h)| <= t*C`.
///
/// El ideal es `fl(C*q_h)` —el producto en binary64— porque es el que DTC-089
/// declara normativo, y la cota se evalúa en copias por el mismo motivo. El
/// oráculo no vuelve a multiplicar en `f64`: convierte ese producto a su
/// racional exacto y opera en `u128`, de modo que acepta pesos que no son
/// diádicos en decimal —`0,3`, `0,7`— sin perder un bit. `tolerance` es la
/// tolerancia **en copias**: `(1, 2)` es el `0,5` exacto predeterminado.
fn assert_largest_remainder(
    counts: &[u64],
    declared: &[&str],
    copies: u64,
    tolerance: (u128, u128),
) {
    let weights = effective_weights(declared);
    let (tolerance_numerator, tolerance_denominator) = tolerance;
    assert_eq!(counts.len(), weights.len(), "un recuento por rango");
    let copies_f64 = copies as f64;
    // Ideal normativo de cada rango: `fl(C*q_h)`, ya como racional exacto.
    let ideals: Vec<(u128, u128)> = weights
        .iter()
        .map(|weight| exact_rational(copies_f64 * *weight))
        .collect();
    let common = ideals
        .iter()
        .map(|(_, denominator)| *denominator)
        .max()
        .expect("al menos un rango");
    let copies = u128::from(copies);
    let mut floors_total = 0_u128;
    let mut increments = 0_u128;
    let mut incremented_remainders: Vec<u128> = Vec::new();
    let mut untouched_remainders: Vec<u128> = Vec::new();
    let mut total = 0_u128;
    for (rank, (numerator, denominator)) in ideals.iter().enumerate() {
        assert_eq!(common % denominator, 0, "denominadores potencia de dos");
        let floor = numerator / denominator;
        // Restos llevados al denominador común, para poder compararlos y
        // sumarlos sin dividir nunca.
        let remainder = (numerator % denominator) * (common / denominator);
        floors_total += floor;
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

        // `|n - num/den| <= tn/td`  ⟺  `|n*den - num| * td <= tn * den`.
        let left = (count * denominator) as i128 - *numerator as i128;
        assert!(
            left.unsigned_abs() * tolerance_denominator <= tolerance_numerator * denominator,
            "rango {rank}: la cuota {count} incumple la cota exacta en copias"
        );
    }
    assert_eq!(total, copies, "la cuota no suma C");
    assert!(floors_total <= copies, "los pisos no pueden superar C");
    assert_eq!(
        increments,
        copies - floors_total,
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
            // Cota efectiva predeterminada: `0,5` exacto **en copias**.
            assert_largest_remainder(quota.counts(), &weights, copies, (1, 2));
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
        assert_largest_remainder(quota.counts(), &fixture.weights, 3, (1, 2));
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
        assert_largest_remainder(quota.counts(), &fixture.weights, 6, (1, 2));
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
            assert_largest_remainder(quota.counts(), &fixture.weights, copies, (1, 2));
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

// ---------------------------------------------------------------------------
// F diploide y F realizado
// ---------------------------------------------------------------------------

/// Rango que el ámbito `founder_sample` de una copia elige por su cuenta.
///
/// Recorre los pesos secuencialmente con la suma compensada de DTC-089, igual
/// que el oráculo del muestreo sin F, y no llama a ninguna función del
/// productor.
fn independent_rank(
    weights: &[f64],
    root_seed: RootSeed,
    scenario_id: &ScenarioId,
    unit: &F5FounderUnit,
    individual_id: &IndividualId,
    copy_index: u32,
) -> u64 {
    let mut rng = derive_f5_subseed(
        root_seed,
        F5RandomDecision::founder_sample(
            scenario_id,
            unit.population_id(),
            unit.block_id(),
            individual_id,
            copy_index,
        ),
    )
    .expect("subsemilla")
    .local_rng();
    let threshold = (rng.next_u64() >> 11) as f64 / (1_u64 << 53) as f64;
    let mut prefix = NeumaierOracle::default();
    let mut chosen = weights
        .iter()
        .rposition(|weight| *weight > 0.0)
        .expect("algún peso positivo");
    for (rank, weight) in weights.iter().enumerate() {
        if threshold < prefix.add(*weight) {
            chosen = rank;
            break;
        }
    }
    chosen as u64
}

/// Decisión de la Bernoulli de F: `bits53 < floor(F*2^53)` sobre el primer
/// U53 del ámbito `founder_f` del individuo.
fn reuses_copy_zero(
    root_seed: RootSeed,
    scenario_id: &ScenarioId,
    unit: &F5FounderUnit,
    individual_id: &IndividualId,
    inbreeding: f64,
) -> bool {
    let mut rng = derive_f5_subseed(
        root_seed,
        F5RandomDecision::founder_f(
            scenario_id,
            unit.population_id(),
            unit.block_id(),
            individual_id,
        ),
    )
    .expect("subsemilla")
    .local_rng();
    (rng.next_u64() >> 11) < (inbreeding * (1_u64 << 53) as f64).floor() as u64
}

/// Individuos de una unidad diploide como parejas `(copia 0, copia 1)`.
fn diploid_pairs(unit: &F5FounderUnit) -> Vec<(IndividualId, u64, u64)> {
    unit.assignments()
        .chunks(2)
        .map(|pair| {
            assert_eq!(
                pair.len(),
                2,
                "un diploide publica dos copias por individuo"
            );
            assert_eq!(
                pair[0].individual_id(),
                pair[1].individual_id(),
                "las copias de un individuo se publican juntas"
            );
            assert_eq!((pair[0].copy_index(), pair[1].copy_index()), (0, 1));
            (
                pair[0].individual_id().clone(),
                pair[0].haplotype_rank(),
                pair[1].haplotype_rank(),
            )
        })
        .collect()
}

/// `F_realized` exacto en racionales de 128 bits, escrito desde su definición
/// y con una expresión distinta de la del productor:
/// `1 - het*C² / (N*(C² - Σ n_h²))`, que es `1 - H_obs/H_exp` sin dividir
/// nunca antes de tiempo. Devuelve `None` cuando `H_expected` es exactamente
/// cero, que es el caso `not_applicable` del contrato.
///
/// `frequencies` decide si se usan las frecuencias **realizadas** o el
/// objetivo `q` del perfil: con el objetivo, el valor es otro, y eso es lo que
/// hace discriminante la comparación.
fn realized_f_oracle(unit: &F5FounderUnit, haplotype_count: usize) -> Option<f64> {
    let mut counts = vec![0_i128; haplotype_count];
    let mut heterozygous = 0_i128;
    let mut individuals = 0_i128;
    for (_, left, right) in diploid_pairs(unit) {
        counts[left as usize] += 1;
        counts[right as usize] += 1;
        individuals += 1;
        if left != right {
            heterozygous += 1;
        }
    }
    let copies = 2 * individuals;
    let squares: i128 = counts.iter().map(|count| count * count).sum();
    let expected_numerator = copies * copies - squares;
    if expected_numerator == 0 {
        return None;
    }
    Some(1.0 - (heterozygous * copies * copies) as f64 / (individuals * expected_numerator) as f64)
}

/// El mismo cociente con el objetivo `q` del perfil en lugar de la frecuencia
/// realizada, que es la variante que DTC-089 prohíbe.
fn target_frequency_variant(unit: &F5FounderUnit, weights: &[f64]) -> f64 {
    let pairs = diploid_pairs(unit);
    let heterozygous = pairs
        .iter()
        .filter(|(_, left, right)| left != right)
        .count();
    let h_observed = heterozygous as f64 / pairs.len() as f64;
    let h_expected = 1.0 - weights.iter().map(|weight| weight * weight).sum::<f64>();
    1.0 - h_observed / h_expected
}

/// `F=0` muestrea cada copia por su cuenta y coincide exactamente con el HWE
/// del mismo escenario sin F declarado.
#[test]
fn a_zero_f_reproduces_hardy_weinberg_copy_by_copy() {
    let root = temporary("f-zero");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8).inbreeding("0.0")],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();

    let mut checked = 0_u64;
    for unit in assignment.units() {
        for copy in unit.assignments() {
            assert_eq!(
                copy.haplotype_rank(),
                independent_rank(
                    &weights,
                    root_seed,
                    &scenario_id,
                    unit,
                    copy.individual_id(),
                    copy.copy_index()
                ),
                "con F=0 la copia {}#{} debe muestrearse por su cuenta",
                copy.individual_id(),
                copy.copy_index()
            );
            checked += 1;
        }
    }
    assert_eq!(checked, 32, "dos bloques × ocho individuos diploides");

    // El mismo escenario sin declarar F produce exactamente las mismas filas.
    let bare_root = temporary("f-absent");
    let bare = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8)],
        ..Fixture::default()
    };
    let bare_chain = build_chain(&bare_root, &bare);
    let bare_assignment = assign(&bare_chain);
    assert_eq!(
        assignment.units(),
        bare_assignment.units(),
        "F=0 coincide con HWE en todas las unidades"
    );
    // Y el escenario debe tener heterocigotos: sin ellos la comparación con
    // F=1 no distinguiría nada.
    assert!(
        assignment
            .units()
            .iter()
            .flat_map(diploid_pairs)
            .any(|(_, left, right)| left != right),
        "el fixture debe producir heterocigotos con F=0"
    );
}

/// El F **de población** de `[f5]` gobierna una subpoblación que no declara el
/// suyo.
///
/// REQ-EVO-46 exige que F pueda definirse para la población completa; sin esta
/// prueba, borrar el respaldo `.or(...)` al F de `[f5]` no rompería nada y la
/// subpoblación pasaría a muestrear en HWE.
#[test]
fn the_population_level_f_governs_a_subpopulation_without_its_own() {
    let root = temporary("scenario-f");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8)],
        scenario_inbreeding: Some("1.0"),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();
    let mut differing = 0_u64;
    for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
        for (individual_id, first, second) in diploid_pairs(unit) {
            assert_eq!(first, second, "{individual_id} debe ser homocigoto con F=1");
            // Con F=0 la copia 1 habría sorteado esta otra fila: la diferencia
            // es lo que distingue «F de población aplicado» de «sin F».
            if independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 1) != first
            {
                differing += 1;
            }
        }
        assert_eq!(
            realized.target_f(),
            Some(1.0),
            "el objetivo publicado es el F de población"
        );
    }
    assert!(
        differing > 0,
        "sin el respaldo al F de población, alguna copia 1 sortearía otra fila"
    );
}

/// El F de la subpoblación **sustituye** al de población, no al revés.
///
/// Con la precedencia invertida, `[f5] inbreeding_f = 1.0` produciría
/// homocigosis completa donde el contrato pide HWE, sin diagnóstico y sin que
/// ninguna salida lo delatara.
#[test]
fn a_subpopulation_f_overrides_the_population_f() {
    let root = temporary("override-f");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8).inbreeding("0.0")],
        scenario_inbreeding: Some("1.0"),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();
    let mut heterozygous = 0_u64;
    for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
        assert_eq!(
            realized.target_f(),
            Some(0.0),
            "gana el F de la subpoblación"
        );
        for copy in unit.assignments() {
            assert_eq!(
                copy.haplotype_rank(),
                independent_rank(
                    &weights,
                    root_seed,
                    &scenario_id,
                    unit,
                    copy.individual_id(),
                    copy.copy_index()
                ),
                "con F=0 efectivo cada copia se muestrea por su cuenta"
            );
        }
        for (_, first, second) in diploid_pairs(unit) {
            if first != second {
                heterozygous += 1;
            }
        }
    }
    assert!(heterozygous > 0, "F=0 efectivo debe dejar heterocigotos");
}

/// Cada subpoblación recibe **su** F, no el de otra.
#[test]
fn two_subpopulations_receive_their_own_f() {
    let root = temporary("two-f");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 8).inbreeding("1.0"),
            Subpopulation::new("P2", 8).inbreeding("0.0"),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();
    let mut heterozygous_in_second = 0_u64;
    let mut homozygous_in_first = 0_u64;
    for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
        let first_population = unit.population_id().as_str() == "P1";
        assert_eq!(
            realized.target_f(),
            Some(if first_population { 1.0 } else { 0.0 }),
            "cada subpoblación publica su propio objetivo"
        );
        for (individual_id, first, second) in diploid_pairs(unit) {
            if first_population {
                assert_eq!(first, second, "{individual_id} de P1 tiene F=1");
                homozygous_in_first += 1;
            } else {
                assert_eq!(
                    second,
                    independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 1),
                    "{individual_id} de P2 tiene F=0"
                );
                if first != second {
                    heterozygous_in_second += 1;
                }
            }
        }
    }
    assert_eq!(homozygous_in_first, 16, "dos bloques × ocho individuos");
    assert!(
        heterozygous_in_second > 0,
        "P2 debe mostrar heterocigotos: sin ellos, tomar el F de P1 no se vería"
    );
}

/// `0<F<1`: la copia 1 reutiliza la fila de la copia 0 exactamente cuando su
/// Bernoulli lo dice, y se muestrea por su cuenta cuando no.
///
/// La batería exige además los dos casos que se confunden con facilidad: un
/// individuo homocigoto **por muestreo independiente**, que no es una
/// reutilización, y uno heterocigoto.
#[test]
fn an_intermediate_f_reuses_exactly_when_its_bernoulli_says_so() {
    let root = temporary("f-intermediate");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 16).inbreeding("0.5")],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();

    let mut reused = 0_u64;
    let mut independent = 0_u64;
    let mut homozygous_by_independent_sampling = 0_u64;
    let mut heterozygous = 0_u64;
    for unit in assignment.units() {
        for (individual_id, first, second) in diploid_pairs(unit) {
            assert_eq!(
                first,
                independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 0),
                "la copia 0 nunca reutiliza"
            );
            if reuses_copy_zero(root_seed, &scenario_id, unit, &individual_id, 0.5) {
                assert_eq!(
                    second, first,
                    "{individual_id}: su Bernoulli reutiliza la fila de la copia 0"
                );
                reused += 1;
            } else {
                let own =
                    independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 1);
                assert_eq!(
                    second, own,
                    "{individual_id}: sin reutilizar, la copia 1 usa su propio ámbito"
                );
                independent += 1;
                if own == first {
                    homozygous_by_independent_sampling += 1;
                }
            }
            if first != second {
                heterozygous += 1;
            }
        }
    }
    assert_eq!(
        reused + independent,
        32,
        "dos bloques × dieciséis individuos"
    );
    assert!(reused > 0, "la batería debe cubrir la reutilización");
    assert!(independent > 0, "y el muestreo independiente");
    assert!(
        homozygous_by_independent_sampling > 0,
        "un homocigoto por muestreo independiente no es una reutilización"
    );
    assert!(heterozygous > 0, "y debe haber heterocigotos");
}

/// `F=1` produce homocigosis completa: cada copia 1 reutiliza la fila de su
/// copia 0 y ninguna abre su propio muestreo.
#[test]
fn a_unit_f_produces_complete_homozygosity() {
    let root = temporary("f-one");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8).inbreeding("1.0")],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let scenario_id = ScenarioId::parse("scenario", "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", SEED).expect("seed");
    let weights = chain.config.profile().haplotype_weights().to_vec();

    let mut differing = 0_u64;
    for unit in assignment.units() {
        let mut ranks = BTreeSet::new();
        for (individual_id, first, second) in diploid_pairs(unit) {
            assert_eq!(
                first,
                independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 0),
                "la copia 0 se muestrea igual con cualquier F"
            );
            assert_eq!(second, first, "{individual_id} debe ser homocigoto con F=1");
            // La fila que la copia 1 habría sorteado por su cuenta: cuando
            // difiere, la reutilización es observable.
            if independent_rank(&weights, root_seed, &scenario_id, unit, &individual_id, 1) != first
            {
                differing += 1;
            }
            ranks.insert(first);
        }
        assert!(
            ranks.len() > 1,
            "el fixture debe realizar más de una categoría"
        );
        let realized = assignment
            .f_realized()
            .iter()
            .find(|row| row.block_id() == unit.block_id())
            .expect("fila de F realizado");
        assert_eq!(realized.status(), F5RealizedFStatus::Evaluated);
        assert_eq!(realized.target_f(), Some(1.0));
        assert_eq!(
            realized.value(),
            Some(1.0),
            "sin heterocigotos, F realizado vale exactamente 1"
        );
    }
    assert!(
        differing > 0,
        "al menos una copia 1 habría sorteado otra fila: la reutilización se observa"
    );
}

/// Las probabilidades genotípicas multialélicas de REQ-EVO-48 se cumplen en la
/// vía productiva completa: `P(ii) = q_i² + F·q_i·(1-q_i)` y
/// `P(ij) = 2·(1-F)·q_i·q_j` con `i<j`.
///
/// Las demás pruebas fijan el **mecanismo** —la copia 1 reutiliza si y solo si
/// su Bernoulli lo dice—; esta comprueba la **distribución** que ese mecanismo
/// induce, que es lo que REQ-VAL-23 pide verificar. Con `F=1/2`,
/// `q = (1/2, 1/4, 1/8, 1/8)` y 8 000 individuos diploides, la desviación
/// máxima observada es 0,0045; la cota de 0,02 deja margen frente al muestreo
/// sin dejar pasar una fórmula equivocada, cuyos errores son un orden mayor.
#[test]
fn the_multiallelic_f_formulas_hold_end_to_end() {
    let root = temporary("formulas");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 4000).inbreeding("0.5")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let weights = chain.config.profile().haplotype_weights().to_vec();
    let inbreeding = 0.5_f64;
    let mut observed = vec![vec![0_u64; weights.len()]; weights.len()];
    let mut total = 0_u64;
    for unit in assignment.units() {
        for (_, first, second) in diploid_pairs(unit) {
            let (low, high) = if first <= second {
                (first as usize, second as usize)
            } else {
                (second as usize, first as usize)
            };
            observed[low][high] += 1;
            total += 1;
        }
    }
    assert_eq!(total, 8000, "dos bloques × cuatro mil individuos");
    let mut categories = 0_u64;
    for i in 0..weights.len() {
        for j in i..weights.len() {
            let expected = if i == j {
                weights[i] * weights[i] + inbreeding * weights[i] * (1.0 - weights[i])
            } else {
                2.0 * (1.0 - inbreeding) * weights[i] * weights[j]
            };
            let seen = observed[i][j] as f64 / total as f64;
            assert!(
                (seen - expected).abs() <= 0.02,
                "P({i},{j}) observada {seen} frente a {expected} (N={total})"
            );
            categories += 1;
        }
    }
    assert_eq!(
        categories, 10,
        "las diez categorías i<=j de un perfil de cuatro"
    );
}

/// `F_realized` se contrasta con su racional exacto en cada unidad, y usar el
/// objetivo `q` en lugar de la frecuencia realizada daría otro valor.
#[test]
fn the_realized_f_matches_its_rational_oracle_in_both_modes() {
    // La tercera fila es una cuota controlada **no exacta**: con `C=6` los
    // ideales `6*q_h` valen 3, 1,5, 0,75 y 0,75, el reparto queda `[3,1,1,1]` y
    // la frecuencia realizada `(1/2, 1/6, 1/6, 1/6)` no es el objetivo. Sin
    // ella, la variante «frecuencia objetivo» no se distinguiría en controlado.
    for (label, mode, declared, individuals) in [
        ("probabilistic", "probabilistic", Some("0.5"), 8_u64),
        ("controlled-exact", "controlled", None, 8),
        ("controlled-inexact", "controlled", None, 3),
    ] {
        let root = temporary(&format!("realized-{label}"));
        let mut subpopulation = Subpopulation::new("P1", individuals).sampling(mode);
        if let Some(value) = declared {
            subpopulation = subpopulation.inbreeding(value);
        }
        let fixture = Fixture {
            subpopulations: vec![subpopulation],
            ..dyadic_fixture()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let weights = chain.config.profile().haplotype_weights().to_vec();

        assert_eq!(
            assignment.f_realized().len(),
            assignment.units().len(),
            "{label}: una fila de F realizado por unidad"
        );
        let mut compared = 0_u64;
        for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
            assert_eq!(realized.population_id(), unit.population_id());
            assert_eq!(realized.block_id(), unit.block_id());
            assert_eq!(
                realized.target_f(),
                Some(declared.map_or(0.0, |value| value.parse::<f64>().expect("F decimal"))),
                "{label}: el objetivo publicado es el F efectivo"
            );
            let expected = realized_f_oracle(unit, weights.len()).expect("H_expected positivo");
            assert_eq!(realized.status(), F5RealizedFStatus::Evaluated);
            let value = realized.value().expect("F realizado evaluado");
            assert!(
                (value - expected).abs() <= 1e-12,
                "{label}: F realizado {value} frente al racional exacto {expected}"
            );
            // Con el objetivo `q` del perfil el cociente es otro: la
            // comparación anterior distingue las dos frecuencias. Solo en la
            // cuota controlada exacta no aplica, porque allí la frecuencia
            // realizada **es** el objetivo.
            if label != "controlled-exact" {
                let with_target = target_frequency_variant(unit, &weights);
                assert!(
                    (with_target - expected).abs() > 1e-9,
                    "{label}: la frecuencia realizada debe diferir del objetivo q"
                );
            }
            compared += 1;
        }
        assert_eq!(compared, 2, "{label}: dos bloques comparados");
    }
}

/// Una subpoblación con reparto **mixto** publica su `F_realized` por bloque,
/// con el modo que cada bloque resolvió.
///
/// La sustitución deja el bloque 0 en controlado y el 1 en el predeterminado
/// probabilístico: la fila de F realizado es por bloque, no por subpoblación, y
/// el modo no la cambia.
#[test]
fn a_mixed_repartition_publishes_a_realized_f_per_block() {
    let root = temporary("mixed-realized");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8)],
        overrides: vec![Substitution {
            population: "P1",
            block_index: 0,
            mode: "controlled",
            tolerance: None,
        }],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let weights = chain.config.profile().haplotype_weights().to_vec();
    let mut modes = Vec::new();
    for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
        assert_eq!(realized.block_id(), unit.block_id());
        assert_eq!(realized.target_f(), Some(0.0));
        assert_eq!(realized.status(), F5RealizedFStatus::Evaluated);
        let expected = realized_f_oracle(unit, weights.len()).expect("H_expected positivo");
        let value = realized.value().expect("F realizado evaluado");
        assert!(
            (value - expected).abs() <= 1e-12,
            "bloque {}: F realizado {value} frente al racional exacto {expected}",
            unit.block_index()
        );
        modes.push((unit.block_index(), unit.mode()));
    }
    modes.sort_unstable();
    assert_eq!(
        modes,
        vec![
            (0, F5SamplingMode::Controlled),
            (1, F5SamplingMode::Probabilistic)
        ],
        "el fixture debe mezclar los dos modos en la misma subpoblación"
    );
}

/// Una sola categoría realizada deja `H_expected = 0` y el estado publicado es
/// `not_applicable` con valor nulo, nunca cero.
#[test]
fn a_null_expected_heterozygosity_is_published_as_not_applicable() {
    let root = temporary("h-expected-zero");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1).inbreeding("1.0")],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
        let pairs = diploid_pairs(unit);
        assert_eq!(pairs.len(), 1, "un solo individuo");
        assert_eq!(pairs[0].1, pairs[0].2, "homocigoto por F=1");
        assert_eq!(realized.status(), F5RealizedFStatus::NotApplicable);
        assert_eq!(realized.value(), None, "no aplicable no es cero");
        assert_eq!(realized.target_f(), Some(1.0), "el objetivo sí se conserva");
        assert_eq!(realized_f_oracle(unit, 2), None, "el oráculo coincide");
    }
}

/// Declarar F con `P=1` o `P>2` se rechaza al admitir el escenario, citando
/// ámbito y ploidía, y sin llegar a ninguna asignación.
#[test]
fn declaring_f_outside_a_diploid_is_rejected_before_founders() {
    for ploidy in [1_u64, 3] {
        for value in ["0.0", "0.25", "1.0"] {
            let root = temporary(&format!("f-ploidy-{ploidy}-{value}"));
            let fixture = Fixture {
                ploidy,
                subpopulations: vec![Subpopulation::new("P1", 2).inbreeding(value)],
                ..Fixture::default()
            };
            let general = fixture.write(&root);
            let failure =
                load_and_validate_f5_config(&general).expect_err("F fuera de P=2 debe rechazarse");
            assert_eq!(failure.code(), "GPS-CONF-101");
            assert_eq!(
                failure.context().get("fields").map(String::as_str),
                Some("inbreeding_f,ploidy"),
                "el rechazo debe citar ámbito y ploidía"
            );
            assert_eq!(
                failure.context().get("detail").map(String::as_str),
                Some("F solo se admite con ploidía 2")
            );
            let _ = fs::remove_dir_all(&root);
        }
    }
}

/// Sin F declarado, una ploidía distinta de 2 sí se ejecuta y publica su F
/// realizado como no aplicable, con objetivo y valor nulos.
///
/// Las dos primeras combinaciones dejan `C` **impar**, que en un diploide no
/// puede ocurrir: es donde `q_realized_h = n_h/C` deja de tener denominador
/// par y la fórmula tampoco debe evaluarse.
#[test]
fn a_non_diploid_publishes_a_not_applicable_realized_f() {
    for (ploidy, individuals, copies) in [(1_u64, 3_u64, 3_usize), (3, 1, 3), (3, 2, 6)] {
        let root = temporary(&format!("f-polyploid-{ploidy}-{individuals}"));
        let fixture = Fixture {
            ploidy,
            subpopulations: vec![Subpopulation::new("P1", individuals)],
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        assert_eq!(assignment.f_realized().len(), assignment.units().len());
        for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
            assert_eq!(
                unit.assignments().len(),
                copies,
                "P={ploidy} con N={individuals} asigna C={copies} copias"
            );
            assert_eq!(realized.status(), F5RealizedFStatus::NotApplicable);
            assert_eq!(
                realized.target_f(),
                None,
                "fuera de P=2 el objetivo es nulo"
            );
            assert_eq!(realized.value(), None);
        }
    }
}

/// El recuento de `F_realized` se suma al coste publicado de la unidad, no se
/// ejecuta fuera de la contabilidad.
///
/// La cota inferior solo usa las partes deterministas del modo controlado: el
/// recorrido de estados (`H`), el multiconjunto (`C`), la permutación (`C-1`),
/// la escritura de las asignaciones (`C`) y el recuento diploide (`C/2 + H`).
/// Sin este último, el total observado queda por debajo.
#[test]
fn the_realized_count_is_added_to_the_published_unit_cost() {
    for individuals in [8_u64, 32] {
        let root = temporary(&format!("realized-cost-{individuals}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", individuals).sampling("controlled")],
            ..dyadic_fixture()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let copies = individuals * 2;
        let haplotypes = chain.config.profile().haplotype_weights().len() as u64;
        let floor = 3 * copies + 2 * haplotypes - 1 + copies / 2;
        for unit in assignment.units() {
            assert!(
                unit.operations() >= floor,
                "C={copies}: {} visitas no incluyen el recuento de F realizado ({floor})",
                unit.operations()
            );
            assert!(
                unit.operations() <= 8 * (copies + 4),
                "C={copies}: el recuento no puede romper la cota lineal"
            );
        }
    }
}

/// F5 no reproduce: el módulo no contiene la fórmula de Ne por contribuciones
/// gaméticas, que es de F7.3, y por tanto no puede publicar un Ne evaluado.
///
/// El registro de «Ne reproductivo no aplicable» que pide DTC-089 es esta
/// ausencia: ningún esquema publicado tiene un campo de Ne y este hito no
/// inventa ninguno.
#[test]
fn f5_never_computes_a_reproductive_ne() {
    let source = include_str!("../src/f5_founders.rs");
    for forbidden in ["k_bar", "V_k", "k_j", "N_p"] {
        assert!(
            !source.contains(forbidden),
            "la fórmula de Ne es de F7.3 y no puede aparecer aquí: {forbidden}"
        );
    }
    // La búsqueda negativa solo vale si el módulo sí nombra lo suyo.
    assert!(source.contains("F_realized"));
}

/// F5.3.2 convierte en positiva la búsqueda negativa de F5.3.1: el módulo abre
/// ahora el ámbito `founder_f`, y lo abre **exactamente una vez**, dentro del
/// único camino que lo consume.
///
/// Las pruebas que fijan cuándo se abre y cuándo no son las unitarias de
/// `draw_copy_rank`, que observan el orden de apertura; esta solo impide que
/// vuelva a haber un segundo consumidor sin pasar por ellas.
#[test]
fn the_founder_f_scope_is_opened_exactly_once_in_the_module() {
    let source = include_str!("../src/f5_founders.rs");
    assert_eq!(
        source.matches("F5RandomDecision::founder_f(").count(),
        1,
        "un solo consumidor de founder_f, el de la copia 1 diploide"
    );
    assert!(
        source.contains("F5RandomDecision::founder_sample"),
        "la búsqueda solo vale si el módulo sí abre los ámbitos suyos"
    );
}

// ---------------------------------------------------------------------------
// F5.3.3. Identidades fundadoras, pertenencia y reconstructor independiente
// ---------------------------------------------------------------------------

/// Contexto BLAKE3 de DTC-039, copiado del encabezado de `id-catalog-v4`.
const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
/// Versión de derivación que el catálogo fija para las cinco clases F5.
const ID_DERIVATION_VERSION: u32 = 1;

/// Un componente DTC-037: longitud del nombre, nombre, longitud del valor y
/// valor, todos UTF-8 y con la longitud en `u32` big-endian.
fn push_component(material: &mut Vec<u8>, name: &str, value: &str) {
    material.extend_from_slice(&u32::try_from(name.len()).expect("nombre").to_be_bytes());
    material.extend_from_slice(name.as_bytes());
    material.extend_from_slice(&u32::try_from(value.len()).expect("valor").to_be_bytes());
    material.extend_from_slice(value.as_bytes());
}

/// Codificación canónica DTC-037 del material de un identificador automático.
fn id_material(root: &[u8; 32], class: &str, components: &[(&str, String)]) -> Vec<u8> {
    let mut material = Vec::new();
    material.extend_from_slice(root);
    material.extend_from_slice(&ID_DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(
        &u32::try_from(components.len() + 1)
            .expect("recuento")
            .to_be_bytes(),
    );
    push_component(&mut material, "class", class);
    for (name, value) in components {
        push_component(&mut material, name, value);
    }
    material
}

/// Base32 sin relleno escrita bit a bit, sin reutilizar la del productor.
fn base32(bytes: &[u8]) -> String {
    const ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
    let bits = bytes.len() * 8;
    (0..bits.div_ceil(5))
        .map(|group| {
            let mut index = 0_usize;
            for offset in 0..5 {
                let bit = group * 5 + offset;
                let set = bit < bits && (bytes[bit / 8] >> (7 - bit % 8)) & 1 == 1;
                index = (index << 1) | usize::from(set);
            }
            char::from(ALPHABET[index])
        })
        .collect()
}

/// Rederiva un identificador automático desde sus componentes publicados.
fn derived_identifier(
    root: &[u8; 32],
    prefix: &str,
    class: &str,
    components: &[(&str, String)],
) -> String {
    let material = id_material(root, class, components);
    let digest = blake3::derive_key(ID_CONTEXT, &material);
    format!("{prefix}{}", base32(&digest[..16]))
}

fn decode_hex(text: &str) -> Vec<u8> {
    assert!(text.len().is_multiple_of(2), "hexadecimal de longitud par");
    text.as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            u8::from_str_radix(std::str::from_utf8(pair).expect("hexadecimal"), 16)
                .expect("hexadecimal")
        })
        .collect()
}

fn encode_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

fn root_seed_bytes(seed: &str) -> [u8; 32] {
    decode_hex(seed).try_into().expect("semilla de 256 bits")
}

/// Componentes que el catálogo da a `FCP_`, en su orden exacto.
fn founder_copy_components(
    scenario_id: &str,
    reference_id: &str,
    copy: &F5FounderCopy,
) -> Vec<(&'static str, String)> {
    vec![
        ("execution_id", EXECUTION.to_owned()),
        ("scenario_id", scenario_id.to_owned()),
        ("reference_id", reference_id.to_owned()),
        ("population_id", copy.population_id().as_str().to_owned()),
        ("individual_id", copy.individual_id().as_str().to_owned()),
        ("sequence_id", copy.sequence_id().as_str().to_owned()),
        ("copy_index_u32", copy.copy_index().to_string()),
    ]
}

/// Rederiva cada `FCP_`, `HAP_` y `FHA_` leyendo solo lo que el repositorio
/// publica más el escenario, la referencia, la ejecución y la semilla
/// efectiva, y sin llamar a ninguna función del productor.
fn assert_every_identity_is_rederivable(chain: &Chain, identities: &Materialized) {
    let root = root_seed_bytes(chain.config.seed());
    let scenario_id = chain.config.scenario_id();
    let reference_id = chain.config.reference_id();

    for copy in identities.copies() {
        let expected = derived_identifier(
            &root,
            "FCP_",
            "f5_founder_copy",
            &founder_copy_components(scenario_id, reference_id, copy),
        );
        assert_eq!(
            copy.founder_copy_id().as_str(),
            expected,
            "la copia {}/{}/{}/{} no rederiva su FCP_",
            copy.population_id(),
            copy.individual_id(),
            copy.sequence_id(),
            copy.copy_index()
        );
    }

    for haplotype in identities.haplotypes() {
        let expected = derived_identifier(
            &root,
            "HAP_",
            "f5_haplotype",
            &[
                ("execution_id", EXECUTION.to_owned()),
                ("scenario_id", scenario_id.to_owned()),
                ("reference_id", reference_id.to_owned()),
                (
                    "population_id",
                    haplotype.population_id().as_str().to_owned(),
                ),
                ("block_index_u64", haplotype.block_index().to_string()),
                ("haplotype_rank_u64", haplotype.haplotype_rank().to_string()),
                ("alleles", haplotype.alleles().to_owned()),
            ],
        );
        assert_eq!(
            haplotype.haplotype_id().as_str(),
            expected,
            "la fila {} del bloque {} no rederiva su HAP_",
            haplotype.haplotype_rank(),
            haplotype.block_id()
        );
    }

    for assignment in identities.assignments() {
        let expected = derived_identifier(
            &root,
            "FHA_",
            "f5_block_haplotype_assignment",
            &[
                ("execution_id", EXECUTION.to_owned()),
                ("scenario_id", scenario_id.to_owned()),
                ("reference_id", reference_id.to_owned()),
                (
                    "population_id",
                    assignment.population_id().as_str().to_owned(),
                ),
                ("block_index_u64", assignment.block_index().to_string()),
                (
                    "founder_copy_id",
                    assignment.founder_copy_id().as_str().to_owned(),
                ),
                (
                    "haplotype_id",
                    assignment.haplotype_id().as_str().to_owned(),
                ),
            ],
        );
        assert_eq!(
            assignment.block_assignment_id().as_str(),
            expected,
            "la asignación de {} en {} no rederiva su FHA_",
            assignment.founder_copy_id(),
            assignment.block_id()
        );

        // La asignación enlaza dos identidades ya publicadas y no inventa
        // ninguna: su copia y su fila existen en las otras dos tablas.
        let copy = identities
            .copy(assignment.founder_copy_id())
            .expect("la asignación resuelve su pertenencia");
        assert_eq!(copy.population_id(), assignment.population_id());
        assert_eq!(copy.individual_id(), assignment.individual_id());
        assert_eq!(copy.sequence_id(), assignment.sequence_id());
        assert_eq!(copy.copy_index(), assignment.copy_index());
        assert!(
            identities.haplotypes().iter().any(|haplotype| {
                haplotype.haplotype_id() == assignment.haplotype_id()
                    && haplotype.block_id() == assignment.block_id()
                    && haplotype.population_id() == assignment.population_id()
                    && haplotype.haplotype_rank() == assignment.haplotype_rank()
            }),
            "la asignación nombra una fila haplotípica que no se publicó"
        );
    }
}

/// El oráculo de identidades debe reproducir los vectores congelados del
/// catálogo antes de poder juzgar al productor.
///
/// Sin esta prueba, un oráculo equivocado y un productor equivocado de la
/// misma forma se darían la razón mutuamente.
#[test]
fn the_independent_id_oracle_reproduces_the_frozen_catalog_vectors() {
    const CATALOG: &str = include_str!("../catalogs/id-catalog-v4.tsv");
    let mut checked = 0_usize;
    for line in CATALOG.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('|').collect::<Vec<_>>();
        assert_eq!(fields.len(), 9, "fila del catálogo: {line}");
        let (class, prefix) = (fields[0], fields[1]);
        if !class.starts_with("f5_") {
            continue;
        }
        let mut components = fields[4]
            .split(';')
            .map(|component| component.split_once('=').expect("componente nombre=valor"));
        let declared_class = components.next().expect("clase");
        assert_eq!(declared_class, ("class", class));
        let components = components
            .map(|(name, value)| (name, value.to_owned()))
            .collect::<Vec<_>>();
        let root: [u8; 32] = decode_hex(fields[5]).try_into().expect("raíz de 256 bits");
        let borrowed = components
            .iter()
            .map(|(name, value)| (*name, value.clone()))
            .collect::<Vec<_>>();
        let material = id_material(&root, class, &borrowed);
        assert_eq!(encode_hex(&material), fields[6], "material de {class}");
        let digest = blake3::derive_key(ID_CONTEXT, &material);
        assert_eq!(encode_hex(&digest), fields[7], "digest de {class}");
        assert_eq!(
            derived_identifier(&root, prefix, class, &borrowed),
            fields[8],
            "identificador de {class}"
        );
        checked += 1;
    }
    assert_eq!(checked, 5, "las cinco clases F5 del catálogo");
}

/// Reconstructor completo: dos subpoblaciones —una con IDs aportados y otra
/// con claves derivadas—, dos secuencias, dos bloques y `P=2`.
#[test]
fn every_published_founder_identity_is_rederivable_by_an_independent_reader() {
    let root = temporary("reconstructor");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 2),
            Subpopulation::new("P2", 2).keys(),
        ],
        sequences: 2,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    assert_eq!(
        identities.copies().len(),
        2 * 2 * 2 * 2,
        "N*P por secuencia"
    );
    assert_eq!(identities.assignments().len(), 2 * 2 * 2 * 2);
    assert_eq!(
        identities.haplotypes().len(),
        2 * 2 * chain.config.profile().haplotype_weights().len(),
        "H filas por bloque×subpoblación"
    );
    assert_every_identity_is_rederivable(&chain, &identities);

    // Idempotencia: materializar dos veces las mismas entradas da lo mismo. Se
    // compara lo conservado —copias y unidades— y lo recogido de los dos
    // recorridos, que es todo lo que el repositorio produce.
    let repeated = materialize(&chain);
    assert_eq!(identities.copies(), repeated.copies());
    assert_eq!(identities.units(), repeated.units());
    assert_eq!(identities.haplotypes(), repeated.haplotypes());
    assert_eq!(identities.assignments(), repeated.assignments());
    for copy in identities.copies() {
        assert_eq!(
            identities.copy(copy.founder_copy_id()),
            repeated.copy(copy.founder_copy_id()),
            "la pertenencia también es idéntica"
        );
    }
}

/// Los IDs aportados por la usuaria y los derivados de una clave conviven en
/// la misma ejecución y producen copias con el mismo contrato.
#[test]
fn provided_and_derived_individual_ids_produce_the_same_copy_contract() {
    let root = temporary("provided-and-derived");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 2).founders(&["given-a", "given-b"]),
            Subpopulation::new("P2", 2).keys(),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let provided = identities
        .copies()
        .iter()
        .filter(|copy| copy.population_id().as_str() == "P1")
        .map(|copy| copy.individual_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        provided,
        BTreeSet::from(["given-a", "given-b"]),
        "un founder_id aportado se usa literalmente"
    );
    let derived = identities
        .copies()
        .iter()
        .filter(|copy| copy.population_id().as_str() == "P2")
        .map(|copy| copy.individual_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(derived.len(), 2);
    assert!(
        derived.iter().all(|id| id.starts_with("IND_")),
        "una founder_key deriva su IndividualId: {derived:?}"
    );
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// Una copia conserva **una** identidad a través de todos los bloques de su
/// secuencia, y cada bloque le da una asignación distinta.
#[test]
fn one_copy_keeps_a_single_identity_across_every_block_of_its_sequence() {
    let root = temporary("copy-across-blocks");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2)],
        block_count: 2,
        sequences: 1,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let blocks = identities
        .assignments()
        .iter()
        .map(|assignment| assignment.block_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(blocks.len(), 2, "el fixture debe tener dos bloques");

    let mut by_copy: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    let mut assignment_ids: BTreeSet<&str> = BTreeSet::new();
    for assignment in identities.assignments() {
        by_copy
            .entry(assignment.founder_copy_id().as_str())
            .or_default()
            .insert(assignment.block_id().as_str());
        assert!(
            assignment_ids.insert(assignment.block_assignment_id().as_str()),
            "dos asignaciones no pueden compartir FHA_"
        );
        assert_ne!(
            assignment.block_assignment_id().as_str(),
            assignment.founder_copy_id().as_str(),
            "FHA_ no puede sustituir a FCP_"
        );
    }
    assert_eq!(by_copy.len(), 4, "N*P copias de una sola secuencia");
    for (copy, seen) in &by_copy {
        assert_eq!(
            seen.len(),
            2,
            "la copia {copy} debe aparecer en los dos bloques con el mismo FCP_"
        );
    }
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// Un individuo tiene una copia propia **por secuencia**: dos secuencias
/// producen dos `FCP_` distintos para el mismo `copy_index`, y cada bloque usa
/// el de su propia secuencia.
#[test]
fn an_individual_has_one_copy_identity_per_sequence() {
    let root = temporary("copy-per-sequence");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1)],
        block_count: 2,
        sequences: 2,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let sequences = identities
        .copies()
        .iter()
        .map(|copy| copy.sequence_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        sequences,
        BTreeSet::from(["chr1", "chr2"]),
        "las dos secuencias pobladas del layout"
    );
    assert_eq!(
        identities.copies().len(),
        4,
        "1 individuo × 2 copias × 2 secuencias"
    );

    for copy_index in 0..2_u32 {
        let ids = identities
            .copies()
            .iter()
            .filter(|copy| copy.copy_index() == copy_index)
            .map(|copy| copy.founder_copy_id().as_str())
            .collect::<BTreeSet<_>>();
        assert_eq!(
            ids.len(),
            2,
            "la copia {copy_index} tiene un FCP_ distinto por secuencia"
        );
    }

    // Cada bloque vive en una secuencia y sus asignaciones usan exactamente
    // las copias de esa secuencia.
    let mut sequence_of_block: BTreeMap<&str, &str> = BTreeMap::new();
    for assignment in identities.assignments() {
        let sequence = sequence_of_block
            .entry(assignment.block_id().as_str())
            .or_insert_with(|| assignment.sequence_id().as_str());
        assert_eq!(
            *sequence,
            assignment.sequence_id().as_str(),
            "un bloque no puede recibir copias de dos secuencias"
        );
        let copy = identities
            .copy(assignment.founder_copy_id())
            .expect("pertenencia");
        assert_eq!(copy.sequence_id(), assignment.sequence_id());
    }
    assert_eq!(
        sequence_of_block.values().collect::<BTreeSet<_>>().len(),
        2,
        "los dos bloques del fixture están en secuencias distintas"
    );
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// La ploidía fija cuántas copias existen y qué índices publican: `P=1` no
/// materializa ninguna copia 1 y `P=2` materializa exactamente dos por
/// secuencia.
#[test]
fn the_copy_identities_follow_the_effective_ploidy() {
    for ploidy in [1_u64, 2, 3] {
        let root = temporary(&format!("ploidy-{ploidy}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", 3)],
            ploidy,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let identities = materialize(&chain);

        let expected = 3 * ploidy as usize;
        assert_eq!(
            identities.copies().len(),
            expected,
            "N*P copias con P={ploidy}"
        );
        let indices = identities
            .copies()
            .iter()
            .map(F5FounderCopy::copy_index)
            .collect::<BTreeSet<_>>();
        assert_eq!(
            indices,
            (0..ploidy as u32).collect::<BTreeSet<_>>(),
            "los índices de copia con P={ploidy}"
        );
        assert_eq!(
            identities.assignments().len(),
            expected * 2,
            "una asignación por copia y bloque"
        );
        assert_every_identity_is_rederivable(&chain, &identities);
    }
}

/// La pertenencia es total, única e inyectiva: cada copia publicada tiene
/// exactamente un `FCP_` y cada `FCP_` resuelve una sola subpoblación.
#[test]
fn every_copy_has_exactly_one_membership_and_every_membership_one_copy() {
    let root = temporary("membership");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2), Subpopulation::new("P2", 3)],
        sequences: 2,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let mut identifiers: BTreeSet<&str> = BTreeSet::new();
    let mut tuples: BTreeSet<(&str, &str, &str, u32)> = BTreeSet::new();
    let mut population_of: BTreeMap<&str, &str> = BTreeMap::new();
    for copy in identities.copies() {
        assert!(
            identifiers.insert(copy.founder_copy_id().as_str()),
            "un FCP_ repetido es una pertenencia duplicada"
        );
        assert!(
            tuples.insert((
                copy.population_id().as_str(),
                copy.individual_id().as_str(),
                copy.sequence_id().as_str(),
                copy.copy_index(),
            )),
            "una tupla de identidad repetida es una pertenencia duplicada"
        );
        assert_eq!(
            identities
                .copy(copy.founder_copy_id())
                .expect("pertenencia")
                .population_id(),
            copy.population_id()
        );
        // Un individuo pertenece a una sola subpoblación en toda la ejecución.
        let population = population_of
            .entry(copy.individual_id().as_str())
            .or_insert_with(|| copy.population_id().as_str());
        assert_eq!(*population, copy.population_id().as_str());
    }
    assert_eq!(
        identifiers.len(),
        5 * 2 * 2,
        "(2+3) individuos × 2 secuencias × P=2"
    );
    assert_eq!(population_of.len(), 5);

    // Ninguna asignación queda sin pertenencia y ninguna copia queda sin
    // asignación en los bloques de su secuencia.
    let mut assigned: BTreeSet<&str> = BTreeSet::new();
    for assignment in identities.assignments() {
        assert!(
            identities.copy(assignment.founder_copy_id()).is_some(),
            "pertenencia ausente para {}",
            assignment.founder_copy_id()
        );
        assigned.insert(assignment.founder_copy_id().as_str());
    }
    assert_eq!(assigned, identifiers, "toda copia recibe asignación");

    // Un `FCP_` que nadie materializó no tiene pertenencia: la consulta no
    // puede devolver una copia cualquiera.
    let unknown = FounderCopyId::parse("FCP_AAAAAAAAAAAAAAAAAAAAAAAAAA", "test").expect("id");
    assert!(identities.copy(&unknown).is_none());
}

/// DTC-089: `IndividualId` identifica al fundador una sola vez en toda la
/// ejecución. El cruce entre dos subpoblaciones ya se rechaza en la admisión,
/// antes de que ninguna identidad se materialice.
///
/// Esta prueba fija **dónde** vive el invariante: `f5_founders` conserva una
/// guarda propia, con pruebas unitarias, porque una pertenencia duplicada no
/// produciría colisión de `FCP_`, pero ninguna configuración admitida puede
/// alcanzarla.
#[test]
fn an_individual_declared_in_two_subpopulations_is_rejected_before_any_identity() {
    let root = temporary("shared-individual");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 2).founders(&["shared", "only-p1"]),
            Subpopulation::new("P2", 2).founders(&["shared", "only-p2"]),
        ],
        ..Fixture::default()
    };
    let general = fixture.write(&root);
    let failure = load_and_validate_f5_config(&general).expect_err("fundador compartido");
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("subpopulations.founder_ids")
    );
    let _ = fs::remove_dir_all(&root);
}

/// La misma `founder_key` en dos subpoblaciones **no** es el mismo fundador:
/// la clave se deriva en el ámbito de su población, de modo que produce dos
/// `IndividualId` distintos, dos pertenencias distintas y dos `FCP_` distintos.
#[test]
fn the_same_founder_key_in_two_subpopulations_yields_two_founders() {
    let root = temporary("shared-key");
    let fixture = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 2).keys().founders(&["k0", "k1"]),
            Subpopulation::new("P2", 2).keys().founders(&["k0", "k1"]),
        ],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let individuals = identities
        .copies()
        .iter()
        .map(|copy| copy.individual_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        individuals.len(),
        4,
        "dos claves por subpoblación, no dos en total"
    );
    let identifiers = identities
        .copies()
        .iter()
        .map(|copy| copy.founder_copy_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(identifiers.len(), identities.copies().len());
    for copy in identities.copies() {
        assert_eq!(
            identities
                .copy(copy.founder_copy_id())
                .expect("pertenencia")
                .population_id(),
            copy.population_id()
        );
    }
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// El layout aporta la secuencia de cada copia: si no lo produjo esta
/// ejecución, sus identificadores de bloque no se rederivan y la
/// materialización se rechaza antes de derivar ninguna identidad.
#[test]
fn a_layout_from_another_execution_is_rejected() {
    let root = temporary("foreign-layout");
    let chain = build_chain(&root, &Fixture::default());
    let assignment = assign(&chain);
    let foreign = ExecutionId::parse("RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &foreign,
        &chain.layout,
        &chain.search,
        &assignment,
    )
    .expect_err("el layout no pertenece a esa ejecución");
    assert_eq!(error.code(), "GPS-CONF-101");
}

/// La admisión sellada de otra configuración no autoriza materializar estas
/// identidades.
#[test]
fn a_sealed_admission_from_another_configuration_is_rejected() {
    let own_root = temporary("own-admission");
    let other_root = temporary("other-admission");
    let own = build_chain(&own_root, &Fixture::default());
    let other = build_chain(
        &other_root,
        &Fixture {
            subpopulations: vec![Subpopulation::new("Q1", 2)],
            ..Fixture::default()
        },
    );
    let assignment = assign(&own);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &own.config,
        &other.sealed,
        &execution,
        &own.layout,
        &own.search,
        &assignment,
    )
    .expect_err("la admisión sellada es de otra configuración");
    assert_eq!(error.code(), "GPS-CONF-101");
}

/// Cada tabla se publica en la clave que DTC-089 le da.
#[test]
fn every_published_table_is_ordered_by_its_own_key() {
    let root = temporary("order");
    // Una sola secuencia con dos bloques: así cada copia tiene **dos** filas y
    // `block_id` es la componente que desempata. Con un bloque por secuencia la
    // tupla anterior ya sería única y la quinta clave nunca se ejercitaría.
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P2", 2), Subpopulation::new("P1", 2)],
        sequences: 1,
        block_count: 2,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let copies = identities
        .copies()
        .iter()
        .map(|copy| {
            (
                copy.population_id().as_str(),
                copy.individual_id().as_str(),
                copy.sequence_id().as_str(),
                copy.copy_index(),
            )
        })
        .collect::<Vec<_>>();
    let mut sorted = copies.clone();
    sorted.sort_unstable();
    assert_eq!(
        copies, sorted,
        "copias por población, individuo, secuencia y copia"
    );

    let haplotypes = identities
        .haplotypes()
        .iter()
        .map(|haplotype| {
            (
                haplotype.population_id().as_str(),
                haplotype.block_id().as_str(),
                haplotype.haplotype_rank(),
            )
        })
        .collect::<Vec<_>>();
    let mut sorted = haplotypes.clone();
    sorted.sort_unstable();
    assert_eq!(haplotypes, sorted, "filas por población, bloque y rango");

    let assignments = identities
        .assignments()
        .iter()
        .map(|assignment| {
            (
                assignment.population_id().as_str(),
                assignment.individual_id().as_str(),
                assignment.sequence_id().as_str(),
                assignment.copy_index(),
                assignment.block_id().as_str(),
            )
        })
        .collect::<Vec<_>>();
    let mut sorted = assignments.clone();
    sorted.sort_unstable();
    assert_eq!(
        assignments, sorted,
        "asignaciones por población, individuo, secuencia, copia y bloque"
    );

    // La quinta componente solo discrimina si ordenar por `block_id` y por
    // `block_index` dan resultados distintos.
    let mut by_index = identities
        .assignments()
        .iter()
        .map(|assignment| (assignment.block_index(), assignment.block_id().as_str()))
        .collect::<Vec<_>>();
    by_index.sort_unstable();
    by_index.dedup();
    assert_eq!(by_index.len(), 2, "el fixture debe tener dos bloques");
    assert_ne!(
        by_index[0].1 < by_index[1].1,
        by_index[0].0 < by_index[1].0,
        "los dos órdenes deben divergir para que la prueba distinga"
    );
    let blocks_per_copy = identities.assignments().len() / identities.copies().len();
    assert_eq!(
        blocks_per_copy, 2,
        "cada copia debe recibir los dos bloques"
    );
}

/// Las filas haplotípicas publican **todos** los rangos de la matriz ganadora
/// con sus alelos literales, que es lo que DTC-089 hace componente del `HAP_`.
#[test]
fn the_haplotype_rows_publish_every_rank_with_the_alleles_of_the_winning_matrix() {
    let root = temporary("haplotype-rows");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 2)],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    for haplotype in identities.haplotypes() {
        let block = chain
            .search
            .block(haplotype.population_id(), haplotype.block_index())
            .expect("bloque ganador");
        let matrix = block.matrix();
        let rank = usize::try_from(haplotype.haplotype_rank()).expect("rango");
        let expected = (0..matrix.site_count())
            .map(|site| {
                if matrix.allele(rank, site).expect("alelo") {
                    '1'
                } else {
                    '0'
                }
            })
            .collect::<String>();
        assert_eq!(haplotype.alleles(), expected);
        assert_eq!(haplotype.alleles().len(), matrix.site_count());
    }

    for unit in assign(&chain).units() {
        let published = identities
            .haplotypes()
            .iter()
            .filter(|haplotype| {
                haplotype.population_id() == unit.population_id()
                    && haplotype.block_index() == unit.block_index()
            })
            .map(PublishedHaplotype::haplotype_rank)
            .collect::<Vec<_>>();
        let haplotype_count = chain
            .search
            .block(unit.population_id(), unit.block_index())
            .expect("bloque")
            .matrix()
            .haplotype_count() as u64;
        assert_eq!(published, (0..haplotype_count).collect::<Vec<_>>());

        // La asignación enlaza exactamente el rango que F5.3.1 sorteó.
        for seat in unit.assignments() {
            let assignment = identities
                .assignments()
                .iter()
                .find(|assignment| {
                    assignment.population_id() == unit.population_id()
                        && assignment.block_index() == unit.block_index()
                        && assignment.individual_id() == seat.individual_id()
                        && assignment.copy_index() == seat.copy_index()
                })
                .expect("la asignación publicada de esa copia");
            assert_eq!(assignment.haplotype_rank(), seat.haplotype_rank());
            let row = identities
                .haplotypes()
                .iter()
                .find(|haplotype| {
                    haplotype.population_id() == unit.population_id()
                        && haplotype.block_index() == unit.block_index()
                        && haplotype.haplotype_rank() == seat.haplotype_rank()
                })
                .expect("la fila publicada de ese rango");
            assert_eq!(
                assignment.haplotype_id(),
                row.haplotype_id(),
                "el FHA_ debe enlazar el HAP_ del rango que se le asignó"
            );
        }
    }
}

/// F5 no rompe ni recompone haplotipos: el módulo no nombra un proceso
/// haplotípico, un cruce, un mosaico ni un mapa de recombinación, y la
/// identidad de la copia no se deriva de una fase.
#[test]
fn the_module_never_names_a_haplotype_process_a_cross_a_mosaic_or_a_map() {
    let source = include_str!("../src/f5_founders.rs");
    for forbidden in [
        "HaplotypeProcess",
        "haplotype_process",
        "PhaseId",
        "phase_id",
        "mosaic",
        "mosaico",
        "recombin",
        "CrossId",
        "cross_id",
        "genetic_map",
        "recombination_map",
        "centimorgan",
    ] {
        assert!(
            !source.contains(forbidden),
            "F5 no crea mosaicos ni proceso haplotípico: {forbidden}"
        );
    }
    // La búsqueda negativa solo vale si el módulo sí nombra lo suyo.
    assert!(source.contains("F5FounderCopyScope"));
    assert!(source.contains("F5BlockHaplotypeAssignmentScope"));
    assert!(source.contains("copy_index"));
}

/// Layout, búsqueda y asignación de la **misma** configuración bajo otra
/// identidad de ejecución.
///
/// Los ámbitos RNG de F5 no llevan `execution_id`, de modo que la geometría
/// sale idéntica y lo único que cambia son los identificadores derivados: es
/// justo la pieza ajena que la materialización debe reconocer.
fn chain_under(
    chain: &Chain,
    execution: &ExecutionId,
) -> (F5PhysicalLayout, F5HaplotypeSearch, F5FounderAssignment) {
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&chain.config, &chain.sealed, execution).expect("solicitud de layout"),
    )
    .expect("layout físico");
    let search = search_f5_haplotypes(&chain.config, &chain.sealed, execution, &layout)
        .expect("búsqueda haplotípica terminada");
    let assignment = assign_f5_founders(&chain.config, &chain.sealed, &search)
        .expect("asignación fundadora completa");
    (layout, search, assignment)
}

/// Un layout de **otra referencia** de la misma ejecución publica copias sobre
/// secuencias que la referencia propia no contiene.
///
/// Escenario, referencia, semilla, ejecución y `block_count` coinciden, así que
/// ni la admisión sellada ni la rederivación de `block_id` lo distinguen: lo
/// distingue que sus bloques no caben en ningún tramo elegible propio.
#[test]
fn a_layout_of_another_reference_in_the_same_execution_is_rejected() {
    let own_root = temporary("own-reference");
    let foreign_root = temporary("foreign-reference");
    let own = build_chain(
        &own_root,
        &Fixture {
            sequences: 1,
            block_count: 2,
            ..Fixture::default()
        },
    );
    let foreign = build_chain(
        &foreign_root,
        &Fixture {
            sequences: 2,
            block_count: 2,
            ..Fixture::default()
        },
    );
    // El layout ajeno coloca un bloque en `chr2`, que la referencia propia no
    // tiene: sin la comprobación de tramos se publicarían `FCP_` fantasma.
    assert!(
        foreign
            .layout
            .blocks()
            .iter()
            .any(|block| block.sequence_id().as_str() == "chr2")
    );
    let assignment = assign(&own);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &own.config,
        &own.sealed,
        &execution,
        &foreign.layout,
        &own.search,
        &assignment,
    )
    .expect_err("el layout es de otra referencia");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("layout.geometry")
    );
}

/// La dirección contraria también: un layout con **menos** secuencias sitúa sus
/// bloques en coordenadas que ningún tramo propio contiene.
#[test]
fn a_layout_with_coordinates_outside_the_own_tracts_is_rejected() {
    let own_root = temporary("own-tracts");
    let foreign_root = temporary("foreign-tracts");
    let own = build_chain(
        &own_root,
        &Fixture {
            sequences: 2,
            block_count: 2,
            ..Fixture::default()
        },
    );
    let foreign = build_chain(
        &foreign_root,
        &Fixture {
            sequences: 1,
            block_count: 2,
            ..Fixture::default()
        },
    );
    let assignment = assign(&own);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &own.config,
        &own.sealed,
        &execution,
        &foreign.layout,
        &own.search,
        &assignment,
    )
    .expect_err("el layout ocupa coordenadas ajenas");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("layout.geometry")
    );
}

/// El recuento de bloques se comprueba, y no lo cubre la contención: los tres
/// bloques del layout ajeno caben de sobra dentro del único tramo propio.
#[test]
fn a_layout_with_another_block_count_is_rejected() {
    let own_root = temporary("own-block-count");
    let foreign_root = temporary("foreign-block-count");
    let own = build_chain(
        &own_root,
        &Fixture {
            block_count: 2,
            slack: 2,
            ..Fixture::default()
        },
    );
    let foreign = build_chain(
        &foreign_root,
        &Fixture {
            block_count: 3,
            ..Fixture::default()
        },
    );
    assert_eq!(foreign.layout.blocks().len(), 3);
    for block in foreign.layout.blocks() {
        assert!(
            block.end_1based() <= 12,
            "el layout ajeno debe caber en el tramo propio para aislar el recuento"
        );
    }
    let assignment = assign(&own);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &own.config,
        &own.sealed,
        &execution,
        &foreign.layout,
        &own.search,
        &assignment,
    )
    .expect_err("el layout tiene otro número de bloques");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("preflight.block_count")
    );
}

/// Una asignación construida sobre la búsqueda de **otra** ejecución nombra
/// bloques con otros `BLK_`.
#[test]
fn an_assignment_from_another_execution_is_rejected() {
    let root = temporary("foreign-assignment");
    let chain = build_chain(&root, &Fixture::default());
    let own = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let other = ExecutionId::parse("RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    let (_, _, foreign_assignment) = chain_under(&chain, &other);
    let error = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &own,
        &chain.layout,
        &chain.search,
        &foreign_assignment,
    )
    .expect_err("la asignación es de otra ejecución");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("assignment.identity")
    );
}

/// `alleles` debe medir los `S` sitios de **su** bloque: una matriz ganadora de
/// otra geometría se rechaza antes de derivar ningún `HAP_`.
#[test]
fn a_winning_matrix_with_another_site_count_is_rejected() {
    let own_root = temporary("own-sites");
    let foreign_root = temporary("foreign-sites");
    let own = build_chain(&own_root, &Fixture::default());
    let foreign = build_chain(
        &foreign_root,
        &Fixture {
            gap: 2,
            bins: vec![(1, 2, 1.0, 0.05)],
            ..Fixture::default()
        },
    );
    let own_sites = own.layout.blocks()[0].sites().len();
    let foreign_sites = foreign.search.blocks()[0].matrix().site_count();
    assert_ne!(
        own_sites, foreign_sites,
        "el fixture ajeno debe tener otra cardinalidad de sitios"
    );
    let assignment = assign(&own);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let error = materialize_f5_founder_identities(
        &own.config,
        &own.sealed,
        &execution,
        &own.layout,
        &foreign.search,
        &assignment,
    )
    .expect_err("la matriz no tiene los sitios de su bloque");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("search.geometry")
    );
}

/// Un `founder_id` literal con el texto exacto de un `IND_` derivado de una
/// `founder_key` de otra subpoblación es el mismo fundador dos veces.
///
/// La admisión no puede verlo —trabaja sobre el documento declarado, donde el
/// identificador derivado todavía no existe—, así que lo rechaza el censo.
#[test]
fn a_literal_founder_id_equal_to_a_derived_one_is_rejected() {
    let derived_root = temporary("derived-identity");
    let colliding_root = temporary("colliding-identity");
    let keys = ["k0", "k1"];
    let derived_fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P2", 2).keys().founders(&keys)],
        ..Fixture::default()
    };
    let general = derived_fixture.write(&derived_root);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let rosters = resolve_f5_founder_roster(&config).expect("censo fundador");
    let derived = rosters[0].individuals()[0]
        .individual_id()
        .as_str()
        .to_owned();
    assert!(derived.starts_with("IND_"), "{derived} debe ser derivado");
    let _ = fs::remove_dir_all(&derived_root);

    let literal = [derived.as_str()];
    let colliding = Fixture {
        subpopulations: vec![
            Subpopulation::new("P1", 1).founders(&literal),
            Subpopulation::new("P2", 2).keys().founders(&keys),
        ],
        ..Fixture::default()
    };
    let general = colliding.write(&colliding_root);
    let config = load_and_validate_f5_config(&general)
        .expect("la admisión no puede distinguir el literal de un derivado");
    let error = resolve_f5_founder_roster(&config).expect_err("el mismo fundador dos veces");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("subpopulations.founder_ids,subpopulations.founder_keys")
    );
    assert!(
        format!("{error:?}").contains(&derived),
        "el diagnóstico nombra al fundador repetido"
    );
    let _ = fs::remove_dir_all(&colliding_root);
}

/// El orden de `sequence_id` es el **lexicográfico** de la columna, no el de
/// aparición en el FASTA.
///
/// Con diez secuencias los dos órdenes divergen: el FASTA da `chr1 … chr10` y
/// la columna da `chr1, chr10, chr2, …`. Es el fixture que distingue las dos
/// lecturas posibles de la clave que DTC-089 escribe como «secuencia».
#[test]
fn the_sequence_key_is_lexicographic_and_not_the_fasta_order() {
    let root = temporary("lexicographic-sequences");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1)],
        sequences: 10,
        block_count: 10,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let identities = materialize(&chain);

    let published = identities
        .copies()
        .iter()
        .map(|copy| copy.sequence_id().as_str())
        .collect::<Vec<_>>();
    let mut lexicographic = published.clone();
    lexicographic.sort_unstable();
    assert_eq!(
        published, lexicographic,
        "la columna ordena su propio valor"
    );

    let mut order = Vec::new();
    for sequence in &published {
        if order.last() != Some(sequence) {
            order.push(*sequence);
        }
    }
    let mut fasta = chain
        .layout
        .blocks()
        .iter()
        .map(|block| (block.sequence_ordinal(), block.sequence_id().as_str()))
        .collect::<Vec<_>>();
    fasta.sort_unstable();
    fasta.dedup();
    let fasta_order = fasta
        .iter()
        .map(|(_, sequence)| *sequence)
        .collect::<Vec<_>>();
    assert_eq!(order.len(), 10);
    assert_eq!(fasta_order.len(), 10);
    assert_eq!(order[1], "chr10", "lexicográfico: chr10 va detrás de chr1");
    assert_eq!(fasta_order[1], "chr2", "FASTA: chr2 va detrás de chr1");
    assert_ne!(
        order, fasta_order,
        "con chr1…chr10 los dos órdenes tienen que divergir"
    );
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// El modo controlado materializa sus enlaces igual que el probabilístico: la
/// permutación de F5.3.1 reparte los rangos y cada copia conserva su `FCP_`.
#[test]
fn the_controlled_mode_materializes_its_links() {
    let root = temporary("controlled-links");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 8).sampling("controlled")],
        ..dyadic_fixture()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let identities = materialize(&chain);

    assert!(
        assignment
            .units()
            .iter()
            .all(|unit| unit.mode() == F5SamplingMode::Controlled && unit.quota().is_some()),
        "el fixture debe resolverse en controlado"
    );
    assert_eq!(identities.copies().len(), 16, "N=8 con P=2");
    for unit in assignment.units() {
        for seat in unit.assignments() {
            let row = identities
                .assignments()
                .iter()
                .find(|row| {
                    row.population_id() == unit.population_id()
                        && row.block_index() == unit.block_index()
                        && row.individual_id() == seat.individual_id()
                        && row.copy_index() == seat.copy_index()
                })
                .expect("la asignación publicada de esa copia");
            assert_eq!(row.haplotype_rank(), seat.haplotype_rank());
            let expected = identities
                .haplotypes()
                .iter()
                .find(|haplotype| {
                    haplotype.population_id() == unit.population_id()
                        && haplotype.block_index() == unit.block_index()
                        && haplotype.haplotype_rank() == seat.haplotype_rank()
                })
                .expect("la fila publicada de ese rango");
            assert_eq!(row.haplotype_id(), expected.haplotype_id());
        }
    }
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// Con `F>0` la copia 1 puede reutilizar la fila de la copia 0, pero **no** su
/// identidad: cada copia conserva su `FCP_` y cada bloque su `FHA_`.
#[test]
fn a_reused_row_never_collapses_two_copies_into_one_identity() {
    let root = temporary("inbreeding-links");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 6).inbreeding("0.5")],
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);
    let identities = materialize(&chain);

    let reused = assignment
        .units()
        .iter()
        .flat_map(|unit| {
            unit.assignments()
                .iter()
                .filter(|seat| seat.copy_index() == 1)
                .filter(|seat| {
                    unit.assignments().iter().any(|other| {
                        other.copy_index() == 0
                            && other.individual_id() == seat.individual_id()
                            && other.haplotype_rank() == seat.haplotype_rank()
                    })
                })
        })
        .count();
    assert!(reused > 0, "con F=1/2 alguna copia debe compartir fila");

    let copy_ids = identities
        .copies()
        .iter()
        .map(|copy| copy.founder_copy_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(copy_ids.len(), identities.copies().len());
    let assignment_ids = identities
        .assignments()
        .iter()
        .map(|row| row.block_assignment_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(assignment_ids.len(), identities.assignments().len());
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// Una secuencia **retenida por la referencia pero sin bloques** no aporta
/// ninguna copia fundadora, y eso no es un error.
///
/// DTC-089 identifica la copia «por secuencia … a través de todos los bloques
/// de esa secuencia», y las tuplas FCP son las de las secuencias que el layout
/// puebla. La consecuencia está escrita en DTC-089 y es deliberada: el fundador
/// no tiene copia de la secuencia vacía en `truth.jsonl`, aunque la referencia
/// la retenga.
///
/// El caso exige holgura: con `slack: 1` cada secuencia recibe exactamente sus
/// bloques y la divergencia es estructuralmente imposible. Con `slack: 3` y
/// esta semilla el layout coloca los dos bloques en `chr1` y deja `chr2` vacía.
#[test]
fn a_retained_sequence_without_blocks_publishes_no_copy() {
    const SEED: &str = "0000000000000000000000000000000000000000000000000000000000000001";
    let root = temporary("unpopulated-sequence");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1)],
        sequences: 2,
        block_count: 2,
        slack: 3,
        seed: SEED,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    assert_eq!(
        chain.sealed.retained_sequence_ids(),
        2,
        "la referencia retiene chr1 y chr2"
    );
    let mut populated = chain
        .layout
        .blocks()
        .iter()
        .map(|block| block.sequence_id().as_str())
        .collect::<Vec<_>>();
    populated.sort_unstable();
    populated.dedup();
    assert_eq!(
        populated,
        vec!["chr1"],
        "el fixture debe dejar chr2 retenida y sin bloques"
    );

    let identities = materialize(&chain);
    let published = identities
        .copies()
        .iter()
        .map(|copy| copy.sequence_id().as_str())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        published,
        BTreeSet::from(["chr1"]),
        "las copias son las de las secuencias pobladas, no las retenidas"
    );
    assert_eq!(
        identities.copies().len(),
        2,
        "un individuo, dos copias, una sola secuencia poblada"
    );
    assert!(
        (identities.copies().len() as u64) < chain.sealed.retained_sequence_ids() * 2,
        "publicar por secuencia retenida daría el doble de copias"
    );
    // Toda copia publicada recibe asignación en los dos bloques de su
    // secuencia: no queda ninguna copia huérfana ni ninguna secuencia a medias.
    assert_eq!(identities.assignments().len(), 4);
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// El predeterminado en copias es `0.5` **exacto**, no `0.5/C` reescalado.
///
/// Desde la resolución de la cota controlada, ese `0,5` exacto **es** el valor
/// normativo que DTC-089 declara, no una conveniencia de implementación: la
/// cota se evalúa en copias —`|n_h - fl(C*q_h)| <= t*C`— y el predeterminado
/// aporta el `0,5` literal, porque `(0.5/C)*C` deja de valer `0,5`
/// para muchos `C`. Con `C=49` vale `0.49999999999999994`, de modo que el
/// reescalado rechaza con `GPS-DOMAIN-103` una cuota que la cota inclusiva del
/// contrato acepta: `q={1/2,1/2}` reparte `[25,24]` y su error es exactamente
/// medio individuo.
///
/// El caso exige `C` impar y mayor que 16; ninguna prueba anterior pasaba de
/// ahí, y por eso el reescalado sobrevivía a la batería entera.
#[test]
fn the_default_copy_tolerance_is_exactly_one_half_and_not_a_rescaled_one() {
    // La divergencia es del `f64`, no del fixture: se comprueba aquí para que
    // la prueba falle por su causa y no por una elección de `C` equivocada.
    let copies = 49.0_f64;
    assert_ne!((0.5 / copies) * copies, 0.5);
    assert!((0.5 / copies) * copies < 0.5);

    let root = temporary("exact-half-tolerance");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 49).sampling("controlled")],
        ploidy: 1,
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);

    for unit in assignment.units() {
        let quota = unit.quota().expect("controlado publica su cuota");
        assert_eq!(quota.copies(), 49);
        let mut counts = quota.counts().to_vec();
        counts.sort_unstable();
        assert_eq!(
            counts,
            vec![24, 25],
            "la cuota entera de q={{1/2,1/2}} con C=49 reparte 24 y 25"
        );
        // La cota se cumple justo en el borde: medio individuo de error.
        assert_eq!(quota.max_abs_error_effective(), 0.5 / 49.0);
        for count in quota.counts() {
            let error = (*count as f64 - 49.0 * 0.5).abs();
            assert_eq!(error, 0.5, "el error en copias es exactamente medio");
        }
    }
    let identities = materialize(&chain);
    assert_eq!(
        identities.copies().len(),
        49,
        "P=1 da una copia por individuo"
    );
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// El modo controlado no es solo diploide: la costura F5.3.1→F5.3.3 se fija
/// también en `P=1` y `P=3`.
///
/// Con el perfil diádico y `C` potencia de dos la cuota es exacta, `copy_index`
/// recorre `0..P-1`, cada copia recibe su `FCP_` y cada bloque su `FHA_`, y
/// `F_realized` es `not_applicable` fuera de `P=2`.
#[test]
fn the_controlled_mode_works_outside_diploids() {
    for (ploidy, individuals) in [(1_u64, 16_u64), (3, 16)] {
        let root = temporary(&format!("controlled-ploidy-{ploidy}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", individuals).sampling("controlled")],
            ploidy,
            ..dyadic_fixture()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let identities = materialize(&chain);
        let copies = individuals * ploidy;

        for (unit, realized) in assignment.units().iter().zip(assignment.f_realized()) {
            assert_eq!(unit.mode(), F5SamplingMode::Controlled);
            let quota = unit.quota().expect("controlado publica su cuota");
            assert_eq!(quota.copies(), copies);
            assert_eq!(quota.counts().iter().sum::<u64>(), copies);
            assert_largest_remainder(quota.counts(), &fixture.weights, copies, (1, 2));
            assert_eq!(
                realized.status(),
                F5RealizedFStatus::NotApplicable,
                "fuera de P=2 el F realizado no aplica"
            );
            assert_eq!(unit.assignments().len() as u64, copies);
        }

        let indices = identities
            .copies()
            .iter()
            .map(F5FounderCopy::copy_index)
            .collect::<BTreeSet<_>>();
        assert_eq!(indices, (0..ploidy as u32).collect::<BTreeSet<_>>());
        assert_eq!(identities.copies().len() as u64, copies);
        assert_eq!(
            identities.assignments().len() as u64,
            copies * fixture.block_count
        );
        assert_every_identity_is_rederivable(&chain, &identities);
    }
}

/// La cuota ideal normativa es `fl(C*q_h)` y la cota se evalúa **en copias**.
///
/// `q={0,7, 0,3}` con cinco fundadores haploides es el caso que lo revela. En
/// decimales exactos —lo que la usuaria escribe— la cuota cumple la cota **con
/// igualdad**: `|3/5-7/10|` y `|2/5-3/10|` valen ambos `1/10`, que es la
/// tolerancia predeterminada `0,5/C`. En binary64, en cambio, `0,3` es
/// `0,29999999999999998890`, el ideal de cada rango cae en un semientero
/// exacto —`fl(5*0,7)=3,5` y `fl(5*0,3)=1,5`— y la desigualdad literal
/// `|n_h/C-q_h|<=t` se excede en un ulp.
///
/// DTC-089 resuelve que la desigualdad normativa es la de copias. Esta prueba
/// la fija por los dos lados: comprueba que la cuota publicada la cumple con
/// igualdad, y que las dos alternativas descartadas —evaluar `|n/C-q|` en
/// binary64 y evaluarla exactamente en racionales sobre esos mismos pesos— la
/// rechazarían.
#[test]
fn a_half_integer_ideal_is_accepted_because_the_bound_is_evaluated_in_copies() {
    let weights = ["0.7", "0.3"];
    let root = temporary("half-integer-ideal");
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 5).sampling("controlled")],
        ploidy: 1,
        weights: weights.to_vec(),
        ..Fixture::default()
    };
    let chain = build_chain(&root, &fixture);
    let assignment = assign(&chain);

    let effective = effective_weights(&weights);
    for unit in assignment.units() {
        let quota = unit.quota().expect("controlado publica su cuota");
        assert_eq!(quota.copies(), 5);
        assert_eq!(quota.counts().iter().sum::<u64>(), 5);
        assert_eq!(
            quota.counts(),
            &[3, 2],
            "el desempate del grupo fronterizo publica esta cuota"
        );
        assert_largest_remainder(quota.counts(), &weights, 5, (1, 2));

        for (rank, count) in quota.counts().iter().enumerate() {
            let ideal = 5.0_f64 * effective[rank];
            // El ideal cae en un semientero exacto y la cota se cumple justo.
            assert_eq!(ideal.fract(), 0.5);
            assert_eq!((*count as f64 - ideal).abs(), 0.5);
        }

        // Alternativa descartada 1: la desigualdad literal en binary64.
        let literal_tolerance = 0.5 / 5.0;
        let literal_rejects =
            quota.counts().iter().enumerate().any(|(rank, count)| {
                (*count as f64 / 5.0 - effective[rank]).abs() > literal_tolerance
            });
        assert!(
            literal_rejects,
            "evaluar |n/C-q| en binary64 rechazaría esta cuota legal"
        );

        // Alternativa descartada 2: la misma desigualdad en racionales
        // exactos sobre los pesos binary64. `|n/C - q| <= 1/10` se comprueba
        // como `|n*den - C*num| * 10 <= C * den`.
        let exact_rejects = quota.counts().iter().enumerate().any(|(rank, count)| {
            let (numerator, denominator) = exact_rational(effective[rank]);
            let left = (u128::from(*count) * denominator) as i128 - (5 * numerator) as i128;
            left.unsigned_abs() * 10 > 5 * denominator
        });
        assert!(
            exact_rejects,
            "evaluar la desigualdad en racionales exactos sobre los pesos binary64 también la rechazaría"
        );
    }

    // Y en decimales exactos, que es lo que la usuaria declaró, la cuota es
    // legal con igualdad: `|3/5-7/10| = |2/5-3/10| = 1/10`.
    for (count, (numerator, denominator)) in
        [3_i128, 2].into_iter().zip([(7_i128, 10_i128), (3, 10)])
    {
        assert_eq!(
            (count * denominator - 5 * numerator).abs() * 10,
            5 * denominator
        );
    }

    let identities = materialize(&chain);
    assert_every_identity_is_rederivable(&chain, &identities);
}

/// Barrido de perfiles **no diádicos** contra el oráculo racional exacto.
///
/// Hasta ahora el oráculo abortaba con «el oráculo exige pesos diádicos», de
/// modo que la clase de entrada que revela la semántica de la cota no podía
/// llegar a la comprobación. Con el ideal normativo `fl(C*q_h)` convertido a su
/// racional exacto, el oráculo valida cualquier perfil cuyo peso efectivo sea
/// el declarado.
///
/// Quince combinaciones —treinta unidades bloque×subpoblación—: cinco perfiles
/// cuyos decimales no son diádicos, con `C ∈ {5, 7, 10}` en haploides. Entre
/// ellas están los ideales semienteros de `0,7/0,3` con `C=5`, los que caen un
/// ulp por debajo de un entero —`fl(7*0,7)=4,8999999999999995`— y los que dan
/// un rango con cero copias —`0,9/0,1` con `C=5`—.
#[test]
fn the_integer_quota_oracle_accepts_non_dyadic_profiles() {
    // Cinco perfiles de dos rangos cuyos decimales no son diádicos. Los de
    // tres rangos quedan fuera de este barrido porque la geometría mínima del
    // fixture —`S=3` con `H=3`— no le da a la búsqueda LD ninguna matriz que
    // cumpla a la vez bins y signos, no por nada de la cuota; el perfil
    // diádico de `H=4` ya cubre el caso multirrango.
    const PROFILES: [&[&str]; 5] = [
        &["0.7", "0.3"],
        &["0.6", "0.4"],
        &["0.9", "0.1"],
        &["0.8", "0.2"],
        &["0.55", "0.45"],
    ];
    let mut checked = 0_usize;
    for (profile_index, weights) in PROFILES.iter().enumerate() {
        for individuals in [5_u64, 7, 10] {
            let root = temporary(&format!("non-dyadic-{profile_index}-{individuals}"));
            // `H=2` solo realiza `r²` en 0 o 1, de modo que su objetivo
            // alcanzable es 1,0. El perfil de pesos es lo único que el barrido
            // varía respecto del fixture por defecto.
            let fixture = Fixture {
                subpopulations: vec![Subpopulation::new("P1", individuals).sampling("controlled")],
                ploidy: 1,
                weights: weights.to_vec(),
                max_haplotype_attempts: 512,
                max_swap_proposals: 2048,
                ..Fixture::default()
            };
            let chain = build_chain(&root, &fixture);
            let assignment = assign(&chain);
            for unit in assignment.units() {
                let quota = unit.quota().expect("controlado publica su cuota");
                assert_eq!(quota.copies(), individuals);
                assert_largest_remainder(quota.counts(), weights, individuals, (1, 2));
                checked += 1;
            }
            // Las identidades se materializan igual sobre una cuota no diádica.
            let identities = materialize(&chain);
            assert_eq!(identities.copies().len() as u64, individuals);
            assert_every_identity_is_rederivable(&chain, &identities);
        }
    }
    assert_eq!(
        checked,
        PROFILES.len() * 3 * 2,
        "cinco perfiles × tres tamaños × dos bloques"
    );
}

// ---------------------------------------------------------------------------
// Lo retenido frente a lo presupuestado
// ---------------------------------------------------------------------------

/// Bytes que el repositorio retiene, contados por el tamaño real de cada tipo.
///
/// No cuenta el desbordamiento de nodo de los `BTreeMap` ni el texto compartido
/// de los identificadores, que es una sola asignación por identidad: es una
/// **cota inferior** de lo retenido, que es el lado conservador para la
/// pregunta que interesa —cómo crece— y el honesto para la que compara contra
/// el presupuesto.
fn retained_bytes(identities: &F5FounderIdentities) -> usize {
    let copies = size_of_val(identities.copies());
    // Una entrada de pertenencia: el identificador y la posición.
    let membership = identities.copies().len() * (size_of::<FounderCopyId>() + size_of::<usize>());
    let units: usize = identities
        .units()
        .iter()
        .map(|unit| {
            size_of::<F5UnitHaplotypes>() + size_of_val(unit.haplotype_ids()) + size_of::<usize>()
        })
        .sum();
    copies + membership + units
}

/// Bytes de la única asignación del texto de una identidad derivada.
///
/// Cada `FCP_`, `HAP_` o `FHA_` es un `Arc<str>`: la cabecera de dos contadores
/// y los treinta caracteres del prefijo más el base32 de 128 bits.
const SHARED_IDENTIFIER_BYTES: usize = 2 * size_of::<usize>() + 30;

/// Bytes de una entrada del registro de colisiones.
///
/// La pareja `(clase, texto)` que vive en el nodo del árbol más la asignación
/// de su texto, sin contar el redondeo del asignador ni el relleno parcial de
/// los nodos.
const REGISTERED_IDENTITY_BYTES: usize = size_of::<(u8, String)>() + 32;

/// Reparto, término a término, de lo que la fase fundadora conserva viva.
///
/// Cada campo se compara contra el término de `final_formula_peak` que le
/// corresponde: agregar los cuatro antes de comparar dejaría que un término
/// sobrado tapase a otro corto, que es exactamente lo que ocurría con `64*H`.
#[derive(Debug)]
struct FounderTerms {
    /// Copias `FCP_` con su pertenencia y su texto: `192*X`.
    copies: u64,
    /// Filas haplotípicas retenidas: `B*N_pop*(160+S_bound)*H`.
    haplotype_rows: u64,
    /// Asientos de `F5FounderUnit::assignments`, 32 B cada uno: `64*B*C_total`.
    seats: u64,
    /// Registro de colisiones: `128*(X+B*N_pop*H+B*C_total)`.
    registry: u64,
}

impl FounderTerms {
    fn total(&self) -> u64 {
        self.copies + self.haplotype_rows + self.seats + self.registry
    }

    /// Ningún término vivo puede exceder el que v8 le concede.
    fn assert_within(&self, budget: &Self, label: &str) {
        for (name, live, granted) in [
            ("copias", self.copies, budget.copies),
            (
                "filas haplotípicas",
                self.haplotype_rows,
                budget.haplotype_rows,
            ),
            ("asientos", self.seats, budget.seats),
            ("registro", self.registry, budget.registry),
        ] {
            assert!(
                live <= granted,
                "{label}: {live} bytes vivos de {name} exceden los {granted} presupuestados"
            );
        }
    }
}

/// Lo que la fase fundadora conserva viva a la vez, contado por el tamaño real
/// de cada tipo más la única asignación de cada texto de identidad.
///
/// Es una **cota inferior** declarada, como [`retained_bytes`]: no cuenta el
/// desbordamiento de nodo de los árboles ni el redondeo del asignador. Ese es
/// el lado honesto para compararla contra un presupuesto que sí debe
/// cubrirlos.
fn founder_footprint(
    identities: &F5FounderIdentities,
    assignment: &F5FounderAssignment,
) -> FounderTerms {
    let copies = identities.copies().len()
        * (size_of::<F5FounderCopy>()
            + size_of::<FounderCopyId>()
            + size_of::<usize>()
            + SHARED_IDENTIFIER_BYTES);
    let haplotype_rows: usize = identities
        .units()
        .iter()
        .map(|unit| {
            size_of::<F5UnitHaplotypes>()
                + size_of::<usize>()
                + unit.haplotype_ids().len()
                    * (size_of::<F5HaplotypeId>() + SHARED_IDENTIFIER_BYTES)
        })
        .sum();
    let seats: usize = assignment
        .units()
        .iter()
        .map(|unit| size_of_val(unit.assignments()))
        .sum();
    FounderTerms {
        copies: copies as u64,
        haplotype_rows: haplotype_rows as u64,
        seats: seats as u64,
        registry: (identities.registered_identities() * REGISTERED_IDENTITY_BYTES) as u64,
    }
}

/// Bytes que `f5-resource-accounting-v9` concede a esas mismas estructuras,
/// escritos desde la fórmula de DTC-089 y no desde el código que la evalúa.
fn founder_budget(
    chain: &Chain,
    populations: u64,
    haplotypes: u64,
    founder_count: u64,
    ploidy: u64,
    sequences: u64,
) -> FounderTerms {
    let blocks = chain.sealed.block_count();
    let site_bound = chain.sealed.site_bound_per_block();
    let c_total = founder_count * ploidy;
    let copy_sequences = c_total * sequences;
    let seats = blocks * c_total;
    FounderTerms {
        copies: 192 * copy_sequences,
        haplotype_rows: blocks * populations * (160 + site_bound) * haplotypes,
        seats: 64 * seats,
        registry: 128 * (copy_sequences + blocks * populations * haplotypes + seats),
    }
}

/// Lo que el repositorio retiene **no crece con el número de bloques**.
///
/// Es la relación que faltaba y que habría detectado el defecto sin revisión:
/// la versión anterior retenía `B*C_total` asignaciones de 136 bytes y
/// `B*N_pop*H` filas haplotípicas con su cadena de alelos, de modo que
/// multiplicar `B` por cuatro multiplicaba lo retenido por casi cuatro mientras
/// la admisión sellaba el mismo presupuesto.
///
/// Ahora el único término que crece con `B` son los `H` identificadores `HAP_`
/// de cada unidad, y la prueba fija su tamaño **exacto**. Las asignaciones
/// `FHA_` se entregan una a una y no se retienen: su número sí crece con `B`.
#[test]
fn the_retained_repository_does_not_grow_with_the_block_count() {
    // Alambre de tropiezo: el repositorio conserva exactamente los campos que
    // no crecen con `B` —contexto de derivación, registro, copias, pertenencia,
    // unidades y su índice—. Retener otra tabla añadiría un campo y cambiaría
    // este tamaño, que `retained_bytes` no podría ver porque solo recorre la
    // superficie pública.
    assert_eq!(
        size_of::<F5FounderIdentities>(),
        232,
        "una tabla retenida nueva cambiaría el tamaño del repositorio"
    );
    let haplotypes = 2; // `H` del fixture por defecto.
    let mut measured = Vec::new();
    for blocks in [2_u64, 8] {
        let root = temporary(&format!("retained-{blocks}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", 2), Subpopulation::new("P2", 2)],
            block_count: blocks,
            sequences: 1,
            // Holgura para que los ocho bloques quepan; la geometría no
            // interviene en lo que se mide, solo `B`.
            slack: 2,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
        let mut identities = materialize_f5_founder_identities(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &assignment,
        )
        .expect("identidades fundadoras materializadas");

        let mut emitted = 0_usize;
        identities
            .for_each_assignment(&assignment, |_| {
                emitted += 1;
                Ok(())
            })
            .expect("recorrido de asignaciones");

        let retained = retained_bytes(&identities);
        // La promesa que la admisión selló para toda la fase.
        assert!(
            retained as u64 <= chain.sealed.working_memory_bytes(),
            "B={blocks}: {retained} bytes retenidos exceden el pico sellado {}",
            chain.sealed.working_memory_bytes()
        );
        // La relación que el paso 2 de B2 fija: todo lo que la fase fundadora
        // conserva viva —tablas, textos, registro y asientos— cabe dentro de
        // los términos que `final_formula_peak` le concede.
        founder_footprint(&identities, &assignment).assert_within(
            &founder_budget(&chain, 2, haplotypes as u64, 4, 2, 1),
            &format!("B={blocks}"),
        );
        measured.push((
            blocks,
            identities.copies().len(),
            identities.units().len(),
            retained,
            emitted,
            identities.registered_identities(),
        ));
    }

    let (_, copies_small, units_small, retained_small, emitted_small, registered_small) =
        measured[0];
    let (_, copies_large, units_large, retained_large, emitted_large, registered_large) =
        measured[1];

    assert_eq!(
        copies_small, copies_large,
        "las copias no dependen del número de bloques"
    );
    assert_eq!(copies_small, 8, "dos subpoblaciones de dos individuos, P=2");
    assert_eq!(units_small, 2 * 2, "B×N_pop unidades");
    assert_eq!(units_large, 8 * 2);

    // El único crecimiento es el de las unidades, y es exactamente el suyo.
    let per_unit = size_of::<F5UnitHaplotypes>()
        + haplotypes * size_of::<F5HaplotypeId>()
        + size_of::<usize>();
    assert_eq!(
        retained_large - retained_small,
        (units_large - units_small) * per_unit,
        "lo retenido solo puede crecer por los HAP_ de cada unidad"
    );

    // Las asignaciones sí crecen como `B*C_total`, y por eso no se retienen.
    assert_eq!(emitted_small, 2 * copies_small, "B=2 por C_total=8");
    assert_eq!(emitted_large, 8 * copies_large, "B=8 por C_total=8");
    assert_eq!(
        emitted_large,
        4 * emitted_small,
        "cuadruplicar B cuadruplica las asignaciones emitidas"
    );

    // El término dominante: el registro de colisiones conserva una entrada por
    // identidad derivada y crece como `B*C_total`.
    // `f5-resource-accounting-v9` lo presupuesta como
    // `128*(X+B*N_pop*H+B*C_total)`; v7 no lo contenía en ninguna fórmula.
    assert_eq!(
        registered_small,
        copies_small + units_small * haplotypes + emitted_small,
        "FCP_ por copia, HAP_ por rango y FHA_ por asignación"
    );
    assert_eq!(
        registered_large,
        copies_large + units_large * haplotypes + emitted_large
    );
    assert!(
        registered_large > 3 * registered_small,
        "el registro crece con B aunque las tablas ya no lo hagan"
    );
}

/// El presupuesto sellado **crece exactamente** con las estructuras
/// fundadoras que la fase retiene, y no por casualidad.
///
/// La ploidía es la única dimensión que mueve `C_total` —y con él `X`, los
/// asientos de asignación y las entradas del registro de identidades— sin
/// tocar ninguna otra magnitud de `final_formula_peak`: bloques, tramos, IDs
/// de secuencia, filas de perfil, normalizaciones, `H`, `S_bound` y bins
/// quedan idénticos. La diferencia entre los dos picos sellados es por tanto
/// la de los tres términos fundadores y de nada más, y esta prueba la escribe
/// desde DTC-089 en lugar de leerla del código que la calcula.
///
/// Con v7 esta igualdad era `128*ΔC_total`: no contaba las copias por
/// secuencia, ni los asientos, ni el registro. Suprimir hoy cualquiera de los
/// tres términos, o devolver `192*X` a `128*C_total`, la mata.
#[test]
fn the_sealed_budget_grows_exactly_with_the_founder_structures_it_retains() {
    let mut sealed = Vec::new();
    let mut footprints = Vec::new();
    for ploidy in [2_u64, 3] {
        let root = temporary(&format!("budget-ploidy-{ploidy}"));
        let fixture = Fixture {
            subpopulations: vec![Subpopulation::new("P1", 2), Subpopulation::new("P2", 2)],
            block_count: 2,
            sequences: 1,
            slack: 2,
            ploidy,
            ..Fixture::default()
        };
        let chain = build_chain(&root, &fixture);
        let assignment = assign(&chain);
        let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
        let mut identities = materialize_f5_founder_identities(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &assignment,
        )
        .expect("identidades fundadoras materializadas");
        identities
            .for_each_assignment(&assignment, |_| Ok(()))
            .expect("recorrido de asignaciones");

        // Lo retenido cabe en lo presupuestado también con ploidía impar, y
        // término a término.
        let footprint = founder_footprint(&identities, &assignment);
        footprint.assert_within(
            &founder_budget(&chain, 2, 2, 4, ploidy, 1),
            &format!("ploidía {ploidy}"),
        );
        footprints.push(footprint.total());
        sealed.push((
            chain.sealed.working_memory_bytes(),
            chain.sealed.block_count(),
            chain.sealed.site_bound_per_block(),
            chain.sealed.retained_sequence_ids(),
        ));
    }

    let (diploid_peak, blocks, site_bound, sequences) = sealed[0];
    let (triploid_peak, blocks_again, site_bound_again, sequences_again) = sealed[1];
    assert_eq!(
        blocks, blocks_again,
        "la ploidía no mueve el número de bloques"
    );
    assert_eq!(site_bound, site_bound_again, "ni la cota de sitios");
    assert_eq!(sequences, sequences_again, "ni las secuencias retenidas");

    // Cuatro individuos, una secuencia: la copia extra de cada individuo añade
    // `founder_count` copias a `C_total`, otras tantas a `X` y una por bloque
    // a los asientos, y cada una de esas identidades entra en el registro.
    let founder_count = 4_u64;
    let extra_copies = founder_count * sequences;
    let extra_seats = blocks * founder_count;
    let expected = 192 * extra_copies + 64 * extra_seats + 128 * (extra_copies + extra_seats);
    assert_eq!(
        triploid_peak - diploid_peak,
        expected,
        "el pico sellado debe crecer exactamente por los términos fundadores"
    );
    assert!(
        footprints[1] > footprints[0],
        "la tercera copia de cada individuo retiene más de lo que retenían dos"
    );
    assert!(
        footprints[1] - footprints[0] <= expected,
        "lo que la ploidía añade de verdad no puede exceder lo que v8 le concede"
    );
}

/// El recorrido de asignaciones es **único**: repetirlo derivaría cada `FHA_`
/// por segunda vez, que es una colisión.
#[test]
fn the_assignment_traversal_can_only_be_consumed_once() {
    let root = temporary("single-traversal");
    let chain = build_chain(&root, &Fixture::default());
    let assignment = assign(&chain);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let mut identities = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &execution,
        &chain.layout,
        &chain.search,
        &assignment,
    )
    .expect("identidades fundadoras materializadas");

    identities
        .for_each_assignment(&assignment, |_| Ok(()))
        .expect("el primer recorrido materializa las asignaciones");
    let error = identities
        .for_each_assignment(&assignment, |_| Ok(()))
        .expect_err("el segundo recorrido no puede repetir la derivación");
    assert_eq!(error.code(), "GPS-INTERNAL-101");
    assert!(format!("{error:?}").contains("ya se materializaron"));
}

/// Una asignación que no es la que materializó el repositorio se rechaza antes
/// de derivar ningún `FHA_`.
#[test]
fn a_traversal_with_another_assignment_is_rejected() {
    let root = temporary("foreign-traversal");
    let chain = build_chain(
        &root,
        &Fixture {
            block_count: 2,
            ..Fixture::default()
        },
    );
    let other_root = temporary("foreign-traversal-other");
    let other = build_chain(
        &other_root,
        &Fixture {
            subpopulations: vec![Subpopulation::new("Q1", 1)],
            block_count: 2,
            ..Fixture::default()
        },
    );
    let assignment = assign(&chain);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let mut identities = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &execution,
        &chain.layout,
        &chain.search,
        &assignment,
    )
    .expect("identidades fundadoras materializadas");

    let foreign = assign(&other);
    let error = identities
        .for_each_assignment(&foreign, |_| {
            panic!("no puede visitarse ninguna fila");
        })
        .expect_err("la asignación es de otra configuración");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(
        error.context().get("fields").map(String::as_str),
        Some("assignment.identity")
    );
    // El recorrido no se consumió: la asignación propia sigue siendo válida.
    let mut emitted = 0_usize;
    identities
        .for_each_assignment(&assignment, |_| {
            emitted += 1;
            Ok(())
        })
        .expect("el recorrido propio sigue disponible");
    assert_eq!(emitted, 2 * 2);
}

/// El error del visitante se propaga sin envolverlo y detiene el recorrido.
#[test]
fn the_visitor_error_stops_the_traversal() {
    let root = temporary("visitor-error");
    let chain = build_chain(&root, &Fixture::default());
    let assignment = assign(&chain);
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let mut identities = materialize_f5_founder_identities(
        &chain.config,
        &chain.sealed,
        &execution,
        &chain.layout,
        &chain.search,
        &assignment,
    )
    .expect("identidades fundadoras materializadas");

    let mut visited = 0_usize;
    let error = identities
        .for_each_assignment(&assignment, |_| {
            visited += 1;
            if visited == 2 {
                return Err(Diagnostic::f5_internal("scenario", "el escritor falló"));
            }
            Ok(())
        })
        .expect_err("el error del visitante llega al llamador");
    assert_eq!(error.code(), "GPS-INTERNAL-101");
    assert!(format!("{error:?}").contains("el escritor falló"));
    assert_eq!(visited, 2, "el recorrido se detiene en la fila que falló");
}
