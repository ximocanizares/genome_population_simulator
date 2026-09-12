//! F5.4.2. Oráculo independiente de aceptación sobre los bytes cerrados.
//!
//! La batería no comprueba el oráculo contra el productor: publica un paquete
//! real, construye la **declaración** que F5.4.3 copiará al manifiesto y
//! comprueba que el oráculo la sostiene o la tumba. Cada mutación —de los
//! bytes o de la declaración— se contrasta contra la causa tipada que le
//! corresponde, que es el criterio literal del hito: no basta con que falle,
//! tiene que fallar **por su causa**.
//!
//! Los oráculos internos —identidades DTC-037/DTC-039 y subsemillas
//! DTC-089— se validan contra los catálogos congelados en las unitarias de
//! `src/f5_oracle.rs`, antes de que este fichero los use para nada.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_artifacts::{
    F5_TSV_CONTRACTS, F5ReferenceSequences, F5SuccessArtifacts, F5SuccessRequest,
    write_f5_success_artifacts,
};
use genome_population_simulator::f5_founders::{
    F5FounderAssignment, F5FounderIdentities, F5RealizedFStatus, F5SamplingMode,
    assign_f5_founders, materialize_f5_founder_identities,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_oracle::{
    F5_RNG_CATALOG_DECISIONS, F5DeclaredAttemptSummary, F5DeclaredBounds, F5DeclaredLdBin,
    F5DeclaredLink, F5DeclaredMode, F5DeclaredRealizedF, F5DeclaredSampling,
    F5DeclaredSpectrumBand, F5DeclaredSubpopulation, F5DeclaredTruthLink, F5GenotypeModel,
    F5OracleCause, F5PackageDeclaration, audit_f5_package, recount_f5_man70,
};
use genome_population_simulator::f5_preflight::{
    F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, F5_DEFAULT_MAX_SWAP_PROPOSALS, F5OtherResources,
    SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_search::{F5HaplotypeSearch, search_f5_haplotypes};
use genome_population_simulator::identity::{ExecutionId, SequenceId};
use genome_population_simulator::publication::ExecutionWorkspace;
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
        "gps-f5-oracle-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

// ---------------------------------------------------------------------------
// Fixture de escenario
// ---------------------------------------------------------------------------

#[derive(Clone)]
struct Subpopulation<'a> {
    id: &'a str,
    individual_count: u64,
    sampling: Option<&'a str>,
    inbreeding_f: Option<&'a str>,
    keys: bool,
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
    const fn keys(mut self) -> Self {
        self.keys = true;
        self
    }
    const fn founders(mut self, founder_ids: &'a [&'a str]) -> Self {
        self.founder_ids = Some(founder_ids);
        self
    }
}

#[derive(Clone)]
struct Substitution<'a> {
    population: &'a str,
    block_index: u64,
    mode: &'a str,
    tolerance: Option<&'a str>,
}

struct Fixture<'a> {
    subpopulations: Vec<Subpopulation<'a>>,
    overrides: Vec<Substitution<'a>>,
    block_count: u64,
    block_length: u64,
    /// Secuencias que reciben bloques.
    sequences: u64,
    /// Secuencias adicionales compuestas solo de `N`: la referencia las retiene
    /// y el layout no las puebla, que es el caso que DTC-089 declara explícito.
    barren_sequences: u64,
    weights: Vec<&'a str>,
    gap: u64,
    spectrum: &'a str,
    bins: Vec<(u64, u64, f64, f64)>,
    max_haplotype_attempts: u64,
    max_swap_proposals: u64,
    ploidy: u64,
    seed: &'a str,
    slack: u64,
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
            barren_sequences: 0,
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
    /// Escribe los documentos y devuelve la ruta del general y las longitudes
    /// de las secuencias, que el escritor necesita para `truth.jsonl`.
    fn write(&self, root: &Path) -> (PathBuf, Vec<(String, u64)>) {
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
            "cada secuencia poblada debe recibir el mismo número entero de bloques"
        );
        let span =
            usize::try_from(self.block_length * (self.block_count / self.sequences) * self.slack)
                .expect("span");
        let mut fasta = String::new();
        let mut lengths = Vec::new();
        for ordinal in 0..self.sequences {
            let sequence = (0..span)
                .map(|index| char::from(bases[(index + ordinal as usize) % 4]))
                .collect::<String>();
            let id = format!("chr{}", ordinal + 1);
            fasta.push_str(&format!(">{id}\n{sequence}\n"));
            lengths.push((id, span as u64));
        }
        for ordinal in 0..self.barren_sequences {
            // Sin base elegible no hay tramo y, por tanto, ningún bloque.
            let id = format!("barren{}", ordinal + 1);
            fasta.push_str(&format!(">{id}\n{}\n", "N".repeat(span)));
            lengths.push((id, span as u64));
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
        (general, lengths)
    }
}

// ---------------------------------------------------------------------------
// Cadena productiva y publicación
// ---------------------------------------------------------------------------

struct Chain {
    config: EffectiveF5Config,
    sealed: SealedF5Admission,
    layout: F5PhysicalLayout,
    search: F5HaplotypeSearch,
    assignment: F5FounderAssignment,
    identities: F5FounderIdentities,
    sequences: F5ReferenceSequences,
    root: PathBuf,
}

impl Drop for Chain {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn build_chain(root: &Path, fixture: &Fixture<'_>) -> Chain {
    let (general, lengths) = fixture.write(root);
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
    let assignment =
        assign_f5_founders(&config, &sealed, &search).expect("asignación fundadora completa");
    let identities = materialize_f5_founder_identities(
        &config,
        &sealed,
        &execution,
        &layout,
        &search,
        &assignment,
    )
    .expect("identidades fundadoras materializadas");
    let mut sequences = F5ReferenceSequences::new();
    for (id, length) in lengths {
        assert!(sequences.insert(SequenceId::parse(id, "test").expect("sequence"), length));
    }
    Chain {
        config,
        sealed,
        layout,
        search,
        assignment,
        identities,
        sequences,
        root: root.to_path_buf(),
    }
}

/// Un paquete F5.4.1 escrito en el provisional, con sus bytes accesibles.
struct Published {
    chain: Chain,
    payload: PathBuf,
    artifacts: F5SuccessArtifacts,
}

impl Published {
    fn bytes(&self, relative: &str) -> Vec<u8> {
        fs::read(self.payload.join(relative)).unwrap_or_else(|_| panic!("artefacto {relative}"))
    }
}

fn payload_of(root: &Path, destination: &Path) -> PathBuf {
    let name = destination.file_name().expect("name").to_string_lossy();
    fs::read_dir(root)
        .expect("root")
        .map(|entry| entry.expect("entry").path())
        .find(|path| {
            path.file_name().is_some_and(|value| {
                value
                    .to_string_lossy()
                    .starts_with(&format!(".{name}.gps-incomplete-"))
            })
        })
        .expect("provisional")
        .join("payload")
}

fn publish(label: &str, fixture: &Fixture<'_>) -> Published {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("root");
    let chain = build_chain(&root, fixture);
    publish_chain(chain)
}

fn publish_chain(mut chain: Chain) -> Published {
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let artifacts = {
        let request = F5SuccessRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
            &chain.sequences,
        )
        .expect("solicitud de publicación");
        write_f5_success_artifacts(&mut workspace, &request, &mut chain.identities)
            .expect("paquete F5.4.1 escrito")
    };
    let payload = payload_of(&chain.root, &destination);
    Published {
        chain,
        payload,
        artifacts,
    }
}

/// Directorio temporal que se borra aunque la prueba falle: sin él, una
/// aserción rota deja restos en `/tmp`.
struct TempDir(PathBuf);

impl TempDir {
    fn new(label: &str) -> Self {
        let root = temporary(label);
        fs::create_dir_all(&root).expect("destino temporal");
        Self(root)
    }
    fn path(&self) -> &Path {
        &self.0
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

// ---------------------------------------------------------------------------
// La declaración que el oráculo juzga
// ---------------------------------------------------------------------------

fn declared_mode(mode: F5SamplingMode) -> F5DeclaredMode {
    match mode {
        F5SamplingMode::Probabilistic => F5DeclaredMode::Probabilistic,
        F5SamplingMode::Controlled => F5DeclaredMode::Controlled,
    }
}

/// Construye la declaración exactamente como F5.4.3 la copiará al manifiesto.
///
/// El oráculo no recibe ninguna de estas estructuras: recibe sus cifras, y las
/// contrasta contra los bytes.
fn declaration_for(published: &Published) -> F5PackageDeclaration {
    let chain = &published.chain;
    let config = &chain.config;
    let sealed = &chain.sealed;
    let profile = config.profile();
    let scenario = config.declared_scenario();
    let mut subpopulations = Vec::new();
    for roster in chain.assignment.rosters() {
        let inbreeding_f = chain
            .assignment
            .f_realized()
            .iter()
            .find(|row| row.population_id() == roster.population_id())
            .and_then(genome_population_simulator::f5_founders::F5RealizedF::target_f);
        subpopulations.push(F5DeclaredSubpopulation {
            population_id: roster.population_id().as_str().to_owned(),
            individual_count: roster.individual_count(),
            inbreeding_f,
        });
    }
    let sampling = chain
        .assignment
        .sampling()
        .iter()
        .map(|row| F5DeclaredSampling {
            population_id: row.population_id().as_str().to_owned(),
            block_id: row.block_id().as_str().to_owned(),
            requested_mode: declared_mode(row.requested_mode()),
            effective_mode: declared_mode(row.effective_mode()),
            controlled_max_abs_error_declared: row.controlled_max_abs_error_declared(),
            controlled_max_abs_error_effective: row.controlled_max_abs_error_effective(),
        })
        .collect();
    let f_realized = chain
        .assignment
        .f_realized()
        .iter()
        .map(|row| F5DeclaredRealizedF {
            population_id: row.population_id().as_str().to_owned(),
            block_id: row.block_id().as_str().to_owned(),
            target_f: row.target_f(),
            status: match row.status() {
                F5RealizedFStatus::Evaluated => "evaluated".to_owned(),
                F5RealizedFStatus::NotApplicable => "not_applicable".to_owned(),
            },
            value: row.value(),
        })
        .collect();
    let mut ld_bins = Vec::new();
    for subpopulation in &subpopulations {
        for bin in profile
            .ld_bins(&subpopulation.population_id)
            .expect("bins LD del perfil")
        {
            ld_bins.push(F5DeclaredLdBin {
                population_id: subpopulation.population_id.clone(),
                distance_start_bp: bin.distance_start_bp(),
                distance_end_bp: bin.distance_end_bp(),
                target_mean_r2: bin.target_mean_r2(),
                max_abs_error: bin.max_abs_error(),
            });
        }
    }
    let artifacts = published
        .artifacts
        .artifacts()
        .iter()
        .map(|link| F5DeclaredLink {
            artifact_id: link.artifact_id().to_owned(),
            artifact_kind: link.artifact_kind().to_owned(),
            schema_version: link.schema_version().to_owned(),
            path: link.path().to_owned(),
            size_bytes: link.size_bytes(),
            sha256: link.sha256().to_owned(),
        })
        .collect();
    let truth = published.artifacts.truth();
    F5PackageDeclaration {
        scenario_id: config.scenario_id().to_owned(),
        execution_id: EXECUTION.to_owned(),
        reference_id: config.reference_id().to_owned(),
        reference_sha256: sealed.reference_sha256().to_owned(),
        seed: config.seed().to_owned(),
        ploidy: config.ploidy(),
        profile_id: sealed.profile_id().to_owned(),
        profile_metadata_sha256: sealed.profile_metadata_sha256().to_owned(),
        haplotype_weights: profile.haplotype_weights().to_vec(),
        snp_spectrum: profile
            .snp_spectrum()
            .iter()
            .map(|band| F5DeclaredSpectrumBand {
                alt_af_min: band.minimum(),
                alt_af_max: band.maximum(),
                weight: band.weight(),
            })
            .collect(),
        subpopulations,
        sampling,
        f_realized,
        ld_bins,
        normalizations: profile.normalizations().to_vec(),
        artifacts,
        truth: F5DeclaredTruthLink {
            artifact_id: truth.artifact_id().to_owned(),
            artifact_kind: truth.artifact_kind().to_owned(),
            schema_version: truth.schema_version().to_owned(),
            projection_version: truth.projection_version().to_owned(),
            path: truth.path().to_owned(),
            size_bytes: truth.size_bytes(),
            sha256: truth.sha256().to_owned(),
        },
        bounds: F5DeclaredBounds {
            block_count: sealed.block_count(),
            eligible_bases: sealed.eligible_bases(),
            site_bound_per_block: sealed.site_bound_per_block(),
            bin_count_total: sealed.bin_count_total(),
            max_haplotype_attempts: scenario
                .max_haplotype_attempts
                .as_ref()
                .map_or(F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, |value| *value.value()),
            max_swap_proposals: scenario
                .max_swap_proposals
                .as_ref()
                .map_or(F5_DEFAULT_MAX_SWAP_PROPOSALS, |value| *value.value()),
            tsv_rows_bound: sealed.tsv_rows_bound(),
            json_rows_bound: sealed.json_rows_bound(),
            output_bytes: sealed.output_bytes(),
            output_limit_bytes: sealed.output_limit_bytes(),
        },
        id_catalog_id: "id-catalog-v4".to_owned(),
        rng_catalog_id: "f5-rng-catalog-v1".to_owned(),
        attempt_summary: F5DeclaredAttemptSummary {
            layout_attempts_used: chain.layout.successful_attempt_index() + 1,
            haplotype_attempts_total: chain.search.haplotype_attempts_total(),
            max_attempt_index: chain.search.max_attempt_index(),
        },
    }
}

/// Copia el paquete publicado a un destino propio y mutable.
struct Staged {
    root: PathBuf,
    declaration: F5PackageDeclaration,
}

impl Drop for Staged {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Staged {
    fn new(published: &Published, label: &str) -> Self {
        let root = temporary(label);
        fs::create_dir_all(root.join("outputs")).expect("destino");
        for contract in &F5_TSV_CONTRACTS {
            fs::write(root.join(contract.path()), published.bytes(contract.path()))
                .expect("copia del artefacto");
        }
        fs::write(
            root.join("outputs/truth.jsonl"),
            published.bytes("outputs/truth.jsonl"),
        )
        .expect("copia de la proyección");
        Self {
            root,
            declaration: declaration_for(published),
        }
    }

    fn text(&self, relative: &str) -> String {
        fs::read_to_string(self.root.join(relative)).expect("artefacto legible")
    }

    /// Reescribe un artefacto **y su enlace**, de modo que la mutación se
    /// juzgue por su semántica y no por la huella que ya la delataría.
    fn rewrite(&mut self, relative: &str, contents: &str) {
        fs::write(self.root.join(relative), contents).expect("artefacto reescrito");
        let size = contents.len() as u64;
        let digest = sha256(contents.as_bytes());
        if relative == "outputs/truth.jsonl" {
            self.declaration.truth.size_bytes = size;
            self.declaration.truth.sha256 = digest;
            return;
        }
        let link = self
            .declaration
            .artifacts
            .iter_mut()
            .find(|link| link.path == relative)
            .expect("enlace del artefacto");
        link.size_bytes = size;
        link.sha256 = digest;
    }

    fn audit(&self) -> Result<(), F5OracleCause> {
        audit_f5_package(&self.root, &self.declaration)
            .map(|_| ())
            .map_err(|rejection| rejection.cause())
    }
}

// ---------------------------------------------------------------------------
// Fixtures de cobertura
// ---------------------------------------------------------------------------

/// Dos subpoblaciones con censos distintos, dos secuencias con dos bloques
/// cada una, una secuencia estéril que la referencia retiene, ploidía 2,
/// modos mixtos y fundadores aportados frente a derivados.
fn wide_fixture<'a>() -> Fixture<'a> {
    Fixture {
        subpopulations: vec![
            Subpopulation::new("PA", 2).founders(&["zeta", "alfa"]),
            Subpopulation::new("PB", 1).keys().sampling("controlled"),
        ],
        overrides: vec![Substitution {
            population: "PA",
            block_index: 2,
            mode: "controlled",
            tolerance: None,
        }],
        block_count: 4,
        sequences: 2,
        barren_sequences: 1,
        block_length: 3,
        weights: vec!["0.5", "0.5"],
        bins: vec![(1, 1, 1.0, 0.05), (2, 2, 1.0, 0.05)],
        max_haplotype_attempts: 256,
        max_swap_proposals: 1024,
        ..Fixture::default()
    }
}

/// Bloques con hueco 2: las distancias realizadas son 2 y 4, de modo que dos
/// bins del perfil quedan `not_applicable` y dos `evaluated`.
fn sparse_fixture<'a>() -> Fixture<'a> {
    Fixture {
        block_length: 5,
        gap: 2,
        block_count: 2,
        bins: vec![
            (1, 1, 1.0, 0.05),
            (2, 2, 1.0, 0.05),
            (3, 3, 1.0, 0.05),
            (4, 4, 1.0, 0.05),
        ],
        max_haplotype_attempts: 512,
        max_swap_proposals: 2048,
        ..Fixture::default()
    }
}

/// Perfil diádico de cuatro pesos distintos: `1/2, 1/4, 1/8, 1/8`.
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

// ---------------------------------------------------------------------------
// Aceptación
// ---------------------------------------------------------------------------

#[test]
fn the_oracle_accepts_a_whole_published_package_and_summarises_it() {
    let published = publish("accept", &wide_fixture());
    let staged = Staged::new(&published, "accept-stage");
    let report = audit_f5_package(&staged.root, &staged.declaration).expect("paquete válido");
    assert_eq!(report.copies(), published.artifacts.copy_count());
    assert_eq!(report.genotypes(), published.artifacts.genotype_count());
    assert!(!report.initialization_summaries().is_empty());
    for summary in report.initialization_summaries() {
        assert!(summary.recount_matches());
        assert_eq!(
            summary.recount_source(),
            "outputs/truth.jsonl#f5-founder-truth-projection-v2"
        );
        assert_eq!(
            summary.realized_genotype_counts(),
            summary.recounted_genotype_counts()
        );
        assert_eq!(
            summary.realized_genotype_counts().len(),
            usize::from(summary.ploidy()) + 1
        );
        assert_eq!(
            summary.realized_genotype_counts().iter().sum::<u64>(),
            summary.individual_count()
        );
        let probabilities: f64 = summary.expected_genotype_probabilities().iter().sum();
        assert!((probabilities - 1.0).abs() < 1e-9, "{probabilities}");
        assert_eq!(
            summary.expected_genotype_frequencies(),
            summary.expected_genotype_probabilities()
        );
        let expected_total: f64 = summary.expected_genotype_counts().iter().sum();
        assert!(
            (expected_total - summary.individual_count() as f64).abs() < 1e-9,
            "{expected_total}"
        );
        let recounted: f64 = summary.recounted_genotype_frequencies().iter().sum();
        assert!((recounted - 1.0).abs() < 1e-12, "{recounted}");
    }
    // Los dos modos del escenario producen los dos modelos que les tocan.
    let models: Vec<F5GenotypeModel> = report
        .initialization_summaries()
        .iter()
        .map(genome_population_simulator::f5_oracle::F5InitializationSummary::model)
        .collect();
    assert!(models.contains(&F5GenotypeModel::Hwe), "{models:?}");
    assert!(
        models.contains(&F5GenotypeModel::ControlledHweConditionedCounts),
        "{models:?}"
    );
    // Las cotas consumidas caben en las selladas.
    assert!(report.resources().tsv_rows() <= staged.declaration.bounds.tsv_rows_bound);
    // Cifra exacta del fixture ancho: 3 fundadores + 2 subpoblaciones + 8
    // filas de muestreo + 8 de F realizado + 24 resúmenes MAN-70
    // (2 subpoblaciones × 4 bloques × 3 sitios) + 2×5 filas de normalización
    // (1+1+1 intervalos y 2 pesos) + 9 objetos fijos.
    assert_eq!(
        report.resources().observable_json_rows(),
        3 + 2 + 8 + 8 + 24 + 10 + 9
    );
    assert!(report.resources().observable_json_rows() <= staged.declaration.bounds.json_rows_bound);
    assert_eq!(
        report.resources().truth_records(),
        1 + report.copies() + report.genotypes()
    );
    assert!(report.resources().output_bytes() <= staged.declaration.bounds.output_bytes);
}

/// El recuento MAN-70 abre `truth.jsonl` y nada más: se ejecuta sobre un
/// directorio donde **no existe ningún TSV ni ningún manifiesto**.
#[test]
fn the_man70_recount_needs_only_the_founder_truth() {
    let published = publish("recount", &wide_fixture());
    let directory = TempDir::new("recount-alone");
    let root = directory.path();
    let truth = root.join("truth.jsonl");
    fs::write(&truth, published.bytes("outputs/truth.jsonl")).expect("copia");
    let recount = recount_f5_man70(&truth).expect("recuento independiente");
    assert_eq!(recount.copies(), published.artifacts.copy_count());
    assert_eq!(recount.genotypes(), published.artifacts.genotype_count());
    assert_eq!(recount.execution_id(), EXECUTION);
    assert!(!recount.loci().is_empty());
    for locus in recount.loci() {
        assert_eq!(locus.counts().len(), usize::from(recount.ploidy()) + 1);
        assert_eq!(locus.counts().iter().sum::<u64>(), locus.individual_count());
    }
    // El recuento no lee ninguna otra ruta: el directorio solo tiene la verdad.
    let entries: Vec<String> = fs::read_dir(root)
        .expect("directorio")
        .map(|entry| {
            entry
                .expect("entrada")
                .file_name()
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    assert_eq!(entries, vec!["truth.jsonl".to_owned()]);
}

// ---------------------------------------------------------------------------
// Mutación aislada de los bytes
// ---------------------------------------------------------------------------

fn replace_once(text: &str, from: &str, to: &str) -> String {
    let at = text
        .find(from)
        .unwrap_or_else(|| panic!("el fixture no contiene {from:?}"));
    format!("{}{to}{}", &text[..at], &text[at + from.len()..])
}

fn replace_after(text: &str, anchor: &str, from: &str, to: &str) -> String {
    // El ancla se consume: buscar `from` desde su inicio acertaría dentro de la
    // propia clave y dejaría el valor intacto, que es el falso negativo que la
    // revisión de F5.4.2 destapó en «case_id ajeno».
    let start = text
        .find(anchor)
        .unwrap_or_else(|| panic!("el fixture no contiene el ancla {anchor:?}"))
        + anchor.len();
    let at = start
        + text[start..]
            .find(from)
            .unwrap_or_else(|| panic!("el fixture no contiene {from:?} tras {anchor:?}"));
    format!("{}{to}{}", &text[..at], &text[at + from.len()..])
}

/// Suma uno al primer entero decimal que aparece tras el ancla.
fn bump_number_after(text: &str, anchor: &str) -> String {
    let start = text.find(anchor).expect("ancla") + anchor.len();
    let end = start
        + text[start..]
            .find(|character: char| !character.is_ascii_digit())
            .expect("fin del entero");
    let value: u64 = text[start..end].parse().expect("entero");
    format!("{}{}{}", &text[..start], value + 1, &text[end..])
}

/// Cambia el dígito binario que sigue al ancla.
fn flip_digit_after(text: &str, anchor: &str) -> String {
    let at = text.find(anchor).expect("ancla") + anchor.len();
    let flipped = if text.as_bytes()[at] == b'0' {
        '1'
    } else {
        '0'
    };
    format!("{}{flipped}{}", &text[..at], &text[at + 1..])
}

/// Rota la base canónica que sigue al ancla.
fn rotate_base_after(text: &str, anchor: &str) -> String {
    let at = text.find(anchor).expect("ancla") + anchor.len();
    let rotated = match text.as_bytes()[at] {
        b'A' => 'C',
        b'C' => 'G',
        b'G' => 'T',
        _ => 'A',
    };
    format!("{}{rotated}{}", &text[..at], &text[at + 1..])
}

fn cell(text: &str, row: usize, column: usize) -> String {
    text.lines()
        .nth(row + 1)
        .expect("fila")
        .split('\t')
        .nth(column)
        .expect("columna")
        .to_owned()
}

fn set_cell(text: &str, row: usize, column: usize, value: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    let mut columns: Vec<String> = lines[row + 1].split('\t').map(str::to_owned).collect();
    columns[column] = value.to_owned();
    lines[row + 1] = columns.join("\t");
    format!("{}\n", lines.join("\n"))
}

fn swap_rows(text: &str, left: usize, right: usize) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    lines.swap(left + 1, right + 1);
    format!("{}\n", lines.join("\n"))
}

/// Muta un carácter Base32 del identificador de una celda.
fn damage_identifier(text: &str, row: usize, column: usize) -> String {
    let value = cell(text, row, column);
    let mut bytes = value.into_bytes();
    let last = bytes.len() - 1;
    bytes[last] = if bytes[last] == b'A' { b'B' } else { b'A' };
    set_cell(text, row, column, &String::from_utf8(bytes).expect("UTF-8"))
}

const LAYOUT_PATH: &str = "outputs/f5-block-layout.tsv";
const SITES_PATH: &str = "outputs/f5-sites.tsv";
const HAPLOTYPES_PATH: &str = "outputs/f5-haplotypes.tsv";
const FREQUENCIES_PATH: &str = "outputs/f5-frequencies.tsv";
const LD_PATH: &str = "outputs/f5-ld-evaluation.tsv";
const FOUNDERS_PATH: &str = "outputs/f5-founders.tsv";
const ATTEMPTS_PATH: &str = "outputs/f5-attempts.tsv";
const TRUTH_PATH: &str = "outputs/truth.jsonl";

/// Una corrupción aislada y la causa **exacta** con la que debe morir.
struct ByteMutation {
    name: &'static str,
    path: &'static str,
    cause: F5OracleCause,
    apply: fn(&str) -> String,
}

/// La primera base canónica que no es ninguna de las dos dadas.
fn other_base(reference: &str, alternate: &str) -> String {
    ["A", "C", "G", "T"]
        .into_iter()
        .find(|base| *base != reference && *base != alternate)
        .expect("una tercera base")
        .to_owned()
}

/// Borra la **última** fila por `site_id` de la primera unidad marginal. Antes
/// del cierre real de unidad nada lo detectaba directamente: moría de rebote
/// por `max_marginal_error` o por el `pair_count` de un bin, con causa ajena.
fn drop_last_site_of_first_unit(text: &str) -> String {
    let rows: Vec<&str> = text.lines().skip(1).collect();
    let first: Vec<&str> = rows[0].split('\t').collect();
    let last = rows
        .iter()
        .rposition(|row| {
            let fields: Vec<&str> = row.split('\t').collect();
            fields[0] == first[0] && fields[1] == first[1]
        })
        .expect("última fila de la unidad");
    assert!(
        last > 0,
        "fixture degenerado: la unidad tiene un solo sitio"
    );
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    lines.remove(last + 1);
    format!("{}\n", lines.join("\n"))
}

/// Mueve `target_alt_af` **y su error** dejando intacto el máximo de su unidad,
/// de modo que ni `abs_error` ni `max_marginal_error` del mejor intento lo
/// delatan: solo rehacer el sorteo `marginal_spectrum` lo hace.
fn move_target_marginal(text: &str) -> String {
    let rows: Vec<Vec<String>> = text
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect();
    let mut maximum: BTreeMap<(String, String), f64> = BTreeMap::new();
    for row in &rows {
        let error: f64 = row[6].parse().expect("error");
        let slot = maximum
            .entry((row[0].clone(), row[1].clone()))
            .or_insert(0.0);
        if error > *slot {
            *slot = error;
        }
    }
    for (index, row) in rows.iter().enumerate() {
        let error: f64 = row[6].parse().expect("error");
        let group = maximum[&(row[0].clone(), row[1].clone())];
        // El objetivo pasa a ser la frecuencia realizada: el error nuevo es
        // cero, estrictamente menor que el máximo de su unidad.
        if error > 0.0 && error < group {
            let text = set_cell(text, index, 4, &row[5]);
            return set_cell(&text, index, 6, &0.0_f64.to_string());
        }
    }
    panic!("fixture degenerado: ninguna unidad tiene un error marginal intermedio");
}

/// Mueve la frecuencia realizada **y su error** de forma coherente: la fila
/// queda internamente consistente y lo único que puede delatarla es rehacer
/// `Σ_h q_h x_hi` sobre la matriz publicada.
fn move_realized_marginal(text: &str) -> String {
    let target: f64 = cell(text, 0, 4).parse().expect("objetivo");
    let realized: f64 = cell(text, 0, 5).parse().expect("realizada");
    let moved = if realized > 0.5 {
        realized - 0.25
    } else {
        realized + 0.25
    };
    let error = (moved - target).abs();
    let text = set_cell(text, 0, 5, &moved.to_string());
    set_cell(&text, 0, 6, &error.to_string())
}

/// Mueve la media `r²` **y su error** dentro de la tolerancia del bin: la fila
/// sigue siendo coherente y solo rehacer `D` y `r²` pareja a pareja la delata.
fn move_mean_r2(text: &str) -> String {
    let row = text
        .lines()
        .skip(1)
        .position(|line| line.split('\t').nth(4) == Some("evaluated"))
        .expect("bin evaluado");
    let target: f64 = cell(text, row, 7).parse().expect("objetivo");
    let tolerance: f64 = cell(text, row, 8).parse().expect("tolerancia");
    let mean = target - tolerance / 2.0;
    let error = (mean - target).abs();
    assert!(
        error <= tolerance,
        "el desplazamiento debe caber en la tolerancia"
    );
    let text = set_cell(text, row, 6, &mean.to_string());
    set_cell(&text, row, 9, &error.to_string())
}

/// Lleva la posición de la **última** modificación de la primera copia fuera
/// de su secuencia. Es el único sitio donde la guarda
/// `position >= sequence_length` es alcanzable: las otras dos vistas se
/// comparan contra las modificaciones y morirían antes.
fn push_last_modification_out_of_the_sequence(text: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    let index = lines
        .iter()
        .position(|line| line.contains("\"record_type\":\"haplotype_copy\""))
        .expect("copia");
    let anchor = "\"position\":";
    let at = lines[index].rfind(anchor).expect("posición") + anchor.len();
    let end = at
        + lines[index][at..]
            .find(|character: char| !character.is_ascii_digit())
            .expect("fin de la posición");
    let position = lines[index][at..end].to_owned();
    // El mismo sitio aparece en las tres vistas con su posición; moverlas las
    // tres deja el registro coherente consigo mismo y solo la cota de la
    // secuencia lo delata. Las posiciones son únicas dentro de una copia y la
    // clave siguiente es siempre `provenance`, de modo que la sustitución no
    // puede acertar el prefijo de otra.
    let from = format!("\"position\":{position},\"provenance\"");
    assert_eq!(lines[index].matches(&from).count(), 3, "las tres vistas");
    lines[index] = lines[index].replace(&from, "\"position\":999999999,\"provenance\"");
    format!("{}\n", lines.join("\n"))
}

/// Duplica un portador de un genotipo heterocigoto y **ajusta sus dosis** para
/// que el registro cuadre consigo mismo.
///
/// Con la guarda de portador repetido en pie, muere por la verdad. Sin ella,
/// la única que lo mata es la igualdad entre lo realizado desde los TSV y lo
/// recontado desde la verdad: ahí es donde esa igualdad tiene dientes.
fn duplicate_one_carrier(text: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    let index = lines
        .iter()
        .position(|line| {
            if !line.contains("\"record_type\":\"genotype\"") {
                return false;
            }
            let record: serde_json::Value = serde_json::from_str(line).expect("genotipo");
            let carriers = record["carriers"].as_array().expect("portadores");
            carriers.len() == 2 && carriers[0]["allele_index"] != carriers[1]["allele_index"]
        })
        .expect("un genotipo heterocigoto");
    let mut record: serde_json::Value = serde_json::from_str(&lines[index]).expect("genotipo");
    let first = record["carriers"][0].clone();
    let allele_index = first["allele_index"].as_u64().expect("índice") as usize;
    record["carriers"][1] = first;
    let alleles = record["alleles"].as_array().expect("alelos").len();
    record["allele_dosages"] = serde_json::Value::Array(
        (0..alleles)
            .map(|position| {
                serde_json::Value::from(if position == allele_index { 2_u64 } else { 0 })
            })
            .collect(),
    );
    lines[index] = serde_json::to_string(&record).expect("serialización");
    format!("{}\n", lines.join("\n"))
}

/// Cambia el alelo alternativo de un SNP de la primera copia **en sus dos
/// vistas que lo publican** —modificaciones y fragmento vigente—, dejando
/// `f5-sites.tsv` intacto.
///
/// El registro queda internamente coherente, de modo que la comprobación de
/// las tres vistas no lo ve: solo el cruce con los bloques que
/// `f5-founders.tsv` asigna lo delata.
fn relabel_one_copy_snp(text: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    let index = lines
        .iter()
        .position(|line| line.contains("\"record_type\":\"haplotype_copy\""))
        .expect("copia");
    let anchor = "\"alternate\":\"";
    let at = lines[index].find(anchor).expect("alelo") + anchor.len();
    let alternate = lines[index][at..at + 1].to_owned();
    let event_anchor = "\"event_id\":\"";
    let event_at = lines[index][at..].find(event_anchor).expect("evento") + at + event_anchor.len();
    let event_end = lines[index][event_at..].find('"').expect("cierre") + event_at;
    let site_id = lines[index][event_at..event_end].to_owned();
    let reference_anchor = "\"reference\":\"";
    let reference_at = lines[index][at..]
        .find(reference_anchor)
        .expect("referencia")
        + at
        + reference_anchor.len();
    let reference = lines[index][reference_at..reference_at + 1].to_owned();
    let replacement = other_base(&reference, &alternate);
    let from = format!("\"alternate\":\"{alternate}\",\"event_id\":\"{site_id}\"");
    let to = format!("\"alternate\":\"{replacement}\",\"event_id\":\"{site_id}\"");
    assert_eq!(
        lines[index].matches(&from).count(),
        2,
        "el alelo debe aparecer en las dos vistas que lo publican"
    );
    lines[index] = lines[index].replace(&from, &to);
    format!("{}\n", lines.join("\n"))
}

/// Retira el último genotipo y ajusta el recuento de la cabecera, de modo que
/// el documento sigue cuadrando consigo mismo y lo único que falla es el
/// recuento por locus contra el censo declarado.
fn drop_last_genotype(text: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    lines.pop().expect("genotipo");
    let anchor = "\"genotypes\":";
    let at = lines[0].find(anchor).expect("recuento") + anchor.len();
    let end = at
        + lines[0][at..]
            .find(|character: char| !character.is_ascii_digit())
            .expect("fin del recuento");
    let value: u64 = lines[0][at..end].parse().expect("recuento");
    lines[0] = format!("{}{}{}", &lines[0][..at], value - 1, &lines[0][end..]);
    format!("{}\n", lines.join("\n"))
}

/// Cambia el `HAP_` de la traza en las **tres** vistas de la primera copia,
/// de modo que el registro queda internamente coherente y lo único que falla
/// es el cruce con la fila que `f5-founders.tsv` asignó.
fn retrace_first_copy(text: &str) -> String {
    let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
    let index = lines
        .iter()
        .position(|line| line.contains("\"record_type\":\"haplotype_copy\""))
        .expect("copia");
    let anchor = "\"matrix_source_id\":\"";
    let at = lines[index].find(anchor).expect("procedencia") + anchor.len();
    let end = lines[index][at..].find('"').expect("cierre") + at;
    let original = lines[index][at..end].to_owned();
    let mut bytes = original.clone().into_bytes();
    let last = bytes.len() - 1;
    bytes[last] = if bytes[last] == b'A' { b'B' } else { b'A' };
    let replacement = String::from_utf8(bytes).expect("UTF-8");
    lines[index] = lines[index].replace(&original, &replacement);
    format!("{}\n", lines.join("\n"))
}

/// Sustituye el haplotipo de una fila por el de otra fila **de su misma
/// unidad**: el `HAP_` existe en su tabla propietaria, de modo que lo que se
/// rompe es el `FHA_` y no la búsqueda del rango.
fn swap_haplotype_within_unit(text: &str) -> String {
    let rows = text.lines().count() - 1;
    for left in 0..rows {
        for right in (left + 1)..rows {
            if cell(text, left, 0) == cell(text, right, 0)
                && cell(text, left, 4) == cell(text, right, 4)
                && cell(text, left, 7) != cell(text, right, 7)
            {
                return set_cell(text, left, 7, &cell(text, right, 7));
            }
        }
    }
    panic!("fixture degenerado: ninguna unidad asigna dos filas haplotípicas distintas");
}

/// Cambia la polaridad publicada de la primera fila que no sea la del sitio 0.
fn flip_polarity(text: &str) -> String {
    let mut rows = 0_usize;
    loop {
        let value = cell(text, rows, 3);
        if value == "negative" {
            return set_cell(text, rows, 3, "positive");
        }
        if rows > 0 {
            return set_cell(text, rows, 3, "negative");
        }
        rows += 1;
    }
}

const BYTE_MUTATIONS: &[ByteMutation] = &[
    // --- contrato de bytes ---
    ByteMutation {
        name: "BOM al principio",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| format!("\u{feff}{text}"),
    },
    ByteMutation {
        name: "final de línea CRLF",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| replace_once(text, "\n", "\r\n"),
    },
    ByteMutation {
        name: "sin LF final",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| text.trim_end_matches('\n').to_owned(),
    },
    ByteMutation {
        name: "línea vacía",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| replace_once(text, "\n", "\n\n"),
    },
    ByteMutation {
        name: "comentario",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| replace_once(text, "\n", "\n# nota\n"),
    },
    ByteMutation {
        name: "columna de más",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| replace_once(text, "\n", "\textra\n"),
    },
    ByteMutation {
        name: "cabecera alterada",
        path: LAYOUT_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| replace_once(text, "block_id\t", "blockid\t"),
    },
    ByteMutation {
        name: "celda entrecomillada",
        path: SITES_PATH,
        cause: F5OracleCause::ByteContract,
        apply: |text| {
            let value = cell(text, 0, 5);
            set_cell(text, 0, 5, &format!("\"{value}\""))
        },
    },
    // --- layout ---
    ByteMutation {
        name: "layout con inicio movido",
        path: LAYOUT_PATH,
        cause: F5OracleCause::BlockLayout,
        apply: |text| {
            let start: u64 = cell(text, 0, 2).parse().expect("inicio");
            set_cell(text, 0, 2, &(start + 1).to_string())
        },
    },
    ByteMutation {
        name: "layout con longitud incoherente",
        path: LAYOUT_PATH,
        cause: F5OracleCause::BlockLayout,
        apply: |text| {
            let length: u64 = cell(text, 0, 4).parse().expect("longitud");
            set_cell(text, 0, 4, &(length + 1).to_string())
        },
    },
    ByteMutation {
        name: "layout desordenado",
        path: LAYOUT_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| swap_rows(text, 0, 1),
    },
    ByteMutation {
        name: "layout con BLK_ ajeno",
        path: LAYOUT_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| damage_identifier(text, 0, 0),
    },
    // --- sitios ---
    ByteMutation {
        name: "site_index no consecutivo",
        path: SITES_PATH,
        cause: F5OracleCause::Site,
        apply: |text| set_cell(text, 0, 1, "1"),
    },
    ByteMutation {
        name: "sitios desordenados",
        path: SITES_PATH,
        cause: F5OracleCause::Site,
        apply: |text| swap_rows(text, 0, 1),
    },
    ByteMutation {
        name: "posición del sitio movida",
        path: SITES_PATH,
        cause: F5OracleCause::Site,
        apply: |text| {
            let position: u64 = cell(text, 0, 4).parse().expect("posición");
            set_cell(text, 0, 4, &(position + 1).to_string())
        },
    },
    ByteMutation {
        name: "SITE_ ajeno",
        path: SITES_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| damage_identifier(text, 0, 2),
    },
    ByteMutation {
        name: "alelo alternativo ajeno al sorteo",
        path: SITES_PATH,
        cause: F5OracleCause::Site,
        apply: |text| {
            set_cell(
                text,
                0,
                6,
                &other_base(&cell(text, 0, 5), &cell(text, 0, 6)),
            )
        },
    },
    ByteMutation {
        name: "base de referencia cambiada",
        path: SITES_PATH,
        cause: F5OracleCause::Site,
        apply: |text| set_cell(text, 0, 5, &cell(text, 0, 6)),
    },
    // --- haplotipos: H×S y q ---
    ByteMutation {
        name: "q sin decimal mínimo",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Haplotype,
        apply: |text| set_cell(text, 0, 4, &format!("{}0", cell(text, 0, 4))),
    },
    ByteMutation {
        name: "q ajena al perfil",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Haplotype,
        apply: |text| set_cell(text, 0, 4, "0.25"),
    },
    ByteMutation {
        name: "rango haplotípico repetido",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Haplotype,
        apply: |text| set_cell(text, 1, 3, "0"),
    },
    ByteMutation {
        name: "cadena de alelos truncada",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Haplotype,
        apply: |text| {
            let alleles = cell(text, 0, 5);
            set_cell(text, 0, 5, &alleles[..alleles.len() - 1])
        },
    },
    ByteMutation {
        name: "alelo volteado",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| {
            let alleles = cell(text, 0, 5);
            let mut bytes = alleles.into_bytes();
            bytes[0] = if bytes[0] == b'0' { b'1' } else { b'0' };
            set_cell(text, 0, 5, &String::from_utf8(bytes).expect("UTF-8"))
        },
    },
    ByteMutation {
        name: "HAP_ ajeno",
        path: HAPLOTYPES_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| damage_identifier(text, 0, 2),
    },
    // --- marginales ---
    ByteMutation {
        name: "frecuencia realizada movida",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: |text| set_cell(text, 0, 5, "0.125"),
    },
    ByteMutation {
        name: "frecuencia realizada movida con su error",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: move_realized_marginal,
    },
    ByteMutation {
        name: "objetivo marginal movido",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: |text| set_cell(text, 0, 4, "0.375"),
    },
    ByteMutation {
        name: "objetivo marginal movido con su error y bajo el máximo",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: move_target_marginal,
    },
    ByteMutation {
        name: "NA fuera de un bin",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: |text| set_cell(text, 0, 4, "NA"),
    },
    ByteMutation {
        name: "sitios marginales desordenados",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: |text| swap_rows(text, 0, 1),
    },
    ByteMutation {
        name: "último sitio de una unidad marginal ausente",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Marginal,
        apply: drop_last_site_of_first_unit,
    },
    // --- signo y RNG ---
    ByteMutation {
        name: "polaridad en mayúsculas",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Polarity,
        apply: |text| set_cell(text, 0, 3, "POSITIVE"),
    },
    ByteMutation {
        name: "polaridad volteada",
        path: FREQUENCIES_PATH,
        cause: F5OracleCause::Polarity,
        apply: flip_polarity,
    },
    // --- D y r² ---
    ByteMutation {
        name: "media r² movida",
        path: LD_PATH,
        cause: F5OracleCause::LinkageDisequilibrium,
        apply: |text| set_cell(text, 0, 6, "0.5"),
    },
    ByteMutation {
        name: "media r² movida con su error y dentro de tolerancia",
        path: LD_PATH,
        cause: F5OracleCause::LinkageDisequilibrium,
        apply: move_mean_r2,
    },
    ByteMutation {
        name: "parejas del bin alteradas",
        path: LD_PATH,
        cause: F5OracleCause::LinkageDisequilibrium,
        apply: |text| bump_number_after(text, &format!("\t{}\t", cell(text, 0, 4))),
    },
    ByteMutation {
        name: "estado no aplicable con parejas",
        path: LD_PATH,
        cause: F5OracleCause::LinkageDisequilibrium,
        apply: |text| set_cell(text, 0, 4, "not_applicable"),
    },
    ByteMutation {
        name: "bins desordenados",
        path: LD_PATH,
        cause: F5OracleCause::LinkageDisequilibrium,
        apply: |text| swap_rows(text, 0, 1),
    },
    // --- intentos ---
    ByteMutation {
        name: "sin mejor intento",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::Attempt,
        apply: |text| replace_once(text, "\ttrue\n", "\tfalse\n"),
    },
    ByteMutation {
        name: "dos mejores intentos",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::Attempt,
        apply: |text| replace_once(text, "\tfalse\n", "\ttrue\n"),
    },
    ByteMutation {
        name: "selected_best no booleano",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::Attempt,
        apply: |text| replace_once(text, "\ttrue\n", "\t1\n"),
    },
    ByteMutation {
        name: "más aceptadas que propuestas",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::Attempt,
        apply: |text| {
            let proposals: u64 = cell(text, 0, 3).parse().expect("propuestas");
            set_cell(text, 0, 4, &(proposals + 1).to_string())
        },
    },
    ByteMutation {
        name: "intento fuera del límite",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::ResourceBound,
        apply: |text| set_cell(text, 0, 2, "4096"),
    },
    ByteMutation {
        name: "signos fallidos alterados",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::SignRule,
        apply: |text| {
            let row = text
                .lines()
                .skip(1)
                .position(|line| line.ends_with("\ttrue"))
                .expect("mejor intento");
            let failed: u64 = cell(text, row, 7).parse().expect("signos");
            set_cell(text, row, 7, &(failed + 1).to_string())
        },
    },
    ByteMutation {
        name: "filas repetidas alteradas",
        path: ATTEMPTS_PATH,
        cause: F5OracleCause::Attempt,
        apply: |text| {
            let row = text
                .lines()
                .skip(1)
                .position(|line| line.ends_with("\ttrue"))
                .expect("mejor intento");
            let duplicates: u64 = cell(text, row, 5).parse().expect("duplicados");
            set_cell(text, row, 5, &(duplicates + 1).to_string())
        },
    },
    // --- asignación fundadora ---
    ByteMutation {
        name: "fundadores desordenados",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Assignment,
        apply: |text| swap_rows(text, 0, 1),
    },
    ByteMutation {
        name: "índice de copia cambiado",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| set_cell(text, 0, 3, "1"),
    },
    ByteMutation {
        name: "FCP_ ajeno",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| damage_identifier(text, 0, 5),
    },
    ByteMutation {
        name: "FHA_ ajeno",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Identity,
        apply: |text| damage_identifier(text, 0, 6),
    },
    ByteMutation {
        name: "haplotipo inexistente asignado",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Assignment,
        apply: |text| damage_identifier(text, 0, 7),
    },
    ByteMutation {
        name: "haplotipo de otra unidad asignado",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Assignment,
        apply: |text| {
            let other = cell(text, 1, 7);
            assert_ne!(other, cell(text, 0, 7), "fixture degenerado");
            set_cell(text, 0, 7, &other)
        },
    },
    ByteMutation {
        name: "otra fila haplotípica de la misma unidad",
        path: FOUNDERS_PATH,
        cause: F5OracleCause::Identity,
        apply: swap_haplotype_within_unit,
    },
    // --- proyección de verdad ---
    ByteMutation {
        name: "truth publicado como v1",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| replace_once(text, "\"schema_version\":\"2\"", "\"schema_version\":\"1\""),
    },
    ByteMutation {
        name: "versión de proyección ajena",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            replace_once(
                text,
                "f5-founder-truth-projection-v2",
                "f5-founder-truth-projection-v1",
            )
        },
    },
    ByteMutation {
        name: "clase de proyección ajena",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| replace_once(text, "\"kind\":\"f5_founder\"", "\"kind\":\"f5_diploid\""),
    },
    ByteMutation {
        name: "censo de la cabecera alterado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| bump_number_after(text, "\"individual_count\":"),
    },
    ByteMutation {
        name: "recuento de copias alterado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| bump_number_after(text, "\"counts\":{\"copies\":"),
    },
    ByteMutation {
        name: "recuento de genotipos alterado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| bump_number_after(text, "\"genotypes\":"),
    },
    ByteMutation {
        name: "huella de la referencia alterada",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| replace_after(text, "\"reference_content_sha256\":\"", "a", "z"),
    },
    ByteMutation {
        name: "pertenencia de una copia cambiada",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            replace_after(
                text,
                "\"founder_membership\":{",
                "\"population_id\":\"PA\"",
                "\"population_id\":\"PB\"",
            )
        },
    },
    ByteMutation {
        name: "fase distinta del índice de copia",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| bump_number_after(text, "\"phase\":"),
    },
    ByteMutation {
        name: "copias desordenadas",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
            lines.swap(1, 2);
            format!("{}\n", lines.join("\n"))
        },
    },
    ByteMutation {
        name: "alelo alternativo de una copia cambiado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| rotate_base_after(text, "\"alternate\":\""),
    },
    ByteMutation {
        name: "sitio del locus cambiado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            replace_after(
                text,
                "\"f5_locus\":{",
                "\"site_id\":\"SITE_",
                "\"site_id\":\"SITX_",
            )
        },
    },
    ByteMutation {
        name: "población del locus cambiada",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            replace_after(
                text,
                "\"f5_locus\":{",
                "\"population_id\":\"PA\"",
                "\"population_id\":\"PB\"",
            )
        },
    },
    ByteMutation {
        name: "portador con otro alelo",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| flip_digit_after(text, "\"carriers\":[{\"allele_index\":"),
    },
    ByteMutation {
        name: "portador repetido con sus dosis ajustadas",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: duplicate_one_carrier,
    },
    ByteMutation {
        name: "case_id ajeno",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| replace_after(text, "\"case_id\":\"CASE_", "a", "b"),
    },
    ByteMutation {
        name: "procedencia incoherente entre las tres vistas",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: |text| {
            let anchor = "\"matrix_source_id\":\"HAP_";
            let at = text.find(anchor).expect("procedencia") + anchor.len();
            let mut bytes = text.as_bytes().to_vec();
            bytes[at] = if bytes[at] == b'A' { b'B' } else { b'A' };
            String::from_utf8(bytes).expect("UTF-8")
        },
    },
    ByteMutation {
        name: "una copia entera traza otra fila haplotípica",
        path: TRUTH_PATH,
        cause: F5OracleCause::Assignment,
        apply: retrace_first_copy,
    },
    ByteMutation {
        name: "un SNP de copia con otro alelo en sus dos vistas",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: relabel_one_copy_snp,
    },
    ByteMutation {
        name: "un SNP de copia fuera de su secuencia",
        path: TRUTH_PATH,
        cause: F5OracleCause::Truth,
        apply: push_last_modification_out_of_the_sequence,
    },
    ByteMutation {
        name: "dosis alterada",
        path: TRUTH_PATH,
        cause: F5OracleCause::Man70,
        apply: |text| bump_number_after(text, "\"allele_dosages\":["),
    },
    ByteMutation {
        name: "un genotipo menos con su recuento ajustado",
        path: TRUTH_PATH,
        cause: F5OracleCause::Man70,
        apply: drop_last_genotype,
    },
];

#[test]
fn every_isolated_byte_mutation_dies_of_its_own_cause() {
    let published = publish("byte-mutations", &wide_fixture());
    let mut detected = Vec::new();
    for mutation in BYTE_MUTATIONS {
        let mut staged = Staged::new(&published, "byte-mutation");
        staged.audit().expect("el paquete íntegro es válido");
        let original = staged.text(mutation.path);
        let corrupted = (mutation.apply)(&original);
        assert_ne!(
            corrupted, original,
            "{}: la mutación no cambió ningún byte",
            mutation.name
        );
        // El enlace se rehace con los bytes mutados: así la corrupción se juzga
        // por su semántica y no por la huella, que la delataría siempre.
        staged.rewrite(mutation.path, &corrupted);
        let cause = staged
            .audit()
            .err()
            .unwrap_or_else(|| panic!("«{}» pasó inadvertida", mutation.name));
        assert_eq!(
            cause, mutation.cause,
            "«{}» murió por {cause:?} y le corresponde {:?}",
            mutation.name, mutation.cause
        );
        detected.push(mutation.name);
    }
    assert_eq!(detected.len(), BYTE_MUTATIONS.len());
}

/// El enlace es la otra mitad: los mismos bytes con una huella o un tamaño
/// declarados que no son los suyos.
#[test]
fn a_link_that_does_not_describe_its_bytes_is_refused() {
    let published = publish("links", &wide_fixture());
    for path in [LAYOUT_PATH, FOUNDERS_PATH] {
        let mut staged = Staged::new(&published, "link-mutation");
        let text = staged.text(path);
        fs::write(staged.root.join(path), text).expect("reescritura");
        let link = staged
            .declaration
            .artifacts
            .iter_mut()
            .find(|link| link.path == path)
            .expect("enlace");
        link.size_bytes += 1;
        assert_eq!(staged.audit(), Err(F5OracleCause::Link), "{path}");
    }
    let mut staged = Staged::new(&published, "link-digest");
    staged.declaration.truth.sha256 = sha256(b"otra cosa");
    assert_eq!(staged.audit(), Err(F5OracleCause::Link));
    let mut staged = Staged::new(&published, "link-id");
    staged.declaration.artifacts[0].artifact_id = "f5-01-block-layout-v2".to_owned();
    assert_eq!(staged.audit(), Err(F5OracleCause::Link));
    // Un byte de más en un artefacto, con el enlace intacto, también muere.
    let staged = Staged::new(&published, "link-byte");
    let text = staged.text(SITES_PATH);
    fs::write(staged.root.join(SITES_PATH), format!("{text}\n")).expect("reescritura");
    assert_eq!(staged.audit(), Err(F5OracleCause::Link));
}

#[test]
fn a_missing_artifact_is_refused() {
    let published = publish("missing", &wide_fixture());
    for path in [LAYOUT_PATH, SITES_PATH, FOUNDERS_PATH, TRUTH_PATH] {
        let staged = Staged::new(&published, "missing-artifact");
        fs::remove_file(staged.root.join(path)).expect("artefacto retirado");
        assert_eq!(
            staged.audit(),
            Err(F5OracleCause::ArtifactMissing),
            "{path}"
        );
    }
}

// ---------------------------------------------------------------------------
// Mutación aislada de la declaración
// ---------------------------------------------------------------------------

type Perturbation = fn(&mut F5PackageDeclaration);

struct DeclarationMutation {
    name: &'static str,
    cause: F5OracleCause,
    apply: Perturbation,
}

const DECLARATION_MUTATIONS: &[DeclarationMutation] = &[
    // --- semilla y catálogos ---
    DeclarationMutation {
        name: "semilla truncada",
        cause: F5OracleCause::Seed,
        apply: |declaration| declaration.seed.truncate(63),
    },
    DeclarationMutation {
        name: "semilla en mayúsculas",
        cause: F5OracleCause::Seed,
        apply: |declaration| declaration.seed = declaration.seed.to_uppercase(),
    },
    DeclarationMutation {
        name: "otra semilla raíz",
        cause: F5OracleCause::Identity,
        apply: |declaration| declaration.seed = "0".repeat(64),
    },
    DeclarationMutation {
        name: "catálogo de identidades ajeno",
        cause: F5OracleCause::Randomness,
        apply: |declaration| declaration.id_catalog_id = "id-catalog-v3".to_owned(),
    },
    DeclarationMutation {
        name: "catálogo RNG ajeno",
        cause: F5OracleCause::Randomness,
        apply: |declaration| declaration.rng_catalog_id = "f5-rng-catalog-v2".to_owned(),
    },
    // --- identidad de la ejecución ---
    DeclarationMutation {
        name: "otra ejecución",
        cause: F5OracleCause::Identity,
        apply: |declaration| {
            declaration.execution_id = "RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4".to_owned();
        },
    },
    DeclarationMutation {
        name: "otra huella de referencia",
        cause: F5OracleCause::Truth,
        apply: |declaration| declaration.reference_sha256 = "0".repeat(64),
    },
    DeclarationMutation {
        name: "otro perfil de matriz",
        cause: F5OracleCause::Truth,
        apply: |declaration| declaration.profile_id = "otro-perfil".to_owned(),
    },
    DeclarationMutation {
        name: "subpoblaciones desordenadas",
        cause: F5OracleCause::Assignment,
        apply: |declaration| declaration.subpopulations.reverse(),
    },
    DeclarationMutation {
        name: "ploidía cero",
        cause: F5OracleCause::Assignment,
        apply: |declaration| declaration.ploidy = 0,
    },
    // --- DTC-049 ---
    DeclarationMutation {
        name: "auditorías DTC-049 desordenadas",
        cause: F5OracleCause::Normalization,
        apply: |declaration| declaration.normalizations.swap(0, 1),
    },
    DeclarationMutation {
        name: "una auditoría DTC-049 ausente",
        cause: F5OracleCause::Normalization,
        apply: |declaration| {
            declaration.normalizations.pop();
        },
    },
    DeclarationMutation {
        name: "pesos efectivos ajenos a su auditoría",
        cause: F5OracleCause::Normalization,
        apply: |declaration| {
            declaration.haplotype_weights = vec![0.75, 0.25];
        },
    },
    // --- modo, cuota y tolerancia ---
    DeclarationMutation {
        name: "modo efectivo degradado",
        cause: F5OracleCause::Sampling,
        apply: |declaration| {
            let row = declaration
                .sampling
                .iter_mut()
                .find(|row| row.effective_mode == F5DeclaredMode::Controlled)
                .expect("unidad controlada");
            row.requested_mode = F5DeclaredMode::Probabilistic;
        },
    },
    DeclarationMutation {
        name: "tolerancia en una unidad probabilística",
        cause: F5OracleCause::Sampling,
        apply: |declaration| {
            let row = declaration
                .sampling
                .iter_mut()
                .find(|row| row.effective_mode == F5DeclaredMode::Probabilistic)
                .expect("unidad probabilística");
            row.controlled_max_abs_error_effective = Some(0.1);
        },
    },
    DeclarationMutation {
        name: "tolerancia controlada por encima de 0.5/C",
        cause: F5OracleCause::Sampling,
        apply: |declaration| {
            let row = declaration
                .sampling
                .iter_mut()
                .find(|row| row.effective_mode == F5DeclaredMode::Controlled)
                .expect("unidad controlada");
            row.controlled_max_abs_error_declared = Some(0.75);
            row.controlled_max_abs_error_effective = Some(0.75);
        },
    },
    DeclarationMutation {
        name: "tolerancia declarada distinta de la efectiva",
        cause: F5OracleCause::Sampling,
        apply: |declaration| {
            let row = declaration
                .sampling
                .iter_mut()
                .find(|row| row.effective_mode == F5DeclaredMode::Controlled)
                .expect("unidad controlada");
            row.controlled_max_abs_error_declared = Some(0.0);
        },
    },
    DeclarationMutation {
        name: "una unidad sin fila de muestreo",
        cause: F5OracleCause::Sampling,
        apply: |declaration| {
            declaration.sampling.pop();
        },
    },
    // --- F ---
    DeclarationMutation {
        name: "F realizado alterado",
        cause: F5OracleCause::Inbreeding,
        apply: |declaration| {
            let row = declaration
                .f_realized
                .iter_mut()
                .find(|row| row.value.is_some())
                .expect("F realizado evaluado");
            row.value = row.value.map(|value| value + 0.25);
        },
    },
    DeclarationMutation {
        name: "F realizado declarado no aplicable",
        cause: F5OracleCause::Inbreeding,
        apply: |declaration| {
            let row = declaration
                .f_realized
                .iter_mut()
                .find(|row| row.value.is_some())
                .expect("F realizado evaluado");
            row.status = "not_applicable".to_owned();
            row.value = None;
        },
    },
    DeclarationMutation {
        name: "objetivo de F ajeno a su subpoblación",
        cause: F5OracleCause::Inbreeding,
        apply: |declaration| {
            declaration.f_realized[0].target_f = Some(0.5);
        },
    },
    DeclarationMutation {
        name: "F sin declarar en una subpoblación diploide",
        cause: F5OracleCause::Inbreeding,
        apply: |declaration| declaration.subpopulations[0].inbreeding_f = None,
    },
    DeclarationMutation {
        name: "F realizado desordenado",
        cause: F5OracleCause::Inbreeding,
        apply: |declaration| declaration.f_realized.reverse(),
    },
    // --- cotas de recurso ---
    DeclarationMutation {
        name: "bloques sellados de más",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.block_count += 1,
    },
    DeclarationMutation {
        name: "cota de sitios por bloque reducida",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.site_bound_per_block = 1,
    },
    DeclarationMutation {
        name: "cota de filas TSV reducida",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.tsv_rows_bound = 1,
    },
    DeclarationMutation {
        name: "cota de registros JSON reducida",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.json_rows_bound = 1,
    },
    DeclarationMutation {
        name: "presupuesto de bytes reducido",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.output_bytes = 1,
    },
    DeclarationMutation {
        name: "presupuesto por encima de su límite",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| {
            declaration.bounds.output_bytes = declaration.bounds.output_limit_bytes + 1;
        },
    },
    DeclarationMutation {
        name: "límite de intentos reducido",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.max_haplotype_attempts = 1,
    },
    DeclarationMutation {
        name: "límite de propuestas reducido",
        cause: F5OracleCause::ResourceBound,
        apply: |declaration| declaration.bounds.max_swap_proposals = 0,
    },
];

#[test]
fn every_isolated_declaration_mutation_dies_of_its_own_cause() {
    let published = publish("declaration-mutations", &wide_fixture());
    for mutation in DECLARATION_MUTATIONS {
        let mut staged = Staged::new(&published, "declaration-mutation");
        staged.audit().expect("el paquete íntegro es válido");
        let original = staged.declaration.clone();
        (mutation.apply)(&mut staged.declaration);
        assert_ne!(
            staged.declaration, original,
            "{}: la mutación no cambió la declaración",
            mutation.name
        );
        let cause = staged
            .audit()
            .err()
            .unwrap_or_else(|| panic!("«{}» pasó inadvertida", mutation.name));
        assert_eq!(
            cause, mutation.cause,
            "«{}» murió por {cause:?} y le corresponde {:?}",
            mutation.name, mutation.cause
        );
    }
}

// ---------------------------------------------------------------------------
// Ámbitos de configuración
// ---------------------------------------------------------------------------

/// Publica y audita un escenario entero, devolviendo sus resúmenes MAN-70.
fn accepted(label: &str, fixture: &Fixture<'_>) -> Vec<(F5GenotypeModel, Vec<u64>, Vec<f64>, u8)> {
    let published = publish(label, fixture);
    let staged = Staged::new(&published, label);
    let report = audit_f5_package(&staged.root, &staged.declaration)
        .unwrap_or_else(|rejection| panic!("{label}: {}", rejection.diagnostic().render_human()));
    report
        .initialization_summaries()
        .iter()
        .map(|summary| {
            (
                summary.model(),
                summary.recounted_genotype_counts().to_vec(),
                summary.expected_genotype_probabilities().to_vec(),
                summary.ploidy(),
            )
        })
        .collect()
}

/// Ploidías 1, 2 y 3: el número de dosis por locus es siempre `P+1` y el
/// recuento cuadra con el censo.
#[test]
fn the_oracle_accepts_every_covered_ploidy() {
    for ploidy in [1_u64, 2, 3] {
        let fixture = Fixture {
            ploidy,
            subpopulations: vec![Subpopulation::new("P1", 2)],
            ..Fixture::default()
        };
        let summaries = accepted(&format!("ploidy-{ploidy}"), &fixture);
        assert!(!summaries.is_empty(), "ploidía {ploidy}");
        for (model, counts, probabilities, published_ploidy) in &summaries {
            assert_eq!(u64::from(*published_ploidy), ploidy);
            assert_eq!(counts.len(), ploidy as usize + 1);
            assert_eq!(probabilities.len(), ploidy as usize + 1);
            assert_eq!(counts.iter().sum::<u64>(), 2);
            assert_eq!(*model, F5GenotypeModel::Hwe);
        }
    }
}

/// `F=0` es HWE y `F=1` es el modelo diploide con F: el oráculo distingue el
/// modelo por la declaración, no por el resultado.
#[test]
fn the_two_extremes_of_f_choose_their_own_model() {
    let hwe = accepted(
        "f-zero",
        &Fixture {
            scenario_inbreeding: Some("0.0"),
            subpopulations: vec![Subpopulation::new("P1", 2)],
            ..Fixture::default()
        },
    );
    assert!(hwe.iter().all(|(model, ..)| *model == F5GenotypeModel::Hwe));
    let inbred = accepted(
        "f-one",
        &Fixture {
            scenario_inbreeding: Some("1.0"),
            subpopulations: vec![Subpopulation::new("P1", 2)],
            ..Fixture::default()
        },
    );
    assert!(
        inbred
            .iter()
            .all(|(model, ..)| *model == F5GenotypeModel::DiploidF)
    );
    // Con `F=1` las dos copias son idénticas por descendencia: ningún locus
    // publica heterocigotos y el modelo no les da probabilidad.
    for (_, counts, probabilities, _) in &inbred {
        assert_eq!(counts[1], 0, "{counts:?}");
        assert_eq!(probabilities[1], 0.0, "{probabilities:?}");
    }
}

/// Fundadores aportados por identificador y derivados de clave publican el
/// mismo contrato.
#[test]
fn provided_and_derived_founder_identifiers_are_both_accepted() {
    for (label, subpopulation) in [
        (
            "provided",
            Subpopulation::new("P1", 2).founders(&["zeta", "alfa"]),
        ),
        ("derived", Subpopulation::new("P1", 2).keys()),
    ] {
        let summaries = accepted(
            label,
            &Fixture {
                subpopulations: vec![subpopulation],
                ..Fixture::default()
            },
        );
        assert!(!summaries.is_empty(), "{label}");
    }
}

/// Bins `evaluated` y `not_applicable` en la misma ejecución.
#[test]
fn both_ld_states_pass_the_oracle() {
    let published = publish("sparse", &sparse_fixture());
    let staged = Staged::new(&published, "sparse-stage");
    audit_f5_package(&staged.root, &staged.declaration).expect("paquete con bins no aplicables");
    let text = staged.text("outputs/f5-ld-evaluation.tsv");
    assert!(text.contains("\tevaluated\t"), "{text}");
    assert!(text.contains("\tnot_applicable\t"), "{text}");
    assert!(text.contains("\tNA\t"), "{text}");
}

/// Cuatro pesos distintos: la columna `q` no es degenerada y el oráculo la
/// contrasta contra el perfil efectivo.
#[test]
fn a_dyadic_profile_publishes_four_distinct_weights() {
    let published = publish("dyadic", &dyadic_fixture());
    let staged = Staged::new(&published, "dyadic-stage");
    audit_f5_package(&staged.root, &staged.declaration).expect("perfil diádico");
    let text = staged.text("outputs/f5-haplotypes.tsv");
    for weight in ["0.5", "0.25", "0.125"] {
        assert!(text.contains(&format!("\t{weight}\t")), "{text}");
    }
}

/// Una cuota controlada que su tolerancia efectiva no sostiene se rechaza por
/// muestreo, no por cualquier otra causa.
#[test]
fn a_controlled_quota_outside_its_tolerance_is_refused() {
    let fixture = Fixture {
        subpopulations: vec![Subpopulation::new("P1", 1).sampling("controlled")],
        weights: vec!["0.5", "0.25", "0.125", "0.125"],
        block_length: 3,
        bins: vec![(1, 2, 0.25, 0.05)],
        max_haplotype_attempts: 512,
        max_swap_proposals: 2048,
        ..Fixture::default()
    };
    let published = publish("quota", &fixture);
    let mut staged = Staged::new(&published, "quota-stage");
    audit_f5_package(&staged.root, &staged.declaration).expect("cuota dentro de su tolerancia");
    // `C=2` y la tolerancia efectiva es `0.5/C=0.25`: con `0.1` la cuota
    // realizada del segundo rango deja de caber.
    for row in &mut staged.declaration.sampling {
        assert_eq!(row.effective_mode, F5DeclaredMode::Controlled);
        row.controlled_max_abs_error_declared = Some(0.1);
        row.controlled_max_abs_error_effective = Some(0.1);
    }
    assert_eq!(staged.audit(), Err(F5OracleCause::Sampling));
}

/// Una secuencia estéril no aporta ninguna copia ni ningún genotipo: el
/// oráculo no espera una copia por secuencia retenida.
#[test]
fn a_barren_sequence_contributes_nothing() {
    let published = publish("barren", &wide_fixture());
    let staged = Staged::new(&published, "barren-stage");
    let report = audit_f5_package(&staged.root, &staged.declaration).expect("paquete válido");
    let founders = staged.text("outputs/f5-founders.tsv");
    assert!(!founders.contains("barren1"), "{founders}");
    // Dos subpoblaciones, dos secuencias pobladas, ploidía 2 y censos 2 y 1.
    assert_eq!(report.copies(), (2 + 1) * 2 * 2);
}

// ---------------------------------------------------------------------------
// Pico real del oráculo, medido con asignador contador
// ---------------------------------------------------------------------------
//
// Gemelo del arnés de `tests/f5_artifacts.rs`, que vigila lo que el escritor
// asigna. Aquí se vigila lo que el **validador** asigna: un oráculo que
// materializase el paquete entero para juzgarlo tendría el mismo defecto que
// F5.4.1 tuvo al materializar los registros JSON.

// ---------------------------------------------------------------------------
// Pico real del escritor, medido con asignador contador
// ---------------------------------------------------------------------------
//
// Gemelo de `the_retained_repository_does_not_grow_with_the_block_count` de
// `tests/f5_founders.rs`, que vigila lo que F5.3 **retiene**. Aquí se vigila
// lo que F5.4.1 **asigna mientras escribe**, que la relación de memoria de
// DTC-089 no contabilizaba: `truth_bound` acota bytes de salida y
// `final_formula_peak` no reserva ni un byte por copia serializada.
//
// El asignador solo cuenta el hilo que mide, de modo que las demás pruebas
// del binario —que corren en paralelo— no contaminan la cifra. Los tres
// contadores son `Cell` de inicialización constante y sin `Drop`: leerlos
// desde `alloc`/`dealloc` no asigna ni registra destructores.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;

thread_local! {
    static TRACKING: Cell<bool> = const { Cell::new(false) };
    static LIVE: Cell<i64> = const { Cell::new(0) };
    static PEAK: Cell<i64> = const { Cell::new(0) };
}

fn note(delta: i64) {
    if TRACKING.try_with(Cell::get).unwrap_or(false) {
        let live = LIVE.with(|cell| {
            let value = cell.get() + delta;
            cell.set(value);
            value
        });
        PEAK.with(|cell| {
            if live > cell.get() {
                cell.set(live);
            }
        });
    }
}

struct CountingAllocator;

// SAFETY: cada método delega en `System`, que satisface el contrato de
// `GlobalAlloc`, y solo añade un contador de enteros que no toca el puntero
// devuelto ni la memoria asignada.
unsafe impl GlobalAlloc for CountingAllocator {
    unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc(layout) };
        if !pointer.is_null() {
            note(layout.size() as i64);
        }
        pointer
    }
    unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
        let pointer = unsafe { System.alloc_zeroed(layout) };
        if !pointer.is_null() {
            note(layout.size() as i64);
        }
        pointer
    }
    unsafe fn dealloc(&self, pointer: *mut u8, layout: Layout) {
        note(-(layout.size() as i64));
        unsafe { System.dealloc(pointer, layout) };
    }
    unsafe fn realloc(&self, pointer: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
        let next = unsafe { System.realloc(pointer, layout, new_size) };
        if !next.is_null() {
            note(new_size as i64 - layout.size() as i64);
        }
        next
    }
}

#[global_allocator]
static ALLOCATOR: CountingAllocator = CountingAllocator;

/// Ejecuta `body` y devuelve su resultado con el pico de bytes vivos que
/// asignó **este** hilo mientras duraba.
fn peak_bytes<T>(body: impl FnOnce() -> T) -> (T, u64) {
    LIVE.with(|cell| cell.set(0));
    PEAK.with(|cell| cell.set(0));
    TRACKING.with(|cell| cell.set(true));
    let value = body();
    TRACKING.with(|cell| cell.set(false));
    let peak = PEAK.with(Cell::get);
    (value, u64::try_from(peak).unwrap_or_default())
}

/// Holgura por byte del registro que el lector tiene en vuelo.
///
/// La razón exceso/registro converge: vale 4,56 con `B=8`, 7,67 con `B=32` y
/// 8,48 con `B=128`, que es ya la asíntota del árbol `serde_json::Value` con
/// `arbitrary_precision`. Diez deja margen sin dejar pasar una regresión que
/// la llevara a doce.
const ORACLE_RECORD_ALLOWANCE: u64 = 10;

/// El pico del oráculo cabe en el presupuesto sellado **más** un múltiplo del
/// registro que tiene en vuelo, y ese múltiplo no crece.
///
/// El oráculo retiene lo mismo que el productor —geometría, matrices
/// ganadoras en bits, un rango por asiento y un vector de dosis por locus— más
/// el índice de sitios por copia que el recuento MAN-70 necesita para resolver
/// portadores sin abrir ningún TSV. Todo eso son términos que
/// `f5-resource-accounting-v8` ya concede a las mismas dimensiones.
///
/// Lo que **no** cabe en `working_memory_bytes` es el registro
/// `haplotype_copy` que el lector parsea: el escritor lo emite por flujo y la
/// fórmula lo contabiliza en `output_bytes` a través de `truth_bound`, de modo
/// que un lector que deba interpretarlo necesita materializarlo una vez.
///
/// **El exceso crece linealmente con `B`**, porque una copia cubre los bloques
/// de su secuencia; lo que lo salva es que el registro más largo crece también
/// linealmente, de modo que la razón entre ambos converge. No es el defecto de
/// F5.4.1 —allí el pico crecía con el **número** de registros, no con su
/// tamaño—, pero decir «acotado por el registro más largo, no por el paquete»
/// sin más sería engañoso: lo acotado es la razón, y por eso la prueba la
/// mide en tres tamaños y la afirma.
#[test]
fn the_oracle_peak_keeps_a_bounded_ratio_to_the_record_it_holds() {
    for blocks in [8_u64, 32, 128] {
        let fixture = Fixture {
            block_count: blocks,
            slack: 2,
            ..Fixture::default()
        };
        let published = publish(&format!("oracle-peak-{blocks}"), &fixture);
        let sealed_peak = published.chain.sealed.working_memory_bytes();
        let staged = Staged::new(&published, &format!("oracle-peak-stage-{blocks}"));
        let longest = staged
            .text("outputs/truth.jsonl")
            .lines()
            .map(str::len)
            .max()
            .expect("registro") as u64;
        let (report, peak) = peak_bytes(|| audit_f5_package(&staged.root, &staged.declaration));
        report.expect("paquete válido");
        let excess = peak.saturating_sub(sealed_peak);
        assert!(
            excess <= ORACLE_RECORD_ALLOWANCE * longest,
            "B={blocks}: el oráculo excede el sellado en {excess} bytes con un registro de \
             {longest}, razón {:.2} por encima de {ORACLE_RECORD_ALLOWANCE}",
            excess as f64 / longest as f64
        );
    }

    // Multiplicar por dieciséis los genotipos no puede multiplicar el pico: lo
    // que crece es un rango por asiento y un índice de sitios por copia, no el
    // paquete. La misma razón vale, y además el pico no llega ni a duplicarse.
    let mut census_peaks = Vec::new();
    for census in [1_u64, 16] {
        let fixture = Fixture {
            block_count: 8,
            slack: 2,
            subpopulations: vec![Subpopulation::new("P1", census)],
            ..Fixture::default()
        };
        let published = publish(&format!("oracle-census-{census}"), &fixture);
        let sealed_peak = published.chain.sealed.working_memory_bytes();
        let staged = Staged::new(&published, &format!("oracle-census-stage-{census}"));
        let longest = staged
            .text("outputs/truth.jsonl")
            .lines()
            .map(str::len)
            .max()
            .expect("registro") as u64;
        let (report, peak) = peak_bytes(|| audit_f5_package(&staged.root, &staged.declaration));
        let genotypes = report.expect("paquete válido").genotypes();
        assert!(
            peak.saturating_sub(sealed_peak) <= ORACLE_RECORD_ALLOWANCE * longest,
            "censo={census}: el oráculo pica {peak} bytes sobre {sealed_peak} sellados"
        );
        census_peaks.push((genotypes, peak));
    }
    let (small_genotypes, small_peak) = census_peaks[0];
    let (large_genotypes, large_peak) = census_peaks[1];
    assert_eq!(large_genotypes, small_genotypes * 16);
    assert!(
        large_peak < 2 * small_peak,
        "dieciséis veces los genotipos ({small_genotypes} a {large_genotypes}) \
         llevaron el pico de {small_peak} a {large_peak}"
    );
}

// ---------------------------------------------------------------------------
// Sustitución coherente del alelo alternativo
// ---------------------------------------------------------------------------

/// `case_id` de DTC-061, reimplementado aquí para poder rehacerlo tras mover
/// el alelo alternativo de un locus.
fn oracle_case_id(
    population_id: &str,
    individual_id: &str,
    sequence_id: &str,
    position_1based: u64,
    alternate: char,
) -> String {
    let mut hasher = Sha256::new();
    for value in [
        population_id,
        individual_id,
        sequence_id,
        &position_1based.to_string(),
        &format!("REF,SNP:{alternate}"),
    ] {
        hasher.update(value.as_bytes());
        hasher.update([0]);
    }
    format!("CASE_{}", &format!("{:x}", hasher.finalize())[..59])
}

/// Sustituye el ALT de un sitio en **todas** sus vistas: la columna de
/// `f5-sites.tsv`, las modificaciones y el fragmento vigente de cada copia que
/// lo lleva, el `f5_locus`, el alelo tipado de cada genotipo y su `case_id`.
///
/// Es la reproducción literal del hallazgo bloqueante de la revisión: antes de
/// rederivar el ALT desde `site_alt`, el paquete resultante era **aceptado**,
/// y la verdad fundadora afirmaba que los individuos portan una base que el
/// contrato no sortea.
fn substitute_alternate_base(staged: &mut Staged) -> (String, String) {
    let sites = staged.text(SITES_PATH);
    let site_id = cell(&sites, 0, 2);
    let reference = cell(&sites, 0, 5);
    let alternate = cell(&sites, 0, 6);
    let replacement = other_base(&reference, &alternate);
    staged.rewrite(SITES_PATH, &set_cell(&sites, 0, 6, &replacement));

    let new_base = replacement.chars().next().expect("base");
    let truth = staged.text(TRUTH_PATH);
    let mut lines = Vec::with_capacity(truth.lines().count());
    for line in truth.lines() {
        if line.contains("\"record_type\":\"haplotype_copy\"") {
            // Las dos vistas que publican el alelo lo escriben pegado a su
            // `event_id`; el historial no lo publica.
            lines.push(line.replace(
                &format!("\"alternate\":\"{alternate}\",\"event_id\":\"{site_id}\""),
                &format!("\"alternate\":\"{replacement}\",\"event_id\":\"{site_id}\""),
            ));
            continue;
        }
        if !line.contains(&format!("\"site_id\":\"{site_id}\"")) {
            lines.push(line.to_owned());
            continue;
        }
        let mut record: serde_json::Value = serde_json::from_str(line).expect("genotipo");
        let object = record.as_object_mut().expect("objeto");
        let population_id = object["f5_locus"]["population_id"]
            .as_str()
            .expect("población")
            .to_owned();
        let individual_id = object["individual_id"]
            .as_str()
            .expect("individuo")
            .to_owned();
        let sequence_id = object["locus"]["sequence_id"]
            .as_str()
            .expect("secuencia")
            .to_owned();
        let start = object["locus"]["start"].as_u64().expect("inicio");
        object["f5_locus"]["alt"] = serde_json::Value::String(replacement.clone());
        for allele in object["alleles"].as_array_mut().expect("alelos") {
            if allele["kind"] == "snp" {
                allele["base"] = serde_json::Value::String(replacement.clone());
            }
        }
        object["case_id"] = serde_json::Value::String(oracle_case_id(
            &population_id,
            &individual_id,
            &sequence_id,
            start + 1,
            new_base,
        ));
        lines.push(serde_json::to_string(&record).expect("serialización"));
    }
    staged.rewrite(TRUTH_PATH, &format!("{}\n", lines.join("\n")));
    (alternate, replacement)
}

/// El paquete coherente debe seguir siendo válido tras un viaje de ida y
/// vuelta por `serde_json`: si no lo fuera, la mutación siguiente estaría
/// midiendo el serializador y no el contrato.
#[test]
fn a_genotype_record_round_trips_through_serde_json_byte_for_byte() {
    let published = publish("roundtrip", &wide_fixture());
    let staged = Staged::new(&published, "roundtrip-stage");
    let truth = staged.text(TRUTH_PATH);
    let mut seen = 0_usize;
    for line in truth
        .lines()
        .filter(|line| line.contains("\"record_type\":\"genotype\""))
    {
        let record: serde_json::Value = serde_json::from_str(line).expect("genotipo");
        assert_eq!(serde_json::to_string(&record).expect("serialización"), line);
        seen += 1;
    }
    assert!(seen > 0);
}

/// Un alelo alternativo sustituido **de forma coherente en todas sus vistas**
/// muere por el sitio: es el sorteo `site_alt` quien lo fija, y el oráculo lo
/// rederiva desde `scenario_id`, `block_id` y `site_id`.
#[test]
fn a_coherently_substituted_alternate_base_is_refused() {
    let published = publish("alt-substitution", &wide_fixture());
    let mut staged = Staged::new(&published, "alt-substitution-stage");
    staged.audit().expect("el paquete íntegro es válido");
    let (original, replacement) = substitute_alternate_base(&mut staged);
    assert_ne!(original, replacement);
    // El paquete mutado es internamente coherente: sus ocho ficheros y sus dos
    // enlaces se sostienen entre sí, y solo la rederivación del sorteo lo
    // delata.
    assert_eq!(staged.audit(), Err(F5OracleCause::Site));
}

// ---------------------------------------------------------------------------
// Guardas que solo son alcanzables por su propia frontera
// ---------------------------------------------------------------------------

/// La cota `position < sequence_length` de una copia no la alcanza
/// `audit_f5_package` —el cruce con `f5-sites.tsv` mata antes cualquier
/// posición movida—, pero sí la frontera pública del recuento, que no abre
/// ningún TSV. Ahí es donde tiene que estar y donde se prueba.
#[test]
fn a_copy_snp_outside_its_sequence_is_refused_by_the_recount() {
    let published = publish("outside", &wide_fixture());
    let directory = TempDir::new("outside-alone");
    let truth = directory.path().join("truth.jsonl");
    fs::write(&truth, published.bytes("outputs/truth.jsonl")).expect("copia");
    recount_f5_man70(&truth).expect("la verdad íntegra se reconta");
    let moved =
        push_last_modification_out_of_the_sequence(&fs::read_to_string(&truth).expect("verdad"));
    fs::write(&truth, moved).expect("verdad movida");
    let rejection = recount_f5_man70(&truth).expect_err("un SNP fuera de su secuencia");
    assert_eq!(rejection.cause(), F5OracleCause::Truth);
    assert!(
        rejection
            .diagnostic()
            .render_human()
            .contains("fuera de su secuencia"),
        "{}",
        rejection.diagnostic().render_human()
    );
}

/// Un bin publicado cuyo error excede su tolerancia se rechaza.
///
/// Es el criterio de aceptación de DTC-089 y solo es alcanzable moviendo a la
/// vez la tolerancia declarada y la columna publicada: cualquier mutación que
/// suba el error muere antes contra la recomputación de `r²`.
#[test]
fn a_bin_whose_error_exceeds_its_tolerance_is_refused() {
    let published = publish("tolerance", &dyadic_fixture());
    let mut staged = Staged::new(&published, "tolerance-stage");
    staged.audit().expect("el paquete íntegro es válido");
    let text = staged.text(LD_PATH);
    let row = text
        .lines()
        .skip(1)
        .position(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            fields[4] == "evaluated" && fields[9].parse::<f64>().is_ok_and(|error| error > 0.0)
        })
        .expect("un bin evaluado con error positivo");
    let error: f64 = cell(&text, row, 9).parse().expect("error");
    let tightened = error / 2.0;
    staged.rewrite(LD_PATH, &set_cell(&text, row, 8, &tightened.to_string()));
    // La declaración acompaña a la columna: el perfil pasa a exigir esa misma
    // tolerancia, de modo que la fila es coherente y lo único que falla es el
    // criterio de aceptación.
    let start: u64 = cell(&text, row, 2).parse().expect("inicio");
    for bin in &mut staged.declaration.ld_bins {
        if bin.distance_start_bp == start {
            bin.max_abs_error = tightened;
        }
    }
    assert_eq!(staged.audit(), Err(F5OracleCause::LinkageDisequilibrium));
}

/// El catálogo RNG publicado tiene las doce decisiones que el oráculo exige.
///
/// `validate_rng_catalog` cerraba con `seen == 0`, de modo que un catálogo
/// recortado a una fila dejaba al oráculo autovalidándose contra una sola
/// decisión mientras su reparto 3+9 —tres rederivadas, nueve declaradas
/// creídas— seguía escrito sobre doce. La constante y el fichero se contrastan
/// aquí: si alguien añade o quita una fila sin tocar la cifra, esto muere.
#[test]
fn the_published_rng_catalogue_has_the_twelve_decisions_the_oracle_demands() {
    const CATALOG: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");
    let rows = CATALOG
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#') && !line.starts_with("decision"))
        .count();
    assert_eq!(rows, F5_RNG_CATALOG_DECISIONS);
    assert_eq!(F5_RNG_CATALOG_DECISIONS, 12);
}
