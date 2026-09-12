//! F5.4.1. Los siete TSV de éxito y la proyección fundadora `truth.jsonl` v2.
//!
//! La batería no comprueba la escritura contra el escritor: reconstruye el
//! paquete publicado con un lector propio que solo abre bytes. Ese lector
//! rederiva con DTC-037/DTC-039 —BLAKE3 y Base32 reimplementados aquí y
//! validados antes contra los vectores congelados del catálogo— los cinco
//! identificadores automáticos, cruza los siete ficheros entre sí y, por
//! separado, reconstruye identidad, pertenencia, fase, alelos, dosis y
//! recuentos por locus×subpoblación **leyendo solo `truth.jsonl` v2**.
//!
//! El oráculo de orden de bloque es el propio identificador: `BLK_` lleva
//! `block_index_u64` en su material y ningún fichero publica esa columna, de
//! modo que el lector la deriva del ordinal de fila de `f5-block-layout.tsv`.
//! Si «bloque» no fuera el orden canónico del layout, la rederivación fallaría.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::f5_artifacts::{
    F5_TSV_CONTRACTS, F5ReferenceSequences, F5SuccessArtifacts, F5SuccessRequest,
    write_f5_success_artifacts,
};
use genome_population_simulator::f5_founders::{
    F5FounderAssignment, F5FounderIdentities, assign_f5_founders, materialize_f5_founder_identities,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
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
        "gps-f5-artifacts-{}-{}-{label}",
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
    const fn inbreeding(mut self, value: &'a str) -> Self {
        self.inbreeding_f = Some(value);
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
    fn text(&self, relative: &str) -> String {
        String::from_utf8(self.bytes(relative)).expect("UTF-8")
    }
    fn package(&self) -> Package {
        let mut files = BTreeMap::new();
        for contract in &F5_TSV_CONTRACTS {
            files.insert(contract.path().to_owned(), self.bytes(contract.path()));
        }
        files.insert(
            self.artifacts.truth().path().to_owned(),
            self.bytes(self.artifacts.truth().path()),
        );
        Package {
            files,
            seed: self.chain.config.seed().to_owned(),
        }
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

/// Publica sustituyendo las longitudes de secuencia declaradas.
fn publish_with_sequences(
    label: &str,
    fixture: &Fixture<'_>,
    sequences: F5ReferenceSequences,
) -> Diagnostic {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("root");
    let mut chain = build_chain(&root, fixture);
    chain.sequences = sequences;
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
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
        .expect_err("el escritor debe rechazar")
}

// ---------------------------------------------------------------------------
// Oráculo independiente de identidades (DTC-037/DTC-039)
// ---------------------------------------------------------------------------

const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
const ID_DERIVATION_VERSION: u32 = 1;

fn push_component(material: &mut Vec<u8>, name: &str, value: &str) {
    material.extend_from_slice(&u32::try_from(name.len()).expect("nombre").to_be_bytes());
    material.extend_from_slice(name.as_bytes());
    material.extend_from_slice(&u32::try_from(value.len()).expect("valor").to_be_bytes());
    material.extend_from_slice(value.as_bytes());
}

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

fn root_seed_bytes(seed: &str) -> [u8; 32] {
    decode_hex(seed).try_into().expect("semilla de 256 bits")
}

/// Contexto común de las cinco clases F5.
struct IdOracle {
    root: [u8; 32],
    execution_id: String,
    scenario_id: String,
    reference_id: String,
}

impl IdOracle {
    fn head(&self) -> Vec<(&'static str, String)> {
        vec![
            ("execution_id", self.execution_id.clone()),
            ("scenario_id", self.scenario_id.clone()),
            ("reference_id", self.reference_id.clone()),
        ]
    }

    fn block(&self, block_index: u64) -> String {
        let mut components = self.head();
        components.push(("block_index_u64", block_index.to_string()));
        derived_identifier(&self.root, "BLK_", "f5_block", &components)
    }

    fn site(&self, block_index: u64, site_index: u64) -> String {
        let mut components = self.head();
        components.push(("block_index_u64", block_index.to_string()));
        components.push(("site_index_u64", site_index.to_string()));
        derived_identifier(&self.root, "SITE_", "f5_site", &components)
    }

    fn haplotype(&self, population: &str, block_index: u64, rank: u64, alleles: &str) -> String {
        let mut components = self.head();
        components.push(("population_id", population.to_owned()));
        components.push(("block_index_u64", block_index.to_string()));
        components.push(("haplotype_rank_u64", rank.to_string()));
        components.push(("alleles", alleles.to_owned()));
        derived_identifier(&self.root, "HAP_", "f5_haplotype", &components)
    }

    fn founder_copy(
        &self,
        population: &str,
        individual: &str,
        sequence: &str,
        copy_index: u32,
    ) -> String {
        let mut components = self.head();
        components.push(("population_id", population.to_owned()));
        components.push(("individual_id", individual.to_owned()));
        components.push(("sequence_id", sequence.to_owned()));
        components.push(("copy_index_u32", copy_index.to_string()));
        derived_identifier(&self.root, "FCP_", "f5_founder_copy", &components)
    }

    fn assignment(
        &self,
        population: &str,
        block_index: u64,
        founder_copy_id: &str,
        haplotype_id: &str,
    ) -> String {
        let mut components = self.head();
        components.push(("population_id", population.to_owned()));
        components.push(("block_index_u64", block_index.to_string()));
        components.push(("founder_copy_id", founder_copy_id.to_owned()));
        components.push(("haplotype_id", haplotype_id.to_owned()));
        derived_identifier(
            &self.root,
            "FHA_",
            "f5_block_haplotype_assignment",
            &components,
        )
    }
}

// ---------------------------------------------------------------------------
// Lector estricto de los siete TSV
// ---------------------------------------------------------------------------

struct Package {
    files: BTreeMap<String, Vec<u8>>,
    seed: String,
}

impl Package {
    fn file(&self, path: &str) -> Result<&[u8], String> {
        self.files
            .get(path)
            .map(Vec::as_slice)
            .ok_or_else(|| format!("falta el artefacto {path}"))
    }
}

/// Aplica el contrato de bytes común a los siete TSV y devuelve sus filas.
fn read_tsv(package: &Package, index: usize) -> Result<Vec<Vec<String>>, String> {
    let contract = &F5_TSV_CONTRACTS[index];
    let bytes = package.file(contract.path())?;
    if bytes.starts_with(&[0xEF, 0xBB, 0xBF]) {
        return Err(format!("{}: BOM prohibido", contract.path()));
    }
    let text =
        std::str::from_utf8(bytes).map_err(|_| format!("{}: no es UTF-8", contract.path()))?;
    if !text.ends_with('\n') {
        return Err(format!("{}: falta el LF final", contract.path()));
    }
    if text.contains('\r') {
        return Err(format!("{}: CR prohibido", contract.path()));
    }
    if text.contains('"') {
        return Err(format!("{}: comillas prohibidas", contract.path()));
    }
    let mut lines = text.strip_suffix('\n').unwrap_or(text).split('\n');
    let header = lines
        .next()
        .ok_or_else(|| format!("{}: sin cabecera", contract.path()))?;
    if header != contract.header() {
        return Err(format!(
            "{}: cabecera inesperada {header:?}",
            contract.path()
        ));
    }
    let columns = contract.header().split('\t').count();
    let mut rows = Vec::new();
    for line in lines {
        if line.is_empty() {
            return Err(format!("{}: línea vacía prohibida", contract.path()));
        }
        if line.starts_with('#') {
            return Err(format!("{}: comentario prohibido", contract.path()));
        }
        let row: Vec<String> = line.split('\t').map(str::to_owned).collect();
        if row.len() != columns {
            return Err(format!(
                "{}: la fila tiene {} columnas y la cabecera {columns}",
                contract.path(),
                row.len()
            ));
        }
        rows.push(row);
    }
    Ok(rows)
}

fn canonical_u64(path: &str, value: &str) -> Result<u64, String> {
    if value.is_empty()
        || (value.len() > 1 && value.starts_with('0'))
        || !value.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(format!(
            "{path}: {value:?} no es un entero decimal canónico"
        ));
    }
    value
        .parse()
        .map_err(|_| format!("{path}: {value:?} no cabe en u64"))
}

/// `f64` con el literal mínimo de ida y vuelta, como exige el contrato.
fn canonical_f64(path: &str, value: &str) -> Result<f64, String> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| format!("{path}: {value:?} no es un f64"))?;
    if !parsed.is_finite() || parsed.to_string() != value {
        return Err(format!(
            "{path}: {value:?} no es el decimal mínimo de ida y vuelta"
        ));
    }
    Ok(parsed)
}

fn probability(path: &str, value: &str) -> Result<f64, String> {
    let parsed = canonical_f64(path, value)?;
    if !(0.0..=1.0).contains(&parsed) {
        return Err(format!("{path}: {value:?} no pertenece a [0,1]"));
    }
    Ok(parsed)
}

// ---------------------------------------------------------------------------
// Lector independiente de `truth.jsonl` v2
// ---------------------------------------------------------------------------

use serde_json::{Map, Value, json};

fn exact_object<'a>(
    value: &'a Value,
    what: &str,
    keys: &[&str],
) -> Result<&'a Map<String, Value>, String> {
    let object = value
        .as_object()
        .ok_or_else(|| format!("{what}: no es un objeto"))?;
    let observed: BTreeSet<&str> = object.keys().map(String::as_str).collect();
    let expected: BTreeSet<&str> = keys.iter().copied().collect();
    if observed != expected {
        return Err(format!(
            "{what}: claves {observed:?} en lugar de {expected:?}"
        ));
    }
    Ok(object)
}

fn text(object: &Map<String, Value>, what: &str, key: &str) -> Result<String, String> {
    object
        .get(key)
        .and_then(Value::as_str)
        .map(str::to_owned)
        .ok_or_else(|| format!("{what}: {key} no es una cadena"))
}

fn number(object: &Map<String, Value>, what: &str, key: &str) -> Result<u64, String> {
    object
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| format!("{what}: {key} no es un entero no negativo"))
}

/// Cabecera v2 ya validada.
struct TruthHeader {
    execution_id: String,
    reference_id: String,
    scenario_id: String,
    /// Huella de contenido de la referencia materializada, que DTC-061
    /// mantiene separada de los digests de artefacto de `source_digests`.
    reference_content_sha256: String,
    /// `source_digests` en el orden publicado: clase, identificador y huella.
    source_digests: Vec<(String, String, String)>,
    populations: Vec<(String, u64, u64)>,
    copies: u64,
    genotypes: u64,
}

fn read_truth_header(value: &Value) -> Result<TruthHeader, String> {
    let object = exact_object(
        value,
        "truth.header",
        &[
            "record_type",
            "document_kind",
            "schema_version",
            "execution_id",
            "reference_id",
            "reference_content_sha256",
            "coordinate_convention",
            "counts",
            "source_digests",
            "group_digests",
            "projection",
            "populations",
        ],
    )?;
    if text(object, "truth.header", "record_type")? != "header" {
        return Err("truth.header: record_type no es header".to_owned());
    }
    if text(object, "truth.header", "document_kind")? != "gps_truth" {
        return Err("truth.header: document_kind no es gps_truth".to_owned());
    }
    let version = text(object, "truth.header", "schema_version")?;
    if version != "2" {
        return Err(format!(
            "truth.header: schema_version {version:?}; una proyección F5 exige la 2 y publicar v1 es error de versión"
        ));
    }
    if text(object, "truth.header", "coordinate_convention")? != "zero_based_half_open_v1" {
        return Err("truth.header: convenio de coordenadas ajeno".to_owned());
    }
    let digest = text(object, "truth.header", "reference_content_sha256")?;
    if digest.len() != 64
        || !digest
            .bytes()
            .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
    {
        return Err("truth.header: reference_content_sha256 no es hex minúscula de 64".to_owned());
    }
    let projection = exact_object(
        object
            .get("projection")
            .ok_or("truth.header: sin projection")?,
        "truth.projection",
        &["kind", "projection_version", "scenario_id"],
    )?;
    if text(projection, "truth.projection", "kind")? != "f5_founder" {
        return Err("truth.projection: kind no es f5_founder".to_owned());
    }
    if text(projection, "truth.projection", "projection_version")?
        != "f5-founder-truth-projection-v2"
    {
        return Err("truth.projection: projection_version ajena".to_owned());
    }
    let counts = exact_object(
        object.get("counts").ok_or("truth.header: sin counts")?,
        "truth.counts",
        &["copies", "genotypes"],
    )?;
    let mut populations = Vec::new();
    let listed = object
        .get("populations")
        .and_then(Value::as_array)
        .ok_or("truth.header: populations no es una lista")?;
    if listed.is_empty() {
        return Err("truth.header: populations vacía".to_owned());
    }
    for entry in listed {
        let population = exact_object(
            entry,
            "truth.populations",
            &["population_id", "individual_count", "ploidy"],
        )?;
        populations.push((
            text(population, "truth.populations", "population_id")?,
            number(population, "truth.populations", "individual_count")?,
            number(population, "truth.populations", "ploidy")?,
        ));
    }
    if populations.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
        return Err("truth.header: populations no está ordenada por population_id".to_owned());
    }
    let listed_digests = object
        .get("source_digests")
        .and_then(Value::as_array)
        .ok_or("truth.header: source_digests no es una lista")?;
    if listed_digests.is_empty() {
        return Err("truth.header: source_digests vacía".to_owned());
    }
    let mut source_digests = Vec::new();
    for entry in listed_digests {
        let entry = exact_object(entry, "truth.source_digests", &["source", "sha256"])?;
        let source = exact_object(
            entry
                .get("source")
                .ok_or("truth.source_digests: sin source")?,
            "truth.source_digests.source",
            &["kind", "source_id"],
        )?;
        let value = text(entry, "truth.source_digests", "sha256")?;
        if value.len() != 64
            || !value
                .bytes()
                .all(|byte| byte.is_ascii_hexdigit() && !byte.is_ascii_uppercase())
        {
            return Err("truth.source_digests: sha256 no es hex minúscula de 64".to_owned());
        }
        source_digests.push((
            text(source, "truth.source_digests.source", "kind")?,
            text(source, "truth.source_digests.source", "source_id")?,
            value,
        ));
    }
    Ok(TruthHeader {
        execution_id: text(object, "truth.header", "execution_id")?,
        reference_id: text(object, "truth.header", "reference_id")?,
        scenario_id: text(projection, "truth.projection", "scenario_id")?,
        reference_content_sha256: digest,
        source_digests,
        populations,
        copies: number(counts, "truth.counts", "copies")?,
        genotypes: number(counts, "truth.counts", "genotypes")?,
    })
}

/// Identidad de copia reconstruida.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CopyKey {
    reference_id: String,
    individual_id: String,
    sequence_id: String,
    copy_index: u64,
}

fn read_copy_identity(value: &Value, what: &str) -> Result<CopyKey, String> {
    let object = exact_object(
        value,
        what,
        &["reference_id", "individual_id", "sequence_id", "copy_index"],
    )?;
    Ok(CopyKey {
        reference_id: text(object, what, "reference_id")?,
        individual_id: text(object, what, "individual_id")?,
        sequence_id: text(object, what, "sequence_id")?,
        copy_index: number(object, what, "copy_index")?,
    })
}

/// Todo lo que el contrato exige recuperar leyendo **solo** `truth.jsonl` v2.
struct TruthWorld {
    header: TruthHeader,
    /// Pertenencia y `FCP_` de cada copia.
    membership: BTreeMap<CopyKey, (String, String)>,
    /// Fase publicada de cada copia.
    phase: BTreeMap<CopyKey, u64>,
    /// Sitios donde la copia lleva el alelo alternativo, con su base.
    alternate_sites: BTreeMap<CopyKey, BTreeMap<String, (u64, char, char)>>,
    /// Procedencia de matriz de cada SNP de cada copia: `matrix_source_id` y
    /// `population_id`, ya comprobados iguales en las tres vistas.
    snp_provenance: BTreeMap<CopyKey, BTreeMap<String, (String, String)>>,
    /// `case_id` de cada genotipo y la tupla que lo determina.
    cases: BTreeMap<String, (String, String, String, u64, char)>,
    /// Dosis alternativa por locus×subpoblación e individuo.
    dosage: BTreeMap<(String, String), BTreeMap<String, u64>>,
    /// Recuento de individuos por dosis y locus×subpoblación.
    counts: BTreeMap<(String, String), BTreeMap<u64, u64>>,
    /// Geometría del locus vista desde el genotipo: bloque, secuencia,
    /// posición 1-based, `ref` y `alt`.
    locus: BTreeMap<String, (String, String, u64, char, char)>,
    /// Orden de aparición de los grupos población×bloque en los genotipos.
    block_order: Vec<(String, String)>,
    copies: u64,
    genotypes: u64,
}

/// Reconstruye el mundo fundador sin abrir ningún TSV ni llamar al productor.
#[allow(clippy::too_many_lines)]
fn reconstruct_from_truth(bytes: &[u8]) -> Result<TruthWorld, String> {
    let text_body = std::str::from_utf8(bytes).map_err(|_| "truth: no es UTF-8".to_owned())?;
    if !text_body.ends_with('\n') {
        return Err("truth: falta el LF final".to_owned());
    }
    let mut records = text_body
        .strip_suffix('\n')
        .unwrap_or(text_body)
        .split('\n');
    let header_line = records.next().ok_or("truth: documento vacío")?;
    let header_value: Value =
        serde_json::from_str(header_line).map_err(|error| format!("truth.header: {error}"))?;
    let header = read_truth_header(&header_value)?;
    let declared: BTreeMap<String, (u64, u64)> = header
        .populations
        .iter()
        .map(|(id, individuals, ploidy)| (id.clone(), (*individuals, *ploidy)))
        .collect();

    let mut world = TruthWorld {
        membership: BTreeMap::new(),
        phase: BTreeMap::new(),
        alternate_sites: BTreeMap::new(),
        snp_provenance: BTreeMap::new(),
        cases: BTreeMap::new(),
        dosage: BTreeMap::new(),
        counts: BTreeMap::new(),
        locus: BTreeMap::new(),
        block_order: Vec::new(),
        copies: 0,
        genotypes: 0,
        header,
    };
    let mut previous_copy: Option<(String, String, String, String, u64)> = None;
    let mut current_group: Option<(String, String)> = None;
    let mut previous_locus: Option<(String, String)> = None;
    let mut seen_genotype = false;
    for line in records {
        if line.is_empty() {
            return Err("truth: línea vacía prohibida".to_owned());
        }
        let value: Value = serde_json::from_str(line)
            .map_err(|error| format!("truth: registro ilegible: {error}"))?;
        let kind = value
            .get("record_type")
            .and_then(Value::as_str)
            .ok_or("truth: registro sin record_type")?;
        match kind {
            "haplotype_copy" => {
                if seen_genotype {
                    return Err(
                        "truth: una copia aparece después de un genotipo; las copias van primero"
                            .to_owned(),
                    );
                }
                let object = exact_object(
                    &value,
                    "truth.copy",
                    &[
                        "record_type",
                        "identity",
                        "phase",
                        "provenance",
                        "segments",
                        "current_pieces",
                        "history",
                        "founder_membership",
                    ],
                )?;
                let identity = read_copy_identity(
                    object.get("identity").ok_or("truth.copy: sin identity")?,
                    "truth.copy.identity",
                )?;
                let phase = number(object, "truth.copy", "phase")?;
                if phase != identity.copy_index {
                    return Err("truth.copy: la fase no es el índice de copia".to_owned());
                }
                let membership = exact_object(
                    object
                        .get("founder_membership")
                        .ok_or("truth.copy: sin founder_membership")?,
                    "truth.copy.founder_membership",
                    &["population_id", "founder_copy_id"],
                )?;
                let population =
                    text(membership, "truth.copy.founder_membership", "population_id")?;
                if !declared.contains_key(&population) {
                    return Err(format!(
                        "truth.copy: la pertenencia nombra la subpoblación {population}, ausente de la cabecera"
                    ));
                }
                let founder_copy_id = text(
                    membership,
                    "truth.copy.founder_membership",
                    "founder_copy_id",
                )?;
                let key = (
                    identity.reference_id.clone(),
                    population.clone(),
                    identity.individual_id.clone(),
                    identity.sequence_id.clone(),
                    phase,
                );
                if previous_copy.as_ref().is_some_and(|last| last >= &key) {
                    return Err(
                        "truth: las copias no llegan por referencia, población, individuo, secuencia y fase"
                            .to_owned(),
                    );
                }
                previous_copy = Some(key);
                let (sites, provenance) =
                    read_copy_modifications(object, &identity, &population, &world.header)?;
                if world
                    .membership
                    .insert(identity.clone(), (population, founder_copy_id))
                    .is_some()
                {
                    return Err("truth: una copia aparece dos veces".to_owned());
                }
                world.phase.insert(identity.clone(), phase);
                world.snp_provenance.insert(identity.clone(), provenance);
                world.alternate_sites.insert(identity, sites);
                world.copies += 1;
            }
            "genotype" => {
                seen_genotype = true;
                let object = exact_object(
                    &value,
                    "truth.genotype",
                    &[
                        "record_type",
                        "case_id",
                        "reference_id",
                        "individual_id",
                        "case_kind",
                        "locus",
                        "alleles",
                        "allele_dosages",
                        "carriers",
                        "f5_locus",
                    ],
                )?;
                let placed = read_genotype(&mut world, object, &declared)?;
                let group = (placed.population_id.clone(), placed.block_id.clone());
                if current_group.as_ref() != Some(&group) {
                    if world.block_order.contains(&group) {
                        return Err(format!(
                            "truth: el grupo {}×{} vuelve a aparecer; los genotipos deben agruparse por bloque",
                            group.0, group.1
                        ));
                    }
                    if world
                        .block_order
                        .last()
                        .is_some_and(|last| last.0 > group.0)
                    {
                        return Err(
                            "truth: los genotipos no llegan ordenados por población".to_owned()
                        );
                    }
                    world.block_order.push(group.clone());
                    current_group = Some(group);
                    previous_locus = None;
                }
                // «Sitio» en la clave de los genotipos es `site_id`: es el
                // único campo de sitio del registro y la misma lectura que
                // `f5-frequencies.tsv`, cuya única columna de sitio también lo
                // es. La posición física no ordena nada.
                let key = (placed.site_id.clone(), placed.individual_id.clone());
                if previous_locus.as_ref().is_some_and(|last| last >= &key) {
                    return Err(
                        "truth: dentro de un bloque los genotipos no llegan por site_id e individuo"
                            .to_owned(),
                    );
                }
                previous_locus = Some(key);
                world.genotypes += 1;
            }
            other => return Err(format!("truth: record_type desconocido {other:?}")),
        }
    }
    if world.copies != world.header.copies || world.genotypes != world.header.genotypes {
        return Err(format!(
            "truth: la cabecera declara {}/{} y el documento trae {}/{}",
            world.header.copies, world.header.genotypes, world.copies, world.genotypes
        ));
    }
    for ((site_id, population), individuals) in &world.dosage {
        let expected = declared
            .get(population)
            .ok_or_else(|| format!("truth: la subpoblación {population} no está declarada"))?;
        if individuals.len() as u64 != expected.0 {
            return Err(format!(
                "truth: el locus {site_id}×{population} trae {} individuos y la cabecera declara {}",
                individuals.len(),
                expected.0
            ));
        }
    }
    Ok(world)
}

/// Comprueba que las tres vistas de una copia —modificaciones del segmento,
/// SNP vigentes del fragmento e historial— describen los mismos alelos, y
/// devuelve los sitios alternativos que la copia lleva.
#[allow(clippy::type_complexity)]
fn read_copy_modifications(
    object: &Map<String, Value>,
    identity: &CopyKey,
    population_id: &str,
    header: &TruthHeader,
) -> Result<
    (
        BTreeMap<String, (u64, char, char)>,
        BTreeMap<String, (String, String)>,
    ),
    String,
> {
    let provenance = exact_object(
        object
            .get("provenance")
            .ok_or("truth.copy: sin provenance")?,
        "truth.copy.provenance",
        &["origin", "sources"],
    )?;
    let sources = provenance
        .get("sources")
        .and_then(Value::as_array)
        .ok_or("truth.copy.provenance: sources no es una lista")?;
    if sources.is_empty() {
        return Err("truth.copy.provenance: sources vacía".to_owned());
    }
    // Las fuentes de la copia son exactamente las que la cabecera declara en
    // `source_digests`: perder una —o inventarla— rompe la traza de la que
    // cuelga toda la procedencia fina.
    let mut declared_sources = Vec::new();
    for source in sources {
        let source = exact_object(
            source,
            "truth.copy.provenance.source",
            &["kind", "source_id"],
        )?;
        declared_sources.push((
            text(source, "truth.copy.provenance.source", "kind")?,
            text(source, "truth.copy.provenance.source", "source_id")?,
        ));
    }
    let expected_sources: Vec<(String, String)> = header
        .source_digests
        .iter()
        .map(|(kind, source_id, _)| (kind.clone(), source_id.clone()))
        .collect();
    if declared_sources != expected_sources {
        return Err(format!(
            "truth.copy.provenance: las fuentes {declared_sources:?} no son las {expected_sources:?} de la cabecera"
        ));
    }
    let matrix_source_id = expected_sources
        .iter()
        .find(|(kind, _)| kind == "matrix")
        .map(|(_, source_id)| source_id.clone())
        .ok_or("truth.copy.provenance: ninguna fuente es la matriz")?;
    let segments = object
        .get("segments")
        .and_then(Value::as_array)
        .ok_or("truth.copy: segments no es una lista")?;
    if segments.len() != 1 {
        return Err("truth.copy: una copia fundadora publica un segmento por secuencia".to_owned());
    }
    let segment = exact_object(
        &segments[0],
        "truth.copy.segment",
        &["source", "modifications"],
    )?;
    let source = read_segment_source(
        segment
            .get("source")
            .ok_or("truth.copy.segment: sin source")?,
        identity,
    )?;
    let mut sites = BTreeMap::new();
    let mut provenances: BTreeMap<String, (String, String)> = BTreeMap::new();
    let mut ordered = Vec::new();
    for value in segment
        .get("modifications")
        .and_then(Value::as_array)
        .ok_or("truth.copy.segment: modifications no es una lista")?
    {
        let modification = exact_object(
            value,
            "truth.copy.modification",
            &[
                "kind",
                "event_id",
                "position",
                "reference",
                "alternate",
                "provenance",
            ],
        )?;
        if text(modification, "truth.copy.modification", "kind")? != "snp" {
            return Err("truth.copy: F5 solo publica modificaciones SNP".to_owned());
        }
        let event_id = text(modification, "truth.copy.modification", "event_id")?;
        let position = number(modification, "truth.copy.modification", "position")?;
        if position >= source {
            return Err("truth.copy: un SNP cae fuera de su secuencia".to_owned());
        }
        let reference = single_base(modification, "truth.copy.modification", "reference")?;
        let alternate = single_base(modification, "truth.copy.modification", "alternate")?;
        if reference == alternate {
            return Err("truth.copy: un SNP publica el mismo alelo dos veces".to_owned());
        }
        let matrix =
            read_snp_provenance(modification, &event_id, population_id, &matrix_source_id)?;
        if sites
            .insert(event_id.clone(), (position + 1, reference, alternate))
            .is_some()
        {
            return Err("truth.copy: un sitio aparece dos veces en la misma copia".to_owned());
        }
        provenances.insert(event_id.clone(), matrix);
        ordered.push((position, event_id));
    }
    if ordered.windows(2).any(|pair| pair[0].0 >= pair[1].0) {
        return Err("truth.copy: las modificaciones no llegan por posición ascendente".to_owned());
    }

    let pieces = object
        .get("current_pieces")
        .and_then(Value::as_array)
        .ok_or("truth.copy: current_pieces no es una lista")?;
    if pieces.len() != 1 {
        return Err("truth.copy: una copia fundadora publica un fragmento vigente".to_owned());
    }
    let piece = exact_object(&pieces[0], "truth.copy.piece", &["kind", "source", "snps"])?;
    if text(piece, "truth.copy.piece", "kind")? != "source" {
        return Err("truth.copy: el fragmento vigente no procede de la referencia".to_owned());
    }
    let piece_source = read_segment_source(
        piece.get("source").ok_or("truth.copy.piece: sin source")?,
        identity,
    )?;
    if piece_source != source {
        return Err("truth.copy: el fragmento vigente no cubre la misma secuencia".to_owned());
    }
    let mut current = BTreeMap::new();
    for value in piece
        .get("snps")
        .and_then(Value::as_array)
        .ok_or("truth.copy.piece: snps no es una lista")?
    {
        let snp = exact_object(
            value,
            "truth.copy.piece.snp",
            &[
                "event_id",
                "position",
                "reference",
                "alternate",
                "provenance",
            ],
        )?;
        let event_id = text(snp, "truth.copy.piece.snp", "event_id")?;
        if provenances.get(&event_id)
            != Some(&read_snp_provenance(
                snp,
                &event_id,
                population_id,
                &matrix_source_id,
            )?)
        {
            return Err(
                "truth.copy: el fragmento vigente no repite la procedencia de su modificación"
                    .to_owned(),
            );
        }
        current.insert(
            event_id,
            (
                number(snp, "truth.copy.piece.snp", "position")? + 1,
                single_base(snp, "truth.copy.piece.snp", "reference")?,
                single_base(snp, "truth.copy.piece.snp", "alternate")?,
            ),
        );
    }
    if current != sites {
        return Err(
            "truth.copy: el fragmento vigente y las modificaciones describen alelos distintos"
                .to_owned(),
        );
    }

    let mut history = BTreeMap::new();
    for value in object
        .get("history")
        .and_then(Value::as_array)
        .ok_or("truth.copy: history no es una lista")?
    {
        let event = exact_object(
            value,
            "truth.copy.history",
            &["kind", "event_id", "position", "provenance"],
        )?;
        if text(event, "truth.copy.history", "kind")? != "snp" {
            return Err("truth.copy: F5 solo publica hitos SNP".to_owned());
        }
        let event_id = text(event, "truth.copy.history", "event_id")?;
        if provenances.get(&event_id)
            != Some(&read_snp_provenance(
                event,
                &event_id,
                population_id,
                &matrix_source_id,
            )?)
        {
            return Err(
                "truth.copy: el historial no repite la procedencia de su modificación".to_owned(),
            );
        }
        history.insert(
            event_id,
            number(event, "truth.copy.history", "position")? + 1,
        );
    }
    let expected: BTreeMap<String, u64> = sites
        .iter()
        .map(|(site, (position, _, _))| (site.clone(), *position))
        .collect();
    if history != expected {
        return Err("truth.copy: el historial no coincide con las modificaciones".to_owned());
    }
    Ok((sites, provenances))
}

/// Valida la fuente del segmento y devuelve la longitud de la secuencia.
fn read_segment_source(value: &Value, identity: &CopyKey) -> Result<u64, String> {
    let source = exact_object(
        value,
        "truth.copy.segment.source",
        &[
            "source",
            "reference_id",
            "sequence_id",
            "start",
            "end",
            "sequence_length",
        ],
    )?;
    exact_object(
        source
            .get("source")
            .ok_or("truth.copy.segment.source: sin source")?,
        "truth.copy.segment.source.source",
        &["kind", "source_id"],
    )?;
    if text(source, "truth.copy.segment.source", "reference_id")? != identity.reference_id
        || text(source, "truth.copy.segment.source", "sequence_id")? != identity.sequence_id
    {
        return Err("truth.copy: el segmento no pertenece a la copia".to_owned());
    }
    let start = number(source, "truth.copy.segment.source", "start")?;
    let end = number(source, "truth.copy.segment.source", "end")?;
    let length = number(source, "truth.copy.segment.source", "sequence_length")?;
    if start != 0 || end != length || length == 0 {
        return Err("truth.copy: el segmento no cubre la secuencia completa".to_owned());
    }
    Ok(length)
}

/// Devuelve `(matrix_source_id, population_id)` de la procedencia de matriz.
fn read_snp_provenance(
    object: &Map<String, Value>,
    event_id: &str,
    population_id: &str,
    matrix_source_id: &str,
) -> Result<(String, String), String> {
    let provenance = exact_object(
        object
            .get("provenance")
            .ok_or("truth: SNP sin procedencia")?,
        "truth.snp.provenance",
        &["kind", "version", "source", "matrix"],
    )?;
    if text(provenance, "truth.snp.provenance", "kind")? != "matrix" {
        return Err("truth: un SNP fundador no procede de una matriz".to_owned());
    }
    if text(provenance, "truth.snp.provenance", "version")? != "1" {
        return Err("truth: versión de procedencia ajena".to_owned());
    }
    let source = exact_object(
        provenance
            .get("source")
            .ok_or("truth.snp.provenance: sin source")?,
        "truth.snp.provenance.source",
        &["kind", "source_id"],
    )?;
    if text(source, "truth.snp.provenance.source", "kind")? != "matrix"
        || text(source, "truth.snp.provenance.source", "source_id")? != matrix_source_id
    {
        return Err("truth: la fuente del SNP no es la matriz que la cabecera declara".to_owned());
    }
    let matrix = exact_object(
        provenance
            .get("matrix")
            .ok_or("truth.snp.provenance: sin matrix")?,
        "truth.snp.provenance.matrix",
        &["matrix_source_id", "population_id", "variant_id"],
    )?;
    if text(matrix, "truth.snp.provenance.matrix", "variant_id")? != event_id {
        return Err("truth: la procedencia de matriz no nombra el sitio del evento".to_owned());
    }
    let declared = text(matrix, "truth.snp.provenance.matrix", "population_id")?;
    if declared != population_id {
        return Err(format!(
            "truth: la procedencia del SNP {event_id} declara {declared} y su copia pertenece a {population_id}"
        ));
    }
    Ok((
        text(matrix, "truth.snp.provenance.matrix", "matrix_source_id")?,
        declared,
    ))
}

fn single_base(object: &Map<String, Value>, what: &str, key: &str) -> Result<char, String> {
    let value = text(object, what, key)?;
    let mut characters = value.chars();
    match (characters.next(), characters.next()) {
        (Some(base @ ('A' | 'C' | 'G' | 'T')), None) => Ok(base),
        _ => Err(format!("{what}: {key} no es una base canónica")),
    }
}

/// Genotipo ya validado, con lo que el orden y el cruce necesitan.
struct PlacedGenotype {
    population_id: String,
    block_id: String,
    individual_id: String,
    site_id: String,
}

fn read_genotype(
    world: &mut TruthWorld,
    object: &Map<String, Value>,
    declared: &BTreeMap<String, (u64, u64)>,
) -> Result<PlacedGenotype, String> {
    if text(object, "truth.genotype", "case_kind")? != "snp" {
        return Err("truth.genotype: F5 solo publica casos SNP".to_owned());
    }
    let reference_id = text(object, "truth.genotype", "reference_id")?;
    let individual_id = text(object, "truth.genotype", "individual_id")?;
    let locus = exact_object(
        object.get("locus").ok_or("truth.genotype: sin locus")?,
        "truth.genotype.locus",
        &["sequence_id", "start", "end"],
    )?;
    let start = number(locus, "truth.genotype.locus", "start")?;
    let end = number(locus, "truth.genotype.locus", "end")?;
    if end != start + 1 {
        return Err("truth.genotype: el locus SNP no mide una base".to_owned());
    }
    let sequence_id = text(locus, "truth.genotype.locus", "sequence_id")?;
    let f5_locus = exact_object(
        object
            .get("f5_locus")
            .ok_or("truth.genotype: sin f5_locus")?,
        "truth.genotype.f5_locus",
        &["population_id", "block_id", "site_id", "ref", "alt"],
    )?;
    let population_id = text(f5_locus, "truth.genotype.f5_locus", "population_id")?;
    let block_id = text(f5_locus, "truth.genotype.f5_locus", "block_id")?;
    let site_id = text(f5_locus, "truth.genotype.f5_locus", "site_id")?;
    let reference_base = single_base(f5_locus, "truth.genotype.f5_locus", "ref")?;
    let alternate_base = single_base(f5_locus, "truth.genotype.f5_locus", "alt")?;
    if reference_base == alternate_base {
        return Err("truth.genotype: f5_locus repite el alelo".to_owned());
    }
    let ploidy = declared
        .get(&population_id)
        .ok_or_else(|| format!("truth.genotype: subpoblación {population_id} no declarada"))?
        .1;
    let geometry = (
        block_id.clone(),
        sequence_id.clone(),
        start + 1,
        reference_base,
        alternate_base,
    );
    match world.locus.get(&site_id) {
        None => {
            world.locus.insert(site_id.clone(), geometry);
        }
        Some(previous) if *previous == geometry => {}
        Some(_) => {
            return Err(format!(
                "truth.genotype: el sitio {site_id} cambia de bloque, secuencia, posición o alelos"
            ));
        }
    }

    let alleles = object
        .get("alleles")
        .and_then(Value::as_array)
        .ok_or("truth.genotype: alleles no es una lista")?;
    let mut alternate_index = None;
    for (index, allele) in alleles.iter().enumerate() {
        let kind = allele
            .get("kind")
            .and_then(Value::as_str)
            .ok_or("truth.genotype: alelo sin kind")?;
        match kind {
            "reference" => {
                exact_object(allele, "truth.genotype.allele", &["kind"])?;
            }
            "snp" => {
                let snp = exact_object(allele, "truth.genotype.allele", &["kind", "base"])?;
                if single_base(snp, "truth.genotype.allele", "base")? != alternate_base {
                    return Err("truth.genotype: el alelo tipado no es el ALT del locus".to_owned());
                }
                if alternate_index.replace(index).is_some() {
                    return Err(
                        "truth.genotype: dos alelos alternativos en un sitio bialélico".to_owned(),
                    );
                }
            }
            other => return Err(format!("truth.genotype: alelo {other:?} ajeno a F5")),
        }
    }
    let alternate_index =
        alternate_index.ok_or("truth.genotype: el sitio no publica su alelo alternativo")?;

    let dosages: Vec<u64> = object
        .get("allele_dosages")
        .and_then(Value::as_array)
        .ok_or("truth.genotype: allele_dosages no es una lista")?
        .iter()
        .map(|value| value.as_u64().ok_or("dosis no entera".to_owned()))
        .collect::<Result<_, _>>()?;
    if dosages.len() != alleles.len() {
        return Err("truth.genotype: las dosis no acompañan a los alelos".to_owned());
    }
    if dosages.iter().sum::<u64>() != ploidy {
        return Err("truth.genotype: las dosis no suman la ploidía declarada".to_owned());
    }

    let mut from_carriers = vec![0_u64; alleles.len()];
    let mut carrier_copies = BTreeSet::new();
    for value in object
        .get("carriers")
        .and_then(Value::as_array)
        .ok_or("truth.genotype: carriers no es una lista")?
    {
        let carrier = exact_object(value, "truth.genotype.carrier", &["allele_index", "copy"])?;
        let index = usize::try_from(number(carrier, "truth.genotype.carrier", "allele_index")?)
            .map_err(|_| "truth.genotype: allele_index desbordado".to_owned())?;
        if index >= alleles.len() {
            return Err("truth.genotype: allele_index fuera de la lista de alelos".to_owned());
        }
        let copy = read_copy_identity(
            carrier
                .get("copy")
                .ok_or("truth.genotype.carrier: sin copy")?,
            "truth.genotype.carrier.copy",
        )?;
        if copy.reference_id != reference_id
            || copy.individual_id != individual_id
            || copy.sequence_id != sequence_id
        {
            return Err("truth.genotype: un portador no es una copia de su individuo".to_owned());
        }
        let membership = world
            .membership
            .get(&copy)
            .ok_or("truth.genotype: un portador no tiene copia publicada")?;
        if membership.0 != population_id {
            return Err(
                "truth.genotype: f5_locus no coincide con la pertenencia de un portador".to_owned(),
            );
        }
        let carries = world
            .alternate_sites
            .get(&copy)
            .is_some_and(|sites| sites.contains_key(&site_id));
        if carries != (index == alternate_index) {
            return Err(format!(
                "truth.genotype: el portador {} de {site_id} contradice los alelos de su copia",
                copy.copy_index
            ));
        }
        if !carrier_copies.insert(copy) {
            return Err("truth.genotype: un portador aparece dos veces".to_owned());
        }
        from_carriers[index] += 1;
    }
    if from_carriers != dosages {
        return Err("truth.genotype: las dosis no son las de sus portadores".to_owned());
    }
    if carrier_copies.len() as u64 != ploidy {
        return Err("truth.genotype: los portadores no son las copias de la ploidía".to_owned());
    }

    // `case_id` es la clave del caso: única por población, individuo,
    // secuencia y posición. Un literal constante colapsaría todas en una.
    let case = text(object, "truth.genotype", "case_id")?;
    if world
        .cases
        .insert(
            case.clone(),
            (
                population_id.clone(),
                individual_id.clone(),
                sequence_id.clone(),
                start + 1,
                alternate_base,
            ),
        )
        .is_some()
    {
        return Err(format!("truth.genotype: el case_id {case} se repite"));
    }

    let locus_key = (site_id, population_id.clone());
    let alternate_dosage = dosages[alternate_index];
    if world
        .dosage
        .entry(locus_key.clone())
        .or_default()
        .insert(individual_id.clone(), alternate_dosage)
        .is_some()
    {
        return Err("truth.genotype: un individuo repite su genotipo en un locus".to_owned());
    }
    *world
        .counts
        .entry(locus_key.clone())
        .or_default()
        .entry(alternate_dosage)
        .or_default() += 1;

    Ok(PlacedGenotype {
        population_id,
        block_id,
        individual_id,
        site_id: locus_key.0,
    })
}

// ---------------------------------------------------------------------------
// Cruce completo de los siete TSV con `truth.jsonl`
// ---------------------------------------------------------------------------

const BLOCK_LAYOUT: usize = 0;
const SITES: usize = 1;
const HAPLOTYPES: usize = 2;
const FREQUENCIES: usize = 3;
const LD_EVALUATION: usize = 4;
const FOUNDERS: usize = 5;
const ATTEMPTS: usize = 6;

struct BlockRow {
    block_id: String,
    sequence_id: String,
    start_1based: u64,
    end_1based: u64,
}

struct SiteRow {
    site_id: String,
    position_1based: u64,
    reference_base: char,
    alternate_base: char,
}

/// Aplica el contrato completo al paquete publicado. Devuelve el mundo
/// reconstruido desde `truth.jsonl` cuando todo cuadra.
#[allow(clippy::too_many_lines)]
fn check_package(package: &Package) -> Result<TruthWorld, String> {
    let truth = reconstruct_from_truth(package.file("outputs/truth.jsonl")?)?;
    let oracle = IdOracle {
        root: root_seed_bytes(&package.seed),
        execution_id: truth.header.execution_id.clone(),
        scenario_id: truth.header.scenario_id.clone(),
        reference_id: truth.header.reference_id.clone(),
    };
    let ploidy = truth.header.populations[0].2;
    let populations: Vec<&str> = truth
        .header
        .populations
        .iter()
        .map(|(id, _, _)| id.as_str())
        .collect();

    // --- f5-block-layout.tsv -------------------------------------------------
    let path = F5_TSV_CONTRACTS[BLOCK_LAYOUT].path();
    let rows = read_tsv(package, BLOCK_LAYOUT)?;
    if rows.is_empty() {
        return Err(format!("{path}: el layout tiene al menos un bloque"));
    }
    let mut blocks = Vec::new();
    let mut block_index_of = BTreeMap::new();
    let mut sequence_order: Vec<String> = Vec::new();
    let mut previous: Option<(usize, u64, u64)> = None;
    for (index, row) in rows.iter().enumerate() {
        let block_index = index as u64;
        if row[0] != oracle.block(block_index) {
            return Err(format!(
                "{path}: la fila {index} no rederiva su BLK_; «bloque» no es el orden canónico del layout"
            ));
        }
        let start = canonical_u64(path, &row[2])?;
        let end = canonical_u64(path, &row[3])?;
        let length = canonical_u64(path, &row[4])?;
        if start == 0 || start > end || length != end - start + 1 {
            return Err(format!("{path}: geometría inválida en la fila {index}"));
        }
        let ordinal = match sequence_order.iter().position(|id| *id == row[1]) {
            Some(position) => position,
            None => {
                sequence_order.push(row[1].clone());
                sequence_order.len() - 1
            }
        };
        let key = (ordinal, start, end);
        if previous.is_some_and(|last| last >= key) {
            return Err(format!(
                "{path}: el orden no es FASTA, inicio y fin en la fila {index}"
            ));
        }
        previous = Some(key);
        block_index_of.insert(row[0].clone(), block_index);
        blocks.push(BlockRow {
            block_id: row[0].clone(),
            sequence_id: row[1].clone(),
            start_1based: start,
            end_1based: end,
        });
    }

    // --- f5-sites.tsv --------------------------------------------------------
    let path = F5_TSV_CONTRACTS[SITES].path();
    let rows = read_tsv(package, SITES)?;
    let mut sites: Vec<Vec<SiteRow>> = (0..blocks.len()).map(|_| Vec::new()).collect();
    let mut previous_block: Option<String> = None;
    for row in &rows {
        let block_index = *block_index_of
            .get(&row[0])
            .ok_or_else(|| format!("{path}: {} no existe en el layout", row[0]))?;
        if previous_block.as_ref() != Some(&row[0]) {
            if previous_block.as_ref().is_some_and(|last| *last >= row[0]) {
                return Err(format!(
                    "{path}: los bloques no llegan por block_id ascendente"
                ));
            }
            previous_block = Some(row[0].clone());
        }
        let block = &blocks[block_index as usize];
        let site_index = canonical_u64(path, &row[1])?;
        let slot = &mut sites[block_index as usize];
        if site_index != slot.len() as u64 {
            return Err(format!("{path}: site_index no es consecutivo desde 0"));
        }
        if row[2] != oracle.site(block_index, site_index) {
            return Err(format!("{path}: el sitio no rederiva su SITE_"));
        }
        if row[3] != block.sequence_id {
            return Err(format!(
                "{path}: el sitio no nombra la secuencia de su bloque"
            ));
        }
        let position = canonical_u64(path, &row[4])?;
        if position < block.start_1based || position > block.end_1based {
            return Err(format!("{path}: el sitio cae fuera de su bloque"));
        }
        let reference_base = one_base(path, &row[5])?;
        let alternate_base = one_base(path, &row[6])?;
        if reference_base == alternate_base {
            return Err(format!("{path}: ref y alt coinciden"));
        }
        if slot
            .last()
            .is_some_and(|last| last.position_1based >= position)
        {
            return Err(format!(
                "{path}: las posiciones no crecen dentro del bloque"
            ));
        }
        slot.push(SiteRow {
            site_id: row[2].clone(),
            position_1based: position,
            reference_base,
            alternate_base,
        });
    }
    if sites.iter().any(Vec::is_empty) {
        return Err(format!("{path}: un bloque no publica ningún sitio"));
    }

    // --- f5-haplotypes.tsv ---------------------------------------------------
    let path = F5_TSV_CONTRACTS[HAPLOTYPES].path();
    let rows = read_tsv(package, HAPLOTYPES)?;
    let mut haplotypes: BTreeMap<(String, String), Vec<(String, String)>> = BTreeMap::new();
    let mut previous: Option<(String, String, u64)> = None;
    let mut weights: Option<Vec<String>> = None;
    for row in &rows {
        if !populations.contains(&row[0].as_str()) {
            return Err(format!(
                "{path}: {} no es una subpoblación publicada",
                row[0]
            ));
        }
        let block_index = *block_index_of
            .get(&row[1])
            .ok_or_else(|| format!("{path}: {} no existe en el layout", row[1]))?;
        let rank = canonical_u64(path, &row[3])?;
        let key = (row[0].clone(), row[1].clone(), rank);
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(format!("{path}: el orden no es población, bloque y rango"));
        }
        previous = Some(key);
        let entry = haplotypes
            .entry((row[0].clone(), row[1].clone()))
            .or_default();
        if rank != entry.len() as u64 {
            return Err(format!("{path}: haplotype_rank no es consecutivo desde 0"));
        }
        let q = canonical_f64(path, &row[4])?;
        if !(q > 0.0 && q <= 1.0) {
            return Err(format!("{path}: q={q} no pertenece a (0,1]"));
        }
        let alleles = &row[5];
        if alleles.len() != sites[block_index as usize].len()
            || !alleles.bytes().all(|byte| byte == b'0' || byte == b'1')
        {
            return Err(format!("{path}: alleles no son S caracteres 0/1"));
        }
        if row[2] != oracle.haplotype(&row[0], block_index, rank, alleles) {
            return Err(format!("{path}: la fila no rederiva su HAP_"));
        }
        entry.push((row[2].clone(), alleles.clone()));
        match &weights {
            None if rank == 0 => weights = Some(vec![row[4].clone()]),
            Some(values) if rank as usize == values.len() => {
                let mut values = values.clone();
                values.push(row[4].clone());
                weights = Some(values);
            }
            Some(values) if values[rank as usize] != row[4] => {
                return Err(format!("{path}: q depende del bloque o de la población"));
            }
            _ => {}
        }
    }
    let weights = weights.ok_or_else(|| format!("{path}: sin filas haplotípicas"))?;
    let total: f64 = weights
        .iter()
        .map(|value| value.parse::<f64>().expect("q ya validado"))
        .sum();
    if (total - 1.0).abs() > 1e-9 {
        return Err(format!(
            "{path}: los pesos no suman 1 por población y bloque"
        ));
    }
    if haplotypes.len() != populations.len() * blocks.len() {
        return Err(format!("{path}: falta alguna unidad bloque×subpoblación"));
    }

    // --- f5-frequencies.tsv --------------------------------------------------
    let path = F5_TSV_CONTRACTS[FREQUENCIES].path();
    let rows = read_tsv(package, FREQUENCIES)?;
    // La clave literal del fichero, sin derivarla de ningún otro productor:
    // `f5-frequencies.tsv` no publica `site_index`, de modo que su «sitio» es
    // la columna `site_id` y la terna debe crecer estrictamente.
    let mut previous: Option<(&str, &str, &str)> = None;
    let mut frequency_units: BTreeSet<(String, String)> = BTreeSet::new();
    let mut sites_per_unit: BTreeMap<(String, String), BTreeSet<String>> = BTreeMap::new();
    for row in &rows {
        let block_index = *block_index_of
            .get(&row[1])
            .ok_or_else(|| format!("{path}: {} no existe en el layout", row[1]))?;
        let unit = (row[0].clone(), row[1].clone());
        frequency_units.insert(unit.clone());
        let block_sites = &sites[block_index as usize];
        if !block_sites.iter().any(|site| site.site_id == row[2]) {
            return Err(format!(
                "{path}: el sitio {} no pertenece a su bloque",
                row[2]
            ));
        }
        if !sites_per_unit
            .entry(unit)
            .or_default()
            .insert(row[2].clone())
        {
            return Err(format!(
                "{path}: el sitio {} se repite en su unidad",
                row[2]
            ));
        }
        let key = (row[0].as_str(), row[1].as_str(), row[2].as_str());
        if previous.is_some_and(|last| last >= key) {
            return Err(format!(
                "{path}: el orden no es population_id, block_id y site_id estrictamente creciente"
            ));
        }
        previous = Some(key);
        if row[3] != "positive" && row[3] != "negative" {
            return Err(format!("{path}: polarity {:?} ajena", row[3]));
        }
        for value in row.iter().take(7).skip(4) {
            probability(path, value)?;
        }
        let target = probability(path, &row[4])?;
        let realized = probability(path, &row[5])?;
        let error = probability(path, &row[6])?;
        if ((realized - target).abs() - error).abs() > 1e-12 {
            return Err(format!("{path}: abs_error no es |realizado-objetivo|"));
        }
    }
    if frequency_units.len() != populations.len() * blocks.len() {
        return Err(format!("{path}: falta alguna unidad bloque×subpoblación"));
    }
    for ((population, block_id), published) in &sites_per_unit {
        let block_index = block_index_of[block_id];
        let expected: BTreeSet<String> = sites[block_index as usize]
            .iter()
            .map(|site| site.site_id.clone())
            .collect();
        if *published != expected {
            return Err(format!(
                "{path}: la unidad {population}×{block_id} no publica todos los sitios de su bloque"
            ));
        }
    }

    // --- f5-ld-evaluation.tsv ------------------------------------------------
    let path = F5_TSV_CONTRACTS[LD_EVALUATION].path();
    let rows = read_tsv(package, LD_EVALUATION)?;
    let mut previous: Option<(String, String, u64)> = None;
    let mut ld_units: BTreeSet<(String, String)> = BTreeSet::new();
    for row in &rows {
        if !block_index_of.contains_key(&row[1]) {
            return Err(format!("{path}: {} no existe en el layout", row[1]));
        }
        let start = canonical_u64(path, &row[2])?;
        let end = canonical_u64(path, &row[3])?;
        if start == 0 || start > end {
            return Err(format!("{path}: bin inválido"));
        }
        let key = (row[0].clone(), row[1].clone(), start);
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(format!(
                "{path}: el orden no es población, bloque e inicio de bin"
            ));
        }
        previous = Some(key);
        ld_units.insert((row[0].clone(), row[1].clone()));
        let pairs = canonical_u64(path, &row[5])?;
        match row[4].as_str() {
            "evaluated" => {
                if pairs == 0 {
                    return Err(format!(
                        "{path}: un bin evaluado no puede tener cero parejas"
                    ));
                }
                let mean = probability(path, &row[6])?;
                let target = probability(path, &row[7])?;
                let tolerance = probability(path, &row[8])?;
                let error = probability(path, &row[9])?;
                if ((mean - target).abs() - error).abs() > 1e-12 {
                    return Err(format!("{path}: abs_error no es |media-objetivo|"));
                }
                if error > tolerance {
                    return Err(format!("{path}: un bin publicado excede su tolerancia"));
                }
            }
            "not_applicable" => {
                if pairs != 0 {
                    return Err(format!("{path}: un bin no aplicable tiene parejas"));
                }
                if row[6] != "NA" || row[9] != "NA" {
                    return Err(format!("{path}: un bin no aplicable debe usar NA"));
                }
                probability(path, &row[7])?;
                probability(path, &row[8])?;
            }
            other => return Err(format!("{path}: status {other:?} ajeno")),
        }
        if row[4] == "evaluated" && (row[6] == "NA" || row[9] == "NA") {
            return Err(format!("{path}: NA fuera de un bin no aplicable"));
        }
    }
    if ld_units.len() != populations.len() * blocks.len() {
        return Err(format!("{path}: falta alguna unidad bloque×subpoblación"));
    }

    // --- f5-attempts.tsv -----------------------------------------------------
    let path = F5_TSV_CONTRACTS[ATTEMPTS].path();
    let rows = read_tsv(package, ATTEMPTS)?;
    let mut previous: Option<(String, String, u64)> = None;
    let mut best: BTreeMap<(String, String), u64> = BTreeMap::new();
    for row in &rows {
        if !block_index_of.contains_key(&row[1]) {
            return Err(format!("{path}: {} no existe en el layout", row[1]));
        }
        let attempt = canonical_u64(path, &row[2])?;
        let key = (row[0].clone(), row[1].clone(), attempt);
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(format!(
                "{path}: el orden no es población, bloque e intento"
            ));
        }
        previous = Some(key);
        for value in row.iter().take(8).skip(3) {
            canonical_u64(path, value)?;
        }
        for value in row.iter().take(11).skip(8) {
            canonical_f64(path, value)?;
        }
        match row[11].as_str() {
            "true" => *best.entry((row[0].clone(), row[1].clone())).or_default() += 1,
            "false" => {
                best.entry((row[0].clone(), row[1].clone())).or_default();
            }
            other => return Err(format!("{path}: selected_best {other:?} ajeno")),
        }
    }
    if best.len() != populations.len() * blocks.len() || best.values().any(|count| *count != 1) {
        return Err(format!(
            "{path}: cada bloque×subpoblación necesita exactamente un mejor intento"
        ));
    }

    // --- f5-founders.tsv -----------------------------------------------------
    let path = F5_TSV_CONTRACTS[FOUNDERS].path();
    let rows = read_tsv(package, FOUNDERS)?;
    let mut previous: Option<(String, String, String, u64, String)> = None;
    let mut alternate_from_tsv: BTreeMap<CopyKey, BTreeMap<String, (u64, char, char)>> =
        BTreeMap::new();
    let mut copy_membership: BTreeMap<CopyKey, (String, String)> = BTreeMap::new();
    for row in &rows {
        if !populations.contains(&row[0].as_str()) {
            return Err(format!(
                "{path}: {} no es una subpoblación publicada",
                row[0]
            ));
        }
        let copy_index = canonical_u64(path, &row[3])?;
        if copy_index >= ploidy {
            return Err(format!("{path}: copy_index fuera de 0..ploidy-1"));
        }
        let block_index = *block_index_of
            .get(&row[4])
            .ok_or_else(|| format!("{path}: {} no existe en el layout", row[4]))?;
        let key = (
            row[0].clone(),
            row[1].clone(),
            row[2].clone(),
            copy_index,
            row[4].clone(),
        );
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(format!(
                "{path}: el orden no es población, individuo, secuencia, copia y bloque"
            ));
        }
        previous = Some(key);
        let founder_copy_id = oracle.founder_copy(
            &row[0],
            &row[1],
            &row[2],
            u32::try_from(copy_index).map_err(|_| "copy_index desbordado".to_owned())?,
        );
        if row[5] != founder_copy_id {
            return Err(format!("{path}: la fila no rederiva su FCP_"));
        }
        if blocks[block_index as usize].sequence_id != row[2] {
            return Err(format!(
                "{path}: la copia no pertenece a la secuencia del bloque"
            ));
        }
        let unit = haplotypes
            .get(&(row[0].clone(), row[4].clone()))
            .ok_or_else(|| format!("{path}: la unidad no publica filas haplotípicas"))?;
        let rank = unit
            .iter()
            .position(|(haplotype_id, _)| *haplotype_id == row[7])
            .ok_or_else(|| format!("{path}: {} no existe en f5-haplotypes.tsv", row[7]))?;
        if row[6] != oracle.assignment(&row[0], block_index, &row[5], &row[7]) {
            return Err(format!("{path}: la fila no rederiva su FHA_"));
        }
        let identity = CopyKey {
            reference_id: truth.header.reference_id.clone(),
            individual_id: row[1].clone(),
            sequence_id: row[2].clone(),
            copy_index,
        };
        copy_membership.insert(identity.clone(), (row[0].clone(), row[5].clone()));
        let carried = alternate_from_tsv.entry(identity).or_default();
        for (site, allele) in sites[block_index as usize].iter().zip(unit[rank].1.bytes()) {
            if allele == b'1' {
                carried.insert(
                    site.site_id.clone(),
                    (
                        site.position_1based,
                        site.reference_base,
                        site.alternate_base,
                    ),
                );
            }
        }
    }

    // --- cruce con `truth.jsonl` --------------------------------------------
    if copy_membership != truth.membership {
        return Err(
            "cruce: las copias y su pertenencia no coinciden entre f5-founders.tsv y truth.jsonl"
                .to_owned(),
        );
    }
    if alternate_from_tsv != truth.alternate_sites {
        return Err(
            "cruce: los alelos reconstruidos desde los TSV no son los de truth.jsonl".to_owned(),
        );
    }
    for (site_id, (block_id, sequence_id, position, reference_base, alternate_base)) in &truth.locus
    {
        let block_index = *block_index_of
            .get(block_id)
            .ok_or_else(|| format!("cruce: el genotipo nombra el bloque {block_id}, ausente"))?;
        let site = sites[block_index as usize]
            .iter()
            .find(|candidate| candidate.site_id == *site_id)
            .ok_or_else(|| format!("cruce: el sitio {site_id} no pertenece a su bloque"))?;
        if site.position_1based != *position
            || site.reference_base != *reference_base
            || site.alternate_base != *alternate_base
            || blocks[block_index as usize].sequence_id != *sequence_id
        {
            return Err(format!("cruce: el locus {site_id} difiere de f5-sites.tsv"));
        }
    }
    let mut identifiers: Vec<String> = blocks.iter().map(|block| block.block_id.clone()).collect();
    identifiers.sort();
    let canonical_groups: Vec<(String, String)> = populations
        .iter()
        .flat_map(|population| {
            identifiers
                .iter()
                .map(move |block_id| ((*population).to_owned(), block_id.clone()))
        })
        .collect();
    if truth.block_order != canonical_groups {
        return Err(
            "cruce: los genotipos no recorren población y bloque en el orden canónico".to_owned(),
        );
    }
    let mut expected_dosage: BTreeMap<(String, String), BTreeMap<String, u64>> = BTreeMap::new();
    for (identity, carried) in &alternate_from_tsv {
        let population = &copy_membership[identity].0;
        for (block_index, block) in sites.iter().enumerate() {
            if blocks[block_index].sequence_id != identity.sequence_id {
                continue;
            }
            for site in block {
                let entry = expected_dosage
                    .entry((site.site_id.clone(), population.clone()))
                    .or_default()
                    .entry(identity.individual_id.clone())
                    .or_default();
                if carried.contains_key(&site.site_id) {
                    *entry += 1;
                }
            }
        }
    }
    if expected_dosage != truth.dosage {
        return Err("cruce: las dosis por locus×subpoblación no cuadran con los TSV".to_owned());
    }
    Ok(truth)
}

fn one_base(path: &str, value: &str) -> Result<char, String> {
    let mut characters = value.chars();
    match (characters.next(), characters.next()) {
        (Some(base @ ('A' | 'C' | 'G' | 'T')), None) => Ok(base),
        _ => Err(format!("{path}: {value:?} no es una base canónica")),
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

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

#[test]
fn the_independent_id_oracle_reproduces_the_frozen_catalog_vectors() {
    let catalog = fs::read_to_string("catalogs/id-catalog-v4.tsv").expect("catálogo de IDs");
    let mut checked = BTreeSet::new();
    for line in catalog.lines() {
        if line.starts_with('#') {
            continue;
        }
        let columns: Vec<&str> = line.split('|').collect();
        if columns.len() < 9 || !columns[0].starts_with("f5_") {
            continue;
        }
        let class = columns[0];
        let pieces: Vec<(&str, String)> = columns[4]
            .split(';')
            .filter(|piece| !piece.is_empty())
            .map(|piece| {
                let (name, value) = piece.split_once('=').expect("componente nombre=valor");
                (name, value.to_owned())
            })
            .collect();
        let (head, components) = pieces.split_first().expect("el primero es class");
        assert_eq!((head.0, head.1.as_str()), ("class", class));
        let root: [u8; 32] = decode_hex(columns[5]).try_into().expect("semilla");
        let material = id_material(&root, class, components);
        assert_eq!(
            encode_hex(&material),
            columns[6],
            "material canónico de {class}"
        );
        let digest = blake3::derive_key(ID_CONTEXT, &material);
        assert_eq!(encode_hex(&digest), columns[7], "digest de {class}");
        let prefix = columns[1];
        assert_eq!(
            derived_identifier(&root, prefix, class, components),
            columns[8],
            "vector congelado de {class}"
        );
        checked.insert(class.to_owned());
    }
    assert_eq!(
        checked.len(),
        5,
        "el catálogo publica un vector por cada una de las cinco clases F5"
    );
}

fn encode_hex(bytes: &[u8]) -> String {
    bytes.iter().map(|byte| format!("{byte:02x}")).collect()
}

#[test]
fn a_success_publishes_exactly_the_seven_tsv_and_the_founder_truth() {
    let published = publish("seven", &wide_fixture());
    let links = published.artifacts.artifacts();
    let identifiers: Vec<&str> = links.iter().map(|link| link.artifact_id()).collect();
    assert_eq!(
        identifiers,
        vec![
            "f5-01-block-layout-v1",
            "f5-02-sites-v1",
            "f5-03-haplotypes-v1",
            "f5-04-frequencies-v1",
            "f5-05-ld-evaluation-v1",
            "f5-06-founders-v1",
            "f5-07-attempts-v1",
        ]
    );
    let mut lexicographic = identifiers.clone();
    lexicographic.sort_unstable();
    assert_eq!(
        identifiers, lexicographic,
        "el orden de la tabla es el lexicográfico de artifact_id"
    );
    let kinds: Vec<&str> = links.iter().map(|link| link.artifact_kind()).collect();
    assert_eq!(
        kinds,
        vec![
            "f5_block_layout",
            "f5_sites",
            "f5_haplotypes",
            "f5_frequencies",
            "f5_ld_evaluation",
            "f5_founders",
            "f5_attempts",
        ]
    );
    for link in links {
        assert_eq!(link.schema_version(), "1");
        assert!(link.path().starts_with("outputs/"));
        assert!(!link.path().contains(".."));
        let bytes = published.bytes(link.path());
        assert_eq!(link.size_bytes(), bytes.len() as u64);
        assert_eq!(link.sha256(), sha256(&bytes));
        assert_eq!(
            link.to_json()
                .as_object()
                .expect("objeto")
                .keys()
                .map(String::as_str)
                .collect::<BTreeSet<_>>(),
            BTreeSet::from([
                "artifact_id",
                "artifact_kind",
                "schema_version",
                "path",
                "size_bytes",
                "sha256"
            ])
        );
    }
    let truth = published.artifacts.truth();
    assert_eq!(truth.artifact_id(), "truth-jsonl");
    assert_eq!(truth.artifact_kind(), "truth_jsonl");
    assert_eq!(truth.schema_version(), "2");
    assert_eq!(truth.projection_version(), "f5-founder-truth-projection-v2");
    assert_eq!(truth.path(), "outputs/truth.jsonl");
    let bytes = published.bytes(truth.path());
    assert_eq!(truth.size_bytes(), bytes.len() as u64);
    assert_eq!(truth.sha256(), sha256(&bytes));
    assert_eq!(
        truth
            .to_json()
            .as_object()
            .expect("objeto")
            .keys()
            .map(String::as_str)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([
            "artifact_id",
            "artifact_kind",
            "schema_version",
            "projection_version",
            "path",
            "size_bytes",
            "sha256"
        ])
    );
    assert!(
        !identifiers.contains(&"truth-jsonl"),
        "el enlace de verdad es un campo separado y no una octava fila"
    );
}

#[test]
fn the_published_package_satisfies_its_whole_contract() {
    let published = publish("contract", &wide_fixture());
    check_package(&published.package()).expect("paquete conforme");
}

#[test]
fn the_block_order_is_the_layout_and_not_the_lexicographic_identifier() {
    let published = publish("block-order", &wide_fixture());
    let package = published.package();
    let rows = read_tsv(&package, BLOCK_LAYOUT).expect("layout legible");
    let written: Vec<String> = rows.iter().map(|row| row[0].clone()).collect();
    let mut sorted = written.clone();
    sorted.sort();
    assert_ne!(
        written, sorted,
        "el fixture sería degenerado si los BLK_ ya estuvieran en orden lexicográfico"
    );
    check_package(&package).expect("paquete conforme");
}

#[test]
fn a_barren_sequence_contributes_no_copy_and_no_genotype() {
    let published = publish("barren", &wide_fixture());
    let package = published.package();
    let world = check_package(&package).expect("paquete conforme");
    assert!(
        world
            .membership
            .keys()
            .all(|copy| !copy.sequence_id.starts_with("barren")),
        "una secuencia que el layout no puebla no aporta ninguna copia"
    );
    let founders = read_tsv(&package, FOUNDERS).expect("fundadores");
    assert!(founders.iter().all(|row| !row[2].starts_with("barren")));
}

#[test]
fn the_truth_header_publishes_the_exact_projection_and_populations() {
    let published = publish("header", &wide_fixture());
    let body = published.text("outputs/truth.jsonl");
    let header: Value = serde_json::from_str(body.lines().next().expect("cabecera")).expect("JSON");
    assert_eq!(header["schema_version"], "2");
    assert_eq!(
        header["projection"],
        serde_json::json!({
            "kind": "f5_founder",
            "projection_version": "f5-founder-truth-projection-v2",
            "scenario_id": "scenario"
        })
    );
    assert_eq!(
        header["populations"],
        serde_json::json!([
            {"population_id": "PA", "individual_count": 2, "ploidy": 2},
            {"population_id": "PB", "individual_count": 1, "ploidy": 2}
        ])
    );
    let world = check_package(&published.package()).expect("paquete conforme");
    assert_eq!(world.copies, header["counts"]["copies"].as_u64().unwrap());
    assert_eq!(
        world.genotypes,
        header["counts"]["genotypes"].as_u64().unwrap()
    );
    // 3 fundadores × 2 secuencias pobladas × ploidía 2.
    assert_eq!(world.copies, 12);
    // 3 fundadores × 4 bloques × 3 sitios.
    assert_eq!(world.genotypes, 36);
}

#[test]
fn every_copy_carries_its_membership_and_every_genotype_its_locus() {
    let published = publish("membership", &wide_fixture());
    let body = published.text("outputs/truth.jsonl");
    let mut copies = 0_usize;
    let mut genotypes = 0_usize;
    for line in body.lines().skip(1) {
        let value: Value = serde_json::from_str(line).expect("JSON");
        match value["record_type"].as_str().expect("record_type") {
            "haplotype_copy" => {
                copies += 1;
                let membership = value["founder_membership"].as_object().expect("objeto");
                assert_eq!(
                    membership
                        .keys()
                        .map(String::as_str)
                        .collect::<BTreeSet<_>>(),
                    BTreeSet::from(["population_id", "founder_copy_id"])
                );
                let identity = value["identity"].as_object().expect("objeto");
                assert_eq!(
                    identity.keys().map(String::as_str).collect::<BTreeSet<_>>(),
                    BTreeSet::from(["reference_id", "individual_id", "sequence_id", "copy_index"])
                );
            }
            "genotype" => {
                genotypes += 1;
                let locus = value["f5_locus"].as_object().expect("objeto");
                assert_eq!(
                    locus.keys().map(String::as_str).collect::<BTreeSet<_>>(),
                    BTreeSet::from(["population_id", "block_id", "site_id", "ref", "alt"])
                );
            }
            other => panic!("registro inesperado {other}"),
        }
    }
    assert_eq!((copies, genotypes), (12, 36));
}

#[test]
fn deriving_the_assignments_twice_is_refused() {
    let root = temporary("twice");
    fs::create_dir_all(&root).expect("root");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let mut first = ExecutionWorkspace::create(&destination).expect("workspace");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
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
    write_f5_success_artifacts(&mut first, &request, &mut chain.identities)
        .expect("primera publicación");
    let mut second =
        ExecutionWorkspace::create(&chain.root.join("republished")).expect("workspace");
    let failure = write_f5_success_artifacts(&mut second, &request, &mut chain.identities)
        .expect_err("la segunda pasada deriva un FHA_ repetido");
    assert_eq!(failure.code(), "GPS-INTERNAL-101");
}

// ---------------------------------------------------------------------------
// Corrupción aislada: cada mutación debe ser detectada por el lector
// ---------------------------------------------------------------------------

fn replace_once(text: &str, from: &str, to: &str) -> String {
    let at = text
        .find(from)
        .unwrap_or_else(|| panic!("el fixture no contiene {from:?}"));
    format!("{}{to}{}", &text[..at], &text[at + from.len()..])
}

fn replace_after(text: &str, anchor: &str, from: &str, to: &str) -> String {
    let start = text
        .find(anchor)
        .unwrap_or_else(|| panic!("el fixture no contiene el ancla {anchor:?}"));
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

struct Mutation {
    name: &'static str,
    path: &'static str,
    apply: fn(&str) -> String,
}

const LAYOUT_PATH: &str = "outputs/f5-block-layout.tsv";
const SITES_PATH: &str = "outputs/f5-sites.tsv";
const HAPLOTYPES_PATH: &str = "outputs/f5-haplotypes.tsv";
const FREQUENCIES_PATH: &str = "outputs/f5-frequencies.tsv";
const LD_PATH: &str = "outputs/f5-ld-evaluation.tsv";
const FOUNDERS_PATH: &str = "outputs/f5-founders.tsv";
const ATTEMPTS_PATH: &str = "outputs/f5-attempts.tsv";
const TRUTH_PATH: &str = "outputs/truth.jsonl";

const MUTATIONS: &[Mutation] = &[
    // --- contrato de bytes ---
    Mutation {
        name: "BOM al principio",
        path: LAYOUT_PATH,
        apply: |text| format!("\u{feff}{text}"),
    },
    Mutation {
        name: "final de línea CRLF",
        path: LAYOUT_PATH,
        apply: |text| replace_once(text, "\n", "\r\n"),
    },
    Mutation {
        name: "sin LF final",
        path: LAYOUT_PATH,
        apply: |text| text.trim_end_matches('\n').to_owned(),
    },
    Mutation {
        name: "línea vacía",
        path: LAYOUT_PATH,
        apply: |text| replace_once(text, "\n", "\n\n"),
    },
    Mutation {
        name: "comentario",
        path: LAYOUT_PATH,
        apply: |text| replace_once(text, "\n", "\n# nota\n"),
    },
    Mutation {
        name: "columna de más",
        path: LAYOUT_PATH,
        apply: |text| replace_once(text, "\n", "\textra\n"),
    },
    Mutation {
        name: "cabecera alterada",
        path: LAYOUT_PATH,
        apply: |text| replace_once(text, "block_id\t", "blockid\t"),
    },
    Mutation {
        name: "celda entrecomillada",
        path: SITES_PATH,
        apply: |text| {
            let value = cell(text, 0, 5);
            set_cell(text, 0, 5, &format!("\"{value}\""))
        },
    },
    // --- layout ---
    Mutation {
        name: "layout desordenado",
        path: LAYOUT_PATH,
        apply: |text| swap_rows(text, 0, 1),
    },
    Mutation {
        name: "layout con BLK_ ajeno",
        path: LAYOUT_PATH,
        apply: |text| damage_identifier(text, 0, 0),
    },
    Mutation {
        name: "layout con inicio movido",
        path: LAYOUT_PATH,
        apply: |text| bump_number_after(text, &format!("{}\t", cell(text, 0, 1))),
    },
    Mutation {
        name: "layout con longitud incoherente",
        path: LAYOUT_PATH,
        apply: |text| {
            let length: u64 = cell(text, 0, 4).parse().expect("longitud");
            set_cell(text, 0, 4, &(length + 1).to_string())
        },
    },
    // --- sitios ---
    Mutation {
        name: "site_index no consecutivo",
        path: SITES_PATH,
        apply: |text| set_cell(text, 0, 1, "1"),
    },
    Mutation {
        name: "SITE_ ajeno",
        path: SITES_PATH,
        apply: |text| damage_identifier(text, 0, 2),
    },
    Mutation {
        name: "base de referencia cambiada",
        path: SITES_PATH,
        apply: |text| {
            let base = cell(text, 0, 5);
            set_cell(text, 0, 5, if base == "A" { "C" } else { "A" })
        },
    },
    Mutation {
        name: "posición del sitio movida",
        path: SITES_PATH,
        apply: |text| {
            let position: u64 = cell(text, 0, 4).parse().expect("posición");
            set_cell(text, 0, 4, &(position + 1).to_string())
        },
    },
    Mutation {
        name: "sitios desordenados",
        path: SITES_PATH,
        apply: |text| swap_rows(text, 0, 1),
    },
    // --- haplotipos ---
    Mutation {
        name: "alelo volteado",
        path: HAPLOTYPES_PATH,
        apply: |text| {
            let alleles = cell(text, 0, 5);
            let mut bytes = alleles.into_bytes();
            bytes[0] = if bytes[0] == b'0' { b'1' } else { b'0' };
            set_cell(text, 0, 5, &String::from_utf8(bytes).expect("UTF-8"))
        },
    },
    Mutation {
        name: "q sin decimal mínimo",
        path: HAPLOTYPES_PATH,
        apply: |text| set_cell(text, 0, 4, &format!("{}0", cell(text, 0, 4))),
    },
    Mutation {
        name: "rango haplotípico repetido",
        path: HAPLOTYPES_PATH,
        apply: |text| set_cell(text, 1, 3, "0"),
    },
    Mutation {
        name: "HAP_ ajeno",
        path: HAPLOTYPES_PATH,
        apply: |text| damage_identifier(text, 0, 2),
    },
    // --- frecuencias ---
    Mutation {
        name: "polaridad ajena",
        path: FREQUENCIES_PATH,
        apply: |text| set_cell(text, 0, 3, "POSITIVE"),
    },
    Mutation {
        name: "frecuencia realizada movida",
        path: FREQUENCIES_PATH,
        apply: |text| set_cell(text, 0, 5, "0.125"),
    },
    Mutation {
        name: "NA fuera de un bin",
        path: FREQUENCIES_PATH,
        apply: |text| set_cell(text, 0, 4, "NA"),
    },
    // --- LD ---
    Mutation {
        name: "estado no aplicable con parejas",
        path: LD_PATH,
        apply: |text| set_cell(text, 0, 4, "not_applicable"),
    },
    Mutation {
        name: "bin evaluado sin parejas",
        path: LD_PATH,
        apply: |text| set_cell(text, 0, 5, "0"),
    },
    Mutation {
        name: "NA en un bin evaluado",
        path: LD_PATH,
        apply: |text| set_cell(text, 0, 6, "NA"),
    },
    Mutation {
        name: "bins desordenados",
        path: LD_PATH,
        apply: |text| swap_rows(text, 0, 1),
    },
    // --- fundadores ---
    Mutation {
        name: "fundadores desordenados",
        path: FOUNDERS_PATH,
        apply: |text| swap_rows(text, 0, 1),
    },
    Mutation {
        name: "FCP_ ajeno",
        path: FOUNDERS_PATH,
        apply: |text| damage_identifier(text, 0, 5),
    },
    Mutation {
        name: "FHA_ ajeno",
        path: FOUNDERS_PATH,
        apply: |text| damage_identifier(text, 0, 6),
    },
    Mutation {
        name: "otra fila haplotípica asignada",
        path: FOUNDERS_PATH,
        apply: |text| {
            // El `HAP_` de otro rango de la misma unidad sí existe en su tabla
            // propietaria: lo que rompe es el `FHA_` y los alelos de la copia.
            let other = cell(text, 1, 7);
            if other == cell(text, 0, 7) {
                panic!("fixture degenerado: las dos primeras filas comparten haplotipo");
            }
            set_cell(text, 0, 7, &other)
        },
    },
    Mutation {
        name: "índice de copia cambiado",
        path: FOUNDERS_PATH,
        apply: |text| set_cell(text, 0, 3, "1"),
    },
    // --- intentos ---
    Mutation {
        name: "sin mejor intento",
        path: ATTEMPTS_PATH,
        apply: |text| replace_once(text, "\ttrue\n", "\tfalse\n"),
    },
    Mutation {
        name: "dos mejores intentos",
        path: ATTEMPTS_PATH,
        apply: |text| replace_once(text, "\tfalse\n", "\ttrue\n"),
    },
    Mutation {
        name: "selected_best no booleano",
        path: ATTEMPTS_PATH,
        apply: |text| replace_once(text, "\ttrue\n", "\t1\n"),
    },
    // --- proyección de verdad ---
    Mutation {
        name: "truth publicado como v1",
        path: TRUTH_PATH,
        apply: |text| replace_once(text, "\"schema_version\":\"2\"", "\"schema_version\":\"1\""),
    },
    Mutation {
        name: "versión de proyección ajena",
        path: TRUTH_PATH,
        apply: |text| {
            replace_once(
                text,
                "f5-founder-truth-projection-v2",
                "f5-founder-truth-projection-v1",
            )
        },
    },
    Mutation {
        name: "clase de proyección ajena",
        path: TRUTH_PATH,
        apply: |text| replace_once(text, "\"kind\":\"f5_founder\"", "\"kind\":\"f5_diploid\""),
    },
    Mutation {
        name: "proyección ausente",
        path: TRUTH_PATH,
        apply: |text| {
            let at = text.find("\"projection\":{").expect("projection");
            let end = text[at..].find('}').expect("fin") + at + 2;
            format!("{}{}", &text[..at], &text[end..])
        },
    },
    Mutation {
        name: "censo de la cabecera alterado",
        path: TRUTH_PATH,
        apply: |text| bump_number_after(text, "\"individual_count\":"),
    },
    Mutation {
        name: "recuento de copias alterado",
        path: TRUTH_PATH,
        apply: |text| bump_number_after(text, "\"counts\":{\"copies\":"),
    },
    Mutation {
        name: "recuento de genotipos alterado",
        path: TRUTH_PATH,
        apply: |text| bump_number_after(text, "\"genotypes\":"),
    },
    Mutation {
        name: "huella de la referencia truncada",
        path: TRUTH_PATH,
        apply: |text| replace_after(text, "\"reference_content_sha256\":\"", "a", "z"),
    },
    Mutation {
        name: "pertenencia de una copia cambiada",
        path: TRUTH_PATH,
        apply: |text| {
            replace_after(
                text,
                "\"founder_membership\":{",
                "\"population_id\":\"PA\"",
                "\"population_id\":\"PB\"",
            )
        },
    },
    Mutation {
        name: "población del locus cambiada",
        path: TRUTH_PATH,
        apply: |text| {
            replace_after(
                text,
                "\"f5_locus\":{",
                "\"population_id\":\"PA\"",
                "\"population_id\":\"PB\"",
            )
        },
    },
    Mutation {
        name: "dosis alterada",
        path: TRUTH_PATH,
        apply: |text| bump_number_after(text, "\"allele_dosages\":["),
    },
    Mutation {
        name: "portador con otro alelo",
        path: TRUTH_PATH,
        apply: |text| flip_digit_after(text, "\"carriers\":[{\"allele_index\":"),
    },
    Mutation {
        name: "fase distinta del índice de copia",
        path: TRUTH_PATH,
        apply: |text| bump_number_after(text, "\"phase\":"),
    },
    Mutation {
        name: "copias desordenadas",
        path: TRUTH_PATH,
        apply: |text| {
            let mut lines: Vec<String> = text.lines().map(str::to_owned).collect();
            lines.swap(1, 2);
            format!("{}\n", lines.join("\n"))
        },
    },
    Mutation {
        name: "alelo alternativo de una copia cambiado",
        path: TRUTH_PATH,
        apply: |text| rotate_base_after(text, "\"alternate\":\""),
    },
    Mutation {
        name: "sitio del locus cambiado",
        path: TRUTH_PATH,
        apply: |text| {
            replace_after(
                text,
                "\"f5_locus\":{",
                "\"site_id\":\"SITE_",
                "\"site_id\":\"SITX_",
            )
        },
    },
];

#[test]
fn every_isolated_corruption_is_detected_by_an_independent_reader() {
    let published = publish("mutations", &wide_fixture());
    let base = published.package();
    check_package(&base).expect("el paquete íntegro es válido");
    let mut detected = Vec::new();
    for mutation in MUTATIONS {
        let original = String::from_utf8(base.files[mutation.path].clone()).expect("UTF-8");
        let corrupted = (mutation.apply)(&original);
        assert_ne!(
            corrupted, original,
            "{}: la mutación no cambió ningún byte",
            mutation.name
        );
        let mut files = base.files.clone();
        files.insert(mutation.path.to_owned(), corrupted.into_bytes());
        let package = Package {
            files,
            seed: base.seed.clone(),
        };
        let error = check_package(&package)
            .err()
            .unwrap_or_else(|| panic!("«{}» pasó inadvertida", mutation.name));
        detected.push((mutation.name, error));
    }
    assert_eq!(detected.len(), MUTATIONS.len());
}

#[test]
fn a_missing_artifact_is_detected() {
    let published = publish("missing", &wide_fixture());
    let base = published.package();
    for path in [LAYOUT_PATH, SITES_PATH, FOUNDERS_PATH, TRUTH_PATH] {
        let mut files = base.files.clone();
        files.remove(path);
        let package = Package {
            files,
            seed: base.seed.clone(),
        };
        let error = check_package(&package)
            .err()
            .unwrap_or_else(|| panic!("falta {path} y el lector no lo notó"));
        assert!(error.contains(path), "{error}");
    }
}

#[test]
fn every_link_declares_the_literal_path_size_and_digest_of_its_bytes() {
    let published = publish("links", &wide_fixture());
    let expected = [
        "outputs/f5-block-layout.tsv",
        "outputs/f5-sites.tsv",
        "outputs/f5-haplotypes.tsv",
        "outputs/f5-frequencies.tsv",
        "outputs/f5-ld-evaluation.tsv",
        "outputs/f5-founders.tsv",
        "outputs/f5-attempts.tsv",
    ];
    for (link, path) in published.artifacts.artifacts().iter().zip(expected) {
        assert_eq!(link.path(), path);
        let bytes = published.bytes(path);
        assert_eq!(link.size_bytes(), bytes.len() as u64);
        assert_eq!(link.sha256(), sha256(&bytes));
        // Un byte distinto rompe tamaño y huella: el manifiesto no puede
        // declarar éxito sobre un fichero alterado.
        let mut altered = bytes.clone();
        altered.push(b'x');
        assert_ne!(link.size_bytes(), altered.len() as u64);
        assert_ne!(link.sha256(), sha256(&altered));
    }
    let truth = published.artifacts.truth();
    assert_eq!(truth.path(), "outputs/truth.jsonl");
    let bytes = published.bytes(truth.path());
    assert_eq!(truth.size_bytes(), bytes.len() as u64);
    assert_eq!(truth.sha256(), sha256(&bytes));
}

// ---------------------------------------------------------------------------
// Ámbitos de configuración: ploidías, modos, F, IDs y bins no aplicables
// ---------------------------------------------------------------------------

/// Bloques con hueco 2: las distancias realizadas son 2 y 4, de modo que los
/// bins 1 y 3 del perfil quedan `not_applicable` y los otros dos `evaluated`.
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

#[test]
fn both_ld_states_appear_and_only_the_inapplicable_one_uses_na() {
    let published = publish("ld-states", &sparse_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let rows = read_tsv(&package, LD_EVALUATION).expect("evaluación LD");
    let evaluated = rows.iter().filter(|row| row[4] == "evaluated").count();
    let inapplicable = rows.iter().filter(|row| row[4] == "not_applicable").count();
    assert!(
        evaluated > 0 && inapplicable > 0,
        "el fixture sería degenerado sin los dos estados: {evaluated}/{inapplicable}"
    );
    for row in &rows {
        let na = row[6] == "NA" || row[9] == "NA";
        assert_eq!(
            na,
            row[4] == "not_applicable",
            "NA solo existe en un bin no aplicable"
        );
        assert_eq!(row[5] == "0", row[4] == "not_applicable");
    }
}

#[test]
fn a_haploid_execution_publishes_one_copy_per_founder_and_sequence() {
    let published = publish(
        "haploid",
        &Fixture {
            ploidy: 1,
            subpopulations: vec![Subpopulation::new("P1", 2)],
            ..Fixture::default()
        },
    );
    let package = published.package();
    let world = check_package(&package).expect("paquete conforme");
    assert_eq!(world.header.populations[0].2, 1);
    assert_eq!(world.copies, 2);
    assert!(
        world.membership.keys().all(|copy| copy.copy_index == 0),
        "una ploidía 1 solo publica la copia 0"
    );
    for individuals in world.dosage.values() {
        assert!(individuals.values().all(|dosage| *dosage <= 1));
    }
}

#[test]
fn a_triploid_execution_publishes_three_phases_per_founder_and_sequence() {
    let published = publish(
        "triploid",
        &Fixture {
            ploidy: 3,
            subpopulations: vec![Subpopulation::new("P1", 1)],
            ..Fixture::default()
        },
    );
    let package = published.package();
    let world = check_package(&package).expect("paquete conforme");
    assert_eq!(world.header.populations[0].2, 3);
    assert_eq!(world.copies, 3);
    assert_eq!(
        world
            .membership
            .keys()
            .map(|copy| copy.copy_index)
            .collect::<BTreeSet<_>>(),
        BTreeSet::from([0, 1, 2])
    );
    let dosages: BTreeSet<u64> = world
        .dosage
        .values()
        .flat_map(BTreeMap::values)
        .copied()
        .collect();
    assert!(
        dosages.iter().all(|dosage| *dosage <= 3),
        "ninguna dosis puede exceder la ploidía: {dosages:?}"
    );
    assert!(
        dosages.contains(&1) || dosages.contains(&2),
        "el fixture sería degenerado sin ningún heterocigoto triploide: {dosages:?}"
    );
}

#[test]
fn a_declared_population_f_changes_the_published_genotypes() {
    let base = publish(
        "f-zero",
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 4)],
            ..Fixture::default()
        },
    );
    let inbred = publish(
        "f-one",
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 4).inbreeding("1.0")],
            ..Fixture::default()
        },
    );
    let with_f = check_package(&inbred.package()).expect("paquete conforme");
    check_package(&base.package()).expect("paquete conforme");
    // `F=1` produce homocigosis completa: ninguna dosis intermedia.
    for individuals in with_f.dosage.values() {
        for dosage in individuals.values() {
            assert!(*dosage == 0 || *dosage == 2, "F=1 excluye el heterocigoto");
        }
    }
    let heterozygous = check_package(&base.package())
        .expect("paquete conforme")
        .dosage
        .values()
        .flat_map(BTreeMap::values)
        .filter(|dosage| **dosage == 1)
        .count();
    assert!(
        heterozygous > 0,
        "el fixture sería degenerado si `F=0` tampoco produjera heterocigotos"
    );
}

#[test]
fn provided_and_derived_founder_identifiers_publish_the_same_contract() {
    let provided = publish(
        "provided",
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 2).founders(&["zeta", "alfa"])],
            ..Fixture::default()
        },
    );
    let derived = publish(
        "derived",
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 2).keys()],
            ..Fixture::default()
        },
    );
    let provided_world = check_package(&provided.package()).expect("paquete conforme");
    let derived_world = check_package(&derived.package()).expect("paquete conforme");
    assert_eq!(provided_world.copies, derived_world.copies);
    assert_eq!(provided_world.genotypes, derived_world.genotypes);
    let literal: Vec<String> = provided_world
        .membership
        .keys()
        .map(|copy| copy.individual_id.clone())
        .collect();
    assert_eq!(literal.first().map(String::as_str), Some("alfa"));
    assert!(
        literal.contains(&"zeta".to_owned()),
        "los identificadores aportados se publican literalmente y ordenados"
    );
    assert!(
        derived_world
            .membership
            .keys()
            .all(|copy| copy.individual_id.starts_with("IND_")),
        "una clave fundadora deriva su identidad"
    );
}

#[test]
fn a_controlled_only_and_a_probabilistic_only_execution_both_publish() {
    for (label, mode) in [
        ("controlled", "controlled"),
        ("probabilistic", "probabilistic"),
    ] {
        let published = publish(
            label,
            &Fixture {
                subpopulations: vec![Subpopulation::new("P1", 2).sampling(mode)],
                ..Fixture::default()
            },
        );
        check_package(&published.package()).unwrap_or_else(|error| panic!("modo {mode}: {error}"));
    }
}

#[test]
fn the_same_seed_publishes_the_same_bytes() {
    let first = publish("repeat-a", &wide_fixture());
    let second = publish("repeat-b", &wide_fixture());
    for (left, right) in first
        .artifacts
        .artifacts()
        .iter()
        .zip(second.artifacts.artifacts())
    {
        assert_eq!(left, right, "{}", left.artifact_id());
        assert_eq!(first.bytes(left.path()), second.bytes(right.path()));
    }
    assert_eq!(first.artifacts.truth(), second.artifacts.truth());
    assert_eq!(
        first.bytes("outputs/truth.jsonl"),
        second.bytes("outputs/truth.jsonl")
    );
}

#[test]
fn a_sequence_without_a_declared_length_cannot_be_published() {
    let failure = publish_with_sequences(
        "no-length",
        &Fixture::default(),
        F5ReferenceSequences::new(),
    );
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("reference.sequence_lengths")
    );
}

#[test]
fn a_block_that_does_not_fit_its_declared_sequence_is_refused() {
    let mut sequences = F5ReferenceSequences::new();
    assert!(sequences.insert(SequenceId::parse("chr1", "test").expect("sequence"), 1));
    let failure = publish_with_sequences("short-sequence", &Fixture::default(), sequences);
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("reference.sequence_lengths")
    );
}

#[test]
fn a_zero_length_or_repeated_sequence_is_not_accepted_by_the_declaration() {
    let mut sequences = F5ReferenceSequences::new();
    let identifier = SequenceId::parse("chr1", "test").expect("sequence");
    assert!(!sequences.insert(identifier.clone(), 0));
    assert!(sequences.insert(identifier.clone(), 10));
    assert!(!sequences.insert(identifier.clone(), 11));
    assert_eq!(sequences.length(&identifier), Some(10));
    assert_eq!(sequences.len(), 1);
    assert!(!sequences.is_empty());
}

// ---------------------------------------------------------------------------
// Fidelidad: el escritor publica los valores que sus productores calcularon
// ---------------------------------------------------------------------------
//
// El lector independiente comprueba el contrato de los bytes, pero hay
// columnas cuya verdad no es rederivable desde los propios ficheros: la
// polaridad de cada sitio, las tres frecuencias, los diez campos LD, los nueve
// contadores de intento y `q`. Sustituir cualquiera de ellas por una constante
// plausible dejaría un paquete internamente coherente. Aquí se contrastan
// contra el valor que F5.2 publicó, que es la verdad del escritor.

#[test]
fn the_polarity_column_is_not_constant_and_is_the_one_each_site_received() {
    let published = publish("polarity", &wide_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let rows = read_tsv(&package, FREQUENCIES).expect("frecuencias");
    let polarities: BTreeSet<&str> = rows.iter().map(|row| row[3].as_str()).collect();
    assert_eq!(
        polarities,
        BTreeSet::from(["negative", "positive"]),
        "el fixture sería degenerado si todos los sitios compartieran polaridad"
    );

    // El `expected` no hereda el orden de emisión del productor: se ordena
    // por la clave literal del fichero —población, `block_id` y `site_id`—,
    // que es la que DTC-089 fija y la única columna de sitio que publica.
    let mut expected = Vec::new();
    for unit in published.chain.identities.units() {
        let block = published
            .chain
            .search
            .block(unit.population_id(), unit.block_index())
            .expect("bloque ganador");
        let mut unit_rows: Vec<(String, String, String, String, String, String, String)> = block
            .columns()
            .iter()
            .map(|column| {
                (
                    unit.population_id().as_str().to_owned(),
                    unit.block_id().as_str().to_owned(),
                    column.site_id().as_str().to_owned(),
                    column.polarity().as_str().to_owned(),
                    column.target_alt_af().to_string(),
                    column.realized_alt_af().to_string(),
                    column.abs_error().to_string(),
                )
            })
            .collect();
        unit_rows.sort_by(|left, right| left.2.cmp(&right.2));
        expected.extend(unit_rows);
    }
    expected
        .sort_by(|left, right| (&left.0, &left.1, &left.2).cmp(&(&right.0, &right.1, &right.2)));
    let observed: Vec<_> = rows
        .iter()
        .map(|row| {
            (
                row[0].clone(),
                row[1].clone(),
                row[2].clone(),
                row[3].clone(),
                row[4].clone(),
                row[5].clone(),
                row[6].clone(),
            )
        })
        .collect();
    assert_eq!(observed, expected);
}

#[test]
fn the_published_ld_rows_are_the_producers_own_bins() {
    let published = publish("ld-fidelity", &sparse_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let rows = read_tsv(&package, LD_EVALUATION).expect("evaluación LD");
    let mut expected = Vec::new();
    for unit in published.chain.identities.units() {
        let evaluation = published
            .chain
            .search
            .ld_evaluation()
            .iter()
            .find(|candidate| {
                candidate.population_id() == unit.population_id()
                    && candidate.block_index() == unit.block_index()
            })
            .expect("evaluación del bloque");
        for bin in evaluation.bins() {
            expected.push(vec![
                unit.population_id().as_str().to_owned(),
                unit.block_id().as_str().to_owned(),
                bin.distance_start_bp().to_string(),
                bin.distance_end_bp().to_string(),
                bin.status().to_owned(),
                bin.pair_count().to_string(),
                bin.mean_r2()
                    .map_or_else(|| "NA".to_owned(), |value| value.to_string()),
                bin.target_mean_r2().to_string(),
                bin.max_abs_error().to_string(),
                bin.abs_error()
                    .map_or_else(|| "NA".to_owned(), |value| value.to_string()),
            ]);
        }
    }
    assert_eq!(rows, expected);
}

#[test]
fn the_published_attempts_are_the_producers_own_records() {
    let published = publish("attempts-fidelity", &wide_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let rows = read_tsv(&package, ATTEMPTS).expect("intentos");
    let mut expected = Vec::new();
    for unit in published.chain.identities.units() {
        for record in published.chain.search.records() {
            if record.population_id() != unit.population_id()
                || record.block_index() != unit.block_index()
            {
                continue;
            }
            expected.push(vec![
                record.population_id().as_str().to_owned(),
                record.block_id().as_str().to_owned(),
                record.attempt_index().to_string(),
                record.proposals().to_string(),
                record.accepted().to_string(),
                record.duplicate_rows().to_string(),
                record.bins_outside().to_string(),
                record.signs_failed().to_string(),
                record.max_r2_error().to_string(),
                record.sum_r2_error().to_string(),
                record.max_marginal_error().to_string(),
                record.selected_best().to_string(),
            ]);
        }
    }
    assert_eq!(rows, expected);
    assert!(
        rows.iter().any(|row| row[8] != "0") || rows.iter().any(|row| row[3] != "0"),
        "el fixture sería degenerado si todos los contadores publicados fueran cero"
    );
}

/// Perfil diádico de cuatro pesos **distintos**: `1/2, 1/4, 1/8, 1/8`.
///
/// Un perfil uniforme haría degenerada la columna `q`, y una `q` sustituida
/// por la del primer rango pasaría inadvertida.
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

#[test]
fn the_published_weights_are_the_effective_profile_weights() {
    let published = publish("weights", &dyadic_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let rows = read_tsv(&package, HAPLOTYPES).expect("haplotipos");
    let weights = published.chain.config.profile().haplotype_weights();
    assert!(
        weights.windows(2).any(|pair| pair[0] != pair[1]),
        "el fixture sería degenerado con pesos uniformes: una `q` constante pasaría inadvertida"
    );
    for row in &rows {
        let rank: usize = row[3].parse().expect("rango");
        assert_eq!(row[4], weights[rank].to_string());
    }
    assert_eq!(
        rows.len(),
        published.chain.identities.units().len() * weights.len()
    );
    assert_eq!(
        rows.iter()
            .map(|row| row[4].clone())
            .collect::<BTreeSet<_>>(),
        weights
            .iter()
            .map(std::string::ToString::to_string)
            .collect::<BTreeSet<_>>()
    );
}

#[test]
fn the_published_layout_and_sites_are_the_producers_own_geometry() {
    let published = publish("geometry", &wide_fixture());
    let package = published.package();
    check_package(&package).expect("paquete conforme");
    let layout_rows = read_tsv(&package, BLOCK_LAYOUT).expect("layout");
    let expected: Vec<Vec<String>> = published
        .chain
        .layout
        .blocks()
        .iter()
        .map(|block| {
            vec![
                block.block_id().as_str().to_owned(),
                block.sequence_id().as_str().to_owned(),
                block.start_1based().to_string(),
                block.end_1based().to_string(),
                block.length_bp().to_string(),
            ]
        })
        .collect();
    assert_eq!(layout_rows, expected);

    let site_rows = read_tsv(&package, SITES).expect("sitios");
    let mut ordered: Vec<&_> = published.chain.layout.blocks().iter().collect();
    ordered.sort_by(|left, right| left.block_id().cmp(right.block_id()));
    let expected: Vec<Vec<String>> = ordered
        .iter()
        .flat_map(|block| {
            block.sites().iter().map(move |site| {
                vec![
                    block.block_id().as_str().to_owned(),
                    site.site_index().to_string(),
                    site.site_id().as_str().to_owned(),
                    block.sequence_id().as_str().to_owned(),
                    site.position_1based().to_string(),
                    char::from(site.reference_base()).to_string(),
                    char::from(site.alternate_base()).to_string(),
                ]
            })
        })
        .collect();
    assert_eq!(site_rows, expected);
}

#[test]
fn the_published_founder_rows_are_the_single_assignment_traversal() {
    let root = temporary("founder-fidelity");
    fs::create_dir_all(&root).expect("root");
    let mut chain = build_chain(&root, &wide_fixture());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let mut visited = Vec::new();
    {
        let request = F5SuccessRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
            &chain.sequences,
        )
        .expect("solicitud");
        write_f5_success_artifacts(&mut workspace, &request, &mut chain.identities)
            .expect("paquete escrito");
    }
    let payload = payload_of(&chain.root, &destination);
    let rows: Vec<Vec<String>> = fs::read_to_string(payload.join(FOUNDERS_PATH))
        .expect("fundadores")
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect();
    // Una segunda cadena idéntica reproduce el recorrido: la pasada del
    // escritor es exactamente la de `for_each_assignment`, sin reordenar.
    let second_root = temporary("founder-fidelity-mirror");
    fs::create_dir_all(&second_root).expect("root");
    let mut mirror = build_chain(&second_root, &wide_fixture());
    mirror
        .identities
        .for_each_assignment(&mirror.assignment, |row| {
            visited.push(vec![
                row.population_id().as_str().to_owned(),
                row.individual_id().as_str().to_owned(),
                row.sequence_id().as_str().to_owned(),
                row.copy_index().to_string(),
                row.block_id().as_str().to_owned(),
                row.founder_copy_id().as_str().to_owned(),
                row.block_assignment_id().as_str().to_owned(),
                row.haplotype_id().as_str().to_owned(),
            ]);
            Ok(())
        })
        .expect("recorrido");
    assert_eq!(rows, visited);
    assert!(rows.len() > 1);
}

#[test]
fn the_recount_by_locus_and_subpopulation_comes_only_from_truth() {
    let published = publish("recount", &wide_fixture());
    let world = check_package(&published.package()).expect("paquete conforme");
    let declared: BTreeMap<&str, (u64, u64)> = world
        .header
        .populations
        .iter()
        .map(|(id, individuals, ploidy)| (id.as_str(), (*individuals, *ploidy)))
        .collect();
    assert!(!world.counts.is_empty());
    for ((site_id, population), histogram) in &world.counts {
        let (individuals, ploidy) = declared[population.as_str()];
        assert_eq!(
            histogram.values().sum::<u64>(),
            individuals,
            "el locus {site_id}×{population} debe contar exactamente sus individuos"
        );
        assert!(histogram.keys().all(|dosage| *dosage <= ploidy));
        let from_histogram: u64 = histogram.iter().map(|(dosage, count)| dosage * count).sum();
        let from_individuals: u64 = world.dosage[&(site_id.clone(), population.clone())]
            .values()
            .sum();
        assert_eq!(from_histogram, from_individuals);
    }
    // Un locus por sitio y subpoblación: la geometría es común y cada
    // subpoblación cubre todos los bloques.
    assert_eq!(world.counts.len(), world.locus.len() * declared.len());
    assert!(
        world.counts.values().any(|histogram| histogram.len() > 1),
        "el fixture sería degenerado si todo locus tuviera una sola dosis"
    );
}

// ---------------------------------------------------------------------------
// Fidelidad de `truth.jsonl`: la sexta, simétrica a las cinco de los TSV
// ---------------------------------------------------------------------------

/// `case_id` reimplementado aquí desde DTC-061: SHA-256 de población,
/// individuo, secuencia, posición 1-based y la lista de claves de alelo, cada
/// campo terminado en `NUL`, con el prefijo `CASE_` y 59 hex.
fn oracle_case_id(
    population: &str,
    individual: &str,
    sequence: &str,
    position_1based: u64,
    alternate: char,
) -> String {
    let mut hasher = Sha256::new();
    for value in [
        population,
        individual,
        sequence,
        &position_1based.to_string(),
        &format!("REF,SNP:{alternate}"),
    ] {
        hasher.update(value.as_bytes());
        hasher.update([0]);
    }
    format!("CASE_{}", &format!("{:x}", hasher.finalize())[..59])
}

/// Los campos de `truth.jsonl` cuya verdad **no** vive en sus propios bytes:
/// las dos huellas de la cabecera, la traza FHA→HAP de cada SNP, la
/// pertenencia que esa traza declara y el `case_id` de cada genotipo.
///
/// Sustituir cualquiera de ellos por una constante plausible o por el valor
/// vecino dejaba hasta ahora un documento internamente coherente: el lector
/// solo comprobaba su forma.
#[test]
fn the_published_truth_carries_the_digests_the_trace_and_the_cases_it_received() {
    let published = publish("truth-fidelity", &wide_fixture());
    let package = published.package();
    let world = check_package(&package).expect("paquete conforme");
    let sealed = &published.chain.sealed;

    // --- las dos huellas de la cabecera -------------------------------------
    assert_eq!(
        world.header.reference_content_sha256,
        sealed.reference_sha256(),
        "la huella de contenido es la del FASTA sellado"
    );
    assert_ne!(
        world.header.reference_content_sha256,
        sealed.profile_metadata_sha256(),
        "el fixture sería degenerado si las dos huellas coincidieran"
    );
    assert_eq!(
        world.header.source_digests,
        vec![
            (
                "reference".to_owned(),
                sealed.reference_id().to_owned(),
                sealed.reference_sha256().to_owned(),
            ),
            (
                "matrix".to_owned(),
                sealed.profile_id().to_owned(),
                sealed.profile_metadata_sha256().to_owned(),
            ),
        ],
        "`source_digests` son exactamente la referencia y el perfil con sus digests"
    );

    // --- la traza FHA→HAP de cada SNP ---------------------------------------
    let site_rows = read_tsv(&package, SITES).expect("sitios");
    let block_of_site: BTreeMap<String, String> = site_rows
        .iter()
        .map(|row| (row[2].clone(), row[0].clone()))
        .collect();
    let founder_rows = read_tsv(&package, FOUNDERS).expect("fundadores");
    let mut assigned: BTreeMap<(String, String, u64, String), String> = BTreeMap::new();
    for row in &founder_rows {
        let copy_index: u64 = row[3].parse().expect("copy_index");
        assert!(
            assigned
                .insert(
                    (row[1].clone(), row[2].clone(), copy_index, row[4].clone()),
                    row[7].clone(),
                )
                .is_none(),
            "f5-founders.tsv no repite copia×bloque"
        );
    }
    let mut checked = 0_usize;
    let mut haplotypes_seen: BTreeSet<String> = BTreeSet::new();
    for (identity, sites) in &world.snp_provenance {
        let population = &world.membership[identity].0;
        for (site_id, (matrix_source_id, declared_population)) in sites {
            let block_id = block_of_site
                .get(site_id)
                .unwrap_or_else(|| panic!("el sitio {site_id} no existe en f5-sites.tsv"));
            let expected = assigned
                .get(&(
                    identity.individual_id.clone(),
                    identity.sequence_id.clone(),
                    identity.copy_index,
                    block_id.clone(),
                ))
                .unwrap_or_else(|| panic!("sin asignación para {site_id}"));
            assert_eq!(
                matrix_source_id, expected,
                "el SNP {site_id} debe citar el HAP_ que f5-founders.tsv asignó a esa copia"
            );
            assert_eq!(declared_population, population);
            haplotypes_seen.insert(matrix_source_id.clone());
            checked += 1;
        }
    }
    assert!(checked > 0, "el fixture sería degenerado sin ningún SNP");
    assert!(
        haplotypes_seen.len() > 1,
        "el fixture sería degenerado si un solo HAP_ cubriera todos los SNP: {haplotypes_seen:?}"
    );

    // --- el `case_id` de cada genotipo --------------------------------------
    assert_eq!(world.cases.len() as u64, world.genotypes);
    let mut tuples = BTreeSet::new();
    for (case, (population, individual, sequence, position, alternate)) in &world.cases {
        assert_eq!(
            *case,
            oracle_case_id(population, individual, sequence, *position, *alternate),
            "el case_id no es el que su tupla determina"
        );
        assert!(
            tuples.insert((
                population.clone(),
                individual.clone(),
                sequence.clone(),
                *position
            )),
            "dos genotipos comparten población, individuo, secuencia y posición"
        );
    }
    assert_eq!(
        tuples.len() as u64,
        world.genotypes,
        "un caso por población, individuo, secuencia y posición"
    );
}

// ---------------------------------------------------------------------------
// Ataduras de la solicitud con la admisión y con las identidades
// ---------------------------------------------------------------------------

#[test]
fn publishing_under_an_execution_other_than_the_one_that_derived_the_identities_is_refused() {
    let root = temporary("alien-execution");
    fs::create_dir_all(&root).expect("root");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
    let alien = ExecutionId::parse("RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("execution");
    assert_ne!(alien.as_str(), EXECUTION);
    let failure = {
        let request = F5SuccessRequest::new(
            &chain.config,
            &chain.sealed,
            &alien,
            &chain.layout,
            &chain.search,
            &chain.assignment,
            &chain.sequences,
        )
        .expect("la solicitud no puede ver todavía el repositorio");
        write_f5_success_artifacts(&mut workspace, &request, &mut chain.identities)
            .expect_err("una ejecución ajena a las identidades no puede publicarse")
    };
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("execution.identity")
    );
}

#[test]
fn a_sequence_declaration_wider_than_the_sealed_reference_is_refused() {
    // Más secuencias que las que la admisión retuvo.
    let mut extra = F5ReferenceSequences::new();
    assert!(extra.insert(SequenceId::parse("chr1", "test").expect("sequence"), 6));
    assert!(extra.insert(SequenceId::parse("chr9", "test").expect("sequence"), 6));
    let failure = publish_with_sequences("too-many-sequences", &Fixture::default(), extra);
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("reference.sequence_lengths")
    );

    // Una longitud inventada: la suma declarada excede las bases recorridas.
    let mut inflated = F5ReferenceSequences::new();
    assert!(inflated.insert(
        SequenceId::parse("chr1", "test").expect("sequence"),
        1_000_000_000
    ));
    let failure = publish_with_sequences("inflated-sequence", &Fixture::default(), inflated);
    assert_eq!(failure.code(), "GPS-CONF-101");
    assert_eq!(
        failure.context().get("fields").map(String::as_str),
        Some("reference.sequence_lengths")
    );
}

// ---------------------------------------------------------------------------
// El `truth.jsonl` **publicado** contra `truth-v2.schema.json`
// ---------------------------------------------------------------------------

/// Copia literal del arnés de `tests/f5_0_1_external_schemas.rs`: aísla el
/// esquema y sus dependencias en un directorio propio y delega la validación
/// en `jsonschema` de Python, que no comparte código con el productor.
fn assert_external_cases(schema: &str, dependencies: &[&str], cases: Vec<(Value, bool)>) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let bundle = temporary(&format!("schema-{schema}"));
    fs::create_dir_all(&bundle).expect("isolated schema bundle");
    for name in std::iter::once(schema).chain(dependencies.iter().copied()) {
        assert!(
            !name.contains('/'),
            "el paquete solo admite nombres hermanos"
        );
        fs::copy(root.join("schemas").join(name), bundle.join(name)).expect("copy schema");
    }
    let case_path = bundle.join("cases.json");
    let case_document = cases
        .into_iter()
        .map(|(instance, expected)| json!({"instance": instance, "expected": expected}))
        .collect::<Vec<_>>();
    fs::write(
        &case_path,
        serde_json::to_vec(&case_document).expect("JSON cases"),
    )
    .expect("write cases");
    let output = std::process::Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(bundle.join(schema))
        .arg(&case_path)
        .arg("--cases")
        .stderr(std::process::Stdio::piped())
        .output()
        .expect("external validator");
    fs::remove_dir_all(&bundle).expect("remove isolated schema bundle");
    assert!(
        output.status.success(),
        "{schema}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

/// Cada línea **publicada** valida contra el esquema externo, en las tres
/// ploidías cubiertas.
///
/// No es lo mismo que validar registros sintéticos: el `case_id` publicado
/// mide exactamente 64 caracteres, el máximo del patrón `$defs/id`, y sin esta
/// prueba nada vigilaría ese cero de margen.
#[test]
fn every_published_truth_line_validates_against_the_external_schema() {
    let mut lengths = BTreeSet::new();
    for ploidy in [1_u64, 2, 3] {
        let published = publish(
            &format!("schema-ploidy-{ploidy}"),
            &Fixture {
                ploidy,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
        );
        check_package(&published.package()).expect("paquete conforme");
        let document = published.text("outputs/truth.jsonl");
        let mut cases = Vec::new();
        for line in document.lines() {
            let value: Value = serde_json::from_str(line).expect("línea JSON");
            if let Some(case) = value.get("case_id").and_then(Value::as_str) {
                lengths.insert(case.chars().count());
            }
            cases.push((value, true));
        }
        assert!(
            cases.len() > 1,
            "una ploidía {ploidy} publica cabecera, copias y genotipos"
        );
        assert_external_cases("truth-v2.schema.json", &[], cases);
    }
    assert_eq!(
        lengths,
        BTreeSet::from([64]),
        "el case_id publicado agota el patrón `$defs/id`: 64 caracteres exactos"
    );
}

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

/// Lo que el escritor asigna **cabe en el pico que la admisión selló**, y no
/// crece con el número de bloques.
///
/// Es la relación que faltaba y que habría detectado el defecto sin revisión:
/// la versión anterior materializaba el registro `haplotype_copy` con tres
/// `collect::<Vec<Value>>()` vivos a la vez y picaba 2,32× el sellado con
/// `B=8`, 4,39× con `B=32` y 6,23× con `B=128`, creciendo linealmente con los
/// bloques mientras la admisión sellaba el mismo presupuesto.
///
/// Ahora los únicos términos que crecen con `B` son los índices de posición
/// —`B` para el layout, `B` ordenados por `block_id`, `B*N_pop` para las
/// unidades y sus tramos, y hasta `S_bound` índices por bloque para el orden
/// de sitio—, y la prueba acota su coste por bloque.
#[test]
fn the_writer_peak_fits_the_sealed_budget_and_does_not_grow_with_the_block_count() {
    let mut measured = Vec::new();
    for blocks in [8_u64, 32] {
        let root = temporary(&format!("writer-peak-{blocks}"));
        fs::create_dir_all(&root).expect("root");
        let mut chain = build_chain(
            &root,
            &Fixture {
                block_count: blocks,
                // Holgura para que los bloques quepan; la geometría no
                // interviene en lo que se mide, solo `B`.
                slack: 2,
                ..Fixture::default()
            },
        );
        let destination = chain.root.join("published");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
        let sealed_peak = chain.sealed.working_memory_bytes();
        let peak = {
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
            let (artifacts, peak) = peak_bytes(|| {
                write_f5_success_artifacts(&mut workspace, &request, &mut chain.identities)
            });
            artifacts.expect("paquete F5.4.1 escrito");
            peak
        };
        // La promesa que la admisión selló para toda la fase.
        assert!(
            peak <= sealed_peak,
            "B={blocks}: el escritor pica {peak} bytes y la admisión selló {sealed_peak}"
        );
        measured.push((blocks, peak, sealed_peak));
    }

    let (small_blocks, small_peak, _) = measured[0];
    let (large_blocks, large_peak, _) = measured[1];
    // Cuadruplicar `B` no puede cuadruplicar el pico: lo único que crece son
    // los índices de posición. Un kibibyte por bloque los cubre con holgura y
    // queda un orden de magnitud por debajo de lo que costaba materializar
    // los registros, que eran unos 25 KiB por bloque en este mismo fixture.
    let allowance = (large_blocks - small_blocks) * 1024;
    assert!(
        large_peak <= small_peak + allowance,
        "el pico creció {} bytes al pasar de {small_blocks} a {large_blocks} bloques, \
         más de los {allowance} que los índices de posición necesitan",
        large_peak.saturating_sub(small_peak)
    );
}
