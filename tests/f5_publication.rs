//! F5.4.3. Configuración efectiva v5, manifiesto v5, fallo y publicación
//! atómica.
//!
//! La batería no juzga el publicador con el publicador. Cada comprobación usa
//! una de estas tres fuentes independientes:
//!
//! 1. **Los bytes publicados**, releídos del destino final y contrastados con
//!    los ficheros que enlazan: tamaño y SHA-256 se recomputan aquí.
//! 2. **El oráculo de F5.4.2**, al que se le entrega una declaración
//!    construida **desde el manifiesto publicado** —no desde el productor— y
//!    desde los TSV de perfil que el fixture escribió. La única pieza que no
//!    puede reconstruirse así es `normalizations`, cuyo tipo no tiene
//!    constructor público; se toma del perfil admitido y se comprueba aparte
//!    contra las cuatro auditorías publicadas.
//! 3. **El validador JSON Schema externo en Python**, ajeno al serializador
//!    Rust, sobre los tres esquemas que este hito instancia.
//!
//! La igualdad exacta entre efectivo y manifiesto se comprueba releyendo los
//! **dos ficheros del destino**, nunca dos salidas del mismo serializador en
//! memoria, y su discriminación se mide mutando un término de una copia y
//! exigiendo que la puerta física lo tumbe citando ese término.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::collections::BTreeSet;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::f5_artifacts::{
    F5_TSV_CONTRACTS, F5ReferenceSequences, F5SuccessRequest, write_f5_success_artifacts,
};
use genome_population_simulator::f5_founders::{
    F5FounderAssignment, F5FounderIdentities, assign_f5_founders, materialize_f5_founder_identities,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_oracle::{
    F5DeclaredAttemptSummary, F5DeclaredBounds, F5DeclaredLdBin, F5DeclaredLink, F5DeclaredMode,
    F5DeclaredRealizedF, F5DeclaredSampling, F5DeclaredSpectrumBand, F5DeclaredSubpopulation,
    F5DeclaredTruthLink, F5FailureDeclaration, F5PackageDeclaration, audit_f5_failure_package,
    audit_f5_package,
};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_publication::{
    F5_EFFECTIVE_CONFIG_PATH, F5_FAILURE_ARTIFACT_ID, F5_FAILURE_PATH, F5_MANIFEST_PATH,
    F5FailureCode, F5FailureReport, F5FailureStage, F5PublicationRequest, publish_f5_failure,
    publish_f5_success, verify_published_f5_package,
};
use genome_population_simulator::f5_search::{
    F5HaplotypeSearch, F5SearchFailure, search_f5_haplotypes,
};
use genome_population_simulator::identity::{ExecutionId, SequenceId};
use genome_population_simulator::publication::ExecutionWorkspace;
use genome_population_simulator::results::ArtifactRecord;
use serde_json::{Value, json};
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
        "gps-f5-publication-{}-{}-{label}",
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
    /// Objetivo por bases cubiertas en lugar de por número de bloques. Fuerza
    /// el reparto residual del plan de longitudes, que es la vía que puede
    /// devolver «residual final fuera del soporte factible» **dentro** del
    /// bucle de intentos.
    covered_bases: Option<u64>,
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
            covered_bases: None,
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
                "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {founder_count}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = {}\nmax_swap_proposals = {}\nmax_haplotype_cells = 65536\n{}\n{}{subpopulations}{overrides}",
                sha256(metadata.as_bytes()),
                self.max_haplotype_attempts,
                self.max_swap_proposals,
                self.covered_bases.map_or_else(
                    || format!("block_count = {}", self.block_count),
                    |bases| format!("covered_bases = {bases}"),
                ),
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
// Publicación completa
// ---------------------------------------------------------------------------

/// Un destino publicado, con su árbol final y las piezas que lo produjeron.
struct Package {
    chain: Chain,
    destination: PathBuf,
}

impl Package {
    fn bytes(&self, relative: &str) -> Vec<u8> {
        fs::read(self.destination.join(relative))
            .unwrap_or_else(|error| panic!("artefacto {relative}: {error}"))
    }
    fn document(&self, relative: &str) -> Value {
        serde_json::from_slice(&self.bytes(relative)).expect("documento JSON")
    }
    fn manifest(&self) -> Value {
        self.document(F5_MANIFEST_PATH)
    }
    fn effective(&self) -> Value {
        self.document(F5_EFFECTIVE_CONFIG_PATH)
    }
    /// Catálogo del destino publicado, reconstruido abriendo cada fichero.
    fn records(&self) -> Vec<ArtifactRecord> {
        let mut records = Vec::new();
        for link in links_of(&self.manifest()) {
            let relative = link["path"].as_str().expect("ruta");
            let bytes = self.bytes(relative);
            records.push(ArtifactRecord {
                id: link["artifact_id"].as_str().expect("id").to_owned(),
                kind: link["artifact_kind"].as_str().expect("clase").to_owned(),
                path: relative.to_owned(),
                size_bytes: bytes.len() as u64,
                digest: sha256(&bytes),
            });
        }
        records
    }
}

/// Los enlaces físicos del manifiesto, éxito o fallo: los siete TSV más la
/// verdad, o el único enlace no biológico.
fn links_of(manifest: &Value) -> Vec<Value> {
    let f5 = &manifest["f5"];
    if manifest["state"] == json!("complete") {
        let mut links = f5["artifacts"].as_array().expect("artefactos").clone();
        links.push(f5["truth"].clone());
        links
    } else {
        f5["f5_failure"]["artifacts"]
            .as_array()
            .expect("artefactos")
            .clone()
    }
}

fn publish(label: &str, fixture: &Fixture<'_>) -> Package {
    try_publish(label, fixture).expect("publicación con éxito")
}

fn try_publish(label: &str, fixture: &Fixture<'_>) -> Result<Package, Diagnostic> {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain(&root, fixture);
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let published = {
        let request = F5PublicationRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
        )?;
        publish_f5_success(
            workspace,
            &request,
            &mut chain.identities,
            |payload, declaration| {
                audit_f5_package(payload, declaration)
                    .map(|report| report.initialization_summaries().to_vec())
                    .map_err(|rejection| rejection.diagnostic().clone())
            },
        )
        .map_err(|failure| failure.into_diagnostic())?
    };
    assert_eq!(published, destination);
    Ok(Package { chain, destination })
}

/// Declaración reconstruida **desde el manifiesto publicado**, salvo el perfil.
///
/// Cada cifra que el manifiesto publica —enlaces, verdad, subpoblaciones,
/// muestreo, `F` realizado, cotas y catálogos— se relee de sus bytes y se le
/// entrega al oráculo para que la juzgue contra el paquete. Lo que **no**
/// viene de ahí son los pesos haplotípicos, el espectro, los bins LD y las
/// cuatro auditorías DTC-049: proceden del perfil admitido **en memoria**,
/// porque ningún documento publicado los enumera y `F5NormalizationAudit` no
/// tiene constructor público. No hay circularidad —el oráculo contrasta esos
/// valores contra los bytes del paquete—, pero la independencia de esa parte
/// es menor que la del resto y queda dicho.
fn declaration_from_published(package: &Package) -> F5PackageDeclaration {
    let manifest = package.manifest();
    let f5 = &manifest["f5"];
    let profile = package.chain.config.profile();
    let text = |value: &Value| value.as_str().expect("cadena").to_owned();
    let number = |value: &Value| value.as_u64().expect("entero");
    let optional = |value: &Value| value.as_f64();
    let mode = |value: &Value| match value.as_str() {
        Some("controlled") => F5DeclaredMode::Controlled,
        Some("probabilistic") => F5DeclaredMode::Probabilistic,
        other => panic!("modo inesperado {other:?}"),
    };
    let subpopulations = f5["subpopulations"]
        .as_array()
        .expect("subpoblaciones")
        .iter()
        .map(|row| F5DeclaredSubpopulation {
            population_id: text(&row["population_id"]),
            individual_count: number(&row["individual_count"]),
            inbreeding_f: optional(&row["inbreeding_f"]),
        })
        .collect::<Vec<_>>();
    let ld_bins = subpopulations
        .iter()
        .flat_map(|subpopulation| {
            profile
                .ld_bins(&subpopulation.population_id)
                .expect("bins del perfil")
                .iter()
                .map(|bin| F5DeclaredLdBin {
                    population_id: subpopulation.population_id.clone(),
                    distance_start_bp: bin.distance_start_bp(),
                    distance_end_bp: bin.distance_end_bp(),
                    target_mean_r2: bin.target_mean_r2(),
                    max_abs_error: bin.max_abs_error(),
                })
        })
        .collect();
    let limits = &f5["limits"];
    let preflight = &f5["resource_preflight"];
    F5PackageDeclaration {
        scenario_id: text(&f5["scenario_id"]),
        execution_id: text(&manifest["execution_id"]),
        reference_id: text(&f5["reference_id"]),
        reference_sha256: text(&f5["reference_sha256"]),
        seed: text(&f5["seed"]),
        ploidy: u8::try_from(number(&f5["ploidy"])).expect("ploidía"),
        profile_id: text(&f5["profile"]["profile_id"]),
        profile_metadata_sha256: text(&f5["profile"]["sha256"]),
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
        sampling: f5["sampling"]
            .as_array()
            .expect("muestreo")
            .iter()
            .map(|row| F5DeclaredSampling {
                population_id: text(&row["population_id"]),
                block_id: text(&row["block_id"]),
                requested_mode: mode(&row["requested_mode"]),
                effective_mode: mode(&row["effective_mode"]),
                controlled_max_abs_error_declared: optional(
                    &row["controlled_max_abs_error_declared"],
                ),
                controlled_max_abs_error_effective: optional(
                    &row["controlled_max_abs_error_effective"],
                ),
            })
            .collect(),
        f_realized: f5["f_realized"]
            .as_array()
            .expect("F realizado")
            .iter()
            .map(|row| F5DeclaredRealizedF {
                population_id: text(&row["population_id"]),
                block_id: text(&row["block_id"]),
                target_f: optional(&row["target_f"]),
                status: text(&row["status"]),
                value: optional(&row["value"]),
            })
            .collect(),
        ld_bins,
        normalizations: profile.normalizations().to_vec(),
        artifacts: f5["artifacts"]
            .as_array()
            .expect("artefactos")
            .iter()
            .map(|link| F5DeclaredLink {
                artifact_id: text(&link["artifact_id"]),
                artifact_kind: text(&link["artifact_kind"]),
                schema_version: text(&link["schema_version"]),
                path: text(&link["path"]),
                size_bytes: number(&link["size_bytes"]),
                sha256: text(&link["sha256"]),
            })
            .collect(),
        truth: F5DeclaredTruthLink {
            artifact_id: text(&f5["truth"]["artifact_id"]),
            artifact_kind: text(&f5["truth"]["artifact_kind"]),
            schema_version: text(&f5["truth"]["schema_version"]),
            projection_version: text(&f5["truth"]["projection_version"]),
            path: text(&f5["truth"]["path"]),
            size_bytes: number(&f5["truth"]["size_bytes"]),
            sha256: text(&f5["truth"]["sha256"]),
        },
        bounds: F5DeclaredBounds {
            // B no viaja en ningún documento: se cuenta en las filas del
            // layout publicado, que es la tabla propietaria de los bloques.
            eligible_bases: number(&f5["eligible_bases"]),
            block_count: (String::from_utf8(package.bytes("outputs/f5-block-layout.tsv"))
                .expect("UTF-8")
                .lines()
                .count()
                - 1) as u64,
            site_bound_per_block: number(&preflight["site_bound_per_block"]),
            bin_count_total: number(&preflight["bin_count_total"]),
            max_haplotype_attempts: number(&limits["max_haplotype_attempts"]),
            max_swap_proposals: number(&limits["max_swap_proposals"]),
            tsv_rows_bound: number(&preflight["tsv_rows_bound"]),
            json_rows_bound: number(&preflight["json_rows_bound"]),
            output_bytes: number(&preflight["output_bytes"]),
            output_limit_bytes: number(&preflight["output_limit_bytes"]),
        },
        id_catalog_id: text(&f5["id_catalog"]["catalog_id"]),
        rng_catalog_id: text(&f5["rng_catalog"]["catalog_id"]),
        attempt_summary: F5DeclaredAttemptSummary {
            layout_attempts_used: number(&f5["attempt_summary"]["layout_attempts_used"]),
            haplotype_attempts_total: number(&f5["attempt_summary"]["haplotype_attempts_total"]),
            max_attempt_index: number(&f5["attempt_summary"]["max_attempt_index"]),
        },
    }
}

// ---------------------------------------------------------------------------
// Validador JSON Schema externo
// ---------------------------------------------------------------------------

fn assert_external_cases(schema: &str, dependencies: &[&str], cases: Vec<(Value, bool)>) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let bundle = temporary(&format!("schema-{schema}"));
    fs::create_dir_all(&bundle).expect("paquete aislado de esquemas");
    for name in std::iter::once(schema).chain(dependencies.iter().copied()) {
        assert!(
            !name.contains('/'),
            "el paquete solo admite nombres hermanos"
        );
        fs::copy(root.join("schemas").join(name), bundle.join(name)).expect("copiar esquema");
    }
    let case_path = bundle.join("cases.json");
    let document = cases
        .into_iter()
        .map(|(instance, expected)| json!({"instance": instance, "expected": expected}))
        .collect::<Vec<_>>();
    fs::write(
        &case_path,
        serde_json::to_vec(&document).expect("casos JSON"),
    )
    .expect("escribir casos");
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(bundle.join(schema))
        .arg(&case_path)
        .arg("--cases")
        .stderr(Stdio::piped())
        .output()
        .expect("validador externo");
    fs::remove_dir_all(&bundle).expect("retirar el paquete aislado");
    assert!(
        output.status.success(),
        "{schema}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

// ---------------------------------------------------------------------------
// Éxito: 7/7 + truth
// ---------------------------------------------------------------------------

/// El destino publicado contiene exactamente los siete TSV en el orden
/// normativo, la verdad en su campo separado, los dos documentos raíz y nada
/// más. Cada enlace se contrasta reabriendo el fichero.
#[test]
fn a_successful_publication_links_seven_artifacts_and_the_truth() {
    let package = publish("success", &wide_fixture());
    let manifest = package.manifest();
    assert_eq!(manifest["manifest_kind"], json!("gps_execution_manifest"));
    assert_eq!(manifest["schema_version"], json!("5"));
    assert_eq!(manifest["state"], json!("complete"));
    assert_eq!(manifest["execution_id"], json!(EXECUTION));

    let artifacts = manifest["f5"]["artifacts"].as_array().expect("artefactos");
    assert_eq!(artifacts.len(), 7);
    for (link, contract) in artifacts.iter().zip(&F5_TSV_CONTRACTS) {
        assert_eq!(link["artifact_id"], json!(contract.artifact_id()));
        assert_eq!(link["artifact_kind"], json!(contract.artifact_kind()));
        assert_eq!(link["schema_version"], json!("1"));
        assert_eq!(link["path"], json!(contract.path()));
        let bytes = package.bytes(contract.path());
        assert_eq!(link["size_bytes"], json!(bytes.len() as u64));
        assert_eq!(link["sha256"], json!(sha256(&bytes)));
    }
    // El orden de la tabla es el lexicográfico de `artifact_id`.
    let ids = artifacts
        .iter()
        .map(|link| link["artifact_id"].as_str().expect("id").to_owned())
        .collect::<Vec<_>>();
    let mut sorted = ids.clone();
    sorted.sort();
    assert_eq!(ids, sorted);

    let truth = &manifest["f5"]["truth"];
    assert_eq!(truth["artifact_id"], json!("truth-jsonl"));
    assert_eq!(truth["artifact_kind"], json!("truth_jsonl"));
    assert_eq!(truth["schema_version"], json!("2"));
    assert_eq!(
        truth["projection_version"],
        json!("f5-founder-truth-projection-v2")
    );
    assert_eq!(truth["path"], json!("outputs/truth.jsonl"));
    let bytes = package.bytes("outputs/truth.jsonl");
    assert_eq!(truth["size_bytes"], json!(bytes.len() as u64));
    assert_eq!(truth["sha256"], json!(sha256(&bytes)));
    // La verdad no se duplica dentro de los siete.
    assert!(!ids.iter().any(|id| id == "truth-jsonl"));

    // El árbol publicado contiene exactamente diez ficheros y ningún temporal.
    let mut published = BTreeSet::new();
    collect_files(&package.destination, &package.destination, &mut published);
    let expected = F5_TSV_CONTRACTS
        .iter()
        .map(|contract| contract.path().to_owned())
        .chain([
            "outputs/truth.jsonl".to_owned(),
            F5_MANIFEST_PATH.to_owned(),
            F5_EFFECTIVE_CONFIG_PATH.to_owned(),
        ])
        .collect::<BTreeSet<_>>();
    assert_eq!(published, expected);
}

fn collect_files(root: &Path, directory: &Path, found: &mut BTreeSet<String>) {
    for entry in fs::read_dir(directory).expect("directorio") {
        let path = entry.expect("entrada").path();
        if path.is_dir() {
            collect_files(root, &path, found);
        } else {
            found.insert(
                path.strip_prefix(root)
                    .expect("relativa")
                    .to_string_lossy()
                    .into_owned(),
            );
        }
    }
}

/// El oráculo independiente de F5.4.2 acepta el paquete publicado cuando se le
/// entrega la declaración reconstruida **desde el manifiesto**.
#[test]
fn the_independent_oracle_accepts_the_manifest_it_reads_from_the_package() {
    for (label, fixture) in [
        ("oracle-wide", wide_fixture()),
        ("oracle-sparse", sparse_fixture()),
        ("oracle-dyadic", dyadic_fixture()),
    ] {
        let package = publish(label, &fixture);
        let declaration = declaration_from_published(&package);
        let report =
            audit_f5_package(&package.destination, &declaration).unwrap_or_else(|rejection| {
                panic!("{label}: {}", rejection.diagnostic().render_human())
            });
        // Los resúmenes del manifiesto son exactamente los que el oráculo
        // reconstruye desde `truth.jsonl`, uno a uno y en el mismo orden.
        let published = package.manifest()["f5"]["initialization_summaries"]
            .as_array()
            .expect("resúmenes")
            .clone();
        assert_eq!(published.len(), report.initialization_summaries().len());
        assert!(!published.is_empty());
        for (row, summary) in published.iter().zip(report.initialization_summaries()) {
            assert_eq!(row["population_id"], json!(summary.population_id()));
            assert_eq!(row["site_id"], json!(summary.site_id()));
            assert_eq!(row["model"], json!(summary.model().as_str()));
            assert_eq!(
                row["recount_source"],
                json!("outputs/truth.jsonl#f5-founder-truth-projection-v2")
            );
            assert_eq!(row["recount_matches"], json!(true));
            let recounted = row["recounted_genotype_counts"]
                .as_array()
                .expect("recuento")
                .iter()
                .map(|entry| entry["value"].as_u64().expect("u64"))
                .collect::<Vec<_>>();
            assert_eq!(recounted, summary.recounted_genotype_counts());
            let realized = row["realized_genotype_counts"]
                .as_array()
                .expect("realizado")
                .iter()
                .map(|entry| entry["value"].as_u64().expect("u64"))
                .collect::<Vec<_>>();
            assert_eq!(realized, recounted, "{label}: MAN-70 exige la igualdad");
        }
    }
}

// ---------------------------------------------------------------------------
// Igualdad exacta entre efectivo y manifiesto
// ---------------------------------------------------------------------------

/// Los nueve objetos que DTC-089 obliga a reutilizar y las tres copias de la
/// semilla, releídos de los **dos ficheros del destino**.
#[test]
fn the_effective_configuration_and_the_manifest_publish_the_same_objects() {
    let package = publish("equality", &wide_fixture());
    let effective = package.effective();
    let manifest = package.manifest();
    assert_eq!(
        effective["effective_config_kind"],
        json!("gps_effective_config")
    );
    assert_eq!(effective["schema_version"], json!("5"));
    assert_eq!(effective["source_path"], json!("general.toml"));

    for key in [
        "profile",
        "mask",
        "target",
        "algorithms",
        "rng_catalog",
        "id_catalog",
        "limits",
        "resource_preflight",
        "normalizations",
        "sampling",
        "scenario_id",
        "reference_id",
        "reference_sha256",
        "ploidy",
        "founder_count",
        "schema_version",
    ] {
        assert_eq!(
            effective["f5"][key], manifest["f5"][key],
            "/f5/{key} difiere entre las dos copias"
        );
    }
    // Cada escalar de `resource_preflight`, no solo los tres totales.
    let preflight = effective["f5"]["resource_preflight"]
        .as_object()
        .expect("preflight");
    assert_eq!(preflight.len(), 29);
    for (key, value) in preflight {
        assert_eq!(
            manifest["f5"]["resource_preflight"][key], *value,
            "resource_preflight/{key}"
        );
    }

    // Las tres copias de la semilla.
    let seed = effective["seed"].as_str().expect("semilla");
    assert_eq!(seed.len(), 64);
    assert!(seed.bytes().all(|byte| byte.is_ascii_hexdigit()));
    assert_eq!(effective["f5"]["seed"], json!(seed));
    assert_eq!(manifest["f5"]["seed"], json!(seed));

    // El manifiesto copia cuatro de los seis campos de cada subpoblación.
    let left = effective["f5"]["subpopulations"]
        .as_array()
        .expect("subpoblaciones efectivas");
    let right = manifest["f5"]["subpopulations"]
        .as_array()
        .expect("subpoblaciones del manifiesto");
    assert_eq!(left.len(), right.len());
    for (left, right) in left.iter().zip(right) {
        assert_eq!(left.as_object().expect("objeto").len(), 6);
        assert_eq!(right.as_object().expect("objeto").len(), 4);
        for key in [
            "population_id",
            "individual_count",
            "inbreeding_f",
            "inbreeding_f_source",
        ] {
            assert_eq!(left[key], right[key], "subpopulations/{key}");
        }
    }
    // El manifiesto no repite claves fundadoras ni el vector de individuos.
    assert!(right.iter().all(|row| row.get("founders").is_none()));
}

/// Mutar un término de una copia publicada tumba la puerta física citándolo.
///
/// Es la prueba de discriminación de la igualdad: sin ella, dos documentos
/// distintos se publicarían. Cada caso altera **un** término y exige que el
/// detalle del rechazo lo nombre.
#[test]
fn every_divergence_between_the_two_copies_is_rejected_by_its_own_term() {
    let package = publish("divergence", &wide_fixture());
    let records = package.records();
    let intact = verify_published_f5_package(
        &package.destination,
        &package.destination.join(F5_MANIFEST_PATH),
        &records,
    );
    assert!(intact.is_ok(), "el paquete íntegro debe pasar la puerta");

    let cases: Vec<DivergenceCase> = vec![
        (
            "seed de la raíz efectiva",
            |effective, _| effective["seed"] = json!("f".repeat(64)),
            "seed",
        ),
        (
            "seed del objeto f5 efectivo",
            |effective, _| effective["f5"]["seed"] = json!("e".repeat(64)),
            "seed",
        ),
        (
            "seed del manifiesto",
            |_, manifest| manifest["f5"]["seed"] = json!("d".repeat(64)),
            "seed",
        ),
        (
            "tolerancia de una normalización",
            |effective, _| effective["f5"]["normalizations"][0]["tolerance"] = json!(1e-7),
            "/f5/normalizations",
        ),
        (
            "un valor efectivo de una normalización",
            |_, manifest| {
                manifest["f5"]["normalizations"][3]["effective_values"][0]["value"] = json!(0.4);
            },
            "/f5/normalizations",
        ),
        (
            "un escalar de resource_preflight",
            |_, manifest| {
                manifest["f5"]["resource_preflight"]["tract_count"] = json!(99_u64);
            },
            "/f5/resource_preflight",
        ),
        (
            "un límite",
            |effective, _| effective["f5"]["limits"]["max_swap_proposals"] = json!(7_u64),
            "/f5/limits",
        ),
        (
            "el hash del perfil",
            |_, manifest| manifest["f5"]["profile"]["sha256"] = json!("c".repeat(64)),
            "/f5/profile",
        ),
        (
            "el catálogo de identidades",
            |effective, _| effective["f5"]["id_catalog"]["catalog_id"] = json!("id-catalog-v3"),
            "/f5/id_catalog",
        ),
        (
            "el objetivo",
            |_, manifest| manifest["f5"]["target"]["target_bases"] = json!(1_u64),
            "/f5/target",
        ),
        (
            "un algoritmo",
            |effective, _| effective["f5"]["algorithms"]["ld"] = json!("phased-r2-neumaier-v2"),
            "/f5/algorithms",
        ),
        (
            "la fuente de F de una subpoblación",
            |_, manifest| {
                manifest["f5"]["subpopulations"][0]["inbreeding_f_source"] = json!("subpopulation");
            },
            "/f5/subpopulations/0/inbreeding_f_source",
        ),
        (
            "el censo de una subpoblación",
            |effective, _| effective["f5"]["subpopulations"][0]["individual_count"] = json!(9_u64),
            "/f5/subpopulations/0/individual_count",
        ),
        (
            "una fila de muestreo",
            |_, manifest| manifest["f5"]["sampling"][0]["effective_mode"] = json!("controlled"),
            "/f5/sampling",
        ),
        (
            "la fuente del recuento MAN-70",
            |_, manifest| {
                manifest["f5"]["initialization_summaries"][0]["recount_source"] =
                    json!("outputs/f5-founders.tsv");
            },
            "recount_source",
        ),
        (
            "la igualdad obligatoria del recuento",
            |_, manifest| {
                manifest["f5"]["initialization_summaries"][0]["recount_matches"] = json!(false);
            },
            "recount_matches",
        ),
        (
            "el tamaño de un enlace",
            |_, manifest| manifest["f5"]["artifacts"][2]["size_bytes"] = json!(3_u64),
            "/f5/artifacts/2",
        ),
        (
            "el hash de la verdad",
            |_, manifest| manifest["f5"]["truth"]["sha256"] = json!("a".repeat(64)),
            "/f5/truth",
        ),
        (
            "la versión de proyección de la verdad",
            |_, manifest| {
                manifest["f5"]["truth"]["projection_version"] =
                    json!("f5-founder-truth-projection-v1");
            },
            "/f5/truth",
        ),
        (
            "el orden de los siete enlaces",
            |_, manifest| {
                let artifacts = manifest["f5"]["artifacts"].as_array_mut().expect("lista");
                artifacts.swap(0, 1);
            },
            "/f5/artifacts/0",
        ),
        (
            "el estado del manifiesto",
            |_, manifest| manifest["state"] = json!("failed"),
            "/f5/f5_failure",
        ),
    ];

    // Vigésimo segunda divergencia, y la única que el inventario **no** puede
    // ver: compara rutas, y esta cambia el identificador dejando la ruta
    // intacta. Es el eje propio de la comparación catálogo/manifiesto, que se
    // quedó sin prueba cuando el inventario adelantó a la que la ejercitaba.
    let mut renamed = records.clone();
    renamed[1].id = "f5-02-sitios-v1".to_owned();
    let rejection = verify_published_f5_package(
        &package.destination,
        &package.destination.join(F5_MANIFEST_PATH),
        &renamed,
    )
    .expect_err("un artifact_id ajeno sobre la misma ruta debe tumbar la puerta");
    let detail = rejection
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default();
    assert!(
        detail.contains("f5-02-sitios-v1") && detail.contains("f5-02-sites-v1"),
        "el rechazo debe citar los dos catálogos y dijo {detail:?}"
    );
    assert_eq!(rejection.code(), "GPS-FMT-006");

    for (label, mutate, expected) in cases {
        let staged = temporary(&format!("divergence-{}", label.replace(' ', "-")));
        copy_tree(&package.destination, &staged);
        let mut effective: Value =
            serde_json::from_slice(&fs::read(staged.join(F5_EFFECTIVE_CONFIG_PATH)).expect("leer"))
                .expect("JSON");
        let mut manifest: Value =
            serde_json::from_slice(&fs::read(staged.join(F5_MANIFEST_PATH)).expect("leer"))
                .expect("JSON");
        mutate(&mut effective, &mut manifest);
        fs::write(
            staged.join(F5_EFFECTIVE_CONFIG_PATH),
            serde_json::to_vec(&effective).expect("JSON"),
        )
        .expect("escribir");
        fs::write(
            staged.join(F5_MANIFEST_PATH),
            serde_json::to_vec(&manifest).expect("JSON"),
        )
        .expect("escribir");
        let rejection =
            verify_published_f5_package(&staged, &staged.join(F5_MANIFEST_PATH), &records)
                .expect_err(&format!("«{label}» debe tumbar la puerta"));
        let detail = rejection
            .context()
            .get("detail")
            .cloned()
            .unwrap_or_default();
        assert!(
            detail.contains(expected),
            "«{label}»: se esperaba que el rechazo citara {expected:?} y dijo {detail:?}"
        );
        assert_eq!(rejection.code(), "GPS-FMT-006");
        fs::remove_dir_all(&staged).expect("retirar copia");
    }
}

/// Etiqueta, mutación de una de las dos copias y término que el rechazo debe
/// citar.
type DivergenceCase = (&'static str, fn(&mut Value, &mut Value), &'static str);

/// Etiqueta, mutación de un documento publicado y término que el rechazo debe
/// citar.
type DocumentCase = (&'static str, fn(&mut Value), &'static str);

fn copy_tree(source: &Path, destination: &Path) {
    fs::create_dir_all(destination).expect("destino");
    for entry in fs::read_dir(source).expect("origen") {
        let path = entry.expect("entrada").path();
        let target = destination.join(path.file_name().expect("nombre"));
        if path.is_dir() {
            copy_tree(&path, &target);
        } else {
            fs::copy(&path, &target).expect("copiar");
        }
    }
}

// ---------------------------------------------------------------------------
// Las cifras del preflight publicadas
// ---------------------------------------------------------------------------

/// Los escalares que el hito enumera viajan a los dos documentos con el valor
/// que una cuenta independiente sobre el fixture predice.
#[test]
fn the_published_preflight_carries_the_scalars_the_contract_enumerates() {
    let fixture = wide_fixture();
    let package = publish("preflight", &fixture);
    let effective = package.effective();
    let preflight = &effective["f5"]["resource_preflight"];
    assert_eq!(preflight["algorithm"], json!("f5-resource-accounting-v9"));

    // El FASTA del fixture: dos secuencias pobladas y una estéril, todas de
    // `span` bases. Ninguna cabecera cuenta como base procesada.
    let span = fixture.block_length * (fixture.block_count / fixture.sequences) * fixture.slack;
    let sequences = fixture.sequences + fixture.barren_sequences;
    assert_eq!(
        preflight["reference_bases_processed"],
        json!(span * sequences)
    );
    // Pico de IDs FASTA retenidos: toda secuencia retiene el suyo, poblada o
    // no.
    assert_eq!(preflight["retained_sequence_ids"], json!(sequences));
    // Un tramo elegible por secuencia poblada; la estéril no aporta ninguno.
    assert_eq!(preflight["tract_count"], json!(fixture.sequences));

    // `S_bound = 1 + (Lmax-1)/gap_min`.
    assert_eq!(
        preflight["site_bound_per_block"],
        json!(1 + (fixture.block_length - 1) / fixture.gap)
    );
    // Bins totales: los del perfil por subpoblación.
    assert_eq!(
        preflight["bin_count_total"],
        json!(fixture.bins.len() as u64 * fixture.subpopulations.len() as u64)
    );
    // Filas de perfil y de normalización del paquete canónico.
    assert!(preflight["profile_rows"].as_u64().expect("u64") > 0);
    assert!(preflight["normalization_rows"].as_u64().expect("u64") > 0);
    // Las filas realmente escritas caben en la cota sellada.
    let tsv_rows: u64 = F5_TSV_CONTRACTS
        .iter()
        .map(|contract| {
            String::from_utf8(package.bytes(contract.path()))
                .expect("UTF-8")
                .lines()
                .count() as u64
                - 1
        })
        .sum();
    assert!(tsv_rows <= preflight["tsv_rows_bound"].as_u64().expect("u64"));
    assert!(preflight["json_rows_bound"].as_u64().expect("u64") > 0);

    // Aportaciones ajenas: el flujo F5 las sella a cero y las publica.
    assert_eq!(preflight["other_operation_units"], json!(0_u64));
    assert_eq!(preflight["other_working_memory_bytes"], json!(0_u64));
    assert_eq!(preflight["other_output_bytes"], json!(0_u64));

    // Las tres cotas y sus tres límites duros.
    for (used, limit, hard) in [
        ("operation_units", "operation_limit", 1_u64 << 32),
        (
            "working_memory_bytes",
            "working_memory_limit_bytes",
            1 << 33,
        ),
        ("output_bytes", "output_limit_bytes", 1 << 31),
    ] {
        let used = preflight[used].as_u64().expect("u64");
        let limit = preflight[limit].as_u64().expect("u64");
        assert!(used <= limit, "{used} excede {limit}");
        assert_eq!(limit, hard);
    }
    // Los bytes realmente escritos caben en `output_bytes`.
    let written: u64 = F5_TSV_CONTRACTS
        .iter()
        .map(|contract| package.bytes(contract.path()).len() as u64)
        .sum::<u64>()
        + package.bytes("outputs/truth.jsonl").len() as u64;
    assert!(written <= preflight["output_bytes"].as_u64().expect("u64"));

    // Los tres presupuestos efectivos del TOML son los mismos que `limits`.
    let limits = &effective["f5"]["limits"];
    assert_eq!(limits["max_operation_units"], preflight["operation_limit"]);
    assert_eq!(
        limits["max_working_memory_bytes"],
        preflight["working_memory_limit_bytes"]
    );
    assert_eq!(limits["max_output_bytes"], preflight["output_limit_bytes"]);
    // Y los cuatro límites algorítmicos son los declarados.
    assert_eq!(limits["max_layout_attempts"], json!(8_u64));
    assert_eq!(
        limits["max_haplotype_attempts"],
        json!(fixture.max_haplotype_attempts)
    );
    assert_eq!(
        limits["max_swap_proposals"],
        json!(fixture.max_swap_proposals)
    );
    assert_eq!(limits["max_haplotype_cells"], json!(65_536_u64));

    // Las diez identidades físicas posibles, aquí nueve sin máscara.
    let identities = preflight["input_identities"]
        .as_array()
        .expect("identidades");
    assert_eq!(identities.len(), 9);
    assert_eq!(identities[0]["role"], json!("general_toml"));
    assert_eq!(identities[8]["role"], json!("reference_fasta"));
    assert_eq!(preflight["phases"].as_array().expect("fases").len(), 9);
}

/// Las cuatro auditorías DTC-049 viajan completas: claves canónicas, valores
/// de entrada y efectivos, sumas y deltas. No se sustituyen por un resumen.
#[test]
fn the_four_normalizations_travel_complete_and_in_their_closed_order() {
    let package = publish("normalizations", &dyadic_fixture());
    let normalizations = package.effective()["f5"]["normalizations"]
        .as_array()
        .expect("normalizaciones")
        .clone();
    assert_eq!(normalizations.len(), 4);
    let expected = [
        ("block_lengths", "length_min_bp", 1),
        ("snp_gaps", "gap_min_bp", 1),
        ("snp_spectrum", "alt_af_min", 1),
        ("haplotype_weights", "haplotype_rank", 4),
    ];
    for (audit, (set_kind, order, rows)) in normalizations.iter().zip(expected) {
        assert_eq!(audit["set_kind"], json!(set_kind));
        assert_eq!(audit["canonical_order"], json!(order));
        assert_eq!(audit["algorithm"], json!("neumaier-f64-v1"));
        assert_eq!(audit["tolerance"], json!(1e-9));
        let input = audit["input_values"].as_array().expect("entrada");
        let effective = audit["effective_values"].as_array().expect("efectivos");
        assert_eq!(input.len(), rows);
        assert_eq!(effective.len(), rows);
        for (left, right) in input.iter().zip(effective) {
            assert_eq!(left["key"], right["key"]);
            assert!(left["value"].as_f64().expect("f64").is_finite());
        }
        // La suma de salida cierra en 1 dentro de la tolerancia y el delta
        // publicado es exactamente |suma-1|.
        let output_sum = audit["output_sum"].as_f64().expect("f64");
        let output_delta = audit["output_delta"].as_f64().expect("f64");
        assert!((output_sum - 1.0).abs() <= 1e-9, "{output_sum}");
        assert!((output_delta - (output_sum - 1.0).abs()).abs() <= f64::EPSILON);
    }
    // El perfil diádico publica sus cuatro pesos exactos, no un digest.
    let weights = normalizations[3]["effective_values"]
        .as_array()
        .expect("pesos")
        .iter()
        .map(|entry| entry["value"].as_f64().expect("f64"))
        .collect::<Vec<_>>();
    assert_eq!(weights, vec![0.5, 0.25, 0.125, 0.125]);
}

/// Las longitudes de secuencia que `truth.jsonl` exige nacen del escaneo de
/// admisión, no de una declaración del llamador.
///
/// Es la deuda que F5.4.1 dejó abierta: sin este cableado la vía productiva no
/// puede publicar la proyección. El contraste es con el FASTA del fixture,
/// incluida la secuencia estéril, que la referencia retiene y el layout no
/// puebla.
#[test]
fn the_sealed_admission_carries_every_sequence_length_the_fasta_declares() {
    let fixture = wide_fixture();
    let root = temporary("sequence-lengths");
    fs::create_dir_all(&root).expect("raíz");
    let (general, declared) = fixture.write(&root);
    let config = load_and_validate_f5_config(&general).expect("configuración admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let sequences = F5ReferenceSequences::from_sealed(&sealed).expect("longitudes selladas");
    assert_eq!(sequences.len(), declared.len());
    for (id, length) in &declared {
        let typed = SequenceId::parse(id.clone(), "test").expect("secuencia");
        assert_eq!(sequences.length(&typed), Some(*length), "{id}");
    }
    assert_eq!(
        sequences.lengths().sum::<u64>(),
        sealed.reference_bases_processed()
    );
    fs::remove_dir_all(&root).expect("retirar");
}

// ---------------------------------------------------------------------------
// Ámbitos de configuración
// ---------------------------------------------------------------------------

/// Etiqueta, fixture y, por subpoblación, `population_id`, `inbreeding_f`
/// esperado y su procedencia.
type InbreedingCase<'a> = (&'a str, Fixture<'a>, Vec<(&'a str, Value, &'a str)>);

/// Las dos procedencias de `inbreeding_f` y su resolución por ploidía.
///
/// `subpopulation` cuando la subpoblación lo declara; `population` cuando lo
/// declara el escenario para todas **o cuando nadie lo declara**, porque el
/// valor efectivo procede entonces del ámbito poblacional. Fuera de `P=2` el
/// campo es `null`, nunca cero.
#[test]
fn both_provenances_of_the_inbreeding_coefficient_reach_the_two_documents() {
    let cases: Vec<InbreedingCase<'_>> = vec![
        (
            "f-absent",
            Fixture::default(),
            vec![("P1", json!(0.0), "population")],
        ),
        (
            "f-from-scenario",
            Fixture {
                scenario_inbreeding: Some("1.0"),
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
            vec![("P1", json!(1.0), "population")],
        ),
        (
            "f-from-subpopulation",
            Fixture {
                subpopulations: vec![
                    Subpopulation::new("PA", 2).inbreeding("0.0"),
                    Subpopulation::new("PB", 2).inbreeding("1.0"),
                ],
                ..Fixture::default()
            },
            vec![
                ("PA", json!(0.0), "subpopulation"),
                ("PB", json!(1.0), "subpopulation"),
            ],
        ),
        (
            "f-mixed",
            Fixture {
                scenario_inbreeding: Some("1.0"),
                subpopulations: vec![
                    Subpopulation::new("PA", 2),
                    Subpopulation::new("PB", 2).inbreeding("0.0"),
                ],
                ..Fixture::default()
            },
            vec![
                ("PA", json!(1.0), "population"),
                ("PB", json!(0.0), "subpopulation"),
            ],
        ),
        (
            "f-haploid",
            Fixture {
                ploidy: 1,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
            vec![("P1", Value::Null, "population")],
        ),
        (
            "f-triploid",
            Fixture {
                ploidy: 3,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
            vec![("P1", Value::Null, "population")],
        ),
    ];
    for (label, fixture, expected) in cases {
        let package = publish(label, &fixture);
        let effective = package.effective();
        let manifest = package.manifest();
        let rows = effective["f5"]["subpopulations"]
            .as_array()
            .expect("subpoblaciones");
        assert_eq!(rows.len(), expected.len(), "{label}");
        for (row, (population, value, source)) in rows.iter().zip(&expected) {
            assert_eq!(row["population_id"], json!(population), "{label}");
            assert_eq!(row["inbreeding_f"], *value, "{label}");
            assert_eq!(row["inbreeding_f_source"], json!(source), "{label}");
        }
        // El manifiesto copia los dos campos sin reinterpretarlos.
        for (left, right) in rows
            .iter()
            .zip(manifest["f5"]["subpopulations"].as_array().expect("copia"))
        {
            assert_eq!(left["inbreeding_f"], right["inbreeding_f"], "{label}");
            assert_eq!(
                left["inbreeding_f_source"], right["inbreeding_f_source"],
                "{label}"
            );
        }
        // `f_realized` publica `null`/`not_applicable` fuera de `P=2`.
        let ploidy = manifest["f5"]["ploidy"].as_u64().expect("ploidía");
        for row in manifest["f5"]["f_realized"]
            .as_array()
            .expect("F realizado")
        {
            if ploidy == 2 {
                assert!(!row["target_f"].is_null(), "{label}");
            } else {
                assert_eq!(row["target_f"], Value::Null, "{label}");
                assert_eq!(row["status"], json!("not_applicable"), "{label}");
                assert_eq!(row["value"], Value::Null, "{label}");
            }
        }
    }
}

/// Ambos modos por bloque×subpoblación, IDs aportados y derivados, secuencia
/// estéril y dos subpoblaciones, en un solo paquete publicado.
#[test]
fn the_published_documents_cover_both_sampling_modes_and_both_identity_origins() {
    let package = publish("scopes", &wide_fixture());
    let effective = package.effective();
    let sampling = effective["f5"]["sampling"]
        .as_array()
        .expect("muestreo")
        .clone();
    // Una fila exacta por bloque×subpoblación, ordenada por población y
    // bloque.
    assert_eq!(sampling.len(), 8);
    let keys = sampling
        .iter()
        .map(|row| {
            (
                row["population_id"].as_str().expect("id").to_owned(),
                row["block_id"].as_str().expect("id").to_owned(),
            )
        })
        .collect::<Vec<_>>();
    let mut sorted = keys.clone();
    sorted.sort();
    assert_eq!(keys, sorted);
    assert_eq!(keys.iter().collect::<BTreeSet<_>>().len(), 8);

    let modes = sampling
        .iter()
        .map(|row| row["effective_mode"].as_str().expect("modo").to_owned())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        modes,
        BTreeSet::from(["controlled".to_owned(), "probabilistic".to_owned()])
    );
    let sources = sampling
        .iter()
        .map(|row| row["source"].as_str().expect("origen").to_owned())
        .collect::<BTreeSet<_>>();
    assert_eq!(
        sources,
        BTreeSet::from([
            "block_override".to_owned(),
            "subpopulation_default".to_owned()
        ])
    );
    for row in &sampling {
        match row["effective_mode"].as_str() {
            Some("probabilistic") => {
                assert_eq!(row["controlled_max_abs_error_declared"], Value::Null);
                assert_eq!(row["controlled_max_abs_error_effective"], Value::Null);
            }
            _ => {
                let effective = row["controlled_max_abs_error_effective"]
                    .as_f64()
                    .expect("tolerancia efectiva");
                assert!((0.0..=0.5).contains(&effective), "{effective}");
            }
        }
    }

    // IDs aportados y derivados en la misma publicación: `PA` declara sus dos
    // identidades y `PB` las deriva de claves fundadoras.
    let subpopulations = effective["f5"]["subpopulations"]
        .as_array()
        .expect("subpoblaciones");
    let founders = |population: &str| {
        subpopulations
            .iter()
            .find(|row| row["population_id"] == json!(population))
            .expect("subpoblación")["founders"]
            .as_array()
            .expect("fundadores")
            .clone()
    };
    let supplied = founders("PA");
    assert_eq!(supplied.len(), 2);
    assert!(supplied.iter().all(|row| row["founder_key"].is_null()));
    // Ordenados por `individual_id`, no por el orden del documento.
    assert_eq!(supplied[0]["individual_id"], json!("alfa"));
    assert_eq!(supplied[1]["individual_id"], json!("zeta"));
    let derived = founders("PB");
    assert_eq!(derived.len(), 1);
    assert!(derived[0]["founder_key"].is_string());
    assert!(
        derived[0]["individual_id"]
            .as_str()
            .expect("id")
            .starts_with("IND_")
    );

    // La secuencia estéril se retiene y no aporta ninguna copia fundadora.
    let layout = String::from_utf8(package.bytes("outputs/f5-block-layout.tsv")).expect("UTF-8");
    assert!(!layout.contains("barren1"));
    let founders_tsv = String::from_utf8(package.bytes("outputs/f5-founders.tsv")).expect("UTF-8");
    assert!(!founders_tsv.contains("barren1"));
}

/// Las tres ploidías publican un paquete conforme.
#[test]
fn the_three_ploidies_publish_a_conforming_package() {
    for ploidy in [1_u64, 2, 3] {
        let package = publish(
            &format!("ploidy-{ploidy}"),
            &Fixture {
                ploidy,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
        );
        let manifest = package.manifest();
        assert_eq!(manifest["f5"]["ploidy"], json!(ploidy));
        let declaration = declaration_from_published(&package);
        let report =
            audit_f5_package(&package.destination, &declaration).unwrap_or_else(|rejection| {
                panic!("P={ploidy}: {}", rejection.diagnostic().render_human())
            });
        for summary in report.initialization_summaries() {
            assert_eq!(
                summary.recounted_genotype_counts().len(),
                usize::try_from(ploidy).expect("ploidía") + 1
            );
        }
    }
}

// ---------------------------------------------------------------------------
// Fallo, rechazo, interrupción y atomicidad
// ---------------------------------------------------------------------------

/// Fixture cuyo objetivo LD es imposible: `r²` vale 1 entre los dos sitios de
/// cada bloque y el perfil pide 0 con tolerancia 0,05.
fn exhausting_fixture<'a>() -> Fixture<'a> {
    Fixture {
        block_length: 2,
        block_count: 1,
        bins: vec![(1, 1, 0.0, 0.05)],
        max_haplotype_attempts: 2,
        max_swap_proposals: 4,
        ..Fixture::default()
    }
}

/// Un fallo controlado publica el diagnóstico y **cero** artefactos
/// biológicos, aunque los siete TSV ya estuvieran escritos en el provisional.
#[test]
fn a_controlled_failure_publishes_the_diagnostic_and_no_biological_artifact() {
    let root = temporary("failure");
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, &wide_fixture());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    // Escribe el paquete entero: el fallo llega después y debe retirarlo.
    let mut chain = chain;
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
    let payload = workspace.payload_root().to_path_buf();
    assert!(payload.join("outputs/f5-sites.tsv").is_file());

    let diagnostic = Diagnostic::f5_controlled_quota_infeasible(
        chain.config.scenario_id(),
        "PA",
        "la cuota entera no cabe en la tolerancia efectiva",
    );
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        Some(&chain.layout),
        F5FailureStage::Founders,
        diagnostic,
    );
    assert_eq!(report.code().as_str(), "GPS-DOMAIN-103");
    assert_eq!(report.stage(), F5FailureStage::Founders);
    let published =
        publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
            .expect("fallo publicado");
    assert_eq!(published, destination);

    let mut files = BTreeSet::new();
    collect_files(&destination, &destination, &mut files);
    assert_eq!(
        files,
        BTreeSet::from([F5_MANIFEST_PATH.to_owned(), F5_FAILURE_PATH.to_owned()])
    );
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_MANIFEST_PATH)).expect("leer"))
            .expect("JSON");
    assert_eq!(manifest["state"], json!("failed"));
    let failure = &manifest["f5"]["f5_failure"];
    assert_eq!(failure["code"], json!("GPS-DOMAIN-103"));
    assert_eq!(failure["stage"], json!("founders"));
    assert_eq!(failure["best_layout"], Value::Null);
    assert_eq!(failure["best_attempt"], Value::Null);
    let links = failure["artifacts"].as_array().expect("artefactos");
    assert_eq!(links.len(), 1);
    assert_eq!(links[0]["artifact_id"], json!(F5_FAILURE_ARTIFACT_ID));
    assert_eq!(links[0]["path"], json!(F5_FAILURE_PATH));
    assert!(manifest["f5"].get("truth").is_none());
    assert!(manifest["f5"].get("artifacts").is_none());

    let document: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_FAILURE_PATH)).expect("leer"))
            .expect("JSON");
    assert_eq!(document["diagnostic_kind"], json!("f5_failure"));
    assert_eq!(document["schema_version"], json!("1"));
    assert_eq!(document["artifacts"], json!([]));
    assert_eq!(
        document["counters"]["layout_attempts_used"],
        json!(chain.layout.successful_attempt_index() + 1)
    );
    // Los nueve términos comunes son literalmente los mismos.
    for key in [
        "code",
        "stage",
        "scenario_id",
        "population_id",
        "block_id",
        "limits",
        "counters",
        "best_layout",
        "best_attempt",
    ] {
        assert_eq!(document[key], failure[key], "{key}");
    }
    drop(chain);
}

/// Una búsqueda agotada conserva su mejor intento y sus contadores.
#[test]
fn an_exhausted_search_publishes_its_best_attempt_and_counters() {
    let root = temporary("exhausted");
    fs::create_dir_all(&root).expect("raíz");
    let fixture = exhausting_fixture();
    let (general, _) = fixture.write(&root);
    let config = load_and_validate_f5_config(&general).expect("configuración admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud"),
    )
    .expect("layout");
    let failure: F5SearchFailure = search_f5_haplotypes(&config, &sealed, &execution, &layout)
        .expect_err("la búsqueda debe agotarse");
    assert_eq!(failure.code(), "GPS-DOMAIN-104");
    let report = F5FailureReport::from_search(&config, &layout, &failure);

    let destination = root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    publish_f5_failure(workspace, &config, &sealed, &execution, &report).expect("fallo publicado");
    let document: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_FAILURE_PATH)).expect("leer"))
            .expect("JSON");
    assert_eq!(document["code"], json!("GPS-DOMAIN-104"));
    assert_eq!(document["stage"], json!("haplotypes"));
    assert_eq!(document["population_id"], json!("P1"));
    assert!(document["block_id"].is_string());
    let counters = &document["counters"];
    assert_eq!(counters["layout_attempts_used"], json!(1_u64));
    assert_eq!(
        counters["haplotype_attempts_used"],
        json!(fixture.max_haplotype_attempts)
    );
    assert!(counters["swap_proposals"].as_u64().expect("u64") > 0);
    let best = &document["best_attempt"];
    assert!(!best.is_null());
    assert_eq!(
        best["matrix_bytes_sha256"].as_str().expect("hash").len(),
        64
    );
    assert!(best["max_r2_error"].as_f64().expect("f64") > 0.05);
    assert_eq!(document["best_layout"], Value::Null);
    fs::remove_dir_all(&root).expect("retirar");
}

/// Un rechazo de preflight no abre ni publica destino: devuelve el documento
/// DTC-044 y el árbol final no llega a existir.
#[test]
fn a_preflight_rejection_opens_no_destination() {
    // Dos rechazos de admisión de causas distintas: la geometría demostrada
    // inviable y el presupuesto excedido.
    //
    // El conteo de bloques del caso «over-budget» se escaló ×16 —4096→65536—
    // cuando el techo de operaciones subió de 2^28 a 2^32 el 2026-09-13. Es el
    // mismo factor que subió el techo, de modo que el fixture sigue en el mismo
    // punto relativo: rechaza por `operation_units` en la fase
    // `reference_and_mask_scan`, con 121 916 035 339 frente a 4 294 967 296.
    // Con 4096 ya no excedía nada y la prueba no discriminaba.
    for (label, block_count, corrupt_digest, code, exit_code) in [
        ("over-budget", 65_536_u64, false, "GPS-RES-101", 2_u64),
        ("profile-digest", 2, true, "GPS-INPUT-002", 2),
    ] {
        let root = temporary(&format!("preflight-rejection-{label}"));
        fs::create_dir_all(&root).expect("raíz");
        let fixture = Fixture {
            block_count,
            ..Fixture::default()
        };
        let (general, _) = fixture.write(&root);
        if corrupt_digest {
            let scenario = root.join("authoring/scenarios/scenario.toml");
            let text = fs::read_to_string(&scenario).expect("escenario");
            let start = text.find("profile_metadata_sha256 = \"").expect("campo")
                + "profile_metadata_sha256 = \"".len();
            let mut corrupted = text.clone();
            corrupted.replace_range(start..start + 64, &"0".repeat(64));
            fs::write(&scenario, corrupted).expect("escenario corrupto");
        }
        let destination = root.join("published");
        let diagnostic = load_and_validate_f5_config(&general)
            .err()
            .unwrap_or_else(|| panic!("{label}: la admisión debe rechazar"));
        assert_eq!(diagnostic.code(), code, "{label}");
        let document: Value =
            serde_json::from_str(&diagnostic.render_json_document()).expect("documento DTC-044");
        assert_eq!(document["diagnostic_kind"], json!("gps_diagnostics"));
        assert_eq!(document["schema_version"], json!("1"));
        assert_eq!(document["exit_code"], json!(exit_code), "{label}");
        assert_eq!(document["diagnostics"][0]["code"], json!(code));
        assert_eq!(document["diagnostics"].as_array().expect("lista").len(), 1);
        assert!(
            !destination.exists(),
            "{label}: un rechazo de preflight no publica"
        );
        // Tampoco deja provisional: nadie llegó a crearlo.
        let residue = fs::read_dir(&root)
            .expect("raíz")
            .filter_map(|entry| entry.ok())
            .any(|entry| {
                entry
                    .file_name()
                    .to_string_lossy()
                    .contains("gps-incomplete-")
            });
        assert!(
            !residue,
            "{label}: un rechazo de preflight no abre provisional"
        );
        fs::remove_dir_all(&root).expect("retirar");
    }
}

/// Una ejecución interrumpida deja el provisional marcado `INCOMPLETE` y no
/// crea el destino.
#[test]
fn an_interrupted_execution_leaves_the_provisional_marked_incomplete() {
    let root = temporary("interrupted");
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    {
        let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
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
        // La ejecución se interrumpe aquí: el workspace se destruye sin
        // publicar.
    }
    assert!(!destination.exists(), "una interrupción no publica destino");
    let provisional = fs::read_dir(&root)
        .expect("raíz")
        .filter_map(|entry| entry.ok())
        .map(|entry| entry.path())
        .find(|path| {
            path.file_name()
                .is_some_and(|name| name.to_string_lossy().contains("gps-incomplete-"))
        })
        .expect("el provisional sobrevive a la interrupción");
    assert_eq!(
        fs::read(provisional.join("INCOMPLETE")).expect("marcador"),
        b"incomplete\n"
    );
    // El árbol renombrable nunca contiene el marcador.
    assert!(!provisional.join("payload/INCOMPLETE").exists());
    assert!(provisional.join("payload/outputs/f5-sites.tsv").is_file());
    // Y no queda ningún temporal a medio escribir.
    let mut files = BTreeSet::new();
    collect_files(&provisional, &provisional, &mut files);
    assert!(
        !files.iter().any(|name| name.contains(".gps-tmp-")),
        "{files:?}"
    );
}

/// Si el recuento MAN-70 rechaza el paquete, no se publica nada.
///
/// Es la atomicidad del hito: el manifiesto no puede existir sin su recuento,
/// y el destino no puede existir sin manifiesto.
#[test]
fn a_rejected_recount_publishes_no_destination() {
    let root = temporary("rejected-recount");
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let rejection = {
        let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
        let request = F5PublicationRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
        )
        .expect("solicitud");
        publish_f5_success(workspace, &request, &mut chain.identities, |_, _| {
            Err(Diagnostic::f5_internal("scenario", "recuento rechazado"))
        })
        .expect_err("un recuento rechazado no publica")
    };
    assert_eq!(rejection.diagnostic().code(), "GPS-INTERNAL-101");
    assert!(!destination.exists());
    assert!(!destination.join(F5_MANIFEST_PATH).exists());
}

// ---------------------------------------------------------------------------
// Esquemas externos
// ---------------------------------------------------------------------------

/// Los tres documentos que este hito instancia validan contra sus esquemas
/// publicados con un validador ajeno al serializador Rust.
#[test]
fn the_published_documents_validate_against_their_external_schemas() {
    let package = publish("external-schema", &wide_fixture());
    assert_external_cases(
        "effective-config-v5.schema.json",
        &[],
        vec![(package.effective(), true)],
    );
    assert_external_cases(
        "manifest-general-v5.schema.json",
        &[
            "effective-config-v5.schema.json",
            "f5-failure-v1.schema.json",
        ],
        vec![(package.manifest(), true)],
    );

    let root = temporary("external-failure");
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        None,
        F5FailureStage::Publication,
        Diagnostic::f5_internal(chain.config.scenario_id(), "fallo de publicación"),
    );
    publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
        .expect("fallo publicado");
    let failure: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_FAILURE_PATH)).expect("leer"))
            .expect("JSON");
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_MANIFEST_PATH)).expect("leer"))
            .expect("JSON");
    assert_external_cases("f5-failure-v1.schema.json", &[], vec![(failure, true)]);
    assert_external_cases(
        "manifest-general-v5.schema.json",
        &[
            "effective-config-v5.schema.json",
            "f5-failure-v1.schema.json",
        ],
        vec![(manifest, true)],
    );
}

/// Un documento raíz superviviente tumba el fallo controlado.
///
/// `discard_all_artifacts` vacía el catálogo, pero `effective-config.json` no
/// está en él: se escribe como documento raíz. Un fallo detectado **después**
/// de publicarlo dejaría un paquete fallido que presenta la configuración
/// efectiva de un éxito que no ocurrió, y por eso la puerta física comprueba
/// su ausencia en lugar de confiar en la purga del catálogo.
#[test]
fn a_failure_cannot_publish_a_surviving_root_document() {
    let root = temporary("failure-root-document");
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    workspace
        .write_root_document(F5_EFFECTIVE_CONFIG_PATH, b"{}\n")
        .expect("documento raíz");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        None,
        F5FailureStage::Publication,
        Diagnostic::f5_internal(
            chain.config.scenario_id(),
            "fallo tras escribir el efectivo",
        ),
    );
    let rejection =
        publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
            .expect_err("el fallo no puede publicar la configuración efectiva");
    assert_eq!(rejection.code(), "GPS-FMT-006");
    assert!(
        rejection
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains(F5_EFFECTIVE_CONFIG_PATH)),
        "{:?}",
        rejection.context()
    );
    assert!(!destination.exists());
}

// ---------------------------------------------------------------------------
// Pico de la puerta física
// ---------------------------------------------------------------------------

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

fn peak_bytes<T>(body: impl FnOnce() -> T) -> (T, u64) {
    LIVE.with(|cell| cell.set(0));
    PEAK.with(|cell| cell.set(0));
    TRACKING.with(|cell| cell.set(true));
    let value = body();
    TRACKING.with(|cell| cell.set(false));
    let peak = PEAK.with(Cell::get);
    (value, u64::try_from(peak).unwrap_or_default())
}

/// La puerta física cabe en el pico que la admisión selló y no crece con el
/// número de bloques.
///
/// Es la prueba que faltaba y que la revisión adversarial exigió: ninguna de
/// las tres que vigilaban el pico —escritor, repositorio y oráculo— cubría
/// esta fase, de modo que la puerta podía picar cualquier cosa sin que nada
/// lo viera. Medida antes de corregir, picaba **9,57×** el sellado con `B=32`
/// y **12,59×** con `B=128`, por materializar los dos documentos enteros como
/// árboles `serde_json::Value` —incluido `initialization_summaries`, que
/// crece con `B*N_pop*S`— y por clonar las dos listas de subpoblaciones.
#[test]
fn the_verification_peak_fits_the_sealed_budget_and_does_not_grow_with_the_block_count() {
    let mut measured = Vec::new();
    for blocks in [32_u64, 128] {
        let package = publish(
            &format!("verify-peak-{blocks}"),
            &Fixture {
                block_count: blocks,
                slack: 2,
                ..Fixture::default()
            },
        );
        let sealed_peak = package.chain.sealed.working_memory_bytes();
        let records = package.records();
        let manifest_path = package.destination.join(F5_MANIFEST_PATH);
        let (verdict, peak) = peak_bytes(|| {
            verify_published_f5_package(&package.destination, &manifest_path, &records)
        });
        verdict.expect("el paquete publicado pasa su propia puerta");
        assert!(
            peak <= sealed_peak,
            "B={blocks}: la puerta pica {peak} bytes y la admisión selló {sealed_peak} \
             (razón {:.2})",
            peak as f64 / sealed_peak as f64
        );
        measured.push((blocks, peak));
    }
    let (small_blocks, small_peak) = measured[0];
    let (large_blocks, large_peak) = measured[1];
    // Cuadruplicar `B` no puede cuadruplicar el pico. Lo único que crece son
    // las dos copias de la lista `sampling` —una fila por bloque×subpoblación
    // en cada documento—, que el contrato obliga a comparar entre sí y que
    // cuestan del orden de 2,4 KiB por bloque medidos. Cuatro kibibytes las
    // cubren y siguen estando un orden de magnitud por debajo de los ~55 KiB
    // por bloque que costaba materializar los dos documentos enteros.
    let allowance = (large_blocks - small_blocks) * 4096;
    assert!(
        large_peak <= small_peak + allowance,
        "el pico creció {} bytes al pasar de {small_blocks} a {large_blocks} bloques, \
         más de los {allowance} que las listas por bloque necesitan",
        large_peak.saturating_sub(small_peak)
    );
}

/// La publicación **completa** queda acotada por el doble del pico sellado, y
/// esa desviación se declara.
///
/// Mide `publish_f5_success` entero —escritura de los ocho artefactos,
/// construcción de los dos documentos, la puerta física dos veces y el
/// renombrado— con el recuento MAN-70 ya calculado, que es la frontera que
/// este hito declara: el oráculo queda fuera del flujo y su pico no entra en
/// esta cuenta.
///
/// **Lo que no cabe y por qué.** `f5-resource-accounting-v9` nunca modeló la
/// construcción de los dos documentos: su `final_formula_peak` presupuesta el
/// escritor por flujo de F5.4.1, que cumple con 0,2-0,4×. Lo que sobra aquí
/// es el objeto `f5` del manifiesto más lo que la puerta obliga a
/// materializar para compararlo con el efectivo, ambos `O(B*N_pop)`. Medido:
/// **1,63×** con `B=32` (314 412 frente a 193 120) y **1,36×** con `B=128`
/// (763 759 frente a 562 432). Las tres reducciones que sí
/// se aplicaron llevan la cifra desde **19,08×**: emitir los resúmenes MAN-70
/// por flujo en vez de como árbol `Value` (19,08 → 5,23), serializar el
/// manifiesto directamente al fichero y leerlo por flujo en la puerta
/// (5,23 → 1,83) y soltar el documento antes de juzgarlo (1,83 → 1,63).
/// Bajar de 1 exigiría no materializar `sampling`, `f_realized` ni
/// `quota_ties`, que el contrato obliga a comparar término a término entre
/// los dos documentos. **Queda elevado como bloqueante de F5.6**, que es
/// quien abre la vía productiva.
///
/// La cota de dos veces no es cosmética: el defecto que la revisión encontró
/// picaba 19 veces y esta prueba lo habría matado.
#[test]
fn the_whole_publication_stays_within_twice_the_sealed_budget_without_the_oracle() {
    for blocks in [32_u64, 128] {
        let fixture = Fixture {
            block_count: blocks,
            slack: 2,
            ..Fixture::default()
        };
        // Primera publicación: sirve para obtener los resúmenes ya
        // calculados, fuera de toda medición.
        let reference = publish(&format!("whole-peak-source-{blocks}"), &fixture);
        let declaration = declaration_from_published(&reference);
        let summaries = audit_f5_package(&reference.destination, &declaration)
            .expect("oráculo")
            .initialization_summaries()
            .to_vec();

        let root = temporary(&format!("whole-peak-{blocks}"));
        fs::create_dir_all(&root).expect("raíz");
        let mut chain = build_chain(&root, &fixture);
        let destination = chain.root.join("published");
        let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
        let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
        let sealed_peak = chain.sealed.working_memory_bytes();
        let request = F5PublicationRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
        )
        .expect("solicitud");
        let (published, peak) = peak_bytes(|| {
            publish_f5_success(workspace, &request, &mut chain.identities, |_, _| {
                Ok(summaries.clone())
            })
        });
        published.expect("paquete publicado");
        assert!(
            peak <= 2 * sealed_peak,
            "B={blocks}: la publicación pica {peak} bytes y la admisión selló {sealed_peak} \
             (razón {:.2}, declarada como máximo 2,00)",
            peak as f64 / sealed_peak as f64
        );
    }
}

/// Las cinco etapas publican un fallo que su propio esquema acepta, incluida
/// una causa que **no pertenece** al enumerado cerrado.
///
/// Es la regresión de G-1: `from_diagnostic` tomaba el código del diagnóstico
/// sin filtrar, y `GPS-IO-001` —causa natural de la etapa `publication`— no
/// está entre los diez literales de `f5-failure-v1.schema.json`. El paquete se
/// publicaba atómicamente, con destino final creado, y el validador externo lo
/// rechazaba después.
#[test]
fn every_failure_stage_publishes_a_document_its_own_schema_accepts() {
    let cases: Vec<(F5FailureStage, Diagnostic, &str)> = vec![
        (
            F5FailureStage::Layout,
            Diagnostic::f5_layout_exhausted("scenario", "agotado"),
            "GPS-DOMAIN-106",
        ),
        (
            F5FailureStage::Haplotypes,
            Diagnostic::f5_haplotype_cardinality("scenario", "P1", "BLK", "H imposible"),
            "GPS-DOMAIN-105",
        ),
        (
            F5FailureStage::Founders,
            Diagnostic::f5_controlled_quota_infeasible("scenario", "P1", "cuota"),
            "GPS-DOMAIN-103",
        ),
        (
            F5FailureStage::Preflight,
            Diagnostic::f5_resource("scenario", "presupuesto"),
            "GPS-RES-101",
        ),
        // La causa ajena al enumerado: se publica como `GPS-INTERNAL-101`.
        (
            F5FailureStage::Publication,
            Diagnostic::internal_invariant("f5_publication", "rename").with("detail", "ENOSPC"),
            "GPS-INTERNAL-101",
        ),
    ];
    let mut documents = Vec::new();
    let mut manifests = Vec::new();
    for (stage, diagnostic, expected) in cases {
        let root = temporary(&format!("stage-{}", stage.as_str()));
        fs::create_dir_all(&root).expect("raíz");
        let chain = build_chain(&root, &Fixture::default());
        let destination = chain.root.join("published");
        let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
        let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
        let observed = diagnostic.code();
        let report = F5FailureReport::from_diagnostic(&chain.config, None, stage, diagnostic);
        assert_eq!(report.code().as_str(), expected, "{}", stage.as_str());
        if expected != observed {
            // La causa real no se pierde: viaja en el contexto.
            assert_eq!(report.diagnostic().code(), observed);
            assert_eq!(
                report.diagnostic().context().get("published_failure_code"),
                Some(&"GPS-INTERNAL-101".to_owned())
            );
        }
        publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
            .unwrap_or_else(|error| panic!("{}: {}", stage.as_str(), error.render_human()));
        let document: Value =
            serde_json::from_slice(&fs::read(destination.join(F5_FAILURE_PATH)).expect("leer"))
                .expect("JSON");
        assert_eq!(document["code"], json!(expected));
        assert_eq!(document["stage"], json!(stage.as_str()));
        documents.push((document, true));
        manifests.push((
            serde_json::from_slice::<Value>(
                &fs::read(destination.join(F5_MANIFEST_PATH)).expect("leer"),
            )
            .expect("JSON"),
            true,
        ));
    }
    assert_external_cases("f5-failure-v1.schema.json", &[], documents);
    assert_external_cases(
        "manifest-general-v5.schema.json",
        &[
            "effective-config-v5.schema.json",
            "f5-failure-v1.schema.json",
        ],
        manifests,
    );
    // Y las tres causas naturales de la etapa `publication` siguen fuera del
    // enumerado cerrado, de modo que ninguna puede publicarse tal cual.
    for outside in ["GPS-IO-001", "GPS-INTERNAL-001", "GPS-FMT-006"] {
        assert!(F5FailureCode::parse(outside).is_none(), "{outside}");
    }
}

/// Un nombre de general no publicable se rechaza en la **admisión**, no al
/// serializar el efectivo.
///
/// Es la regresión de G-2: con `mi general.toml` se ejecutaban admisión,
/// layout, búsqueda, fundadores, identidades, los ocho artefactos y el
/// recuento MAN-70 completos antes del rechazo, y el diagnóstico culpaba al
/// esquema del documento producido en lugar de decir qué renombrar.
#[test]
fn a_general_whose_name_is_not_publishable_is_rejected_by_admission() {
    for name in ["mi general.toml", "general (copia).toml", "generación.toml"] {
        let root = temporary("unpublishable-general");
        fs::create_dir_all(&root).expect("raíz");
        let fixture = Fixture::default();
        let (general, _) = fixture.write(&root);
        let renamed = general.with_file_name(name);
        fs::rename(&general, &renamed).expect("renombrar");
        let rejection =
            load_and_validate_f5_config(&renamed).expect_err(&format!("{name} no es publicable"));
        assert_eq!(rejection.code(), "GPS-CONF-101", "{name}");
        assert_eq!(
            rejection.context().get("fields"),
            Some(&"general_path".to_owned()),
            "{name}"
        );
        assert!(
            rejection
                .context()
                .get("detail")
                .is_some_and(|detail| detail.contains("renómbralo")),
            "{name}: el diagnóstico debe decir qué hacer"
        );
        fs::remove_dir_all(&root).expect("retirar");
    }
    // Y el nombre ordinario sigue admitiéndose.
    let root = temporary("publishable-general");
    fs::create_dir_all(&root).expect("raíz");
    let (general, _) = Fixture::default().write(&root);
    load_and_validate_f5_config(&general).expect("general.toml es publicable");
    fs::remove_dir_all(&root).expect("retirar");
}

// ---------------------------------------------------------------------------
// El paquete de fallo, auditado desde sus bytes
// ---------------------------------------------------------------------------

/// Publica un fallo controlado real y devuelve su destino y su cadena.
fn publish_failure(label: &str, stage: F5FailureStage, diagnostic: Diagnostic) -> (Chain, PathBuf) {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let report = F5FailureReport::from_diagnostic(&chain.config, None, stage, diagnostic);
    publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
        .expect("fallo publicado");
    (chain, destination)
}

/// Declaración del paquete de fallo, reconstruida **desde el manifiesto
/// publicado**.
fn failure_declaration_from_published(destination: &Path) -> F5FailureDeclaration {
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_MANIFEST_PATH)).expect("leer"))
            .expect("JSON");
    let failure = &manifest["f5"]["f5_failure"];
    let link = &failure["artifacts"][0];
    let text = |value: &Value| value.as_str().expect("cadena").to_owned();
    let optional = |value: &Value| value.as_str().map(str::to_owned);
    let numbers = |object: &Value, keys: &[&str]| {
        keys.iter()
            .map(|key| object[*key].as_u64().expect("u64"))
            .collect::<Vec<_>>()
    };
    let limits = numbers(
        &failure["limits"],
        &[
            "max_layout_attempts",
            "max_haplotype_attempts",
            "max_swap_proposals",
            "max_haplotype_cells",
            "max_operation_units",
            "max_working_memory_bytes",
            "max_output_bytes",
        ],
    );
    let counters = numbers(
        &failure["counters"],
        &[
            "layout_attempts_used",
            "haplotype_attempts_used",
            "swap_proposals",
            "swap_accepted",
        ],
    );
    F5FailureDeclaration {
        code: text(&failure["code"]),
        stage: text(&failure["stage"]),
        scenario_id: text(&failure["scenario_id"]),
        population_id: optional(&failure["population_id"]),
        block_id: optional(&failure["block_id"]),
        limits: limits.try_into().expect("siete límites"),
        counters: counters.try_into().expect("cuatro contadores"),
        best_layout: failure["best_layout"].clone(),
        best_attempt: failure["best_attempt"].clone(),
        link: F5DeclaredLink {
            artifact_id: text(&link["artifact_id"]),
            artifact_kind: text(&link["artifact_kind"]),
            schema_version: text(&link["schema_version"]),
            path: text(&link["path"]),
            size_bytes: link["size_bytes"].as_u64().expect("u64"),
            sha256: text(&link["sha256"]),
        },
    }
}

/// El oráculo independiente acepta un paquete de fallo íntegro.
#[test]
fn the_independent_oracle_accepts_a_whole_failed_package() {
    let (chain, destination) = publish_failure(
        "failure-oracle",
        F5FailureStage::Founders,
        Diagnostic::f5_controlled_quota_infeasible("scenario", "P1", "cuota inviable"),
    );
    let declaration = failure_declaration_from_published(&destination);
    assert_eq!(declaration.code, "GPS-DOMAIN-103");
    assert_eq!(declaration.stage, "founders");
    audit_f5_failure_package(&destination, &declaration)
        .unwrap_or_else(|rejection| panic!("{}", rejection.diagnostic().render_human()));
    drop(chain);
}

/// Cada corrupción del paquete de fallo muere por su causa.
///
/// Es el eje que la Revisión F5.4 exige y que no existía: la puerta leía
/// `code`, `stage` y el enlace **del manifiesto** y nunca abría
/// `outputs/f5-failure.json`, de modo que un diagnóstico que contradijera a su
/// propio manifiesto en causa, etapa, contadores y artefactos pasaba —incluso
/// rehaciendo el enlace con el tamaño y el SHA-256 nuevos—. La batería aplica
/// la misma disciplina que la del éxito: el paquete íntegro pasa primero,
/// `assert_ne!` sobre los bytes, enlace rehecho y causa contrastada.
#[test]
fn every_byte_mutation_of_a_failed_package_dies_of_its_own_cause() {
    let (chain, source) = publish_failure(
        "failure-mutations",
        F5FailureStage::Founders,
        Diagnostic::f5_controlled_quota_infeasible("scenario", "P1", "cuota inviable"),
    );
    let intact = failure_declaration_from_published(&source);
    audit_f5_failure_package(&source, &intact).expect("el paquete íntegro pasa primero");
    let records = failure_records(&source);
    verify_published_f5_package(&source, &source.join(F5_MANIFEST_PATH), &records)
        .expect("y su puerta física también");

    let cases: Vec<DocumentCase> = vec![
        (
            "causa ajena al manifiesto",
            |d| d["code"] = json!("GPS-DOMAIN-101"),
            "code",
        ),
        (
            "etapa ajena al manifiesto",
            |d| d["stage"] = json!("layout"),
            "stage",
        ),
        (
            "escenario ajeno",
            |d| d["scenario_id"] = json!("otro"),
            "scenario_id",
        ),
        (
            "población inventada",
            |d| d["population_id"] = json!("PZ"),
            "population_id",
        ),
        (
            "bloque inventado",
            |d| d["block_id"] = json!("BLK_X"),
            "block_id",
        ),
        (
            "un límite movido",
            |d| d["limits"]["max_swap_proposals"] = json!(7_u64),
            "limits",
        ),
        (
            "contador de layout inventado",
            |d| d["counters"]["layout_attempts_used"] = json!(999_u64),
            "counters",
        ),
        (
            "mejor layout inventado",
            |d| {
                d["best_layout"] = json!({
                    "attempt_index": 0,
                    "placed_blocks": 1,
                    "unplaced_blocks": 0,
                    "placed_bases": 3,
                    "remaining_target_bases": 0,
                    "layout_bytes_sha256": "a".repeat(64),
                });
            },
            "best_layout",
        ),
        (
            "artefacto biológico en el diagnóstico",
            |d| d["artifacts"] = json!([{"path": "outputs/f5-sites.tsv"}]),
            "artifacts",
        ),
        (
            "clase del diagnóstico",
            |d| d["diagnostic_kind"] = json!("gps_diagnostics"),
            "diagnostic_kind",
        ),
        ("clave de más", |d| d["extra"] = json!(1), "object"),
    ];

    for (label, mutate, expected) in cases {
        let staged = temporary(&format!("failure-mutation-{}", label.replace(' ', "-")));
        copy_tree(&source, &staged);
        let document_path = staged.join(F5_FAILURE_PATH);
        let original = fs::read(&document_path).expect("leer");
        let mut document: Value = serde_json::from_slice(&original).expect("JSON");
        mutate(&mut document);
        let mutated = serde_json::to_vec_pretty(&document).expect("JSON");
        assert_ne!(original, mutated, "«{label}» no tocó ningún byte");
        fs::write(&document_path, &mutated).expect("escribir");

        // El enlace se rehace con los bytes nuevos: sin esto toda mutación
        // moriría por el SHA-256 y taparía su causa.
        let mut declaration = intact.clone();
        declaration.link.size_bytes = mutated.len() as u64;
        declaration.link.sha256 = sha256(&mutated);
        let rejection = audit_f5_failure_package(&staged, &declaration)
            .expect_err(&format!("«{label}» debe morir"));
        let detail = rejection.diagnostic().render_human();
        assert!(
            detail.contains(expected),
            "«{label}»: se esperaba que el rechazo citara {expected:?} y dijo {detail:?}"
        );
        fs::remove_dir_all(&staged).expect("retirar");
    }
    drop(chain);
}

fn failure_records(destination: &Path) -> Vec<ArtifactRecord> {
    let bytes = fs::read(destination.join(F5_FAILURE_PATH)).expect("leer");
    vec![ArtifactRecord {
        id: F5_FAILURE_ARTIFACT_ID.to_owned(),
        kind: "f5_failure".to_owned(),
        path: F5_FAILURE_PATH.to_owned(),
        size_bytes: bytes.len() as u64,
        digest: sha256(&bytes),
    }]
}

/// La puerta física también contrasta el diagnóstico contra el manifiesto.
#[test]
fn the_gate_rejects_a_diagnostic_that_contradicts_its_own_manifest() {
    let (chain, source) = publish_failure(
        "failure-gate",
        F5FailureStage::Founders,
        Diagnostic::f5_controlled_quota_infeasible("scenario", "P1", "cuota inviable"),
    );
    let records = failure_records(&source);
    verify_published_f5_package(&source, &source.join(F5_MANIFEST_PATH), &records)
        .expect("el paquete íntegro pasa");

    let cases: Vec<DocumentCase> = vec![
        ("causa", |d| d["code"] = json!("GPS-DOMAIN-101"), "code"),
        ("etapa", |d| d["stage"] = json!("layout"), "stage"),
        (
            "contadores",
            |d| d["counters"]["layout_attempts_used"] = json!(999_u64),
            "counters",
        ),
    ];
    for (label, mutate, expected) in cases {
        let staged = temporary(&format!("failure-gate-{label}"));
        copy_tree(&source, &staged);
        let document_path = staged.join(F5_FAILURE_PATH);
        let original = fs::read(&document_path).expect("leer");
        let mut document: Value = serde_json::from_slice(&original).expect("JSON");
        mutate(&mut document);
        let mutated = serde_json::to_vec_pretty(&document).expect("JSON");
        assert_ne!(original, mutated, "«{label}» no tocó ningún byte");
        fs::write(&document_path, &mutated).expect("escribir");
        // El enlace del manifiesto se rehace: la mutación no puede morir por
        // el hash.
        let manifest_path = staged.join(F5_MANIFEST_PATH);
        let mut manifest: Value =
            serde_json::from_slice(&fs::read(&manifest_path).expect("leer")).expect("JSON");
        manifest["f5"]["f5_failure"]["artifacts"][0]["size_bytes"] = json!(mutated.len() as u64);
        manifest["f5"]["f5_failure"]["artifacts"][0]["sha256"] = json!(sha256(&mutated));
        fs::write(&manifest_path, serde_json::to_vec(&manifest).expect("JSON")).expect("escribir");
        let staged_records = failure_records(&staged);
        let rejection = verify_published_f5_package(&staged, &manifest_path, &staged_records)
            .expect_err(&format!("«{label}» debe tumbar la puerta"));
        let detail = rejection
            .context()
            .get("detail")
            .cloned()
            .unwrap_or_default();
        assert!(
            detail.contains(expected),
            "«{label}»: se esperaba {expected:?} y dijo {detail:?}"
        );
        fs::remove_dir_all(&staged).expect("retirar");
    }
    drop(chain);
}

/// Un fichero que ningún manifiesto declara **no** llega al destino, en los
/// dos estados.
///
/// Ninguna de las dos ramas de la puerta recorría el árbol: la de éxito no
/// miraba, y la de fallo comprobaba `exists()` sobre nueve rutas literales. Un
/// documento raíz escrito por el llamador —que es la única forma que v5 tiene
/// de conservar la copia TOML que DTC-035 exige, porque `artifacts` es
/// `prefixItems`+`items:false`— se publicaba sin ruta, tamaño ni digest.
#[test]
fn an_undeclared_file_never_reaches_the_destination_in_either_state() {
    // Éxito.
    let root = temporary("orphan-success");
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    workspace
        .write_root_document("orphan.json", b"{}\n")
        .expect("documento huérfano");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let request = F5PublicationRequest::new(
        &chain.config,
        &chain.sealed,
        &execution,
        &chain.layout,
        &chain.search,
        &chain.assignment,
    )
    .expect("solicitud");
    let failure = publish_f5_success(
        workspace,
        &request,
        &mut chain.identities,
        |payload, decl| {
            audit_f5_package(payload, decl)
                .map(|report| report.initialization_summaries().to_vec())
                .map_err(|rejection| rejection.diagnostic().clone())
        },
    )
    .expect_err("un fichero sin declarar no puede publicarse");
    let (_, diagnostic) = failure.into_parts();
    assert_eq!(diagnostic.code(), "GPS-FMT-006");
    assert!(
        diagnostic
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("orphan.json")),
        "{:?}",
        diagnostic.context()
    );
    assert!(!destination.exists());
    drop(chain);

    // Fallo.
    let root = temporary("orphan-failure");
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let mut workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    workspace
        .write_root_document("orphan.json", b"{}\n")
        .expect("documento huérfano");
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        None,
        F5FailureStage::Founders,
        Diagnostic::f5_controlled_quota_infeasible("scenario", "P1", "cuota"),
    );
    let rejection =
        publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
            .expect_err("tampoco en un fallo");
    assert!(
        rejection
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("orphan.json")),
        "{:?}",
        rejection.context()
    );
    assert!(!destination.exists());
    drop(chain);
}

/// Un recuento MAN-70 rechazado es un **fallo controlado** y publica su
/// paquete `failed`.
///
/// DTC-045 reserva el provisional `INCOMPLETE` a la interrupción no
/// controlable y obliga a que un fallo controlado posterior al comienzo
/// «publique atómicamente el directorio final como fallido». Con la firma
/// anterior era imposible: el `ExecutionWorkspace` entraba por valor y al
/// devolver `Err` ya estaba destruido, de modo que el llamador se quedaba sin
/// con qué publicar. Ahora lo recibe de vuelta.
#[test]
fn a_rejected_recount_publishes_a_failed_package_instead_of_an_incomplete_provisional() {
    let root = temporary("recount-failure");
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain(&root, &Fixture::default());
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let failure = {
        let request = F5PublicationRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
        )
        .expect("solicitud");
        publish_f5_success(workspace, &request, &mut chain.identities, |_, _| {
            Err(Diagnostic::f5_internal("scenario", "el recuento rechaza"))
        })
        .expect_err("el recuento rechaza")
    };
    let (workspace, diagnostic) = failure.into_parts();
    let workspace = workspace.expect("el provisional vuelve vivo al llamador");
    assert_eq!(diagnostic.code(), "GPS-INTERNAL-101");

    // El llamador cierra el paquete como fallido, que es lo que DTC-045 pide.
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        Some(&chain.layout),
        F5FailureStage::Publication,
        diagnostic,
    );
    assert_eq!(report.stage(), F5FailureStage::Publication);
    let published =
        publish_f5_failure(workspace, &chain.config, &chain.sealed, &execution, &report)
            .expect("el fallo controlado se publica");
    assert_eq!(published, destination);

    // Destino `failed`, con su diagnóstico y sin artefactos biológicos.
    let mut files = BTreeSet::new();
    collect_files(&destination, &destination, &mut files);
    assert_eq!(
        files,
        BTreeSet::from([F5_MANIFEST_PATH.to_owned(), F5_FAILURE_PATH.to_owned()])
    );
    let manifest: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_MANIFEST_PATH)).expect("leer"))
            .expect("JSON");
    assert_eq!(manifest["state"], json!("failed"));
    assert_eq!(manifest["f5"]["f5_failure"]["stage"], json!("publication"));
    // Y el oráculo independiente lo acepta.
    let declaration = failure_declaration_from_published(&destination);
    audit_f5_failure_package(&destination, &declaration).expect("oráculo del fallo");
    // Ya no queda ningún provisional `INCOMPLETE`: esa marca es de la
    // interrupción, no de un fallo controlado.
    let residue = fs::read_dir(&chain.root)
        .expect("raíz")
        .filter_map(Result::ok)
        .any(|entry| {
            entry
                .file_name()
                .to_string_lossy()
                .contains("gps-incomplete-")
        });
    assert!(!residue, "un fallo controlado no deja provisional");
    drop(chain);
}

/// Los tres términos que no tenían juez ya lo tienen.
#[test]
fn eligible_bases_and_the_attempt_summary_are_now_contrasted() {
    let package = publish("unjudged-terms", &wide_fixture());
    let declaration = declaration_from_published(&package);
    audit_f5_package(&package.destination, &declaration).expect("paquete íntegro");
    let records = package.records();
    let manifest_path = package.destination.join(F5_MANIFEST_PATH);

    // `eligible_bases` contra la cota superior que el propio preflight publica.
    let staged = temporary("unjudged-eligible");
    copy_tree(&package.destination, &staged);
    let mut manifest: Value =
        serde_json::from_slice(&fs::read(staged.join(F5_MANIFEST_PATH)).expect("leer"))
            .expect("JSON");
    manifest["f5"]["eligible_bases"] = json!(999_999_u64);
    fs::write(
        staged.join(F5_MANIFEST_PATH),
        serde_json::to_vec(&manifest).expect("JSON"),
    )
    .expect("escribir");
    let rejection = verify_published_f5_package(&staged, &staged.join(F5_MANIFEST_PATH), &records)
        .expect_err("eligible_bases inventado debe morir");
    assert!(
        rejection
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("/f5/eligible_bases")),
        "{:?}",
        rejection.context()
    );
    fs::remove_dir_all(&staged).expect("retirar");

    // Y `attempt_summary` contra `f5-attempts.tsv`, que el oráculo recorre.
    for (label, mutate, expected) in [
        (
            "total de intentos",
            "haplotype_attempts_total",
            "haplotype_attempts_total",
        ),
        ("índice mayor", "max_attempt_index", "max_attempt_index"),
    ] {
        let mut mutated = declaration.clone();
        match mutate {
            "haplotype_attempts_total" => {
                mutated.attempt_summary.haplotype_attempts_total += 1;
            }
            _ => mutated.attempt_summary.max_attempt_index += 1,
        }
        let rejection = audit_f5_package(&package.destination, &mutated)
            .expect_err(&format!("«{label}» debe morir"));
        let detail = rejection.diagnostic().render_human();
        assert!(
            detail.contains(expected),
            "«{label}»: se esperaba {expected:?} y dijo {detail:?}"
        );
    }
    // El paquete íntegro sigue pasando su puerta.
    verify_published_f5_package(&package.destination, &manifest_path, &records)
        .expect("el íntegro pasa");
}

/// Un layout agotado publica los intentos que **realmente** consumió.
///
/// La regla anterior deducía el contador de la presencia de un mejor
/// candidato, y era falsa: `consider_best` corre antes de `materialize_attempt`
/// en la misma iteración, de modo que un fallo del intento *k* ya tiene
/// candidato, y `From<Diagnostic>` lo descartaba dejando el contador a cero.
/// Un fallo en el intento 3 de 10 publicaba `layout_attempts_used = 0`: un
/// dato falso en un manifiesto publicado. La prueba que lo cubría eran cuatro
/// `assert_eq!` sobre una función pura de dos argumentos, que no podía fallar.
#[test]
fn an_exhausted_layout_publishes_the_attempts_it_really_consumed() {
    let root = temporary("layout-exhausted");
    fs::create_dir_all(&root).expect("raíz");
    // Hueco mayor que el bloque: ningún intento coloca dos sitios, de modo que
    // todos caen por `TooFewSites` y la búsqueda de layout se agota.
    let fixture = Fixture {
        block_length: 3,
        gap: 5,
        block_count: 1,
        bins: vec![(1, 2, 1.0, 0.05)],
        ..Fixture::default()
    };
    let (general, _) = fixture.write(&root);
    let config = load_and_validate_f5_config(&general).expect("configuración admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let failure = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud"),
    )
    .expect_err("el layout debe agotarse");
    assert_eq!(failure.code(), "GPS-DOMAIN-106");
    // Los ocho intentos que el escenario declara, no cero.
    assert_eq!(failure.attempts_used(), sealed.max_layout_attempts());
    assert_eq!(failure.attempts_used(), 8);
    assert!(failure.best_layout().is_some());

    let report = F5FailureReport::from_layout(&config, &failure);
    let destination = root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    publish_f5_failure(workspace, &config, &sealed, &execution, &report)
        .expect("fallo de layout publicado");
    let document: Value =
        serde_json::from_slice(&fs::read(destination.join(F5_FAILURE_PATH)).expect("leer"))
            .expect("JSON");
    assert_eq!(document["stage"], json!("layout"));
    assert_eq!(document["code"], json!("GPS-DOMAIN-106"));
    assert_eq!(document["counters"]["layout_attempts_used"], json!(8_u64));
    assert!(!document["best_layout"].is_null());
    // Y el oráculo del fallo lo acepta.
    let declaration = failure_declaration_from_published(&destination);
    audit_f5_failure_package(&destination, &declaration).expect("oráculo del fallo de layout");
    fs::remove_dir_all(&root).expect("retirar");
}

/// El inventario no admite un enlace simbólico en el árbol publicado.
///
/// Un paquete que declara cada fichero por su huella no puede contener un
/// enlace: el destino queda fuera de la carga útil y su contenido no es el que
/// se hasheó. Hoy es inalcanzable en producción —el provisional es `0o700`, lo
/// crea el propio proceso y ninguna pieza F5 crea enlaces—, pero decidir con
/// `is_dir()`, que **sí** sigue enlaces, dejaba además una recursión sin tope:
/// un enlace a un ancestro la colgaría.
///
/// El caso elegido es un enlace a **directorio**, que es donde `is_dir()` y
/// `file_type()` difieren de verdad. Con la guarda, el rechazo cita el enlace;
/// sin ella, el recorrido lo atraviesa y acaba citando el fichero de dentro,
/// que es un veredicto correcto por la razón equivocada.
#[test]
fn the_payload_inventory_refuses_a_symbolic_link() {
    let package = publish("symlink", &Fixture::default());
    let staged = temporary("symlink-staged");
    copy_tree(&package.destination, &staged);
    let records = package.records();
    verify_published_f5_package(&staged, &staged.join(F5_MANIFEST_PATH), &records)
        .expect("la copia íntegra pasa");

    let aside = temporary("symlink-aside");
    fs::create_dir_all(&aside).expect("directorio ajeno");
    fs::write(aside.join("dentro.txt"), b"x\n").expect("fichero ajeno");
    std::os::unix::fs::symlink(&aside, staged.join("enlace")).expect("enlace a directorio");

    let rejection = verify_published_f5_package(&staged, &staged.join(F5_MANIFEST_PATH), &records)
        .expect_err("un enlace no puede publicarse");
    let detail = rejection
        .context()
        .get("detail")
        .cloned()
        .unwrap_or_default();
    assert!(
        detail.contains("enlace") && !detail.contains("dentro.txt"),
        "el rechazo debe citar el enlace y no lo que hay detrás, y dijo {detail:?}"
    );
    assert_eq!(rejection.code(), "GPS-FMT-006");
    fs::remove_dir_all(&staged).expect("retirar");
    fs::remove_dir_all(&aside).expect("retirar");
}
