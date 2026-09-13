//! F5.5.1. Clausura de identidad: el paquete solo lee lo que F5.3 cerró.
//!
//! ## Qué demuestra esta batería, y qué deja demostrado a otra
//!
//! El hito no produce estado biológico ni toca el paquete atómico: produce
//! **evidencia de una propiedad temporal y de una ausencia**. Por eso ninguna
//! prueba de aquí vuelve a juzgar valores —eso es F5.4.2, con su oráculo
//! independiente— y todas juzgan la relación entre el repositorio que F5.3
//! cerró y los bytes que F5.4.3 publicó:
//!
//! 1. **Precedencia.** El sello se toma antes de publicar; la auditoría, sobre
//!    el destino ya renombrado. Toda identidad publicada tiene que estar en el
//!    repositorio sellado.
//! 2. **Conservación.** Entre el sello y la auditoría el registro solo crece
//!    en una `FHA_` por fila de `f5-founders.tsv`. Se comprueba contra un
//!    recuento de filas hecho con un lector propio de esta batería, no con el
//!    número que la auditoría devuelve.
//! 3. **Mutación posterior prohibida.** El árbol publicado se mide byte a byte
//!    antes y después de auditarlo. El instrumento se valida aparte: una
//!    huella que no distinguiera un byte cambiado no probaría nada.
//!
//! ## Cómo discrimina cada prueba
//!
//! Cada rechazo se provoca mutando **el valor** de un campo publicado —nunca
//! su clave, que era el falso negativo que F5.4.2 pagó— sobre una **copia** del
//! árbol, y toda mutación se hace con una sustitución textual que se
//! **verifica**: si el texto que se iba a sustituir no aparecía, la prueba
//! aborta en vez de pasar sin haber cambiado nada. La causa tipada del rechazo
//! se contrasta una a una: no basta con que haya error.
//!
//! Las dos búsquedas negativas —«ausencia de migración, selección o
//! genealogía» y «la auditoría no puede derivar una identidad»— son **guardas
//! léxicas, no semánticas**: leen el texto del módulo. Por eso ambas llevan
//! control positivo, que es lo único que impide que pasen por estar mal
//! escritas.
//!
//! ## Estudio de mutación sobre la implementación
//!
//! Que una prueba falle contra un paquete corrupto no demuestra que fije la
//! comprobación que anuncia. Cada comprobación de `src/f5_closure.rs` se
//! desactiva por separado y se mide qué prueba muere. Treinta mutaciones,
//! **una sola superviviente** y declarada:
//!
//! | Comprobación desactivada | Prueba que muere |
//! |---|---|
//! | acreditación del destino por el testigo | `a_repository_rebuilt_after_publication_cannot_accredit_the_package` |
//! | identidad de ejecución del sello | `a_witness_from_another_execution_is_rejected` |
//! | copias y unidades del sello | `two_repositories_with_different_censuses_do_not_share_a_witness` |
//! | ley de conservación del éxito | `an_unconsumed_equivalent_repository_breaks_the_conservation_law` |
//! | ley de conservación del fallo | `deriving_identities_after_a_failure_breaks_the_conservation_law` |
//! | ausencia de artefacto biológico en el fallo | `a_failed_package_publishes_no_founder_identity` |
//! | cabecera exacta de `f5-founders.tsv` | `the_exact_header_of_the_founders_file_is_required` |
//! | resolución de la copia en el repositorio | `a_copy_the_repository_never_closed_is_rejected` |
//! | pertenencia en `f5-founders.tsv` | `a_copy_moved_to_another_subpopulation_is_rejected` |
//! | individuo contra la copia cerrada | `a_founder_renamed_to_another_closed_founder_is_rejected` |
//! | tupla `sequence_id`/`copy_index` | `a_copy_published_on_another_sequence_is_rejected` y `..._with_another_copy_index_...` |
//! | clase del identificador | las cuatro pruebas de clase `FCP_`, `FHA_`, `HAP_` y la de la verdad |
//! | unidad bloque×subpoblación cerrada | `a_unit_the_repository_never_closed_is_rejected` |
//! | `HAP_` dentro de su unidad cerrada | `a_haplotype_the_repository_never_closed_is_rejected` |
//! | orden de las copias en el TSV | `a_copy_repeated_in_the_founders_file_is_rejected` |
//! | avance de bloques dentro de una copia | `a_block_repeated_within_a_copy_is_rejected` |
//! | columnas de menos | `a_row_with_a_missing_column_is_rejected` |
//! | columnas de más | `a_row_with_an_extra_column_is_rejected` |
//! | LF final | `a_founders_file_without_its_final_lf_is_rejected` |
//! | conjunto de copias publicadas | `a_closed_copy_missing_from_the_package_is_rejected` |
//! | documento de verdad vacío | `an_empty_truth_document_is_rejected` |
//! | cabecera de `truth.jsonl` | `a_truth_document_without_its_header_is_rejected` |
//! | `founder_membership` ausente | `a_copy_record_without_membership_is_rejected` |
//! | `identity` ausente | `a_copy_record_without_identity_is_rejected` |
//! | copia de la verdad no cerrada | `a_truth_copy_the_repository_never_closed_is_rejected` |
//! | pertenencia en `truth.jsonl` | `a_membership_record_moved_to_another_subpopulation_is_rejected` |
//! | identidad en `truth.jsonl` | `a_copy_record_with_another_identity_is_rejected` |
//! | orden de la pertenencia | `a_duplicated_membership_record_is_rejected` |
//! | recuento de pertenencias | `a_missing_membership_record_is_rejected` |
//! | corte en el primer genotipo | ocho pruebas, entre ellas la de ámbitos y la de pico |
//!
//! La única superviviente declarada es `check_repository_order`: desactivarla
//! no mata ninguna prueba porque **ninguna configuración admitida la alcanza**.
//! No es una comprobación sobre el paquete sino una precondición de las dos
//! búsquedas binarias de la propia auditoría, y se conserva por eso. Queda
//! escrito en vez de disimulado.
//!
//! El mismo estudio, en sus dos rondas, corrigió cuatro defectos de esta
//! batería: cuatro pruebas que pasaban por una comprobación distinta de la que
//! anunciaban —de ahí que el orden se juzgue **después** de las comparaciones
//! de identidad y que cada prueba fije campo y detalle—, y cuatro aserciones
//! que no podían fallar, retiradas.

use std::alloc::{GlobalAlloc, Layout, System};
use std::cell::Cell;
use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::f5_artifacts::{F5_TRUTH_PATH, F5_TSV_CONTRACTS};
use genome_population_simulator::f5_closure::{
    F5ClosureCause, F5ClosureRejection, F5ClosureReport, F5ClosureWitness,
    audit_f5_failure_closure, audit_f5_identity_closure, seal_and_publish_f5_failure,
    seal_and_publish_f5_success,
};
use genome_population_simulator::f5_founders::{
    F5FounderAssignment, F5FounderIdentities, assign_f5_founders, materialize_f5_founder_identities,
};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_oracle::audit_f5_package;
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::f5_publication::{
    F5FailureReport, F5FailureStage, F5PublicationRequest,
};
use genome_population_simulator::f5_search::{F5HaplotypeSearch, search_f5_haplotypes};
use genome_population_simulator::identity::ExecutionId;
use genome_population_simulator::publication::ExecutionWorkspace;
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
/// Otra ejecución, para que un sello ajeno sea reconocible por su identidad y
/// no solo por sus recuentos.
const OTHER_EXECUTION: &str = "RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";
const CONTINUOUS_SPECTRUM: &str = "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n";

/// Posición contratada de `f5-founders.tsv` dentro de los siete TSV.
const FOUNDERS: usize = 5;

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-closure-{}-{}-{label}",
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
    root: PathBuf,
}

impl Drop for Chain {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

fn build_chain(root: &Path, fixture: &Fixture<'_>) -> Chain {
    build_chain_as(root, fixture, EXECUTION)
}

fn build_chain_as(root: &Path, fixture: &Fixture<'_>, execution_id: &str) -> Chain {
    let (general, _lengths) = fixture.write(root);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(execution_id, "test").expect("execution");
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
    Chain {
        config,
        sealed,
        layout,
        search,
        assignment,
        identities,
        root: root.to_path_buf(),
    }
}

// ---------------------------------------------------------------------------
// Publicación, con el sello tomado **antes** de publicar
// ---------------------------------------------------------------------------

/// Un destino publicado con el sello que precede a su publicación.
struct Package {
    chain: Chain,
    destination: PathBuf,
    witness: F5ClosureWitness,
}

impl Package {
    fn audit(&self) -> Result<F5ClosureReport, F5ClosureRejection> {
        audit_f5_identity_closure(
            &self.destination,
            &self.witness,
            &self.chain.identities,
            &self.chain.assignment,
        )
    }
    /// Audita el destino ya publicado **después** de alterarlo en su sitio.
    ///
    /// La corrupción se aplica sobre el paquete real, no sobre una copia: el
    /// testigo acredita un destino concreto y auditar otro directorio saldría
    /// por «el testigo no acredita este paquete» antes de llegar a lo que la
    /// prueba mide. Alterar el destino publicado es además, literalmente, la
    /// mutación posterior que este hito prohíbe.
    fn audit_after(&self, tamper: impl FnOnce(&Path)) -> F5ClosureRejection {
        tamper(&self.destination);
        self.audit()
            .expect_err("el paquete alterado no puede cerrar la clausura")
    }
}

/// Fixture ancho: dos subpoblaciones, identidades aportadas y derivadas, dos
/// secuencias pobladas, una estéril y los dos modos de muestreo.
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

fn publish(label: &str, fixture: &Fixture<'_>) -> Package {
    publish_as(label, fixture, EXECUTION)
}

fn publish_as(label: &str, fixture: &Fixture<'_>, execution_id: &str) -> Package {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("raíz");
    let mut chain = build_chain_as(&root, fixture, execution_id);
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(execution_id, "test").expect("ejecución");
    // `seal_and_publish_f5_success` sella y **luego** publica: el testigo que
    // devuelve es lo único que acredita ese orden, y no hay otra forma de
    // obtenerlo.
    let witness = {
        let request = F5PublicationRequest::new(
            &chain.config,
            &chain.sealed,
            &execution,
            &chain.layout,
            &chain.search,
            &chain.assignment,
        )
        .expect("solicitud de publicación");
        seal_and_publish_f5_success(
            workspace,
            &request,
            &mut chain.identities,
            |payload, declaration| {
                audit_f5_package(payload, declaration)
                    .map(|report| report.initialization_summaries().to_vec())
                    .map_err(|rejection| rejection.diagnostic().clone())
            },
        )
        .map_err(|failure| failure.into_diagnostic())
        .expect("publicación con éxito")
    };
    assert_eq!(witness.destination(), destination);
    Package {
        chain,
        destination,
        witness,
    }
}

// ---------------------------------------------------------------------------
// Instrumentos: lectura independiente, copia del árbol y mutación verificada
// ---------------------------------------------------------------------------

/// Filas de `f5-founders.tsv`, leídas con un lector propio de esta batería.
///
/// El recuento que contrasta el informe **no** puede salir del informe.
fn founders_rows(root: &Path) -> Vec<Vec<String>> {
    let path = root.join(F5_TSV_CONTRACTS[FOUNDERS].path());
    let text = String::from_utf8(fs::read(&path).expect("leer fundadores")).expect("UTF-8");
    let mut lines = text.lines();
    assert_eq!(
        lines.next(),
        Some(F5_TSV_CONTRACTS[FOUNDERS].header()),
        "la cabecera publicada no es la del contrato"
    );
    lines
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn truth_lines(root: &Path) -> Vec<String> {
    let text =
        String::from_utf8(fs::read(root.join(F5_TRUTH_PATH)).expect("leer verdad")).expect("UTF-8");
    text.lines().map(str::to_owned).collect()
}

/// Reescribe `f5-founders.tsv` con su cabecera contratada y las filas dadas.
fn write_founders(root: &Path, rows: &[String]) {
    let mut text = String::from(F5_TSV_CONTRACTS[FOUNDERS].header());
    text.push('\n');
    for row in rows {
        text.push_str(row);
        text.push('\n');
    }
    fs::write(root.join(F5_TSV_CONTRACTS[FOUNDERS].path()), text).expect("escribir fundadores");
}

fn write_truth(root: &Path, lines: &[String]) {
    let mut text = lines.join("\n");
    text.push('\n');
    fs::write(root.join(F5_TRUTH_PATH), text).expect("escribir verdad");
}

/// Índices de las líneas `haplotype_copy` de `truth.jsonl`.
fn copy_record_indices(lines: &[String]) -> Vec<usize> {
    lines
        .iter()
        .enumerate()
        .filter(|(_, line)| line.contains(r#""record_type":"haplotype_copy""#))
        .map(|(index, _)| index)
        .collect()
}

/// Copia recursiva del árbol publicado, para no mutar el destino original.
fn clone_tree(source: &Path, label: &str) -> PathBuf {
    let target = temporary(label);
    copy_into(source, &target);
    target
}

fn copy_into(source: &Path, target: &Path) {
    fs::create_dir_all(target).expect("directorio");
    for entry in fs::read_dir(source).expect("leer directorio") {
        let entry = entry.expect("entrada");
        let kind = entry.file_type().expect("tipo");
        let destination = target.join(entry.file_name());
        if kind.is_dir() {
            copy_into(&entry.path(), &destination);
        } else {
            assert!(kind.is_file(), "el árbol publicado solo tiene ficheros");
            fs::copy(entry.path(), &destination).expect("copiar");
        }
    }
}

/// Sustituye el valor de una columna de `f5-founders.tsv` y devuelve el
/// anterior.
///
/// Muta **el valor**, nunca la clave ni la cabecera: el falso negativo que
/// F5.4.2 pagó fue exactamente una mutación que corrompía la clave y medía,
/// sin saberlo, otro contrato.
fn rewrite_founder_field(root: &Path, row: usize, column: usize, value: &str) -> String {
    let path = root.join(F5_TSV_CONTRACTS[FOUNDERS].path());
    let text = String::from_utf8(fs::read(&path).expect("leer")).expect("UTF-8");
    let mut lines = text.lines().map(str::to_owned).collect::<Vec<_>>();
    let index = row + 1;
    assert!(index < lines.len(), "la fila {row} no existe");
    let mut fields = lines[index]
        .split('\t')
        .map(str::to_owned)
        .collect::<Vec<_>>();
    assert_eq!(fields.len(), 8, "la fila no tiene las ocho columnas");
    let previous = std::mem::replace(&mut fields[column], value.to_owned());
    assert_ne!(
        previous, value,
        "la mutación no cambia nada: la prueba no discriminaría"
    );
    lines[index] = fields.join("\t");
    let mut out = lines.join("\n");
    out.push('\n');
    fs::write(&path, out).expect("escribir");
    previous
}

/// Fragmento literal que el escritor emite para la pertenencia de una copia.
fn membership_fragment(founder_copy_id: &str, population_id: &str) -> String {
    format!(
        r#""founder_membership":{{"founder_copy_id":"{founder_copy_id}","population_id":"{population_id}"}}"#
    )
}

/// Sustitución textual **verificada**: si el texto buscado no aparecía el
/// número esperado de veces, la prueba aborta en lugar de pasar sin mutar.
fn substitute(root: &Path, relative: &str, from: &str, to: &str, expected: usize) {
    let path = root.join(relative);
    let text = String::from_utf8(fs::read(&path).expect("leer")).expect("UTF-8");
    assert_eq!(
        text.matches(from).count(),
        expected,
        "la mutación de {relative} no encontró {expected} apariciones de {from}"
    );
    fs::write(&path, text.replace(from, to)).expect("escribir");
}

/// Contrasta un rechazo contra **el contrato concreto** que debía romperse.
///
/// No basta con la causa tipada: varias comprobaciones comparten causa, y un
/// estudio de mutación sobre la propia implementación demostró que cuatro de
/// estas pruebas pasaban por una comprobación distinta de la que anunciaban.
/// Por eso se fija también el campo y un fragmento del detalle.
fn assert_rejection(
    rejection: &F5ClosureRejection,
    cause: F5ClosureCause,
    field: &str,
    fragment: &str,
) {
    assert_eq!(rejection.cause(), cause, "causa");
    let context = rejection.diagnostic().context();
    assert_eq!(
        context.get("field").map(String::as_str),
        Some(field),
        "campo"
    );
    let detail = context.get("detail").expect("todo rechazo lleva detalle");
    assert!(
        detail.contains(fragment),
        "el detalle {detail:?} no es el de la comprobación anunciada ({fragment:?})"
    );
}

/// Sustitución verificada **dentro de una sola línea** de un JSONL.
///
/// Un mismo literal aparece en varios registros —`"individual_id":"alfa"` sale
/// una vez por copia de ese fundador—, de modo que una sustitución global
/// mutaría más de lo que la prueba anuncia. Se exige exactamente una aparición
/// en la línea elegida.
fn substitute_in_record(root: &Path, relative: &str, index: usize, from: &str, to: &str) {
    let path = root.join(relative);
    let text = String::from_utf8(fs::read(&path).expect("leer")).expect("UTF-8");
    let mut lines = text.lines().map(str::to_owned).collect::<Vec<_>>();
    assert!(index < lines.len(), "el registro {index} no existe");
    assert_eq!(
        lines[index].matches(from).count(),
        1,
        "la mutación esperaba una sola aparición de {from} en el registro {index}"
    );
    lines[index] = lines[index].replace(from, to);
    let mut out = lines.join("\n");
    out.push('\n');
    fs::write(&path, out).expect("escribir");
}

/// Huella del árbol completo: ruta relativa, tamaño y SHA-256 de cada fichero.
fn tree_fingerprint(root: &Path) -> BTreeMap<String, (u64, String)> {
    let mut fingerprint = BTreeMap::new();
    collect_fingerprint(root, root, &mut fingerprint);
    fingerprint
}

fn collect_fingerprint(root: &Path, directory: &Path, into: &mut BTreeMap<String, (u64, String)>) {
    for entry in fs::read_dir(directory).expect("leer directorio") {
        let entry = entry.expect("entrada");
        let kind = entry.file_type().expect("tipo");
        if kind.is_dir() {
            collect_fingerprint(root, &entry.path(), into);
        } else {
            let bytes = fs::read(entry.path()).expect("leer fichero");
            let relative = entry
                .path()
                .strip_prefix(root)
                .expect("relativa")
                .to_string_lossy()
                .into_owned();
            into.insert(relative, (bytes.len() as u64, sha256(&bytes)));
        }
    }
}

// ---------------------------------------------------------------------------
// 1. Precedencia, atestiguada
// ---------------------------------------------------------------------------

/// El paquete publicado nombra **exactamente** lo que F5.3 cerró, y el
/// registro creció en una `FHA_` por fila.
#[test]
fn the_published_package_only_names_identities_closed_in_f5_3() {
    let package = publish("closed", &wide_fixture());
    let report = package.audit().expect("clausura auditada");

    // Recuento independiente: lo hace esta batería sobre los bytes, no el
    // informe que se está juzgando.
    let rows = founders_rows(&package.destination);
    let mut distinct_copies = rows.iter().map(|row| row[5].clone()).collect::<Vec<_>>();
    distinct_copies.sort_unstable();
    distinct_copies.dedup();
    let mut distinct_individuals = rows
        .iter()
        .map(|row| (row[0].clone(), row[1].clone()))
        .collect::<Vec<_>>();
    distinct_individuals.sort_unstable();
    distinct_individuals.dedup();

    assert_eq!(report.assignments(), rows.len() as u64);
    assert_eq!(report.copies(), distinct_copies.len() as u64);
    assert_eq!(report.individuals(), distinct_individuals.len() as u64);

    // La ley de conservación, escrita aquí con los dos extremos medidos aparte.
    assert_eq!(
        package.chain.identities.registered_identities(),
        package.witness.seal().registered_identities() + rows.len(),
        "el registro tiene que crecer en exactamente una FHA_ por fila publicada"
    );

    // Identidades aportadas y derivadas: `PA` declara dos `founder_ids`
    // literales y `PB` una `founder_key`, que deriva un `IND_`.
    assert_eq!(report.supplied_individuals(), 2);
    assert_eq!(report.derived_individuals(), 1);
    let derived = distinct_individuals
        .iter()
        .filter(|(_, individual)| individual.starts_with("IND_"))
        .count();
    assert_eq!(
        derived as u64,
        report.derived_individuals(),
        "los derivados publicados son los que llevan el prefijo de clase"
    );
}

/// **La sonda de la revisión adversarial, invertida.**
///
/// Con un sello suelto, rematerializar un repositorio equivalente **después**
/// de que el paquete existiera, sellarlo entonces y consumir su recorrido
/// bastaba para avalar el paquete viejo: las cuatro cifras del sello son una
/// función determinista del repositorio sin consumir. El testigo lo cierra —
/// publicar el repositorio rematerializado exige un destino nuevo, y el
/// testigo que devuelve acredita ese destino, no el paquete que se quería
/// avalar.
#[test]
fn a_repository_rebuilt_after_publication_cannot_accredit_the_package() {
    let package = publish("precedence", &wide_fixture());
    // Control positivo: el testigo verdadero sí acredita su paquete.
    package.audit().expect("clausura auditada");

    // El repositorio rematerializado desde cero, con la misma ejecución, y
    // publicado —porque no hay otra forma de obtener un testigo— a su propio
    // destino, ya existiendo el paquete que se quiere avalar.
    let rebuilt = publish("precedence-rebuilt", &wide_fixture());
    assert_eq!(
        rebuilt.witness.seal().registered_identities(),
        package.witness.seal().registered_identities(),
        "las cifras del sello son idénticas: por eso un sello suelto no bastaba"
    );
    assert_eq!(
        rebuilt.witness.seal().execution_id(),
        package.witness.seal().execution_id()
    );

    let rejection = audit_f5_identity_closure(
        &package.destination,
        &rebuilt.witness,
        &rebuilt.chain.identities,
        &rebuilt.chain.assignment,
    )
    .expect_err("un testigo de otro paquete no acredita este");
    assert_rejection(
        &rejection,
        F5ClosureCause::Seal,
        "destination",
        "el testigo acredita",
    );
}

/// Un repositorio equivalente **sin consumir** rompe la ley de conservación.
///
/// Es la otra mitad de la sonda: aunque el testigo sea el verdadero, el
/// repositorio que se presenta tiene que ser el que la publicación consumió.
/// Uno recién materializado tiene las mismas cifras de sello pero no ha
/// derivado ninguna `FHA_`.
#[test]
fn an_unconsumed_equivalent_repository_breaks_the_conservation_law() {
    let package = publish("conservation", &wide_fixture());
    let root = temporary("conservation-unconsumed");
    fs::create_dir_all(&root).expect("raíz");
    let fresh = build_chain(&root, &wide_fixture());
    assert_eq!(
        fresh.identities.registered_identities(),
        package.witness.seal().registered_identities(),
        "sin consumir, el repositorio nuevo registra lo mismo que el sello"
    );

    let rejection = audit_f5_identity_closure(
        &package.destination,
        &package.witness,
        &fresh.identities,
        &fresh.assignment,
    )
    .expect_err("un repositorio sin consumir no publicó estas asignaciones");
    assert_rejection(
        &rejection,
        F5ClosureCause::Conservation,
        "registered_identities",
        "no derivó exactamente una FHA_ por fila",
    );
}

/// Un testigo de otra ejecución no avala este repositorio.
#[test]
fn a_witness_from_another_execution_is_rejected() {
    let package = publish("own-execution", &Fixture::default());
    let other = publish_as("other-execution", &Fixture::default(), OTHER_EXECUTION);
    assert_ne!(
        package.witness.seal().execution_id(),
        other.witness.seal().execution_id(),
        "los dos testigos tienen que ser de ejecuciones distintas"
    );

    // Se audita el destino que el testigo ajeno sí acredita, para que lo que
    // falle sea la identidad de ejecución y no la acreditación del destino.
    let rejection = audit_f5_identity_closure(
        &other.destination,
        &other.witness,
        &package.chain.identities,
        &package.chain.assignment,
    )
    .expect_err("un testigo ajeno no avala este repositorio");
    assert_rejection(
        &rejection,
        F5ClosureCause::Seal,
        "execution_id",
        "el sello es de la ejecución",
    );
}

/// Dos repositorios de censo distinto no comparten testigo.
///
/// Es la rama `copies,units` de `check_seal`, que la revisión señaló como
/// alcanzable y sin prueba.
#[test]
fn two_repositories_with_different_censuses_do_not_share_a_witness() {
    let small = publish(
        "census-small",
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 1)],
            ..Fixture::default()
        },
    );
    let large_root = temporary("census-large");
    fs::create_dir_all(&large_root).expect("raíz");
    let large = build_chain(
        &large_root,
        &Fixture {
            subpopulations: vec![Subpopulation::new("P1", 3)],
            ..Fixture::default()
        },
    );
    assert_ne!(
        large.identities.copies().len(),
        small.witness.seal().copies(),
        "los dos censos tienen que diferir"
    );

    let rejection = audit_f5_identity_closure(
        &small.destination,
        &small.witness,
        &large.identities,
        &large.assignment,
    )
    .expect_err("el testigo no selló este repositorio");
    assert_rejection(
        &rejection,
        F5ClosureCause::Seal,
        "copies,units",
        "el sello conserva",
    );
}

// ---------------------------------------------------------------------------
// 2. Identidad, pertenencia y censo: cada mutación con su causa
// ---------------------------------------------------------------------------

/// Una copia que el repositorio nunca cerró.
#[test]
fn a_copy_the_repository_never_closed_is_rejected() {
    let package = publish("unclosed-copy", &wide_fixture());
    let rejection = package.audit_after(|root| {
        // Otro `FCP_` sintácticamente válido que ninguna derivación produjo.
        let previous = rewrite_founder_field(root, 0, 5, "FCP_AAAAAAAAAAAAAAAAAAAAAAAAAA");
        assert!(previous.starts_with("FCP_"));
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "founder_copy_id",
        "no es ninguna de las copias que F5.3 cerró",
    );
}

/// Una copia cerrada publicada en otra subpoblación.
#[test]
fn a_copy_moved_to_another_subpopulation_is_rejected() {
    let package = publish("moved-copy", &wide_fixture());
    let rejection = package.audit_after(|root| {
        // La última fila es de `PB`; se la publica como si fuera de `PA`.
        let rows = founders_rows(root);
        let last = rows.len() - 1;
        assert_eq!(rows[last][0], "PB");
        assert_eq!(rewrite_founder_field(root, last, 0, "PA"), "PB");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "population_id",
        "pertenece a PB y el paquete la publica en PA",
    );
}

/// Un fundador renombrado **a otro fundador que sí está en el censo**.
///
/// Es la variante que discrimina: con un nombre ajeno al censo el rechazo lo
/// produciría la búsqueda en el censo, no la comparación contra la copia
/// cerrada, y la prueba mediría otra cosa. `alfa` y `zeta` son los dos
/// `founder_ids` literales de `PA`.
#[test]
fn a_founder_renamed_to_another_closed_founder_is_rejected() {
    let package = publish("renamed-founder", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        assert_eq!(
            rows[0][1], "alfa",
            "la primera fila es la del fundador alfa"
        );
        assert_eq!(rewrite_founder_field(root, 0, 1, "zeta"), "alfa");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Roster,
        "individual_id",
        "es del fundador alfa y el paquete la publica como zeta",
    );
}

/// Una copia publicada sobre otra secuencia de la misma referencia.
///
/// Es la tupla `sequence_id`/`copy_index`, que el plan nombra como «su tupla
/// tiene que ser la cerrada» y que no tenía prueba.
#[test]
fn a_copy_published_on_another_sequence_is_rejected() {
    let package = publish("moved-sequence", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        assert_eq!(rows[0][2], "chr1");
        assert_eq!(rewrite_founder_field(root, 0, 2, "chr2"), "chr1");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "sequence_id,copy_index",
        "y el paquete la publica sobre chr2/0",
    );
}

/// Una copia publicada con otro índice de copia.
#[test]
fn a_copy_published_with_another_copy_index_is_rejected() {
    let package = publish("moved-index", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        assert_eq!(rows[0][3], "0");
        assert_eq!(rewrite_founder_field(root, 0, 3, "1"), "0");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "sequence_id,copy_index",
        "y el paquete la publica sobre chr1/1",
    );
}

/// Una fila haplotípica que el repositorio nunca cerró para esa unidad.
#[test]
fn a_haplotype_the_repository_never_closed_is_rejected() {
    let package = publish("unclosed-haplotype", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let previous = rewrite_founder_field(root, 0, 7, "HAP_AAAAAAAAAAAAAAAAAAAAAAAAAA");
        assert!(previous.starts_with("HAP_"));
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "haplotype_id",
        "no es ninguna de las filas haplotípicas que F5.3 cerró",
    );
}

/// `FHA_` no sustituye a `FCP_`: DTC-089 lo dice y la clausura lo comprueba.
#[test]
fn an_assignment_identifier_cannot_take_the_place_of_a_copy() {
    let package = publish("substitution-fha", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let copy_id = rows[0][5].clone();
        assert!(rewrite_founder_field(root, 0, 6, &copy_id).starts_with("FHA_"));
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Substitution,
        "block_assignment_id",
        "no es un identificador de clase FHA_",
    );
}

/// Y `FCP_` tampoco ocupa el lugar de otra clase.
#[test]
fn a_copy_identifier_of_another_class_is_rejected() {
    let package = publish("substitution-fcp", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let assignment_id = rows[0][6].clone();
        assert!(rewrite_founder_field(root, 0, 5, &assignment_id).starts_with("FCP_"));
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Substitution,
        "founder_copy_id",
        "no es un identificador de clase FCP_",
    );
}

/// Y `HAP_` tampoco.
#[test]
fn a_haplotype_identifier_of_another_class_is_rejected() {
    let package = publish("substitution-hap", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let copy_id = rows[0][5].clone();
        assert!(rewrite_founder_field(root, 0, 7, &copy_id).starts_with("HAP_"));
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Substitution,
        "haplotype_id",
        "no es un identificador de clase HAP_",
    );
}

/// La cabecera exacta es el **plano de bytes** de «ausencia de migración,
/// selección o genealogía»: la mitad que no depende de una guarda léxica.
#[test]
fn the_exact_header_of_the_founders_file_is_required() {
    let package = publish("header", &Fixture::default());
    let rejection = package.audit_after(|root| {
        let path = root.join(F5_TSV_CONTRACTS[FOUNDERS].path());
        let text = String::from_utf8(fs::read(&path).expect("leer")).expect("UTF-8");
        // Una columna de genealogía añadida a la cabecera contratada.
        let (header, body) = text.split_once('\n').expect("cabecera");
        fs::write(&path, format!("{header}\tparent_id\n{body}")).expect("escribir");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "header",
        "no es exactamente la de su contrato",
    );
}

/// Una copia cerrada que el paquete no publica: pertenencia ausente en el TSV.
#[test]
fn a_closed_copy_missing_from_the_package_is_rejected() {
    let package = publish("missing-copy", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let orphan = rows[rows.len() - 1][5].clone();
        let kept = rows
            .iter()
            .filter(|row| row[5] != orphan)
            .map(|row| row.join("\t"))
            .collect::<Vec<_>>();
        assert!(
            kept.len() < rows.len(),
            "la mutación tiene que retirar al menos una fila"
        );
        write_founders(root, &kept);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_copy_id",
        "copias y el paquete publica",
    );
}

/// Una copia publicada dos veces en el TSV.
#[test]
fn a_copy_repeated_in_the_founders_file_is_rejected() {
    let package = publish("repeated-copy", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        // La primera fila se repite al final: su clave vuelve hacia atrás.
        let mut lines = rows.iter().map(|row| row.join("\t")).collect::<Vec<_>>();
        lines.push(rows[0].join("\t"));
        write_founders(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_copy_id",
        "se publica fuera de su orden o más de una vez",
    );
}

/// Dentro de una copia, los bloques tienen que avanzar.
#[test]
fn a_block_repeated_within_a_copy_is_rejected() {
    let package = publish("repeated-block", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        assert_eq!(
            rows[0][5], rows[1][5],
            "las dos primeras filas son bloques de la misma copia"
        );
        let mut lines = rows.iter().map(|row| row.join("\t")).collect::<Vec<_>>();
        // El segundo bloque repite el `block_id` del primero.
        lines[1] = rows[0].join("\t");
        write_founders(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "block_id",
        "los bloques de una copia no avanzan",
    );
}

/// Sin LF final la última fila está truncada, en los dos recorridos.
#[test]
fn a_founders_file_without_its_final_lf_is_rejected() {
    let package = publish("founders-eol", &Fixture::default());
    let rejection = package.audit_after(|root| {
        let path = root.join(F5_TSV_CONTRACTS[FOUNDERS].path());
        let mut bytes = fs::read(&path).expect("leer");
        assert_eq!(bytes.pop(), Some(b'\n'), "el contrato exige LF final");
        fs::write(&path, &bytes).expect("escribir");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "eol",
        "falta el LF final",
    );
}

/// Una fila con una columna de menos.
#[test]
fn a_row_with_a_missing_column_is_rejected() {
    let package = publish("short-row", &Fixture::default());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let mut lines = rows.iter().map(|row| row.join("\t")).collect::<Vec<_>>();
        lines[0] = rows[0][..7].join("\t");
        write_founders(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "columns",
        "la fila tiene 7 columnas y el contrato fija 8",
    );
}

/// Una fila con una columna de más.
#[test]
fn a_row_with_an_extra_column_is_rejected() {
    let package = publish("long-row", &Fixture::default());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let mut lines = rows.iter().map(|row| row.join("\t")).collect::<Vec<_>>();
        lines[0] = format!("{}\tPADRE", rows[0].join("\t"));
        write_founders(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "columns",
        "más de 8 columnas",
    );
}

/// La pertenencia publicada en `truth.jsonl` no es la cerrada.
#[test]
fn a_membership_record_moved_to_another_subpopulation_is_rejected() {
    let package = publish("membership-moved", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let lines = truth_lines(root);
        let indices = copy_record_indices(&lines);
        let last = *indices.last().expect("alguna copia");
        let rows = founders_rows(root);
        let copy_id = rows[rows.len() - 1][5].clone();
        let from = membership_fragment(&copy_id, "PB");
        assert!(lines[last].contains(&from), "la última copia es la de PB");
        substitute(
            root,
            F5_TRUTH_PATH,
            &from,
            &membership_fragment(&copy_id, "PA"),
            1,
        );
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_membership.population_id",
        "pertenece a PB y la verdad la publica en PA",
    );
}

/// Pertenencia duplicada.
#[test]
fn a_duplicated_membership_record_is_rejected() {
    let package = publish("membership-duplicated", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let mut lines = truth_lines(root);
        let indices = copy_record_indices(&lines);
        let first = indices[0];
        lines.insert(first + 1, lines[first].clone());
        write_truth(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_membership",
        "fuera de orden o más de una vez",
    );
}

/// Pertenencia ausente en el recuento.
#[test]
fn a_missing_membership_record_is_rejected() {
    let package = publish("membership-missing", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let mut lines = truth_lines(root);
        let indices = copy_record_indices(&lines);
        let last = *indices.last().expect("alguna copia");
        lines.remove(last);
        write_truth(root, &lines);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_membership",
        "la verdad publica la pertenencia de",
    );
}

/// Una copia publicada **sin** su objeto `founder_membership`.
///
/// Es «pertenencia ausente» en el registro, no en el recuento, y el plan la
/// nombra como evidencia del hito.
#[test]
fn a_copy_record_without_membership_is_rejected() {
    let package = publish("membership-absent", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let lines = truth_lines(root);
        let first = copy_record_indices(&lines)[0];
        let rows = founders_rows(root);
        let from = format!("{},", membership_fragment(&rows[0][5], &rows[0][0]));
        assert!(lines[first].contains(&from));
        substitute(root, F5_TRUTH_PATH, &from, "", 1);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Membership,
        "founder_membership",
        "no declara su pertenencia",
    );
}

/// Una copia publicada **sin** su objeto `identity`.
#[test]
fn a_copy_record_without_identity_is_rejected() {
    let package = publish("identity-absent", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let lines = truth_lines(root);
        let first = copy_record_indices(&lines)[0];
        let rows = founders_rows(root);
        let from = format!(
            r#""identity":{{"copy_index":{},"individual_id":"{}","#,
            rows[0][3], rows[0][1]
        );
        assert!(
            lines[first].contains(&from),
            "el prefijo de identity cambió"
        );
        substitute_in_record(root, F5_TRUTH_PATH, first, &from, r#""identity_absent":{"#);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "identity",
        "no declara su identidad",
    );
}

/// La identidad publicada en `truth.jsonl` no es la que sostiene esa copia.
#[test]
fn a_copy_record_with_another_identity_is_rejected() {
    let package = publish("identity-other", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        assert_eq!(rows[0][1], "alfa");
        let first = copy_record_indices(&truth_lines(root))[0];
        substitute_in_record(
            root,
            F5_TRUTH_PATH,
            first,
            r#""individual_id":"alfa""#,
            r#""individual_id":"zeta""#,
        );
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "identity",
        "y la verdad la publica sobre zeta/",
    );
}

/// En `truth.jsonl` la pertenencia también nombra su propia clase.
#[test]
fn a_membership_identifier_of_another_class_is_rejected() {
    let package = publish("membership-class", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let copy_id = rows[0][5].clone();
        substitute(
            root,
            F5_TRUTH_PATH,
            &format!(r#""founder_copy_id":"{copy_id}""#),
            r#""founder_copy_id":"FHA_AAAAAAAAAAAAAAAAAAAAAAAAAA""#,
            1,
        );
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Substitution,
        "founder_copy_id",
        "no es un identificador de clase FCP_",
    );
}

/// Un bloque×subpoblación que el repositorio nunca cerró.
#[test]
fn a_unit_the_repository_never_closed_is_rejected() {
    let package = publish("unclosed-unit", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let block = rows[0][4].clone();
        // Se le cambia el último carácter: sigue siendo un `BLK_` sintáctico y
        // no es ninguno de los bloques publicados.
        let mut absent = block.clone();
        let last = absent.pop().expect("identificador no vacío");
        absent.push(if last == 'A' { 'B' } else { 'A' });
        assert!(
            rows.iter().all(|row| row[4] != absent),
            "el bloque inventado no puede existir en el paquete"
        );
        assert_eq!(rewrite_founder_field(root, 0, 4, &absent), block);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "block_id",
        "no es ninguna de las unidades que F5.3 cerró",
    );
}

/// La verdad nombra una copia que el repositorio nunca cerró.
///
/// Se distingue de la sustitución de clase: el identificador es un `FCP_`
/// sintáctico y de su propia clase, solo que ninguna derivación lo produjo.
#[test]
fn a_truth_copy_the_repository_never_closed_is_rejected() {
    let package = publish("truth-unclosed", &wide_fixture());
    let rejection = package.audit_after(|root| {
        let rows = founders_rows(root);
        let copy_id = rows[0][5].clone();
        let first = copy_record_indices(&truth_lines(root))[0];
        substitute_in_record(
            root,
            F5_TRUTH_PATH,
            first,
            &format!(r#""founder_copy_id":"{copy_id}""#),
            r#""founder_copy_id":"FCP_AAAAAAAAAAAAAAAAAAAAAAAAAA""#,
        );
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "founder_copy_id",
        "no es ninguna de las copias que F5.3 cerró",
    );
}

/// Un documento de verdad vacío.
#[test]
fn an_empty_truth_document_is_rejected() {
    let package = publish("truth-empty", &Fixture::default());
    let rejection = package.audit_after(|root| {
        fs::write(root.join(F5_TRUTH_PATH), "").expect("vaciar");
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "header",
        "documento vacío",
    );
}

/// Un documento de verdad **sin cabecera**.
///
/// Sin esta comprobación la primera copia se consumiría como cabecera y el
/// rechazo saldría con una causa falsa: exactamente el defecto que el otro
/// recorrido ya había corregido.
#[test]
fn a_truth_document_without_its_header_is_rejected() {
    let package = publish("truth-headless", &Fixture::default());
    let rejection = package.audit_after(|root| {
        let lines = truth_lines(root);
        assert!(lines[0].contains(r#""record_type":"header""#));
        write_truth(root, &lines[1..]);
    });
    assert_rejection(
        &rejection,
        F5ClosureCause::ByteContract,
        "header",
        "no es la cabecera v2",
    );
}

// ---------------------------------------------------------------------------
// 3. Mutación posterior prohibida
// ---------------------------------------------------------------------------

/// La auditoría no escribe un solo byte del paquete publicado.
///
/// El instrumento se valida en la misma prueba: una huella que no distinguiera
/// un byte cambiado no probaría nada.
#[test]
fn the_closure_audit_never_writes_a_byte_of_the_published_package() {
    let package = publish("immutable", &wide_fixture());
    let before = tree_fingerprint(&package.destination);
    package.audit().expect("clausura auditada");
    // Y otra vez: la auditoría es idempotente y sigue sin escribir.
    package.audit().expect("clausura auditada de nuevo");
    let after = tree_fingerprint(&package.destination);
    assert_eq!(before, after, "la auditoría alteró el paquete publicado");

    // Discriminación del instrumento, sobre una copia.
    let corrupt = clone_tree(&package.destination, "immutable-instrument");
    let path = corrupt.join(F5_TSV_CONTRACTS[FOUNDERS].path());
    let mut bytes = fs::read(&path).expect("leer");
    let last = bytes.len() - 2;
    bytes[last] ^= 0x20;
    fs::write(&path, &bytes).expect("escribir");
    assert_ne!(
        before,
        tree_fingerprint(&corrupt),
        "la huella no distingue un byte cambiado y no podría probar nada"
    );
    fs::remove_dir_all(&corrupt).expect("retirar");
}

/// Después de publicar, ninguna vía puede derivar otra identidad.
///
/// El recorrido de asignaciones es de un solo uso: el segundo intento muere
/// **sin** registrar nada, y la clausura sigue cerrando.
#[test]
fn no_identity_can_be_derived_after_publication() {
    let mut package = publish("single-use", &wide_fixture());
    let registered = package.chain.identities.registered_identities();

    let error = package
        .chain
        .identities
        .for_each_assignment(&package.chain.assignment, |_| Ok(()))
        .expect_err("derivar un FHA_ dos veces es una colisión");
    assert_eq!(error.code(), "GPS-INTERNAL-101");
    assert_eq!(
        package.chain.identities.registered_identities(),
        registered,
        "un intento rechazado no puede registrar nada"
    );
    package
        .audit()
        .expect("la clausura sigue cerrando tras el intento rechazado");
}

// ---------------------------------------------------------------------------
// 4. El paquete de fallo
// ---------------------------------------------------------------------------

/// Publica un fallo controlado real y devuelve su cadena, su destino y el
/// testigo que acredita el orden sello→publicación.
fn publish_failure(label: &str, fixture: &Fixture<'_>) -> (Chain, PathBuf, F5ClosureWitness) {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("raíz");
    let chain = build_chain(&root, fixture);
    let destination = chain.root.join("published");
    let workspace = ExecutionWorkspace::create(&destination).expect("provisional");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("ejecución");
    let report = F5FailureReport::from_diagnostic(
        &chain.config,
        None,
        F5FailureStage::Founders,
        Diagnostic::f5_controlled_quota_infeasible(
            chain.config.scenario_id(),
            "P1",
            "la cuota entera no cabe en la tolerancia efectiva",
        ),
    );
    let witness = seal_and_publish_f5_failure(
        workspace,
        &chain.config,
        &chain.sealed,
        &execution,
        &report,
        &chain.identities,
    )
    .expect("fallo publicado");
    assert_eq!(witness.destination(), destination);
    (chain, destination, witness)
}

/// Un fallo controlado no publica ninguna identidad fundadora.
#[test]
fn a_failed_package_publishes_no_founder_identity() {
    let (chain, destination, witness) = publish_failure("failure", &Fixture::default());
    let closure = audit_f5_failure_closure(&destination, &witness, &chain.identities)
        .expect("un fallo no expone identidad");
    // Lo que se comprueba no es que un `default()` sea `default()`, sino que
    // el informe de fallo **no** se construya con cifras de un paquete de
    // éxito: el mismo fixture, publicado con éxito, sí publica copias,
    // asignaciones e individuos.
    let success = publish("failure-contrast", &Fixture::default());
    let successful = success.audit().expect("clausura auditada");
    assert!(
        successful.copies() > 0 && successful.assignments() > 0 && successful.individuals() > 0,
        "el éxito del mismo fixture sí publica identidad"
    );
    assert_eq!(
        (
            closure.copies(),
            closure.assignments(),
            closure.individuals()
        ),
        (0, 0, 0),
        "un fallo no publica ninguna copia, asignación ni individuo"
    );

    // Discriminación: un artefacto biológico plantado en el fallo se rechaza.
    let rejection = {
        fs::create_dir_all(destination.join("outputs")).expect("outputs");
        fs::write(
            destination.join(F5_TSV_CONTRACTS[FOUNDERS].path()),
            format!("{}\n", F5_TSV_CONTRACTS[FOUNDERS].header()),
        )
        .expect("plantar");
        audit_f5_failure_closure(&destination, &witness, &chain.identities)
            .expect_err("un fallo con artefacto biológico expondría identidad")
    };
    assert_rejection(
        &rejection,
        F5ClosureCause::Unclosed,
        "path",
        "un paquete fallido publica un artefacto biológico",
    );
}

/// **La mitad «fallo» de la ley de conservación.**
///
/// Publicar un fallo no deriva ninguna identidad. La revisión adversarial
/// demostró que la versión anterior no lo comprobaba: solo miraba que el
/// registro no encogiese —y el registro nunca encoge—, de modo que derivar el
/// recorrido completo **después** de publicar el fallo pasaba la auditoría.
#[test]
fn deriving_identities_after_a_failure_breaks_the_conservation_law() {
    let (mut chain, destination, witness) = publish_failure("failure-derive", &Fixture::default());
    let sealed_count = witness.seal().registered_identities();
    assert_eq!(chain.identities.registered_identities(), sealed_count);
    // Control positivo: antes de derivar nada, la clausura del fallo cierra.
    audit_f5_failure_closure(&destination, &witness, &chain.identities)
        .expect("un fallo no expone identidad");

    let mut derived = 0_usize;
    chain
        .identities
        .for_each_assignment(&chain.assignment, |_| {
            derived += 1;
            Ok(())
        })
        .expect("el recorrido todavía no se había consumido");
    assert!(derived > 0, "el recorrido tiene que derivar algo");
    assert_eq!(
        chain.identities.registered_identities(),
        sealed_count + derived
    );

    let rejection = audit_f5_failure_closure(&destination, &witness, &chain.identities)
        .expect_err("publicar un fallo no puede derivar ninguna identidad");
    assert_rejection(
        &rejection,
        F5ClosureCause::Conservation,
        "registered_identities",
        "publicar un fallo no deriva ninguna",
    );
}

// ---------------------------------------------------------------------------
// 5. Los ámbitos enteros
// ---------------------------------------------------------------------------

/// Ploidías 1, 2 y 3, los dos modos, `F` en sus extremos, identidades
/// aportadas y derivadas, secuencia estéril y dos subpoblaciones.
#[test]
fn every_ploidy_mode_and_f_extreme_keeps_the_closure() {
    let cases: Vec<(&str, Fixture<'_>, u64, u64)> = vec![
        (
            "haploid",
            Fixture {
                ploidy: 1,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
            2,
            2,
        ),
        (
            "triploid",
            Fixture {
                ploidy: 3,
                subpopulations: vec![Subpopulation::new("P1", 2)],
                ..Fixture::default()
            },
            2,
            6,
        ),
        (
            "f-zero",
            Fixture {
                subpopulations: vec![Subpopulation::new("P1", 2).inbreeding("0.0")],
                ..Fixture::default()
            },
            2,
            4,
        ),
        (
            "f-one",
            Fixture {
                subpopulations: vec![Subpopulation::new("P1", 2).inbreeding("1.0")],
                ..Fixture::default()
            },
            2,
            4,
        ),
        (
            "controlled",
            Fixture {
                subpopulations: vec![Subpopulation::new("P1", 2).sampling("controlled")],
                ..Fixture::default()
            },
            2,
            4,
        ),
        ("wide", wide_fixture(), 3, 12),
    ];
    for (label, fixture, individuals, copies) in cases {
        let package = publish(label, &fixture);
        let report = package
            .audit()
            .unwrap_or_else(|rejection| panic!("{label}: {rejection:?}"));
        let rows = founders_rows(&package.destination);
        assert_eq!(report.assignments(), rows.len() as u64, "{label}: filas");
        assert_eq!(report.individuals(), individuals, "{label}: individuos");
        // `C_total` se calcula aquí desde el fixture, no desde el informe: son
        // `individual_count * ploidy` copias por secuencia poblada.
        assert_eq!(report.copies(), copies, "{label}: copias");
        assert_eq!(
            package.chain.identities.registered_identities(),
            package.witness.seal().registered_identities() + rows.len(),
            "{label}: conservación"
        );
    }
}

// ---------------------------------------------------------------------------
// 6. Búsquedas negativas, declaradas como léxicas
// ---------------------------------------------------------------------------

/// F5 no simula migración, selección ni genealogía, y la clausura tampoco las
/// nombra.
///
/// **Es una guarda léxica, no semántica**: lee el texto del módulo. Su valor
/// entero está en el control positivo, que impide que pase por estar mal
/// escrita.
#[test]
fn the_closure_module_names_no_migration_selection_or_genealogy() {
    let source = include_str!("../src/f5_closure.rs");
    for forbidden in [
        "migration",
        "migrant",
        "selection",
        "fitness",
        "pedigree",
        "genealogy",
        "lineage",
        "ancestor",
        "offspring",
        "parent_id",
        "recombination",
    ] {
        assert!(
            !source.contains(forbidden),
            "F5 no simula {forbidden} y la clausura no puede nombrarlo"
        );
    }
    // Control positivo: la búsqueda solo vale si el módulo sí nombra lo suyo.
    assert!(source.contains("founder_membership"));
    assert!(source.contains("founder_copy_id"));
}

/// La auditoría no puede derivar una identidad, y no por promesa: por firma.
///
/// Derivar exige `&mut F5FounderIdentities`; la clausura solo recibe `&`. La
/// guarda es **léxica** y la refuerza el comprobador de préstamos, que es
/// quien realmente lo impide.
#[test]
fn the_closure_audit_cannot_derive_an_identity_by_its_signature() {
    let source = include_str!("../src/f5_closure.rs");
    // Las dos auditorías reciben el repositorio en préstamo **compartido**.
    for signature in [
        "pub fn audit_f5_identity_closure(",
        "pub fn audit_f5_failure_closure(",
    ] {
        let body = source.split(signature).nth(1).expect("la firma existe");
        let parameters = body.split(')').next().expect("los parámetros");
        assert!(
            parameters.contains("identities: &F5FounderIdentities"),
            "{signature} tiene que recibir el repositorio en préstamo compartido"
        );
        assert!(
            !parameters.contains("&mut"),
            "{signature} no puede recibir nada en préstamo exclusivo"
        );
    }
    // El único préstamo exclusivo del módulo es el que la publicación exige, y
    // lo atraviesa sin tocarlo.
    assert_eq!(
        source.matches("&mut F5FounderIdentities").count(),
        1,
        "el único préstamo exclusivo es el que seal_and_publish_f5_success delega"
    );
    assert!(source.contains("identities: &mut F5FounderIdentities,\n    recount: R,"));
    assert!(
        !source.contains(".for_each_assignment("),
        "el único camino que deriva FHA_ no puede llamarse desde este módulo"
    );
    // Y no escribe: ninguna vía de escritura del sistema de ficheros.
    for writer in [
        "File::create",
        "fs::write",
        "fs::create_dir",
        "fs::remove",
        "fs::rename",
        "OpenOptions",
    ] {
        assert!(
            !source.contains(writer),
            "la clausura no puede escribir: {writer}"
        );
    }
    assert!(source.contains("File::open"), "pero sí abre para leer");
}

/// La posición que la auditoría indexa es la de `f5-founders.tsv`.
///
/// Sin esta prueba, reordenar el arreglo de contratos convertiría la clausura
/// en la auditoría silenciosa de otro fichero.
#[test]
fn the_founders_contract_is_the_one_the_audit_indexes() {
    let contract = &F5_TSV_CONTRACTS[FOUNDERS];
    assert_eq!(contract.artifact_kind(), "f5_founders");
    assert_eq!(contract.path(), "outputs/f5-founders.tsv");
    assert_eq!(
        contract.header(),
        "population_id\tindividual_id\tsequence_id\tcopy_index\tblock_id\tfounder_copy_id\tblock_assignment_id\thaplotype_id"
    );
}

// ---------------------------------------------------------------------------
// 7. Memoria
// ---------------------------------------------------------------------------

thread_local! {
    static LIVE: Cell<i64> = const { Cell::new(0) };
    static PEAK: Cell<i64> = const { Cell::new(0) };
    static TRACKING: Cell<bool> = const { Cell::new(false) };
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

/// La clausura **no** añade una tercera pieza fuera de
/// `f5-resource-accounting-v8` junto al oráculo y la publicación: su pico cabe
/// dentro del presupuesto de memoria de trabajo que la admisión selló.
///
/// Lo que el pico sigue es el **mayor registro publicado**, porque la
/// auditoría sostiene un registro cada vez: la línea que acaba de leer. No
/// sigue a `B*C_total`, ni al censo, ni al número de identidades, y la prueba
/// mueve **las dos** magnitudes por separado: bloques con censo fijo y censo
/// con bloques fijos. La cota que se comprueba es mecánica: el doble del
/// registro más largo —`read_line` duplica su búfer al crecer— más 8 KiB para
/// el `BufReader` y las cuatro estructuras deserializadas.
#[test]
fn the_closure_audit_peak_fits_the_sealed_budget_and_holds_one_record_at_a_time() {
    let cases: Vec<(&str, Fixture<'_>)> = vec![
        // Bloques variables, censo fijo en un individuo diploide.
        (
            "blocks-8",
            Fixture {
                block_count: 8,
                slack: 2,
                ..Fixture::default()
            },
        ),
        (
            "blocks-32",
            Fixture {
                block_count: 32,
                slack: 2,
                ..Fixture::default()
            },
        ),
        // Censo variable, bloques fijos: de dos copias a dieciséis.
        (
            "census-1",
            Fixture {
                block_count: 8,
                slack: 2,
                subpopulations: vec![Subpopulation::new("P1", 1)],
                ..Fixture::default()
            },
        ),
        (
            "census-8",
            Fixture {
                block_count: 8,
                slack: 2,
                subpopulations: vec![
                    Subpopulation::new("PA", 4),
                    Subpopulation::new("PB", 4).keys(),
                ],
                ..Fixture::default()
            },
        ),
    ];
    let mut measured = Vec::new();
    for (label, fixture) in cases {
        let package = publish(&format!("peak-{label}"), &fixture);
        let sealed_peak = package.chain.sealed.working_memory_bytes();
        let (report, peak) = peak_bytes(|| package.audit());
        let report = report.unwrap_or_else(|rejection| panic!("{label}: {rejection:?}"));
        assert!(
            peak <= sealed_peak,
            "{label}: la clausura pica {peak} bytes y la admisión selló {sealed_peak}"
        );
        // El mayor registro publicado, medido aparte sobre los bytes.
        let longest = truth_lines(&package.destination)
            .iter()
            .map(String::len)
            .max()
            .expect("alguna línea") as u64;
        let bound = 2 * longest + 8192;
        assert!(
            peak <= bound,
            "{label}: la clausura pica {peak} bytes y un registro cada vez cabe en {bound}"
        );
        measured.push((report.copies(), peak, longest));
    }
    // El censo sí se movió: si las dos últimas ejecuciones publicaran las
    // mismas copias, la mitad «no crece con el censo» no se habría ejercido.
    let (small_copies, small_peak, small_longest) = measured[2];
    let (large_copies, large_peak, large_longest) = measured[3];
    assert!(
        large_copies >= 8 * small_copies,
        "el censo tiene que multiplicarse: {large_copies} frente a {small_copies} copias"
    );
    // Y lo que explica el crecimiento del pico **no** es el censo sino el
    // registro: multiplicar por ocho las copias solo mueve el pico lo que se
    // movió el mayor registro publicado.
    let explained = 2 * large_longest.saturating_sub(small_longest) + 8192;
    assert!(
        large_peak.saturating_sub(small_peak) <= explained,
        "el pico creció {} bytes con el censo y el mayor registro solo explica {explained}",
        large_peak.saturating_sub(small_peak)
    );
}
