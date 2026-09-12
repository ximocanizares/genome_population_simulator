//! F5.4.1. Los siete TSV de éxito y la proyección fundadora `truth.jsonl` v2.
//!
//! Este módulo **no decide nada**: escribe lo que F5.1, F5.2 y F5.3 ya
//! cerraron, con los identificadores, clases, claves y orden que DTC-089 fija
//! en «Formatos» y en «IDs, RNG y formatos de salida». No extiende ni modifica
//! ningún esquema publicado.
//!
//! ## Tres reglas de escritura
//!
//! 1. **Streaming.** Ninguna tabla completa se materializa en memoria para
//!    ordenarla: cada fuente llega ya en la clave de su fichero y se recorre
//!    una vez. El único búfer vivo es el de la copia que se está serializando,
//!    acotado por los sitios de su secuencia. Ese búfer es de `CopySnp`, no
//!    de JSON: el registro `haplotype_copy` se serializa **por flujo** y sus
//!    tres vistas —modificaciones, fragmento vigente e historial— recorren el
//!    mismo búfer tres veces sin construir ningún `Value` intermedio. Un
//!    registro completo materializado costaría del orden de un kilobyte por
//!    SNP y por vista, y `final_formula_peak` no reserva ni un byte por copia
//!    serializada: `truth_bound` solo entra en `output_bytes`.
//! 2. **Una sola pasada de asignaciones.** `for_each_assignment` deriva cada
//!    `FHA_` exactamente una vez; derivarlo dos veces es una colisión. Esa
//!    pasada alimenta a la vez las filas de `f5-founders.tsv` y los registros
//!    `haplotype_copy` de `truth.jsonl`, que comparten su clave —población,
//!    individuo, secuencia y fase—.
//! 3. **Sin estructuras nuevas retenidas.** Los índices auxiliares son dos
//!    mapas de posición de tamaño `B` y `B*N_pop`, y `B` vectores de hasta
//!    `S_bound` índices para ordenar los sitios por `site_id`: muy por debajo
//!    de lo que `final_formula_peak` ya concede a esas dimensiones, que
//!    incluyen un término `256*S_bound` por bloque×subpoblación.
//!
//! ## Interpretaciones del contrato, declaradas
//!
//! - **Cada término de la clave de un fichero nombra una columna de ese
//!   fichero.** La regla es completa y cierra sola las siete claves de
//!   DTC-089: «bloque» es `block_id` porque esa es la columna; «índice» de
//!   `f5-sites.tsv` es `site_index`; «rango» de `f5-haplotypes.tsv`,
//!   `haplotype_rank`; «intento» de `f5-attempts.tsv`, `attempt_index`; y
//!   «sitio» de `f5-frequencies.tsv` es `site_id`, porque ese fichero **no
//!   publica** `site_index` y su única columna de sitio es `site_id`. El
//!   séptimo, `f5-block-layout.tsv`, escribe su clave entera —«FASTA, inicio,
//!   fin, ID»— porque «FASTA» no es ninguna de sus columnas. Los genotipos de
//!   `truth.jsonl` siguen la misma lectura: su clave «referencia, población,
//!   bloque, sitio e individuo» usa el `site_id` de `f5_locus`, que es el
//!   único campo de sitio del registro.
//! - **«bloque» como clave de orden es `block_id` ascendente.** Es la lectura
//!   que F5.3 ya cerró y que el plan prohíbe reabrir: DTC-089 ordena la lista
//!   `sampling` del efectivo por `population_id`+`block_id`, `assign_f5_founders`
//!   ordena por esa misma pareja y F5.3.3 declara que su recorrido de
//!   asignaciones **es** el orden de `f5-founders.tsv` y que sus unidades
//!   **son** el de `f5-haplotypes.tsv`. Los seis ficheros que dicen «bloque»
//!   usan por tanto `block_id`. El séptimo, `f5-block-layout.tsv`, no dice
//!   «bloque»: su clave está escrita entera —«FASTA, inicio, fin, ID»— y es la
//!   geometría, es decir `block_index` ascendente. Esa distinción es lo que
//!   permite a un lector recuperar `block_index` —componente del `BLK_` que
//!   ninguna columna publica— desde el ordinal de fila del layout.
//! - **La copia fundadora es la referencia con los SNP de sus bloques
//!   aplicados.** Por eso cada `haplotype_copy` publica un segmento por
//!   secuencia completa —`[0, sequence_length)`— cuyas modificaciones son los
//!   sitios donde su fila haplotípica lleva el alelo alternativo. La longitud
//!   de la secuencia no vive en ninguna pieza F5 previa y por eso es una
//!   entrada explícita: ver [`F5ReferenceSequences`].
//! - **Un bit a 1 de la matriz es el alelo alternativo**, con independencia de
//!   la polaridad: DTC-089 usa la polaridad solo en la regla de signos de LD y
//!   la publica en `f5-frequencies.tsv`.
//! - **`reference_content_sha256` es el SHA-256 del FASTA sellado.** F5 no
//!   transforma la referencia, de modo que la huella de su contenido y la del
//!   artefacto declarado coinciden por construcción.
//! - **Cada genotipo publica siempre los dos alelos del sitio** —referencia y
//!   alternativo—, también cuando la dosis de uno es cero. Así el recuento
//!   MAN-70 por dosis `0..ploidy` lee un índice de alelo estable en todos los
//!   portadores en lugar de deducirlo de los que resultaron presentes.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;

use serde_json::{Value, json};

use crate::application::EffectiveF5Config;
use crate::diagnostics::Diagnostic;
use crate::f5_founders::{
    F5BlockAssignment, F5FounderAssignment, F5FounderIdentities, f5_haplotype_alleles,
};
use crate::f5_haplotypes::F5AlleleMatrix;
use crate::f5_layout::{F5LayoutBlock, F5PhysicalLayout};
use crate::f5_preflight::{
    F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, F5_DEFAULT_MAX_SWAP_PROPOSALS, F5OtherResources,
    SealedF5Admission,
};
use crate::f5_search::F5HaplotypeSearch;
use crate::identity::{
    ExecutionId, F5BlockId, F5HaplotypeId, F5SiteId, FounderCopyId, IndividualId, PopulationId,
    SequenceId,
};
use crate::publication::ExecutionWorkspace;
use crate::truth_projection::case_id;
use crate::truth_schema::{
    TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TruthAllele, TruthSourceKind,
};

// ---------------------------------------------------------------------------
// Contratos literales
// ---------------------------------------------------------------------------

/// Versión de esquema de todo enlace TSV F5.
pub const F5_TSV_SCHEMA_VERSION: &str = "1";
/// Enlace obligatorio de la proyección fundadora.
pub const F5_TRUTH_ARTIFACT_ID: &str = "truth-jsonl";
pub const F5_TRUTH_ARTIFACT_KIND: &str = "truth_jsonl";
pub const F5_TRUTH_PATH: &str = "outputs/truth.jsonl";
pub const F5_TRUTH_SCHEMA_VERSION: &str = "2";
pub const F5_TRUTH_PROJECTION_VERSION: &str = "f5-founder-truth-projection-v2";
/// `projection.kind` de la cabecera v2.
pub const F5_TRUTH_PROJECTION_KIND: &str = "f5_founder";

/// Un TSV del contrato: identificador, clase, ruta y cabecera literales.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5TsvContract {
    artifact_id: &'static str,
    artifact_kind: &'static str,
    path: &'static str,
    header: &'static str,
}

impl F5TsvContract {
    #[must_use]
    pub const fn artifact_id(&self) -> &'static str {
        self.artifact_id
    }
    #[must_use]
    pub const fn artifact_kind(&self) -> &'static str {
        self.artifact_kind
    }
    #[must_use]
    pub const fn path(&self) -> &'static str {
        self.path
    }
    #[must_use]
    pub const fn header(&self) -> &'static str {
        self.header
    }
}

/// Los siete TSV exactos de un éxito, en el orden de la tabla de DTC-089, que
/// es también el orden lexicográfico de `artifact_id`.
pub const F5_TSV_CONTRACTS: [F5TsvContract; 7] = [
    F5TsvContract {
        artifact_id: "f5-01-block-layout-v1",
        artifact_kind: "f5_block_layout",
        path: "outputs/f5-block-layout.tsv",
        header: "block_id\tsequence_id\tstart_1based\tend_1based\tlength_bp",
    },
    F5TsvContract {
        artifact_id: "f5-02-sites-v1",
        artifact_kind: "f5_sites",
        path: "outputs/f5-sites.tsv",
        header: "block_id\tsite_index\tsite_id\tsequence_id\tposition_1based\tref\talt",
    },
    F5TsvContract {
        artifact_id: "f5-03-haplotypes-v1",
        artifact_kind: "f5_haplotypes",
        path: "outputs/f5-haplotypes.tsv",
        header: "population_id\tblock_id\thaplotype_id\thaplotype_rank\tq\talleles",
    },
    F5TsvContract {
        artifact_id: "f5-04-frequencies-v1",
        artifact_kind: "f5_frequencies",
        path: "outputs/f5-frequencies.tsv",
        header: "population_id\tblock_id\tsite_id\tpolarity\ttarget_alt_af\trealized_alt_af\tabs_error",
    },
    F5TsvContract {
        artifact_id: "f5-05-ld-evaluation-v1",
        artifact_kind: "f5_ld_evaluation",
        path: "outputs/f5-ld-evaluation.tsv",
        header: "population_id\tblock_id\tdistance_start_bp\tdistance_end_bp\tstatus\tpair_count\tmean_r2\ttarget_mean_r2\tmax_abs_error\tabs_error",
    },
    F5TsvContract {
        artifact_id: "f5-06-founders-v1",
        artifact_kind: "f5_founders",
        path: "outputs/f5-founders.tsv",
        header: "population_id\tindividual_id\tsequence_id\tcopy_index\tblock_id\tfounder_copy_id\tblock_assignment_id\thaplotype_id",
    },
    F5TsvContract {
        artifact_id: "f5-07-attempts-v1",
        artifact_kind: "f5_attempts",
        path: "outputs/f5-attempts.tsv",
        header: "population_id\tblock_id\tattempt_index\tproposals\taccepted\tduplicate_rows\tbins_outside\tsigns_failed\tmax_r2_error\tsum_r2_error\tmax_marginal_error\tselected_best",
    },
];

const BLOCK_LAYOUT: usize = 0;
const SITES: usize = 1;
const HAPLOTYPES: usize = 2;
const FREQUENCIES: usize = 3;
const LD_EVALUATION: usize = 4;
const FOUNDERS: usize = 5;
const ATTEMPTS: usize = 6;

/// Literal que DTC-089 reserva a `mean_r2` y `abs_error` de un bin no
/// aplicable. Ningún otro campo queda vacío ni usa este literal.
const NOT_AVAILABLE: &str = "NA";

// ---------------------------------------------------------------------------
// Enlaces publicados
// ---------------------------------------------------------------------------

/// Enlace cerrado de un TSV F5: `artifact_id`, `artifact_kind`,
/// `schema_version="1"`, `path`, `size_bytes:u64` y `sha256`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5ArtifactLink {
    artifact_id: &'static str,
    artifact_kind: &'static str,
    path: &'static str,
    size_bytes: u64,
    sha256: String,
}

impl F5ArtifactLink {
    #[must_use]
    pub const fn artifact_id(&self) -> &'static str {
        self.artifact_id
    }
    #[must_use]
    pub const fn artifact_kind(&self) -> &'static str {
        self.artifact_kind
    }
    #[must_use]
    pub const fn schema_version(&self) -> &'static str {
        F5_TSV_SCHEMA_VERSION
    }
    #[must_use]
    pub const fn path(&self) -> &'static str {
        self.path
    }
    #[must_use]
    pub const fn size_bytes(&self) -> u64 {
        self.size_bytes
    }
    #[must_use]
    pub fn sha256(&self) -> &str {
        &self.sha256
    }

    /// El objeto exacto que el manifiesto y la configuración efectiva copian.
    #[must_use]
    pub fn to_json(&self) -> Value {
        json!({
            "artifact_id": self.artifact_id,
            "artifact_kind": self.artifact_kind,
            "schema_version": F5_TSV_SCHEMA_VERSION,
            "path": self.path,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        })
    }
}

/// Enlace `truth` obligatorio: separado de los siete, con su propia versión de
/// esquema y su versión de proyección.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5TruthLink {
    size_bytes: u64,
    sha256: String,
}

impl F5TruthLink {
    #[must_use]
    pub const fn artifact_id(&self) -> &'static str {
        F5_TRUTH_ARTIFACT_ID
    }
    #[must_use]
    pub const fn artifact_kind(&self) -> &'static str {
        F5_TRUTH_ARTIFACT_KIND
    }
    #[must_use]
    pub const fn schema_version(&self) -> &'static str {
        F5_TRUTH_SCHEMA_VERSION
    }
    #[must_use]
    pub const fn projection_version(&self) -> &'static str {
        F5_TRUTH_PROJECTION_VERSION
    }
    #[must_use]
    pub const fn path(&self) -> &'static str {
        F5_TRUTH_PATH
    }
    #[must_use]
    pub const fn size_bytes(&self) -> u64 {
        self.size_bytes
    }
    #[must_use]
    pub fn sha256(&self) -> &str {
        &self.sha256
    }

    #[must_use]
    pub fn to_json(&self) -> Value {
        json!({
            "artifact_id": F5_TRUTH_ARTIFACT_ID,
            "artifact_kind": F5_TRUTH_ARTIFACT_KIND,
            "schema_version": F5_TRUTH_SCHEMA_VERSION,
            "projection_version": F5_TRUTH_PROJECTION_VERSION,
            "path": F5_TRUTH_PATH,
            "size_bytes": self.size_bytes,
            "sha256": self.sha256,
        })
    }
}

/// Todo lo que una ejecución con éxito publica en F5.4.1.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5SuccessArtifacts {
    artifacts: [F5ArtifactLink; 7],
    truth: F5TruthLink,
    copy_count: u64,
    genotype_count: u64,
}

impl F5SuccessArtifacts {
    /// Los siete enlaces, una vez y en el orden de la tabla.
    #[must_use]
    pub const fn artifacts(&self) -> &[F5ArtifactLink; 7] {
        &self.artifacts
    }
    #[must_use]
    pub const fn truth(&self) -> &F5TruthLink {
        &self.truth
    }
    /// `counts.copies` de la cabecera v2.
    #[must_use]
    pub const fn copy_count(&self) -> u64 {
        self.copy_count
    }
    /// `counts.genotypes` de la cabecera v2.
    #[must_use]
    pub const fn genotype_count(&self) -> u64 {
        self.genotype_count
    }
}

// ---------------------------------------------------------------------------
// Longitudes de secuencia
// ---------------------------------------------------------------------------

/// Longitudes de las secuencias de la referencia sellada.
///
/// `truth.jsonl` exige `sequence_length` en cada fuente de segmento y ninguna
/// pieza F5 anterior la conserva: el layout retiene bloques, no secuencias, y
/// la admisión retiene tramos elegibles, cuyo último extremo no es la longitud.
/// Por eso es una entrada explícita del escritor y no un valor inventado.
/// Una entrada por secuencia retenida es exactamente lo que
/// `256*retained_sequence_ids` de `final_formula_peak` ya concede.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct F5ReferenceSequences {
    lengths: BTreeMap<SequenceId, u64>,
}

impl F5ReferenceSequences {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Declara la longitud de una secuencia. Devuelve `false` —y no cambia
    /// nada— si la longitud es cero o si la secuencia ya estaba declarada.
    ///
    /// El resultado no puede ignorarse: una declaración rechazada en silencio
    /// publicaría después una copia sin longitud o con la de otra secuencia.
    #[must_use = "una declaración rechazada deja la secuencia sin longitud publicable"]
    pub fn insert(&mut self, sequence_id: SequenceId, length_bp: u64) -> bool {
        if length_bp == 0 || self.lengths.contains_key(&sequence_id) {
            return false;
        }
        self.lengths.insert(sequence_id, length_bp);
        true
    }

    #[must_use]
    pub fn length(&self, sequence_id: &SequenceId) -> Option<u64> {
        self.lengths.get(sequence_id).copied()
    }

    /// Longitudes declaradas, en el orden canónico de `sequence_id`.
    pub fn lengths(&self) -> impl ExactSizeIterator<Item = u64> + '_ {
        self.lengths.values().copied()
    }

    #[must_use]
    pub fn len(&self) -> usize {
        self.lengths.len()
    }

    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.lengths.is_empty()
    }
}

// ---------------------------------------------------------------------------
// Frontera pública
// ---------------------------------------------------------------------------

/// Las piezas cerradas que un éxito F5 publica, enlazadas entre sí.
///
/// El constructor rechaza antes de abrir ningún fichero una admisión que no
/// sea la de esta configuración, de modo que ningún artefacto parcial puede
/// nacer de una combinación de ejecuciones.
#[derive(Clone, Copy, Debug)]
pub struct F5SuccessRequest<'a> {
    config: &'a EffectiveF5Config,
    execution_id: &'a ExecutionId,
    layout: &'a F5PhysicalLayout,
    search: &'a F5HaplotypeSearch,
    assignment: &'a F5FounderAssignment,
    sequences: &'a F5ReferenceSequences,
    reference_sha256: &'a str,
    profile_id: &'a str,
    profile_metadata_sha256: &'a str,
    /// Pico de IDs FASTA que la admisión selló: ninguna declaración de
    /// longitudes puede nombrar más secuencias que la referencia admitida.
    retained_sequence_ids: u64,
    /// Bases de la referencia que el escáner de admisión recorrió: la suma de
    /// las longitudes declaradas no puede excederlas.
    reference_bases_processed: u64,
}

impl<'a> F5SuccessRequest<'a> {
    /// # Errors
    ///
    /// `GPS-CONF-101` si la admisión sellada no corresponde a esta
    /// configuración.
    ///
    /// El enlace de `sequences` con la admisión no se cierra aquí, sino en
    /// [`write_f5_success_artifacts`], que exige además que la declaración no
    /// nombre más secuencias que las selladas, que ninguna secuencia poblada
    /// falte, que todo bloque quepa en la suya y que las longitudes declaradas
    /// no sumen más bases que las que el escáner de admisión recorrió.
    pub fn new(
        config: &'a EffectiveF5Config,
        sealed: &'a SealedF5Admission,
        execution_id: &'a ExecutionId,
        layout: &'a F5PhysicalLayout,
        search: &'a F5HaplotypeSearch,
        assignment: &'a F5FounderAssignment,
        sequences: &'a F5ReferenceSequences,
    ) -> Result<Self, Diagnostic> {
        if !sealed.matches_admission(config, F5OtherResources::default()) {
            return Err(configuration(
                config,
                "preflight.binding",
                "la admisión sellada no corresponde a esta configuración",
            ));
        }
        Ok(Self {
            config,
            execution_id,
            layout,
            search,
            assignment,
            sequences,
            reference_sha256: sealed.reference_sha256(),
            profile_id: sealed.profile_id(),
            profile_metadata_sha256: sealed.profile_metadata_sha256(),
            retained_sequence_ids: sealed.retained_sequence_ids(),
            reference_bases_processed: sealed.reference_bases_processed(),
        })
    }
}

/// Escribe los siete TSV de éxito y `outputs/truth.jsonl` v2.
///
/// Consume la pasada única de asignaciones del repositorio de identidades: al
/// volver, `identities` ya no puede recorrerse otra vez, porque derivar dos
/// veces un `FHA_` es una colisión.
///
/// # Errors
///
/// - `GPS-CONF-101` si las piezas recibidas no son las de una misma ejecución
///   —incluido un `execution_id` distinto del que derivó las identidades—, si
///   una secuencia poblada no declara su longitud, si un bloque no cabe en
///   ella, si la declaración nombra más secuencias que las que la admisión
///   selló o si sus longitudes suman más bases que las que el escáner de
///   admisión recorrió.
/// - `GPS-INTERNAL-101` si una pieza publicada incumple una invariante que sus
///   productores ya demostraron —orden, cardinalidad, dominio numérico o
///   correspondencia posicional—.
/// - `GPS-DOMAIN-002` si la derivación de un `FHA_` repite un identificador.
/// - `GPS-IO-001` ante un fallo de escritura del provisional.
pub fn write_f5_success_artifacts(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    identities: &mut F5FounderIdentities,
) -> Result<F5SuccessArtifacts, Diagnostic> {
    let plan = WritePlan::prepare(request, identities)?;

    let block_layout = write_block_layout(workspace, request)?;
    let sites = write_sites(workspace, request, &plan)?;
    let haplotypes = write_haplotypes(workspace, request, identities, &plan)?;
    let frequencies = write_frequencies(workspace, request, identities, &plan)?;
    let ld_evaluation = write_ld_evaluation(workspace, request, &plan)?;
    let attempts = write_attempts(workspace, request, identities, &plan)?;
    // La pasada irreversible va la última: cualquier incumplimiento anterior
    // ya ha abortado sin haber derivado ningún `FHA_`.
    let (founders, truth) = write_founders_and_truth(workspace, request, identities, &plan)?;

    Ok(F5SuccessArtifacts {
        artifacts: [
            block_layout,
            sites,
            haplotypes,
            frequencies,
            ld_evaluation,
            founders,
            attempts,
        ],
        truth,
        copy_count: plan.copy_count,
        genotype_count: plan.genotype_count,
    })
}

// ---------------------------------------------------------------------------
// Plan de escritura: correspondencias verificadas y recuentos de cabecera
// ---------------------------------------------------------------------------

/// Índices y recuentos que las ocho salidas comparten.
///
/// No es una tabla de filas: tiene una entrada por bloque y otra por
/// bloque×subpoblación, y ninguna crece con las copias ni con los sitios.
struct WritePlan {
    /// Posición en `layout.blocks()` de cada `block_index`.
    block_positions: BTreeMap<u64, usize>,
    /// Posiciones de `layout.blocks()` ordenadas por `block_id`: la clave de
    /// `f5-sites.tsv`. Ordena `B` índices, no `B*S_bound` filas.
    blocks_by_id: Vec<usize>,
    /// Para cada posición de `layout.blocks()`, las posiciones de sus sitios
    /// ordenadas por `site_id`: la clave de `f5-frequencies.tsv` —cuya única
    /// columna de sitio es `site_id`— y la de los genotipos de `truth.jsonl`.
    /// Ordena hasta `S_bound` índices por bloque, no filas.
    sites_by_id: Vec<Vec<usize>>,
    /// Posición en `search.blocks()` —y en `search.ld_evaluation()`— de cada
    /// bloque×subpoblación, que llega en la clave de la tabla de marginales.
    search_positions: BTreeMap<(PopulationId, u64), usize>,
    /// Para cada posición de `identities.units()`, la suya en `search.blocks()`.
    unit_to_search: Vec<usize>,
    /// Tramo de `search.records()` de cada posición de `search.blocks()`.
    attempt_ranges: Vec<(usize, usize)>,
    /// `individual_count` por subpoblación, en el orden del censo.
    individual_counts: BTreeMap<PopulationId, u64>,
    copy_count: u64,
    genotype_count: u64,
}

impl WritePlan {
    fn prepare(
        request: &F5SuccessRequest<'_>,
        identities: &F5FounderIdentities,
    ) -> Result<Self, Diagnostic> {
        let config = request.config;
        // Los identificadores los derivó el repositorio con **su** ejecución.
        // Publicar otra en la cabecera v2 dejaría un paquete cuyos `FCP_`,
        // `HAP_` y `FHA_` ningún reconstructor podría rederivar.
        if identities.execution_id() != request.execution_id {
            return Err(configuration(
                config,
                "execution.identity",
                format!(
                    "la cabecera declararía la ejecución {} y las identidades se derivaron con {}",
                    request.execution_id.as_str(),
                    identities.execution_id().as_str()
                ),
            ));
        }
        // Atadura de las longitudes declaradas con la admisión sellada: ni más
        // secuencias que las retenidas ni más bases que las recorridas.
        let declared = u64::try_from(request.sequences.len())
            .map_err(|_| internal(config, "el recuento de secuencias no cabe en u64"))?;
        if declared > request.retained_sequence_ids {
            return Err(configuration(
                config,
                "reference.sequence_lengths",
                format!(
                    "se declaran {declared} secuencias y la admisión selló {}",
                    request.retained_sequence_ids
                ),
            ));
        }
        let mut declared_bases = 0_u64;
        for length in request.sequences.lengths() {
            declared_bases = declared_bases
                .checked_add(length)
                .ok_or_else(|| internal(config, "las longitudes declaradas desbordan u64"))?;
        }
        if declared_bases > request.reference_bases_processed {
            return Err(configuration(
                config,
                "reference.sequence_lengths",
                format!(
                    "las longitudes declaradas suman {declared_bases} bases y la admisión recorrió {}",
                    request.reference_bases_processed
                ),
            ));
        }
        let mut block_positions = BTreeMap::new();
        let mut sites_by_id = Vec::with_capacity(request.layout.blocks().len());
        let mut previous: Option<(u64, u64, u64)> = None;
        for (position, block) in request.layout.blocks().iter().enumerate() {
            if block.length_bp() != block.end_1based().saturating_sub(block.start_1based()) + 1
                || block.start_1based() == 0
                || block.start_1based() > block.end_1based()
            {
                return Err(internal(
                    config,
                    "un bloque colocado no satisface length_bp=end-start+1 con coordenadas positivas",
                ));
            }
            let key = (
                block.sequence_ordinal(),
                block.start_1based(),
                block.end_1based(),
            );
            if previous.is_some_and(|last| last >= key) {
                return Err(internal(
                    config,
                    "el layout no llega en el orden canónico FASTA, inicio y fin",
                ));
            }
            previous = Some(key);
            let length = request
                .sequences
                .length(block.sequence_id())
                .ok_or_else(|| {
                    configuration(
                        config,
                        "reference.sequence_lengths",
                        format!(
                            "la secuencia {} recibe bloques y no declara su longitud",
                            block.sequence_id().as_str()
                        ),
                    )
                })?;
            if block.end_1based() > length {
                return Err(configuration(
                    config,
                    "reference.sequence_lengths",
                    format!(
                        "el bloque {} termina en {} y su secuencia {} declara {length} bases",
                        block.block_index(),
                        block.end_1based(),
                        block.sequence_id().as_str()
                    ),
                ));
            }
            if block_positions
                .insert(block.block_index(), position)
                .is_some()
            {
                return Err(internal(config, "el layout repite un índice de bloque"));
            }
            let mut ordered: Vec<usize> = (0..block.sites().len()).collect();
            ordered.sort_by(|left, right| {
                block.sites()[*left]
                    .site_id()
                    .cmp(block.sites()[*right].site_id())
            });
            if ordered
                .windows(2)
                .any(|pair| block.sites()[pair[0]].site_id() == block.sites()[pair[1]].site_id())
            {
                return Err(internal(
                    config,
                    "un bloque repite un identificador de sitio",
                ));
            }
            sites_by_id.push(ordered);
        }

        let mut blocks_by_id: Vec<usize> = (0..request.layout.blocks().len()).collect();
        blocks_by_id.sort_by(|left, right| {
            request.layout.blocks()[*left]
                .block_id()
                .cmp(request.layout.blocks()[*right].block_id())
        });

        let units = identities.units();
        if units.len() != request.assignment.units().len()
            || units.len() != request.search.blocks().len()
            || units.len() != request.search.ld_evaluation().len()
        {
            return Err(configuration(
                config,
                "units.identity",
                "las unidades materializadas, asignadas, buscadas y evaluadas no son las mismas",
            ));
        }
        // `search.blocks()` llega en la clave de la tabla de marginales
        // —población y orden canónico de bloque— y las unidades materializadas
        // en la de `f5-founders.tsv` —población y `block_id`—. Ninguna de las
        // dos se reordena: se indexan.
        let mut search_positions = BTreeMap::new();
        for (position, winner) in request.search.blocks().iter().enumerate() {
            let evaluation = &request.search.ld_evaluation()[position];
            if evaluation.population_id() != winner.population_id()
                || evaluation.block_index() != winner.block_index()
            {
                return Err(configuration(
                    config,
                    "units.identity",
                    "la evaluación LD no acompaña posicionalmente a su bloque ganador",
                ));
            }
            if search_positions
                .insert(
                    (winner.population_id().clone(), winner.block_index()),
                    position,
                )
                .is_some()
            {
                return Err(internal(
                    config,
                    "la búsqueda repite un bloque×subpoblación",
                ));
            }
        }
        let mut unit_to_search = Vec::with_capacity(units.len());
        let mut previous_unit: Option<(&PopulationId, &F5BlockId)> = None;
        for (position, unit) in units.iter().enumerate() {
            let assigned = &request.assignment.units()[position];
            if assigned.population_id() != unit.population_id()
                || assigned.block_id() != unit.block_id()
                || assigned.block_index() != unit.block_index()
            {
                return Err(configuration(
                    config,
                    "units.identity",
                    "una unidad bloque×subpoblación no ocupa la misma posición en sus dos fuentes",
                ));
            }
            let key = (unit.population_id(), unit.block_id());
            if previous_unit.is_some_and(|last| last >= key) {
                return Err(internal(
                    config,
                    "las unidades no llegan por población y block_id, que es la clave publicada",
                ));
            }
            previous_unit = Some(key);
            let found = *search_positions
                .get(&(unit.population_id().clone(), unit.block_index()))
                .ok_or_else(|| {
                    internal(config, "una unidad no tiene bloque ganador en la búsqueda")
                })?;
            unit_to_search.push(found);
        }
        // Tramos del historial: llega agrupado por unidad en la clave de la
        // búsqueda, de modo que basta reconocer los cambios de grupo.
        let mut attempt_ranges = vec![(0_usize, 0_usize); units.len()];
        let mut current: Option<(usize, usize)> = None;
        for (index, record) in request.search.records().iter().enumerate() {
            let position = *search_positions
                .get(&(record.population_id().clone(), record.block_index()))
                .ok_or_else(|| {
                    internal(
                        config,
                        "un intento nombra un bloque×subpoblación desconocido",
                    )
                })?;
            match current {
                Some((open, start)) if open == position => {
                    attempt_ranges[open] = (start, index + 1);
                }
                _ => {
                    if attempt_ranges[position] != (0, 0) {
                        return Err(internal(
                            config,
                            "el historial de intentos no llega agrupado por bloque×subpoblación",
                        ));
                    }
                    attempt_ranges[position] = (index, index + 1);
                    current = Some((position, index));
                }
            }
        }

        let mut individual_counts = BTreeMap::new();
        for roster in request.assignment.rosters() {
            individual_counts.insert(roster.population_id().clone(), roster.individual_count());
        }

        let copy_count = u64::try_from(identities.copies().len())
            .map_err(|_| internal(config, "el recuento de copias no cabe en u64"))?;
        let mut genotype_count = 0_u64;
        for unit in units {
            let position = *block_positions.get(&unit.block_index()).ok_or_else(|| {
                internal(
                    config,
                    "una unidad nombra un bloque que el layout no colocó",
                )
            })?;
            let sites = u64::try_from(request.layout.blocks()[position].sites().len())
                .map_err(|_| internal(config, "el recuento de sitios no cabe en u64"))?;
            let individuals = *individual_counts
                .get(unit.population_id())
                .ok_or_else(|| internal(config, "una unidad nombra una subpoblación sin censo"))?;
            genotype_count = individuals
                .checked_mul(sites)
                .and_then(|product| genotype_count.checked_add(product))
                .ok_or_else(|| internal(config, "el recuento de genotipos desborda u64"))?;
        }

        Ok(Self {
            block_positions,
            blocks_by_id,
            sites_by_id,
            search_positions,
            unit_to_search,
            attempt_ranges,
            individual_counts,
            copy_count,
            genotype_count,
        })
    }

    /// El bloque colocado de un `block_index`, sin recorrer el layout.
    fn block<'b>(
        &self,
        layout: &'b F5PhysicalLayout,
        block_index: u64,
    ) -> Option<&'b F5LayoutBlock> {
        let position = *self.block_positions.get(&block_index)?;
        layout.blocks().get(position)
    }

    /// El bloque colocado y las posiciones de sus sitios en orden de
    /// `site_id`, que es la clave de `f5-frequencies.tsv` y la de los
    /// genotipos.
    fn block_sites_by_id<'b>(
        &self,
        layout: &'b F5PhysicalLayout,
        block_index: u64,
    ) -> Option<(&'b F5LayoutBlock, &[usize])> {
        let position = *self.block_positions.get(&block_index)?;
        Some((
            layout.blocks().get(position)?,
            self.sites_by_id.get(position)?,
        ))
    }
}

// ---------------------------------------------------------------------------
// Los cinco TSV que se escriben de una sola fuente
// ---------------------------------------------------------------------------

/// Abre el artefacto, escribe su cabecera y delega las filas, todo por flujo.
fn stream_tsv(
    workspace: &mut ExecutionWorkspace,
    contract: &F5TsvContract,
    rows: impl FnOnce(&mut BufWriter<&mut File>) -> Result<(), Diagnostic>,
) -> Result<F5ArtifactLink, Diagnostic> {
    let record = workspace.write_streamed_artifact(
        contract.artifact_id,
        contract.artifact_kind,
        Path::new(contract.path),
        |file| {
            let mut writer = BufWriter::new(file);
            writeln!(writer, "{}", contract.header)
                .map_err(|error| io(contract.path, "write_header", &error))?;
            rows(&mut writer)?;
            writer
                .flush()
                .map_err(|error| io(contract.path, "flush", &error))
        },
    )?;
    Ok(F5ArtifactLink {
        artifact_id: contract.artifact_id,
        artifact_kind: contract.artifact_kind,
        path: contract.path,
        size_bytes: record.size_bytes,
        sha256: record.digest,
    })
}

fn write_block_layout(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[BLOCK_LAYOUT];
    let layout = request.layout;
    let config = request.config;
    stream_tsv(workspace, contract, |writer| {
        for (ordinal, block) in layout.blocks().iter().enumerate() {
            // El ordinal de fila es la **única** vía por la que un lector
            // recupera `block_index`, componente del material del `BLK_` que
            // ninguna columna publica. Publicar el layout en otro orden lo
            // haría irrecuperable, de modo que la igualdad se afirma aquí.
            if block.block_index() != ordinal as u64 {
                return Err(internal(
                    config,
                    "el ordinal de fila del layout no es el block_index que el BLK_ deriva",
                ));
            }
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}",
                block.block_id().as_str(),
                block.sequence_id().as_str(),
                block.start_1based(),
                block.end_1based(),
                block.length_bp(),
            )
            .map_err(|error| io(contract.path, "write_row", &error))?;
        }
        Ok(())
    })
}

fn write_sites(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    plan: &WritePlan,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[SITES];
    let config = request.config;
    let layout = request.layout;
    stream_tsv(workspace, contract, |writer| {
        for position in &plan.blocks_by_id {
            let block = &layout.blocks()[*position];
            for (expected, site) in block.sites().iter().enumerate() {
                if site.site_index() != expected as u64 {
                    return Err(internal(
                        config,
                        "los índices de sitio de un bloque no son consecutivos desde 0",
                    ));
                }
                if site.position_1based() < block.start_1based()
                    || site.position_1based() > block.end_1based()
                {
                    return Err(internal(config, "un sitio cae fuera de su bloque"));
                }
                let reference = canonical_base(config, site.reference_base())?;
                let alternate = canonical_base(config, site.alternate_base())?;
                if reference == alternate {
                    return Err(internal(
                        config,
                        "un sitio bialélico publica el mismo alelo como referencia y alternativo",
                    ));
                }
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{reference}\t{alternate}",
                    block.block_id().as_str(),
                    site.site_index(),
                    site.site_id().as_str(),
                    block.sequence_id().as_str(),
                    site.position_1based(),
                )
                .map_err(|error| io(contract.path, "write_row", &error))?;
            }
        }
        Ok(())
    })
}

fn write_haplotypes(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    identities: &F5FounderIdentities,
    plan: &WritePlan,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[HAPLOTYPES];
    let config = request.config;
    let weights = config.profile().haplotype_weights();
    let mut total = 0.0_f64;
    for weight in weights {
        if !weight.is_finite() || *weight <= 0.0 || *weight > 1.0 {
            return Err(internal(
                config,
                "un peso haplotípico efectivo no es finito ni pertenece a (0,1]",
            ));
        }
        total += *weight;
    }
    // Los pesos son los mismos en todas las unidades, de modo que comprobar
    // una vez que suman 1 lo comprueba por población×bloque, que es como
    // DTC-089 lo exige.
    if (total - 1.0).abs() > 1e-9 {
        return Err(internal(
            config,
            format!(
                "los pesos haplotípicos publicados suman {total} y no 1 por población y bloque"
            ),
        ));
    }
    let units = identities.units();
    let search = request.search;
    stream_tsv(workspace, contract, |writer| {
        for (position, unit) in units.iter().enumerate() {
            let matrix = search.blocks()[plan.unit_to_search[position]].matrix();
            if matrix.haplotype_count() != weights.len()
                || unit.haplotype_ids().len() != weights.len()
            {
                return Err(internal(
                    config,
                    "la matriz ganadora y el perfil no publican el mismo número de filas",
                ));
            }
            for (rank, haplotype_id) in unit.haplotype_ids().iter().enumerate() {
                let alleles = f5_haplotype_alleles(matrix, rank).ok_or_else(|| {
                    internal(config, "la matriz ganadora no conserva una de sus filas")
                })?;
                if alleles.len() != matrix.site_count() {
                    return Err(internal(
                        config,
                        "la cadena de alelos no mide exactamente S caracteres",
                    ));
                }
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{rank}\t{}\t{alleles}",
                    unit.population_id().as_str(),
                    unit.block_id().as_str(),
                    haplotype_id.as_str(),
                    decimal(weights[rank]),
                )
                .map_err(|error| io(contract.path, "write_row", &error))?;
            }
        }
        Ok(())
    })
}

fn write_frequencies(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    identities: &F5FounderIdentities,
    plan: &WritePlan,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[FREQUENCIES];
    let config = request.config;
    let search = request.search;
    let layout = request.layout;
    let units = identities.units();
    stream_tsv(workspace, contract, |writer| {
        for (position, unit) in units.iter().enumerate() {
            let block = &search.blocks()[plan.unit_to_search[position]];
            let (placed, by_id) = plan
                .block_sites_by_id(layout, unit.block_index())
                .ok_or_else(|| internal(config, "una unidad nombra un bloque sin geometría"))?;
            if placed.sites().len() != block.columns().len() {
                return Err(internal(
                    config,
                    "las columnas realizadas y los sitios del bloque no coinciden",
                ));
            }
            // Las columnas realizadas llegan en el orden posicional del
            // bloque: primero se comprueba esa correspondencia y después se
            // publican en la clave del fichero, que es `site_id`.
            for (expected, column) in block.columns().iter().enumerate() {
                if column.site_index() != expected as u64
                    || column.site_id() != placed.sites()[expected].site_id()
                {
                    return Err(internal(
                        config,
                        "una columna realizada no nombra el sitio de su propia posición",
                    ));
                }
            }
            for site_position in by_id {
                let column = &block.columns()[*site_position];
                for value in [
                    column.target_alt_af(),
                    column.realized_alt_af(),
                    column.abs_error(),
                ] {
                    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
                        return Err(internal(
                            config,
                            "una frecuencia o error publicado no pertenece a [0,1]",
                        ));
                    }
                }
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    block.population_id().as_str(),
                    block.block_id().as_str(),
                    column.site_id().as_str(),
                    column.polarity().as_str(),
                    decimal(column.target_alt_af()),
                    decimal(column.realized_alt_af()),
                    decimal(column.abs_error()),
                )
                .map_err(|error| io(contract.path, "write_row", &error))?;
            }
        }
        Ok(())
    })
}

fn write_ld_evaluation(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    plan: &WritePlan,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[LD_EVALUATION];
    let config = request.config;
    let search = request.search;
    // `unit_to_search` ya está en la clave publicada —población y `block_id`—
    // y apunta a la evaluación que acompaña a cada bloque ganador.
    let order = &plan.unit_to_search;
    stream_tsv(workspace, contract, |writer| {
        for position in order {
            let evaluation = &search.ld_evaluation()[*position];
            let mut previous: Option<u64> = None;
            for bin in evaluation.bins() {
                if previous.is_some_and(|last| last >= bin.distance_start_bp()) {
                    return Err(internal(
                        config,
                        "los bins LD no llegan ordenados por inicio de bin",
                    ));
                }
                previous = Some(bin.distance_start_bp());
                if bin.is_evaluated() == (bin.pair_count() == 0) {
                    return Err(internal(
                        config,
                        "pair_count vale cero si y solo si el bin es not_applicable",
                    ));
                }
                if bin.is_evaluated() == (bin.mean_r2().is_none() || bin.abs_error().is_none()) {
                    return Err(internal(
                        config,
                        "mean_r2 y abs_error faltan si y solo si el bin es not_applicable",
                    ));
                }
                let mean = match bin.mean_r2() {
                    None => NOT_AVAILABLE.to_owned(),
                    Some(value) => bounded_decimal(config, value)?,
                };
                let error = match bin.abs_error() {
                    None => NOT_AVAILABLE.to_owned(),
                    Some(value) => bounded_decimal(config, value)?,
                };
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{mean}\t{}\t{}\t{error}",
                    evaluation.population_id().as_str(),
                    evaluation.block_id().as_str(),
                    bin.distance_start_bp(),
                    bin.distance_end_bp(),
                    bin.status(),
                    bin.pair_count(),
                    bounded_decimal(config, bin.target_mean_r2())?,
                    bounded_decimal(config, bin.max_abs_error())?,
                )
                .map_err(|error| io(contract.path, "write_row", &error))?;
            }
        }
        Ok(())
    })
}

fn write_attempts(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    identities: &F5FounderIdentities,
    plan: &WritePlan,
) -> Result<F5ArtifactLink, Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[ATTEMPTS];
    let config = request.config;
    let search = request.search;
    let units = identities.units();
    // Los límites efectivos que la admisión resolvió para la búsqueda: cada
    // contador publicado debe caber en ellos, como DTC-089 exige.
    let scenario = config.declared_scenario();
    let max_attempts = scenario
        .max_haplotype_attempts
        .as_ref()
        .map_or(F5_DEFAULT_MAX_HAPLOTYPE_ATTEMPTS, |value| *value.value());
    let max_proposals = scenario
        .max_swap_proposals
        .as_ref()
        .map_or(F5_DEFAULT_MAX_SWAP_PROPOSALS, |value| *value.value());
    stream_tsv(workspace, contract, |writer| {
        // El historial no se copia ni se reordena: cada unidad publica su
        // tramo contiguo, ya indexado, en la clave de su fichero.
        for (position, unit) in units.iter().enumerate() {
            let (start, end) = plan.attempt_ranges[plan.unit_to_search[position]];
            if start == end {
                return Err(internal(
                    config,
                    "un bloque×subpoblación no publica ningún intento ejecutado",
                ));
            }
            let mut best = 0_u64;
            let mut expected_attempt = None;
            for record in &search.records()[start..end] {
                if expected_attempt.is_some_and(|last| last >= record.attempt_index()) {
                    return Err(internal(
                        config,
                        "los intentos de una unidad no llegan por índice ascendente",
                    ));
                }
                expected_attempt = Some(record.attempt_index());
                if record.attempt_index() >= max_attempts || record.proposals() > max_proposals {
                    return Err(internal(
                        config,
                        "un contador de intento publicado no cabe en los límites efectivos",
                    ));
                }
                if record.selected_best() {
                    best += 1;
                }
                writeln!(
                    writer,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    record.population_id().as_str(),
                    record.block_id().as_str(),
                    record.attempt_index(),
                    record.proposals(),
                    record.accepted(),
                    record.duplicate_rows(),
                    record.bins_outside(),
                    record.signs_failed(),
                    finite_decimal(config, record.max_r2_error())?,
                    finite_decimal(config, record.sum_r2_error())?,
                    finite_decimal(config, record.max_marginal_error())?,
                    record.selected_best(),
                )
                .map_err(|error| io(contract.path, "write_row", &error))?;
            }
            if best != 1 {
                return Err(internal(
                    config,
                    format!(
                        "la unidad {}×{} marca {best} intentos como mejor",
                        unit.population_id().as_str(),
                        unit.block_id().as_str()
                    ),
                ));
            }
        }
        Ok(())
    })
}

// ---------------------------------------------------------------------------
// `f5-founders.tsv` y `truth.jsonl` v2: una sola pasada de asignaciones
// ---------------------------------------------------------------------------

/// Un SNP que una copia lleva, tal como la proyección lo publica tres veces
/// —modificación del segmento, SNP vigente del fragmento e hito del historial—.
struct CopySnp {
    position_1based: u64,
    reference_base: &'static str,
    alternate_base: &'static str,
    site_id: F5SiteId,
    haplotype_id: F5HaplotypeId,
}

/// La copia cuyo registro se está componiendo.
struct OpenCopy {
    population_id: PopulationId,
    individual_id: IndividualId,
    sequence_id: SequenceId,
    copy_index: u32,
    founder_copy_id: FounderCopyId,
}

/// Estado vivo de la pasada única. El único vector es el de los SNP de la
/// copia en curso, que se vacía en cada cambio de copia.
struct CopyStream {
    open: Option<OpenCopy>,
    snps: Vec<CopySnp>,
    written: u64,
}

fn write_founders_and_truth(
    workspace: &mut ExecutionWorkspace,
    request: &F5SuccessRequest<'_>,
    identities: &mut F5FounderIdentities,
    plan: &WritePlan,
) -> Result<(F5ArtifactLink, F5TruthLink), Diagnostic> {
    let contract = &F5_TSV_CONTRACTS[FOUNDERS];
    let mut founders = workspace.begin_streamed_artifact(
        contract.artifact_id,
        contract.artifact_kind,
        Path::new(contract.path),
    )?;
    let mut truth = workspace.begin_streamed_artifact(
        F5_TRUTH_ARTIFACT_ID,
        F5_TRUTH_ARTIFACT_KIND,
        Path::new(F5_TRUTH_PATH),
    )?;
    let mut founders_file = founders.take_file()?;
    let mut truth_file = truth.take_file()?;
    let streamed = {
        let mut tsv = BufWriter::new(&mut founders_file);
        let mut jsonl = BufWriter::new(&mut truth_file);
        (|| -> Result<(), Diagnostic> {
            writeln!(tsv, "{}", contract.header)
                .map_err(|error| io(contract.path, "write_header", &error))?;
            write_truth_header(&mut jsonl, request, plan)?;
            stream_copies(&mut tsv, &mut jsonl, request, identities, plan)?;
            tsv.flush()
                .map_err(|error| io(contract.path, "flush", &error))?;
            stream_genotypes(&mut jsonl, request, identities, plan)?;
            jsonl
                .flush()
                .map_err(|error| io(F5_TRUTH_PATH, "flush", &error))
        })()
    };
    streamed?;
    founders.restore_file(founders_file);
    truth.restore_file(truth_file);
    let founders_record = workspace.finish_streamed_artifact(founders, |_| Ok(()))?;
    let truth_record = workspace.finish_streamed_artifact(truth, |_| Ok(()))?;
    Ok((
        F5ArtifactLink {
            artifact_id: contract.artifact_id,
            artifact_kind: contract.artifact_kind,
            path: contract.path,
            size_bytes: founders_record.size_bytes,
            sha256: founders_record.digest,
        },
        F5TruthLink {
            size_bytes: truth_record.size_bytes,
            sha256: truth_record.digest,
        },
    ))
}

/// Cabecera v2: los campos de DTC-061 más `projection` y `populations`.
fn write_truth_header(
    writer: &mut impl Write,
    request: &F5SuccessRequest<'_>,
    plan: &WritePlan,
) -> Result<(), Diagnostic> {
    let config = request.config;
    let populations = plan
        .individual_counts
        .iter()
        .map(|(population_id, individual_count)| {
            json!({
                "population_id": population_id.as_str(),
                "individual_count": individual_count,
                "ploidy": config.ploidy(),
            })
        })
        .collect::<Vec<_>>();
    let header = json!({
        "record_type": "header",
        "document_kind": TRUTH_DOCUMENT_KIND,
        "schema_version": F5_TRUTH_SCHEMA_VERSION,
        "execution_id": request.execution_id.as_str(),
        "reference_id": config.reference_id(),
        "reference_content_sha256": request.reference_sha256,
        "coordinate_convention": TRUTH_COORDINATE_CONVENTION,
        "counts": {"copies": plan.copy_count, "genotypes": plan.genotype_count},
        "source_digests": [
            {
                "source": source_ref_json(TruthSourceKind::Reference, config.reference_id()),
                "sha256": request.reference_sha256,
            },
            {
                "source": source_ref_json(TruthSourceKind::Matrix, request.profile_id),
                "sha256": request.profile_metadata_sha256,
            },
        ],
        "group_digests": [],
        "projection": {
            "kind": F5_TRUTH_PROJECTION_KIND,
            "projection_version": F5_TRUTH_PROJECTION_VERSION,
            "scenario_id": config.scenario_id(),
        },
        "populations": populations,
    });
    write_record(writer, &header)
}

/// La pasada única: filas de `f5-founders.tsv` y registros `haplotype_copy`.
fn stream_copies(
    tsv: &mut impl Write,
    jsonl: &mut impl Write,
    request: &F5SuccessRequest<'_>,
    identities: &mut F5FounderIdentities,
    plan: &WritePlan,
) -> Result<(), Diagnostic> {
    let config = request.config;
    let layout = request.layout;
    let search = request.search;
    let sequences = request.sequences;
    let reference_id = config.reference_id();
    let profile_id = request.profile_id;
    let mut stream = CopyStream {
        open: None,
        snps: Vec::new(),
        written: 0,
    };
    identities.for_each_assignment(request.assignment, |row: &F5BlockAssignment| {
        if stream
            .open
            .as_ref()
            .is_none_or(|open| &open.founder_copy_id != row.founder_copy_id())
        {
            flush_copy(
                jsonl,
                &mut stream,
                reference_id,
                profile_id,
                sequences,
                config,
            )?;
            stream.open = Some(OpenCopy {
                population_id: row.population_id().clone(),
                individual_id: row.individual_id().clone(),
                sequence_id: row.sequence_id().clone(),
                copy_index: row.copy_index(),
                founder_copy_id: row.founder_copy_id().clone(),
            });
        }
        if u64::from(row.copy_index()) >= u64::from(config.ploidy()) {
            return Err(internal(
                config,
                "una asignación publica un copy_index fuera de 0..ploidy-1",
            ));
        }
        writeln!(
            tsv,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.population_id().as_str(),
            row.individual_id().as_str(),
            row.sequence_id().as_str(),
            row.copy_index(),
            row.block_id().as_str(),
            row.founder_copy_id().as_str(),
            row.block_assignment_id().as_str(),
            row.haplotype_id().as_str(),
        )
        .map_err(|error| io(F5_TSV_CONTRACTS[FOUNDERS].path, "write_row", &error))?;

        let position = *plan
            .search_positions
            .get(&(row.population_id().clone(), row.block_index()))
            .ok_or_else(|| internal(config, "una asignación nombra una unidad desconocida"))?;
        let block = plan
            .block(layout, row.block_index())
            .ok_or_else(|| internal(config, "una asignación nombra un bloque sin geometría"))?;
        let matrix = search.blocks()[position].matrix();
        let rank = usize::try_from(row.haplotype_rank())
            .map_err(|_| internal(config, "un rango haplotípico no cabe en usize"))?;
        collect_copy_snps(
            config,
            matrix,
            block,
            rank,
            row.haplotype_id(),
            &mut stream.snps,
        )
    })?;
    flush_copy(
        jsonl,
        &mut stream,
        reference_id,
        profile_id,
        sequences,
        config,
    )?;
    if stream.written != plan.copy_count {
        return Err(internal(
            config,
            format!(
                "la proyección publicó {} copias y el repositorio materializó {}",
                stream.written, plan.copy_count
            ),
        ));
    }
    Ok(())
}

/// Añade los sitios del bloque donde la fila asignada lleva el alelo
/// alternativo. Un bit a 1 es siempre ALT; la polaridad no lo reinterpreta.
fn collect_copy_snps(
    config: &EffectiveF5Config,
    matrix: &F5AlleleMatrix,
    block: &F5LayoutBlock,
    rank: usize,
    haplotype_id: &F5HaplotypeId,
    snps: &mut Vec<CopySnp>,
) -> Result<(), Diagnostic> {
    if matrix.site_count() != block.sites().len() {
        return Err(internal(
            config,
            "la matriz ganadora no tiene los sitios de su propio bloque",
        ));
    }
    for (site_index, site) in block.sites().iter().enumerate() {
        let alternate = matrix.allele(rank, site_index).ok_or_else(|| {
            internal(config, "una asignación nombra una celda fuera de la matriz")
        })?;
        if alternate {
            snps.push(CopySnp {
                position_1based: site.position_1based(),
                reference_base: canonical_base(config, site.reference_base())?,
                alternate_base: canonical_base(config, site.alternate_base())?,
                site_id: site.site_id().clone(),
                haplotype_id: haplotype_id.clone(),
            });
        }
    }
    Ok(())
}

/// Serializa la copia abierta **por flujo** y vacía su búfer.
///
/// El registro `haplotype_copy` no se materializa: se escriben sus claves —en
/// el orden canónico de `serde_json::Map`, que es un `BTreeMap`— y sus tres
/// vistas recorren el búfer de [`CopySnp`] tres veces emitiendo cada SNP
/// directamente al destino. Ninguna de las tres construye un `Value`, de modo
/// que lo vivo sigue siendo únicamente el búfer que el módulo promete y el
/// pico del escritor no crece con `B`.
fn flush_copy(
    writer: &mut impl Write,
    stream: &mut CopyStream,
    reference_id: &str,
    profile_id: &str,
    sequences: &F5ReferenceSequences,
    config: &EffectiveF5Config,
) -> Result<(), Diagnostic> {
    let Some(copy) = stream.open.take() else {
        return Ok(());
    };
    let length = sequences.length(&copy.sequence_id).ok_or_else(|| {
        configuration(
            config,
            "reference.sequence_lengths",
            format!(
                "la copia de la secuencia {} no puede publicarse sin su longitud",
                copy.sequence_id.as_str()
            ),
        )
    })?;
    // Los bloques llegan en la clave publicada —`block_id`—, que es un digest
    // y no el orden físico. La secuencia de la copia sí es física: sus SNP se
    // ordenan por coordenada antes de serializarse. Es una ordenación del
    // búfer de **una** copia, acotado por los sitios de su secuencia, no de
    // ninguna tabla de filas.
    stream.snps.sort_by_key(|snp| snp.position_1based);
    if stream
        .snps
        .windows(2)
        .any(|pair| pair[0].position_1based == pair[1].position_1based)
    {
        return Err(internal(
            config,
            "dos sitios de una misma copia ocupan la misma coordenada",
        ));
    }
    let phase = u8::try_from(copy.copy_index)
        .map_err(|_| internal(config, "el índice de copia no cabe en la fase publicada"))?;
    let sequence_id = copy.sequence_id.as_str();
    let population_id = copy.population_id.as_str();

    put(writer, br#"{"current_pieces":[{"kind":"source","snps":["#)?;
    for (index, snp) in stream.snps.iter().enumerate() {
        if index > 0 {
            put(writer, b",")?;
        }
        put(writer, br#"{"alternate":"#)?;
        put_string(writer, snp.alternate_base)?;
        put(writer, br#","event_id":"#)?;
        put_string(writer, snp.site_id.as_str())?;
        put(writer, br#","position":"#)?;
        put_number(writer, snp.position_1based - 1)?;
        put(writer, br#","provenance":"#)?;
        put_snp_provenance(writer, snp, population_id, profile_id)?;
        put(writer, br#","reference":"#)?;
        put_string(writer, snp.reference_base)?;
        put(writer, b"}")?;
    }
    put(writer, br#"],"source":"#)?;
    put_segment_source(writer, reference_id, sequence_id, length)?;
    put(writer, br#"}],"founder_membership":{"founder_copy_id":"#)?;
    put_string(writer, copy.founder_copy_id.as_str())?;
    put(writer, br#","population_id":"#)?;
    put_string(writer, population_id)?;
    put(writer, br#"},"history":["#)?;
    for (index, snp) in stream.snps.iter().enumerate() {
        if index > 0 {
            put(writer, b",")?;
        }
        put(writer, br#"{"event_id":"#)?;
        put_string(writer, snp.site_id.as_str())?;
        put(writer, br#","kind":"snp","position":"#)?;
        put_number(writer, snp.position_1based - 1)?;
        put(writer, br#","provenance":"#)?;
        put_snp_provenance(writer, snp, population_id, profile_id)?;
        put(writer, b"}")?;
    }
    put(writer, br#"],"identity":{"copy_index":"#)?;
    put_number(writer, u64::from(copy.copy_index))?;
    put(writer, br#","individual_id":"#)?;
    put_string(writer, copy.individual_id.as_str())?;
    put(writer, br#","reference_id":"#)?;
    put_string(writer, reference_id)?;
    put(writer, br#","sequence_id":"#)?;
    put_string(writer, sequence_id)?;
    put(writer, br#"},"phase":"#)?;
    put_number(writer, u64::from(phase))?;
    put(
        writer,
        br#","provenance":{"origin":"f5_founder","sources":["#,
    )?;
    put_source_ref(writer, TruthSourceKind::Reference, reference_id)?;
    put(writer, b",")?;
    put_source_ref(writer, TruthSourceKind::Matrix, profile_id)?;
    put(
        writer,
        br#"]},"record_type":"haplotype_copy","segments":[{"modifications":["#,
    )?;
    for (index, snp) in stream.snps.iter().enumerate() {
        if index > 0 {
            put(writer, b",")?;
        }
        put(writer, br#"{"alternate":"#)?;
        put_string(writer, snp.alternate_base)?;
        put(writer, br#","event_id":"#)?;
        put_string(writer, snp.site_id.as_str())?;
        put(writer, br#","kind":"snp","position":"#)?;
        put_number(writer, snp.position_1based - 1)?;
        put(writer, br#","provenance":"#)?;
        put_snp_provenance(writer, snp, population_id, profile_id)?;
        put(writer, br#","reference":"#)?;
        put_string(writer, snp.reference_base)?;
        put(writer, b"}")?;
    }
    put(writer, br#"],"source":"#)?;
    put_segment_source(writer, reference_id, sequence_id, length)?;
    put(writer, b"}]}\n")?;

    stream.snps.clear();
    stream.written += 1;
    Ok(())
}

/// Fuente del segmento: la secuencia completa, `[0, sequence_length)`.
fn put_segment_source(
    writer: &mut impl Write,
    reference_id: &str,
    sequence_id: &str,
    length: u64,
) -> Result<(), Diagnostic> {
    put(writer, br#"{"end":"#)?;
    put_number(writer, length)?;
    put(writer, br#","reference_id":"#)?;
    put_string(writer, reference_id)?;
    put(writer, br#","sequence_id":"#)?;
    put_string(writer, sequence_id)?;
    put(writer, br#","sequence_length":"#)?;
    put_number(writer, length)?;
    put(writer, br#","source":"#)?;
    put_source_ref(writer, TruthSourceKind::Reference, reference_id)?;
    put(writer, br#","start":0}"#)
}

/// Procedencia de un SNP fundador: la fila haplotípica de la que salió el
/// alelo, la subpoblación y el `SITE_`.
fn put_snp_provenance(
    writer: &mut impl Write,
    snp: &CopySnp,
    population_id: &str,
    profile_id: &str,
) -> Result<(), Diagnostic> {
    put(writer, br#"{"kind":"matrix","matrix":{"matrix_source_id":"#)?;
    put_string(writer, snp.haplotype_id.as_str())?;
    put(writer, br#","population_id":"#)?;
    put_string(writer, population_id)?;
    put(writer, br#","variant_id":"#)?;
    put_string(writer, snp.site_id.as_str())?;
    put(writer, br#"},"source":"#)?;
    put_source_ref(writer, TruthSourceKind::Matrix, profile_id)?;
    put(writer, br#","version":"1"}"#)
}

fn put_source_ref(
    writer: &mut impl Write,
    kind: TruthSourceKind,
    source_id: &str,
) -> Result<(), Diagnostic> {
    put(writer, br#"{"kind":"#)?;
    put_string(writer, source_kind_str(kind))?;
    put(writer, br#","source_id":"#)?;
    put_string(writer, source_id)?;
    put(writer, b"}")
}

fn put(writer: &mut impl Write, bytes: &[u8]) -> Result<(), Diagnostic> {
    writer
        .write_all(bytes)
        .map_err(|error| io(F5_TRUTH_PATH, "write_record", &error))
}

fn put_number(writer: &mut impl Write, value: u64) -> Result<(), Diagnostic> {
    write!(writer, "{value}").map_err(|error| io(F5_TRUTH_PATH, "write_record", &error))
}

/// Escribe una cadena JSON con el escapado canónico de `serde_json`, sin
/// materializar ninguna copia intermedia.
fn put_string(writer: &mut impl Write, value: &str) -> Result<(), Diagnostic> {
    serde_json::to_writer(&mut *writer, value)
        .map_err(|error| io(F5_TRUTH_PATH, "write_record", &std::io::Error::other(error)))
}

/// Los genotipos: uno por fundador y sitio, ordenados por referencia,
/// población, bloque, `site_id` e individuo.
fn stream_genotypes(
    writer: &mut impl Write,
    request: &F5SuccessRequest<'_>,
    identities: &F5FounderIdentities,
    plan: &WritePlan,
) -> Result<(), Diagnostic> {
    let config = request.config;
    let reference_id = config.reference_id();
    let ploidy = usize::from(config.ploidy());
    let mut written = 0_u64;
    for (position, unit) in identities.units().iter().enumerate() {
        let (block, by_id) = plan
            .block_sites_by_id(request.layout, unit.block_index())
            .ok_or_else(|| internal(config, "una unidad nombra un bloque sin geometría"))?;
        let matrix = request.search.blocks()[plan.unit_to_search[position]].matrix();
        let seats = request.assignment.units()[position].assignments();
        if matrix.site_count() != block.sites().len() {
            return Err(internal(
                config,
                "la matriz ganadora no tiene los sitios de su propio bloque",
            ));
        }
        // «Sitio» en la clave de los genotipos es el `site_id` de `f5_locus`,
        // el único campo de sitio del registro.
        for site_index in by_id.iter().copied() {
            let site = &block.sites()[site_index];
            let reference_base = canonical_base(config, site.reference_base())?;
            let alternate_base = canonical_base(config, site.alternate_base())?;
            let alleles = [
                TruthAllele::Reference,
                TruthAllele::Snp {
                    base: alternate_base.to_owned(),
                },
            ];
            let mut cursor = 0_usize;
            while cursor < seats.len() {
                let individual = seats[cursor].individual_id();
                let mut carriers = Vec::with_capacity(ploidy);
                let mut dosages = [0_u64, 0_u64];
                while cursor < seats.len() && seats[cursor].individual_id() == individual {
                    let seat = &seats[cursor];
                    let rank = usize::try_from(seat.haplotype_rank())
                        .map_err(|_| internal(config, "un rango haplotípico no cabe en usize"))?;
                    let alternate = matrix.allele(rank, site_index).ok_or_else(|| {
                        internal(config, "una asignación nombra una celda fuera de la matriz")
                    })?;
                    let allele_index = usize::from(alternate);
                    dosages[allele_index] += 1;
                    carriers.push(json!({
                        "allele_index": allele_index,
                        "copy": {
                            "reference_id": reference_id,
                            "individual_id": individual.as_str(),
                            "sequence_id": unit.sequence_id().as_str(),
                            "copy_index": seat.copy_index(),
                        },
                    }));
                    cursor += 1;
                }
                if carriers.len() != ploidy {
                    return Err(internal(
                        config,
                        "un individuo no aporta exactamente `ploidy` copias a su bloque",
                    ));
                }
                let record = json!({
                    "record_type": "genotype",
                    "case_id": case_id(
                        unit.population_id().as_str(),
                        individual.as_str(),
                        unit.sequence_id().as_str(),
                        site.position_1based(),
                        &alleles,
                    ),
                    "reference_id": reference_id,
                    "individual_id": individual.as_str(),
                    "case_kind": "snp",
                    "locus": {
                        "sequence_id": unit.sequence_id().as_str(),
                        "start": site.position_1based() - 1,
                        "end": site.position_1based(),
                    },
                    "alleles": [
                        {"kind": "reference"},
                        {"kind": "snp", "base": alternate_base},
                    ],
                    "allele_dosages": dosages,
                    "carriers": carriers,
                    "f5_locus": {
                        "population_id": unit.population_id().as_str(),
                        "block_id": unit.block_id().as_str(),
                        "site_id": site.site_id().as_str(),
                        "ref": reference_base,
                        "alt": alternate_base,
                    },
                });
                write_record(writer, &record)?;
                written += 1;
            }
        }
    }
    if written != plan.genotype_count {
        return Err(internal(
            config,
            format!(
                "la proyección publicó {written} genotipos y la cabecera declara {}",
                plan.genotype_count
            ),
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Utilidades
// ---------------------------------------------------------------------------

fn source_kind_str(kind: TruthSourceKind) -> &'static str {
    match kind {
        TruthSourceKind::Reference => "reference",
        TruthSourceKind::Founder => "founder",
        TruthSourceKind::Explicit => "explicit",
        TruthSourceKind::Random => "random",
        TruthSourceKind::Matrix => "matrix",
    }
}

fn source_ref_json(kind: TruthSourceKind, source_id: &str) -> Value {
    json!({"kind": source_kind_str(kind), "source_id": source_id})
}

/// Escribe un registro JSON Lines con LF final obligatorio.
fn write_record(writer: &mut impl Write, value: &Value) -> Result<(), Diagnostic> {
    serde_json::to_writer(&mut *writer, value)
        .map_err(|error| io(F5_TRUTH_PATH, "write_record", &std::io::Error::other(error)))?;
    writer
        .write_all(b"\n")
        .map_err(|error| io(F5_TRUTH_PATH, "write_record", &error))
}

/// Literal decimal mínimo de ida y vuelta, el mismo que la admisión exige a
/// los perfiles.
fn decimal(value: f64) -> String {
    value.to_string()
}

fn finite_decimal(config: &EffectiveF5Config, value: f64) -> Result<String, Diagnostic> {
    if !value.is_finite() {
        return Err(internal(config, "un f64 publicado no es finito"));
    }
    Ok(decimal(value))
}

fn bounded_decimal(config: &EffectiveF5Config, value: f64) -> Result<String, Diagnostic> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(internal(
            config,
            "una media o error LD publicado no pertenece a [0,1]",
        ));
    }
    Ok(decimal(value))
}

/// Base canónica como literal estático: la proyección la escribe tres veces
/// por SNP y dos por genotipo, y ninguna de esas escrituras necesita asignar.
fn canonical_base(config: &EffectiveF5Config, base: u8) -> Result<&'static str, Diagnostic> {
    match base {
        b'A' => Ok("A"),
        b'C' => Ok("C"),
        b'G' => Ok("G"),
        b'T' => Ok("T"),
        _ => Err(internal(
            config,
            "un sitio publica una base que no es A, C, G ni T",
        )),
    }
}

fn io(path: &str, operation: &str, error: &std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path)
        .with("operation", operation)
        .with("detail", error.to_string())
}

fn configuration(
    config: &EffectiveF5Config,
    fields: &str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::f5_configuration(config.scenario_path().display().to_string(), fields, detail)
}

fn internal(config: &EffectiveF5Config, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::f5_internal(config.scenario_id(), detail)
}
