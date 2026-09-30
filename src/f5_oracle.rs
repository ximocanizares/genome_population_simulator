//! F5.4.2. Oráculo independiente de aceptación sobre los bytes cerrados.
//!
//! Este módulo **no lee ninguna estructura del productor**. Recibe la ruta del
//! paquete publicado —los siete TSV y `outputs/truth.jsonl` v2— y una
//! *declaración* de lo que la ejecución afirma, y comprueba la segunda contra
//! los primeros. No recibe `F5PhysicalLayout`, `F5HaplotypeSearch`,
//! `F5FounderAssignment`, `F5FounderIdentities` ni `F5SuccessArtifacts`: ni un
//! contador, ni una asignación, ni un acumulador del productor entra aquí.
//!
//! ## Qué rederiva por su cuenta
//!
//! - **Identidades DTC-037/DTC-039.** Reimplementa el material canónico,
//!   BLAKE3 `derive_key` y Base32, y **se valida antes de usarse** contra los
//!   cinco vectores congelados de `catalogs/id-catalog-v4.tsv`. Con eso
//!   rederiva `BLK_`, `SITE_`, `HAP_`, `FCP_` y `FHA_` desde las columnas
//!   publicadas.
//! - **Subsemillas DTC-089.** Reimplementa el material de subsemilla, el
//!   ChaCha12 local, el rechazo por muestreo de `f5_uniform_index`, el U53 de
//!   `f5_uniform_53` y el muestreador de recuentos `poisson-count-v1`; se
//!   valida contra `catalogs/f5-rng-catalog-v1.tsv` y
//!   `catalogs/f5-rng-catalog-v2.tsv`, y rederiva **las cuatro decisiones cuyos
//!   componentes los bytes publican**: el alelo alternativo de cada sitio
//!   (`site_alt`), su marginal objetivo (`marginal_spectrum`), su polaridad
//!   (`polarity`) y el recuento haplotípico de cada bloque×subpoblación
//!   (`block_haplotype_count`, solo cuando el perfil declara su distribución).
//!   Las demás —longitudes, posiciones y huecos del layout, columnas
//!   iniciales, propuestas de intercambio, muestreo fundador, F, desempates de
//!   cuota y permutación controlada— no publican alguno de sus componentes y
//!   quedan fuera.
//! - **Aritmética DTC-089.** Neumaier propio para `p_hat`, `p_hat_ij`, `D`,
//!   `r²`, medias por bin y errores, con la misma agrupación literal que fija
//!   la decisión, de modo que la comparación es bit a bit.
//! - **DTC-049.** Recalcula suma, delta, normalización, suma y delta de
//!   salida de los cuatro conjuntos declarados.
//! - **MAN-70.** [`recount_f5_man70`] abre **solo** `truth.jsonl`, agrupa por
//!   `site_id`×`population_id`, resuelve cada portador contra la copia
//!   publicada y recalcula las dosis. No abre ningún TSV ni consulta ningún
//!   manifiesto.
//!
//! ## Qué retiene, y por qué cabe
//!
//! El oráculo no materializa el paquete: recorre cada fichero por flujo y
//! cruza `f5-founders.tsv` con las copias de `truth.jsonl` **en paso
//! acompasado**, porque ambas llegan en la misma clave. Lo único que retiene
//! es lo que también retiene el productor: la geometría (`B` bloques y sus
//! `B*S` sitios), las matrices ganadoras en bits (`B*N_pop*H*S`), un rango
//! haplotípico por asiento (`B*C_total`) y un vector de dosis por locus
//! (`B*N_pop*S*(P+1)`). El recuento MAN-70, que por contrato no puede
//! compartir nada con lo anterior, retiene además una entrada por copia con
//! su subpoblación y sus sitios alternativos —`C_total`, el mismo orden que el
//! repositorio de identidades del productor—.
//!
//! La única concesión sobre el presupuesto sellado es el **registro que tiene
//! en vuelo**: el escritor emite cada `haplotype_copy` por flujo y la fórmula
//! lo contabiliza en `output_bytes` a través de `truth_bound`, pero un lector
//! que deba interpretarlo tiene que materializarlo una vez. Esa concesión está
//! acotada por el registro más largo, no por el paquete, y una prueba mide el
//! pico real contra `working_memory_bytes` más dieciséis veces ese registro,
//! además de comprobar que dieciséis veces los genotipos no duplican el pico.
//!
//! ## Interpretaciones del contrato, declaradas
//!
//! - **La clave de cada fichero se comprueba literal.** El oráculo no deriva
//!   el orden esperado de otro fichero: exige que la tupla de columnas que
//!   DTC-089 nombra crezca estrictamente. `block_index` es la única excepción
//!   posible, porque ninguna columna lo publica: se recupera del ordinal de
//!   fila de `f5-block-layout.tsv`, y de él cuelga la rederivación de `BLK_`.
//! - **`expected_*` de MAN-70 usa la frecuencia realizada.** `realized_alt_af`
//!   es `Σ_h q_h x_hi`, es decir la probabilidad exacta de que una copia
//!   muestreada probabilísticamente lleve el alelo alternativo. El modelo
//!   esperado es por tanto binomial en `hwe`, la mezcla de identidad por
//!   descendencia en `diploid_f` y la hipergeométrica inducida por permutar
//!   sin reemplazo las `C` copias en `controlled_hwe_conditioned_counts`. En
//!   este último el recuento alélico fijado `K` no es `realized_alt_af*C`
//!   —la cuota reparte `n_h` enteros—, sino el número de copias de la unidad
//!   cuya fila lleva el alelo alternativo, que los bytes publican.
//! - **`f5-attempts.tsv` publica el intento del que salió la matriz
//!   ganadora.** Su fila `selected_best` da el `attempt_index` con el que se
//!   derivó la polaridad de cada sitio —el único componente RNG que no vive en
//!   una columna—, y sin él esa rederivación sería imposible desde los bytes. Sus seis contadores —filas repetidas, bins
//!   fuera, signos fallidos y los tres errores— son funciones puras de esa
//!   matriz, de sus marginales y de sus bins, y el oráculo los recalcula.
//! - **`json_rows_bound` no acota las líneas de `truth.jsonl`.** Sus términos
//!   son los del efectivo y el manifiesto —censo, muestreo, F realizado,
//!   resúmenes MAN-70 y filas de normalización—, y su término de locus es
//!   `B*N_pop*S_bound`, sin el censo. Las líneas de la proyección crecen con
//!   `individual_count` y viven en `output_bytes`: con censo 4 y ocho bloques
//!   la proyección publica 105 líneas y el bound vale 104. El oráculo acota
//!   por tanto las filas JSON que el manifiesto publicará, y los bytes de la
//!   proyección por `output_bytes`.
//! - **El orden FASTA del layout se aproxima por primera aparición.** La clave
//!   de `f5-block-layout.tsv` empieza por «FASTA», y ningún documento
//!   publicado lleva el diccionario ordenado de secuencias: la admisión sella
//!   `retained_sequence_ids`, que es un recuento. El oráculo numera las
//!   secuencias por el orden en que aparecen en el layout, de modo que un
//!   layout que las enumerase en otro orden —pero consistente consigo mismo—
//!   pasaría. Es una frontera, no una comprobación, y cerrarla exige llevar la
//!   lista ordenada de `sequence_id` a la declaración y atarla al FASTA
//!   sellado.
//! - **La auditoría DTC-049, el espectro, los bins LD y los pesos se
//!   comprueban contra sí mismos.** El oráculo rehace toda la aritmética desde
//!   los `input_values` que la propia declaración trae; no vuelve a abrir los
//!   cinco TSV del perfil. Lo que ata esas cifras al perfil real es
//!   `profile_metadata_sha256`, que el oráculo sí cruza contra
//!   `source_digests` de la verdad publicada. Una declaración con un perfil
//!   inventado **y** su huella inventada pasaría la aritmética y moriría en el
//!   cruce; una con la huella correcta y cifras inventadas pasaría. Cerrarla
//!   exige releer los cinco TSV.
//! - **Comparte con el productor los literales de contrato.**
//!   `F5_TSV_CONTRACTS` —identificadores, clases, rutas y cabeceras— y los
//!   `F5_TRUTH_*` vienen de `f5_artifacts`. Son el contrato que se juzga, no
//!   estado: si el escritor los cambiara, el oráculo cambiaría con él y
//!   ninguno de los dos lo notaría. Se declara porque es la única dependencia
//!   textual entre los dos módulos.
//! - **Dos bloques de una secuencia no pueden solaparse.** DTC-089 escribe la
//!   clave del layout como «FASTA, inicio, fin» y el orden estricto no
//!   prohíbe por sí solo que un bloque termine dentro del siguiente. El
//!   oráculo lo prohíbe porque de ello depende que `case_id` identifique un
//!   caso: dos loci de un mismo individuo en la misma coordenada darían el
//!   mismo `case_id`.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader, Read};
use std::path::Path;

use rand_chacha::ChaCha12Rng;
use rand_chacha::rand_core::{Rng, SeedableRng};
use serde_json::{Map, Value};
use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;
use crate::f5_artifacts::{
    F5_TRUTH_ARTIFACT_ID, F5_TRUTH_ARTIFACT_KIND, F5_TRUTH_PATH, F5_TRUTH_PROJECTION_KIND,
    F5_TRUTH_PROJECTION_VERSION, F5_TRUTH_SCHEMA_VERSION, F5_TSV_CONTRACTS, F5_TSV_SCHEMA_VERSION,
    f5_canonical_decimal,
};
use crate::f5_profile::{
    F5_PROFILE_NORMALIZATION_ALGORITHM, F5_PROFILE_NORMALIZATION_TOLERANCE, F5NamedProbability,
    F5NormalizationAudit, F5NormalizationSet,
};

// ---------------------------------------------------------------------------
// Contratos literales del oráculo
// ---------------------------------------------------------------------------

/// Algoritmo del oráculo, tal como el manifiesto lo nombrará en F5.4.3.
pub const F5_ORACLE_ALGORITHM: &str = "f5-acceptance-oracle-v1";
/// `recount_source` exacto de REQ-MAN-70.
pub const F5_MAN70_RECOUNT_SOURCE: &str = "outputs/truth.jsonl#f5-founder-truth-projection-v2";
/// Catálogo de identidades del que salen los vectores congelados.
pub const F5_ORACLE_ID_CATALOG: &str = "id-catalog-v4";
/// Catálogos de subsemillas de los que salen los vectores congelados.
///
/// Un paquete declara v1 o v2, y **cuál de los dos no es opcional**: v2
/// exactamente cuando el perfil declara la distribución del recuento
/// haplotípico, porque solo entonces puede abrirse `block_haplotype_count`.
/// El oráculo valida los dos ficheros congelados antes de usar su propia
/// criptografía, de modo que rederiva los trece vectores aunque el paquete
/// que juzga use solo doce.
pub const F5_ORACLE_RNG_CATALOG: &str = "f5-rng-catalog-v1";
pub const F5_ORACLE_RNG_CATALOG_V2: &str = "f5-rng-catalog-v2";
/// Decisiones exactas de `f5-rng-catalog-v1`: las tres que este módulo
/// rederiva y las nueve que declara creídas.
pub const F5_RNG_CATALOG_DECISIONS: usize = 12;
/// Decisiones exactas de `f5-rng-catalog-v2`: las doce de v1 y
/// `block_haplotype_count`, que este módulo también rederiva.
pub const F5_RNG_CATALOG_V2_DECISIONS: usize = 13;
/// Los nueve objetos fijos que `json_rows_bound` reserva al manifiesto.
const F5_MANIFEST_FIXED_OBJECTS: u64 = 9;

const ID_CONTEXT: &str = "genome_population_simulator/v1/id";
const SUBSEED_CONTEXT: &str = "genome_population_simulator/v1/subseed";
const DERIVATION_VERSION: u32 = 1;
const BASE32_ALPHABET: &[u8; 32] = b"ABCDEFGHIJKLMNOPQRSTUVWXYZ234567";
const FIFTY_THREE_BIT_BOUND: u64 = 1 << 53;
const NOT_AVAILABLE: &str = "NA";

const ID_CATALOG_TSV: &str = include_str!("../catalogs/id-catalog-v4.tsv");
const RNG_CATALOG_TSV: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");
const RNG_CATALOG_V2_TSV: &str = include_str!("../catalogs/f5-rng-catalog-v2.tsv");

const BLOCK_LAYOUT: usize = 0;
const SITES: usize = 1;
const HAPLOTYPES: usize = 2;
const FREQUENCIES: usize = 3;
const LD_EVALUATION: usize = 4;
const FOUNDERS: usize = 5;
const ATTEMPTS: usize = 6;

// ---------------------------------------------------------------------------
// Causas de rechazo
// ---------------------------------------------------------------------------

/// Causa exacta por la que el oráculo rechaza un paquete.
///
/// El contrato de F5.4.2 no se conforma con que una corrupción falle: exige
/// que falle **por su causa**. Por eso el rechazo lleva una causa tipada y no
/// solo un diagnóstico, y por eso la batería contrasta la causa de cada
/// mutación en lugar de limitarse a comprobar que hubo error.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum F5OracleCause {
    /// Falta un artefacto del paquete o no puede abrirse.
    ArtifactMissing,
    /// Codificación, cabecera, columnas, BOM, CR, comillas o LF final.
    ByteContract,
    /// Ruta, tamaño o SHA-256 declarados que no son los de los bytes.
    Link,
    /// La semilla declarada no es una raíz de 256 bits en hex minúscula.
    Seed,
    /// Geometría u orden de `f5-block-layout.tsv`.
    BlockLayout,
    /// Índice, identidad, coordenada, alelos u orden de `f5-sites.tsv`.
    Site,
    /// Rango, `q`, cadena de alelos u orden de `f5-haplotypes.tsv`.
    Haplotype,
    /// `target_alt_af`, `realized_alt_af`, `abs_error` u orden marginal.
    Marginal,
    /// La polaridad publicada no es la que el RNG normativo deriva.
    Polarity,
    /// Bins, `pair_count`, `mean_r2`, `abs_error`, `NA` u orden LD.
    LinkageDisequilibrium,
    /// La regla de signos de DTC-089 sobre las parejas evaluables.
    SignRule,
    /// Contadores, límites, orden o mejor intento de `f5-attempts.tsv`.
    Attempt,
    /// Modo, cuota controlada o tolerancia efectiva de una unidad.
    Sampling,
    /// F declarado, su ámbito de ploidía o el F realizado publicado.
    Inbreeding,
    /// Un identificador automático no se rederiva desde sus componentes.
    Identity,
    /// Una fila fundadora no describe la asignación que los bytes soportan.
    Assignment,
    /// El oráculo de identidades o de subsemillas no reproduce su catálogo.
    Randomness,
    /// Estructura, orden, procedencia o recuentos de `truth.jsonl` v2.
    Truth,
    /// El recuento MAN-70 desde la verdad no cuadra con lo realizado.
    Man70,
    /// La auditoría DTC-049 declarada no es la que sus valores producen.
    Normalization,
    /// Una cota o presupuesto sellado que el paquete publicado excede.
    ResourceBound,
}

impl F5OracleCause {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::ArtifactMissing => "artifact_missing",
            Self::ByteContract => "byte_contract",
            Self::Link => "link",
            Self::Seed => "seed",
            Self::BlockLayout => "block_layout",
            Self::Site => "site",
            Self::Haplotype => "haplotype",
            Self::Marginal => "marginal",
            Self::Polarity => "polarity",
            Self::LinkageDisequilibrium => "linkage_disequilibrium",
            Self::SignRule => "sign_rule",
            Self::Attempt => "attempt",
            Self::Sampling => "sampling",
            Self::Inbreeding => "inbreeding",
            Self::Identity => "identity",
            Self::Assignment => "assignment",
            Self::Randomness => "randomness",
            Self::Truth => "truth",
            Self::Man70 => "man70",
            Self::Normalization => "normalization",
            Self::ResourceBound => "resource_bound",
        }
    }
}

/// Rechazo del oráculo: la causa tipada y el diagnóstico público.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5OracleRejection {
    cause: F5OracleCause,
    diagnostic: Diagnostic,
}

impl F5OracleRejection {
    #[must_use]
    pub const fn cause(&self) -> F5OracleCause {
        self.cause
    }
    #[must_use]
    pub const fn diagnostic(&self) -> &Diagnostic {
        &self.diagnostic
    }
    #[must_use]
    pub fn into_diagnostic(self) -> Diagnostic {
        self.diagnostic
    }
}

impl From<F5OracleRejection> for Diagnostic {
    fn from(rejection: F5OracleRejection) -> Self {
        rejection.diagnostic
    }
}

/// La anotación `#[track_caller]` y la llamada a `record` son el barrido de
/// guardas descrito en [`crate::guard_trace`]: apagado no cuesta nada, y
/// encendido anota **desde qué línea** se rechazó, no esta.
#[track_caller]
fn reject(
    cause: F5OracleCause,
    path: &str,
    field: &str,
    detail: impl Into<String>,
) -> F5OracleRejection {
    crate::guard_trace::record();
    let diagnostic = Diagnostic::format(
        "GPS-FMT-101",
        "El paquete de perfiles F5 no satisface su formato, versión, orden o identidad declarada.",
    )
    .with("path", path)
    .with("field", field)
    .with("cause", cause.as_str())
    .with("detail", detail.into());
    F5OracleRejection { cause, diagnostic }
}

#[track_caller]
fn resource_reject(path: &str, detail: impl Into<String>) -> F5OracleRejection {
    crate::guard_trace::record();
    let diagnostic = Diagnostic::resource(
        "GPS-RES-101",
        "Las dimensiones o presupuestos F5 exceden el límite declarado o duro.",
    )
    .with("scenario_id", path)
    .with("cause", F5OracleCause::ResourceBound.as_str())
    .with("detail", detail.into());
    F5OracleRejection {
        cause: F5OracleCause::ResourceBound,
        diagnostic,
    }
}

#[track_caller]
fn catalog_reject(detail: impl Into<String>) -> F5OracleRejection {
    crate::guard_trace::record();
    let diagnostic = Diagnostic::internal(
        "GPS-INTERNAL-101",
        "Una operación F5 posterior a un preflight válido produjo overflow o un resultado numérico no finito pese a sus cotas verificadas.",
    )
    .with("scenario_id", F5_ORACLE_ALGORITHM)
    .with("cause", F5OracleCause::Randomness.as_str())
    .with("detail", detail.into());
    F5OracleRejection {
        cause: F5OracleCause::Randomness,
        diagnostic,
    }
}

type Checked<T> = Result<T, F5OracleRejection>;

// ---------------------------------------------------------------------------
// Declaración bajo examen
// ---------------------------------------------------------------------------

/// Modo de muestreo declarado para una unidad bloque×subpoblación.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5DeclaredMode {
    Probabilistic,
    Controlled,
}

impl F5DeclaredMode {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Probabilistic => "probabilistic",
            Self::Controlled => "controlled",
        }
    }
}

/// Subpoblación declarada: censo, ploidía efectiva y F.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredSubpopulation {
    pub population_id: String,
    pub individual_count: u64,
    /// F **efectivo**: `None` fuera de ploidía 2 y `Some` —posiblemente `0`—
    /// en diploides.
    pub inbreeding_f: Option<f64>,
}

/// Fila declarada de `sampling`: modo, origen y tolerancias de una unidad.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredSampling {
    pub population_id: String,
    pub block_id: String,
    pub requested_mode: F5DeclaredMode,
    pub effective_mode: F5DeclaredMode,
    pub controlled_max_abs_error_declared: Option<f64>,
    pub controlled_max_abs_error_effective: Option<f64>,
}

/// Fila declarada de `f_realized`.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredRealizedF {
    pub population_id: String,
    pub block_id: String,
    pub target_f: Option<f64>,
    /// `evaluated` o `not_applicable`.
    pub status: String,
    pub value: Option<f64>,
}

/// Banda declarada del espectro `f5_snp_spectrum`, con su peso ya
/// normalizado por DTC-049.
///
/// `alt_af_min` y `weight` se contrastan contra la auditoría DTC-049 del
/// efectivo, que publica la primera como clave y el segundo como valor
/// efectivo. `alt_af_max` **no** aparece en ningún documento publicado: viene
/// del perfil, que la declaración ata por `profile_metadata_sha256` y el
/// oráculo cruza contra `source_digests` de la verdad publicada. Es la misma
/// frontera que ya tienen `ld_bins` y `haplotype_weights`, y queda declarada.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5DeclaredSpectrumBand {
    pub alt_af_min: f64,
    pub alt_af_max: f64,
    pub weight: f64,
}

/// Bin LD declarado por el perfil para una subpoblación.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredLdBin {
    pub population_id: String,
    pub distance_start_bp: u64,
    pub distance_end_bp: u64,
    pub target_mean_r2: f64,
    pub max_abs_error: f64,
}

/// Enlace declarado de uno de los siete TSV.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredLink {
    pub artifact_id: String,
    pub artifact_kind: String,
    pub schema_version: String,
    pub path: String,
    pub size_bytes: u64,
    pub sha256: String,
}

/// Enlace declarado de la proyección fundadora.
#[derive(Clone, Debug, PartialEq)]
pub struct F5DeclaredTruthLink {
    pub artifact_id: String,
    pub artifact_kind: String,
    pub schema_version: String,
    pub projection_version: String,
    pub path: String,
    pub size_bytes: u64,
    pub sha256: String,
}

/// Cotas y presupuestos que la admisión selló y el paquete debe respetar.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5DeclaredBounds {
    pub block_count: u64,
    /// Bases elegibles que la admisión selló. El oráculo no abre la
    /// referencia, de modo que no puede recomputarlas; sí puede exigir que no
    /// sean menos que las que el layout colocó, que es la cota inferior
    /// observable desde `f5-block-layout.tsv`.
    pub eligible_bases: u64,
    pub site_bound_per_block: u64,
    pub bin_count_total: u64,
    pub max_haplotype_attempts: u64,
    pub max_swap_proposals: u64,
    pub tsv_rows_bound: u64,
    pub json_rows_bound: u64,
    pub output_bytes: u64,
    pub output_limit_bytes: u64,
}

/// Todo lo que la ejecución declara y el oráculo contrasta contra los bytes.
///
/// Ninguno de estos campos es un acumulador del productor: son las cifras que
/// el manifiesto publicará y que, si el paquete no las sostiene, deben ser
/// rechazadas. Construirlas es responsabilidad del llamador.
#[derive(Clone, Debug, PartialEq)]
pub struct F5PackageDeclaration {
    pub scenario_id: String,
    pub execution_id: String,
    pub reference_id: String,
    pub reference_sha256: String,
    pub seed: String,
    pub ploidy: u8,
    pub profile_id: String,
    pub profile_metadata_sha256: String,
    pub haplotype_weights: Vec<f64>,
    /// Distribución del recuento haplotípico que el perfil declara, o `None`
    /// cuando H es el número de pesos en todas las unidades.
    pub haplotype_count_distribution: Option<F5DeclaredHaplotypeCountDistribution>,
    /// H publicado de cada bloque×subpoblación, en el mismo orden canónico que
    /// `sampling`. Solo se declara cuando hay distribución: el oráculo **no**
    /// se lo cree, lo rederiva del ámbito `block_haplotype_count` y lo compara
    /// contra esta lista y contra las filas de `f5-haplotypes.tsv`.
    pub haplotype_counts: Vec<F5DeclaredHaplotypeCount>,
    pub snp_spectrum: Vec<F5DeclaredSpectrumBand>,
    pub subpopulations: Vec<F5DeclaredSubpopulation>,
    pub sampling: Vec<F5DeclaredSampling>,
    pub f_realized: Vec<F5DeclaredRealizedF>,
    pub ld_bins: Vec<F5DeclaredLdBin>,
    pub normalizations: Vec<F5NormalizationAudit>,
    pub artifacts: Vec<F5DeclaredLink>,
    pub truth: F5DeclaredTruthLink,
    pub bounds: F5DeclaredBounds,
    pub id_catalog_id: String,
    pub rng_catalog_id: String,
    /// `attempt_summary` del manifiesto. Dos de sus tres términos son
    /// funciones puras de `f5-attempts.tsv` y aquí se contrastan con él;
    /// `layout_attempts_used` no es observable desde ningún TSV publicado y
    /// queda declarado.
    pub attempt_summary: F5DeclaredAttemptSummary,
}

/// Distribución declarada del recuento haplotípico por bloque.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F5DeclaredHaplotypeCountDistribution {
    pub lambda: f64,
    pub minimum: u64,
    pub maximum: u64,
}

/// Recuento haplotípico publicado de un bloque×subpoblación.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5DeclaredHaplotypeCount {
    pub population_id: String,
    pub block_id: String,
    pub haplotype_count: u64,
}

/// `attempt_summary` declarado por el manifiesto.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5DeclaredAttemptSummary {
    pub layout_attempts_used: u64,
    pub haplotype_attempts_total: u64,
    pub max_attempt_index: u64,
}

// ---------------------------------------------------------------------------
// Informe
// ---------------------------------------------------------------------------

/// Modelo genotípico de un locus, con los tres literales de DTC-089.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5GenotypeModel {
    Hwe,
    DiploidF,
    ControlledHweConditionedCounts,
}

impl F5GenotypeModel {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Hwe => "hwe",
            Self::DiploidF => "diploid_f",
            Self::ControlledHweConditionedCounts => "controlled_hwe_conditioned_counts",
        }
    }
}

/// Resumen MAN-70 de un locus×subpoblación, con los siete vectores por dosis.
///
/// Siete, no nueve: esperados —probabilidades, frecuencias y recuentos—,
/// realizados —frecuencias y recuentos— y recontados —frecuencias y recuentos—.
/// Cada uno tiene `ploidy+1` entradas, que es lo que hace que el sobre de este
/// registro en `f5-resource-accounting-v9` crezca con la ploidía.
#[derive(Clone, Debug, PartialEq)]
pub struct F5InitializationSummary {
    population_id: String,
    block_id: String,
    site_id: String,
    ploidy: u8,
    individual_count: u64,
    model: F5GenotypeModel,
    requested_mode: F5DeclaredMode,
    effective_mode: F5DeclaredMode,
    input_alt_af: f64,
    realized_alt_af: f64,
    expected_genotype_probabilities: Vec<f64>,
    expected_genotype_frequencies: Vec<f64>,
    expected_genotype_counts: Vec<f64>,
    realized_genotype_frequencies: Vec<f64>,
    realized_genotype_counts: Vec<u64>,
    recounted_genotype_frequencies: Vec<f64>,
    recounted_genotype_counts: Vec<u64>,
    recount_matches: bool,
}

impl F5InitializationSummary {
    #[must_use]
    pub fn population_id(&self) -> &str {
        &self.population_id
    }
    #[must_use]
    pub fn block_id(&self) -> &str {
        &self.block_id
    }
    #[must_use]
    pub fn site_id(&self) -> &str {
        &self.site_id
    }
    #[must_use]
    pub const fn ploidy(&self) -> u8 {
        self.ploidy
    }
    #[must_use]
    pub const fn individual_count(&self) -> u64 {
        self.individual_count
    }
    #[must_use]
    pub const fn model(&self) -> F5GenotypeModel {
        self.model
    }
    #[must_use]
    pub const fn requested_mode(&self) -> F5DeclaredMode {
        self.requested_mode
    }
    #[must_use]
    pub const fn effective_mode(&self) -> F5DeclaredMode {
        self.effective_mode
    }
    #[must_use]
    pub const fn input_alt_af(&self) -> f64 {
        self.input_alt_af
    }
    #[must_use]
    pub const fn realized_alt_af(&self) -> f64 {
        self.realized_alt_af
    }
    #[must_use]
    pub fn expected_genotype_probabilities(&self) -> &[f64] {
        &self.expected_genotype_probabilities
    }
    #[must_use]
    pub fn expected_genotype_frequencies(&self) -> &[f64] {
        &self.expected_genotype_frequencies
    }
    #[must_use]
    pub fn expected_genotype_counts(&self) -> &[f64] {
        &self.expected_genotype_counts
    }
    #[must_use]
    pub fn realized_genotype_frequencies(&self) -> &[f64] {
        &self.realized_genotype_frequencies
    }
    #[must_use]
    pub fn realized_genotype_counts(&self) -> &[u64] {
        &self.realized_genotype_counts
    }
    #[must_use]
    pub fn recounted_genotype_frequencies(&self) -> &[f64] {
        &self.recounted_genotype_frequencies
    }
    #[must_use]
    pub fn recounted_genotype_counts(&self) -> &[u64] {
        &self.recounted_genotype_counts
    }
    /// `recount_source` obligatorio de REQ-MAN-70.
    #[must_use]
    pub const fn recount_source(&self) -> &'static str {
        F5_MAN70_RECOUNT_SOURCE
    }
    /// Siempre `true` en un informe devuelto: la desigualdad es un rechazo.
    #[must_use]
    pub const fn recount_matches(&self) -> bool {
        self.recount_matches
    }
}

/// Cotas que el paquete publicado consumió, ya contrastadas con las selladas.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5ResourceAudit {
    tsv_rows: u64,
    observable_json_rows: u64,
    truth_records: u64,
    output_bytes: u64,
}

impl F5ResourceAudit {
    #[must_use]
    pub const fn tsv_rows(&self) -> u64 {
        self.tsv_rows
    }
    /// Filas JSON del efectivo y del manifiesto que **son observables desde
    /// los bytes publicados**: fundadores, subpoblaciones, muestreo, F
    /// realizado, resúmenes MAN-70, filas de normalización y los nueve objetos
    /// fijos. No incluye los `2*B*N_pop*H` rangos de `quota_ties`, que ningún
    /// TSV enumera, de modo que es una **cota inferior** de lo que
    /// `json_rows_bound` acota.
    #[must_use]
    pub const fn observable_json_rows(&self) -> u64 {
        self.observable_json_rows
    }
    /// Registros de `truth.jsonl`, que `output_bytes` acota por bytes y
    /// `json_rows_bound` **no** acota por filas.
    #[must_use]
    pub const fn truth_records(&self) -> u64 {
        self.truth_records
    }
    #[must_use]
    pub const fn output_bytes(&self) -> u64 {
        self.output_bytes
    }
}

/// Lo que el oráculo devuelve cuando el paquete satisface su contrato entero.
#[derive(Clone, Debug, PartialEq)]
pub struct F5OracleReport {
    initialization_summaries: Vec<F5InitializationSummary>,
    resources: F5ResourceAudit,
    copies: u64,
    genotypes: u64,
}

impl F5OracleReport {
    /// Resúmenes MAN-70 ordenados por subpoblación, bloque y sitio.
    #[must_use]
    pub fn initialization_summaries(&self) -> &[F5InitializationSummary] {
        &self.initialization_summaries
    }
    #[must_use]
    pub const fn resources(&self) -> &F5ResourceAudit {
        &self.resources
    }
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
    #[must_use]
    pub const fn genotypes(&self) -> u64 {
        self.genotypes
    }
}

// ---------------------------------------------------------------------------
// Suma de Neumaier independiente
// ---------------------------------------------------------------------------

/// Suma compensada `neumaier-f64-v1`, reimplementada aquí.
#[derive(Clone, Copy, Debug, Default)]
struct Neumaier {
    sum: f64,
    compensation: f64,
}

impl Neumaier {
    fn add(&mut self, value: f64) {
        let next = self.sum + value;
        self.compensation += if self.sum.abs() >= value.abs() {
            (self.sum - next) + value
        } else {
            (value - next) + self.sum
        };
        self.sum = next;
    }

    const fn total(self) -> f64 {
        self.sum + self.compensation
    }
}

fn neumaier_total(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut accumulator = Neumaier::default();
    for value in values {
        accumulator.add(value);
    }
    accumulator.total()
}

// ---------------------------------------------------------------------------
// Oráculo de identidades DTC-037/DTC-039
// ---------------------------------------------------------------------------

fn push_component(material: &mut Vec<u8>, name: &str, value: &str) -> Checked<()> {
    // Los dos rechazos siguientes son inalcanzables: exigen un nombre o un valor
    // de más de 4 GiB, que es lo que cabe en el `u32` del formato del material.
    // Los nombres son literales de este módulo y los valores son identificadores,
    // cifras o bases que los lectores ya han materializado línea a línea; ninguna
    // prueba construye una cadena así. Se conservan porque una conversión que
    // truncara en silencio cambiaría el identificador derivado.
    let name_length =
        u32::try_from(name.len()).map_err(|_| catalog_reject("un nombre de componente no cabe"))?;
    let value_length =
        u32::try_from(value.len()).map_err(|_| catalog_reject("un valor de componente no cabe"))?;
    material.extend_from_slice(&name_length.to_be_bytes());
    material.extend_from_slice(name.as_bytes());
    material.extend_from_slice(&value_length.to_be_bytes());
    material.extend_from_slice(value.as_bytes());
    Ok(())
}

fn canonical_material(root: &[u8; 32], components: &[(&str, String)]) -> Checked<Vec<u8>> {
    // Inalcanzable: `components` lo arma cada llamador de este módulo con menos de
    // diez entradas, o lo lee de una fila del catálogo compilado, y solo una
    // cuenta por encima de `u32::MAX` haría fallar la conversión. Se conserva por
    // el mismo motivo que los dos rechazos de `push_component`.
    let count = u32::try_from(components.len())
        .map_err(|_| catalog_reject("el recuento de componentes no cabe"))?;
    let mut material = Vec::new();
    material.extend_from_slice(root);
    material.extend_from_slice(&DERIVATION_VERSION.to_be_bytes());
    material.extend_from_slice(&count.to_be_bytes());
    for (name, value) in components {
        push_component(&mut material, name, value)?;
    }
    Ok(material)
}

fn base32(bytes: &[u8]) -> String {
    let bits = bytes.len() * 8;
    (0..bits.div_ceil(5))
        .map(|group| {
            let mut index = 0_usize;
            for offset in 0..5 {
                let bit = group * 5 + offset;
                let set = bit < bits && (bytes[bit / 8] >> (7 - bit % 8)) & 1 == 1;
                index = (index << 1) | usize::from(set);
            }
            char::from(BASE32_ALPHABET[index])
        })
        .collect()
}

fn decode_hex(text: &str) -> Option<Vec<u8>> {
    if !text.len().is_multiple_of(2) {
        return None;
    }
    text.as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            let text = std::str::from_utf8(pair).ok()?;
            u8::from_str_radix(text, 16).ok()
        })
        .collect()
}

fn is_lowercase_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

/// Rederivador de los cinco identificadores automáticos F5.
struct IdOracle {
    root: [u8; 32],
    execution_id: String,
    scenario_id: String,
    reference_id: String,
}

impl IdOracle {
    /// Construye el oráculo **después** de reproducir los cinco vectores
    /// congelados del catálogo: si la reimplementación de material, BLAKE3 o
    /// Base32 no los reproduce, no puede usarse para juzgar a nadie.
    fn new(
        root: [u8; 32],
        execution_id: &str,
        scenario_id: &str,
        reference_id: &str,
    ) -> Checked<Self> {
        validate_id_catalog()?;
        Ok(Self {
            root,
            execution_id: execution_id.to_owned(),
            scenario_id: scenario_id.to_owned(),
            reference_id: reference_id.to_owned(),
        })
    }

    fn head(&self) -> Vec<(&'static str, String)> {
        vec![
            ("execution_id", self.execution_id.clone()),
            ("scenario_id", self.scenario_id.clone()),
            ("reference_id", self.reference_id.clone()),
        ]
    }

    fn derive(
        &self,
        prefix: &str,
        class: &str,
        tail: Vec<(&'static str, String)>,
    ) -> Checked<String> {
        let mut components = vec![("class", class.to_owned())];
        components.extend(self.head());
        components.extend(tail);
        let material = canonical_material(&self.root, &components)?;
        let digest = blake3::derive_key(ID_CONTEXT, &material);
        Ok(format!("{prefix}{}", base32(&digest[..16])))
    }

    fn block(&self, block_index: u64) -> Checked<String> {
        self.derive(
            "BLK_",
            "f5_block",
            vec![("block_index_u64", block_index.to_string())],
        )
    }

    fn site(&self, block_index: u64, site_index: u64) -> Checked<String> {
        self.derive(
            "SITE_",
            "f5_site",
            vec![
                ("block_index_u64", block_index.to_string()),
                ("site_index_u64", site_index.to_string()),
            ],
        )
    }

    fn haplotype(
        &self,
        population_id: &str,
        block_index: u64,
        rank: u64,
        alleles: &str,
    ) -> Checked<String> {
        self.derive(
            "HAP_",
            "f5_haplotype",
            vec![
                ("population_id", population_id.to_owned()),
                ("block_index_u64", block_index.to_string()),
                ("haplotype_rank_u64", rank.to_string()),
                ("alleles", alleles.to_owned()),
            ],
        )
    }

    fn founder_copy(
        &self,
        population_id: &str,
        individual_id: &str,
        sequence_id: &str,
        copy_index: u64,
    ) -> Checked<String> {
        self.derive(
            "FCP_",
            "f5_founder_copy",
            vec![
                ("population_id", population_id.to_owned()),
                ("individual_id", individual_id.to_owned()),
                ("sequence_id", sequence_id.to_owned()),
                ("copy_index_u32", copy_index.to_string()),
            ],
        )
    }

    fn assignment(
        &self,
        population_id: &str,
        block_index: u64,
        founder_copy_id: &str,
        haplotype_id: &str,
    ) -> Checked<String> {
        self.derive(
            "FHA_",
            "f5_block_haplotype_assignment",
            vec![
                ("population_id", population_id.to_owned()),
                ("block_index_u64", block_index.to_string()),
                ("founder_copy_id", founder_copy_id.to_owned()),
                ("haplotype_id", haplotype_id.to_owned()),
            ],
        )
    }
}

/// Pares `name=value` de una columna de ejemplo del catálogo, en su orden.
fn catalog_components(example: &str) -> Checked<Vec<(&str, String)>> {
    example
        .split(';')
        .map(|pair| {
            pair.split_once('=')
                .map(|(name, value)| (name, value.to_owned()))
                .ok_or_else(|| catalog_reject(format!("componente sin `=`: {pair:?}")))
        })
        .collect()
}

fn catalog_material(root: &[u8; 32], components: &[(&str, String)]) -> Checked<Vec<u8>> {
    let owned: Vec<(&str, String)> = components.to_vec();
    canonical_material(root, &owned)
}

/// Reproduce los cinco vectores F5 congelados de `id-catalog-v4.tsv`.
///
/// Sus seis rechazos se ejercitan con catálogos dañados a través de
/// `validate_id_catalog_text`.
fn validate_id_catalog() -> Checked<()> {
    validate_id_catalog_text(ID_CATALOG_TSV)
}

/// El cuerpo de [`validate_id_catalog`], sobre el texto que se le dé: es lo que
/// permite probar cada rechazo con un catálogo dañado, como ya se hace con
/// `validate_frozen_rng_catalog`.
fn validate_id_catalog_text(catalog: &str) -> Checked<()> {
    let mut seen = 0_usize;
    for line in catalog.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('|').collect();
        if fields.len() != 9 {
            return Err(catalog_reject(
                "una fila del catálogo de IDs no tiene 9 campos",
            ));
        }
        if !fields[0].starts_with("f5_") {
            continue;
        }
        seen += 1;
        let root: [u8; 32] = decode_hex(fields[5])
            .and_then(|bytes| bytes.try_into().ok())
            .ok_or_else(|| catalog_reject("raíz del catálogo de IDs ilegible"))?;
        let components = catalog_components(fields[4])?;
        let material = catalog_material(&root, &components)?;
        if hex(&material) != fields[6] {
            return Err(catalog_reject(format!(
                "el material rederivado de {} no es el congelado",
                fields[0]
            )));
        }
        let digest = blake3::derive_key(ID_CONTEXT, &material);
        if hex(&digest) != fields[7] {
            return Err(catalog_reject(format!(
                "el digest rederivado de {} no es el congelado",
                fields[0]
            )));
        }
        let identifier = format!("{}{}", fields[1], base32(&digest[..16]));
        if identifier != fields[8] {
            return Err(catalog_reject(format!(
                "el identificador rederivado de {} no es el congelado",
                fields[0]
            )));
        }
    }
    if seen != 5 {
        return Err(catalog_reject(format!(
            "el catálogo de IDs declara {seen} clases F5 y el contrato fija cinco"
        )));
    }
    Ok(())
}

fn hex(bytes: &[u8]) -> String {
    let mut text = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        text.push_str(&format!("{byte:02x}"));
    }
    text
}

// ---------------------------------------------------------------------------
// Oráculo de subsemillas DTC-089
// ---------------------------------------------------------------------------

/// Rederivador de la única decisión RNG que los bytes publican: la polaridad.
struct SubseedOracle {
    root: [u8; 32],
    scenario_id: String,
}

impl SubseedOracle {
    fn new(root: [u8; 32], scenario_id: &str) -> Checked<Self> {
        validate_rng_catalog()?;
        Ok(Self {
            root,
            scenario_id: scenario_id.to_owned(),
        })
    }

    /// Marginal objetivo de un sitio×subpoblación: la mecánica cerrada de
    /// DTC-089 —fila y valor, en ese orden, cada uno con su propio U53 y solo
    /// cuando tiene más de una alternativa— sobre la subsemilla
    /// `marginal_spectrum`, cuyos cuatro componentes publica
    /// `f5-frequencies.tsv` columna a columna.
    fn target_alt_af(
        &self,
        population_id: &str,
        block_id: &str,
        site_id: &str,
        bands: &[F5DeclaredSpectrumBand],
    ) -> Checked<f64> {
        let components = vec![
            ("scenario_id", self.scenario_id.clone()),
            ("population_id", population_id.to_owned()),
            ("block_id", block_id.to_owned()),
            ("site_id", site_id.to_owned()),
            ("stage", "marginal_spectrum".to_owned()),
        ];
        let material = canonical_material(&self.root, &components)?;
        let subseed = blake3::derive_key(SUBSEED_CONTEXT, &material);
        let mut rng = ChaCha12Rng::from_seed(subseed);
        let (first, rest) = bands
            .split_first()
            .ok_or_else(|| catalog_reject("el espectro declarado no trae ninguna banda"))?;
        // Una elección con una sola alternativa no consume RNG.
        let band = if rest.is_empty() {
            *first
        } else {
            // Se recorre el peso en orden de fila y se elige la primera suma
            // acumulada estrictamente mayor que U; U **no** se reescala por el
            // total, porque los pesos ya son los efectivos de DTC-049.
            let threshold = uniform_53(&mut rng);
            let mut prefix = Neumaier::default();
            let mut chosen = None;
            for candidate in bands {
                prefix.add(candidate.weight);
                if threshold < prefix.total() {
                    chosen = Some(*candidate);
                    break;
                }
            }
            // La última fila absorbe en exclusiva el residuo de redondeo.
            chosen.unwrap_or(*bands.last().unwrap_or(first))
        };
        if band.alt_af_min == band.alt_af_max {
            return Ok(band.alt_af_min);
        }
        let unit = uniform_53(&mut rng);
        Ok(band.alt_af_min + (band.alt_af_max - band.alt_af_min) * unit)
    }

    /// Alelo alternativo de un sitio: `f5_uniform_index(rng, 3)` sobre la
    /// subsemilla `site_alt` elige entre las tres bases distintas de la
    /// referencia, **en el orden `ACGT`**.
    ///
    /// Los tres componentes de la decisión están publicados: `block_id` y
    /// `site_id` en `f5-sites.tsv` y `scenario_id` en la declaración. Sin esta
    /// rederivación el ALT solo se comprobaba contra `ACGT` y contra REF, y
    /// una sustitución coherente en todas sus vistas —columna, modificaciones,
    /// fragmento vigente, `f5_locus` y alelos tipados— pasaba inadvertida,
    /// publicando una verdad fundadora cuyos individuos portan una base que el
    /// contrato no sortea.
    fn alternate_base(&self, block_id: &str, site_id: &str, reference_base: char) -> Checked<char> {
        let components = vec![
            ("scenario_id", self.scenario_id.clone()),
            ("block_id", block_id.to_owned()),
            ("site_id", site_id.to_owned()),
            ("stage", "site_alt".to_owned()),
        ];
        let material = canonical_material(&self.root, &components)?;
        let subseed = blake3::derive_key(SUBSEED_CONTEXT, &material);
        let mut rng = ChaCha12Rng::from_seed(subseed);
        let alternatives: Vec<char> = ['A', 'C', 'G', 'T']
            .into_iter()
            .filter(|base| *base != reference_base)
            .collect();
        let choice = uniform_index(&mut rng, 3);
        // Los dos rechazos siguientes son inalcanzables: `uniform_index(rng, 3)`
        // devuelve un valor menor que 3, que cabe en cualquier `usize`, y
        // `alternatives` conserva al menos tres de las cuatro bases (solo se quita
        // la de referencia), de modo que `get` nunca da `None`.
        let position = usize::try_from(choice)
            .map_err(|_| catalog_reject("el índice ALT no cabe en usize"))?;
        alternatives
            .get(position)
            .copied()
            .ok_or_else(|| catalog_reject("el selector ALT salió de las tres alternativas"))
    }

    /// Polaridad de un sitio: `positive` sin consumir entropía en el sitio 0 y
    /// Bernoulli(0,5) sobre su propio ChaCha12 en cualquier otro.
    fn polarity(
        &self,
        population_id: &str,
        block_id: &str,
        attempt_index: u64,
        site_index: u64,
    ) -> Checked<&'static str> {
        if site_index == 0 {
            return Ok("positive");
        }
        let components = vec![
            ("scenario_id", self.scenario_id.clone()),
            ("population_id", population_id.to_owned()),
            ("block_id", block_id.to_owned()),
            ("stage", "polarity".to_owned()),
            ("attempt_index_u64", attempt_index.to_string()),
            ("salt.site_index_u64", site_index.to_string()),
        ];
        let material = canonical_material(&self.root, &components)?;
        let subseed = blake3::derive_key(SUBSEED_CONTEXT, &material);
        let mut rng = ChaCha12Rng::from_seed(subseed);
        let threshold = (0.5_f64 * FIFTY_THREE_BIT_BOUND as f64).floor() as u64;
        let bits = rng.next_u64() >> 11;
        Ok(if bits < threshold {
            "positive"
        } else {
            "negative"
        })
    }

    /// Recuento haplotípico de un bloque×subpoblación, rederivado.
    ///
    /// Reimplementa `poisson-count-v1` —inversión por producto de uniformes
    /// bajo `lambda < 10` y el rechazo transformado de Hörmann por encima— y
    /// pinza el resultado en `[h_min, h_max]`. Los tres componentes de la
    /// decisión están publicados: `population_id` y `block_id` los trae cada
    /// fila de `f5-haplotypes.tsv` y `scenario_id` la declaración.
    ///
    /// **Lo que esta rederivación ata son los recuentos, no la lambda.** El
    /// oráculo no abre el perfil: recibe `lambda`, `h_min` y `h_max` de la
    /// declaración y comprueba que los recuentos publicados salgan de ellos.
    /// Como el pinzamiento destruye información —toda la masa fuera del
    /// intervalo colapsa en sus extremos—, muchas lambdas distintas producen
    /// los mismos recuentos y el oráculo las acepta todas. Cerrar esa frontera
    /// exigiría releer el metadata del perfil, que es la misma concesión ya
    /// declarada para la auditoría DTC-049, el espectro, los bins LD y los
    /// pesos.
    fn block_haplotype_count(
        &self,
        population_id: &str,
        block_id: &str,
        distribution: F5DeclaredHaplotypeCountDistribution,
    ) -> Checked<u64> {
        let components = vec![
            ("scenario_id", self.scenario_id.clone()),
            ("population_id", population_id.to_owned()),
            ("block_id", block_id.to_owned()),
            ("stage", "block_haplotype_count".to_owned()),
        ];
        let material = canonical_material(&self.root, &components)?;
        let subseed = blake3::derive_key(SUBSEED_CONTEXT, &material);
        let mut rng = ChaCha12Rng::from_seed(subseed);
        let drawn = poisson_count(&mut rng, distribution.lambda).ok_or_else(|| {
            catalog_reject("la lambda declarada no admite el muestreador de recuentos")
        })?;
        Ok(drawn.clamp(distribution.minimum, distribution.maximum))
    }
}

/// Uniforme abierto de 52 bits que consume `poisson-count-v1`.
fn open_unit(rng: &mut ChaCha12Rng) -> f64 {
    (((rng.next_u64() >> 12) as f64) + 0.5) * (1.0 / 4_503_599_627_370_496.0)
}

/// Reimplementación de `poisson-count-v1` sin llamar al productor.
///
/// Devuelve `None` cuando la lambda no pertenece al dominio del muestreador,
/// cuando un candidato no es numérico o pasa de 2^53, o cuando el contador de
/// la rama de medias pequeñas desborda. Un candidato negativo no da `None`: se
/// repite.
fn poisson_count(rng: &mut ChaCha12Rng, lambda: f64) -> Option<u64> {
    const MAX_LAMBDA: f64 = 4_503_599_627_370_496.0;
    if !lambda.is_finite() || !(0.0..=MAX_LAMBDA).contains(&lambda) {
        return None;
    }
    if lambda == 0.0 {
        return Some(0);
    }
    if lambda < 10.0 {
        let threshold = libm::exp(-lambda);
        let mut product = 1.0_f64;
        let mut count = 0_u64;
        loop {
            product *= open_unit(rng);
            if product <= threshold {
                return Some(count);
            }
            count = count.checked_add(1)?;
        }
    }
    let sqrt_lambda = libm::sqrt(lambda);
    let b = 0.931 + 2.53 * sqrt_lambda;
    let a = -0.059 + 0.02483 * b;
    let inverse_alpha = 1.1239 + 1.1328 / (b - 3.4);
    let vr = 0.9277 - 3.6224 / (b - 2.0);
    let log_lambda = libm::log(lambda);
    loop {
        let u = open_unit(rng) - 0.5;
        let v = open_unit(rng);
        let us = 0.5 - u.abs();
        let candidate = libm::floor((2.0 * a / us + b) * u + lambda + 0.43);
        // Misma regla que el productor (DT-010): un candidato no numérico o
        // por encima de 2^53 no tiene valor que repetir y se rechaza; uno
        // **negativo** se descarta y se vuelve a tirar.
        if candidate.is_nan() || candidate > 9_007_199_254_740_991.0 {
            return None;
        }
        if candidate < 0.0 {
            continue;
        }
        if us >= 0.07 && v <= vr {
            return Some(candidate as u64);
        }
        if us < 0.013 && v > us {
            continue;
        }
        let left = libm::log(v * inverse_alpha / (a / (us * us) + b));
        let right = -lambda + candidate * log_lambda - libm::lgamma(candidate + 1.0);
        if left <= right {
            return Some(candidate as u64);
        }
    }
}

/// Entero uniforme en `0..n` con el rechazo por muestreo de DTC-089.
///
/// Reimplementa `f5_uniform_index` sin llamarlo: el oráculo no importa nada del
/// productor, tampoco su aritmética.
fn uniform_index(rng: &mut ChaCha12Rng, n: u64) -> u64 {
    if n <= 1 {
        return 0;
    }
    let remainder = (u64::MAX % n + 1) % n;
    let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
    loop {
        let value = rng.next_u64();
        if limit.is_none_or(|bound| value < bound) {
            return value % n;
        }
    }
}

/// U53 normativo de DTC-089: `bits53 / 2^53`, cerrado en 0 y abierto en 1.
fn uniform_53(rng: &mut ChaCha12Rng) -> f64 {
    (rng.next_u64() >> 11) as f64 / FIFTY_THREE_BIT_BOUND as f64
}

/// Reproduce las subsemillas y las cuatro primeras palabras de cada decisión
/// de los **dos** catálogos congelados.
fn validate_rng_catalog() -> Checked<()> {
    validate_frozen_rng_catalog(RNG_CATALOG_TSV, F5_RNG_CATALOG_DECISIONS)?;
    validate_frozen_rng_catalog(RNG_CATALOG_V2_TSV, F5_RNG_CATALOG_V2_DECISIONS)
}

fn validate_frozen_rng_catalog(catalog: &str, expected_decisions: usize) -> Checked<()> {
    let mut seen = 0_usize;
    for line in catalog.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('|').collect();
        if fields.len() != 9 {
            return Err(catalog_reject(
                "una fila del catálogo RNG no tiene 9 campos",
            ));
        }
        seen += 1;
        let root: [u8; 32] = decode_hex(fields[2])
            .and_then(|bytes| bytes.try_into().ok())
            .ok_or_else(|| catalog_reject("raíz del catálogo RNG ilegible"))?;
        let components = catalog_components(fields[1])?;
        let material = catalog_material(&root, &components)?;
        if hex(&material) != fields[3] {
            return Err(catalog_reject(format!(
                "el material rederivado de {} no es el congelado",
                fields[0]
            )));
        }
        let subseed = blake3::derive_key(SUBSEED_CONTEXT, &material);
        if hex(&subseed) != fields[4] {
            return Err(catalog_reject(format!(
                "la subsemilla rederivada de {} no es la congelada",
                fields[0]
            )));
        }
        let mut rng = ChaCha12Rng::from_seed(subseed);
        for offset in 0..4 {
            let expected = u64::from_str_radix(fields[5 + offset], 16)
                .map_err(|_| catalog_reject("palabra del catálogo RNG ilegible"))?;
            if rng.next_u64() != expected {
                return Err(catalog_reject(format!(
                    "la palabra {offset} de {} no es la congelada",
                    fields[0]
                )));
            }
        }
    }
    // El reparto de este módulo —las rederivadas y las creídas— se sostiene
    // sobre un número exacto de filas. Exigir solo que no esté vacío dejaba al
    // oráculo autovalidándose contra una sola decisión si el catálogo se
    // recortaba.
    if seen != expected_decisions {
        return Err(catalog_reject(format!(
            "el catálogo RNG trae {seen} decisiones y el contrato fija {expected_decisions}"
        )));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Lector estricto de los siete TSV
// ---------------------------------------------------------------------------

/// Cursor de filas de un TSV publicado que aplica el contrato de bytes sin
/// materializar el fichero ni asignar por celda.
struct TsvCursor<R: BufRead> {
    reader: R,
    path: &'static str,
    columns: usize,
    line: String,
    fields: Vec<(usize, usize)>,
    rows: u64,
    closed_with_lf: bool,
}

impl<R: BufRead> TsvCursor<R> {
    fn open(mut reader: R, index: usize) -> Checked<Self> {
        let contract = &F5_TSV_CONTRACTS[index];
        let path = contract.path();
        let mut header = String::new();
        let read = reader.read_line(&mut header).map_err(|error| {
            reject(
                F5OracleCause::ByteContract,
                path,
                "encoding",
                error.to_string(),
            )
        })?;
        if read == 0 {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                "header",
                "el artefacto está vacío",
            ));
        }
        if header.starts_with('\u{feff}') {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                "bom",
                "BOM prohibido",
            ));
        }
        let closed_with_lf = header.ends_with('\n');
        let body = header.strip_suffix('\n').unwrap_or(&header);
        if body.contains('\r') {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                "eol",
                "CR prohibido",
            ));
        }
        if body != contract.header() {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                "header",
                format!("cabecera {body:?} en lugar de la del contrato"),
            ));
        }
        Ok(Self {
            reader,
            path,
            columns: contract.header().split('\t').count(),
            line: String::new(),
            fields: Vec::new(),
            rows: 0,
            closed_with_lf,
        })
    }

    /// Avanza a la siguiente fila. Devuelve `false` al terminar el fichero, y
    /// solo entonces exige el LF final.
    fn advance(&mut self) -> Checked<bool> {
        self.line.clear();
        let read = self.reader.read_line(&mut self.line).map_err(|error| {
            reject(
                F5OracleCause::ByteContract,
                self.path,
                "encoding",
                error.to_string(),
            )
        })?;
        if read == 0 {
            if !self.closed_with_lf {
                return Err(reject(
                    F5OracleCause::ByteContract,
                    self.path,
                    "eol",
                    "falta el LF final",
                ));
            }
            return Ok(false);
        }
        self.closed_with_lf = self.line.ends_with('\n');
        let body_length = self.line.len() - usize::from(self.closed_with_lf);
        let body = &self.line[..body_length];
        if body.contains('\r') {
            return Err(reject(
                F5OracleCause::ByteContract,
                self.path,
                "eol",
                "CR prohibido",
            ));
        }
        if body.contains('"') {
            return Err(reject(
                F5OracleCause::ByteContract,
                self.path,
                "cell",
                "comillas prohibidas",
            ));
        }
        if body.is_empty() {
            return Err(reject(
                F5OracleCause::ByteContract,
                self.path,
                "row",
                "línea vacía prohibida",
            ));
        }
        if body.starts_with('#') {
            return Err(reject(
                F5OracleCause::ByteContract,
                self.path,
                "row",
                "comentario prohibido",
            ));
        }
        self.fields.clear();
        let mut start = 0_usize;
        for (offset, byte) in body.bytes().enumerate() {
            if byte == b'\t' {
                self.fields.push((start, offset));
                start = offset + 1;
            }
        }
        self.fields.push((start, body_length));
        if self.fields.len() != self.columns {
            return Err(reject(
                F5OracleCause::ByteContract,
                self.path,
                "row",
                format!(
                    "la fila tiene {} columnas y la cabecera {}",
                    self.fields.len(),
                    self.columns
                ),
            ));
        }
        self.rows += 1;
        Ok(true)
    }

    fn field(&self, index: usize) -> &str {
        let (start, end) = self.fields[index];
        &self.line[start..end]
    }

    fn u64(&self, cause: F5OracleCause, index: usize) -> Checked<u64> {
        canonical_u64(cause, self.path, self.field(index))
    }

    fn f64(&self, cause: F5OracleCause, index: usize) -> Checked<f64> {
        canonical_f64(cause, self.path, self.field(index))
    }

    fn probability(&self, cause: F5OracleCause, index: usize) -> Checked<f64> {
        let value = self.f64(cause, index)?;
        if !(0.0..=1.0).contains(&value) {
            return Err(reject(
                cause,
                self.path,
                "probability",
                format!("{value} no pertenece a [0,1]"),
            ));
        }
        Ok(value)
    }

    fn base(&self, cause: F5OracleCause, index: usize) -> Checked<char> {
        let value = self.field(index);
        let mut characters = value.chars();
        match (characters.next(), characters.next()) {
            (Some(base @ ('A' | 'C' | 'G' | 'T')), None) => Ok(base),
            _ => Err(reject(
                cause,
                self.path,
                "base",
                format!("{value:?} no es una base canónica"),
            )),
        }
    }
}

fn canonical_u64(cause: F5OracleCause, path: &str, value: &str) -> Checked<u64> {
    if value.is_empty()
        || (value.len() > 1 && value.starts_with('0'))
        || !value.bytes().all(|byte| byte.is_ascii_digit())
    {
        return Err(reject(
            cause,
            path,
            "integer",
            format!("{value:?} no es un entero decimal canónico"),
        ));
    }
    value
        .parse()
        .map_err(|_| reject(cause, path, "integer", format!("{value:?} no cabe en u64")))
}

/// Reexige sobre los bytes publicados el literal decimal canónico de F5. Como
/// el resto de contratos literales que este módulo comparte —rutas, clases y
/// versiones—, la forma canónica es contrato y no estado del productor.
fn canonical_f64(cause: F5OracleCause, path: &str, value: &str) -> Checked<f64> {
    let parsed: f64 = value
        .parse()
        .map_err(|_| reject(cause, path, "float", format!("{value:?} no es un f64")))?;
    if !parsed.is_finite() || f5_canonical_decimal(parsed) != value {
        return Err(reject(
            cause,
            path,
            "float",
            format!("{value:?} no es el decimal canónico de anchura acotada"),
        ));
    }
    Ok(parsed)
}

fn open_artifact(root: &Path, path: &str) -> Checked<BufReader<File>> {
    let full = root.join(path);
    File::open(&full).map(BufReader::new).map_err(|error| {
        reject(
            F5OracleCause::ArtifactMissing,
            path,
            "path",
            format!("{}: {error}", full.display()),
        )
    })
}

/// Tamaño y SHA-256 de un artefacto, leídos por flujo.
fn measure_artifact(root: &Path, path: &str) -> Checked<(u64, String)> {
    let mut reader = open_artifact(root, path)?;
    let mut hasher = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    let mut size = 0_u64;
    loop {
        let read = reader.read(&mut buffer).map_err(|error| {
            reject(
                F5OracleCause::ArtifactMissing,
                path,
                "read",
                error.to_string(),
            )
        })?;
        if read == 0 {
            break;
        }
        hasher.update(&buffer[..read]);
        size += read as u64;
    }
    Ok((size, format!("{:x}", hasher.finalize())))
}

// ---------------------------------------------------------------------------
// Lector estricto de `truth.jsonl` v2
// ---------------------------------------------------------------------------

/// Envoltorio puro sobre [`reject`]. Lleva `#[track_caller]` para que el
/// barrido anote el punto de rechazo real y no esta línea.
#[track_caller]
fn truth_reject(field: &str, detail: impl Into<String>) -> F5OracleRejection {
    reject(F5OracleCause::Truth, F5_TRUTH_PATH, field, detail)
}

fn exact_object<'a>(
    value: &'a Value,
    what: &str,
    keys: &[&str],
) -> Checked<&'a Map<String, Value>> {
    let object = value
        .as_object()
        .ok_or_else(|| truth_reject(what, "no es un objeto"))?;
    if object.len() != keys.len() || !keys.iter().all(|key| object.contains_key(*key)) {
        let observed: Vec<&str> = object.keys().map(String::as_str).collect();
        return Err(truth_reject(
            what,
            format!("claves {observed:?} en lugar de {keys:?}"),
        ));
    }
    Ok(object)
}

fn text(object: &Map<String, Value>, what: &str, key: &str) -> Checked<String> {
    object
        .get(key)
        .and_then(Value::as_str)
        .map(str::to_owned)
        .ok_or_else(|| truth_reject(what, format!("{key} no es una cadena")))
}

fn number(object: &Map<String, Value>, what: &str, key: &str) -> Checked<u64> {
    object
        .get(key)
        .and_then(Value::as_u64)
        .ok_or_else(|| truth_reject(what, format!("{key} no es un entero no negativo")))
}

fn single_base(object: &Map<String, Value>, what: &str, key: &str) -> Checked<char> {
    let value = text(object, what, key)?;
    let mut characters = value.chars();
    match (characters.next(), characters.next()) {
        (Some(base @ ('A' | 'C' | 'G' | 'T')), None) => Ok(base),
        _ => Err(truth_reject(what, format!("{key} no es una base canónica"))),
    }
}

/// Subpoblación declarada por la cabecera v2.
#[derive(Clone, Debug, Eq, PartialEq)]
struct TruthPopulation {
    population_id: String,
    individual_count: u64,
    ploidy: u64,
}

/// Cabecera v2 ya validada.
#[derive(Clone, Debug, Eq, PartialEq)]
struct TruthHeader {
    execution_id: String,
    reference_id: String,
    scenario_id: String,
    reference_content_sha256: String,
    /// `(kind, source_id, sha256)` en el orden publicado.
    source_digests: Vec<(String, String, String)>,
    populations: Vec<TruthPopulation>,
    copies: u64,
    genotypes: u64,
}

impl TruthHeader {
    fn population(&self, population_id: &str) -> Option<&TruthPopulation> {
        self.populations
            .iter()
            .find(|entry| entry.population_id == population_id)
    }

    fn matrix_source_id(&self) -> Checked<&str> {
        self.source_digests
            .iter()
            .find(|(kind, _, _)| kind == "matrix")
            .map(|(_, source_id, _)| source_id.as_str())
            .ok_or_else(|| truth_reject("source_digests", "ninguna fuente es la matriz"))
    }
}

fn read_truth_header(line: &str) -> Checked<TruthHeader> {
    let value: Value = serde_json::from_str(line)
        .map_err(|error| truth_reject("header", format!("cabecera ilegible: {error}")))?;
    let object = exact_object(
        &value,
        "header",
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
    if text(object, "header", "record_type")? != "header" {
        return Err(truth_reject("header", "record_type no es header"));
    }
    if text(object, "header", "document_kind")? != "gps_truth" {
        return Err(truth_reject("header", "document_kind no es gps_truth"));
    }
    let version = text(object, "header", "schema_version")?;
    if version != F5_TRUTH_SCHEMA_VERSION {
        return Err(truth_reject(
            "header",
            format!("schema_version {version:?}; una proyección F5 exige la 2"),
        ));
    }
    if text(object, "header", "coordinate_convention")? != "zero_based_half_open_v1" {
        return Err(truth_reject("header", "convenio de coordenadas ajeno"));
    }
    let reference_content_sha256 = text(object, "header", "reference_content_sha256")?;
    if !is_lowercase_sha256(&reference_content_sha256) {
        return Err(truth_reject(
            "header",
            "reference_content_sha256 no es hex minúscula de 64",
        ));
    }
    let projection = exact_object(
        object
            .get("projection")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("header", "sin projection"))?,
        "projection",
        &["kind", "projection_version", "scenario_id"],
    )?;
    if text(projection, "projection", "kind")? != F5_TRUTH_PROJECTION_KIND {
        return Err(truth_reject("projection", "kind no es f5_founder"));
    }
    if text(projection, "projection", "projection_version")? != F5_TRUTH_PROJECTION_VERSION {
        return Err(truth_reject("projection", "projection_version ajena"));
    }
    let counts = exact_object(
        object
            .get("counts")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("header", "sin counts"))?,
        "counts",
        &["copies", "genotypes"],
    )?;
    let listed = object
        .get("populations")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("header", "populations no es una lista"))?;
    if listed.is_empty() {
        return Err(truth_reject("header", "populations vacía"));
    }
    let mut populations = Vec::with_capacity(listed.len());
    for entry in listed {
        let population = exact_object(
            entry,
            "populations",
            &["population_id", "individual_count", "ploidy"],
        )?;
        populations.push(TruthPopulation {
            population_id: text(population, "populations", "population_id")?,
            individual_count: number(population, "populations", "individual_count")?,
            ploidy: number(population, "populations", "ploidy")?,
        });
    }
    if populations
        .windows(2)
        .any(|pair| pair[0].population_id >= pair[1].population_id)
    {
        return Err(truth_reject(
            "header",
            "populations no está ordenada por population_id",
        ));
    }
    let listed_digests = object
        .get("source_digests")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("header", "source_digests no es una lista"))?;
    if listed_digests.is_empty() {
        return Err(truth_reject("header", "source_digests vacía"));
    }
    let mut source_digests = Vec::with_capacity(listed_digests.len());
    for entry in listed_digests {
        let entry = exact_object(entry, "source_digests", &["source", "sha256"])?;
        let source = exact_object(
            entry
                .get("source")
                // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
                // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
                .ok_or_else(|| truth_reject("source_digests", "sin source"))?,
            "source_digests.source",
            &["kind", "source_id"],
        )?;
        let digest = text(entry, "source_digests", "sha256")?;
        if !is_lowercase_sha256(&digest) {
            return Err(truth_reject(
                "source_digests",
                "sha256 no es hex minúscula de 64",
            ));
        }
        source_digests.push((
            text(source, "source_digests.source", "kind")?,
            text(source, "source_digests.source", "source_id")?,
            digest,
        ));
    }
    Ok(TruthHeader {
        execution_id: text(object, "header", "execution_id")?,
        reference_id: text(object, "header", "reference_id")?,
        scenario_id: text(projection, "projection", "scenario_id")?,
        reference_content_sha256,
        source_digests,
        populations,
        copies: number(counts, "counts", "copies")?,
        genotypes: number(counts, "counts", "genotypes")?,
    })
}

/// Identidad de copia publicada.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CopyKey {
    reference_id: String,
    individual_id: String,
    sequence_id: String,
    copy_index: u64,
}

fn read_copy_identity(value: &Value, what: &str) -> Checked<CopyKey> {
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

/// Un SNP de copia ya validado en sus tres vistas.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct CopySnp {
    site_id: String,
    position_1based: u64,
    reference_base: char,
    alternate_base: char,
}

/// Una `haplotype_copy` completa, ya comprobada contra sí misma.
#[derive(Clone, Debug, PartialEq)]
struct TruthCopy {
    identity: CopyKey,
    population_id: String,
    founder_copy_id: String,
    sequence_length: u64,
    /// Ordenados por posición ascendente, que es como la copia los publica.
    snps: Vec<CopySnp>,
    /// `matrix_source_id` de cada SNP, en el mismo orden que `snps`.
    haplotype_ids: Vec<String>,
}

fn read_snp_provenance(
    object: &Map<String, Value>,
    event_id: &str,
    population_id: &str,
    matrix_source_id: &str,
) -> Checked<String> {
    let provenance = exact_object(
        object
            .get("provenance")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.snp", "SNP sin procedencia"))?,
        "copy.snp.provenance",
        &["kind", "version", "source", "matrix"],
    )?;
    if text(provenance, "copy.snp.provenance", "kind")? != "matrix" {
        return Err(truth_reject(
            "copy.snp.provenance",
            "un SNP fundador no procede de una matriz",
        ));
    }
    if text(provenance, "copy.snp.provenance", "version")? != "1" {
        return Err(truth_reject(
            "copy.snp.provenance",
            "versión de procedencia ajena",
        ));
    }
    let source = exact_object(
        provenance
            .get("source")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.snp.provenance", "sin source"))?,
        "copy.snp.provenance.source",
        &["kind", "source_id"],
    )?;
    if text(source, "copy.snp.provenance.source", "kind")? != "matrix"
        || text(source, "copy.snp.provenance.source", "source_id")? != matrix_source_id
    {
        return Err(truth_reject(
            "copy.snp.provenance.source",
            "la fuente del SNP no es la matriz que la cabecera declara",
        ));
    }
    let matrix = exact_object(
        provenance
            .get("matrix")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.snp.provenance", "sin matrix"))?,
        "copy.snp.provenance.matrix",
        &["matrix_source_id", "population_id", "variant_id"],
    )?;
    if text(matrix, "copy.snp.provenance.matrix", "variant_id")? != event_id {
        return Err(truth_reject(
            "copy.snp.provenance.matrix",
            "la procedencia no nombra el sitio de su evento",
        ));
    }
    let declared = text(matrix, "copy.snp.provenance.matrix", "population_id")?;
    if declared != population_id {
        return Err(truth_reject(
            "copy.snp.provenance.matrix",
            format!(
                "la procedencia del SNP {event_id} declara {declared} y su copia pertenece a {population_id}"
            ),
        ));
    }
    text(matrix, "copy.snp.provenance.matrix", "matrix_source_id")
}

fn read_segment_source(value: &Value, identity: &CopyKey) -> Checked<u64> {
    let source = exact_object(
        value,
        "copy.segment.source",
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
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.segment.source", "sin source"))?,
        "copy.segment.source.source",
        &["kind", "source_id"],
    )?;
    if text(source, "copy.segment.source", "reference_id")? != identity.reference_id
        || text(source, "copy.segment.source", "sequence_id")? != identity.sequence_id
    {
        return Err(truth_reject(
            "copy.segment.source",
            "el segmento no pertenece a la copia",
        ));
    }
    let start = number(source, "copy.segment.source", "start")?;
    let end = number(source, "copy.segment.source", "end")?;
    let length = number(source, "copy.segment.source", "sequence_length")?;
    if start != 0 || end != length || length == 0 {
        return Err(truth_reject(
            "copy.segment.source",
            "el segmento no cubre la secuencia completa",
        ));
    }
    Ok(length)
}

#[allow(clippy::too_many_lines)]
fn read_truth_copy(value: &Value, header: &TruthHeader) -> Checked<TruthCopy> {
    let object = exact_object(
        value,
        "copy",
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
        object
            .get("identity")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy", "sin identity"))?,
        "copy.identity",
    )?;
    let phase = number(object, "copy", "phase")?;
    if phase != identity.copy_index {
        return Err(truth_reject("copy", "la fase no es el índice de copia"));
    }
    let membership = exact_object(
        object
            .get("founder_membership")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy", "sin founder_membership"))?,
        "copy.founder_membership",
        &["population_id", "founder_copy_id"],
    )?;
    let population_id = text(membership, "copy.founder_membership", "population_id")?;
    if header.population(&population_id).is_none() {
        return Err(truth_reject(
            "copy.founder_membership",
            format!("la pertenencia nombra {population_id}, ausente de la cabecera"),
        ));
    }
    let founder_copy_id = text(membership, "copy.founder_membership", "founder_copy_id")?;

    let provenance = exact_object(
        object
            .get("provenance")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy", "sin provenance"))?,
        "copy.provenance",
        &["origin", "sources"],
    )?;
    if text(provenance, "copy.provenance", "origin")? != "f5_founder" {
        return Err(truth_reject("copy.provenance", "origin ajeno a F5"));
    }
    let sources = provenance
        .get("sources")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy.provenance", "sources no es una lista"))?;
    let mut declared_sources = Vec::with_capacity(sources.len());
    for source in sources {
        let source = exact_object(source, "copy.provenance.source", &["kind", "source_id"])?;
        declared_sources.push((
            text(source, "copy.provenance.source", "kind")?,
            text(source, "copy.provenance.source", "source_id")?,
        ));
    }
    let expected_sources: Vec<(String, String)> = header
        .source_digests
        .iter()
        .map(|(kind, source_id, _)| (kind.clone(), source_id.clone()))
        .collect();
    if declared_sources != expected_sources {
        return Err(truth_reject(
            "copy.provenance",
            format!(
                "las fuentes {declared_sources:?} no son las {expected_sources:?} de la cabecera"
            ),
        ));
    }
    let matrix_source_id = header.matrix_source_id()?.to_owned();

    let segments = object
        .get("segments")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy", "segments no es una lista"))?;
    if segments.len() != 1 {
        return Err(truth_reject(
            "copy.segments",
            "una copia fundadora publica un segmento por secuencia",
        ));
    }
    let segment = exact_object(&segments[0], "copy.segment", &["source", "modifications"])?;
    let sequence_length = read_segment_source(
        segment
            .get("source")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.segment", "sin source"))?,
        &identity,
    )?;

    let mut snps: Vec<CopySnp> = Vec::new();
    let mut haplotype_ids: Vec<String> = Vec::new();
    for value in segment
        .get("modifications")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy.segment", "modifications no es una lista"))?
    {
        let modification = exact_object(
            value,
            "copy.modification",
            &[
                "kind",
                "event_id",
                "position",
                "reference",
                "alternate",
                "provenance",
            ],
        )?;
        if text(modification, "copy.modification", "kind")? != "snp" {
            return Err(truth_reject(
                "copy.modification",
                "F5 solo publica modificaciones SNP",
            ));
        }
        let site_id = text(modification, "copy.modification", "event_id")?;
        let position = number(modification, "copy.modification", "position")?;
        if position >= sequence_length {
            return Err(truth_reject(
                "copy.modification",
                "un SNP cae fuera de su secuencia",
            ));
        }
        let reference_base = single_base(modification, "copy.modification", "reference")?;
        let alternate_base = single_base(modification, "copy.modification", "alternate")?;
        if reference_base == alternate_base {
            return Err(truth_reject(
                "copy.modification",
                "un SNP publica el mismo alelo dos veces",
            ));
        }
        // Orden **estrictamente** creciente: `last.position_1based` es
        // `anterior_0based + 1` y `position` es la coordenada 0-based de esta
        // modificación, de modo que la comparación es `anterior >= actual` y
        // dos modificaciones en la misma posición se rechazan aquí, no aguas
        // abajo.
        if snps
            .last()
            .is_some_and(|last| last.position_1based > position)
        {
            return Err(truth_reject(
                "copy.modification",
                "las modificaciones no llegan por posición ascendente",
            ));
        }
        haplotype_ids.push(read_snp_provenance(
            modification,
            &site_id,
            &population_id,
            &matrix_source_id,
        )?);
        snps.push(CopySnp {
            site_id,
            position_1based: position + 1,
            reference_base,
            alternate_base,
        });
    }

    let pieces = object
        .get("current_pieces")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy", "current_pieces no es una lista"))?;
    if pieces.len() != 1 {
        return Err(truth_reject(
            "copy.current_pieces",
            "una copia fundadora publica un fragmento vigente",
        ));
    }
    let piece = exact_object(&pieces[0], "copy.piece", &["kind", "source", "snps"])?;
    if text(piece, "copy.piece", "kind")? != "source" {
        return Err(truth_reject(
            "copy.piece",
            "el fragmento vigente no procede de la referencia",
        ));
    }
    if read_segment_source(
        piece
            .get("source")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("copy.piece", "sin source"))?,
        &identity,
    )? != sequence_length
    {
        return Err(truth_reject(
            "copy.piece",
            "el fragmento vigente no cubre la misma secuencia",
        ));
    }
    let listed = piece
        .get("snps")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy.piece", "snps no es una lista"))?;
    if listed.len() != snps.len() {
        return Err(truth_reject(
            "copy.piece",
            "el fragmento vigente y las modificaciones no traen los mismos SNP",
        ));
    }
    for (position, value) in listed.iter().enumerate() {
        let snp = exact_object(
            value,
            "copy.piece.snp",
            &[
                "event_id",
                "position",
                "reference",
                "alternate",
                "provenance",
            ],
        )?;
        let site_id = text(snp, "copy.piece.snp", "event_id")?;
        let expected = &snps[position];
        if site_id != expected.site_id
            || number(snp, "copy.piece.snp", "position")? + 1 != expected.position_1based
            || single_base(snp, "copy.piece.snp", "reference")? != expected.reference_base
            || single_base(snp, "copy.piece.snp", "alternate")? != expected.alternate_base
        {
            return Err(truth_reject(
                "copy.piece.snp",
                "el fragmento vigente y las modificaciones describen alelos distintos",
            ));
        }
        if read_snp_provenance(snp, &site_id, &population_id, &matrix_source_id)?
            != haplotype_ids[position]
        {
            return Err(truth_reject(
                "copy.piece.snp",
                "el fragmento vigente no repite la procedencia de su modificación",
            ));
        }
    }

    let history = object
        .get("history")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("copy", "history no es una lista"))?;
    if history.len() != snps.len() {
        return Err(truth_reject(
            "copy.history",
            "el historial no coincide con las modificaciones",
        ));
    }
    for (position, value) in history.iter().enumerate() {
        let event = exact_object(
            value,
            "copy.history",
            &["kind", "event_id", "position", "provenance"],
        )?;
        if text(event, "copy.history", "kind")? != "snp" {
            return Err(truth_reject("copy.history", "F5 solo publica hitos SNP"));
        }
        let site_id = text(event, "copy.history", "event_id")?;
        let expected = &snps[position];
        if site_id != expected.site_id
            || number(event, "copy.history", "position")? + 1 != expected.position_1based
        {
            return Err(truth_reject(
                "copy.history",
                "el historial no coincide con las modificaciones",
            ));
        }
        if read_snp_provenance(event, &site_id, &population_id, &matrix_source_id)?
            != haplotype_ids[position]
        {
            return Err(truth_reject(
                "copy.history",
                "el historial no repite la procedencia de su modificación",
            ));
        }
    }

    Ok(TruthCopy {
        identity,
        population_id,
        founder_copy_id,
        sequence_length,
        snps,
        haplotype_ids,
    })
}

// ---------------------------------------------------------------------------
// Recuento MAN-70 desde `truth.jsonl` y nada más
// ---------------------------------------------------------------------------

/// Un locus×subpoblación recontado exclusivamente desde la proyección v2.
#[derive(Clone, Debug, PartialEq)]
pub struct F5RecountedLocus {
    population_id: String,
    block_id: String,
    site_id: String,
    sequence_id: String,
    position_1based: u64,
    reference_base: char,
    alternate_base: char,
    individual_count: u64,
    counts: Vec<u64>,
    frequencies: Vec<f64>,
}

impl F5RecountedLocus {
    #[must_use]
    pub fn population_id(&self) -> &str {
        &self.population_id
    }
    #[must_use]
    pub fn block_id(&self) -> &str {
        &self.block_id
    }
    #[must_use]
    pub fn site_id(&self) -> &str {
        &self.site_id
    }
    #[must_use]
    pub fn sequence_id(&self) -> &str {
        &self.sequence_id
    }
    #[must_use]
    pub const fn position_1based(&self) -> u64 {
        self.position_1based
    }
    #[must_use]
    pub const fn reference_base(&self) -> char {
        self.reference_base
    }
    #[must_use]
    pub const fn alternate_base(&self) -> char {
        self.alternate_base
    }
    #[must_use]
    pub const fn individual_count(&self) -> u64 {
        self.individual_count
    }
    /// Individuos por dosis alternativa `0..=ploidy`.
    #[must_use]
    pub fn counts(&self) -> &[u64] {
        &self.counts
    }
    #[must_use]
    pub fn frequencies(&self) -> &[f64] {
        &self.frequencies
    }
}

/// Recuento MAN-70 independiente: todo lo que sale de `truth.jsonl` y solo de
/// `truth.jsonl`.
#[derive(Clone, Debug, PartialEq)]
pub struct F5Man70Recount {
    execution_id: String,
    scenario_id: String,
    reference_id: String,
    ploidy: u8,
    copies: u64,
    genotypes: u64,
    loci: Vec<F5RecountedLocus>,
}

impl F5Man70Recount {
    #[must_use]
    pub fn execution_id(&self) -> &str {
        &self.execution_id
    }
    #[must_use]
    pub fn scenario_id(&self) -> &str {
        &self.scenario_id
    }
    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }
    #[must_use]
    pub const fn ploidy(&self) -> u8 {
        self.ploidy
    }
    #[must_use]
    pub const fn copies(&self) -> u64 {
        self.copies
    }
    #[must_use]
    pub const fn genotypes(&self) -> u64 {
        self.genotypes
    }
    /// Loci ordenados por subpoblación, bloque y `site_id`.
    #[must_use]
    pub fn loci(&self) -> &[F5RecountedLocus] {
        &self.loci
    }
    /// `recount_source` que el manifiesto debe enlazar.
    #[must_use]
    pub const fn recount_source(&self) -> &'static str {
        F5_MAN70_RECOUNT_SOURCE
    }
}

fn case_id_oracle(
    population_id: &str,
    individual_id: &str,
    sequence_id: &str,
    position_1based: u64,
    alternate_base: char,
) -> String {
    let mut hasher = Sha256::new();
    for value in [
        population_id,
        individual_id,
        sequence_id,
        &position_1based.to_string(),
        &format!("REF,SNP:{alternate_base}"),
    ] {
        hasher.update(value.as_bytes());
        hasher.update([0]);
    }
    format!("CASE_{}", &format!("{:x}", hasher.finalize())[..59])
}

/// Estado vivo del recuento: lo mínimo para resolver portadores.
struct RecountState {
    /// Subpoblación y sitios alternativos de cada copia, como índices globales
    /// ordenados. Una sola entrada por copia: es el término proporcional a
    /// `C_total` que el recuento necesita para resolver portadores, y el mismo
    /// que el productor retiene en su repositorio de identidades.
    copies: BTreeMap<CopyKey, (String, Vec<u32>)>,
    site_index: BTreeMap<String, u32>,
    founder_copy_ids: BTreeSet<String>,
}

impl RecountState {
    fn intern(&mut self, site_id: &str) -> u32 {
        if let Some(index) = self.site_index.get(site_id) {
            return *index;
        }
        let index = self.site_index.len() as u32;
        self.site_index.insert(site_id.to_owned(), index);
        index
    }
}

/// Recuenta los genotipos fundadores **leyendo solo** `truth.jsonl` v2.
///
/// Agrupa por el par estructurado `site_id`×`population_id`, exige el censo y
/// la ploidía que la cabecera declara, resuelve cada portador contra la copia
/// publicada y recalcula la dosis desde los portadores, nunca desde
/// `allele_dosages`. No abre ningún TSV, no consulta el manifiesto y no recibe
/// ningún contador, asignación ni acumulador del productor.
///
/// # Errors
///
/// [`F5OracleRejection`] con causa [`F5OracleCause::Truth`] si la proyección
/// no satisface su contrato, [`F5OracleCause::Man70`] si el recuento por locus
/// contradice el censo declarado y [`F5OracleCause::ArtifactMissing`] si el
/// fichero no puede abrirse.
#[allow(clippy::too_many_lines)]
pub fn recount_f5_man70(truth_path: &Path) -> Checked<F5Man70Recount> {
    let file = File::open(truth_path).map_err(|error| {
        reject(
            F5OracleCause::ArtifactMissing,
            F5_TRUTH_PATH,
            "path",
            format!("{}: {error}", truth_path.display()),
        )
    })?;
    let mut reader = BufReader::new(file);
    let mut line = String::new();
    let read = reader
        .read_line(&mut line)
        .map_err(|error| truth_reject("encoding", error.to_string()))?;
    if read == 0 {
        return Err(truth_reject("header", "documento vacío"));
    }
    if !line.ends_with('\n') {
        return Err(truth_reject("eol", "falta el LF final"));
    }
    let header = read_truth_header(line.trim_end_matches('\n'))?;
    let ploidy_u64 = header.populations[0].ploidy;
    if header
        .populations
        .iter()
        .any(|population| population.ploidy != ploidy_u64)
    {
        return Err(truth_reject(
            "populations",
            "las subpoblaciones no declaran la misma ploidía",
        ));
    }
    let ploidy = u8::try_from(ploidy_u64)
        .map_err(|_| truth_reject("populations", "la ploidía declarada no cabe en u8"))?;
    if ploidy == 0 {
        return Err(truth_reject("populations", "la ploidía declarada es cero"));
    }

    let mut state = RecountState {
        copies: BTreeMap::new(),
        site_index: BTreeMap::new(),
        founder_copy_ids: BTreeSet::new(),
    };
    let mut loci: Vec<F5RecountedLocus> = Vec::new();
    // Una entrada por locus, no por genotipo: dos loci de una subpoblación no
    // pueden compartir secuencia y coordenada, porque entonces sus casos serían
    // el mismo y `case_id` dejaría de identificarlos.
    let mut placed_sites: BTreeSet<(String, String, u64)> = BTreeSet::new();
    let mut copies = 0_u64;
    let mut genotypes = 0_u64;
    let mut seen_genotype = false;
    let mut previous_copy: Option<(String, String, String, String, u64)> = None;
    let mut groups: Vec<(String, String)> = Vec::new();
    let mut current_group: Option<(String, String)> = None;
    let mut previous_locus: Option<(String, String)> = None;
    let mut closed_with_lf = true;

    loop {
        line.clear();
        let read = reader
            .read_line(&mut line)
            .map_err(|error| truth_reject("encoding", error.to_string()))?;
        if read == 0 {
            break;
        }
        if !closed_with_lf {
            // Inalcanzable: una línea sin LF solo puede ser la última del fichero,
            // y entonces la lectura siguiente devuelve cero y sale por el `break`. El
            // documento sin LF final lo caza la comprobación homónima del cierre.
            return Err(truth_reject("eol", "falta el LF final"));
        }
        closed_with_lf = line.ends_with('\n');
        let body = line.trim_end_matches('\n');
        if body.is_empty() {
            return Err(truth_reject("row", "línea vacía prohibida"));
        }
        let value: Value = serde_json::from_str(body)
            .map_err(|error| truth_reject("row", format!("registro ilegible: {error}")))?;
        let kind = value
            .get("record_type")
            .and_then(Value::as_str)
            .ok_or_else(|| truth_reject("row", "registro sin record_type"))?;
        match kind {
            "haplotype_copy" => {
                if seen_genotype {
                    return Err(truth_reject(
                        "order",
                        "una copia aparece después de un genotipo",
                    ));
                }
                let copy = read_truth_copy(&value, &header)?;
                let key = (
                    copy.identity.reference_id.clone(),
                    copy.population_id.clone(),
                    copy.identity.individual_id.clone(),
                    copy.identity.sequence_id.clone(),
                    copy.identity.copy_index,
                );
                if previous_copy.as_ref().is_some_and(|last| *last >= key) {
                    return Err(truth_reject(
                        "order",
                        "las copias no llegan por referencia, población, individuo, secuencia y fase",
                    ));
                }
                previous_copy = Some(key);
                if u64::from(ploidy) <= copy.identity.copy_index {
                    return Err(truth_reject(
                        "copy.identity",
                        "copy_index fuera de 0..ploidy-1",
                    ));
                }
                if !state.founder_copy_ids.insert(copy.founder_copy_id.clone()) {
                    return Err(truth_reject("copy", "un founder_copy_id se repite"));
                }
                let mut indices: Vec<u32> = copy
                    .snps
                    .iter()
                    .map(|snp| state.intern(&snp.site_id))
                    .collect();
                indices.sort_unstable();
                if indices.windows(2).any(|pair| pair[0] == pair[1]) {
                    return Err(truth_reject(
                        "copy",
                        "un sitio aparece dos veces en la misma copia",
                    ));
                }
                // **Alcanzable**, y es la única defensa de este mapa: la clave que
                // se compara con `previous_copy` tiene cinco campos, con la
                // subpoblación en segundo lugar, y la identidad que entra en
                // `state.copies` tiene cuatro, **sin** subpoblación. Dos copias de
                // subpoblaciones distintas con el mismo individuo, secuencia y fase
                // hacen crecer la primera y repiten la segunda, y nada ata el
                // identificador de individuo a su subpoblación. Sin esta guarda el
                // mapa sobrescribiría una entrada en silencio.
                if state
                    .copies
                    .insert(copy.identity, (copy.population_id, indices))
                    .is_some()
                {
                    return Err(truth_reject("copy", "una copia aparece dos veces"));
                }
                copies += 1;
            }
            "genotype" => {
                seen_genotype = true;
                let genotype = recount_genotype(&value, &header, &state, ploidy)?;
                let placed = &genotype;
                let group = (placed.population_id.clone(), placed.block_id.clone());
                if current_group.as_ref() != Some(&group) {
                    if groups.contains(&group) {
                        return Err(truth_reject(
                            "order",
                            format!("el grupo {}×{} vuelve a aparecer", group.0, group.1),
                        ));
                    }
                    if groups.last().is_some_and(|last| *last >= group) {
                        return Err(truth_reject(
                            "order",
                            "los genotipos no llegan por subpoblación y bloque ascendentes",
                        ));
                    }
                    groups.push(group.clone());
                    current_group = Some(group);
                    previous_locus = None;
                }
                let key = (placed.site_id.clone(), placed.individual_id.clone());
                if previous_locus.as_ref().is_some_and(|last| *last >= key) {
                    return Err(truth_reject(
                        "order",
                        "dentro de un bloque los genotipos no llegan por site_id e individuo",
                    ));
                }
                let new_site = previous_locus
                    .as_ref()
                    .is_none_or(|last| last.0 != placed.site_id);
                previous_locus = Some(key);
                let population = header
                    .population(&placed.population_id)
                    // Inalcanzable: `recount_genotype` acaba de resolver esta misma
                    // subpoblación contra la misma cabecera y rechaza si falta. Ninguna
                    // prueba la ejercita.
                    .ok_or_else(|| truth_reject("f5_locus", "subpoblación no declarada"))?;
                if new_site {
                    if !placed_key(&mut placed_sites, placed) {
                        return Err(truth_reject(
                            "f5_locus",
                            format!(
                                "la subpoblación {} publica dos loci en {}:{}",
                                placed.population_id, placed.sequence_id, placed.position_1based
                            ),
                        ));
                    }
                    loci.push(F5RecountedLocus {
                        population_id: placed.population_id.clone(),
                        block_id: placed.block_id.clone(),
                        site_id: placed.site_id.clone(),
                        sequence_id: placed.sequence_id.clone(),
                        position_1based: placed.position_1based,
                        reference_base: placed.reference_base,
                        alternate_base: placed.alternate_base,
                        individual_count: population.individual_count,
                        counts: vec![0; usize::from(ploidy) + 1],
                        frequencies: Vec::new(),
                    });
                }
                let locus = loci
                    .last_mut()
                    // Inalcanzable: `loci` solo está vacío antes del primer genotipo, y el
                    // primero siempre abre locus porque `previous_locus` es `None`.
                    .ok_or_else(|| truth_reject("order", "un genotipo sin locus abierto"))?;
                // La comparación de `block_id` de abajo es redundante: el bloque
                // define el grupo, y un genotipo de otro bloque abre un grupo nuevo,
                // reinicia `previous_locus` y por tanto nunca llega aquí con un locus
                // abierto de otro bloque. Ninguna prueba puede aislarla.
                if locus.site_id != placed.site_id
                    || locus.sequence_id != placed.sequence_id
                    || locus.position_1based != placed.position_1based
                    || locus.reference_base != placed.reference_base
                    || locus.alternate_base != placed.alternate_base
                    || locus.block_id != placed.block_id
                {
                    return Err(truth_reject(
                        "f5_locus",
                        format!(
                            "el sitio {} cambia de bloque, secuencia, posición o alelos",
                            placed.site_id
                        ),
                    ));
                }
                locus.counts[usize::try_from(placed.alternate_dosage).map_err(|_| {
                    // Inalcanzable en 64 bits: la dosis está acotada por la ploidía, que ya
                    // cabe en un byte. Se conserva por el tipo del documento.
                    truth_reject("carriers", "la dosis alternativa no cabe en usize")
                })?] += 1;
                genotypes += 1;
            }
            other => {
                return Err(truth_reject(
                    "row",
                    format!("record_type desconocido {other:?}"),
                ));
            }
        }
    }
    if !closed_with_lf {
        return Err(truth_reject("eol", "falta el LF final"));
    }
    if copies != header.copies || genotypes != header.genotypes {
        return Err(truth_reject(
            "counts",
            format!(
                "la cabecera declara {}/{} y el documento trae {copies}/{genotypes}",
                header.copies, header.genotypes
            ),
        ));
    }
    for locus in &mut loci {
        let total: u64 = locus.counts.iter().sum();
        if total != locus.individual_count {
            return Err(reject(
                F5OracleCause::Man70,
                F5_TRUTH_PATH,
                "recount",
                format!(
                    "el locus {}×{} reconta {total} individuos y la cabecera declara {}",
                    locus.site_id, locus.population_id, locus.individual_count
                ),
            ));
        }
        let denominator = locus.individual_count as f64;
        locus.frequencies = locus
            .counts
            .iter()
            .map(|count| *count as f64 / denominator)
            .collect();
    }
    // Un sitio no puede repetirse dentro de una subpoblación: los grupos no se
    // reabren y dentro de cada bloque el `site_id` crece estrictamente, de modo
    // que basta comprobar la clave publicada entre loci consecutivos.
    if loci.windows(2).any(|pair| {
        (&pair[0].population_id, &pair[0].block_id, &pair[0].site_id)
            >= (&pair[1].population_id, &pair[1].block_id, &pair[1].site_id)
    }) {
        // Inalcanzable: lo implican las dos comprobaciones de orden de arriba —los
        // grupos no se reabren y dentro de cada bloque el `site_id` crece—. Ninguna
        // prueba la ejercita.
        return Err(truth_reject(
            "order",
            "los loci recontados no crecen por subpoblación, bloque y sitio",
        ));
    }
    Ok(F5Man70Recount {
        execution_id: header.execution_id,
        scenario_id: header.scenario_id,
        reference_id: header.reference_id,
        ploidy,
        copies,
        genotypes,
        loci,
    })
}

/// Genotipo ya validado y con su dosis alternativa recalculada.
struct PlacedGenotype {
    population_id: String,
    block_id: String,
    site_id: String,
    sequence_id: String,
    individual_id: String,
    position_1based: u64,
    reference_base: char,
    alternate_base: char,
    alternate_dosage: u64,
}

#[allow(clippy::too_many_lines)]
fn placed_key(placed_sites: &mut BTreeSet<(String, String, u64)>, placed: &PlacedGenotype) -> bool {
    placed_sites.insert((
        placed.population_id.clone(),
        placed.sequence_id.clone(),
        placed.position_1based,
    ))
}

fn recount_genotype(
    value: &Value,
    header: &TruthHeader,
    state: &RecountState,
    ploidy: u8,
) -> Checked<PlacedGenotype> {
    let object = exact_object(
        value,
        "genotype",
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
    if text(object, "genotype", "case_kind")? != "snp" {
        return Err(truth_reject("genotype", "F5 solo publica casos SNP"));
    }
    let reference_id = text(object, "genotype", "reference_id")?;
    if reference_id != header.reference_id {
        return Err(truth_reject(
            "genotype",
            "un genotipo nombra otra referencia que la cabecera",
        ));
    }
    let individual_id = text(object, "genotype", "individual_id")?;
    let locus = exact_object(
        object
            .get("locus")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("genotype", "sin locus"))?,
        "genotype.locus",
        &["sequence_id", "start", "end"],
    )?;
    let start = number(locus, "genotype.locus", "start")?;
    let end = number(locus, "genotype.locus", "end")?;
    if end != start + 1 {
        return Err(truth_reject(
            "genotype.locus",
            "el locus SNP no mide una base",
        ));
    }
    let sequence_id = text(locus, "genotype.locus", "sequence_id")?;
    let f5_locus = exact_object(
        object
            .get("f5_locus")
            // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
            // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
            .ok_or_else(|| truth_reject("genotype", "sin f5_locus"))?,
        "genotype.f5_locus",
        &["population_id", "block_id", "site_id", "ref", "alt"],
    )?;
    let population_id = text(f5_locus, "genotype.f5_locus", "population_id")?;
    let block_id = text(f5_locus, "genotype.f5_locus", "block_id")?;
    let site_id = text(f5_locus, "genotype.f5_locus", "site_id")?;
    let reference_base = single_base(f5_locus, "genotype.f5_locus", "ref")?;
    let alternate_base = single_base(f5_locus, "genotype.f5_locus", "alt")?;
    if reference_base == alternate_base {
        return Err(truth_reject(
            "genotype.f5_locus",
            "f5_locus repite el alelo",
        ));
    }
    let population = header.population(&population_id).ok_or_else(|| {
        truth_reject(
            "genotype.f5_locus",
            format!("subpoblación {population_id} no declarada"),
        )
    })?;
    let site_slot = state.site_index.get(&site_id).copied();

    let alleles = object
        .get("alleles")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("genotype", "alleles no es una lista"))?;
    let mut alternate_index = None;
    for (index, allele) in alleles.iter().enumerate() {
        let kind = allele
            .get("kind")
            .and_then(Value::as_str)
            .ok_or_else(|| truth_reject("genotype.allele", "alelo sin kind"))?;
        match kind {
            "reference" => {
                exact_object(allele, "genotype.allele", &["kind"])?;
            }
            "snp" => {
                let snp = exact_object(allele, "genotype.allele", &["kind", "base"])?;
                if single_base(snp, "genotype.allele", "base")? != alternate_base {
                    return Err(truth_reject(
                        "genotype.allele",
                        "el alelo tipado no es el ALT del locus",
                    ));
                }
                if alternate_index.replace(index).is_some() {
                    return Err(truth_reject(
                        "genotype.allele",
                        "dos alelos alternativos en un sitio bialélico",
                    ));
                }
            }
            other => {
                return Err(truth_reject(
                    "genotype.allele",
                    format!("alelo {other:?} ajeno a F5"),
                ));
            }
        }
    }
    let alternate_index = alternate_index.ok_or_else(|| {
        truth_reject(
            "genotype.allele",
            "el sitio no publica su alelo alternativo",
        )
    })?;

    let declared_dosages: Vec<u64> = object
        .get("allele_dosages")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("genotype", "allele_dosages no es una lista"))?
        .iter()
        .map(|value| {
            value
                .as_u64()
                .ok_or_else(|| truth_reject("genotype", "dosis no entera"))
        })
        .collect::<Checked<_>>()?;
    if declared_dosages.len() != alleles.len() {
        return Err(truth_reject(
            "genotype",
            "las dosis no acompañan a los alelos",
        ));
    }

    // La dosis se recalcula **desde los portadores**, nunca desde
    // `allele_dosages`: ese campo es lo que se juzga, no la fuente.
    let mut from_carriers = vec![0_u64; alleles.len()];
    let mut carrier_copies: Vec<CopyKey> = Vec::with_capacity(usize::from(ploidy));
    for value in object
        .get("carriers")
        .and_then(Value::as_array)
        .ok_or_else(|| truth_reject("genotype", "carriers no es una lista"))?
    {
        let carrier = exact_object(value, "genotype.carrier", &["allele_index", "copy"])?;
        let index = usize::try_from(number(carrier, "genotype.carrier", "allele_index")?)
            // Inalcanzable en 64 bits: `allele_index` lo acota la lista de alelos, que
            // ya está en memoria. Se conserva por el tipo del documento.
            .map_err(|_| truth_reject("genotype.carrier", "allele_index desbordado"))?;
        if index >= alleles.len() {
            return Err(truth_reject(
                "genotype.carrier",
                "allele_index fuera de la lista de alelos",
            ));
        }
        let copy = read_copy_identity(
            carrier
                .get("copy")
                // Inalcanzable: el `exact_object` de arriba ya exigió esta clave, de modo
                // que `get` nunca devuelve `None`. Ninguna prueba la ejercita.
                .ok_or_else(|| truth_reject("genotype.carrier", "sin copy"))?,
            "genotype.carrier.copy",
        )?;
        if copy.reference_id != reference_id
            || copy.individual_id != individual_id
            || copy.sequence_id != sequence_id
        {
            return Err(truth_reject(
                "genotype.carrier",
                "un portador no es una copia de su individuo",
            ));
        }
        let (membership, sites) = state.copies.get(&copy).ok_or_else(|| {
            truth_reject("genotype.carrier", "un portador no tiene copia publicada")
        })?;
        if *membership != population_id {
            return Err(truth_reject(
                "genotype.carrier",
                "f5_locus no coincide con la pertenencia de un portador",
            ));
        }
        let carried = site_slot.is_some_and(|slot| sites.binary_search(&slot).is_ok());
        if carried != (index == alternate_index) {
            return Err(truth_reject(
                "genotype.carrier",
                format!(
                    "el portador {} de {site_id} contradice los alelos de su copia",
                    copy.copy_index
                ),
            ));
        }
        if carrier_copies.contains(&copy) {
            return Err(truth_reject(
                "genotype.carrier",
                "un portador aparece dos veces",
            ));
        }
        carrier_copies.push(copy);
        from_carriers[index] += 1;
    }
    if carrier_copies.len() != usize::from(ploidy) {
        return Err(truth_reject(
            "genotype.carrier",
            "los portadores no son las copias de la ploidía",
        ));
    }
    if from_carriers != declared_dosages {
        return Err(reject(
            F5OracleCause::Man70,
            F5_TRUTH_PATH,
            "allele_dosages",
            "las dosis publicadas no son las de sus portadores",
        ));
    }
    // Defensa **inalcanzable**: las dosis publicadas ya son las de sus
    // portadores y los portadores ya son exactamente `ploidy`, que la cabecera
    // declara igual en todas las subpoblaciones. Se conserva porque es la
    // afirmación que DTC-089 escribe sobre `allele_dosages`.
    if declared_dosages.iter().sum::<u64>() != population.ploidy {
        return Err(truth_reject(
            "genotype",
            "las dosis no suman la ploidía declarada",
        ));
    }

    let declared_case = text(object, "genotype", "case_id")?;
    let expected_case = case_id_oracle(
        &population_id,
        &individual_id,
        &sequence_id,
        start + 1,
        alternate_base,
    );
    if declared_case != expected_case {
        return Err(truth_reject(
            "genotype.case_id",
            "el case_id publicado no es el que DTC-061 deriva de su caso",
        ));
    }

    Ok(PlacedGenotype {
        population_id,
        block_id,
        site_id,
        sequence_id,
        individual_id,
        position_1based: start + 1,
        reference_base,
        alternate_base,
        alternate_dosage: from_carriers[alternate_index],
    })
}

// ---------------------------------------------------------------------------
// Reconstrucción de los siete TSV
// ---------------------------------------------------------------------------

/// Los seis contadores del mejor intento, que la matriz publicada determina.
#[derive(Clone, Copy, Debug, PartialEq)]
struct BestCounters {
    duplicate_rows: u64,
    bins_outside: u64,
    signs_failed: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
    max_marginal_error: f64,
}

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

/// Geometría reconstruida desde `f5-block-layout.tsv` y `f5-sites.tsv`.
struct Geometry {
    blocks: Vec<BlockRow>,
    /// `block_id` -> `block_index`, que solo el ordinal de fila publica.
    block_index_of: BTreeMap<String, u64>,
    /// `block_id` en orden ascendente: la clave de los seis ficheros que dicen
    /// «bloque».
    sorted_block_ids: Vec<String>,
    sites: Vec<Vec<SiteRow>>,
    layout_rows: u64,
    site_rows: u64,
}

/// Unidad bloque×subpoblación reconstruida desde los bytes.
struct Unit {
    population_id: String,
    block_id: String,
    block_index: u64,
    haplotype_ids: Vec<String>,
    /// Matriz `H*S` en bits, con `words` palabras por fila.
    bits: Vec<u64>,
    words: usize,
    haplotypes: usize,
    /// Pesos efectivos **de esta unidad**: el prefijo de rango del vector
    /// declarado con la altura de su matriz, renormalizado por el propio
    /// oráculo. Sin distribución declarada es el vector completo, literal.
    weights: Vec<f64>,
    sites: usize,
    /// `attempt_index` de la fila `selected_best`.
    best_attempt: Option<u64>,
    /// Los seis contadores que la fila `selected_best` publica y que son
    /// funciones puras de la matriz ganadora.
    best_counters: Option<BestCounters>,
    /// Polaridad publicada por sitio, en el orden de `site_index`.
    polarity: Vec<String>,
    /// `target_alt_af` y `realized_alt_af` por `site_index`.
    target_alt_af: Vec<f64>,
    realized_alt_af: Vec<f64>,
    /// Rango haplotípico de cada asiento, en el orden `(individuo, copia)`.
    ranks: Vec<u32>,
    /// Recuento de copias con el alelo alternativo, por `site_index`.
    alt_copies: Vec<u64>,
    /// Recuento de copias asignadas a cada rango: la cuota realizada.
    rank_counts: Vec<u64>,
    /// Individuos con las dos copias en rangos distintos, para F realizado.
    heterozygous: u64,
    individuals: u64,
}

impl Unit {
    fn bit(&self, haplotype: usize, site: usize) -> bool {
        let position = haplotype * self.words * 64 + site;
        self.bits[position / 64] >> (position % 64) & 1 == 1
    }

    fn set_bit(&mut self, haplotype: usize, site: usize) {
        let position = haplotype * self.words * 64 + site;
        self.bits[position / 64] |= 1 << (position % 64);
    }
}

fn read_layout(
    payload_root: &Path,
    oracle: &IdOracle,
    subseeds: &SubseedOracle,
    bounds: &F5DeclaredBounds,
) -> Checked<Geometry> {
    let path = F5_TSV_CONTRACTS[BLOCK_LAYOUT].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, BLOCK_LAYOUT)?;
    let mut blocks: Vec<BlockRow> = Vec::new();
    let mut block_index_of = BTreeMap::new();
    let mut sequence_order: Vec<String> = Vec::new();
    let mut previous: Option<(usize, u64, u64)> = None;
    let mut placed_bases = 0_u64;
    while cursor.advance()? {
        let block_index = cursor.rows - 1;
        let start = cursor.u64(F5OracleCause::BlockLayout, 2)?;
        let end = cursor.u64(F5OracleCause::BlockLayout, 3)?;
        let length = cursor.u64(F5OracleCause::BlockLayout, 4)?;
        placed_bases = placed_bases.checked_add(length).ok_or_else(|| {
            resource_reject(path, "la suma de longitudes del layout desborda u64")
        })?;
        if start == 0 || start > end || length != end - start + 1 {
            return Err(reject(
                F5OracleCause::BlockLayout,
                path,
                "geometry",
                format!(
                    "la fila {block_index} no satisface length_bp=end-start+1 con coordenadas positivas"
                ),
            ));
        }
        let sequence_id = cursor.field(1).to_owned();
        let ordinal = match sequence_order.iter().position(|id| *id == sequence_id) {
            Some(position) => position,
            None => {
                sequence_order.push(sequence_id.clone());
                sequence_order.len() - 1
            }
        };
        let key = (ordinal, start, end);
        // Orden canónico y, dentro de una secuencia, ausencia de solape: dos
        // bloques que compartieran una base darían a un mismo individuo dos
        // loci en la misma coordenada, y con ellos dos casos de verdad
        // indistinguibles.
        if previous.is_some_and(|last| last >= key)
            || previous.is_some_and(|(last_ordinal, _, last_end)| {
                last_ordinal == ordinal && last_end >= start
            })
        {
            return Err(reject(
                F5OracleCause::BlockLayout,
                path,
                "order",
                format!("la fila {block_index} rompe el orden FASTA, inicio y fin sin solape"),
            ));
        }
        previous = Some(key);
        // El ordinal de fila es la única vía por la que un lector recupera
        // `block_index`, componente del material del `BLK_`. La geometría y su
        // orden se comprueban antes para que una fila movida acuse al layout y
        // un identificador dañado acuse a la identidad.
        if cursor.field(0) != oracle.block(block_index)? {
            return Err(reject(
                F5OracleCause::Identity,
                path,
                "block_id",
                format!("la fila {block_index} no rederiva su BLK_ desde el ordinal del layout"),
            ));
        }
        let block_id = cursor.field(0).to_owned();
        // Inalcanzable: el `BLK_` de cada fila se acaba de exigir igual al que se
        // deriva de su ordinal, y el ordinal es distinto en cada fila, de modo que
        // dos filas con el mismo identificador solo podrían darse por una colisión
        // de 128 bits de BLAKE3. Se conserva porque `block_index_of` no debe
        // sobrescribir una entrada en silencio.
        if block_index_of
            .insert(block_id.clone(), block_index)
            .is_some()
        {
            return Err(reject(
                F5OracleCause::BlockLayout,
                path,
                "block_id",
                "el layout repite un identificador de bloque",
            ));
        }
        blocks.push(BlockRow {
            block_id,
            sequence_id,
            start_1based: start,
            end_1based: end,
        });
    }
    if blocks.is_empty() {
        return Err(reject(
            F5OracleCause::BlockLayout,
            path,
            "row",
            "el layout no coloca ningún bloque",
        ));
    }
    if blocks.len() as u64 != bounds.block_count {
        return Err(resource_reject(
            path,
            format!(
                "el layout publica {} bloques y la admisión selló {}",
                blocks.len(),
                bounds.block_count
            ),
        ));
    }
    // `eligible_bases` no se puede recomputar sin abrir la referencia, pero sí
    // acotar: todo bloque colocado ocupa bases elegibles, de modo que la suma
    // de sus longitudes es una cota inferior observable. La superior —las
    // bases que el escaneo recorrió— la comprueba la puerta física.
    if placed_bases > bounds.eligible_bases {
        return Err(resource_reject(
            path,
            format!(
                "el layout coloca {placed_bases} bases y el manifiesto declara {} elegibles",
                bounds.eligible_bases
            ),
        ));
    }
    let layout_rows = cursor.rows;

    let path = F5_TSV_CONTRACTS[SITES].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, SITES)?;
    let mut sites: Vec<Vec<SiteRow>> = (0..blocks.len()).map(|_| Vec::new()).collect();
    let mut previous_block: Option<String> = None;
    while cursor.advance()? {
        let block_id = cursor.field(0);
        let block_index = *block_index_of.get(block_id).ok_or_else(|| {
            reject(
                F5OracleCause::Site,
                path,
                "block_id",
                format!("{block_id} no existe en el layout"),
            )
        })?;
        if previous_block.as_deref() != Some(block_id) {
            if previous_block
                .as_deref()
                .is_some_and(|last| last >= block_id)
            {
                return Err(reject(
                    F5OracleCause::Site,
                    path,
                    "order",
                    "los bloques no llegan por block_id ascendente",
                ));
            }
            previous_block = Some(block_id.to_owned());
        }
        let slot = &mut sites[block_index as usize];
        let site_index = cursor.u64(F5OracleCause::Site, 1)?;
        if site_index != slot.len() as u64 {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "site_index",
                "site_index no es consecutivo desde 0",
            ));
        }
        if cursor.field(2) != oracle.site(block_index, site_index)? {
            return Err(reject(
                F5OracleCause::Identity,
                path,
                "site_id",
                "la fila no rederiva su SITE_",
            ));
        }
        let block = &blocks[block_index as usize];
        if cursor.field(3) != block.sequence_id {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "sequence_id",
                "el sitio no nombra la secuencia de su bloque",
            ));
        }
        let position = cursor.u64(F5OracleCause::Site, 4)?;
        if position < block.start_1based || position > block.end_1based {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "position_1based",
                "el sitio cae fuera de su bloque",
            ));
        }
        if slot
            .last()
            .is_some_and(|last| last.position_1based >= position)
        {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "order",
                "las posiciones no crecen dentro del bloque",
            ));
        }
        let reference_base = cursor.base(F5OracleCause::Site, 5)?;
        let alternate_base = cursor.base(F5OracleCause::Site, 6)?;
        if reference_base == alternate_base {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "alleles",
                "un sitio bialélico repite su alelo",
            ));
        }
        // El ALT no se cree: lo sortea `site_alt` entre las tres bases
        // distintas de REF y sus tres componentes están publicados.
        let expected_alternate =
            subseeds.alternate_base(&block.block_id, cursor.field(2), reference_base)?;
        if alternate_base != expected_alternate {
            return Err(reject(
                F5OracleCause::Site,
                path,
                "alt",
                format!(
                    "el sitio {} publica {alternate_base} y su selector normativo sortea {expected_alternate}",
                    cursor.field(2)
                ),
            ));
        }
        if slot.len() as u64 >= bounds.site_bound_per_block {
            return Err(resource_reject(
                path,
                format!(
                    "el bloque {} publica más de {} sitios",
                    block.block_id, bounds.site_bound_per_block
                ),
            ));
        }
        slot.push(SiteRow {
            site_id: cursor.field(2).to_owned(),
            position_1based: position,
            reference_base,
            alternate_base,
        });
    }
    if sites.iter().any(Vec::is_empty) {
        return Err(reject(
            F5OracleCause::Site,
            path,
            "row",
            "un bloque no publica ningún sitio",
        ));
    }
    let site_rows = cursor.rows;
    let mut sorted_block_ids: Vec<String> =
        blocks.iter().map(|block| block.block_id.clone()).collect();
    sorted_block_ids.sort();
    Ok(Geometry {
        blocks,
        block_index_of,
        sorted_block_ids,
        sites,
        layout_rows,
        site_rows,
    })
}

/// Crea las `N_pop * B` unidades vacías, en la clave publicada.
fn build_units(declaration: &F5PackageDeclaration, geometry: &Geometry) -> Vec<Unit> {
    let mut units = Vec::with_capacity(declaration.subpopulations.len() * geometry.blocks.len());
    for subpopulation in &declaration.subpopulations {
        for block_id in &geometry.sorted_block_ids {
            let block_index = geometry.block_index_of[block_id];
            let sites = geometry.sites[block_index as usize].len();
            units.push(Unit {
                population_id: subpopulation.population_id.clone(),
                block_id: block_id.clone(),
                block_index,
                haplotype_ids: Vec::new(),
                bits: Vec::new(),
                words: 0,
                haplotypes: 0,
                weights: Vec::new(),
                sites,
                best_attempt: None,
                best_counters: None,
                polarity: Vec::new(),
                target_alt_af: Vec::new(),
                realized_alt_af: Vec::new(),
                ranks: Vec::new(),
                alt_copies: vec![0; sites],
                rank_counts: Vec::new(),
                heterozygous: 0,
                individuals: 0,
            });
        }
    }
    units
}

/// Recorre un TSV cuya clave empieza por `population_id` y `block_id` y
/// devuelve el índice de unidad de cada fila, exigiendo que las unidades
/// lleguen en la clave publicada y sin repetirse.
struct UnitWalk {
    current: Option<usize>,
    visited: Vec<bool>,
}

impl UnitWalk {
    fn new(units: usize) -> Self {
        Self {
            current: None,
            visited: vec![false; units],
        }
    }

    fn visit(&mut self, cause: F5OracleCause, path: &'static str, index: usize) -> Checked<bool> {
        if self.current == Some(index) {
            return Ok(false);
        }
        if self.current.is_some_and(|last| last >= index) || self.visited[index] {
            return Err(reject(
                cause,
                path,
                "order",
                "las unidades no llegan por subpoblación y block_id ascendentes",
            ));
        }
        self.visited[index] = true;
        self.current = Some(index);
        Ok(true)
    }

    fn finish(&self, cause: F5OracleCause, path: &'static str) -> Checked<()> {
        if self.visited.iter().any(|seen| !seen) {
            return Err(reject(
                cause,
                path,
                "row",
                "falta alguna unidad bloque×subpoblación",
            ));
        }
        Ok(())
    }
}

fn unit_index(
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    cause: F5OracleCause,
    path: &'static str,
    population_id: &str,
    block_id: &str,
) -> Checked<usize> {
    let population = declaration
        .subpopulations
        .iter()
        .position(|subpopulation| subpopulation.population_id == population_id)
        .ok_or_else(|| {
            reject(
                cause,
                path,
                "population_id",
                format!("{population_id} no es una subpoblación declarada"),
            )
        })?;
    let block = geometry
        .sorted_block_ids
        .iter()
        .position(|candidate| candidate == block_id)
        .ok_or_else(|| {
            reject(
                cause,
                path,
                "block_id",
                format!("{block_id} no existe en el layout"),
            )
        })?;
    Ok(population * geometry.sorted_block_ids.len() + block)
}

/// `f5-attempts.tsv`: contadores acotados, orden y un único mejor intento.
fn read_attempts(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    units: &mut [Unit],
) -> Checked<u64> {
    let path = F5_TSV_CONTRACTS[ATTEMPTS].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, ATTEMPTS)?;
    let mut walk = UnitWalk::new(units.len());
    let mut previous_attempt: Option<u64> = None;
    let mut best_count = 0_u64;
    let mut index = usize::MAX;
    let mut max_attempt_index = 0_u64;
    while cursor.advance()? {
        let next = unit_index(
            declaration,
            geometry,
            F5OracleCause::Attempt,
            path,
            cursor.field(0),
            cursor.field(1),
        )?;
        if walk.visit(F5OracleCause::Attempt, path, next)? {
            if index != usize::MAX && best_count != 1 {
                return Err(reject(
                    F5OracleCause::Attempt,
                    path,
                    "selected_best",
                    format!("una unidad marca {best_count} intentos como mejor"),
                ));
            }
            best_count = 0;
            previous_attempt = None;
            index = next;
        }
        let attempt = cursor.u64(F5OracleCause::Attempt, 2)?;
        if previous_attempt.is_some_and(|last| last >= attempt) {
            return Err(reject(
                F5OracleCause::Attempt,
                path,
                "order",
                "los intentos de una unidad no llegan por índice ascendente",
            ));
        }
        previous_attempt = Some(attempt);
        max_attempt_index = max_attempt_index.max(attempt);
        if attempt >= declaration.bounds.max_haplotype_attempts {
            return Err(resource_reject(
                path,
                format!(
                    "attempt_index={attempt} no cabe en max_haplotype_attempts={}",
                    declaration.bounds.max_haplotype_attempts
                ),
            ));
        }
        let proposals = cursor.u64(F5OracleCause::Attempt, 3)?;
        if proposals > declaration.bounds.max_swap_proposals {
            return Err(resource_reject(
                path,
                format!(
                    "proposals={proposals} no cabe en max_swap_proposals={}",
                    declaration.bounds.max_swap_proposals
                ),
            ));
        }
        let accepted = cursor.u64(F5OracleCause::Attempt, 4)?;
        if accepted > proposals {
            return Err(reject(
                F5OracleCause::Attempt,
                path,
                "accepted",
                "un intento acepta más propuestas de las que hizo",
            ));
        }
        let duplicate_rows = cursor.u64(F5OracleCause::Attempt, 5)?;
        let bins_outside = cursor.u64(F5OracleCause::Attempt, 6)?;
        let signs_failed = cursor.u64(F5OracleCause::Attempt, 7)?;
        let max_r2_error = cursor.f64(F5OracleCause::Attempt, 8)?;
        let sum_r2_error = cursor.f64(F5OracleCause::Attempt, 9)?;
        let max_marginal_error = cursor.f64(F5OracleCause::Attempt, 10)?;
        match cursor.field(11) {
            "true" => {
                best_count += 1;
                units[next].best_attempt = Some(attempt);
                units[next].best_counters = Some(BestCounters {
                    duplicate_rows,
                    bins_outside,
                    signs_failed,
                    max_r2_error,
                    sum_r2_error,
                    max_marginal_error,
                });
            }
            "false" => {}
            other => {
                return Err(reject(
                    F5OracleCause::Attempt,
                    path,
                    "selected_best",
                    format!("selected_best {other:?} ajeno"),
                ));
            }
        }
    }
    if index != usize::MAX && best_count != 1 {
        return Err(reject(
            F5OracleCause::Attempt,
            path,
            "selected_best",
            format!("una unidad marca {best_count} intentos como mejor"),
        ));
    }
    walk.finish(F5OracleCause::Attempt, path)?;
    // `attempt_summary` del manifiesto: dos de sus tres términos son funciones
    // puras de estas filas y aquí se contrastan con ellas. El tercero,
    // `layout_attempts_used`, no es observable desde ningún TSV publicado y
    // queda declarado; se dice para que no se confunda con un olvido.
    let summary = declaration.attempt_summary;
    if summary.haplotype_attempts_total != cursor.rows {
        return Err(reject(
            F5OracleCause::Attempt,
            path,
            "attempt_summary.haplotype_attempts_total",
            format!(
                "el manifiesto declara {} intentos y el fichero trae {}",
                summary.haplotype_attempts_total, cursor.rows
            ),
        ));
    }
    if summary.max_attempt_index != max_attempt_index {
        return Err(reject(
            F5OracleCause::Attempt,
            path,
            "attempt_summary.max_attempt_index",
            format!(
                "el manifiesto declara {} y el mayor índice publicado es {max_attempt_index}",
                summary.max_attempt_index
            ),
        ));
    }
    Ok(cursor.rows)
}

/// `f5-haplotypes.tsv`: rangos, `q`, `H*S` y la rederivación de cada `HAP_`.
/// Pesos efectivos de un bloque con `haplotype_count` filas.
///
/// El oráculo los rehace con su propio Neumaier: toma el prefijo de rango del
/// vector declarado y lo divide por su suma compensada. Con el vector
/// completo no divide —el vector ya está normalizado y dividir movería sus
/// últimos bits—, que es exactamente la regla que el productor aplica.
fn oracle_block_weights(declared: &[f64], haplotype_count: usize) -> Checked<Vec<f64>> {
    let path = F5_TSV_CONTRACTS[HAPLOTYPES].path();
    if haplotype_count < 2 || haplotype_count > declared.len() {
        return Err(reject(
            F5OracleCause::Haplotype,
            path,
            "haplotype_rank",
            format!(
                "un bloque publica {haplotype_count} filas y el perfil declara {} pesos",
                declared.len()
            ),
        ));
    }
    let prefix = &declared[..haplotype_count];
    if haplotype_count == declared.len() {
        return Ok(prefix.to_vec());
    }
    let total = neumaier_total(prefix.iter().copied());
    if !total.is_finite() || total <= 0.0 {
        return Err(reject(
            F5OracleCause::Haplotype,
            path,
            "q",
            format!("el prefijo de {haplotype_count} pesos suma {total}"),
        ));
    }
    Ok(prefix.iter().map(|weight| *weight / total).collect())
}

/// Recuento haplotípico que **cada** unidad debe publicar.
///
/// Sin distribución declarada es el número de pesos del perfil, como siempre.
/// Con distribución declarada se rederiva del ámbito `block_haplotype_count`
/// y se contrasta contra la lista publicada: el oráculo no lee el H del
/// paquete y se lo cree, lo sortea él y exige que coincidan los dos.
fn expected_haplotype_counts(
    declaration: &F5PackageDeclaration,
    subseeds: &SubseedOracle,
    units: &[Unit],
) -> Checked<Vec<usize>> {
    let path = F5_TSV_CONTRACTS[HAPLOTYPES].path();
    let Some(distribution) = declaration.haplotype_count_distribution else {
        if !declaration.haplotype_counts.is_empty() {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_count",
                "se publican recuentos haplotípicos sin declarar la distribución que los produce",
            ));
        }
        return Ok(vec![declaration.haplotype_weights.len(); units.len()]);
    };
    if declaration.haplotype_counts.len() != units.len() {
        return Err(reject(
            F5OracleCause::Haplotype,
            path,
            "haplotype_count",
            format!(
                "se declaran {} recuentos haplotípicos y el paquete tiene {} unidades",
                declaration.haplotype_counts.len(),
                units.len()
            ),
        ));
    }
    let declared_by_unit: BTreeMap<(&str, &str), u64> = declaration
        .haplotype_counts
        .iter()
        .map(|row| {
            (
                (row.population_id.as_str(), row.block_id.as_str()),
                row.haplotype_count,
            )
        })
        .collect();
    let mut expected = Vec::with_capacity(units.len());
    for unit in units {
        let redrawn =
            subseeds.block_haplotype_count(&unit.population_id, &unit.block_id, distribution)?;
        let declared = *declared_by_unit
            .get(&(unit.population_id.as_str(), unit.block_id.as_str()))
            .ok_or_else(|| {
                reject(
                    F5OracleCause::Haplotype,
                    path,
                    "haplotype_count",
                    format!(
                        "no se declara recuento haplotípico para {}×{}",
                        unit.population_id, unit.block_id
                    ),
                )
            })?;
        if declared != redrawn {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_count",
                format!(
                    "{}×{} declara H={declared} y el ámbito block_haplotype_count sortea {redrawn}",
                    unit.population_id, unit.block_id
                ),
            ));
        }
        // Inalcanzable en 64 bits: `redrawn` está pinzado en `[h_min, h_max]` y
        // `usize` tiene 64 bits en los objetivos del proyecto. Se conserva por el
        // tipo de la conversión.
        expected.push(usize::try_from(redrawn).map_err(|_| {
            reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_count",
                "el recuento haplotípico no cabe en usize",
            )
        })?);
    }
    Ok(expected)
}

fn read_haplotypes(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    oracle: &IdOracle,
    subseeds: &SubseedOracle,
    units: &mut [Unit],
) -> Checked<u64> {
    let path = F5_TSV_CONTRACTS[HAPLOTYPES].path();
    let declared_weights = &declaration.haplotype_weights;
    // Los dos rechazos siguientes son inalcanzables. `audit_normalizations` corre
    // antes que esta función y exige que `declaration.haplotype_weights` sea,
    // valor a valor, `effective_values` de la auditoría DTC-049 de los pesos; y
    // `check_normalization` exige de esa auditoría que no esté vacía y que su
    // suma de salida esté a distancia de uno menor o igual que la misma
    // tolerancia que se usa abajo. Se conservan porque esta función no se fía del
    // orden de sus llamadores.
    if declared_weights.is_empty() {
        return Err(reject(
            F5OracleCause::Haplotype,
            path,
            "q",
            "la declaración no trae ningún peso haplotípico",
        ));
    }
    let weights_total = neumaier_total(declared_weights.iter().copied());
    if (weights_total - 1.0).abs() > F5_PROFILE_NORMALIZATION_TOLERANCE {
        return Err(reject(
            F5OracleCause::Haplotype,
            path,
            "q",
            format!("los pesos declarados suman {weights_total} y no 1"),
        ));
    }
    let expected_counts = expected_haplotype_counts(declaration, subseeds, units)?;
    for (unit, count) in units.iter_mut().zip(expected_counts.iter().copied()) {
        unit.weights = oracle_block_weights(declared_weights, count)?;
    }
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, HAPLOTYPES)?;
    let mut walk = UnitWalk::new(units.len());
    let mut index = usize::MAX;
    while cursor.advance()? {
        let next = unit_index(
            declaration,
            geometry,
            F5OracleCause::Haplotype,
            path,
            cursor.field(0),
            cursor.field(1),
        )?;
        if walk.visit(F5OracleCause::Haplotype, path, next)? {
            index = next;
            let sites = units[index].sites;
            let words = sites.div_ceil(64).max(1);
            let haplotypes = units[index].weights.len();
            units[index].words = words;
            units[index].bits = vec![0; words * haplotypes];
            units[index].haplotypes = haplotypes;
        }
        let rank = cursor.u64(F5OracleCause::Haplotype, 3)?;
        let unit = &units[index];
        if rank != unit.haplotype_ids.len() as u64 {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_rank",
                "haplotype_rank no es consecutivo desde 0",
            ));
        }
        // Inalcanzable: `rank` acaba de compararse con `haplotype_ids.len()`, que
        // es un `usize`, de modo que cabe. Se conserva por el tipo de la conversión.
        let position = usize::try_from(rank).map_err(|_| {
            reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_rank",
                "el rango no cabe en usize",
            )
        })?;
        if position >= unit.weights.len() {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "haplotype_rank",
                "la unidad publica más filas que el recuento rederivado admite",
            ));
        }
        let q = cursor.f64(F5OracleCause::Haplotype, 4)?;
        if !(q > 0.0 && q <= 1.0) {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "q",
                format!("q={q} no pertenece a (0,1]"),
            ));
        }
        if q != unit.weights[position] {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "q",
                format!(
                    "la unidad publica q={q} en el rango {rank} y el prefijo renormalizado vale {}",
                    unit.weights[position]
                ),
            ));
        }
        let alleles = cursor.field(5);
        if alleles.len() != unit.sites || !alleles.bytes().all(|byte| byte == b'0' || byte == b'1')
        {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "alleles",
                "alleles no son S caracteres 0/1",
            ));
        }
        if cursor.field(2) != oracle.haplotype(cursor.field(0), unit.block_index, rank, alleles)? {
            return Err(reject(
                F5OracleCause::Identity,
                path,
                "haplotype_id",
                "la fila no rederiva su HAP_",
            ));
        }
        let haplotype_id = cursor.field(2).to_owned();
        let carried: Vec<usize> = alleles
            .bytes()
            .enumerate()
            .filter_map(|(site, byte)| (byte == b'1').then_some(site))
            .collect();
        let unit = &mut units[index];
        unit.haplotype_ids.push(haplotype_id);
        for site in carried {
            unit.set_bit(position, site);
        }
    }
    walk.finish(F5OracleCause::Haplotype, path)?;
    for unit in units.iter() {
        if unit.haplotype_ids.len() != unit.weights.len() {
            return Err(reject(
                F5OracleCause::Haplotype,
                path,
                "row",
                format!(
                    "la unidad {}×{} publica {} filas y el recuento rederivado es {}",
                    unit.population_id,
                    unit.block_id,
                    unit.haplotype_ids.len(),
                    unit.weights.len()
                ),
            ));
        }
    }
    Ok(cursor.rows)
}

/// Posiciones de los sitios de cada bloque ordenadas por `site_id`, que es la
/// clave literal de `f5-frequencies.tsv` y la de los genotipos.
fn sites_by_id(geometry: &Geometry) -> Vec<Vec<usize>> {
    geometry
        .sites
        .iter()
        .map(|sites| {
            let mut order: Vec<usize> = (0..sites.len()).collect();
            order.sort_by(|left, right| sites[*left].site_id.cmp(&sites[*right].site_id));
            order
        })
        .collect()
}

/// `Σ_h q_h x_hi` con Neumaier en orden canónico `h`, como fija DTC-089.
fn realized_marginal(unit: &Unit, weights: &[f64], site: usize) -> f64 {
    let mut accumulator = Neumaier::default();
    for (haplotype, weight) in weights.iter().enumerate() {
        if unit.bit(haplotype, site) {
            accumulator.add(*weight);
        }
    }
    accumulator.total()
}

/// `Σ_h q_h x_hi x_hj` con la misma aritmética.
fn joint_marginal(unit: &Unit, weights: &[f64], left: usize, right: usize) -> f64 {
    let mut accumulator = Neumaier::default();
    for (haplotype, weight) in weights.iter().enumerate() {
        if unit.bit(haplotype, left) && unit.bit(haplotype, right) {
            accumulator.add(*weight);
        }
    }
    accumulator.total()
}

/// `f5-frequencies.tsv`: marginal realizada bit a bit, error, y la polaridad
/// rederivada desde el RNG normativo.
fn read_frequencies(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    by_id: &[Vec<usize>],
    subseeds: &SubseedOracle,
    units: &mut [Unit],
) -> Checked<u64> {
    let path = F5_TSV_CONTRACTS[FREQUENCIES].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, FREQUENCIES)?;
    let mut walk = UnitWalk::new(units.len());
    let mut index = usize::MAX;
    let mut expected_site = 0_usize;
    let mut previous_site: Option<String> = None;
    while cursor.advance()? {
        let next = unit_index(
            declaration,
            geometry,
            F5OracleCause::Marginal,
            path,
            cursor.field(0),
            cursor.field(1),
        )?;
        if walk.visit(F5OracleCause::Marginal, path, next)? {
            // Cierre real de la unidad anterior: la publicada debe traer
            // exactamente los sitios de su bloque, ni uno menos.
            if index != usize::MAX {
                ensure_all_sites(
                    path,
                    expected_site,
                    &by_id[units[index].block_index as usize],
                    &units[index],
                )?;
            }
            index = next;
            expected_site = 0;
            previous_site = None;
            units[index].polarity = vec![String::new(); units[index].sites];
            units[index].target_alt_af = vec![0.0; units[index].sites];
            units[index].realized_alt_af = vec![0.0; units[index].sites];
        }
        let unit = &units[index];
        let order = &by_id[unit.block_index as usize];
        if expected_site >= order.len() {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "row",
                "la unidad publica más sitios que su bloque",
            ));
        }
        let site_index = order[expected_site];
        let site = &geometry.sites[unit.block_index as usize][site_index];
        // La clave literal es `(population_id, block_id, site_id)`: no se
        // deriva de `f5-sites.tsv`, se exige creciente.
        if previous_site
            .as_deref()
            .is_some_and(|last| last >= cursor.field(2))
        {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "order",
                "el orden no es population_id, block_id y site_id estrictamente creciente",
            ));
        }
        previous_site = Some(cursor.field(2).to_owned());
        if cursor.field(2) != site.site_id {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "site_id",
                format!(
                    "la unidad publica {} donde su bloque tiene {}",
                    cursor.field(2),
                    site.site_id
                ),
            ));
        }
        let target = cursor.probability(F5OracleCause::Marginal, 4)?;
        let realized = cursor.probability(F5OracleCause::Marginal, 5)?;
        let error = cursor.probability(F5OracleCause::Marginal, 6)?;
        // El objetivo tampoco se cree: lo sortea `marginal_spectrum` y sus
        // cuatro componentes están en esta misma fila.
        let expected_target = subseeds.target_alt_af(
            &units[index].population_id,
            &units[index].block_id,
            &site.site_id,
            &declaration.snp_spectrum,
        )?;
        if target != expected_target {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "target_alt_af",
                format!(
                    "el sitio {} publica un objetivo {target} y su espectro normativo sortea {expected_target}",
                    site.site_id
                ),
            ));
        }
        // Los pesos son los de **esta** unidad: el prefijo renormalizado que
        // `read_haplotypes` fijó al rederivar su recuento haplotípico.
        let recomputed = realized_marginal(&units[index], &units[index].weights, site_index);
        if realized != recomputed {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "realized_alt_af",
                format!(
                    "la marginal publicada {realized} no es la suma canónica {recomputed} de su matriz"
                ),
            ));
        }
        if error != (realized - target).abs() {
            return Err(reject(
                F5OracleCause::Marginal,
                path,
                "abs_error",
                "abs_error no es |realizado-objetivo|",
            ));
        }
        // Inalcanzable: `read_attempts` corre antes, exige que cada unidad que
        // recorre marque exactamente un intento como mejor y `walk.finish` que las
        // recorra todas, de modo que todas tienen mejor intento.
        let best_attempt = units[index].best_attempt.ok_or_else(|| {
            reject(
                F5OracleCause::Attempt,
                path,
                "selected_best",
                "la unidad no publica ningún mejor intento",
            )
        })?;
        let expected_polarity = subseeds.polarity(
            &units[index].population_id,
            &units[index].block_id,
            best_attempt,
            site_index as u64,
        )?;
        if cursor.field(3) != expected_polarity {
            return Err(reject(
                F5OracleCause::Polarity,
                path,
                "polarity",
                format!(
                    "el sitio {} publica {:?} y su Bernoulli normativa da {expected_polarity:?}",
                    site.site_id,
                    cursor.field(3)
                ),
            ));
        }
        let unit = &mut units[index];
        unit.polarity[site_index] = expected_polarity.to_owned();
        unit.target_alt_af[site_index] = target;
        unit.realized_alt_af[site_index] = realized;
        expected_site += 1;
    }
    if index != usize::MAX {
        ensure_all_sites(
            path,
            expected_site,
            &by_id[units[index].block_index as usize],
            &units[index],
        )?;
    }
    walk.finish(F5OracleCause::Marginal, path)?;
    Ok(cursor.rows)
}

/// Una unidad marginal publica exactamente los sitios de su bloque.
///
/// La comprobación anterior —`realized_alt_af.len() != sites`— era una guarda
/// muerta: el vector se dimensiona al visitar la unidad y su longitud nunca
/// podía discrepar. Con ella, borrar la última fila por `site_id` de una unidad
/// solo moría de rebote, y por una causa ajena.
fn ensure_all_sites(
    path: &'static str,
    published: usize,
    order: &[usize],
    unit: &Unit,
) -> Checked<()> {
    if published != order.len() {
        return Err(reject(
            F5OracleCause::Marginal,
            path,
            "row",
            format!(
                "la unidad {}×{} publica {published} sitios y su bloque tiene {}",
                unit.population_id,
                unit.block_id,
                order.len()
            ),
        ));
    }
    Ok(())
}

/// Bins LD recalculados desde la matriz, los pesos, las posiciones y las
/// polaridades, con los tres agregados que `f5-attempts.tsv` publica.
struct LdRecomputation {
    pair_count: Vec<u64>,
    mean_r2: Vec<Option<f64>>,
    abs_error: Vec<Option<f64>>,
    bins_outside: u64,
    signs_failed: u64,
    max_r2_error: f64,
    sum_r2_error: f64,
}

fn recompute_ld(
    unit: &Unit,
    weights: &[f64],
    sites: &[SiteRow],
    bins: &[&F5DeclaredLdBin],
) -> LdRecomputation {
    let mut totals = vec![Neumaier::default(); bins.len()];
    let mut pair_count = vec![0_u64; bins.len()];
    let mut signs_failed = 0_u64;
    for left in 0..unit.sites {
        for right in (left + 1)..unit.sites {
            let p_i = unit.realized_alt_af[left];
            let p_j = unit.realized_alt_af[right];
            let p_ij = joint_marginal(unit, weights, left, right);
            let d = p_ij - p_i * p_j;
            let factor_i = p_i * (1.0 - p_i);
            let factor_j = p_j * (1.0 - p_j);
            if !(factor_i > 0.0 && factor_j > 0.0) {
                continue;
            }
            let r_squared = (d * d / (factor_i * factor_j)).min(1.0);
            let polarities_match = unit.polarity[left] == unit.polarity[right];
            let satisfied = if polarities_match { d > 0.0 } else { d < 0.0 };
            if !satisfied {
                signs_failed += 1;
            }
            let distance = sites[right].position_1based - sites[left].position_1based;
            if let Some(bin) = bins.iter().position(|bin| {
                bin.distance_start_bp <= distance && distance <= bin.distance_end_bp
            }) {
                totals[bin].add(r_squared);
                pair_count[bin] += 1;
            }
        }
    }
    let mut mean_r2 = Vec::with_capacity(bins.len());
    let mut abs_error = Vec::with_capacity(bins.len());
    let mut bins_outside = 0_u64;
    let mut max_r2_error = 0.0_f64;
    let mut error_sum = Neumaier::default();
    for (bin, definition) in bins.iter().enumerate() {
        if pair_count[bin] == 0 {
            mean_r2.push(None);
            abs_error.push(None);
            continue;
        }
        let mean = totals[bin].total() / pair_count[bin] as f64;
        let error = (mean - definition.target_mean_r2).abs();
        if error > definition.max_abs_error {
            bins_outside += 1;
        }
        if error > max_r2_error {
            max_r2_error = error;
        }
        error_sum.add(error);
        mean_r2.push(Some(mean));
        abs_error.push(Some(error));
    }
    LdRecomputation {
        pair_count,
        mean_r2,
        abs_error,
        bins_outside,
        signs_failed,
        max_r2_error,
        sum_r2_error: error_sum.total(),
    }
}

/// `f5-ld-evaluation.tsv`: bins, `NA`, medias y errores contra la
/// recomputación propia, y los agregados del mejor intento.
fn read_ld_evaluation(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    units: &[Unit],
) -> Checked<u64> {
    let path = F5_TSV_CONTRACTS[LD_EVALUATION].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, LD_EVALUATION)?;
    let mut walk = UnitWalk::new(units.len());
    let mut index = usize::MAX;
    let mut bins: Vec<&F5DeclaredLdBin> = Vec::new();
    let mut recomputed = None;
    let mut bin_position = 0_usize;
    let mut previous_start: Option<u64> = None;
    while cursor.advance()? {
        let next = unit_index(
            declaration,
            geometry,
            F5OracleCause::LinkageDisequilibrium,
            path,
            cursor.field(0),
            cursor.field(1),
        )?;
        if walk.visit(F5OracleCause::LinkageDisequilibrium, path, next)? {
            if index != usize::MAX {
                ensure_all_bins(path, bin_position, &bins, &units[index])?;
            }
            index = next;
            bins = declaration
                .ld_bins
                .iter()
                .filter(|bin| bin.population_id == units[index].population_id)
                .collect();
            if bins.is_empty() {
                return Err(reject(
                    F5OracleCause::LinkageDisequilibrium,
                    path,
                    "bins",
                    format!(
                        "el perfil no declara bins LD para {}",
                        units[index].population_id
                    ),
                ));
            }
            let block_index = units[index].block_index as usize;
            recomputed = Some(recompute_ld(
                &units[index],
                &units[index].weights,
                &geometry.sites[block_index],
                &bins,
            ));
            bin_position = 0;
            previous_start = None;
        }
        // Inalcanzable: la primera fila abre siempre una unidad —`walk.visit`
        // devuelve `true` cuando no hay unidad actual— y `recomputed` se asigna al
        // abrirla; las siguientes filas heredan ese valor.
        let evaluation = recomputed.as_ref().ok_or_else(|| {
            reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "row",
                "una fila LD sin unidad abierta",
            )
        })?;
        if bin_position >= bins.len() {
            return Err(reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "bins",
                "la unidad publica más bins que su perfil",
            ));
        }
        let definition = bins[bin_position];
        let start = cursor.u64(F5OracleCause::LinkageDisequilibrium, 2)?;
        let end = cursor.u64(F5OracleCause::LinkageDisequilibrium, 3)?;
        if previous_start.is_some_and(|last| last >= start) {
            return Err(reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "order",
                "los bins no llegan por inicio ascendente",
            ));
        }
        previous_start = Some(start);
        if start != definition.distance_start_bp || end != definition.distance_end_bp {
            return Err(reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "bins",
                "el bin publicado no es el que el perfil declara en esa posición",
            ));
        }
        let pairs = cursor.u64(F5OracleCause::LinkageDisequilibrium, 5)?;
        if pairs != evaluation.pair_count[bin_position] {
            return Err(reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "pair_count",
                format!(
                    "el bin publica {pairs} parejas y la recomputación halla {}",
                    evaluation.pair_count[bin_position]
                ),
            ));
        }
        let target = cursor.probability(F5OracleCause::LinkageDisequilibrium, 7)?;
        let tolerance = cursor.probability(F5OracleCause::LinkageDisequilibrium, 8)?;
        if target != definition.target_mean_r2 || tolerance != definition.max_abs_error {
            return Err(reject(
                F5OracleCause::LinkageDisequilibrium,
                path,
                "bins",
                "el objetivo o la tolerancia publicados no son los del perfil",
            ));
        }
        match cursor.field(4) {
            "evaluated" => {
                if pairs == 0 {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "status",
                        "un bin evaluado no puede tener cero parejas",
                    ));
                }
                if cursor.field(6) == NOT_AVAILABLE || cursor.field(9) == NOT_AVAILABLE {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "status",
                        "NA fuera de un bin no aplicable",
                    ));
                }
                let mean = cursor.probability(F5OracleCause::LinkageDisequilibrium, 6)?;
                let error = cursor.probability(F5OracleCause::LinkageDisequilibrium, 9)?;
                // Inalcanzable: `pairs` acaba de igualarse a las parejas de la
                // recomputación y arriba se exigió que no fueran cero; con parejas, la
                // recomputación tiene media.
                let expected_mean = evaluation.mean_r2[bin_position].ok_or_else(|| {
                    reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "mean_r2",
                        "el bin publica media y la recomputación no halla ninguna pareja",
                    )
                })?;
                if mean != expected_mean {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "mean_r2",
                        format!("el bin publica {mean} y la recomputación da {expected_mean}"),
                    ));
                }
                let expected_error = evaluation.abs_error[bin_position].unwrap_or(f64::NAN);
                if error != (mean - target).abs() || error != expected_error {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "abs_error",
                        "abs_error no es |media-objetivo| recalculado",
                    ));
                }
                if error > tolerance {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "abs_error",
                        "un bin publicado excede su tolerancia",
                    ));
                }
            }
            "not_applicable" => {
                if pairs != 0 {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "status",
                        "un bin no aplicable tiene parejas",
                    ));
                }
                if cursor.field(6) != NOT_AVAILABLE || cursor.field(9) != NOT_AVAILABLE {
                    return Err(reject(
                        F5OracleCause::LinkageDisequilibrium,
                        path,
                        "status",
                        "un bin no aplicable debe usar NA",
                    ));
                }
            }
            other => {
                return Err(reject(
                    F5OracleCause::LinkageDisequilibrium,
                    path,
                    "status",
                    format!("status {other:?} ajeno"),
                ));
            }
        }
        bin_position += 1;
        if bin_position == bins.len() {
            // Inalcanzable: por la misma razón que en `read_frequencies`, todas las
            // unidades tienen mejor intento cuando llegan aquí.
            let counters = units[index].best_counters.ok_or_else(|| {
                reject(
                    F5OracleCause::Attempt,
                    path,
                    "selected_best",
                    "la unidad no publica ningún mejor intento",
                )
            })?;
            if counters.bins_outside != evaluation.bins_outside
                || counters.max_r2_error != evaluation.max_r2_error
                || counters.sum_r2_error != evaluation.sum_r2_error
            {
                return Err(reject(
                    F5OracleCause::Attempt,
                    F5_TSV_CONTRACTS[ATTEMPTS].path(),
                    "score",
                    format!(
                        "el mejor intento publica bins_outside={} max_r2_error={} sum_r2_error={} y la matriz ganadora da {} {} {}",
                        counters.bins_outside,
                        counters.max_r2_error,
                        counters.sum_r2_error,
                        evaluation.bins_outside,
                        evaluation.max_r2_error,
                        evaluation.sum_r2_error
                    ),
                ));
            }
            if counters.signs_failed != evaluation.signs_failed {
                return Err(reject(
                    F5OracleCause::SignRule,
                    F5_TSV_CONTRACTS[ATTEMPTS].path(),
                    "signs_failed",
                    format!(
                        "el mejor intento publica {} signos fallidos y la regla de DTC-089 halla {}",
                        counters.signs_failed, evaluation.signs_failed
                    ),
                ));
            }
        }
    }
    if index != usize::MAX {
        ensure_all_bins(path, bin_position, &bins, &units[index])?;
    }
    walk.finish(F5OracleCause::LinkageDisequilibrium, path)?;
    Ok(cursor.rows)
}

fn ensure_all_bins(
    path: &'static str,
    published: usize,
    bins: &[&F5DeclaredLdBin],
    unit: &Unit,
) -> Checked<()> {
    if published != bins.len() {
        return Err(reject(
            F5OracleCause::LinkageDisequilibrium,
            path,
            "bins",
            format!(
                "la unidad {}×{} publica {published} bins y su perfil declara {}",
                unit.population_id,
                unit.block_id,
                bins.len()
            ),
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// `f5-founders.tsv` y las copias de `truth.jsonl`, en paso acompasado
// ---------------------------------------------------------------------------

/// Fila fundadora ya validada en su forma.
struct FounderRow {
    population_id: String,
    individual_id: String,
    sequence_id: String,
    copy_index: u64,
    founder_copy_id: String,
    haplotype_id: String,
    unit: usize,
    rank: u32,
}

/// Lector de `truth.jsonl` que solo entrega las copias, en su orden.
struct TruthCopies {
    reader: BufReader<File>,
    line: String,
    header: TruthHeader,
    closed_with_lf: bool,
    /// Primer registro de genotipo ya leído, que cierra la sección de copias.
    exhausted: bool,
}

impl TruthCopies {
    fn open(payload_root: &Path) -> Checked<Self> {
        let mut reader = open_artifact(payload_root, F5_TRUTH_PATH)?;
        let mut line = String::new();
        let read = reader
            .read_line(&mut line)
            .map_err(|error| truth_reject("encoding", error.to_string()))?;
        if read == 0 {
            return Err(truth_reject("header", "documento vacío"));
        }
        if !line.ends_with('\n') {
            return Err(truth_reject("eol", "falta el LF final"));
        }
        let header = read_truth_header(line.trim_end_matches('\n'))?;
        Ok(Self {
            reader,
            line: String::new(),
            header,
            closed_with_lf: true,
            exhausted: false,
        })
    }

    fn next_copy(&mut self) -> Checked<Option<TruthCopy>> {
        if self.exhausted {
            return Ok(None);
        }
        self.line.clear();
        let read = self
            .reader
            .read_line(&mut self.line)
            .map_err(|error| truth_reject("encoding", error.to_string()))?;
        if read == 0 {
            self.exhausted = true;
            // Sombreada: cuando esta comprobación se dispara, el recuento MAN-70
            // —que corre después sobre el mismo fichero— rechaza el mismo defecto
            // con el mismo campo y el mismo texto, de modo que desde
            // `audit_f5_package` no se puede distinguir cuál de las dos ha
            // hablado. La prueba la dispara (consta en el barrido) pero no puede
            // demostrar que sea ella la que rechaza; quitarla no cambia ningún
            // resultado observable.
            if !self.closed_with_lf {
                return Err(truth_reject("eol", "falta el LF final"));
            }
            return Ok(None);
        }
        self.closed_with_lf = self.line.ends_with('\n');
        let body = self.line.trim_end_matches('\n');
        if body.is_empty() {
            return Err(truth_reject("row", "línea vacía prohibida"));
        }
        let value: Value = serde_json::from_str(body)
            .map_err(|error| truth_reject("row", format!("registro ilegible: {error}")))?;
        match value.get("record_type").and_then(Value::as_str) {
            Some("haplotype_copy") => Ok(Some(read_truth_copy(&value, &self.header)?)),
            Some("genotype") => {
                self.exhausted = true;
                Ok(None)
            }
            _ => Err(truth_reject("row", "registro sin record_type reconocible")),
        }
    }
}

/// Recorre `f5-founders.tsv` y las copias de `truth.jsonl` a la vez.
///
/// Ambas llegan en la misma clave —población, individuo, secuencia y copia—,
/// de modo que un solo acompasamiento basta y nada proporcional a
/// `B*C_total` de SNP llega a retenerse.
#[allow(clippy::too_many_lines)]
fn walk_founders_and_copies(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    oracle: &IdOracle,
    units: &mut [Unit],
) -> Checked<(u64, u64, u64)> {
    let path = F5_TSV_CONTRACTS[FOUNDERS].path();
    let mut cursor = TsvCursor::open(open_artifact(payload_root, path)?, FOUNDERS)?;
    let mut copies = TruthCopies::open(payload_root)?;
    if copies.header.execution_id != declaration.execution_id
        || copies.header.reference_id != declaration.reference_id
        || copies.header.scenario_id != declaration.scenario_id
        || copies.header.reference_content_sha256 != declaration.reference_sha256
    {
        return Err(truth_reject(
            "header",
            "la cabecera v2 no declara la ejecución, referencia o escenario de la declaración",
        ));
    }
    let expected_digests: Vec<(String, String, String)> = vec![
        (
            "reference".to_owned(),
            declaration.reference_id.clone(),
            declaration.reference_sha256.clone(),
        ),
        (
            "matrix".to_owned(),
            declaration.profile_id.clone(),
            declaration.profile_metadata_sha256.clone(),
        ),
    ];
    if copies.header.source_digests != expected_digests {
        return Err(truth_reject(
            "source_digests",
            "las huellas de fuente no son las de la referencia y el perfil declarados",
        ));
    }
    for (position, subpopulation) in declaration.subpopulations.iter().enumerate() {
        let published = copies.header.populations.get(position).ok_or_else(|| {
            truth_reject("populations", "la cabecera declara menos subpoblaciones")
        })?;
        if published.population_id != subpopulation.population_id
            || published.individual_count != subpopulation.individual_count
            || published.ploidy != u64::from(declaration.ploidy)
        {
            return Err(truth_reject(
                "populations",
                format!(
                    "la subpoblación {} no publica el censo o la ploidía declarados",
                    subpopulation.population_id
                ),
            ));
        }
    }
    if copies.header.populations.len() != declaration.subpopulations.len() {
        return Err(truth_reject(
            "populations",
            "la cabecera publica más subpoblaciones que la declaración",
        ));
    }

    let mut previous: Option<(String, String, String, u64, String)> = None;
    let mut group: Vec<FounderRow> = Vec::new();
    let mut copy_count = 0_u64;
    // Individuos distintos: las filas llegan por población e individuo, de
    // modo que basta reconocer los cambios de grupo y no retener ningún
    // conjunto.
    let mut individuals = 0_u64;
    let mut current_individual: Option<(String, String)> = None;
    let mut pending = cursor.advance()?;
    while pending {
        let population_id = cursor.field(0).to_owned();
        let individual_id = cursor.field(1).to_owned();
        let sequence_id = cursor.field(2).to_owned();
        let copy_index = cursor.u64(F5OracleCause::Assignment, 3)?;
        let block_id = cursor.field(4).to_owned();
        if copy_index >= u64::from(declaration.ploidy) {
            return Err(reject(
                F5OracleCause::Assignment,
                path,
                "copy_index",
                "copy_index fuera de 0..ploidy-1",
            ));
        }
        let key = (
            population_id.clone(),
            individual_id.clone(),
            sequence_id.clone(),
            copy_index,
            block_id.clone(),
        );
        if previous.as_ref().is_some_and(|last| *last >= key) {
            return Err(reject(
                F5OracleCause::Assignment,
                path,
                "order",
                "el orden no es población, individuo, secuencia, copia y bloque",
            ));
        }
        previous = Some(key);
        let individual = (population_id.clone(), individual_id.clone());
        if current_individual.as_ref() != Some(&individual) {
            individuals += 1;
            current_individual = Some(individual);
        }
        let unit = unit_index(
            declaration,
            geometry,
            F5OracleCause::Assignment,
            path,
            &population_id,
            &block_id,
        )?;
        let block_index = units[unit].block_index;
        if geometry.blocks[block_index as usize].sequence_id != sequence_id {
            return Err(reject(
                F5OracleCause::Assignment,
                path,
                "sequence_id",
                "la copia no pertenece a la secuencia de su bloque",
            ));
        }
        let founder_copy_id = cursor.field(5).to_owned();
        if founder_copy_id
            != oracle.founder_copy(&population_id, &individual_id, &sequence_id, copy_index)?
        {
            return Err(reject(
                F5OracleCause::Identity,
                path,
                "founder_copy_id",
                "la fila no rederiva su FCP_",
            ));
        }
        let haplotype_id = cursor.field(7).to_owned();
        let rank = units[unit]
            .haplotype_ids
            .iter()
            .position(|candidate| *candidate == haplotype_id)
            .ok_or_else(|| {
                reject(
                    F5OracleCause::Assignment,
                    path,
                    "haplotype_id",
                    format!("{haplotype_id} no existe en f5-haplotypes.tsv"),
                )
            })?;
        if cursor.field(6)
            != oracle.assignment(&population_id, block_index, &founder_copy_id, &haplotype_id)?
        {
            return Err(reject(
                F5OracleCause::Identity,
                path,
                "block_assignment_id",
                "la fila no rederiva su FHA_",
            ));
        }
        let row = FounderRow {
            population_id,
            individual_id,
            sequence_id,
            copy_index,
            founder_copy_id,
            haplotype_id,
            unit,
            // Inalcanzable: `rank` es la posición del haplotipo dentro de los de su
            // unidad, y `f5-haplotypes.tsv` lleva una fila de texto por rango.
            rank: u32::try_from(rank).map_err(|_| {
                reject(
                    F5OracleCause::Assignment,
                    path,
                    "haplotype_id",
                    "el rango no cabe en u32",
                )
            })?,
        };
        let closes = group
            .first()
            .is_some_and(|first| first.founder_copy_id != row.founder_copy_id);
        if closes {
            close_copy_group(&mut copies, declaration, geometry, units, &group)?;
            copy_count += 1;
            group.clear();
        }
        group.push(row);
        pending = cursor.advance()?;
    }
    if !group.is_empty() {
        close_copy_group(&mut copies, declaration, geometry, units, &group)?;
        copy_count += 1;
    }
    if copies.next_copy()?.is_some() {
        return Err(truth_reject(
            "copy",
            "truth.jsonl publica una copia que f5-founders.tsv no asigna",
        ));
    }
    if copy_count != copies.header.copies {
        return Err(truth_reject(
            "counts",
            format!(
                "la cabecera declara {} copias y el paquete asigna {copy_count}",
                copies.header.copies
            ),
        ));
    }
    Ok((cursor.rows, copy_count, individuals))
}

/// Cierra una copia: contrasta el grupo de filas fundadoras con su registro
/// `haplotype_copy` y acumula cuota, alelos y heterocigosis.
fn close_copy_group(
    copies: &mut TruthCopies,
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    units: &mut [Unit],
    group: &[FounderRow],
) -> Checked<()> {
    let path = F5_TSV_CONTRACTS[FOUNDERS].path();
    let first = &group[0];
    let copy = copies.next_copy()?.ok_or_else(|| {
        truth_reject(
            "copy",
            "f5-founders.tsv asigna una copia que truth.jsonl no publica",
        )
    })?;
    if copy.identity.reference_id != declaration.reference_id
        || copy.identity.individual_id != first.individual_id
        || copy.identity.sequence_id != first.sequence_id
        || copy.identity.copy_index != first.copy_index
        || copy.population_id != first.population_id
        || copy.founder_copy_id != first.founder_copy_id
    {
        return Err(truth_reject(
            "copy.identity",
            format!(
                "la copia publicada no es la que f5-founders.tsv asigna a {}",
                first.founder_copy_id
            ),
        ));
    }
    // Los alelos que la copia debe llevar salen de sus bloques asignados, no
    // de la propia proyección: cada bloque aporta los sitios donde su fila
    // haplotípica lleva el alelo alternativo.
    let mut expected: Vec<(CopySnp, String)> = Vec::new();
    for row in group {
        let unit = &units[row.unit];
        let block_index = unit.block_index as usize;
        for (site_index, site) in geometry.sites[block_index].iter().enumerate() {
            if unit.bit(row.rank as usize, site_index) {
                expected.push((
                    CopySnp {
                        site_id: site.site_id.clone(),
                        position_1based: site.position_1based,
                        reference_base: site.reference_base,
                        alternate_base: site.alternate_base,
                    },
                    row.haplotype_id.clone(),
                ));
            }
        }
    }
    expected.sort_by_key(|(snp, _)| snp.position_1based);
    if expected.len() != copy.snps.len() {
        return Err(truth_reject(
            "copy.modification",
            format!(
                "la copia {} publica {} SNP y sus bloques asignan {}",
                copy.founder_copy_id,
                copy.snps.len(),
                expected.len()
            ),
        ));
    }
    for (position, (snp, haplotype_id)) in expected.iter().enumerate() {
        if copy.snps[position] != *snp {
            return Err(truth_reject(
                "copy.modification",
                format!(
                    "la copia {} publica un SNP que sus bloques no asignan",
                    copy.founder_copy_id
                ),
            ));
        }
        if copy.haplotype_ids[position] != *haplotype_id {
            return Err(reject(
                F5OracleCause::Assignment,
                F5_TRUTH_PATH,
                "copy.snp.provenance.matrix",
                format!(
                    "la procedencia del SNP {} nombra {} y f5-founders.tsv asigna {haplotype_id}",
                    snp.site_id, copy.haplotype_ids[position]
                ),
            ));
        }
    }
    // Cuota realizada, alelos por sitio y heterocigosis, unidad a unidad.
    for row in group {
        let ploidy = usize::from(declaration.ploidy);
        let unit = &mut units[row.unit];
        if unit.rank_counts.is_empty() {
            unit.rank_counts = vec![0; unit.haplotypes];
        }
        unit.rank_counts[row.rank as usize] += 1;
        unit.ranks.push(row.rank);
        for site in 0..unit.sites {
            if unit.bit(row.rank as usize, site) {
                unit.alt_copies[site] += 1;
            }
        }
        if unit.ranks.len().is_multiple_of(ploidy) {
            unit.individuals += 1;
            let seats = &unit.ranks[unit.ranks.len() - ploidy..];
            if seats.iter().any(|rank| *rank != seats[0]) {
                unit.heterozygous += 1;
            }
        }
    }
    // Defensa **inalcanzable**: las filas de una copia llegan por `block_id`
    // estrictamente creciente y cada una existe en el layout, de modo que no
    // pueden ser más que los bloques colocados. Se conserva porque acota el
    // vector `expected` sin depender de esa cadena.
    if group.len() as u64 > geometry.blocks.len() as u64 {
        return Err(reject(
            F5OracleCause::Assignment,
            path,
            "block_id",
            "una copia recibe más bloques de los que el layout coloca",
        ));
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Auditoría DTC-049
// ---------------------------------------------------------------------------

const NORMALIZATION_ORDER: [F5NormalizationSet; 4] = [
    F5NormalizationSet::BlockLengths,
    F5NormalizationSet::SnpGaps,
    F5NormalizationSet::SnpSpectrum,
    F5NormalizationSet::HaplotypeWeights,
];

/// Recalcula una auditoría DTC-049 a partir de sus propias cifras declaradas.
///
/// Se expone como función de escalares —y no de la estructura del perfil—
/// para que la batería pueda perturbar cada una de las seis por separado y
/// comprobar que ninguna sobrevive.
#[allow(clippy::too_many_arguments)]
fn check_normalization(
    set: F5NormalizationSet,
    algorithm: &str,
    tolerance: f64,
    input_keys: &[&str],
    input_values: &[f64],
    input_sum: f64,
    input_delta: f64,
    effective_keys: &[&str],
    effective_values: &[f64],
    output_sum: f64,
    output_delta: f64,
) -> Checked<()> {
    let path = "effective-config.json#/f5/normalizations";
    if algorithm != F5_PROFILE_NORMALIZATION_ALGORITHM
        || tolerance != F5_PROFILE_NORMALIZATION_TOLERANCE
    {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "algorithm",
            "el algoritmo o la tolerancia declarados no son los de DTC-049",
        ));
    }
    if input_values.is_empty() || input_keys.len() != input_values.len() {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "input_values",
            format!(
                "la auditoría {} no declara sus valores de entrada",
                set.as_str()
            ),
        ));
    }
    let recomputed_sum = neumaier_total(input_values.iter().copied());
    if recomputed_sum != input_sum {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "input_sum",
            format!(
                "{} declara input_sum={input_sum} y Neumaier da {recomputed_sum}",
                set.as_str()
            ),
        ));
    }
    if (input_sum - 1.0).abs() != input_delta {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "input_delta",
            format!("{} declara un delta de entrada ajeno", set.as_str()),
        ));
    }
    if !(input_sum.is_finite() && input_sum > 0.0 && input_delta <= tolerance) {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "input_delta",
            format!("{} no se normaliza bajo DTC-049", set.as_str()),
        ));
    }
    if effective_keys.len() != input_keys.len() || effective_values.len() != input_values.len() {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "effective_values",
            format!("{} pierde o inventa claves al normalizar", set.as_str()),
        ));
    }
    for position in 0..input_keys.len() {
        if input_keys[position] != effective_keys[position] {
            return Err(reject(
                F5OracleCause::Normalization,
                path,
                "canonical_order",
                format!(
                    "{} no conserva el orden canónico de sus claves",
                    set.as_str()
                ),
            ));
        }
        if effective_values[position] != input_values[position] / input_sum {
            return Err(reject(
                F5OracleCause::Normalization,
                path,
                "effective_values",
                format!(
                    "{}: la clave {} declara {} y la normalización da {}",
                    set.as_str(),
                    input_keys[position],
                    effective_values[position],
                    input_values[position] / input_sum
                ),
            ));
        }
    }
    let recomputed_output = neumaier_total(effective_values.iter().copied());
    if recomputed_output != output_sum {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "output_sum",
            format!(
                "{} declara output_sum={output_sum} y la normalización da {recomputed_output}",
                set.as_str()
            ),
        ));
    }
    if (output_sum - 1.0).abs() != output_delta || output_delta > tolerance {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "output_delta",
            format!(
                "{} declara un delta de salida ajeno o fuera de tolerancia",
                set.as_str()
            ),
        ));
    }
    Ok(())
}

/// Recalcula las cuatro auditorías DTC-049 declaradas.
///
/// No cree ni una de sus seis cifras: rehace la suma de Neumaier de las
/// entradas, el delta, la división por la suma, la suma de salida y su delta,
/// y exige el orden canónico, el algoritmo y la tolerancia literales.
fn audit_normalizations(declaration: &F5PackageDeclaration) -> Checked<()> {
    let path = "effective-config.json#/f5/normalizations";
    if declaration.normalizations.len() != NORMALIZATION_ORDER.len() {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "normalizations",
            format!(
                "se declaran {} auditorías y DTC-089 exige cuatro",
                declaration.normalizations.len()
            ),
        ));
    }
    for (audit, expected) in declaration.normalizations.iter().zip(NORMALIZATION_ORDER) {
        if audit.set() != expected {
            return Err(reject(
                F5OracleCause::Normalization,
                path,
                "set_kind",
                format!(
                    "la auditoría {} aparece donde DTC-089 pone {}",
                    audit.set().as_str(),
                    expected.as_str()
                ),
            ));
        }
        let input_keys: Vec<&str> = audit
            .input_values()
            .iter()
            .map(F5NamedProbability::key)
            .collect();
        let input_values: Vec<f64> = audit
            .input_values()
            .iter()
            .map(F5NamedProbability::value)
            .collect();
        let effective_keys: Vec<&str> = audit
            .effective_values()
            .iter()
            .map(F5NamedProbability::key)
            .collect();
        let effective_values: Vec<f64> = audit
            .effective_values()
            .iter()
            .map(F5NamedProbability::value)
            .collect();
        check_normalization(
            audit.set(),
            audit.algorithm(),
            audit.tolerance(),
            &input_keys,
            &input_values,
            audit.input_sum(),
            audit.input_delta(),
            &effective_keys,
            &effective_values,
            audit.output_sum(),
            audit.output_delta(),
        )?;
    }
    let weights: Vec<f64> = declaration.normalizations[3]
        .effective_values()
        .iter()
        .map(F5NamedProbability::value)
        .collect();
    if weights != declaration.haplotype_weights {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "haplotype_weights",
            "los pesos efectivos declarados no son los de la auditoría DTC-049",
        ));
    }
    // Las bandas del espectro se atan a la misma auditoría en todo lo que el
    // efectivo publica: su clave es `alt_af_min` y su valor efectivo el peso.
    // `alt_af_max` no aparece en ningún documento y queda declarada como
    // frontera del perfil, atado por `profile_metadata_sha256`.
    let spectrum = &declaration.normalizations[2];
    if declaration.snp_spectrum.len() != spectrum.effective_values().len() {
        return Err(reject(
            F5OracleCause::Normalization,
            path,
            "snp_spectrum",
            format!(
                "se declaran {} bandas de espectro y la auditoría DTC-049 trae {}",
                declaration.snp_spectrum.len(),
                spectrum.effective_values().len()
            ),
        ));
    }
    let mut previous: Option<f64> = None;
    for (band, effective) in declaration
        .snp_spectrum
        .iter()
        .zip(spectrum.effective_values())
    {
        // Inalcanzable: `F5NormalizationAudit` no tiene constructor público. La
        // única forma de obtener una es el cargador del perfil, que escribe cada
        // clave del espectro como el texto de un `f64` que él mismo leyó, de modo
        // que una prueba no puede fabricar una auditoría con una clave ilegible.
        let key: f64 = effective.key().parse().map_err(|_| {
            reject(
                F5OracleCause::Normalization,
                path,
                "snp_spectrum",
                format!("la clave {:?} no es un alt_af_min", effective.key()),
            )
        })?;
        if band.alt_af_min != key || band.weight != effective.value() {
            return Err(reject(
                F5OracleCause::Normalization,
                path,
                "snp_spectrum",
                format!(
                    "la banda [{}, {}] con peso {} no es la fila {:?}={} de la auditoría",
                    band.alt_af_min,
                    band.alt_af_max,
                    band.weight,
                    effective.key(),
                    effective.value()
                ),
            ));
        }
        if !(0.0..=1.0).contains(&band.alt_af_min)
            || !(0.0..=1.0).contains(&band.alt_af_max)
            || band.alt_af_min > band.alt_af_max
            || previous.is_some_and(|last| band.alt_af_min <= last)
        {
            return Err(reject(
                F5OracleCause::Normalization,
                path,
                "snp_spectrum",
                "las bandas del espectro no son intervalos disjuntos y crecientes de [0,1]",
            ));
        }
        previous = Some(band.alt_af_max);
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Cuota, modo, tolerancia y F
// ---------------------------------------------------------------------------

fn declared_sampling<'a>(
    declaration: &'a F5PackageDeclaration,
    unit: &Unit,
) -> Checked<&'a F5DeclaredSampling> {
    let mut found = None;
    for row in &declaration.sampling {
        if row.population_id == unit.population_id && row.block_id == unit.block_id {
            if found.is_some() {
                return Err(reject(
                    F5OracleCause::Sampling,
                    "effective-config.json#/f5/sampling",
                    "sampling",
                    format!(
                        "la unidad {}×{} aparece dos veces",
                        unit.population_id, unit.block_id
                    ),
                ));
            }
            found = Some(row);
        }
    }
    found.ok_or_else(|| {
        reject(
            F5OracleCause::Sampling,
            "effective-config.json#/f5/sampling",
            "sampling",
            format!(
                "la unidad {}×{} no declara modo de muestreo",
                unit.population_id, unit.block_id
            ),
        )
    })
}

/// Modo, cuota y tolerancia de una unidad, contrastados con lo que sus copias
/// realizaron.
fn audit_sampling(
    declaration: &F5PackageDeclaration,
    unit: &Unit,
    inbreeding: Option<f64>,
) -> Checked<F5DeclaredMode> {
    let path = "effective-config.json#/f5/sampling";
    let row = declared_sampling(declaration, unit)?;
    if row.requested_mode != row.effective_mode {
        return Err(reject(
            F5OracleCause::Sampling,
            path,
            "effective_mode",
            "DTC-089 prohíbe degradar una cuota: el modo efectivo debe ser el solicitado",
        ));
    }
    let copies = unit.individuals * u64::from(declaration.ploidy);
    match row.effective_mode {
        F5DeclaredMode::Probabilistic => {
            if row.controlled_max_abs_error_declared.is_some()
                || row.controlled_max_abs_error_effective.is_some()
            {
                return Err(reject(
                    F5OracleCause::Sampling,
                    path,
                    "controlled_max_abs_error_effective",
                    "una unidad probabilística no puede declarar tolerancia controlada",
                ));
            }
        }
        F5DeclaredMode::Controlled => {
            let effective = row.controlled_max_abs_error_effective.ok_or_else(|| {
                reject(
                    F5OracleCause::Sampling,
                    path,
                    "controlled_max_abs_error_effective",
                    "una unidad controlada debe declarar su tolerancia efectiva",
                )
            })?;
            if copies == 0 {
                return Err(reject(
                    F5OracleCause::Sampling,
                    path,
                    "copies",
                    "una unidad controlada sin copias no tiene cuota",
                ));
            }
            let ceiling = 0.5 / copies as f64;
            if !effective.is_finite() || !(0.0..=ceiling).contains(&effective) {
                return Err(reject(
                    F5OracleCause::Sampling,
                    path,
                    "controlled_max_abs_error_effective",
                    format!("la tolerancia efectiva {effective} no pertenece a [0,0.5/C]"),
                ));
            }
            if let Some(declared) = row.controlled_max_abs_error_declared {
                if declared != effective {
                    return Err(reject(
                        F5OracleCause::Sampling,
                        path,
                        "controlled_max_abs_error_declared",
                        "la tolerancia declarada y la efectiva difieren",
                    ));
                }
            } else if effective != ceiling {
                return Err(reject(
                    F5OracleCause::Sampling,
                    path,
                    "controlled_max_abs_error_effective",
                    "sin tolerancia declarada la efectiva debe ser 0.5/C",
                ));
            }
            if inbreeding.is_some_and(|value| value != 0.0) {
                return Err(reject(
                    F5OracleCause::Inbreeding,
                    path,
                    "inbreeding_f",
                    format!(
                        "la unidad controlada {}×{} exige F efectivo igual a cero",
                        unit.population_id, unit.block_id
                    ),
                ));
            }
            // La cuota realizada: cada rango recibe `n_h` copias y la
            // desviación por rango cabe en la tolerancia efectiva.
            //
            // El reparto total es una defensa **inalcanzable**: `rank_counts`
            // se incrementa una vez por asiento recorrido y el censo ya exigió
            // `ranks.len() == individual_count*ploidy`. Se conserva porque es
            // la invariante que la desviación por rango presupone.
            let total: u64 = unit.rank_counts.iter().sum();
            if total != copies {
                return Err(reject(
                    F5OracleCause::Sampling,
                    path,
                    "quota",
                    format!("la cuota reparte {total} copias y la unidad tiene {copies}"),
                ));
            }
            let tolerance_copies = effective * copies as f64;
            for (rank, count) in unit.rank_counts.iter().enumerate() {
                // Inalcanzable: `rank_counts` se crea con `unit.haplotypes` entradas y
                // `read_haplotypes` fija `unit.weights` con el mismo número, de modo
                // que todo rango tiene peso.
                let target = unit.weights.get(rank).copied().ok_or_else(|| {
                    reject(
                        F5OracleCause::Sampling,
                        path,
                        "quota",
                        "la cuota nombra un rango que la unidad no publica",
                    )
                })? * copies as f64;
                if (*count as f64 - target).abs() > tolerance_copies {
                    return Err(reject(
                        F5OracleCause::Sampling,
                        path,
                        "quota",
                        format!(
                            "el rango {rank} de {}×{} recibe {count} copias y su cuota es {target} con tolerancia {tolerance_copies}",
                            unit.population_id, unit.block_id
                        ),
                    ));
                }
            }
        }
    }
    Ok(row.effective_mode)
}

/// `F_realized = 1 - H_observed/H_expected`, reconstruido desde los rangos que
/// `f5-founders.tsv` publica.
///
/// **Frontera declarada:** esto es una *transcripción* de la fórmula, no una
/// rederivación independiente de ella. La suma de homocigosis se escribe en el
/// mismo orden de rango y con la misma aritmética binary64 no compensada que
/// DTC-089 fija, porque el valor publicado es bit a bit el de ese orden y
/// cualquier otra asociatividad daría otro último bit. Lo que sí es
/// independiente es su **entrada**: los rangos salen de `f5-founders.tsv` y de
/// `f5-haplotypes.tsv`, no del productor. Un cambio de la fórmula en DTC-089
/// obligaría a cambiar las dos copias, y por eso se dice aquí.
fn recompute_realized_f(unit: &Unit) -> Option<f64> {
    if unit.individuals == 0 {
        return None;
    }
    let copies = unit.ranks.len() as f64;
    let h_observed = unit.heterozygous as f64 / unit.individuals as f64;
    let mut homozygosity = 0.0_f64;
    for count in &unit.rank_counts {
        let frequency = *count as f64 / copies;
        homozygosity += frequency * frequency;
    }
    let h_expected = 1.0 - homozygosity;
    if h_expected.is_nan() || h_expected <= 0.0 {
        return None;
    }
    let value = 1.0 - h_observed / h_expected;
    value.is_finite().then_some(value)
}

fn audit_inbreeding(declaration: &F5PackageDeclaration, units: &[Unit]) -> Checked<()> {
    let path = "manifest.json#/f5/f_realized";
    for subpopulation in &declaration.subpopulations {
        match (declaration.ploidy, subpopulation.inbreeding_f) {
            (2, None) => {
                return Err(reject(
                    F5OracleCause::Inbreeding,
                    path,
                    "inbreeding_f",
                    format!(
                        "la subpoblación diploide {} no declara su F efectivo",
                        subpopulation.population_id
                    ),
                ));
            }
            // `!value.is_finite()` es redundante: `contains` ya es falso para NaN y
            // para los infinitos, de modo que ninguna prueba puede aislarla.
            (2, Some(value)) if !value.is_finite() || !(0.0..=1.0).contains(&value) => {
                return Err(reject(
                    F5OracleCause::Inbreeding,
                    path,
                    "inbreeding_f",
                    format!("F={value} no pertenece a [0,1]"),
                ));
            }
            (ploidy, Some(_)) if ploidy != 2 => {
                return Err(reject(
                    F5OracleCause::Inbreeding,
                    path,
                    "inbreeding_f",
                    format!("F solo se admite con ploidía 2 y la efectiva es {ploidy}"),
                ));
            }
            _ => {}
        }
    }
    if declaration.f_realized.len() != units.len() {
        return Err(reject(
            F5OracleCause::Inbreeding,
            path,
            "f_realized",
            format!(
                "se declaran {} filas de F realizado y hay {} unidades",
                declaration.f_realized.len(),
                units.len()
            ),
        ));
    }
    for (row, unit) in declaration.f_realized.iter().zip(units) {
        if row.population_id != unit.population_id || row.block_id != unit.block_id {
            return Err(reject(
                F5OracleCause::Inbreeding,
                path,
                "f_realized",
                "las filas de F realizado no llegan por subpoblación y bloque",
            ));
        }
        let declared_f = declaration
            .subpopulations
            .iter()
            .find(|subpopulation| subpopulation.population_id == unit.population_id)
            .and_then(|subpopulation| subpopulation.inbreeding_f);
        if row.target_f != declared_f {
            return Err(reject(
                F5OracleCause::Inbreeding,
                path,
                "target_f",
                format!(
                    "la unidad {}×{} publica target_f={:?} y su subpoblación declara {declared_f:?}",
                    unit.population_id, unit.block_id, row.target_f
                ),
            ));
        }
        let expected = if declaration.ploidy == 2 {
            recompute_realized_f(unit)
        } else {
            None
        };
        let (status, value) = match expected {
            Some(value) => ("evaluated", Some(value)),
            None => ("not_applicable", None),
        };
        if row.status != status || row.value != value {
            return Err(reject(
                F5OracleCause::Inbreeding,
                path,
                "value",
                format!(
                    "la unidad {}×{} publica {}/{:?} y la reconstrucción da {status}/{value:?}",
                    unit.population_id, unit.block_id, row.status, row.value
                ),
            ));
        }
    }
    Ok(())
}

// ---------------------------------------------------------------------------
// Vectores MAN-70 por dosis
// ---------------------------------------------------------------------------

/// `C(n, r)` en `f64`, escrito como producto para no desbordar con `C` grande
/// y `r` pequeño, que es el único caso que la ploidía produce.
fn choose(n: u64, r: u64) -> f64 {
    if r > n {
        return 0.0;
    }
    let mut value = 1.0_f64;
    for step in 0..r {
        value = value * (n - step) as f64 / (step + 1) as f64;
    }
    value
}

/// HWE sin condicionar: `P` copias independientes con probabilidad `p`.
fn binomial_probabilities(ploidy: u8, p: f64) -> Vec<f64> {
    (0..=u64::from(ploidy))
        .map(|dosage| {
            choose(u64::from(ploidy), dosage)
                * p.powi(dosage as i32)
                * (1.0 - p).powi((u64::from(ploidy) - dosage) as i32)
        })
        .collect()
}

/// F diploide: con probabilidad `F` las dos copias son idénticas por
/// descendencia y con `1-F` se muestrean por separado.
fn diploid_f_probabilities(p: f64, f: f64) -> Vec<f64> {
    let heterozygous = 2.0 * p * (1.0 - p) * (1.0 - f);
    vec![
        (1.0 - p) * (1.0 - p) + f * p * (1.0 - p),
        heterozygous,
        p * p + f * p * (1.0 - p),
    ]
}

/// Hipergeométrica inducida por permutar sin reemplazo las `C` copias con el
/// recuento ALT fijado `K`, que es lo que el modo controlado realiza.
fn hypergeometric_probabilities(copies: u64, alt: u64, ploidy: u8) -> Vec<f64> {
    let denominator = choose(copies, u64::from(ploidy));
    (0..=u64::from(ploidy))
        .map(|dosage| {
            if denominator == 0.0 {
                return 0.0;
            }
            choose(alt, dosage) * choose(copies - alt, u64::from(ploidy) - dosage) / denominator
        })
        .collect()
}

/// Individuos por dosis alternativa, reconstruidos desde los rangos asignados.
fn realized_counts(unit: &Unit, site: usize, ploidy: usize) -> Vec<u64> {
    let mut counts = vec![0_u64; ploidy + 1];
    for seats in unit.ranks.chunks(ploidy) {
        let dosage = seats
            .iter()
            .filter(|rank| unit.bit(**rank as usize, site))
            .count();
        counts[dosage] += 1;
    }
    counts
}

// ---------------------------------------------------------------------------
// Auditoría completa
// ---------------------------------------------------------------------------

fn root_seed_bytes(seed: &str) -> Checked<[u8; 32]> {
    if !is_lowercase_sha256(seed) {
        return Err(reject(
            F5OracleCause::Seed,
            "effective-config.json#/f5/seed",
            "seed",
            "la semilla declarada no es hex minúscula de 64",
        ));
    }
    // Inalcanzable: `is_lowercase_sha256` acaba de exigir 64 dígitos hexadecimales
    // minúsculos, que `decode_hex` siempre decodifica en 32 bytes.
    decode_hex(seed)
        .and_then(|bytes| bytes.try_into().ok())
        .ok_or_else(|| {
            reject(
                F5OracleCause::Seed,
                "effective-config.json#/f5/seed",
                "seed",
                "la semilla declarada no es una raíz de 256 bits",
            )
        })
}

/// Contrasta los ocho enlaces declarados con los bytes que hay en el destino.
fn audit_links(payload_root: &Path, declaration: &F5PackageDeclaration) -> Checked<u64> {
    if declaration.artifacts.len() != F5_TSV_CONTRACTS.len() {
        return Err(reject(
            F5OracleCause::Link,
            "manifest.json#/f5/artifacts",
            "artifacts",
            format!(
                "se declaran {} enlaces y un éxito publica siete",
                declaration.artifacts.len()
            ),
        ));
    }
    let mut output_bytes = 0_u64;
    for (link, contract) in declaration.artifacts.iter().zip(F5_TSV_CONTRACTS) {
        if link.artifact_id != contract.artifact_id()
            || link.artifact_kind != contract.artifact_kind()
            || link.path != contract.path()
            || link.schema_version != F5_TSV_SCHEMA_VERSION
        {
            return Err(reject(
                F5OracleCause::Link,
                contract.path(),
                "artifact_id",
                format!(
                    "el enlace declara {}/{}/{}/{} y el contrato fija {}/{}/{}/{}",
                    link.artifact_id,
                    link.artifact_kind,
                    link.schema_version,
                    link.path,
                    contract.artifact_id(),
                    contract.artifact_kind(),
                    F5_TSV_SCHEMA_VERSION,
                    contract.path()
                ),
            ));
        }
        let (size, digest) = measure_artifact(payload_root, contract.path())?;
        if link.size_bytes != size || link.sha256 != digest {
            return Err(reject(
                F5OracleCause::Link,
                contract.path(),
                "sha256",
                format!(
                    "el enlace declara {} bytes y {} y los bytes miden {size} y {digest}",
                    link.size_bytes, link.sha256
                ),
            ));
        }
        output_bytes += size;
    }
    let truth = &declaration.truth;
    if truth.artifact_id != F5_TRUTH_ARTIFACT_ID
        || truth.artifact_kind != F5_TRUTH_ARTIFACT_KIND
        || truth.schema_version != F5_TRUTH_SCHEMA_VERSION
        || truth.projection_version != F5_TRUTH_PROJECTION_VERSION
        || truth.path != F5_TRUTH_PATH
    {
        return Err(reject(
            F5OracleCause::Link,
            F5_TRUTH_PATH,
            "artifact_id",
            "el enlace de la proyección fundadora no es el literal de DTC-089",
        ));
    }
    let (size, digest) = measure_artifact(payload_root, F5_TRUTH_PATH)?;
    if truth.size_bytes != size || truth.sha256 != digest {
        return Err(reject(
            F5OracleCause::Link,
            F5_TRUTH_PATH,
            "sha256",
            format!(
                "el enlace declara {} bytes y {} y los bytes miden {size} y {digest}",
                truth.size_bytes, truth.sha256
            ),
        ));
    }
    Ok(output_bytes + size)
}

/// Valida el paquete publicado contra su declaración y devuelve los resúmenes
/// MAN-70 que el manifiesto de F5.4.3 copiará.
///
/// `payload_root` es la raíz bajo la que viven `outputs/...`. La función no
/// recibe ninguna estructura del productor: todo lo que juzga sale de los
/// bytes o de la declaración que se juzga con ellos.
///
/// # Errors
///
/// [`F5OracleRejection`], cuya [`F5OracleCause`] nombra el contrato incumplido.
#[allow(clippy::too_many_lines)]
pub fn audit_f5_package(
    payload_root: &Path,
    declaration: &F5PackageDeclaration,
) -> Checked<F5OracleReport> {
    // El catálogo RNG no es una preferencia del publicador: v2 exactamente
    // cuando hay distribución declarada, v1 exactamente cuando no la hay. Un
    // paquete que sortease recuentos y publicase v1 —o al revés— describiría
    // mal el espacio de decisiones del que salen sus bytes.
    let expected_rng_catalog = if declaration.haplotype_count_distribution.is_some() {
        F5_ORACLE_RNG_CATALOG_V2
    } else {
        F5_ORACLE_RNG_CATALOG
    };
    if declaration.id_catalog_id != F5_ORACLE_ID_CATALOG
        || declaration.rng_catalog_id != expected_rng_catalog
    {
        return Err(reject(
            F5OracleCause::Randomness,
            "effective-config.json#/f5",
            "id_catalog,rng_catalog",
            format!(
                "se declaran {}/{} y los catálogos vigentes son {F5_ORACLE_ID_CATALOG}/{expected_rng_catalog}",
                declaration.id_catalog_id, declaration.rng_catalog_id
            ),
        ));
    }
    if declaration.ploidy == 0 {
        return Err(reject(
            F5OracleCause::Assignment,
            "effective-config.json#/f5/ploidy",
            "ploidy",
            "la ploidía efectiva es cero",
        ));
    }
    if declaration.subpopulations.is_empty() {
        return Err(reject(
            F5OracleCause::Assignment,
            "effective-config.json#/f5/subpopulations",
            "subpopulations",
            "no se declara ninguna subpoblación",
        ));
    }
    if declaration
        .subpopulations
        .windows(2)
        .any(|pair| pair[0].population_id >= pair[1].population_id)
    {
        return Err(reject(
            F5OracleCause::Assignment,
            "effective-config.json#/f5/subpopulations",
            "subpopulations",
            "las subpoblaciones no llegan ordenadas por population_id",
        ));
    }
    let root = root_seed_bytes(&declaration.seed)?;
    let ids = IdOracle::new(
        root,
        &declaration.execution_id,
        &declaration.scenario_id,
        &declaration.reference_id,
    )?;
    let subseeds = SubseedOracle::new(root, &declaration.scenario_id)?;
    audit_normalizations(declaration)?;
    let output_bytes = audit_links(payload_root, declaration)?;

    let geometry = read_layout(payload_root, &ids, &subseeds, &declaration.bounds)?;
    let by_id = sites_by_id(&geometry);
    let mut units = build_units(declaration, &geometry);
    let attempt_rows = read_attempts(payload_root, declaration, &geometry, &mut units)?;
    let haplotype_rows = read_haplotypes(
        payload_root,
        declaration,
        &geometry,
        &ids,
        &subseeds,
        &mut units,
    )?;
    // Los contadores del mejor intento que solo dependen de la matriz se
    // comprueban aquí; los que dependen de los bins, al cerrar cada unidad LD.
    for unit in &units {
        // Inalcanzable: `read_attempts` ya exigió que toda unidad tenga un único
        // mejor intento.
        let counters = unit.best_counters.ok_or_else(|| {
            reject(
                F5OracleCause::Attempt,
                F5_TSV_CONTRACTS[ATTEMPTS].path(),
                "selected_best",
                format!(
                    "la unidad {}×{} no publica ningún mejor intento",
                    unit.population_id, unit.block_id
                ),
            )
        })?;
        let mut rows: Vec<Vec<bool>> = Vec::with_capacity(unit.haplotypes);
        for haplotype in 0..unit.haplotypes {
            rows.push(
                (0..unit.sites)
                    .map(|site| unit.bit(haplotype, site))
                    .collect(),
            );
        }
        rows.sort();
        let duplicates = rows.windows(2).filter(|pair| pair[0] == pair[1]).count() as u64;
        if counters.duplicate_rows != duplicates {
            return Err(reject(
                F5OracleCause::Attempt,
                F5_TSV_CONTRACTS[ATTEMPTS].path(),
                "duplicate_rows",
                format!(
                    "el mejor intento de {}×{} publica {} filas repetidas y la matriz trae {duplicates}",
                    unit.population_id, unit.block_id, counters.duplicate_rows
                ),
            ));
        }
    }
    let frequency_rows = read_frequencies(
        payload_root,
        declaration,
        &geometry,
        &by_id,
        &subseeds,
        &mut units,
    )?;
    // `max_marginal_error` del mejor intento es el máximo de los errores que
    // `f5-frequencies.tsv` acaba de publicar.
    for unit in &units {
        // Inalcanzable, por la misma razón que el bucle anterior.
        let counters = unit.best_counters.ok_or_else(|| {
            reject(
                F5OracleCause::Attempt,
                F5_TSV_CONTRACTS[ATTEMPTS].path(),
                "selected_best",
                "sin mejor intento",
            )
        })?;
        let mut observed = 0.0_f64;
        for site in 0..unit.sites {
            let error = (unit.realized_alt_af[site] - unit.target_alt_af[site]).abs();
            if error > observed {
                observed = error;
            }
        }
        if counters.max_marginal_error != observed {
            return Err(reject(
                F5OracleCause::Attempt,
                F5_TSV_CONTRACTS[ATTEMPTS].path(),
                "max_marginal_error",
                format!(
                    "el mejor intento de {}×{} publica {} y los marginales dan {observed}",
                    unit.population_id, unit.block_id, counters.max_marginal_error
                ),
            ));
        }
    }
    let ld_rows = read_ld_evaluation(payload_root, declaration, &geometry, &units)?;
    let (founder_rows, copies, founder_individuals) =
        walk_founders_and_copies(payload_root, declaration, &geometry, &ids, &mut units)?;

    // Censo y cuota por unidad, ya con los asientos recorridos.
    for unit in &units {
        // Inalcanzable: `units` se construye desde `declaration.subpopulations` y
        // cada unidad lleva el `population_id` de una de ellas.
        let census = declaration
            .subpopulations
            .iter()
            .find(|subpopulation| subpopulation.population_id == unit.population_id)
            .ok_or_else(|| {
                reject(
                    F5OracleCause::Assignment,
                    F5_TSV_CONTRACTS[FOUNDERS].path(),
                    "population_id",
                    "una unidad nombra una subpoblación sin censo",
                )
            })?;
        if unit.individuals != census.individual_count
            || unit.ranks.len() as u64 != census.individual_count * u64::from(declaration.ploidy)
        {
            return Err(reject(
                F5OracleCause::Assignment,
                F5_TSV_CONTRACTS[FOUNDERS].path(),
                "row",
                format!(
                    "la unidad {}×{} recibe {} individuos y {} copias, y el censo declara {} con ploidía {}",
                    unit.population_id,
                    unit.block_id,
                    unit.individuals,
                    unit.ranks.len(),
                    census.individual_count,
                    declaration.ploidy
                ),
            ));
        }
    }
    if declaration.sampling.len() != units.len() {
        return Err(reject(
            F5OracleCause::Sampling,
            "effective-config.json#/f5/sampling",
            "sampling",
            format!(
                "se declaran {} filas de muestreo y hay {} unidades",
                declaration.sampling.len(),
                units.len()
            ),
        ));
    }
    let mut modes = Vec::with_capacity(units.len());
    for unit in &units {
        let inbreeding = declaration
            .subpopulations
            .iter()
            .find(|subpopulation| subpopulation.population_id == unit.population_id)
            .and_then(|subpopulation| subpopulation.inbreeding_f);
        modes.push(audit_sampling(declaration, unit, inbreeding)?);
    }
    audit_inbreeding(declaration, &units)?;

    // El recuento MAN-70 abre `truth.jsonl` **por su cuenta** y no recibe nada
    // de lo anterior.
    let recount = recount_f5_man70(&payload_root.join(F5_TRUTH_PATH))?;
    // Inalcanzable: `recount.copies` es el número de líneas `haplotype_copy` que
    // `recount_f5_man70` acaba de igualar a `counts.copies` de la cabecera, y
    // `walk_founders_and_copies` igualó `copies` a ese mismo campo de esa misma
    // cabecera: los dos son iguales a él. Se conserva porque el recuento es una
    // frontera pública propia.
    if recount.copies != copies {
        return Err(reject(
            F5OracleCause::Man70,
            F5_TRUTH_PATH,
            "counts",
            format!(
                "el recuento independiente halla {} copias y el cruce con los TSV {copies}",
                recount.copies
            ),
        ));
    }
    // Defensa **inalcanzable**: la ploidía del recuento sale de la cabecera v2,
    // que `walk_founders_and_copies` ya contrastó contra la declarada. Se
    // conserva porque el recuento es una frontera pública propia y nada obliga
    // a que sus dos llamadores lleguen en ese orden.
    if recount.ploidy != declaration.ploidy {
        return Err(reject(
            F5OracleCause::Man70,
            F5_TRUTH_PATH,
            "ploidy",
            "la ploidía recontada no es la declarada",
        ));
    }
    let summaries = build_summaries(declaration, &geometry, &by_id, &units, &modes, &recount)?;

    let tsv_rows = geometry.layout_rows
        + geometry.site_rows
        + haplotype_rows
        + frequency_rows
        + ld_rows
        + founder_rows
        + attempt_rows;
    if tsv_rows > declaration.bounds.tsv_rows_bound {
        return Err(resource_reject(
            "outputs",
            format!(
                "el paquete publica {tsv_rows} filas TSV y la admisión selló {}",
                declaration.bounds.tsv_rows_bound
            ),
        ));
    }
    // `json_rows_bound` acota las **filas JSON del efectivo y del manifiesto**
    // —fundadores, subpoblaciones, muestreo, F realizado, rangos de empate,
    // resúmenes MAN-70, filas de normalización y nueve objetos fijos—, no las
    // líneas de `truth.jsonl`, que viven en `output_bytes` a través de
    // `truth_bound`. La confusión es fácil y el oráculo la declara: con censo 4
    // y ocho bloques, `truth.jsonl` publica 105 líneas y el bound vale 104,
    // porque su término de locus es `B*N_pop*S_bound` y no
    // `B*N_pop*S_bound*individual_count`.
    //
    // La cifra es una **cota inferior observable**: los `2*B*N_pop*H` rangos
    // de `quota_ties` no se pueden leer desde los bytes publicados —ningún TSV
    // los enumera— y por eso no se suman. Comparar una cota inferior contra el
    // bound sigue siendo correcto; lo que no puede hacerse es leerla como el
    // recuento exacto del manifiesto, y por eso el nombre lo dice.
    let observable_json_rows = founder_individuals
        + declaration.subpopulations.len() as u64
        + declaration.sampling.len() as u64
        + declaration.f_realized.len() as u64
        + summaries.len() as u64
        + 2 * declaration
            .normalizations
            .iter()
            .map(|audit| audit.input_values().len() as u64)
            .sum::<u64>()
        + F5_MANIFEST_FIXED_OBJECTS;
    if observable_json_rows > declaration.bounds.json_rows_bound {
        return Err(resource_reject(
            "manifest.json#/f5",
            format!(
                "el manifiesto publicaría al menos {observable_json_rows} filas JSON y la admisión selló {}",
                declaration.bounds.json_rows_bound
            ),
        ));
    }
    let truth_records = 1 + recount.copies + recount.genotypes;
    if output_bytes > declaration.bounds.output_bytes
        || declaration.bounds.output_bytes > declaration.bounds.output_limit_bytes
    {
        return Err(resource_reject(
            "outputs",
            format!(
                "el paquete ocupa {output_bytes} bytes, la admisión reservó {} y el límite es {}",
                declaration.bounds.output_bytes, declaration.bounds.output_limit_bytes
            ),
        ));
    }

    Ok(F5OracleReport {
        initialization_summaries: summaries,
        resources: F5ResourceAudit {
            tsv_rows,
            observable_json_rows,
            truth_records,
            output_bytes,
        },
        copies,
        genotypes: recount.genotypes,
    })
}

/// Construye los resúmenes MAN-70 cruzando lo realizado desde los TSV con lo
/// recontado desde la verdad, y exige la igualdad que REQ-MAN-70 obliga.
fn build_summaries(
    declaration: &F5PackageDeclaration,
    geometry: &Geometry,
    by_id: &[Vec<usize>],
    units: &[Unit],
    modes: &[F5DeclaredMode],
    recount: &F5Man70Recount,
) -> Checked<Vec<F5InitializationSummary>> {
    let ploidy = usize::from(declaration.ploidy);
    let mut summaries = Vec::new();
    let mut cursor = 0_usize;
    for (position, unit) in units.iter().enumerate() {
        // Inalcanzable, por la misma razón que el censo de `audit_f5_package`: las
        // unidades salen de `declaration.subpopulations`.
        let census = declaration
            .subpopulations
            .iter()
            .find(|subpopulation| subpopulation.population_id == unit.population_id)
            .ok_or_else(|| {
                reject(
                    F5OracleCause::Man70,
                    F5_TRUTH_PATH,
                    "population_id",
                    "una unidad nombra una subpoblación sin censo",
                )
            })?;
        let inbreeding = census.inbreeding_f.unwrap_or(0.0);
        let effective_mode = modes[position];
        let sampling = declared_sampling(declaration, unit)?;
        for site_index in by_id[unit.block_index as usize].iter().copied() {
            let site = &geometry.sites[unit.block_index as usize][site_index];
            let locus = recount.loci.get(cursor).ok_or_else(|| {
                reject(
                    F5OracleCause::Man70,
                    F5_TRUTH_PATH,
                    "recount",
                    format!(
                        "la verdad no reconta el locus {}×{}",
                        site.site_id, unit.population_id
                    ),
                )
            })?;
            cursor += 1;
            if locus.population_id != unit.population_id
                || locus.block_id != unit.block_id
                || locus.site_id != site.site_id
                || locus.sequence_id != geometry.blocks[unit.block_index as usize].sequence_id
                || locus.position_1based != site.position_1based
                || locus.reference_base != site.reference_base
                || locus.alternate_base != site.alternate_base
                || locus.individual_count != census.individual_count
            {
                return Err(reject(
                    F5OracleCause::Man70,
                    F5_TRUTH_PATH,
                    "recount",
                    format!(
                        "el locus recontado {}×{} no es el de f5-sites.tsv",
                        locus.site_id, locus.population_id
                    ),
                ));
            }
            // La igualdad que REQ-MAN-70 obliga. Es **inalcanzable**: los
            // portadores de cada genotipo son las `ploidy` copias del individuo
            // en esa secuencia, cada copia lleva exactamente los sitios que
            // `close_copy_group` la obligó a llevar según los TSV, y el locus ya
            // coincide en sitio y secuencia con el de `f5-sites.tsv` (comprobación
            // anterior), de modo que las dosis recontadas desde la verdad son las
            // realizadas desde los TSV. Ninguna mutación aislada puede separarlas.
            // Se conserva a propósito: es el disparador que queda si alguna vez se
            // debilita un eslabón de cualquiera de las dos cadenas. Al sustituir el
            // recuento por `allele_dosages` en una variante defectuosa, la mutación
            // «dosis alterada» dejó de morir por MAN-70.
            let realized_counts = realized_counts(unit, site_index, ploidy);
            if realized_counts != locus.counts {
                return Err(reject(
                    F5OracleCause::Man70,
                    F5_TRUTH_PATH,
                    "recounted_genotype_counts",
                    format!(
                        "el locus {}×{} realiza {realized_counts:?} y la verdad reconta {:?}",
                        site.site_id, unit.population_id, locus.counts
                    ),
                ));
            }
            let denominator = census.individual_count as f64;
            let realized_frequencies: Vec<f64> = realized_counts
                .iter()
                .map(|count| *count as f64 / denominator)
                .collect();
            let copies = census.individual_count * u64::from(declaration.ploidy);
            let (model, probabilities) = match effective_mode {
                F5DeclaredMode::Controlled => (
                    F5GenotypeModel::ControlledHweConditionedCounts,
                    hypergeometric_probabilities(
                        copies,
                        unit.alt_copies[site_index],
                        declaration.ploidy,
                    ),
                ),
                F5DeclaredMode::Probabilistic if declaration.ploidy == 2 && inbreeding > 0.0 => (
                    F5GenotypeModel::DiploidF,
                    diploid_f_probabilities(unit.realized_alt_af[site_index], inbreeding),
                ),
                F5DeclaredMode::Probabilistic => (
                    F5GenotypeModel::Hwe,
                    binomial_probabilities(declaration.ploidy, unit.realized_alt_af[site_index]),
                ),
            };
            let expected_counts: Vec<f64> = probabilities
                .iter()
                .map(|probability| denominator * probability)
                .collect();
            summaries.push(F5InitializationSummary {
                population_id: unit.population_id.clone(),
                block_id: unit.block_id.clone(),
                site_id: site.site_id.clone(),
                ploidy: declaration.ploidy,
                individual_count: census.individual_count,
                model,
                requested_mode: sampling.requested_mode,
                effective_mode,
                input_alt_af: unit.target_alt_af[site_index],
                realized_alt_af: unit.realized_alt_af[site_index],
                expected_genotype_frequencies: probabilities.clone(),
                expected_genotype_probabilities: probabilities,
                expected_genotype_counts: expected_counts,
                realized_genotype_frequencies: realized_frequencies,
                realized_genotype_counts: realized_counts,
                recounted_genotype_frequencies: locus.frequencies.clone(),
                recounted_genotype_counts: locus.counts.clone(),
                recount_matches: true,
            });
        }
    }
    if cursor != recount.loci.len() {
        return Err(reject(
            F5OracleCause::Man70,
            F5_TRUTH_PATH,
            "recount",
            format!(
                "la verdad reconta {} loci y los TSV describen {cursor}",
                recount.loci.len()
            ),
        ));
    }
    Ok(summaries)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn catalog_row(source: &str, decision: &str) -> Vec<String> {
        source
            .lines()
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .find(|line| line.starts_with(&format!("{decision}|")))
            .expect("fila del catálogo")
            .split('|')
            .map(str::to_owned)
            .collect()
    }

    /// El oráculo se valida **antes** de juzgar a nadie: reproduce los cinco
    /// vectores de identidad y todas las subsemillas congeladas.
    #[test]
    fn the_oracle_reproduces_the_two_frozen_catalogs() {
        validate_id_catalog().expect("catálogo de identidades reproducido");
        validate_rng_catalog().expect("catálogo de subsemillas reproducido");
    }

    /// Si la reimplementación de Base32 o del material fuese otra, la
    /// validación anterior no podría pasar: este caso fija el vector literal.
    #[test]
    fn the_identity_oracle_reproduces_the_frozen_block_identifier() {
        let row = catalog_row(ID_CATALOG_TSV, "f5_block");
        let root: [u8; 32] = decode_hex(&row[5])
            .expect("raíz")
            .try_into()
            .expect("32 bytes");
        let oracle = IdOracle::new(
            root,
            "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4",
            "scenario-f5",
            "reference-main",
        )
        .expect("oráculo validado");
        assert_eq!(oracle.block(7).expect("BLK_"), row[8]);
        assert_eq!(
            oracle.site(7, 5).expect("SITE_"),
            catalog_row(ID_CATALOG_TSV, "f5_site")[8]
        );
        assert_eq!(
            oracle
                .haplotype("population-A", 7, 3, "101001")
                .expect("HAP_"),
            catalog_row(ID_CATALOG_TSV, "f5_haplotype")[8]
        );
        assert_eq!(
            oracle
                .founder_copy("population-A", "individual-A", "sequence-1", 0)
                .expect("FCP_"),
            catalog_row(ID_CATALOG_TSV, "f5_founder_copy")[8]
        );
        assert_eq!(
            oracle
                .assignment(
                    "population-A",
                    7,
                    "FCP_3MYOBVN22OO7JPHBBXMUFWFJQQ",
                    "HAP_LWRZFIHLKXNC3GHSPCIH6SJV7A"
                )
                .expect("FHA_"),
            catalog_row(ID_CATALOG_TSV, "f5_block_haplotype_assignment")[8]
        );
    }

    /// La polaridad rederivada sale del vector congelado del catálogo RNG, no
    /// de volver a llamar al productor.
    #[test]
    fn the_subseed_oracle_draws_the_frozen_polarity() {
        let row = catalog_row(RNG_CATALOG_TSV, "polarity");
        let root: [u8; 32] = decode_hex(&row[2])
            .expect("raíz")
            .try_into()
            .expect("32 bytes");
        let first = u64::from_str_radix(&row[5], 16).expect("primera palabra");
        let expected = if (first >> 11) < (1_u64 << 52) {
            "positive"
        } else {
            "negative"
        };
        let oracle = SubseedOracle::new(root, "scenario-f5").expect("oráculo validado");
        assert_eq!(
            oracle
                .polarity("population-A", "BLK_blockA", 3, 5)
                .expect("polaridad"),
            expected
        );
        // El sitio 0 no consume entropía y es siempre `positive`.
        assert_eq!(
            oracle
                .polarity("population-A", "BLK_blockA", 3, 0)
                .expect("polaridad"),
            "positive"
        );
    }

    /// El ALT rederivado sale del vector congelado del catálogo RNG: con la
    /// primera palabra de `site_alt` y el rechazo por muestreo de `n=3`, la
    /// base elegida es la que el productor sortea.
    #[test]
    fn the_subseed_oracle_draws_the_frozen_alternate_base() {
        let row = catalog_row(RNG_CATALOG_TSV, "site_alt");
        let root: [u8; 32] = decode_hex(&row[2])
            .expect("raíz")
            .try_into()
            .expect("32 bytes");
        let first = u64::from_str_radix(&row[5], 16).expect("primera palabra");
        // `n=3` deja `limit = u64::MAX`, de modo que solo la palabra máxima se
        // rechaza y el vector congelado no la trae.
        assert!(first < u64::MAX);
        let oracle = SubseedOracle::new(root, "scenario-f5").expect("oráculo validado");
        for (reference, alternatives) in [
            ('A', ['C', 'G', 'T']),
            ('C', ['A', 'G', 'T']),
            ('G', ['A', 'C', 'T']),
            ('T', ['A', 'C', 'G']),
        ] {
            let expected = alternatives[usize::try_from(first % 3).expect("índice")];
            assert_eq!(
                oracle
                    .alternate_base("BLK_blockA", "SITE_siteA", reference)
                    .expect("ALT"),
                expected,
                "REF={reference}"
            );
        }
    }

    /// El rechazo por muestreo se escribe entero: con `n` que no divide el
    /// rango, la palabra de la cola se descarta en vez de sesgar el resto.
    #[test]
    fn the_rejection_sampling_bound_is_the_one_dtc_089_fixes() {
        for n in [2_u64, 3, 5, 7] {
            let remainder = (u64::MAX % n + 1) % n;
            let limit = (remainder != 0).then(|| u64::MAX - (remainder - 1));
            assert!(limit.is_none_or(|bound| bound % n == 0 || bound == u64::MAX));
        }
        let mut rng = ChaCha12Rng::from_seed([7_u8; 32]);
        for _ in 0..64 {
            assert!(uniform_index(&mut rng, 3) < 3);
        }
        assert_eq!(uniform_index(&mut rng, 1), 0);
    }

    /// Unidad sintética con `H` rangos y `P=2`, para fijar `F_realized` con
    /// casos racionales exactos.
    fn diploid_unit(ranks: &[u32], haplotypes: usize) -> Unit {
        let mut unit = Unit {
            population_id: "P1".to_owned(),
            block_id: "BLK_x".to_owned(),
            block_index: 0,
            haplotype_ids: Vec::new(),
            bits: Vec::new(),
            words: 1,
            haplotypes,
            weights: vec![1.0 / haplotypes as f64; haplotypes],
            sites: 0,
            best_attempt: None,
            best_counters: None,
            polarity: Vec::new(),
            target_alt_af: Vec::new(),
            realized_alt_af: Vec::new(),
            ranks: ranks.to_vec(),
            alt_copies: Vec::new(),
            rank_counts: vec![0; haplotypes],
            heterozygous: 0,
            individuals: 0,
        };
        for pair in ranks.chunks(2) {
            unit.individuals += 1;
            if pair[0] != pair[1] {
                unit.heterozygous += 1;
            }
            for rank in pair {
                unit.rank_counts[*rank as usize] += 1;
            }
        }
        unit
    }

    /// `F_realized` en tres casos racionales exactos, incluida la estimación
    /// negativa que REQ-EVO-52 obliga a conservar.
    #[test]
    fn the_realized_f_reproduces_three_exact_rational_cases() {
        // Cuatro individuos, dos rangos equifrecuentes y la mitad
        // heterocigotos: `H_obs = H_exp = 1/2` y `F = 0`.
        assert_eq!(
            recompute_realized_f(&diploid_unit(&[0, 0, 0, 1, 1, 0, 1, 1], 2)),
            Some(0.0)
        );
        // Ningún heterocigoto con las mismas frecuencias: `F = 1`.
        assert_eq!(
            recompute_realized_f(&diploid_unit(&[0, 0, 0, 0, 1, 1, 1, 1], 2)),
            Some(1.0)
        );
        // Todos heterocigotos: `H_obs = 1`, `H_exp = 1/2` y `F = -1`, que se
        // conserva en vez de recortarse a cero.
        assert_eq!(
            recompute_realized_f(&diploid_unit(&[0, 1, 0, 1, 0, 1, 0, 1], 2)),
            Some(-1.0)
        );
        // Un solo rango deja `H_exp = 0`: no hay F realizado que publicar.
        assert_eq!(recompute_realized_f(&diploid_unit(&[0, 0, 0, 0], 1)), None);
    }

    fn normalization_case() -> (Vec<&'static str>, Vec<f64>, f64, f64, Vec<f64>, f64, f64) {
        let keys = vec!["0", "1", "2", "3"];
        let input = vec![0.5_f64, 0.25, 0.125, 0.125];
        let input_sum = neumaier_total(input.iter().copied());
        let input_delta = (input_sum - 1.0).abs();
        let effective: Vec<f64> = input.iter().map(|value| value / input_sum).collect();
        let output_sum = neumaier_total(effective.iter().copied());
        let output_delta = (output_sum - 1.0).abs();
        (
            keys,
            input,
            input_sum,
            input_delta,
            effective,
            output_sum,
            output_delta,
        )
    }

    /// Las seis cifras de una auditoría DTC-049 se recalculan, no se creen:
    /// perturbar cualquiera de ellas la tumba por su causa.
    #[test]
    fn every_perturbed_dtc_049_figure_is_rejected() {
        let (keys, input, input_sum, input_delta, effective, output_sum, output_delta) =
            normalization_case();
        let run = |algorithm: &str,
                   tolerance: f64,
                   input_values: &[f64],
                   input_sum: f64,
                   input_delta: f64,
                   effective_keys: &[&str],
                   effective_values: &[f64],
                   output_sum: f64,
                   output_delta: f64| {
            check_normalization(
                F5NormalizationSet::HaplotypeWeights,
                algorithm,
                tolerance,
                &keys,
                input_values,
                input_sum,
                input_delta,
                effective_keys,
                effective_values,
                output_sum,
                output_delta,
            )
        };
        run(
            F5_PROFILE_NORMALIZATION_ALGORITHM,
            F5_PROFILE_NORMALIZATION_TOLERANCE,
            &input,
            input_sum,
            input_delta,
            &keys,
            &effective,
            output_sum,
            output_delta,
        )
        .expect("la auditoría íntegra se acepta");

        type Perturbation<'a> = Box<dyn Fn() -> Checked<()> + 'a>;
        let perturbations: Vec<(&str, Perturbation<'_>)> = vec![
            (
                "algoritmo ajeno",
                Box::new(|| {
                    run(
                        "kahan-f64-v1",
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta,
                        &keys,
                        &effective,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "tolerancia ajena",
                Box::new(|| {
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        1e-6,
                        &input,
                        input_sum,
                        input_delta,
                        &keys,
                        &effective,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "input_sum alterado",
                Box::new(|| {
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum + 1e-12,
                        input_delta,
                        &keys,
                        &effective,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "input_delta alterado",
                Box::new(|| {
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta + 1e-12,
                        &keys,
                        &effective,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "valor efectivo alterado",
                Box::new(|| {
                    let mut broken = effective.clone();
                    broken[0] = 0.25;
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta,
                        &keys,
                        &broken,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "orden canónico roto",
                Box::new(|| {
                    let reordered = vec!["1", "0", "2", "3"];
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta,
                        &reordered,
                        &effective,
                        output_sum,
                        output_delta,
                    )
                }),
            ),
            (
                "output_sum alterado",
                Box::new(|| {
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta,
                        &keys,
                        &effective,
                        output_sum + 1e-12,
                        output_delta,
                    )
                }),
            ),
            (
                "output_delta alterado",
                Box::new(|| {
                    run(
                        F5_PROFILE_NORMALIZATION_ALGORITHM,
                        F5_PROFILE_NORMALIZATION_TOLERANCE,
                        &input,
                        input_sum,
                        input_delta,
                        &keys,
                        &effective,
                        output_sum,
                        output_delta + 1e-12,
                    )
                }),
            ),
        ];
        for (name, perturbation) in perturbations {
            let rejection = perturbation()
                .err()
                .unwrap_or_else(|| panic!("«{name}» pasó inadvertida"));
            assert_eq!(
                rejection.cause(),
                F5OracleCause::Normalization,
                "«{name}» falló por otra causa: {}",
                rejection.diagnostic().render_human()
            );
        }
    }

    /// La hipergeométrica del modo controlado suma 1 y respeta sus extremos.
    #[test]
    fn the_conditioned_hypergeometric_is_a_distribution() {
        let probabilities = hypergeometric_probabilities(10, 4, 2);
        assert_eq!(probabilities.len(), 3);
        let total: f64 = probabilities.iter().sum();
        assert!((total - 1.0).abs() < 1e-12, "{probabilities:?}");
        // Con `K=0` la dosis 0 es segura y ninguna otra es posible.
        assert_eq!(hypergeometric_probabilities(10, 0, 2), vec![1.0, 0.0, 0.0]);
        // Con `K=C` lo es la dosis `P`.
        assert_eq!(hypergeometric_probabilities(10, 10, 2), vec![0.0, 0.0, 1.0]);
    }

    /// F diploide en sus dos extremos: `F=0` es HWE y `F=1` no deja
    /// heterocigotos.
    #[test]
    fn the_diploid_f_model_spans_its_extremes() {
        assert_eq!(
            diploid_f_probabilities(0.5, 0.0),
            binomial_probabilities(2, 0.5)
        );
        let fixed = diploid_f_probabilities(0.5, 1.0);
        assert_eq!(fixed[1], 0.0);
        assert!((fixed[0] + fixed[2] - 1.0).abs() < 1e-12);
    }
}

// ---------------------------------------------------------------------------
// Oráculo del paquete de fallo
// ---------------------------------------------------------------------------

/// Lo que una ejecución declara de su fallo controlado y el oráculo contrasta
/// contra los bytes publicados.
///
/// Es el gemelo de [`F5PackageDeclaration`] para el otro estado. Ninguno de
/// sus campos es un acumulador del productor: son las cifras que el manifiesto
/// fallido publica y que, si el paquete no las sostiene, deben rechazarse.
#[derive(Clone, Debug, PartialEq)]
pub struct F5FailureDeclaration {
    pub code: String,
    pub stage: String,
    pub scenario_id: String,
    pub population_id: Option<String>,
    pub block_id: Option<String>,
    /// Los cuatro límites algorítmicos y los tres presupuestos efectivos, en
    /// el orden del esquema.
    pub limits: [u64; 7],
    /// `layout_attempts_used`, `haplotype_attempts_used`, `swap_proposals` y
    /// `swap_accepted`.
    pub counters: [u64; 4],
    /// Los dos objetos opcionales, tal y como el manifiesto los publica.
    ///
    /// Viajan como JSON porque el contrato los define como «`null` u objeto
    /// cerrado» y lo que el oráculo tiene que comprobar es que el diagnóstico
    /// publique **exactamente** lo mismo que su manifiesto, sin reinterpretar
    /// ninguno de sus catorce campos.
    pub best_layout: serde_json::Value,
    pub best_attempt: serde_json::Value,
    pub link: F5DeclaredLink,
}

/// Los diez literales exactos de `f5-failure-v1.schema.json#/$defs/code`.
const F5_FAILURE_CODES: [&str; 10] = [
    "GPS-CONF-101",
    "GPS-FMT-101",
    "GPS-DOMAIN-101",
    "GPS-DOMAIN-102",
    "GPS-DOMAIN-103",
    "GPS-DOMAIN-104",
    "GPS-DOMAIN-105",
    "GPS-DOMAIN-106",
    "GPS-RES-101",
    "GPS-INTERNAL-101",
];

/// Las cinco etapas exactas del contrato.
const F5_FAILURE_STAGES: [&str; 5] = [
    "preflight",
    "layout",
    "haplotypes",
    "founders",
    "publication",
];

/// Las siete claves de `limits`, en el orden del esquema.
const F5_FAILURE_LIMIT_KEYS: [&str; 7] = [
    "max_layout_attempts",
    "max_haplotype_attempts",
    "max_swap_proposals",
    "max_haplotype_cells",
    "max_operation_units",
    "max_working_memory_bytes",
    "max_output_bytes",
];

/// Las cuatro claves de `counters`, en el orden del esquema.
const F5_FAILURE_COUNTER_KEYS: [&str; 4] = [
    "layout_attempts_used",
    "haplotype_attempts_used",
    "swap_proposals",
    "swap_accepted",
];

/// Ruta publicada del diagnóstico de un fallo controlado.
pub const F5_FAILURE_DOCUMENT_PATH: &str = "outputs/f5-failure.json";

/// Audita un paquete **fallido** desde sus bytes.
///
/// El oráculo de F5.4.2 exige los siete TSV y `truth.jsonl`: no tiene modo de
/// fallo, de modo que la mitad «fallo» del criterio de la Revisión F5.4 no
/// estaba demostrada porque no estaba implementada. Esta función la cubre y no
/// comparte ninguna estructura con el productor:
///
/// 1. relee `outputs/f5-failure.json` y comprueba sus doce claves, el
///    enumerado cerrado de `code`, las cinco etapas y `artifacts=[]`;
/// 2. contrasta los **nueve términos comunes** contra la declaración;
/// 3. exige que los contadores quepan en la cota que cada uno tiene: los dos
///    de intentos en su límite homónimo y las propuestas en
///    `max_swap_proposals * haplotype_attempts_used`, porque
///    `max_swap_proposals` presupuesta **un intento** y `swap_proposals` es
///    acumulado (DTC-089 l. 1196 y l. 1847-1848);
/// 4. **enumera el árbol** y exige que no exista ninguno de los ocho ficheros
///    biológicos ni ningún fichero sin declarar;
/// 5. reabre el enlace y recomputa su tamaño y SHA-256.
///
/// # Errors
///
/// [`F5OracleRejection`], cuya [`F5OracleCause`] nombra el contrato incumplido.
pub fn audit_f5_failure_package(
    payload_root: &Path,
    declaration: &F5FailureDeclaration,
) -> Checked<()> {
    let path = F5_FAILURE_DOCUMENT_PATH;
    if !F5_FAILURE_CODES.contains(&declaration.code.as_str()) {
        return Err(reject(
            F5OracleCause::ByteContract,
            path,
            "code",
            format!("{} no pertenece al catálogo cerrado", declaration.code),
        ));
    }
    if !F5_FAILURE_STAGES.contains(&declaration.stage.as_str()) {
        return Err(reject(
            F5OracleCause::ByteContract,
            path,
            "stage",
            format!("{} no es una de las cinco etapas", declaration.stage),
        ));
    }
    let document = read_failure_document(payload_root, path)?;
    let text = |key: &str| -> Checked<String> {
        document
            .get(key)
            .and_then(serde_json::Value::as_str)
            .map(str::to_owned)
            .ok_or_else(|| {
                reject(
                    F5OracleCause::ByteContract,
                    path,
                    key,
                    "falta o no es cadena",
                )
            })
    };
    if text("diagnostic_kind")? != "f5_failure" || text("schema_version")? != "1" {
        return Err(reject(
            F5OracleCause::ByteContract,
            path,
            "diagnostic_kind",
            "clase o versión incompatible con f5-failure-v1",
        ));
    }
    const KEYS: [&str; 12] = [
        "diagnostic_kind",
        "schema_version",
        "code",
        "stage",
        "scenario_id",
        "population_id",
        "block_id",
        "limits",
        "counters",
        "best_layout",
        "best_attempt",
        "artifacts",
    ];
    if document.len() != KEYS.len() || KEYS.iter().any(|key| !document.contains_key(*key)) {
        return Err(reject(
            F5OracleCause::ByteContract,
            path,
            "object",
            format!("el objeto cerrado exige {KEYS:?}"),
        ));
    }
    if document
        .get("artifacts")
        .and_then(serde_json::Value::as_array)
        .is_none_or(|artifacts| !artifacts.is_empty())
    {
        return Err(reject(
            F5OracleCause::Link,
            path,
            "artifacts",
            "el diagnóstico publica una lista de artefactos no vacía",
        ));
    }
    // Los nueve términos comunes, contra la declaración.
    for (key, declared) in [
        ("code", declaration.code.clone()),
        ("stage", declaration.stage.clone()),
        ("scenario_id", declaration.scenario_id.clone()),
    ] {
        if text(key)? != declared {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                key,
                format!(
                    "el diagnóstico publica {} y se declara {declared}",
                    text(key)?
                ),
            ));
        }
    }
    for (key, declared) in [
        ("population_id", declaration.population_id.as_deref()),
        ("block_id", declaration.block_id.as_deref()),
    ] {
        let published = document.get(key).and_then(|value| match value {
            serde_json::Value::Null => None,
            other => other.as_str(),
        });
        if published != declared {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                key,
                format!("el diagnóstico publica {published:?} y se declara {declared:?}"),
            ));
        }
    }
    let numbers = |key: &str, keys: &[&str]| -> Checked<Vec<u64>> {
        let object = document
            .get(key)
            .and_then(serde_json::Value::as_object)
            .ok_or_else(|| {
                reject(
                    F5OracleCause::ByteContract,
                    path,
                    key,
                    "falta o no es objeto",
                )
            })?;
        if object.len() != keys.len() || keys.iter().any(|name| !object.contains_key(*name)) {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                key,
                format!("el objeto cerrado exige {keys:?}"),
            ));
        }
        keys.iter()
            .map(|name| {
                object
                    .get(*name)
                    .and_then(serde_json::Value::as_u64)
                    .ok_or_else(|| reject(F5OracleCause::ByteContract, path, name, "no es un u64"))
            })
            .collect()
    };
    let limits = numbers("limits", &F5_FAILURE_LIMIT_KEYS)?;
    if limits != declaration.limits {
        return Err(reject(
            F5OracleCause::ResourceBound,
            path,
            "limits",
            format!(
                "el diagnóstico publica {limits:?} y se declaran {:?}",
                declaration.limits
            ),
        ));
    }
    let counters = numbers("counters", &F5_FAILURE_COUNTER_KEYS)?;
    if counters != declaration.counters {
        return Err(reject(
            F5OracleCause::Attempt,
            path,
            "counters",
            format!(
                "el diagnóstico publica {counters:?} y se declaran {:?}",
                declaration.counters
            ),
        ));
    }
    for (key, declared) in [
        ("best_layout", &declaration.best_layout),
        ("best_attempt", &declaration.best_attempt),
    ] {
        let published = document.get(key).unwrap_or(&serde_json::Value::Null);
        if published != declared {
            return Err(reject(
                F5OracleCause::Attempt,
                path,
                key,
                "el diagnóstico y el manifiesto publican objetos distintos",
            ));
        }
        if !matches!(
            published,
            serde_json::Value::Null | serde_json::Value::Object(_)
        ) {
            return Err(reject(
                F5OracleCause::ByteContract,
                path,
                key,
                "el contrato admite null u objeto cerrado",
            ));
        }
    }
    // Los dos contadores de intentos caben cada uno en **su propio** límite:
    // los de layout en `max_layout_attempts` y los de haplotipos en
    // `max_haplotype_attempts`. Son los únicos dos cuya correspondencia con
    // `limits` es contador contra su límite homónimo.
    for (index, name) in [
        (0_usize, "layout_attempts_used"),
        (1, "haplotype_attempts_used"),
    ] {
        if counters[index] > limits[index] {
            return Err(resource_reject(
                path,
                format!(
                    "{name}={} excede su límite efectivo {}",
                    counters[index], limits[index]
                ),
            ));
        }
    }
    // `swap_proposals` **no** se compara contra `max_swap_proposals`, y por eso
    // sale del bucle anterior: que los tres compartieran recorrido posicional
    // hacía pasar por simetría lo que era una coincidencia de índices.
    //
    // DTC-089 l. 1196 fija `max_swap_proposals` como presupuesto de **un
    // intento** —«Un intento termina … o al consumir `max_swap_proposals`»— y
    // l. 1847-1848 enumera `swap_proposals` junto a `haplotype_attempts_used`,
    // que es **acumulado** sobre los intentos; leer el tercero como «de un
    // intento» rompería la coherencia de su propio grupo. La cota que el
    // contrato sostiene es por tanto el producto: cada uno de los
    // `haplotype_attempts_used` intentos consume como mucho
    // `max_swap_proposals` propuestas. Sin el producto, todo agotamiento de
    // búsqueda —donde el acumulado vale exactamente
    // `max_haplotype_attempts * max_swap_proposals`— se rechazaría con una
    // causa de recurso falsa.
    //
    // La cota es **ajustada**, no holgada: el productor acumula los dos
    // contadores en el mismo recorrido —un intento sumado y sus propuestas
    // sumadas por cada registro de la unidad que el diagnóstico nombra—, de
    // modo que `swap_proposals` es la suma de `haplotype_attempts_used`
    // sumandos y ninguno puede pasar de `max_swap_proposals`. Un fallo que no
    // nace de la búsqueda deja los dos en cero y la cota sigue valiendo.
    //
    // El producto se hace con `checked_mul`: el techo contractual es
    // `10^7 * 10^4 = 10^11` y cabe de sobra en `u64`, pero un desbordamiento
    // silencioso aquí abriría la cota en lugar de cerrarla.
    let proposal_bound = limits[2].checked_mul(counters[1]).ok_or_else(|| {
        resource_reject(
            path,
            format!(
                "la cota de propuestas desborda: max_swap_proposals={} por haplotype_attempts_used={}",
                limits[2], counters[1]
            ),
        )
    })?;
    if counters[2] > proposal_bound {
        return Err(resource_reject(
            path,
            format!(
                "swap_proposals={} excede max_swap_proposals={} por haplotype_attempts_used={}",
                counters[2], limits[2], counters[1]
            ),
        ));
    }
    if counters[3] > counters[2] {
        return Err(reject(
            F5OracleCause::Attempt,
            path,
            "swap_accepted",
            "se aceptan más propuestas de las que se hicieron",
        ));
    }
    // El árbol publicado, enumerado: ni un artefacto biológico ni un fichero
    // sin declarar.
    let mut found = BTreeSet::new();
    collect_failure_files(payload_root, payload_root, path, &mut found)?;
    let declared = BTreeSet::from(["manifest.json".to_owned(), path.to_owned()]);
    if found != declared {
        return Err(reject(
            F5OracleCause::Link,
            path,
            "payload",
            format!("el árbol publicado es {found:?} y un fallo publica {declared:?}"),
        ));
    }
    // Y el enlace declarado corresponde a los bytes.
    let link = &declaration.link;
    if link.artifact_id != "f5-00-failure-v1"
        || link.artifact_kind != "f5_failure"
        || link.schema_version != "1"
        || link.path != path
    {
        return Err(reject(
            F5OracleCause::Link,
            path,
            "artifact_id",
            format!(
                "el enlace declara {}/{}/{}/{}",
                link.artifact_id, link.artifact_kind, link.schema_version, link.path
            ),
        ));
    }
    let (size_bytes, sha256) = measure_artifact(payload_root, path)?;
    if size_bytes != link.size_bytes || sha256 != link.sha256 {
        return Err(reject(
            F5OracleCause::Link,
            path,
            "sha256",
            format!(
                "el enlace declara {} bytes y {} y el fichero mide {size_bytes} con {sha256}",
                link.size_bytes, link.sha256
            ),
        ));
    }
    Ok(())
}

fn read_failure_document(
    payload_root: &Path,
    path: &str,
) -> Checked<serde_json::Map<String, serde_json::Value>> {
    let reader = open_artifact(payload_root, path)?;
    let value: serde_json::Value = serde_json::from_reader(reader).map_err(|error| {
        reject(
            F5OracleCause::ByteContract,
            path,
            "json",
            format!("no es un documento JSON legible: {error}"),
        )
    })?;
    match value {
        serde_json::Value::Object(object) => Ok(object),
        _ => Err(reject(
            F5OracleCause::ByteContract,
            path,
            "json",
            "la raíz no es un objeto",
        )),
    }
}

fn collect_failure_files(
    payload_root: &Path,
    directory: &Path,
    path: &str,
    found: &mut BTreeSet<String>,
) -> Checked<()> {
    let entries = std::fs::read_dir(directory).map_err(|error| {
        reject(
            F5OracleCause::ArtifactMissing,
            path,
            "payload",
            format!("no se pudo recorrer {}: {error}", directory.display()),
        )
    })?;
    for entry in entries {
        // Inalcanzable en la práctica: `ReadDir` solo da error de iteración ante un
        // fallo de entrada y salida del sistema de ficheros a mitad de un listado,
        // que una prueba no puede provocar de forma determinista.
        let entry = entry.map_err(|error| {
            reject(
                F5OracleCause::ArtifactMissing,
                path,
                "payload",
                format!("entrada ilegible en {}: {error}", directory.display()),
            )
        })?;
        let child = entry.path();
        // `file_type` no sigue enlaces; `is_dir` sí. Un paquete declarado por
        // huella no admite un enlace en su árbol.
        // Inalcanzable en la práctica: en Linux `file_type` sale del tipo que el
        // propio listado ya trae y solo consulta el disco si el sistema de ficheros
        // no lo rellena; falla si el fichero desaparece entre el listado y la
        // consulta, una carrera que una prueba no puede provocar de forma
        // determinista.
        let kind = entry.file_type().map_err(|error| {
            reject(
                F5OracleCause::ArtifactMissing,
                path,
                "payload",
                format!("no se pudo tipar {}: {error}", child.display()),
            )
        })?;
        if kind.is_dir() {
            collect_failure_files(payload_root, &child, path, found)?;
            continue;
        }
        if !kind.is_file() {
            return Err(reject(
                F5OracleCause::Link,
                path,
                "payload",
                format!(
                    "{} no es un fichero regular ni un directorio",
                    child.display()
                ),
            ));
        }
        // Inalcanzable: `child` sale de un listado de `directory`, y `directory` es
        // `payload_root` o uno de sus descendientes, de modo que siempre cuelga de
        // él.
        let relative = child.strip_prefix(payload_root).map_err(|_| {
            reject(
                F5OracleCause::ArtifactMissing,
                path,
                "payload",
                format!("{} no cuelga de la carga útil", child.display()),
            )
        })?;
        found.insert(relative.to_string_lossy().replace('\\', "/"));
    }
    Ok(())
}

/// Pruebas de la deuda de cobertura de los jueces (DT-001): cada guarda del
/// oráculo que ninguna prueba disparaba, con su motivo exacto.
#[cfg(test)]
mod guard_tests {
    use super::*;

    /// El texto del detalle de un rechazo de catálogo.
    fn detail_of(rejection: &F5OracleRejection) -> String {
        rejection
            .diagnostic()
            .context()
            .get("detail")
            .expect("todo rechazo lleva detalle")
            .clone()
    }

    /// El catálogo v1 con una fila reescrita: `edit` recibe los nueve campos de
    /// la primera fila de decisión (`layout_length`) y los deja como quiera.
    fn catalog_with_first_row_edited(edit: impl FnOnce(&mut Vec<String>)) -> String {
        let mut edit = Some(edit);
        let mut lines = Vec::new();
        for line in RNG_CATALOG_TSV.lines() {
            match edit.take() {
                Some(edit) if !line.is_empty() && !line.starts_with('#') => {
                    let mut fields: Vec<String> = line.split('|').map(str::to_owned).collect();
                    edit(&mut fields);
                    lines.push(fields.join("|"));
                }
                other => {
                    edit = other;
                    lines.push(line.to_owned());
                }
            }
        }
        assert!(edit.is_none(), "el catálogo tiene alguna fila de decisión");
        lines.join("\n") + "\n"
    }

    /// **La autovalidación del catálogo RNG, guarda por guarda.**
    ///
    /// El oráculo se prueba contra el catálogo congelado antes de juzgar a
    /// nadie, pero con el catálogo íntegro todas sus comprobaciones pasan, de
    /// modo que ninguna había dicho nunca «no». Aquí cada una recibe un
    /// catálogo dañado en **exactamente** el campo que ella vigila, y la
    /// prueba exige el motivo exacto: quitar la comprobación o cambiar su
    /// condición deja pasar el catálogo o lo rechaza por otra razón.
    ///
    /// Las siete pasan por el mismo camino porque la función es la misma que el
    /// oráculo usa con el catálogo real: no se prueba una réplica.
    #[test]
    fn the_frozen_rng_catalog_check_refuses_each_damaged_field() {
        validate_frozen_rng_catalog(RNG_CATALOG_TSV, F5_RNG_CATALOG_DECISIONS)
            .expect("el catálogo íntegro se reproduce");
        let flip = |field: &mut String| {
            // Cambia el primer carácter hexadecimal por otro distinto.
            let first = field.remove(0);
            field.insert(0, if first == '0' { '1' } else { '0' });
        };
        type RowEdit = Box<dyn Fn(&mut Vec<String>)>;
        let cases: [(&str, RowEdit, &str); 7] = [
            (
                "una fila con ocho campos",
                Box::new(|fields| {
                    fields.pop();
                }),
                "no tiene 9 campos",
            ),
            (
                "una raíz que no es hexadecimal",
                Box::new(|fields| fields[2] = "zz".to_owned()),
                "raíz del catálogo RNG ilegible",
            ),
            (
                "una raíz hexadecimal de 31 bytes",
                Box::new(|fields| {
                    fields[2].truncate(62);
                }),
                "raíz del catálogo RNG ilegible",
            ),
            (
                "un material que no es el de sus componentes",
                Box::new(move |fields| flip(&mut fields[3])),
                "el material rederivado",
            ),
            (
                "una subsemilla que no es la del material",
                Box::new(move |fields| flip(&mut fields[4])),
                "la subsemilla rederivada",
            ),
            (
                "una palabra que no es hexadecimal",
                Box::new(|fields| fields[6] = "no-es-hex".to_owned()),
                "palabra del catálogo RNG ilegible",
            ),
            (
                "una palabra que no es la que el ChaCha12 produce",
                Box::new(move |fields| flip(&mut fields[7])),
                "la palabra 2 de",
            ),
        ];
        for (name, edit, expected) in cases {
            let catalog = catalog_with_first_row_edited(edit);
            assert_ne!(catalog, RNG_CATALOG_TSV, "«{name}» no cambió el catálogo");
            let rejection = validate_frozen_rng_catalog(&catalog, F5_RNG_CATALOG_DECISIONS)
                .err()
                .unwrap_or_else(|| panic!("«{name}» pasó inadvertido"));
            assert_eq!(rejection.cause(), F5OracleCause::Randomness, "{name}");
            let detail = detail_of(&rejection);
            assert!(
                detail.contains(expected),
                "«{name}» se rechazó por otra razón: {detail}"
            );
        }
        // El número de decisiones se exige exacto en las dos direcciones: un
        // catálogo recortado y uno ampliado se rechazan por su propio recuento.
        for declared in [F5_RNG_CATALOG_DECISIONS - 1, F5_RNG_CATALOG_DECISIONS + 1] {
            let rejection = validate_frozen_rng_catalog(RNG_CATALOG_TSV, declared)
                .expect_err("un recuento de decisiones que no es el del catálogo");
            assert!(
                detail_of(&rejection).contains("el contrato fija"),
                "{}",
                detail_of(&rejection)
            );
        }
    }

    /// Los componentes de una fila de catálogo son pares `nombre=valor`: un par
    /// sin `=` no es una fila legible.
    ///
    /// La comprobación pasa por `catalog_components`, la misma que lee las
    /// filas del catálogo real, y se ejercita también a través de la validación
    /// completa para que la prueba muera si esa validación dejara de llamarla.
    #[test]
    fn a_catalog_component_without_an_equals_sign_is_refused() {
        let rejection = catalog_components("scenario_id=s;stage")
            .expect_err("un componente sin `=` no es legible");
        assert!(
            detail_of(&rejection).contains("componente sin `=`"),
            "{}",
            detail_of(&rejection)
        );
        assert_eq!(
            catalog_components("a=1;b=2").expect("dos pares legibles"),
            vec![("a", "1".to_owned()), ("b", "2".to_owned())]
        );
        let catalog =
            catalog_with_first_row_edited(|fields| fields[1] = "scenario_id;stage".to_owned());
        let rejection = validate_frozen_rng_catalog(&catalog, F5_RNG_CATALOG_DECISIONS)
            .expect_err("una fila con un componente sin `=`");
        assert!(
            detail_of(&rejection).contains("componente sin `=`"),
            "{}",
            detail_of(&rejection)
        );
    }

    /// Un espectro sin bandas no tiene marginal objetivo que sortear.
    #[test]
    fn a_marginal_target_needs_at_least_one_spectrum_band() {
        let oracle = SubseedOracle::new([7_u8; 32], "scenario-f5").expect("oráculo validado");
        let band = F5DeclaredSpectrumBand {
            alt_af_min: 0.25,
            alt_af_max: 0.75,
            weight: 1.0,
        };
        oracle
            .target_alt_af("P1", "BLK_a", "SITE_a", &[band])
            .expect("una banda basta");
        let rejection = oracle
            .target_alt_af("P1", "BLK_a", "SITE_a", &[])
            .expect_err("sin bandas no hay marginal que sortear");
        assert!(
            detail_of(&rejection).contains("no trae ninguna banda"),
            "{}",
            detail_of(&rejection)
        );
    }

    /// Una lambda fuera del dominio del muestreador de recuentos no se sortea:
    /// se rechaza en vez de pinzar un resultado inventado.
    ///
    /// Ejercita solo el rechazo **por la lambda** (negativa, no numérica,
    /// infinita). El otro camino que acaba en este mismo mensaje, un candidato
    /// por encima de 2^53, se prueba llamando a `poisson_count` directamente
    /// en `poisson_count_rejects_a_candidate_above_two_to_the_53` (desde aquí
    /// haría falta encontrar un identificador de bloque cuya subsemilla lo
    /// produzca). Desde DT-010 un candidato negativo no llega a este mensaje:
    /// se repite.
    #[test]
    fn a_haplotype_count_needs_a_lambda_the_sampler_admits() {
        let oracle = SubseedOracle::new([7_u8; 32], "scenario-f5").expect("oráculo validado");
        let with_lambda = |lambda: f64| F5DeclaredHaplotypeCountDistribution {
            lambda,
            minimum: 2,
            maximum: 8,
        };
        oracle
            .block_haplotype_count("P1", "BLK_a", with_lambda(3.0))
            .expect("una lambda del dominio");
        for lambda in [-1.0, f64::NAN, f64::INFINITY] {
            let rejection = oracle
                .block_haplotype_count("P1", "BLK_a", with_lambda(lambda))
                .expect_err("una lambda fuera del dominio");
            assert!(
                detail_of(&rejection).contains("no admite el muestreador"),
                "lambda={lambda}: {}",
                detail_of(&rejection)
            );
        }
    }

    /// Fuente de palabras que apunta lo que entrega: alimenta al productor y
    /// deja la cinta para saber cuántas vueltas dio.
    struct RecordingTape {
        rng: ChaCha12Rng,
        words: Vec<u64>,
    }

    impl crate::randomness::RandomWordSource for RecordingTape {
        fn next_word(&mut self) -> u64 {
            let word = self.rng.next_u64();
            self.words.push(word);
            word
        }
    }

    /// Cuántas de las vueltas que dio el muestreador (dos palabras cada una)
    /// propusieron un candidato negativo. Rederiva el candidato de la cinta con
    /// la fórmula del método, sin pasar por ninguno de los dos muestreadores.
    fn negative_candidates_on(words: &[u64], lambda: f64) -> usize {
        let open_unit = |word: u64| (((word >> 12) as f64) + 0.5) * (1.0 / 4_503_599_627_370_496.0);
        let b = 0.931 + 2.53 * libm::sqrt(lambda);
        let a = -0.059 + 0.02483 * b;
        words
            .chunks_exact(2)
            .filter(|pair| {
                let u = open_unit(pair[0]) - 0.5;
                let us = 0.5 - u.abs();
                libm::floor((2.0 * a / us + b) * u + lambda + 0.43) < 0.0
            })
            .count()
    }

    /// DT-010: un candidato negativo del método de rechazo se **repite**, igual
    /// que hace el productor. El oráculo lo devolvía como `None` y convertía un
    /// escenario correcto en un rechazo de la lambda.
    ///
    /// La misma cinta (derivada de la semilla) alimenta al productor y al
    /// oráculo, con la media 12 y la 20 (la de DT-010), y deben coincidir en
    /// todas. La prueba exige además que **haya habido** casos con candidato
    /// negativo: sin ellos la coincidencia no probaría nada sobre la regla.
    /// Las dos medias no son comprobaciones independientes: con las 256
    /// semillas, las cintas con candidato negativo resultan ser las mismas en
    /// ambas (69, 86, 89, 113, 160, 193 y 233).
    #[test]
    fn poisson_count_repeats_a_negative_candidate_like_the_producer() {
        for lambda in [12.0_f64, 20.0] {
            let mut with_negative_candidate = 0_usize;
            for seed in 0_u64..256 {
                let mut tape = RecordingTape {
                    rng: ChaCha12Rng::seed_from_u64(seed),
                    words: Vec::new(),
                };
                let produced =
                    crate::randomness::poisson_count_with_source(&mut tape, lambda, None)
                        .expect("el productor repite el candidato negativo, no falla");
                assert_eq!(
                    tape.words.len() % 2,
                    0,
                    "lambda={lambda} seed={seed}: la rama de rejilla consume palabras por pares"
                );
                if negative_candidates_on(&tape.words, lambda) > 0 {
                    with_negative_candidate += 1;
                }
                let mut rng = ChaCha12Rng::seed_from_u64(seed);
                assert_eq!(
                    poisson_count(&mut rng, lambda),
                    Some(produced),
                    "lambda={lambda} seed={seed}: el oráculo debe coincidir con el productor"
                );
                // Mismo consumo de palabras: la repetición no cambia la cinta.
                assert_eq!(
                    rng.next_u64(),
                    tape.rng.next_u64(),
                    "lambda={lambda} seed={seed}: el oráculo consumió otra cantidad de palabras"
                );
            }
            assert!(
                with_negative_candidate >= 1,
                "lambda={lambda}: ninguna cinta propuso un candidato negativo; el vector no ejercita la regla"
            );
            eprintln!(
                "lambda={lambda}: {with_negative_candidate}/256 cintas con candidato negativo"
            );
        }
    }

    /// El límite superior del candidato (2^53) se aplica igual en el productor
    /// y en el oráculo. Hace falta la media máxima, 2^52: con ella la
    /// probabilidad de pasar de 2^53 por vuelta es de orden 2^-30, y con medias
    /// pequeñas es inalcanzable en la práctica.
    ///
    /// Dos grupos de semillas, encontrados por búsqueda fuera de línea. Las
    /// semillas 0 a 7 **coinciden** en productor y oráculo con un candidato
    /// dentro del límite; la semilla 0 lo da por encima de 2^52 (mata una cota
    /// bajada a 2^52). La semilla 1_392_399_649 propone en su primera vuelta un
    /// candidato por encima de 2^53: el oráculo debe dar `None` y el productor
    /// `Err(NumericFailure)` (mata un oráculo sin cota superior).
    #[test]
    fn poisson_count_rejects_a_candidate_above_two_to_the_53() {
        const MAX_LAMBDA: f64 = 4_503_599_627_370_496.0;
        for seed in 0_u64..8 {
            let mut tape = RecordingTape {
                rng: ChaCha12Rng::seed_from_u64(seed),
                words: Vec::new(),
            };
            let produced =
                crate::randomness::poisson_count_with_source(&mut tape, MAX_LAMBDA, None)
                    .expect("estas semillas no proponen un candidato por encima de 2^53");
            let mut rng = ChaCha12Rng::seed_from_u64(seed);
            assert_eq!(
                poisson_count(&mut rng, MAX_LAMBDA),
                Some(produced),
                "seed={seed}: el oráculo debe coincidir con el productor"
            );
        }
        for seed in [1_392_399_649_u64, 2_723_730_605, 3_350_934_628] {
            let mut tape = RecordingTape {
                rng: ChaCha12Rng::seed_from_u64(seed),
                words: Vec::new(),
            };
            assert_eq!(
                crate::randomness::poisson_count_with_source(&mut tape, MAX_LAMBDA, None),
                Err(crate::randomness::RandomSamplingError::NumericFailure),
                "seed={seed}: el productor rechaza el candidato por encima de 2^53"
            );
            assert_eq!(
                tape.words.len(),
                2,
                "seed={seed}: falla en la primera vuelta"
            );
            let mut rng = ChaCha12Rng::seed_from_u64(seed);
            assert_eq!(
                poisson_count(&mut rng, MAX_LAMBDA),
                None,
                "seed={seed}: el oráculo rechaza el candidato por encima de 2^53"
            );
        }
    }

    /// La auditoría DTC-049 con cifras internamente coherentes, para que cada
    /// caso rompa una sola cosa.
    fn audit_of(
        keys: &[&str],
        input: &[f64],
        effective_keys: &[&str],
        effective: &[f64],
    ) -> Checked<()> {
        let input_sum = neumaier_total(input.iter().copied());
        let output_sum = neumaier_total(effective.iter().copied());
        check_normalization(
            F5NormalizationSet::HaplotypeWeights,
            F5_PROFILE_NORMALIZATION_ALGORITHM,
            F5_PROFILE_NORMALIZATION_TOLERANCE,
            keys,
            input,
            input_sum,
            (input_sum - 1.0).abs(),
            effective_keys,
            effective,
            output_sum,
            (output_sum - 1.0).abs(),
        )
    }

    /// **La auditoría DTC-049, guardas de forma y de dominio.**
    ///
    /// Las seis cifras perturbadas ya morían; estas tres guardas miran otra
    /// cosa: que haya valores de entrada, que la entrada se pueda normalizar
    /// y que la salida conserve todas las claves. Cada caso deja coherentes
    /// las cifras que no toca, de modo que solo su guarda puede cazarlo.
    #[test]
    fn the_normalization_audit_refuses_a_malformed_shape_by_its_own_reason() {
        let keys = ["0", "1"];
        audit_of(&keys, &[0.5, 0.5], &keys, &[0.5, 0.5]).expect("la auditoría íntegra se acepta");
        let cases: [(&str, Checked<()>, &str); 5] = [
            (
                "sin valores de entrada",
                audit_of(&[], &[], &[], &[]),
                "no declara sus valores de entrada",
            ),
            (
                "menos claves de entrada que valores",
                audit_of(&["0"], &[0.5, 0.5], &keys, &[0.5, 0.5]),
                "no declara sus valores de entrada",
            ),
            (
                "una entrada que no suma uno dentro de la tolerancia",
                audit_of(&keys, &[0.5, 0.25], &keys, &[2.0 / 3.0, 1.0 / 3.0]),
                "no se normaliza bajo DTC-049",
            ),
            (
                "una salida con una clave de menos",
                audit_of(&keys, &[0.5, 0.5], &["0"], &[0.5, 0.5]),
                "pierde o inventa claves",
            ),
            (
                "una salida con un valor de menos",
                audit_of(&keys, &[0.5, 0.5], &keys, &[0.5]),
                "pierde o inventa claves",
            ),
        ];
        for (name, outcome, expected) in cases {
            let rejection = outcome
                .err()
                .unwrap_or_else(|| panic!("«{name}» pasó inadvertida"));
            assert_eq!(rejection.cause(), F5OracleCause::Normalization, "{name}");
            let detail = detail_of(&rejection);
            assert!(
                detail.contains(expected),
                "«{name}» se rechazó por otra razón: {detail}"
            );
        }
    }

    /// **Los pesos de un bloque: cuántos rangos admite el perfil y qué suman.**
    ///
    /// Un bloque publica entre dos filas y tantas como pesos declara el perfil,
    /// y el prefijo que se renormaliza tiene que sumar algo positivo y finito.
    #[test]
    fn the_block_weights_refuse_a_rank_count_and_a_prefix_the_profile_cannot_give() {
        // Un prefijo válido se renormaliza y el vector completo se devuelve tal cual.
        assert_eq!(
            oracle_block_weights(&[0.5, 0.25, 0.25], 3).expect("el vector completo"),
            vec![0.5, 0.25, 0.25]
        );
        let renormalized =
            oracle_block_weights(&[0.5, 0.25, 0.25], 2).expect("un prefijo renormalizado");
        assert_eq!(renormalized, vec![0.5 / 0.75, 0.25 / 0.75]);
        let cases: [(&str, Checked<Vec<f64>>, &str); 4] = [
            (
                "un solo rango",
                oracle_block_weights(&[0.5, 0.5], 1),
                "publica 1 filas",
            ),
            (
                "más rangos que pesos",
                oracle_block_weights(&[0.5, 0.5], 3),
                "publica 3 filas",
            ),
            (
                "un prefijo que suma cero",
                oracle_block_weights(&[0.0, 0.0, 1.0], 2),
                "el prefijo de 2 pesos suma 0",
            ),
            (
                "un prefijo que no es finito",
                oracle_block_weights(&[f64::INFINITY, 0.5, 0.5], 2),
                "el prefijo de 2 pesos suma",
            ),
        ];
        for (name, outcome, expected) in cases {
            let rejection = outcome
                .err()
                .unwrap_or_else(|| panic!("«{name}» pasó inadvertida"));
            assert_eq!(rejection.cause(), F5OracleCause::Haplotype, "{name}");
            let detail = detail_of(&rejection);
            assert!(
                detail.contains(expected),
                "«{name}» se rechazó por otra razón: {detail}"
            );
        }
    }

    /// El catálogo de identidades con una fila F5 reescrita (`edit` recibe los
    /// nueve campos de la primera) o con una fila F5 menos si `edit` es `None`.
    type RowEditRef<'a> = Option<&'a dyn Fn(&mut Vec<String>)>;
    fn id_catalog_with_first_f5_row(edit: RowEditRef<'_>) -> String {
        let mut done = false;
        let mut lines = Vec::new();
        for line in ID_CATALOG_TSV.lines() {
            if !done && line.starts_with("f5_") {
                done = true;
                if let Some(edit) = edit {
                    let mut fields: Vec<String> = line.split('|').map(str::to_owned).collect();
                    edit(&mut fields);
                    lines.push(fields.join("|"));
                }
            } else {
                lines.push(line.to_owned());
            }
        }
        assert!(done, "el catálogo tiene alguna fila F5");
        lines.join("\n") + "\n"
    }

    /// **La autovalidación del catálogo de identidades, guarda por guarda.**
    ///
    /// Cada rechazo recibe un catálogo dañado en exactamente el campo que vigila
    /// y la prueba exige el motivo exacto, como con el catálogo RNG.
    #[test]
    fn the_id_catalog_check_refuses_each_damaged_field() {
        validate_id_catalog_text(ID_CATALOG_TSV).expect("el catálogo íntegro se reproduce");
        let flip = |field: &mut String| {
            let first = field.remove(0);
            field.insert(0, if first == '0' { '1' } else { '0' });
        };
        type RowEdit = Box<dyn Fn(&mut Vec<String>)>;
        let cases: [(&str, RowEdit, &str); 6] = [
            (
                "una fila con ocho campos",
                Box::new(|fields| {
                    fields.pop();
                }),
                "no tiene 9 campos",
            ),
            (
                "una raíz que no es hexadecimal",
                Box::new(|fields| fields[5] = "zz".to_owned()),
                "raíz del catálogo de IDs ilegible",
            ),
            (
                "un material que no es el de sus componentes",
                Box::new(move |fields| flip(&mut fields[6])),
                "el material rederivado",
            ),
            (
                "un digest que no es el del material",
                Box::new(move |fields| flip(&mut fields[7])),
                "el digest rederivado",
            ),
            (
                "un identificador que no es el del digest",
                Box::new(|fields| fields[8].push('X')),
                "el identificador rederivado",
            ),
            (
                "unos componentes sin `=`",
                Box::new(|fields| fields[4] = "class;otra".to_owned()),
                "componente sin `=`",
            ),
        ];
        for (name, edit, expected) in cases {
            let catalog = id_catalog_with_first_f5_row(Some(&*edit));
            assert_ne!(catalog, ID_CATALOG_TSV, "«{name}» no cambió el catálogo");
            let rejection = validate_id_catalog_text(&catalog)
                .err()
                .unwrap_or_else(|| panic!("«{name}» pasó inadvertido"));
            assert_eq!(rejection.cause(), F5OracleCause::Randomness, "{name}");
            assert!(
                detail_of(&rejection).contains(expected),
                "«{name}» se rechazó por otra razón: {}",
                detail_of(&rejection)
            );
        }
        let short = id_catalog_with_first_f5_row(None);
        let rejection = validate_id_catalog_text(&short).expect_err("una clase F5 de menos");
        assert!(
            detail_of(&rejection).contains("declara 4 clases F5"),
            "{}",
            detail_of(&rejection)
        );
    }
}
