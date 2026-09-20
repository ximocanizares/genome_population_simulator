//! Admisión y configuración efectiva mínima de F1.

use std::collections::BTreeSet;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::adapters::{
    admit_toml_header, admitted_is_f5_general_v3, parse_declared_config, parse_declared_f5_general,
    parse_declared_f5_scenario, validate_declared_f5_pair,
};
use crate::biological_contamination::BiologicalContaminationSettings;
use crate::chimera::ChimeraSettings;
use crate::configuration::{
    DeclaredF5GeneralConfig, DeclaredF5Scenario, DeclaredGeneralConfig, DeclaredPopulationSnvs,
    DeclaredReductionTarget, DeclaredReferenceSource, DeclaredTechnicalEventTarget,
    ExpectedConfigKind, SequenceClass,
};
use crate::coordinates::{ClosedInterval1, GenomicInterval};
use crate::diagnostics::Diagnostic;
use crate::error_quality::SequencingErrorsMode;
use crate::exact_decimal::ExactDecimal;
use crate::f5_admission::{
    AdmissionPhase, AdmittedInputIdentity, BoundedInput, F5AdmissionReceipt, F5ResourceLedger,
    InputRole, LedgerCharge, LedgerLimits, OperationCost, OutputReservation, PhysicalInputIdentity,
    WorkingMemoryPeak,
};
use crate::f5_preflight::{
    F5_HARD_OPERATION_LIMIT, F5_HARD_OUTPUT_LIMIT, F5_HARD_WORKING_MEMORY_LIMIT,
};
use crate::f5_profile::{F5Profile, load_f5_profile_under_budget};
use crate::haplotype::{HaplotypePhase, Ploidy};
use crate::identity::{
    AutomaticIdRequest, EventId, ExecutionId, FounderKey, IdClass, IdentityRegistry, IndividualId,
    PopulationId, ReferenceId, RepeatLocusId, SequenceId,
};
use crate::population_snvs::{
    FrequencyAssignmentMode, FrequencySampling, FrequencySpectrumComponent,
    GeneratedFrequencySpectrum, GeneratedLocusCount, GeneratedPopulationFrequency,
    GeneratedSpatialMode, GeneratedSubstitutionProfile, GeneratedWeightedRegion,
    PopulationSnvsGenerationConfig, validate_generation_config,
};
use crate::randomness::{
    EntropySource, ReproducibilityMetadata, RuntimeSubseedScope, Subseed, SystemEntropy,
    derive_subseed, neumaier_sum, resolve_root_seed,
};
use crate::reference_generation::normalize_composition_probabilities;
use crate::tandem_repeat::{RepeatMotif, TandemRepeatLocus, TandemRepeatTract};
use crate::technical_profile::GPS_HISEQX_PCRFREE_PE150_V1;
use crate::vcf_contract::VcfOutputMode;

/// Antes de poder resolver el límite declarado por el escenario, el TOML
/// general debe poder leerse sin reservar una cantidad arbitraria de memoria.
/// La frontera CLI aplica esta misma cota antes de clasificar un general v2 o
/// v3: la clase todavía no es fiable hasta alcanzar EOF.
const F5_PRE_ADMISSION_TOML_LIMIT_BYTES: u64 = 65_536;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ValidatedGeneralSummary {
    pub total_genome_size: u64,
}

/// Tope por omisión del tamaño del tracto, **en unidades completas**
/// (pendiente P-5, confirmado por el responsable el 20/09/2026).
///
/// Un microsatélite de trabajo rara vez pasa de treinta o cuarenta unidades, de
/// modo que doscientas no limitan ningún uso previsto y actúan solo como red
/// contra un error de escritura. Quien quiera simular expansiones patológicas
/// —cientos o miles de unidades— debe subirlo explícitamente.
pub const DEFAULT_REPEAT_MAX_UNITS: u64 = 200;

#[derive(Clone, Debug, Eq, PartialEq)]

pub struct EffectiveSequence {
    id: String,
    class: SequenceClass,
    length: u64,
}

impl EffectiveSequence {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

/// Versión del documento general que dio origen a la configuración efectiva.
///
/// Manda qué versión de configuración efectiva se publica: un documento v2
/// sigue publicando `effective-config-v4` **byte a byte** como antes, y solo el
/// v4, que es el que puede declarar repeticiones, publica v6.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DeclaredDocumentVersion {
    V2,
    V4,
}

impl DeclaredDocumentVersion {
    /// Versión del documento de configuración efectiva que le corresponde.
    #[must_use]
    pub const fn effective_config_schema_version(self) -> &'static str {
        match self {
            Self::V2 => "4",
            Self::V4 => "6",
        }
    }

    #[must_use]
    pub const fn effective_config_schema_file(self) -> &'static str {
        match self {
            Self::V2 => "effective-config-v4.schema.json",
            Self::V4 => "effective-config-v6.schema.json",
        }
    }
}

/// Un alelo de repetición ya resuelto: su **secuencia materializada** y la
/// frecuencia con la que se sortea entre las copias haplotípicas.
///
/// No guarda número de unidades ni pureza: los dos se derivan del tracto, que
/// es la única autoridad (decisión S-02).
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveRepeatAllele {
    tract: TandemRepeatTract,
    declared_frequency: f64,
    frequency: f64,
}

impl EffectiveRepeatAllele {
    #[must_use]
    pub const fn tract(&self) -> &TandemRepeatTract {
        &self.tract
    }

    /// Frecuencia **tal y como el usuario la escribió**, sin normalizar.
    ///
    /// Se conserva y se publica aparte porque, si la suma declarada no era
    /// exactamente uno, el valor normalizado **no es el que el usuario
    /// escribió** y no habría dónde consultarlo.
    #[must_use]
    pub const fn declared_frequency(&self) -> f64 {
        self.declared_frequency
    }

    /// Frecuencia **normalizada**, que es la que consume el sorteo.
    #[must_use]
    pub const fn frequency(&self) -> f64 {
        self.frequency
    }
}

/// Auditoría de la normalización de las frecuencias de un locus.
///
/// DTC-049 no se conforma con normalizar: exige **recomputar** la suma y su
/// desviación después de normalizar y publicar las dos, junto con la tolerancia
/// efectiva y el algoritmo. Es lo que permite comprobar desde fuera que la
/// normalización cerró y con qué margen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RepeatFrequencyNormalization {
    tolerance: f64,
    original_sum: f64,
    original_delta: f64,
    effective_sum: f64,
    effective_delta: f64,
}

impl RepeatFrequencyNormalization {
    #[must_use]
    pub const fn tolerance(self) -> f64 {
        self.tolerance
    }

    #[must_use]
    pub const fn original_sum(self) -> f64 {
        self.original_sum
    }

    #[must_use]
    pub const fn original_delta(self) -> f64 {
        self.original_delta
    }

    #[must_use]
    pub const fn effective_sum(self) -> f64 {
        self.effective_sum
    }

    #[must_use]
    pub const fn effective_delta(self) -> f64 {
        self.effective_delta
    }
}

/// Un locus de repetición ya resuelto.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveRepeatLocus {
    /// El locus del dominio, con su identificador tipado, su intervalo en el
    /// servicio común de coordenadas y su motivo validado.
    ///
    /// La frontera de configuración **no** mantiene una representación paralela
    /// con campos crudos: había dos del mismo concepto y una de ellas, la
    /// tipada, no tenía usuarios.
    declaration: TandemRepeatLocus,
    max_units: u64,
    admission: EffectiveRepeatAdmission,
    alleles: Vec<EffectiveRepeatAllele>,
    normalization: RepeatFrequencyNormalization,
}

impl EffectiveRepeatLocus {
    #[must_use]
    pub const fn declaration(&self) -> &TandemRepeatLocus {
        &self.declaration
    }

    #[must_use]
    pub fn id(&self) -> &str {
        self.declaration.id().as_str()
    }

    #[must_use]
    pub fn sequence_id(&self) -> &str {
        self.declaration.interval().sequence_id().as_str()
    }

    /// Inicio del intervalo en **base uno, inclusivo**, convertido por el
    /// servicio común de coordenadas y no por una suma a mano.
    pub fn start(&self) -> Result<u64, Diagnostic> {
        Ok(self.closed_interval()?.start())
    }

    /// Final del intervalo en **base uno, inclusivo**.
    pub fn end(&self) -> Result<u64, Diagnostic> {
        Ok(self.closed_interval()?.end())
    }

    fn closed_interval(&self) -> Result<ClosedInterval1, Diagnostic> {
        self.declaration.interval().to_closed_1().map_err(|error| {
            Diagnostic::internal_invariant("application", "repeat_locus_interval")
                .with("detail", format!("{error:?}"))
        })
    }

    /// Los motivos del locus, en orden.
    #[must_use]
    pub fn motifs(&self) -> &[RepeatMotif] {
        self.declaration.motifs()
    }

    /// El motivo cuando el locus tiene uno solo.
    #[must_use]
    pub fn single_motif(&self) -> Option<&RepeatMotif> {
        self.declaration.single_motif()
    }

    /// Mínimos de admisión **efectivos de este locus**: los suyos si los declara
    /// y, si no, los generales.
    #[must_use]
    pub const fn admission(&self) -> EffectiveRepeatAdmission {
        self.admission
    }

    /// Tope propio del tamaño del tracto, en **unidades completas**.
    #[must_use]
    pub const fn max_units(&self) -> u64 {
        self.max_units
    }

    #[must_use]
    pub fn alleles(&self) -> &[EffectiveRepeatAllele] {
        &self.alleles
    }

    #[must_use]
    pub const fn normalization(&self) -> RepeatFrequencyNormalization {
        self.normalization
    }
}

/// Mínimos de la guarda de admisión con referencia real, ya resueltos.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffectiveRepeatAdmission {
    min_purity: f64,
    min_longest_run_units: u64,
}

impl EffectiveRepeatAdmission {
    /// Valores por omisión: se publican siempre, se declaren o no.
    pub const DEFAULT: Self = Self {
        min_purity: 0.8,
        min_longest_run_units: 3,
    };

    #[must_use]
    pub const fn min_purity(self) -> f64 {
        self.min_purity
    }

    #[must_use]
    pub const fn min_longest_run_units(self) -> u64 {
        self.min_longest_run_units
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveComposition {
    values: [f64; 4],
    unknown_base_fraction: f64,
}

/// Parámetro efectivo compartido por las normalizaciones DTC-049 consumidas en F2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffectiveNormalization {
    tolerance: f64,
}

impl EffectiveNormalization {
    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        self.tolerance
    }
}

impl EffectiveComposition {
    #[must_use]
    pub const fn values(&self) -> &[f64; 4] {
        &self.values
    }

    /// Fracción independiente de posiciones ambiguas (`N`) de la referencia generada.
    #[must_use]
    pub const fn unknown_base_fraction(&self) -> f64 {
        self.unknown_base_fraction
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReferenceSource {
    Simulated,
    Supplied { original_path: PathBuf },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveReductionTarget {
    Bases { declared: u64 },
    ReductionPercent { literal: String },
}

impl EffectiveReductionTarget {
    #[must_use]
    pub const fn mode(&self) -> &'static str {
        match self {
            Self::Bases { .. } => "target_bases",
            Self::ReductionPercent { .. } => "reduction_percent",
        }
    }

    #[must_use]
    pub fn declared_value(&self) -> String {
        match self {
            Self::Bases { declared } => declared.to_string(),
            Self::ReductionPercent { literal } => literal.clone(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveReducedReference {
    id: String,
    source_reference_id: String,
    region_size_bases: u64,
    target: EffectiveReductionTarget,
    target_bases: u64,
    selection_subseed: Subseed,
}

impl EffectiveReducedReference {
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub fn source_reference_id(&self) -> &str {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn region_size_bases(&self) -> u64 {
        self.region_size_bases
    }

    #[must_use]
    pub const fn target(&self) -> &EffectiveReductionTarget {
        &self.target
    }

    #[must_use]
    pub const fn target_bases(&self) -> u64 {
        self.target_bases
    }

    #[must_use]
    pub const fn selection_subseed(&self) -> Subseed {
        self.selection_subseed
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveConfigurationInput {
    original_path: String,
    bytes: Arc<[u8]>,
}

impl EffectiveConfigurationInput {
    #[must_use]
    pub fn original_path(&self) -> &str {
        &self.original_path
    }

    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveConfig {
    source_path: String,
    reference_id: String,
    ploidy: Ploidy,
    normalization: EffectiveNormalization,
    reference_source: EffectiveReferenceSource,
    composition: Option<EffectiveComposition>,
    reduced_references: Vec<EffectiveReducedReference>,
    total_genome_size: u64,
    sequences: Vec<EffectiveSequence>,
    configuration_input: EffectiveConfigurationInput,
    reproducibility: ReproducibilityMetadata,
    population_snvs: Option<EffectivePopulationSnvs>,
    initial_populations: Vec<EffectiveInitialPopulation>,
    direct_events: Vec<EffectiveDirectEvent>,
    technical_event: Option<EffectiveTechnicalEvent>,
    vcf_output: Option<VcfOutputMode>,
    document_version: DeclaredDocumentVersion,
    repeat_loci: Vec<EffectiveRepeatLocus>,
    repeat_admission: EffectiveRepeatAdmission,
}

/// Configuración F5 ya resuelta entre sus dos documentos de autoría y el
/// paquete de perfil F5. La máscara y la referencia pertenecen a F5.0.4.
#[derive(Clone, Debug)]
pub struct EffectiveF5Config {
    general_path: PathBuf,
    scenario_path: PathBuf,
    profile_metadata_path: PathBuf,
    excluded_regions_path: Option<PathBuf>,
    reference_path: PathBuf,
    seed: String,
    seed_origin: String,
    ploidy: u8,
    reference_id: String,
    reference_sha256: String,
    scenario_id: String,
    profile: F5Profile,
    bootstrap_documents: F5BootstrapDocuments,
    input_documents: F5InputDocuments,
    profile_input_identities: Vec<AdmittedInputIdentity>,
    /// Estado interno de la admisión v9 después de bootstrap y perfil, que la
    /// primera visita FASTA/máscara continúa. No es contabilidad publicable
    /// mientras la fórmula final no lo selle.
    admission_ledger: F5ResourceLedger,
    /// Admisión sellada: la única pasada FASTA/máscara, sus aportaciones
    /// externas y el receipt que une identidades y fases. Es lo único que abre
    /// layout, y no recalcula ni duplica la geometría retenida.
    admitted: Option<crate::f5_preflight::SealedF5Admission>,
    admitted_ledger: Option<F5ResourceLedger>,
    general: DeclaredF5GeneralConfig,
    scenario: DeclaredF5Scenario,
}

/// Observación física de uno de los dos TOML que forman el bootstrap F5.
///
/// El contador y SHA-256 se calculan mientras se lee el mismo handle y se
/// conservan ya como `PhysicalInputIdentity`. Aún no es un recibo de admisión.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5BootstrapDocument {
    identity: PhysicalInputIdentity,
}

impl F5BootstrapDocument {
    #[must_use]
    pub const fn size_bytes(&self) -> u64 {
        self.identity.observed_size_bytes()
    }

    #[must_use]
    pub fn sha256(&self) -> &str {
        self.identity.observed_sha256()
    }

    #[must_use]
    pub const fn identity(&self) -> &PhysicalInputIdentity {
        &self.identity
    }
}

/// Los dos documentos ya leídos, medidos y hasheados por el bootstrap.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5BootstrapDocuments {
    general: F5BootstrapDocument,
    scenario: F5BootstrapDocument,
}

/// Resultado sellado del primer estado de admisión F5. No ha abierto ningún
/// recurso posterior al escenario ni ha resuelto una semilla.
#[derive(Clone, Debug)]
pub struct F5BootstrapConfig {
    general_path: PathBuf,
    scenario_path: PathBuf,
    profile_metadata_path: PathBuf,
    excluded_regions_path: Option<PathBuf>,
    reference_path: PathBuf,
    documents: F5BootstrapDocuments,
    general: DeclaredF5GeneralConfig,
    scenario: DeclaredF5Scenario,
}

impl F5BootstrapConfig {
    #[must_use]
    pub fn documents(&self) -> &F5BootstrapDocuments {
        &self.documents
    }
}

impl F5BootstrapDocuments {
    #[must_use]
    pub const fn general(&self) -> &F5BootstrapDocument {
        &self.general
    }

    #[must_use]
    pub const fn scenario(&self) -> &F5BootstrapDocument {
        &self.scenario
    }
}

/// Tamaños físicos de todos los documentos que el lector F5 retuvo y parseó
/// antes de iniciar el preflight geométrico. Se derivan exclusivamente de las
/// seis identidades observadas, tras la admisión secuencial del metadata y sus
/// cinco TSV.
///
/// No contienen identidad semántica: sus bytes de relleno no cambian el RNG ni
/// el layout, pero sí forman parte de la reserva transitoria de admisión.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct F5InputDocuments {
    general_toml_bytes: u64,
    scenario_toml_bytes: u64,
    profile_metadata_bytes: u64,
    profile_tsv_bytes: u64,
    profile_tsv_peak_bytes: u64,
}

impl F5InputDocuments {
    fn from_observed(
        general_toml_bytes: u64,
        scenario_toml_bytes: u64,
        identities: &[AdmittedInputIdentity],
    ) -> Result<Self, Diagnostic> {
        let (metadata, tables) = identities.split_first().ok_or_else(|| {
            Diagnostic::f5_resource("<unresolved>", "faltan identidades físicas del perfil")
        })?;
        let expected_table_roles = [
            InputRole::ProfileBlockLengths,
            InputRole::ProfileSnpGaps,
            InputRole::ProfileSnpSpectrum,
            InputRole::ProfileHaplotypeWeights,
            InputRole::ProfileLd,
        ];
        if metadata.role() != InputRole::ProfileMetadata
            || tables.len() != expected_table_roles.len()
            || tables
                .iter()
                .zip(expected_table_roles)
                .any(|(identity, role)| identity.role() != role)
        {
            return Err(Diagnostic::f5_resource(
                "<unresolved>",
                "las identidades físicas del perfil no tienen el orden contractual",
            ));
        }
        let profile_metadata_bytes = metadata.observed_size_bytes();
        let (profile_tsv_bytes, profile_tsv_peak_bytes) =
            tables
                .iter()
                .try_fold((0_u64, 0_u64), |(total, peak), identity| {
                    total
                        .checked_add(identity.observed_size_bytes())
                        .map(|total| (total, peak.max(identity.observed_size_bytes())))
                        .ok_or_else(|| {
                            Diagnostic::f5_resource(
                                "<unresolved>",
                                "overflow al sumar los bytes observados del paquete de perfil",
                            )
                        })
                })?;
        Ok(Self {
            general_toml_bytes,
            scenario_toml_bytes,
            profile_metadata_bytes,
            profile_tsv_bytes,
            profile_tsv_peak_bytes,
        })
    }

    #[must_use]
    pub const fn general_toml_bytes(self) -> u64 {
        self.general_toml_bytes
    }
    #[must_use]
    pub const fn scenario_toml_bytes(self) -> u64 {
        self.scenario_toml_bytes
    }
    #[must_use]
    pub const fn profile_metadata_bytes(self) -> u64 {
        self.profile_metadata_bytes
    }
    #[must_use]
    pub const fn profile_tsv_bytes(self) -> u64 {
        self.profile_tsv_bytes
    }

    #[must_use]
    pub const fn profile_tsv_peak_bytes(self) -> u64 {
        self.profile_tsv_peak_bytes
    }

    /// Cada parser conserva lectura y representación temporal de su propia
    /// entrada. Los TSV son fases secuenciales, por lo que solo contribuye la
    /// mayor tabla, no la suma de cinco búferes que no coexisten.
    pub fn parser_peak_bytes(self, scenario_id: &str) -> Result<u64, Diagnostic> {
        [
            self.general_toml_bytes
                .checked_add(self.scenario_toml_bytes)
                .and_then(|total| total.checked_mul(2)),
            self.profile_metadata_bytes.checked_mul(2),
            self.profile_tsv_peak_bytes.checked_mul(2),
        ]
        .into_iter()
        .try_fold(0_u64, |peak, value| value.map(|value| peak.max(value)))
        .ok_or_else(|| {
            Diagnostic::f5_resource(
                scenario_id,
                "overflow en la reserva de lectura, retención y parser de documentos F5",
            )
        })
    }
}

#[derive(Debug)]
struct BoundedF5Toml {
    bytes: Vec<u8>,
    document: F5BootstrapDocument,
}

/// Lee exactamente un documento TOML de bootstrap desde un único handle.
///
/// Tras retener como máximo 65.536 bytes, se lee un único byte adicional para
/// distinguir EOF del exceso. Ese byte no entra ni en el búfer ni en el hash.
fn read_bounded_f5_toml(path: &Path) -> Result<BoundedF5Toml, Diagnostic> {
    let file = File::open(path).map_err(|error| input_error(path, error))?;
    read_bounded_f5_toml_from_reader(path, InputRole::GeneralToml, file)
}

fn read_bounded_f5_toml_with_opener(
    path: &Path,
    role: InputRole,
    opener: &mut impl FnMut(&Path) -> std::io::Result<File>,
) -> Result<BoundedF5Toml, Diagnostic> {
    let file = opener(path).map_err(|error| input_error(path, error))?;
    read_bounded_f5_toml_from_reader(path, role, file)
}

/// Variante inyectable del lector acotado. Mantiene la lectura `límite + 1`
/// separada del parser para que la cota se pueda demostrar sin materializar el
/// resto de la entrada.
fn read_bounded_f5_toml_from_reader(
    path: &Path,
    role: InputRole,
    reader: impl std::io::Read,
) -> Result<BoundedF5Toml, Diagnostic> {
    let bounded = BoundedInput::new(
        role,
        AdmissionPhase::Bootstrap,
        path.to_path_buf(),
        F5_PRE_ADMISSION_TOML_LIMIT_BYTES,
        reader,
    )
    .read_to_end()
    .map_err(|diagnostic| {
        // Antes de EOF no existe una versión ni un escenario fiables. Se
        // conserva, por tanto, el diagnóstico neutral de F5.1.2.1 aunque el
        // lector común ya haya contado el byte de exceso físicamente.
        if diagnostic.code() == "GPS-RES-101"
            && diagnostic.context().get("dimension").is_some_and(|dimension| dimension == "physical_input_bytes")
        {
            Diagnostic::input(
                "GPS-INPUT-001",
                "No se puede abrir o leer una entrada declarada.",
            )
            .with("path", path.display().to_string())
            .with(
                "detail",
                format!(
                    "toml_bytes>{F5_PRE_ADMISSION_TOML_LIMIT_BYTES}/{F5_PRE_ADMISSION_TOML_LIMIT_BYTES}; el byte 65537 confirma que el documento excede el límite físico antes de clasificarlo"
                ),
            )
        } else {
            diagnostic
        }
    })?;
    let (bytes, identity) = bounded.into_parts();
    Ok(BoundedF5Toml {
        bytes,
        document: F5BootstrapDocument { identity },
    })
}

/// Coste mínimo y temporal del bootstrap, contrastado antes de cualquier
/// apertura posterior.
fn admit_f5_bootstrap(
    general: &F5BootstrapDocument,
    scenario: &F5BootstrapDocument,
    declared: &DeclaredF5Scenario,
) -> Result<(), Diagnostic> {
    let document_bytes = general
        .size_bytes()
        .checked_add(scenario.size_bytes())
        .ok_or_else(|| {
            Diagnostic::f5_resource(
                declared.scenario_id.value(),
                "overflow al sumar los documentos del bootstrap F5",
            )
        })?;
    // Cada byte se conserva como documento admitido y como representación
    // transitoria del parser. No se presupuestan perfil ni geometría aquí.
    let bootstrap_bytes = document_bytes.checked_mul(2).ok_or_else(|| {
        Diagnostic::f5_resource(
            declared.scenario_id.value(),
            "overflow al calcular la memoria transitoria del bootstrap F5",
        )
    })?;
    let working_limit = declared
        .max_working_memory_bytes
        .as_ref()
        .map_or(F5_HARD_WORKING_MEMORY_LIMIT, |limit| *limit.value());
    if bootstrap_bytes > working_limit {
        return Err(Diagnostic::f5_resource(
            declared.scenario_id.value(),
            format!(
                "working_memory_bytes={bootstrap_bytes}/{working_limit}, phase=bootstrap; general y escenario se rechazan antes de abrir perfil, referencia, máscara, destino o RNG"
            ),
        ));
    }
    Ok(())
}

/// Resultado de la frontera CLI.
///
/// La rama F5 lleva la configuración **ya admitida**: el dispatcher completa
/// la admisión física .4 una sola vez y `run` la reutiliza, igual que hace la
/// rama histórica con su `EffectiveConfig`. Repetirla abriría referencia,
/// perfil y máscara por segunda vez, ensancharía la ventana TOCTOU y
/// desmentiría la contabilidad publicada, que sella **una** pasada.
#[derive(Clone, Debug)]
pub enum CliValidatedConfig {
    Legacy(Box<EffectiveConfig>),
    F5(Box<EffectiveF5Config>),
}

impl EffectiveF5Config {
    #[must_use]
    pub fn general_path(&self) -> &Path {
        &self.general_path
    }
    #[must_use]
    pub fn scenario_path(&self) -> &Path {
        &self.scenario_path
    }
    #[must_use]
    pub fn profile_metadata_path(&self) -> &Path {
        &self.profile_metadata_path
    }
    #[must_use]
    pub fn excluded_regions_path(&self) -> Option<&Path> {
        self.excluded_regions_path.as_deref()
    }
    #[must_use]
    pub fn reference_path(&self) -> &Path {
        &self.reference_path
    }
    #[must_use]
    pub fn seed(&self) -> &str {
        &self.seed
    }
    #[must_use]
    pub fn seed_origin(&self) -> &str {
        &self.seed_origin
    }
    #[must_use]
    pub const fn ploidy(&self) -> u8 {
        self.ploidy
    }
    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }
    #[must_use]
    pub fn reference_sha256(&self) -> &str {
        &self.reference_sha256
    }
    #[must_use]
    pub fn scenario_id(&self) -> &str {
        &self.scenario_id
    }
    #[must_use]
    pub const fn profile(&self) -> &F5Profile {
        &self.profile
    }
    #[must_use]
    pub const fn input_documents(&self) -> F5InputDocuments {
        self.input_documents
    }
    #[must_use]
    pub fn profile_input_identities(&self) -> &[AdmittedInputIdentity] {
        &self.profile_input_identities
    }
    /// Las ocho identidades físicas que esta configuración ya observó: los dos
    /// TOML de autoría, el metadata y los cinco TSV, en orden de rol.
    ///
    /// Referencia y máscara no aparecen aquí porque pertenecen al escaneo, no
    /// a la carga de documentos.
    pub(crate) fn admitted_document_identities(&self) -> Vec<AdmittedInputIdentity> {
        let mut identities = Vec::with_capacity(2 + self.profile_input_identities.len());
        identities.push(AdmittedInputIdentity::undeclared(
            self.bootstrap_documents.general().identity().clone(),
        ));
        identities.push(AdmittedInputIdentity::undeclared(
            self.bootstrap_documents.scenario().identity().clone(),
        ));
        identities.extend_from_slice(&self.profile_input_identities);
        identities
    }
    #[must_use]
    pub(crate) fn admission_ledger(&self) -> F5ResourceLedger {
        self.admission_ledger.clone()
    }
    #[must_use]
    pub(crate) fn admitted(&self) -> Option<&crate::f5_preflight::SealedF5Admission> {
        self.admitted.as_ref()
    }
    #[must_use]
    pub const fn bootstrap_documents(&self) -> &F5BootstrapDocuments {
        &self.bootstrap_documents
    }
    #[must_use]
    pub const fn declared_general(&self) -> &DeclaredF5GeneralConfig {
        &self.general
    }
    #[must_use]
    pub const fn declared_scenario(&self) -> &DeclaredF5Scenario {
        &self.scenario
    }
}

/// Perfil sintético único que F3 genera antes de que F4 incorpore perfiles
/// cuantitativos o tecnologías alternativas.
pub const ILLUMINA_SE_MINIMAL_V1: &str = "illumina-se-minimal-v1";
pub const ILLUMINA_SE_MINIMAL_READ_LENGTH: u16 = 150;
pub const ILLUMINA_SE_MINIMAL_ACGT_PHRED: u8 = 40;
pub const ILLUMINA_SE_MINIMAL_N_PHRED: u8 = 0;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EffectiveTechnicalEventTarget {
    CoverageWgs(f64),
    ReadCount(u64),
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveTechnicalEvent {
    event_id: EventId,
    individual_ids: Vec<IndividualId>,
    target: EffectiveTechnicalEventTarget,
    sequencing_errors: SequencingErrorsMode,
    read_length: u16,
    chimera_rate: f64,
    index_hopping_rate: f64,
    biological_contamination: BiologicalContaminationSettings,
    libraries: Vec<EffectiveTechnicalLibrary>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveTechnicalLibrary {
    individual_id: IndividualId,
    i7: String,
    i5: String,
}

impl EffectiveTechnicalLibrary {
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }
    #[must_use]
    pub fn i7(&self) -> &str {
        &self.i7
    }
    #[must_use]
    pub fn i5(&self) -> &str {
        &self.i5
    }
}

impl EffectiveTechnicalEvent {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub fn individual_ids(&self) -> &[IndividualId] {
        &self.individual_ids
    }

    #[must_use]
    pub const fn target(&self) -> EffectiveTechnicalEventTarget {
        self.target
    }

    #[must_use]
    pub const fn sequencing_errors(&self) -> SequencingErrorsMode {
        self.sequencing_errors
    }

    #[must_use]
    pub fn libraries(&self) -> &[EffectiveTechnicalLibrary] {
        &self.libraries
    }

    #[must_use]
    pub const fn profile(&self) -> &'static str {
        GPS_HISEQX_PCRFREE_PE150_V1
    }

    #[must_use]
    pub const fn read_length(&self) -> u16 {
        self.read_length
    }

    /// Tasa Bernoulli ya validada y materializada; por omisión es el neutro
    /// del perfil (`0`).
    #[must_use]
    pub const fn chimera_rate(&self) -> f64 {
        self.chimera_rate
    }

    /// Tasa Bernoulli posterior a contaminación biológica; por omisión usa el
    /// valor neutro del perfil (`0`).
    #[must_use]
    pub const fn index_hopping_rate(&self) -> f64 {
        self.index_hopping_rate
    }

    #[must_use]
    pub const fn biological_contamination(&self) -> &BiologicalContaminationSettings {
        &self.biological_contamination
    }

    /// Ajuste tipado que se entrega al muestreador PE sin volver a interpretar
    /// la configuración.
    #[must_use]
    pub fn chimera_settings(&self) -> ChimeraSettings {
        ChimeraSettings::new(self.chimera_rate).expect("effective chimera rate is validated")
    }

    #[must_use]
    pub const fn acgt_phred(&self) -> u8 {
        ILLUMINA_SE_MINIMAL_ACGT_PHRED
    }

    #[must_use]
    pub const fn n_phred(&self) -> u8 {
        ILLUMINA_SE_MINIMAL_N_PHRED
    }

    /// El perfil inicial no usa amplificación; por tanto no existe una tasa
    /// PCR efectiva que pueda configurar o materializar F4.
    #[must_use]
    pub const fn pcr_duplicates_applicable(&self) -> bool {
        false
    }

    /// El perfil inicial tampoco modela duplicados ópticos. El recuento es un
    /// hecho efectivo, no una tasa de entrada ni el resultado de un sorteo.
    #[must_use]
    pub const fn realized_optical_duplicate_pairs(&self) -> u64 {
        0
    }
}

/// Origen efectivo del recurso canónico F2.2.7.
#[derive(Clone, Debug, PartialEq)]
pub enum EffectivePopulationSnvs {
    Supplied {
        input_path: PathBuf,
    },
    Reuse {
        metadata_path: PathBuf,
    },
    Generate {
        parameters: PopulationSnvsGenerationConfig,
    },
}

/// Modelo de genotipo solicitado para una cohorte directa.
///
/// F2.3.1 solo conserva esta elección validada. El sorteo de composición y la
/// fase se incorporan en F2.3.2 y F2.3.3, respectivamente.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EffectiveInitialGenotypeModel {
    HardyWeinberg,
    FDiploid { inbreeding_coefficient: f64 },
}

/// Cambio directo de F2.3.4 ya ligado a una identidad individual y a una
/// fase existente. Las coordenadas de indels mantienen el convenio actual
/// del compositor, pues se evalúan en orden declarado tras los SNP de matriz.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum EffectiveDirectEvent {
    Snp {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        position_1based: u64,
        reference: u8,
        alternate: u8,
    },
    Insertion {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        point: u64,
        bases: Arc<[u8]>,
    },
    Deletion {
        event_id: EventId,
        individual_id: IndividualId,
        phase: HaplotypePhase,
        sequence_id: SequenceId,
        start: u64,
        end: u64,
    },
}

impl EffectiveDirectEvent {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        match self {
            Self::Snp { event_id, .. }
            | Self::Insertion { event_id, .. }
            | Self::Deletion { event_id, .. } => event_id,
        }
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        match self {
            Self::Snp { individual_id, .. }
            | Self::Insertion { individual_id, .. }
            | Self::Deletion { individual_id, .. } => individual_id,
        }
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        match self {
            Self::Snp { phase, .. }
            | Self::Insertion { phase, .. }
            | Self::Deletion { phase, .. } => *phase,
        }
    }

    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        match self {
            Self::Snp { sequence_id, .. }
            | Self::Insertion { sequence_id, .. }
            | Self::Deletion { sequence_id, .. } => sequence_id,
        }
    }
}

/// Identidades y parámetros ya resueltos de una cohorte directa.
///
/// La población de origen es un perfil de frecuencias de la matriz, no una
/// población evolutiva ni una estructura reproductiva.
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveInitialPopulation {
    population_id: PopulationId,
    individual_ids: Vec<IndividualId>,
    founder_keys: Option<Vec<FounderKey>>,
    ploidy: Ploidy,
    genotype_model: EffectiveInitialGenotypeModel,
}

impl EffectiveInitialPopulation {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub fn individual_ids(&self) -> &[IndividualId] {
        &self.individual_ids
    }

    #[must_use]
    pub fn founder_keys(&self) -> Option<&[FounderKey]> {
        self.founder_keys.as_deref()
    }

    #[must_use]
    pub const fn individual_count(&self) -> usize {
        self.individual_ids.len()
    }

    #[must_use]
    pub const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    #[must_use]
    pub const fn genotype_model(&self) -> EffectiveInitialGenotypeModel {
        self.genotype_model
    }
}

impl EffectiveConfig {
    #[must_use]
    pub fn source_path(&self) -> &str {
        &self.source_path
    }

    /// Versión del documento general que la originó, que es la que decide qué
    /// versión de configuración efectiva se publica.
    #[must_use]
    pub const fn document_version(&self) -> DeclaredDocumentVersion {
        self.document_version
    }

    #[must_use]
    pub fn repeat_loci(&self) -> &[EffectiveRepeatLocus] {
        &self.repeat_loci
    }

    #[must_use]
    pub const fn repeat_admission(&self) -> EffectiveRepeatAdmission {
        self.repeat_admission
    }

    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }

    #[must_use]
    pub const fn ploidy(&self) -> Ploidy {
        self.ploidy
    }

    #[must_use]
    pub const fn normalization(&self) -> EffectiveNormalization {
        self.normalization
    }

    #[must_use]
    pub const fn reference_source(&self) -> &EffectiveReferenceSource {
        &self.reference_source
    }

    #[must_use]
    pub const fn composition(&self) -> Option<&EffectiveComposition> {
        self.composition.as_ref()
    }

    #[must_use]
    pub fn reduced_references(&self) -> &[EffectiveReducedReference] {
        &self.reduced_references
    }

    #[must_use]
    pub const fn total_genome_size(&self) -> u64 {
        self.total_genome_size
    }

    #[must_use]
    pub fn sequences(&self) -> &[EffectiveSequence] {
        &self.sequences
    }

    #[must_use]
    pub const fn configuration_input(&self) -> &EffectiveConfigurationInput {
        &self.configuration_input
    }

    #[must_use]
    pub const fn reproducibility(&self) -> &ReproducibilityMetadata {
        &self.reproducibility
    }
    #[must_use]
    pub fn population_snvs_generation(&self) -> Option<&PopulationSnvsGenerationConfig> {
        match self.population_snvs.as_ref() {
            Some(EffectivePopulationSnvs::Generate { parameters }) => Some(parameters),
            _ => None,
        }
    }

    #[must_use]
    pub fn population_snvs(&self) -> Option<&EffectivePopulationSnvs> {
        self.population_snvs.as_ref()
    }

    /// Cohortes directas listas para que F2.3.2 consuma sus identidades.
    #[must_use]
    pub fn initial_populations(&self) -> &[EffectiveInitialPopulation] {
        &self.initial_populations
    }

    /// Eventos directos tipados y en su orden de declaración global.
    #[must_use]
    pub fn direct_events(&self) -> &[EffectiveDirectEvent] {
        &self.direct_events
    }

    #[must_use]
    pub const fn technical_event(&self) -> Option<&EffectiveTechnicalEvent> {
        self.technical_event.as_ref()
    }

    #[must_use]
    pub const fn vcf_output(&self) -> Option<VcfOutputMode> {
        self.vcf_output
    }
}

pub fn load_and_validate_config(path: &Path) -> Result<EffectiveConfig, Diagnostic> {
    load_and_validate_config_with_entropy(path, &mut SystemEntropy)
}

pub fn load_and_validate_config_with_entropy(
    path: &Path,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    // La API histórica sigue devolviendo exclusivamente el general v2, pero
    // antes de EOF cualquier TOML puede ser el general v3 de F5. Comparte por
    // ello la frontera física de 65.536+1 bytes con el dispatcher.
    let BoundedF5Toml { bytes, .. } = read_bounded_f5_toml(path)?;
    let admitted = admit_toml_header(path, bytes, ExpectedConfigKind::General)?;
    let general = parse_declared_config(admitted)?;
    build_effective_config_with_entropy(&general, entropy)
}

/// Despacha general v2 al flujo histórico y general v3 al bootstrap F5.
///
/// La apertura física es única: primero se lee el documento completo bajo la
/// cota común y solo entonces se clasifica su versión. El flujo histórico no
/// cambia una vez admitido un general v2.
pub fn load_and_validate_cli_config(path: &Path) -> Result<CliValidatedConfig, Diagnostic> {
    load_and_validate_cli_config_with_toml_opener(path, &mut |opened_path: &Path| {
        File::open(opened_path)
    })
}

/// Variante instrumental del dispatcher para demostrar su frontera de apertura.
///
/// No altera el contrato de producción: el observador se invoca justo antes de
/// cada apertura física de TOML del bootstrap.
#[doc(hidden)]
pub fn load_and_validate_cli_config_with_bootstrap_open_observer(
    path: &Path,
    observer: &mut impl FnMut(&Path),
) -> Result<CliValidatedConfig, Diagnostic> {
    load_and_validate_cli_config_with_toml_opener(path, &mut |opened_path| {
        observer(opened_path);
        File::open(opened_path)
    })
}

fn load_and_validate_cli_config_with_toml_opener(
    path: &Path,
    opener: &mut impl FnMut(&Path) -> std::io::Result<File>,
) -> Result<CliValidatedConfig, Diagnostic> {
    let general_path = absolute_normalized_path(path)?;
    let BoundedF5Toml { bytes, document } =
        read_bounded_f5_toml_with_opener(&general_path, InputRole::GeneralToml, opener)?;
    // `general_path` es solo la ruta física para abrir y resolver enlaces F5.
    // La ruta de autoría que se publica en el flujo v2 es exactamente la que
    // solicitó el llamador.
    let admitted = admit_toml_header(path, bytes, ExpectedConfigKind::GeneralDispatch)?;
    if admitted_is_f5_general_v3(&admitted) {
        // El nombre del general se publica en `effective-config.json`: la
        // regla vive aquí, de modo que `validate` y `run` la aplican igual y
        // ningún otro lenguaje necesita repetirla (DTC-051).
        require_publishable_general_name(&general_path)?;
        let bootstrap = bootstrap_f5_config_from_admitted_with_toml_opener(
            general_path,
            document,
            admitted,
            opener,
        )?;
        // La admisión .4 corre **una** vez y su resultado viaja: `validate` la
        // usa y la descarta; `run` la reutiliza sin reabrir nada.
        let effective = materialize_f5_config_after_bootstrap(
            bootstrap,
            &mut SystemEntropy,
            crate::f5_preflight::F5OtherResources::default(),
        )?;
        return Ok(CliValidatedConfig::F5(Box::new(effective)));
    }
    let general = parse_declared_config(admitted)?;
    build_effective_config_with_entropy(&general, &mut SystemEntropy)
        .map(Box::new)
        .map(CliValidatedConfig::Legacy)
}

/// Construye la única transición bootstrap→perfil usada por la CLI y la
/// biblioteca. El pico del bootstrap queda en el histórico, mientras sus
/// búferes se declaran liberados (`retained_after=0`) antes del metadata.
fn profile_ledger_after_bootstrap(
    scenario: &DeclaredF5Scenario,
    documents: &F5BootstrapDocuments,
) -> Result<F5ResourceLedger, Diagnostic> {
    // Los tres predeterminados son exactamente los tres topes duros: no se
    // reescriben aquí como literales porque una copia silenciosa es la que
    // sobrevive a un cambio de contrato y deja el bootstrap con el techo
    // viejo mientras el resto del flujo usa el nuevo.
    let operation_limit = scenario
        .max_operation_units
        .as_ref()
        .map_or(F5_HARD_OPERATION_LIMIT, |limit| *limit.value());
    let working_memory_limit = scenario
        .max_working_memory_bytes
        .as_ref()
        .map_or(F5_HARD_WORKING_MEMORY_LIMIT, |limit| *limit.value());
    let output_limit = scenario
        .max_output_bytes
        .as_ref()
        .map_or(F5_HARD_OUTPUT_LIMIT, |limit| *limit.value());
    let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
        OperationCost::new(operation_limit),
        WorkingMemoryPeak::new(working_memory_limit),
        OutputReservation::new(output_limit),
    ));
    let bootstrap_bytes = documents
        .general()
        .size_bytes()
        .checked_add(documents.scenario().size_bytes())
        .ok_or_else(|| {
            Diagnostic::f5_resource(scenario.scenario_id.value(), "overflow de bootstrap")
        })?;
    let bootstrap_cost = bootstrap_bytes.checked_mul(2).ok_or_else(|| {
        Diagnostic::f5_resource(scenario.scenario_id.value(), "overflow de coste bootstrap")
    })?;
    ledger.charge(
        InputRole::GeneralToml,
        LedgerCharge::Operation(OperationCost::new(bootstrap_cost)),
    )?;
    ledger.charge_memory_phase(InputRole::GeneralToml, 0, bootstrap_cost, 0)?;
    ledger.transition_to(AdmissionPhase::ProfileMetadata)?;
    Ok(ledger)
}

/// Carga el par F5 desde la ruta de CLI. La ruta de CLI se ancla una sola vez
/// al directorio actual; el escenario se ancla al directorio del general.
pub fn load_and_validate_f5_config(path: &Path) -> Result<EffectiveF5Config, Diagnostic> {
    load_and_validate_f5_config_with_entropy(path, &mut SystemEntropy)
}

pub fn load_and_validate_f5_config_with_entropy(
    path: &Path,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveF5Config, Diagnostic> {
    load_and_validate_f5_config_with_entropy_and_resources(
        path,
        entropy,
        crate::f5_preflight::F5OtherResources::default(),
    )
}

/// Variante instrumental que incorpora contribuciones externas antes de abrir
/// referencia o máscara. Producción usa actualmente el valor nulo.
#[doc(hidden)]
pub fn load_and_validate_f5_config_with_entropy_and_resources(
    path: &Path,
    entropy: &mut impl EntropySource,
    other: crate::f5_preflight::F5OtherResources,
) -> Result<EffectiveF5Config, Diagnostic> {
    let general_path = absolute_normalized_path(path)?;
    require_publishable_general_name(&general_path)?;
    let BoundedF5Toml { bytes, document } = read_bounded_f5_toml(&general_path)?;
    let admitted_general = admit_toml_header(path, bytes, ExpectedConfigKind::F5General)
        .map_err(Diagnostic::into_f5_configuration)?;
    let bootstrap = bootstrap_f5_config_from_admitted(general_path, document, admitted_general)?;
    materialize_f5_config_after_bootstrap(bootstrap, entropy, other)
}

/// Continúa la admisión F5 a partir del general ya leído por el dispatcher.
/// No vuelve a abrirlo ni vuelve a calcular su tamaño o hash en otro handle.
fn bootstrap_f5_config_from_admitted(
    general_path: PathBuf,
    general_document: F5BootstrapDocument,
    admitted_general: crate::configuration::AdmittedToml,
) -> Result<F5BootstrapConfig, Diagnostic> {
    bootstrap_f5_config_from_admitted_with_toml_opener(
        general_path,
        general_document,
        admitted_general,
        &mut |opened_path: &Path| File::open(opened_path),
    )
}

fn bootstrap_f5_config_from_admitted_with_toml_opener(
    general_path: PathBuf,
    general_document: F5BootstrapDocument,
    admitted_general: crate::configuration::AdmittedToml,
    opener: &mut impl FnMut(&Path) -> std::io::Result<File>,
) -> Result<F5BootstrapConfig, Diagnostic> {
    let general =
        parse_declared_f5_general(admitted_general).map_err(Diagnostic::into_f5_configuration)?;
    let parent = general_path.parent().ok_or_else(|| {
        config_error_f5(
            &general.source,
            "path",
            "el general no tiene directorio padre",
        )
    })?;
    let scenario_path = lexical_join(parent, Path::new(general.scenario_path.value()))?;
    let reference_path = lexical_join(parent, Path::new(general.reference_path.value()))?;
    let BoundedF5Toml {
        bytes: scenario_bytes,
        document: scenario_document,
    } = read_bounded_f5_toml_with_opener(&scenario_path, InputRole::ScenarioToml, opener)?;
    let admitted_scenario = admit_toml_header(
        &scenario_path,
        scenario_bytes,
        ExpectedConfigKind::F5Scenario,
    )
    .map_err(Diagnostic::into_f5_configuration)?;
    let scenario =
        parse_declared_f5_scenario(admitted_scenario).map_err(Diagnostic::into_f5_configuration)?;
    validate_declared_f5_pair(&general, &scenario).map_err(Diagnostic::into_f5_configuration)?;
    admit_f5_bootstrap(&general_document, &scenario_document, &scenario)?;
    let scenario_parent = scenario_path.parent().ok_or_else(|| {
        config_error_f5(
            &scenario.source,
            "path",
            "el escenario no tiene directorio padre",
        )
    })?;
    let profile_metadata_path = lexical_join(
        scenario_parent,
        Path::new(scenario.profile_metadata_path.value()),
    )?;
    let excluded_regions_path = scenario
        .excluded_regions_path
        .as_ref()
        .map(|path| lexical_join(scenario_parent, Path::new(path.value())))
        .transpose()?;
    Ok(F5BootstrapConfig {
        general_path,
        scenario_path,
        profile_metadata_path,
        excluded_regions_path,
        reference_path,
        documents: F5BootstrapDocuments {
            general: general_document,
            scenario: scenario_document,
        },
        general,
        scenario,
    })
}

/// Materializa el perfil inmediatamente después del bootstrap usando la misma
/// contabilidad v9. La semilla se resuelve únicamente si las seis entradas
/// físicas han sido admitidas.
fn materialize_f5_config_after_bootstrap(
    bootstrap: F5BootstrapConfig,
    entropy: &mut impl EntropySource,
    other: crate::f5_preflight::F5OtherResources,
) -> Result<EffectiveF5Config, Diagnostic> {
    let mut effective = prepare_f5_config_after_bootstrap(bootstrap)?;
    let mut ledger = effective.admission_ledger();
    let mut preflight =
        crate::f5_preflight::admit_f5_resources_under_ledger(&effective, other, &mut ledger)?;
    // Solo un resultado completamente admitido puede pedir entropía. Las
    // semillas declaradas conservan exactamente la misma resolución.
    let (root_seed, origin) = resolve_root_seed(
        effective
            .general
            .seed
            .as_ref()
            .map(|seed| seed.value().as_str()),
        entropy,
    )
    .map_err(|error| error.with("path", effective.general_path.display().to_string()))?;
    effective.seed = root_seed.to_hex();
    effective.seed_origin = origin.as_str().to_owned();
    preflight.bind_configuration(&effective)?;
    // El sellado cierra la contabilidad y une el receipt: identidades de las
    // diez entradas posibles en orden de rol y fases en su orden cerrado.
    let phases = ledger.seal_phases()?;
    let identities = crate::f5_preflight::sealed_input_identities(
        &effective.admitted_document_identities(),
        &preflight,
    );
    let receipt = F5AdmissionReceipt::seal(identities, phases)?;
    effective.admitted = Some(preflight.seal(receipt));
    effective.admitted_ledger = Some(ledger);
    Ok(effective)
}

/// Construye el estado tipado previo a referencia sin resolver entropía.
fn prepare_f5_config_after_bootstrap(
    bootstrap: F5BootstrapConfig,
) -> Result<EffectiveF5Config, Diagnostic> {
    let F5BootstrapConfig {
        general_path,
        scenario_path,
        profile_metadata_path,
        excluded_regions_path,
        reference_path,
        documents,
        general,
        scenario,
    } = bootstrap;
    let population_ids = scenario
        .subpopulations
        .iter()
        .map(|population| population.population_id.value().clone())
        .collect::<BTreeSet<_>>();
    let mut ledger = profile_ledger_after_bootstrap(&scenario, &documents)?;
    let loaded_profile = load_f5_profile_under_budget(
        &profile_metadata_path,
        scenario.profile_metadata_sha256.value(),
        &population_ids,
        &mut ledger,
    )?;
    ledger.transition_to(AdmissionPhase::ReferenceAndMaskScan)?;
    let (profile, profile_input_identities) = loaded_profile.into_parts();
    let input_documents = F5InputDocuments::from_observed(
        documents.general().size_bytes(),
        documents.scenario().size_bytes(),
        &profile_input_identities,
    )?;
    Ok(EffectiveF5Config {
        general_path,
        scenario_path,
        profile_metadata_path,
        excluded_regions_path,
        reference_path,
        seed: String::new(),
        seed_origin: String::new(),
        ploidy: u8::try_from(*general.ploidy.value())
            .map_err(|_| config_error_f5(&general.source, "ploidy", "ploidía fuera de rango"))?,
        reference_id: general.reference_id.value().clone(),
        reference_sha256: general.reference_sha256.value().clone(),
        scenario_id: general.scenario_id.value().clone(),
        profile,
        bootstrap_documents: documents,
        input_documents,
        profile_input_identities,
        admission_ledger: ledger,
        admitted: None,
        admitted_ledger: None,
        general,
        scenario,
    })
}

/// El nombre del general tiene que ser publicable.
///
/// Los cinco enlaces que los documentos declaran —escenario, referencia,
/// metadata, máscara y los cinco TSV de perfil— ya pasan por
/// `require_relative_path`, pero el nombre del **propio** general lo elige el
/// usuario en la línea de órdenes y nadie lo validaba. `effective-config.json`
/// lo publica como `source_path`, cuyo `$defs/path` admite exactamente
/// `[A-Za-z0-9._-]`, de modo que un `mi general.toml` o un
/// `general (copia).toml` ejecutaba admisión, layout, búsqueda, fundadores,
/// identidades, los ocho artefactos y el recuento MAN-70 **completos** antes
/// de morir al serializar el efectivo. La frontera correcta es la admisión: se
/// rechaza aquí, antes de abrir el fichero, y el diagnóstico dice qué
/// renombrar en lugar de culpar al esquema del documento producido.
fn require_publishable_general_name(general_path: &Path) -> Result<(), Diagnostic> {
    let publishable = general_path
        .file_name()
        .and_then(|name| name.to_str())
        .is_some_and(|name| {
            !name.is_empty()
                && !matches!(name, "." | "..")
                && name
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
        });
    if publishable {
        return Ok(());
    }
    Err(Diagnostic::f5_configuration(
        general_path.display().to_string(),
        "general_path",
        "el nombre del documento general se publica en effective-config.json y solo admite [A-Za-z0-9._-]; renómbralo antes de ejecutar",
    ))
}

fn absolute_normalized_path(path: &Path) -> Result<PathBuf, Diagnostic> {
    use std::path::Component;

    let absolute = if path.is_absolute() {
        path.to_path_buf()
    } else {
        std::env::current_dir()
            .map_err(|error| {
                Diagnostic::io(
                    "GPS-IO-001",
                    "No se pudo resolver el directorio de trabajo.",
                )
                .with("detail", error.to_string())
            })?
            .join(path)
    };
    let mut normalized = PathBuf::new();
    for component in absolute.components() {
        match component {
            Component::Prefix(_) | Component::RootDir => normalized.push(component.as_os_str()),
            Component::CurDir => {}
            Component::Normal(value) => normalized.push(value),
            // La ruta inicial viene de la CLI y se ancla al cwd; `..` solo se
            // prohíbe en los enlaces declarados, donde impediría la auditoría
            // de la procedencia general→escenario→recurso.
            Component::ParentDir => {
                normalized.pop();
            }
        }
    }
    Ok(normalized)
}

fn lexical_join(base: &Path, child: &Path) -> Result<PathBuf, Diagnostic> {
    use std::path::Component;
    let mut result = PathBuf::new();
    for component in base.components().chain(child.components()) {
        match component {
            Component::Prefix(_) | Component::RootDir => result.push(component.as_os_str()),
            Component::CurDir => {}
            Component::Normal(value) => result.push(value),
            Component::ParentDir => {
                return Err(Diagnostic::f5_configuration(
                    child.display().to_string(),
                    "path",
                    "las rutas F5 no admiten ..",
                ));
            }
        }
    }
    Ok(result)
}

fn config_error_f5(
    source: &crate::configuration::AdmittedToml,
    field: &str,
    detail: &str,
) -> Diagnostic {
    Diagnostic::f5_configuration(source.path().display().to_string(), field, detail)
}

pub fn validate_general_config(
    general: &DeclaredGeneralConfig,
) -> Result<ValidatedGeneralSummary, Diagnostic> {
    effective_ploidy(general)?;
    ReferenceId::parse(general.reference().id.value(), "configuration")?;
    if general.reference().sequences.is_empty() {
        return Err(config_error(
            general,
            "reference.sequences",
            "se exige al menos una secuencia",
        ));
    }
    let mut ids = BTreeSet::new();
    let mut total = 0_u64;
    for sequence in &general.reference().sequences {
        SequenceId::parse(sequence.id.value(), general.reference().id.value())?;
        if !ids.insert(sequence.id.value()) {
            return Err(config_error(
                general,
                "reference.sequences.id",
                format!("identificador repetido: {}", sequence.id.value()),
            ));
        }
        if *sequence.length.value() == 0 {
            return Err(config_error(
                general,
                "reference.sequences.length",
                format!("{} tiene longitud cero", sequence.id.value()),
            ));
        }
        total = total.checked_add(*sequence.length.value()).ok_or_else(|| {
            config_error(
                general,
                "reference.sequences.length",
                "el tamaño total del genoma desborda u64",
            )
        })?;
    }
    if let DeclaredReferenceSource::Simulated {
        composition,
        unknown_base_fraction,
    } = &general.reference().source
    {
        normalize_composition(general, composition)?;
        parse_unknown_base_fraction(general, unknown_base_fraction)?;
    }
    if !general.reduced_references().is_empty()
        && !matches!(
            &general.reference().source,
            DeclaredReferenceSource::Supplied { .. }
        )
    {
        return Err(config_error(
            general,
            "reduced_references",
            "F1 solo permite crear referencias reducidas desde una referencia suministrada",
        ));
    }
    let mut reference_ids = BTreeSet::from([general.reference().id.value().as_str()]);
    for reduced in general.reduced_references() {
        ReferenceId::parse(reduced.id.value(), "configuration")?;
        if !reference_ids.insert(reduced.id.value()) {
            return Err(config_error(
                general,
                "reduced_references.id",
                format!("identificador repetido: {}", reduced.id.value()),
            ));
        }
        if reduced.source_reference_id.value() != general.reference().id.value() {
            return Err(config_error(
                general,
                "reduced_references.source_reference_id",
                "F1 solo permite reducir la referencia inicial declarada",
            ));
        }
        if *reduced.region_size_bases.value() == 0 {
            return Err(config_error(
                general,
                "reduced_references.region_size_bases",
                "el tamaño de región debe ser positivo",
            ));
        }
    }
    Ok(ValidatedGeneralSummary {
        total_genome_size: total,
    })
}

pub fn build_effective_config(
    general: &DeclaredGeneralConfig,
) -> Result<EffectiveConfig, Diagnostic> {
    build_effective_config_with_entropy(general, &mut SystemEntropy)
}

pub fn build_effective_config_with_entropy(
    general: &DeclaredGeneralConfig,
    entropy: &mut impl EntropySource,
) -> Result<EffectiveConfig, Diagnostic> {
    let summary = validate_general_config(general)?;
    let ploidy = effective_ploidy(general)?;
    let normalization = EffectiveNormalization {
        tolerance: effective_normalization_tolerance(general)?,
    };
    let (root_seed, seed_origin) = resolve_root_seed(
        general.seed().map(|declared| declared.value().as_str()),
        entropy,
    )
    .map_err(|diagnostic| {
        let mut diagnostic = diagnostic.with("path", general.source().path().display().to_string());
        if let Some(seed) = general.seed() {
            diagnostic = diagnostic.with("location", seed.location().to_string());
        }
        diagnostic
    })?;
    let composition = match &general.reference().source {
        DeclaredReferenceSource::Simulated {
            composition,
            unknown_base_fraction,
        } => Some(EffectiveComposition {
            values: normalize_composition(general, composition)?,
            unknown_base_fraction: parse_unknown_base_fraction(general, unknown_base_fraction)?,
        }),
        DeclaredReferenceSource::Supplied { .. } => None,
    };
    let sequences: Vec<EffectiveSequence> = general
        .reference()
        .sequences
        .iter()
        .map(|sequence| EffectiveSequence {
            id: sequence.id.value().clone(),
            class: *sequence.class.value(),
            length: *sequence.length.value(),
        })
        .collect();
    let mut reduced_references = Vec::with_capacity(general.reduced_references().len());
    for reduction in general.reduced_references() {
        let (target, target_bases) = match &reduction.target {
            DeclaredReductionTarget::Bases(value) => (
                EffectiveReductionTarget::Bases {
                    declared: *value.value(),
                },
                *value.value(),
            ),
            DeclaredReductionTarget::ReductionPercent(value) => {
                let exact = ExactDecimal::parse("reduction_percent", value.literal()).map_err(
                    |diagnostic| {
                        diagnostic
                            .with("path", general.source().path().display().to_string())
                            .with("location", value.location().to_string())
                    },
                )?;
                (
                    EffectiveReductionTarget::ReductionPercent {
                        literal: value.literal().to_owned(),
                    },
                    exact
                        .retained_after_percent(summary.total_genome_size)
                        .map_err(|diagnostic| {
                            diagnostic
                                .with("path", general.source().path().display().to_string())
                                .with("location", value.location().to_string())
                        })?,
                )
            }
        };
        if target_bases == 0 || target_bases >= summary.total_genome_size {
            return Err(Diagnostic::configuration(
                "GPS-CONF-021",
                "El objetivo de la referencia reducida no es menor que el espacio elegible.",
            )
            .with("path", general.source().path().display().to_string())
            .with("reference", reduction.id.value())
            .with("target_bases", target_bases.to_string())
            .with("eligible_bases", summary.total_genome_size.to_string()));
        }
        reduced_references.push(EffectiveReducedReference {
            id: reduction.id.value().clone(),
            source_reference_id: reduction.source_reference_id.value().clone(),
            region_size_bases: *reduction.region_size_bases.value(),
            target,
            target_bases,
            selection_subseed: derive_subseed(
                root_seed,
                RuntimeSubseedScope::ReferenceReduction {
                    source_reference_id: reduction.source_reference_id.value(),
                    derived_reference_id: reduction.id.value(),
                },
            )?,
        });
    }
    let reference_source = match &general.reference().source {
        DeclaredReferenceSource::Simulated { .. } => EffectiveReferenceSource::Simulated,
        DeclaredReferenceSource::Supplied { fasta } => EffectiveReferenceSource::Supplied {
            original_path: resolve_link(general.source().path(), fasta.value()),
        },
    };
    let repeat_admission = effective_repeat_admission(general)?;
    let repeat_loci = effective_repeat_loci(
        general,
        &sequences,
        normalization.tolerance,
        repeat_admission,
    )?;
    let reproducibility = ReproducibilityMetadata::new(root_seed, seed_origin);
    let initial_populations = effective_initial_populations(general, root_seed)?;
    let population_snvs = effective_population_snvs(general)?;
    let direct_events = effective_direct_events(general, &initial_populations, &sequences)?;
    let technical_event = effective_technical_event(general, &initial_populations, &direct_events)?;
    let vcf_output = general
        .vcf_output()
        .map(|output| match output.mode.value().as_str() {
            "text" => VcfOutputMode::Text,
            "bgzf" => VcfOutputMode::Bgzf,
            _ => unreachable!("the TOML adapter validates vcf_output.mode"),
        });
    if !direct_events.is_empty() && population_snvs.is_none() {
        return Err(config_error(
            general,
            "direct_events,population_snvs",
            "direct_events exige una matriz SNV poblacional para ligar genotipo y fase",
        ));
    }
    if let Some(EffectivePopulationSnvs::Generate {
        parameters: generation,
    }) = &population_snvs
    {
        validate_generation_config(generation).map_err(|error| {
            config_error(general, "population_snvs.generate", error.to_string())
        })?;
    }
    if let Some(PopulationSnvsGenerationConfig {
        spatial: GeneratedSpatialMode::WeightedRegions { regions, .. },
        ..
    }) = population_snvs_generation_from_effective(population_snvs.as_ref())
    {
        for region in regions {
            let valid = sequences
                .iter()
                .find(|sequence| sequence.id() == region.contig)
                .is_some_and(|sequence| {
                    region.start_1based > 0 && region.end_1based <= sequence.length()
                });
            if !valid {
                return Err(config_error(
                    general,
                    "population_snvs.generate.regions",
                    format!(
                        "región fuera de la referencia efectiva: {}:{}..{}",
                        region.contig, region.start_1based, region.end_1based
                    ),
                ));
            }
        }
    }
    Ok(EffectiveConfig {
        source_path: general.source().path().display().to_string(),
        reference_id: general.reference().id.value().clone(),
        ploidy,
        normalization,
        reference_source,
        composition,
        reduced_references,
        total_genome_size: summary.total_genome_size,
        sequences,
        configuration_input: EffectiveConfigurationInput {
            original_path: general.source().path().display().to_string(),
            bytes: Arc::from(general.source().bytes()),
        },
        reproducibility,
        population_snvs,
        initial_populations,
        direct_events,
        technical_event,
        vcf_output,
        document_version: match general.document_version() {
            "4" => DeclaredDocumentVersion::V4,
            _ => DeclaredDocumentVersion::V2,
        },
        repeat_loci,
        repeat_admission,
    })
}

fn effective_technical_event(
    general: &DeclaredGeneralConfig,
    populations: &[EffectiveInitialPopulation],
    direct_events: &[EffectiveDirectEvent],
) -> Result<Option<EffectiveTechnicalEvent>, Diagnostic> {
    let Some(declared) = general.technical_event() else {
        return Ok(None);
    };
    let event_id = EventId::parse(declared.event_id.value().clone(), "technical_event")
        .map_err(|error| with_initial_population_path(general, error))?;
    if direct_events
        .iter()
        .any(|event| event.event_id() == &event_id)
    {
        return Err(config_error(
            general,
            "technical_event.event_id",
            format!("identificador de evento repetido en el ámbito de ejecución: {event_id}"),
        ));
    }
    let declared_individuals = declared.individual_ids.as_slice();
    let mut known_individuals = BTreeSet::new();
    for population in populations {
        known_individuals.extend(population.individual_ids().iter().cloned());
    }
    let mut individual_ids = Vec::with_capacity(declared_individuals.len());
    let mut selected = BTreeSet::new();
    for declared_id in declared_individuals {
        let individual_id = IndividualId::parse(declared_id.value().clone(), "technical_event")
            .map_err(|error| with_initial_population_path(general, error))?;
        if !known_individuals.contains(&individual_id) {
            return Err(config_error(
                general,
                "technical_event.individual_ids",
                format!("individuo no declarado en initial_populations: {individual_id}"),
            ));
        }
        if !selected.insert(individual_id.clone()) {
            return Err(config_error(
                general,
                "technical_event.individual_ids",
                format!("individuo repetido: {individual_id}"),
            ));
        }
        individual_ids.push(individual_id);
    }
    let target = match &declared.target {
        DeclaredTechnicalEventTarget::CoverageWgs(value) => {
            let coverage = value.literal().parse::<f64>().map_err(|_| {
                config_error(
                    general,
                    "technical_event.target_coverage_wgs",
                    "debe ser un f64 válido",
                )
            })?;
            if !coverage.is_finite() || coverage < 0.0 {
                return Err(config_error(
                    general,
                    "technical_event.target_coverage_wgs",
                    "debe ser finita y mayor o igual que cero",
                ));
            }
            EffectiveTechnicalEventTarget::CoverageWgs(coverage)
        }
        DeclaredTechnicalEventTarget::ReadCount(value) => {
            EffectiveTechnicalEventTarget::ReadCount(*value.value())
        }
    };
    let sequencing_errors = match declared
        .sequencing_errors
        .as_ref()
        .map(|value| value.value().as_str())
    {
        None | Some("enabled") => SequencingErrorsMode::Enabled,
        Some("disabled") => SequencingErrorsMode::Disabled,
        Some(_) => {
            return Err(config_error(
                general,
                "technical_event.sequencing_errors",
                "debe ser enabled o disabled",
            ));
        }
    };
    let read_length = declared
        .read_length_bases
        .as_ref()
        .map_or(u64::from(ILLUMINA_SE_MINIMAL_READ_LENGTH), |value| {
            *value.value()
        });
    if !(1..=150).contains(&read_length) {
        return Err(config_error(
            general,
            "technical_event.read_length_bases",
            "debe estar entre 1 y 150 y no superar el inserto mínimo del perfil",
        ));
    }
    let read_length = u16::try_from(read_length).map_err(|_| {
        config_error(
            general,
            "technical_event.read_length_bases",
            "no cabe en la longitud de lectura del perfil",
        )
    })?;
    let chimera_rate = declared.chimera_rate.as_ref().map_or(Ok(0.0), |value| {
        value.literal().parse::<f64>().map_err(|_| {
            config_error(
                general,
                "technical_event.chimera_rate",
                "debe ser un f64 válido",
            )
        })
    })?;
    if !chimera_rate.is_finite() || !(0.0..=1.0).contains(&chimera_rate) {
        return Err(config_error(
            general,
            "technical_event.chimera_rate",
            "debe ser una fracción finita entre 0 y 1",
        ));
    }
    let index_hopping_rate = declared
        .index_hopping_rate
        .as_ref()
        .map_or(Ok(0.0), |value| {
            value.literal().parse::<f64>().map_err(|_| {
                config_error(
                    general,
                    "technical_event.index_hopping_rate",
                    "debe ser un f64 válido",
                )
            })
        })?;
    if !index_hopping_rate.is_finite() || !(0.0..=1.0).contains(&index_hopping_rate) {
        return Err(config_error(
            general,
            "technical_event.index_hopping_rate",
            "debe ser una fracción finita entre 0 y 1",
        ));
    }
    let biological_contamination = {
        let minimum_percent =
            declared
                .contamination_min_percent
                .as_ref()
                .map_or(Ok(0.0), |value| {
                    value.literal().parse::<f64>().map_err(|_| {
                        config_error(
                            general,
                            "technical_event.contamination_min_percent",
                            "debe ser un f64 válido",
                        )
                    })
                })?;
        let maximum_percent =
            declared
                .contamination_max_percent
                .as_ref()
                .map_or(Ok(0.0), |value| {
                    value.literal().parse::<f64>().map_err(|_| {
                        config_error(
                            general,
                            "technical_event.contamination_max_percent",
                            "debe ser un f64 válido",
                        )
                    })
                })?;
        let recipient_ids = declared
            .contamination_recipient_ids
            .as_ref()
            .map(|ids| {
                ids.iter()
                    .map(|id| {
                        IndividualId::parse(
                            id.value().clone(),
                            "technical_event.contamination_recipient_ids",
                        )
                        .map_err(|error| with_initial_population_path(general, error))
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()
            })
            .transpose()?
            .map(|mut ids| {
                ids.sort();
                ids
            });
        BiologicalContaminationSettings::new(minimum_percent, maximum_percent, recipient_ids)
            .map_err(|error| {
                config_error(
                    general,
                    "technical_event.contamination_min_percent/contamination_max_percent",
                    error.to_string(),
                )
            })?
    };
    for recipient in biological_contamination
        .configured_recipients()
        .unwrap_or(&[])
    {
        if !selected.contains(recipient) {
            return Err(config_error(
                general,
                "technical_event.contamination_recipient_ids",
                format!("receptora fuera de individual_ids: {recipient}"),
            ));
        }
    }
    let mut libraries = Vec::with_capacity(declared.libraries.len());
    let mut library_individuals = BTreeSet::new();
    let mut i7_values = BTreeSet::new();
    let mut i5_values = BTreeSet::new();
    for library in &declared.libraries {
        let individual_id = IndividualId::parse(
            library.individual_id.value().clone(),
            "technical_event.libraries",
        )
        .map_err(|error| with_initial_population_path(general, error))?;
        if !selected.contains(&individual_id) || !library_individuals.insert(individual_id.clone())
        {
            return Err(config_error(
                general,
                "technical_event.libraries.individual_id",
                "debe coincidir exactamente y sin repetidos con individual_ids",
            ));
        }
        for (field, value) in [("i7", library.i7.value()), ("i5", library.i5.value())] {
            if value.len() != 8
                || !value
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
            {
                return Err(config_error(
                    general,
                    "technical_event.libraries",
                    format!("{field} debe contener exactamente ocho bases A/C/G/T"),
                ));
            }
        }
        if !i7_values.insert(library.i7.value().clone())
            || !i5_values.insert(library.i5.value().clone())
        {
            return Err(config_error(
                general,
                "technical_event.libraries",
                "cada componente UDI i7 e i5 debe ser único en la tanda",
            ));
        }
        libraries.push(EffectiveTechnicalLibrary {
            individual_id,
            i7: library.i7.value().clone(),
            i5: library.i5.value().clone(),
        });
    }
    if library_individuals != selected {
        return Err(config_error(
            general,
            "technical_event.libraries",
            "debe existir exactamente una biblioteca por individuo seleccionado",
        ));
    }
    if index_hopping_rate > 0.0 && libraries.len() < 2 {
        return Err(config_error(
            general,
            "technical_event.index_hopping_rate",
            "una tasa positiva exige al menos dos bibliotecas declaradas",
        ));
    }
    libraries.sort_by(|left, right| left.individual_id.cmp(&right.individual_id));
    Ok(Some(EffectiveTechnicalEvent {
        event_id,
        individual_ids,
        target,
        sequencing_errors,
        read_length,
        chimera_rate,
        index_hopping_rate,
        biological_contamination,
        libraries,
    }))
}

fn effective_direct_events(
    general: &DeclaredGeneralConfig,
    populations: &[EffectiveInitialPopulation],
    sequences: &[EffectiveSequence],
) -> Result<Vec<EffectiveDirectEvent>, Diagnostic> {
    let mut individuals = std::collections::BTreeMap::new();
    for population in populations {
        for individual in population.individual_ids() {
            individuals.insert(individual.clone(), population.ploidy());
        }
    }
    let mut registry = IdentityRegistry::default();
    general
        .direct_events()
        .iter()
        .map(|declared| {
            let event_id = EventId::parse(declared.event_id.value().clone(), "direct_events")
                .map_err(|error| with_initial_population_path(general, error))?;
            let event_id = registry
                .register_event("execution", event_id)
                .map_err(|error| with_initial_population_path(general, error))?;
            let individual_id =
                IndividualId::parse(declared.individual_id.value().clone(), "direct_events")
                    .map_err(|error| with_initial_population_path(general, error))?;
            let ploidy = individuals.get(&individual_id).ok_or_else(|| {
                config_error(
                    general,
                    "direct_events.individual_id",
                    format!("individuo no declarado en initial_populations: {individual_id}"),
                )
            })?;
            let phase_index = u8::try_from(*declared.phase.value()).map_err(|_| {
                config_error(general, "direct_events.phase", "no se puede representar")
            })?;
            if phase_index >= ploidy.copy_count() {
                return Err(config_error(
                    general,
                    "direct_events.phase",
                    format!(
                        "la fase {phase_index} no existe para ploidía {}",
                        ploidy.copy_count()
                    ),
                ));
            }
            let phase = HaplotypePhase::from_copy_index(phase_index)
                .map_err(|_| config_error(general, "direct_events.phase", "fase inválida"))?;
            let sequence_id = SequenceId::parse(declared.contig.value().clone(), "direct_events")
                .map_err(|error| with_initial_population_path(general, error))?;
            if !sequences
                .iter()
                .any(|sequence| sequence.id() == sequence_id.as_str())
            {
                return Err(config_error(
                    general,
                    "direct_events.contig",
                    format!("contig ausente de la referencia efectiva: {sequence_id}"),
                ));
            }
            match declared.kind.value().as_str() {
                "snp" => {
                    let position_1based = *declared
                        .position_1based
                        .as_ref()
                        .expect("parser checked")
                        .value();
                    let reference = single_canonical_base(
                        general,
                        "direct_events.ref",
                        declared.reference.as_ref().expect("parser checked").value(),
                    )?;
                    let alternate = single_canonical_base(
                        general,
                        "direct_events.alt",
                        declared.alternate.as_ref().expect("parser checked").value(),
                    )?;
                    if reference == alternate {
                        return Err(config_error(
                            general,
                            "direct_events.ref,direct_events.alt",
                            "deben ser bases canónicas distintas",
                        ));
                    }
                    Ok(EffectiveDirectEvent::Snp {
                        event_id,
                        individual_id,
                        phase,
                        sequence_id,
                        position_1based,
                        reference,
                        alternate,
                    })
                }
                "insertion" => {
                    let bases = declared
                        .bases
                        .as_ref()
                        .expect("parser checked")
                        .value()
                        .as_bytes();
                    if bases.is_empty()
                        || bases.iter().any(|base| {
                            !matches!(base.to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T')
                        })
                    {
                        return Err(config_error(
                            general,
                            "direct_events.bases",
                            "debe contener una o más bases A/C/G/T",
                        ));
                    }
                    Ok(EffectiveDirectEvent::Insertion {
                        event_id,
                        individual_id,
                        phase,
                        sequence_id,
                        point: *declared.point.as_ref().expect("parser checked").value(),
                        bases: Arc::from(
                            bases.iter().map(u8::to_ascii_uppercase).collect::<Vec<_>>(),
                        ),
                    })
                }
                "deletion" => Ok(EffectiveDirectEvent::Deletion {
                    event_id,
                    individual_id,
                    phase,
                    sequence_id,
                    start: *declared.start.as_ref().expect("parser checked").value(),
                    end: *declared.end.as_ref().expect("parser checked").value(),
                }),
                _ => unreachable!("parser validates direct event kind"),
            }
        })
        .collect()
}

fn single_canonical_base(
    general: &DeclaredGeneralConfig,
    field: &str,
    value: &str,
) -> Result<u8, Diagnostic> {
    let bytes = value.as_bytes();
    if bytes.len() != 1 || !matches!(bytes[0].to_ascii_uppercase(), b'A' | b'C' | b'G' | b'T') {
        return Err(config_error(
            general,
            field,
            "debe ser exactamente una base A/C/G/T",
        ));
    }
    Ok(bytes[0].to_ascii_uppercase())
}

fn effective_initial_populations(
    general: &DeclaredGeneralConfig,
    root_seed: crate::randomness::RootSeed,
) -> Result<Vec<EffectiveInitialPopulation>, Diagnostic> {
    let mut population_ids = BTreeSet::new();
    let mut derivation_registry = IdentityRegistry::default();
    let mut execution_registry = IdentityRegistry::default();
    let mut effective = Vec::with_capacity(general.initial_populations().len());

    for declared in general.initial_populations() {
        let population_id = PopulationId::parse(
            declared.population_id.value().clone(),
            "initial_populations",
        )
        .map_err(|error| with_initial_population_path(general, error))?;
        if !population_ids.insert(population_id.clone()) {
            return Err(config_error(
                general,
                "initial_populations.population_id",
                format!("identificador repetido: {population_id}"),
            ));
        }

        let individual_count = *declared.individual_count.value();
        if individual_count == 0 {
            return Err(config_error(
                general,
                "initial_populations.individual_count",
                "debe ser mayor que cero",
            ));
        }
        let individual_count = usize::try_from(individual_count).map_err(|_| {
            config_error(
                general,
                "initial_populations.individual_count",
                "no se puede representar en esta plataforma",
            )
        })?;
        let ploidy = u8::try_from(*declared.ploidy.value())
            .ok()
            .and_then(|value| Ploidy::new(value).ok())
            .ok_or_else(|| {
                config_error(
                    general,
                    "initial_populations.ploidy",
                    format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
                )
            })?;
        let genotype_model = match declared.genotype_model.value().as_str() {
            "hwe" => EffectiveInitialGenotypeModel::HardyWeinberg,
            "f_diploid" => {
                if ploidy != Ploidy::DIPLOID {
                    return Err(config_error(
                        general,
                        "initial_populations.genotype_model,initial_populations.ploidy",
                        "f_diploid solo se admite con ploidía 2",
                    ));
                }
                let coefficient = declared
                    .inbreeding_coefficient
                    .as_ref()
                    .ok_or_else(|| {
                        config_error(
                            general,
                            "initial_populations.inbreeding_coefficient",
                            "f_diploid exige inbreeding_coefficient",
                        )
                    })?
                    .literal()
                    .parse::<f64>()
                    .map_err(|_| {
                        config_error(
                            general,
                            "initial_populations.inbreeding_coefficient",
                            "debe ser un f64 válido",
                        )
                    })?;
                if !coefficient.is_finite() || !(0.0..=1.0).contains(&coefficient) {
                    return Err(config_error(
                        general,
                        "initial_populations.inbreeding_coefficient",
                        "debe ser finito y pertenecer a [0,1]",
                    ));
                }
                EffectiveInitialGenotypeModel::FDiploid {
                    inbreeding_coefficient: coefficient,
                }
            }
            _ => unreachable!("parser validates initial genotype model"),
        };

        let (individual_ids, founder_keys) = match (
            &declared.individual_ids,
            &declared.founder_keys,
        ) {
            (Some(provided), None) => {
                if provided.len() != individual_count {
                    return Err(config_error(
                        general,
                        "initial_populations.individual_ids,initial_populations.individual_count",
                        "individual_ids debe contener exactamente individual_count identificadores",
                    ));
                }
                let individual_ids = provided
                    .iter()
                    .map(|id| {
                        let id = IndividualId::parse(id.value().clone(), "initial_populations")
                            .map_err(|error| with_initial_population_path(general, error))?;
                        execution_registry
                            .register_individual("execution", id)
                            .map_err(|error| with_initial_population_path(general, error))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                (individual_ids, None)
            }
            (None, Some(provided)) => {
                if provided.len() != individual_count {
                    return Err(config_error(
                        general,
                        "initial_populations.founder_keys,initial_populations.individual_count",
                        "founder_keys debe contener exactamente individual_count claves",
                    ));
                }
                let mut founder_key_registry = BTreeSet::new();
                let founder_keys = provided
                    .iter()
                    .map(|key| {
                        let key = FounderKey::parse(key.value().clone(), "initial_populations")
                            .map_err(|error| with_initial_population_path(general, error))?;
                        if !founder_key_registry.insert(key.clone()) {
                            return Err(config_error(
                                general,
                                "initial_populations.founder_keys",
                                format!("founder_key repetida: {key}"),
                            ));
                        }
                        Ok(key)
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()?;
                let individual_ids = founder_keys
                    .iter()
                    .map(|founder_key| {
                        let request = AutomaticIdRequest::new(
                            IdClass::Individual,
                            population_id.as_str(),
                            &[("founder_key", founder_key.as_str())],
                        )
                        .map_err(|error| with_initial_population_path(general, error))?;
                        let id = derivation_registry
                            .derive_individual(root_seed, &request)
                            .map_err(|error| with_initial_population_path(general, error))?;
                        execution_registry
                            .register_individual("execution", id)
                            .map_err(|error| with_initial_population_path(general, error))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                (individual_ids, Some(founder_keys))
            }
            (Some(_), Some(_)) => {
                return Err(config_error(
                    general,
                    "initial_populations.individual_ids,initial_populations.founder_keys",
                    "individual_ids y founder_keys son mutuamente excluyentes",
                ));
            }
            (None, None) => {
                return Err(config_error(
                    general,
                    "initial_populations.founder_keys",
                    "founder_keys es obligatorio cuando se omite individual_ids",
                ));
            }
        };
        effective.push(EffectiveInitialPopulation {
            population_id,
            individual_ids,
            founder_keys,
            ploidy,
            genotype_model,
        });
    }
    Ok(effective)
}

fn with_initial_population_path(
    general: &DeclaredGeneralConfig,
    diagnostic: Diagnostic,
) -> Diagnostic {
    diagnostic.with("path", general.source().path().display().to_string())
}

fn population_snvs_generation_from_effective(
    value: Option<&EffectivePopulationSnvs>,
) -> Option<&PopulationSnvsGenerationConfig> {
    match value {
        Some(EffectivePopulationSnvs::Generate { parameters }) => Some(parameters),
        _ => None,
    }
}

fn effective_population_snvs(
    general: &DeclaredGeneralConfig,
) -> Result<Option<EffectivePopulationSnvs>, Diagnostic> {
    let Some(declared) = general.population_snvs() else {
        return Ok(None);
    };
    match declared {
        DeclaredPopulationSnvs::Supplied { input_path } => {
            return Ok(Some(EffectivePopulationSnvs::Supplied {
                input_path: resolve_link(general.source().path(), input_path.value()),
            }));
        }
        DeclaredPopulationSnvs::Reuse { metadata_path } => {
            return Ok(Some(EffectivePopulationSnvs::Reuse {
                metadata_path: resolve_link(general.source().path(), metadata_path.value()),
            }));
        }
        DeclaredPopulationSnvs::Generate { parameters: _ } => {}
    }
    let DeclaredPopulationSnvs::Generate { parameters } = declared else {
        unreachable!("handled non-generation population_snvs modes")
    };
    let decimal = |field: &str, value: &crate::configuration::DeclaredNumber| {
        value
            .literal()
            .parse::<f64>()
            .map_err(|_| config_error(general, field, "debe ser un f64 válido"))
    };
    let count = if let Some(value) = &parameters.locus_count {
        GeneratedLocusCount::Exact(*value.value())
    } else if let Some(value) = &parameters.mean_loci {
        GeneratedLocusCount::PoissonMean(decimal("population_snvs.generate.mean_loci", value)?)
    } else {
        GeneratedLocusCount::PoissonDensityPerBase(decimal(
            "population_snvs.generate.mean_density_per_base",
            parameters.mean_density_per_base.as_ref().ok_or_else(|| {
                config_error(general, "population_snvs.generate", "recuento ausente")
            })?,
        )?)
    };
    let spatial = if parameters.spatial_mode.value() == "uniform" {
        GeneratedSpatialMode::Uniform
    } else {
        GeneratedSpatialMode::WeightedRegions {
            outside_weight: decimal(
                "population_snvs.generate.outside_weight",
                parameters.outside_weight.as_ref().ok_or_else(|| {
                    config_error(
                        general,
                        "population_snvs.generate.outside_weight",
                        "falta outside_weight",
                    )
                })?,
            )?,
            regions: parameters
                .regions
                .iter()
                .map(|region| {
                    Ok(GeneratedWeightedRegion {
                        contig: region.contig.value().clone(),
                        start_1based: *region.start.value(),
                        end_1based: *region.end.value(),
                        weight: decimal("population_snvs.generate.regions.weight", &region.weight)?,
                    })
                })
                .collect::<Result<Vec<_>, Diagnostic>>()?,
        }
    };
    Ok(Some(EffectivePopulationSnvs::Generate { parameters: PopulationSnvsGenerationConfig {
        count,
        spatial,
        substitution_profile: match parameters.substitution_profile.value().as_str() {
            "uniform_snv12_v1" => GeneratedSubstitutionProfile::UniformSnv12V1,
            "human_harris2017_snv12_approx_v1" => {
                GeneratedSubstitutionProfile::HumanHarris2017Snv12ApproxV1
            }
            _ => return Err(config_error(general, "population_snvs.generate.substitution_profile", "perfil de sustitución desconocido")),
        },
        populations: parameters
            .populations
            .iter()
            .map(|population| {
                let spectrum = match population.frequency_model.value().as_str() {
                    "uniform" => None,
                    "human_1000g_global_sfs_approx_v1" => Some(human_1000g_spectrum(
                        population.frequency_mode.as_ref().map_or("stochastic", |mode| mode.value()),
                    )?),
                    "spectrum" => Some(GeneratedFrequencySpectrum {
                        mode: frequency_assignment_mode(
                            population.frequency_mode.as_ref().map_or("stochastic", |mode| mode.value()),
                        )?,
                        components: population.spectrum.iter().map(|component| {
                            let weight = decimal("population_snvs.generate.populations.spectrum.weight", &component.weight)?;
                            match component.kind.value().as_str() {
                                "fixed" => Ok(FrequencySpectrumComponent::Fixed {
                                    value: decimal("population_snvs.generate.populations.spectrum.fixed_alt_af", component.fixed_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.fixed_alt_af", "fixed exige fixed_alt_af"))?)?,
                                    weight,
                                }),
                                "interval" => Ok(FrequencySpectrumComponent::Interval {
                                    min: decimal("population_snvs.generate.populations.spectrum.min_alt_af", component.min_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.min_alt_af", "interval exige min_alt_af"))?)?,
                                    max: decimal("population_snvs.generate.populations.spectrum.max_alt_af", component.max_alt_af.as_ref().ok_or_else(|| config_error(general, "population_snvs.generate.populations.spectrum.max_alt_af", "interval exige max_alt_af"))?)?,
                                    weight,
                                    sampling: match component.sampling.as_ref().map_or("", |sampling| sampling.value().as_str()) {
                                        "linear_uniform" => FrequencySampling::LinearUniform,
                                        "log_uniform" => FrequencySampling::LogUniform,
                                        _ => return Err(config_error(general, "population_snvs.generate.populations.spectrum.sampling", "interval exige linear_uniform o log_uniform")),
                                    },
                                }),
                                _ => Err(config_error(general, "population_snvs.generate.populations.spectrum.kind", "debe ser fixed o interval")),
                            }
                        }).collect::<Result<Vec<_>, Diagnostic>>()?,
                    }),
                    _ => return Err(config_error(general, "population_snvs.generate.populations.frequency_model", "modelo de frecuencia desconocido")),
                };
                Ok(GeneratedPopulationFrequency {
                    population_id: population.population_id.value().clone(),
                    min_alt_af: population.min_alt_af.as_ref().map(|value| decimal(
                        "population_snvs.generate.populations.min_alt_af",
                        value,
                    )).transpose()?.unwrap_or(0.0),
                    max_alt_af: population.max_alt_af.as_ref().map(|value| decimal(
                        "population_snvs.generate.populations.max_alt_af",
                        value,
                    )).transpose()?.unwrap_or(1.0),
                    spectrum,
                })
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?,
    }}))
}

fn frequency_assignment_mode(value: &str) -> Result<FrequencyAssignmentMode, Diagnostic> {
    match value {
        "stochastic" => Ok(FrequencyAssignmentMode::Stochastic),
        "exact" => Ok(FrequencyAssignmentMode::Exact),
        _ => Err(Diagnostic::configuration(
            "GPS-CONF-008",
            "El modo de espectro de frecuencias no es válido.",
        )),
    }
}

fn human_1000g_spectrum(mode: &str) -> Result<GeneratedFrequencySpectrum, Diagnostic> {
    Ok(GeneratedFrequencySpectrum {
        mode: frequency_assignment_mode(mode)?,
        components: vec![
            FrequencySpectrumComponent::Interval {
                min: 1.0 / 5008.0,
                max: 0.005,
                weight: 64.0,
                sampling: FrequencySampling::LogUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.005,
                max: 0.05,
                weight: 12.0,
                sampling: FrequencySampling::LogUniform,
            },
            FrequencySpectrumComponent::Interval {
                min: 0.05,
                max: 1.0,
                weight: 8.0,
                sampling: FrequencySampling::LinearUniform,
            },
        ],
    })
}

fn effective_normalization_tolerance(general: &DeclaredGeneralConfig) -> Result<f64, Diagnostic> {
    const DEFAULT: f64 = 1e-9;
    const MAXIMUM: f64 = 1e-6;
    let Some(declared) = general.normalization_tolerance() else {
        return Ok(DEFAULT);
    };
    let value = declared
        .literal()
        .parse::<f64>()
        .map_err(|_| config_error(general, "normalization_tolerance", "debe ser un f64 válido"))?;
    if !value.is_finite() || !(0.0..=MAXIMUM).contains(&value) {
        return Err(config_error(
            general,
            "normalization_tolerance",
            "debe ser finita y pertenecer a [0,1e-6]",
        ));
    }
    Ok(value)
}

fn effective_ploidy(general: &DeclaredGeneralConfig) -> Result<Ploidy, Diagnostic> {
    let Some(declared) = general.ploidy() else {
        return Ok(Ploidy::DIPLOID);
    };
    let value = *declared.value();
    let value = u8::try_from(value).map_err(|_| {
        config_error(
            general,
            "ploidy",
            format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
        )
    })?;
    Ploidy::new(value).map_err(|_| {
        config_error(
            general,
            "ploidy",
            format!("debe pertenecer a {}..={}", Ploidy::MIN, Ploidy::MAX),
        )
    })
}

fn normalize_composition(
    general: &DeclaredGeneralConfig,
    declared: &[crate::configuration::DeclaredCompositionValue],
) -> Result<[f64; 4], Diagnostic> {
    if declared.len() != 4 {
        return Err(config_error(
            general,
            "reference.composition",
            "se exigen A, C, G y T exactamente una vez",
        ));
    }
    let mut parsed = [0.0_f64; 4];
    let mut seen = [false; 4];
    for item in declared {
        let index = match item.base.value().as_str() {
            "A" => 0,
            "C" => 1,
            "G" => 2,
            "T" => 3,
            _ => {
                return Err(config_error(
                    general,
                    "reference.composition",
                    format!("base no admitida: {}", item.base.value()),
                ));
            }
        };
        if std::mem::replace(&mut seen[index], true) {
            return Err(config_error(
                general,
                "reference.composition",
                format!("base repetida: {}", item.base.value()),
            ));
        }
        let value = item.value.literal().parse::<f64>().map_err(|_| {
            config_error(
                general,
                "reference.composition",
                format!("valor decimal inválido para {}", item.base.value()),
            )
        })?;
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(config_error(
                general,
                "reference.composition",
                format!("{} debe ser finito y pertenecer a [0,1]", item.base.value()),
            ));
        }
        parsed[index] = value;
    }
    normalize_composition_probabilities(parsed)
        .map_err(|detail| config_error(general, "reference.composition", detail))
}

fn parse_unknown_base_fraction(
    general: &DeclaredGeneralConfig,
    declared: &Option<crate::configuration::DeclaredNumber>,
) -> Result<f64, Diagnostic> {
    let Some(declared) = declared else {
        return Ok(0.0);
    };
    let value = declared.literal().parse::<f64>().map_err(|_| {
        config_error(
            general,
            "reference.composition.N",
            "valor decimal inválido para N",
        )
    })?;
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(config_error(
            general,
            "reference.composition.N",
            "N debe ser finito y pertenecer a [0,1]",
        ));
    }
    Ok(value)
}

/// Resuelve los loci de repetición declarados.
///
/// # El tope se comprueba **antes** de materializar nada
///
/// Un alelo que rebasa el tope se rechaza **sin construir su secuencia**. La
/// versión anterior materializaba primero y comparaba después, y un error de
/// dedo de nueve cifras consumía gigabytes y decenas de segundos antes de
/// emitir el rechazo correcto —o no llegaba a emitirlo—. El tope existe
/// justamente «para que un error de dedo no genere un tracto desmesurado», de
/// modo que materializarlo primero anulaba su razón de ser.
///
/// El diagnóstico nombra **el locus y el alelo**, no el individuo: con
/// frecuencias (decisión S-04) el alelo se conoce al leer la configuración,
/// antes de que exista ningún individuo al que asignárselo.
///
/// # Los mínimos de admisión se aplican también al alelo declarado
///
/// La guarda de admisión sobre la referencia real —que compara el tracto que
/// hay en el genoma contra estos mismos mínimos— **todavía no existe**. Lo que
/// sí existe, y se aplica aquí, es la misma comparación sobre el alelo que el
/// usuario declara. Publicar los mínimos sin aplicarlos a nada sería publicar
/// el cartel de una guarda inexistente: quien leyera la configuración efectiva
/// concluiría que un alelo de pureza cero pasó un umbral de 0,8.
///
/// Y es la comprobación correcta por sí misma: un alelo que no contiene ni una
/// copia de su propio motivo no es un alelo de ese locus venga de donde venga, y
/// además haría fallar más adelante al paso de unidades sobre la tirada más
/// larga.
fn effective_repeat_loci(
    general: &DeclaredGeneralConfig,
    sequences: &[EffectiveSequence],
    tolerance: f64,
    admission: EffectiveRepeatAdmission,
) -> Result<Vec<EffectiveRepeatLocus>, Diagnostic> {
    let reference_id = ReferenceId::parse(general.reference().id.value().as_str(), "execution")
        .map_err(|diagnostic| {
            diagnostic.with("path", general.source().path().display().to_string())
        })?;
    let mut resolved: Vec<EffectiveRepeatLocus> =
        Vec::with_capacity(general.reference().repeat_loci.len());
    for locus in &general.reference().repeat_loci {
        // El identificador del locus es **tipado**: lo declara el usuario y se
        // valida con la misma sintaxis que el resto de identificadores
        // aportados. Sin esto, un identificador con espacios o acentos se
        // admitía, se ejecutaba y se publicaba, y solo el validador externo de
        // la configuración efectiva lo cazaba después, contra su propio patrón.
        let locus_id = RepeatLocusId::parse(locus.id.value().as_str(), "repeat_locus").map_err(
            |diagnostic| {
                diagnostic
                    .with("path", general.source().path().display().to_string())
                    .with("location", locus.id.location().to_string())
            },
        )?;
        let mut motifs = Vec::with_capacity(locus.motifs.len());
        for declared in &locus.motifs {
            motifs.push(
                RepeatMotif::new(declared.value().as_bytes()).map_err(|error| {
                    repeat_error(
                        general,
                        "reference.repeat_loci.motif",
                        locus.id.value(),
                        error.detail(),
                    )
                })?,
            );
        }
        let sequence = sequences
            .iter()
            .find(|sequence| sequence.id() == locus.sequence.value())
            .ok_or_else(|| {
                repeat_error(
                    general,
                    "reference.repeat_loci.sequence",
                    locus.id.value(),
                    format!("secuencia desconocida: {}", locus.sequence.value()),
                )
            })?;
        let sequence_id =
            SequenceId::parse(sequence.id(), reference_id.as_str()).map_err(|diagnostic| {
                diagnostic.with("path", general.source().path().display().to_string())
            })?;
        // El intervalo pasa por el servicio común de coordenadas, que es quien
        // convierte de base uno cerrada a base cero semiabierta y quien
        // comprueba que cabe en la secuencia. Aquí no se resta ni se suma uno.
        let closed =
            ClosedInterval1::new(*locus.start.value(), *locus.end.value()).map_err(|error| {
                repeat_error(
                    general,
                    "reference.repeat_loci.start,reference.repeat_loci.end",
                    locus.id.value(),
                    format!("intervalo declarado inválido: {error:?}"),
                )
            })?;
        let interval = GenomicInterval::from_closed_1(
            reference_id.clone(),
            sequence_id,
            closed,
            sequence.length(),
        )
        .map_err(|error| {
            repeat_error(
                general,
                "reference.repeat_loci.start,reference.repeat_loci.end",
                locus.id.value(),
                format!(
                    "el intervalo {}..{} no cabe en {} ({} bases): {error:?}",
                    locus.start.value(),
                    locus.end.value(),
                    sequence.id(),
                    sequence.length()
                ),
            )
        })?;
        // Dos loci que se solapan describirían el mismo tracto dos veces y
        // producirían dos alelos sobre las mismas bases.
        if let Some(previous) = resolved.iter().find(|other| {
            other.declaration().interval().sequence_id() == interval.sequence_id()
                && other.declaration().interval().start() < interval.end()
                && interval.start() < other.declaration().interval().end()
        }) {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.start,reference.repeat_loci.end",
                locus.id.value(),
                format!("su intervalo se solapa con el del locus {}", previous.id()),
            ));
        }
        // La guarda de intervalo contra motivo vive en el tipo del dominio, no
        // reescrita aquí: una sola representación del mismo concepto.
        let declaration =
            TandemRepeatLocus::new(locus_id, interval, motifs.clone()).map_err(|error| {
                repeat_error(
                    general,
                    "reference.repeat_loci.start,reference.repeat_loci.end",
                    locus.id.value(),
                    error.detail(),
                )
            })?;
        // Los mínimos del locus sustituyen a los generales, que quedan como
        // valor de partida. El mismo número no significa lo mismo según el
        // motivo, y un umbral único o deja pasar basura en los motivos cortos o
        // rechaza loci legítimos en los largos.
        let locus_admission = resolve_admission(general, locus.admission.as_ref(), admission)
            .map_err(|diagnostic| diagnostic.with("fields", "reference.repeat_loci.min_purity"))?;
        let max_units = locus
            .max_units
            .as_ref()
            .map_or(DEFAULT_REPEAT_MAX_UNITS, |declared| *declared.value());
        if max_units == 0 {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.max_units",
                locus.id.value(),
                "el tope debe ser un entero positivo en unidades",
            ));
        }

        let mut declared_weights = Vec::with_capacity(locus.alleles.len());
        let mut tracts = Vec::with_capacity(locus.alleles.len());
        for (index, allele) in locus.alleles.iter().enumerate() {
            // ANTES de materializar: el tope se decide sobre lo declarado.
            if let Some(units) = declared_unit_count(&motifs, allele)
                && units > max_units as f64
            {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles",
                    locus.id.value(),
                    format!(
                        "el alelo {index} mide {units} unidades y rebasa el tope de {max_units}"
                    ),
                ));
            }
            let tract = materialize_declared_allele(general, locus, &motifs, index, allele)?;
            let annotation = tract.annotation();
            // Red por si la forma declarada no permitía contar las unidades de
            // antemano. Hoy siempre permite, así que esta rama no se ejercita;
            // se conserva porque el tope es un contrato y no puede depender de
            // que la forma declarada sea una de las dos de hoy.
            if annotation.units_as_f64() > max_units as f64 {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles",
                    locus.id.value(),
                    format!(
                        "el alelo {index} mide {} unidades y rebasa el tope de {max_units}",
                        annotation.units_as_f64()
                    ),
                ));
            }
            if annotation.purity_as_f64() < locus_admission.min_purity() {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles",
                    locus.id.value(),
                    format!(
                        "el alelo {index} tiene pureza {} y el mínimo declarado es {}",
                        annotation.purity_as_f64(),
                        locus_admission.min_purity()
                    ),
                ));
            }
            if annotation.longest_run_units() < locus_admission.min_longest_run_units() {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles",
                    locus.id.value(),
                    format!(
                        "el alelo {index} tiene una tirada más larga de {} unidades y el mínimo declarado es {}",
                        annotation.longest_run_units(),
                        locus_admission.min_longest_run_units()
                    ),
                ));
            }
            let frequency = allele.frequency.literal().parse::<f64>().map_err(|_| {
                repeat_error(
                    general,
                    "reference.repeat_loci.alleles.frequency",
                    locus.id.value(),
                    format!("el alelo {index} declara una frecuencia decimal inválida"),
                )
            })?;
            if !frequency.is_finite() || !(0.0..=1.0).contains(&frequency) {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles.frequency",
                    locus.id.value(),
                    format!("el alelo {index} declara una frecuencia fuera de [0,1]"),
                ));
            }
            declared_weights.push(frequency);
            tracts.push(tract);
        }
        // DTC-049: suma compensada con **el mismo** algoritmo que consume el
        // sorteo categórico, y auditoría recomputada después de normalizar.
        let original_sum = neumaier_sum(declared_weights.iter().copied());
        if !original_sum.is_finite()
            || original_sum <= 0.0
            || (original_sum - 1.0).abs() > tolerance
        {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.alleles.frequency",
                locus.id.value(),
                format!(
                    "las frecuencias suman {original_sum} y se apartan de 1 más que {tolerance}"
                ),
            ));
        }
        let effective_weights: Vec<f64> = declared_weights
            .iter()
            .map(|weight| weight / original_sum)
            .collect();
        // Recomputar sobre lo normalizado, no dar por hecho que suma 1: es lo
        // que el resto del proyecto hace y lo que hace publicable la auditoría.
        let effective_sum = neumaier_sum(effective_weights.iter().copied());
        let normalization = RepeatFrequencyNormalization {
            tolerance,
            original_sum,
            original_delta: (original_sum - 1.0).abs(),
            effective_sum,
            effective_delta: (effective_sum - 1.0).abs(),
        };
        let alleles = tracts
            .into_iter()
            .zip(declared_weights.iter().copied())
            .zip(effective_weights)
            .map(
                |((tract, declared_frequency), frequency)| EffectiveRepeatAllele {
                    tract,
                    declared_frequency,
                    frequency,
                },
            )
            .collect();
        resolved.push(EffectiveRepeatLocus {
            declaration,
            max_units,
            admission: locus_admission,
            alleles,
            normalization,
        });
    }
    Ok(resolved)
}

/// Número **total** de unidades que la declaración implica, sin materializar
/// nada. Es lo que el tope compara, y sumar los bloques es lo correcto: lo que
/// el tope protege es el tamaño del tracto entero.
///
/// Devuelve `None` cuando la declaración no es una de las dos formas válidas;
/// en ese caso el diagnóstico correcto lo produce la materialización.
fn declared_unit_count(
    motifs: &[RepeatMotif],
    allele: &crate::configuration::DeclaredRepeatAllele,
) -> Option<f64> {
    match (allele.units.is_empty(), &allele.sequence) {
        (false, None) if allele.units.len() == motifs.len() => {
            Some(allele.units.iter().map(|units| *units.value() as f64).sum())
        }
        (true, Some(sequence)) if motifs.len() == 1 => {
            Some(sequence.value().len() as f64 / motifs[0].length() as f64)
        }
        _ => None,
    }
}

/// Materializa la secuencia del alelo desde su declaración compacta.
///
/// Las dos formas son excluyentes: **o** el número de unidades de cada bloque
/// con sus interrupciones, **o** la secuencia exacta. Una interrupción
/// **sustituye** las bases que ocupan su desplazamiento —absoluto dentro del
/// tracto— y conserva la longitud, de modo que la estructura en bloques no se
/// mueve.
///
/// La secuencia exacta solo se admite con **un** motivo: unas bases sueltas no
/// pueden decir dónde acaba un bloque y empieza el siguiente, y esa es
/// justamente la información que la estructura aporta.
fn materialize_declared_allele(
    general: &DeclaredGeneralConfig,
    locus: &crate::configuration::DeclaredRepeatLocus,
    motifs: &[RepeatMotif],
    index: usize,
    allele: &crate::configuration::DeclaredRepeatAllele,
) -> Result<TandemRepeatTract, Diagnostic> {
    match (allele.units.is_empty(), &allele.sequence) {
        (false, Some(_)) | (true, None) => Err(repeat_error(
            general,
            "reference.repeat_loci.alleles",
            locus.id.value(),
            format!("el alelo {index} debe declarar units o sequence, exactamente uno"),
        )),
        (true, Some(sequence)) => {
            if motifs.len() != 1 {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles.sequence",
                    locus.id.value(),
                    format!(
                        "el alelo {index} declara su secuencia exacta y el locus tiene {} motivos: una secuencia no puede decir dónde acaba un bloque y empieza el siguiente",
                        motifs.len()
                    ),
                ));
            }
            if !allele.interruptions.is_empty() {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.alleles.interruptions",
                    locus.id.value(),
                    format!(
                        "el alelo {index} declara su secuencia exacta: las interrupciones ya están en ella"
                    ),
                ));
            }
            TandemRepeatTract::new(motifs[0].clone(), sequence.value().as_bytes()).map_err(
                |error| {
                    repeat_error(
                        general,
                        "reference.repeat_loci.alleles.sequence",
                        locus.id.value(),
                        format!("el alelo {index}: {}", error.detail()),
                    )
                },
            )
        }
        (false, None) => {
            let specification: Vec<(RepeatMotif, u64)> = motifs
                .iter()
                .cloned()
                .zip(allele.units.iter().map(|units| *units.value()))
                .collect();
            let tract = TandemRepeatTract::pure_blocks(&specification).map_err(|error| {
                repeat_error(
                    general,
                    "reference.repeat_loci.alleles.units",
                    locus.id.value(),
                    format!("el alelo {index}: {}", error.detail()),
                )
            })?;
            if allele.interruptions.is_empty() {
                return Ok(tract);
            }
            let mut bases = tract.bases().to_vec();
            for interruption in &allele.interruptions {
                let offset = usize::try_from(*interruption.offset.value()).map_err(|_| {
                    repeat_error(
                        general,
                        "reference.repeat_loci.alleles.interruptions.offset",
                        locus.id.value(),
                        format!("el alelo {index} declara un desplazamiento inalcanzable"),
                    )
                })?;
                let replacement = interruption.bases.value().as_bytes();
                if replacement.is_empty() || offset + replacement.len() > bases.len() {
                    return Err(repeat_error(
                        general,
                        "reference.repeat_loci.alleles.interruptions",
                        locus.id.value(),
                        format!(
                            "el alelo {index} tiene una interrupción en {offset} que no cabe en sus {} bases",
                            bases.len()
                        ),
                    ));
                }
                bases[offset..offset + replacement.len()].copy_from_slice(replacement);
            }
            // La sustitución conserva la longitud, así que los bloques del
            // tracto puro siguen valiendo tal cual.
            TandemRepeatTract::from_blocks(tract.blocks().to_vec(), bases).map_err(|error| {
                repeat_error(
                    general,
                    "reference.repeat_loci.alleles.interruptions",
                    locus.id.value(),
                    format!("el alelo {index}: {}", error.detail()),
                )
            })
        }
    }
}

/// Resuelve unos mínimos declarados contra los que sirven de valor de partida.
///
/// Es la misma función para los generales —cuyo valor de partida son los de
/// omisión— y para los de cada locus, cuyo valor de partida son los generales.
fn resolve_admission(
    general: &DeclaredGeneralConfig,
    declared: Option<&crate::configuration::DeclaredRepeatAdmission>,
    fallback: EffectiveRepeatAdmission,
) -> Result<EffectiveRepeatAdmission, Diagnostic> {
    let Some(declared) = declared else {
        return Ok(fallback);
    };
    let min_purity = match &declared.min_purity {
        None => fallback.min_purity(),
        Some(value) => {
            let parsed = value
                .literal()
                .parse::<f64>()
                .map_err(|_| config_error(general, "min_purity", "debe ser un decimal válido"))?;
            if !parsed.is_finite() || !(0.0..=1.0).contains(&parsed) {
                return Err(config_error(
                    general,
                    "min_purity",
                    "debe ser finito y pertenecer a [0,1]",
                ));
            }
            parsed
        }
    };
    let min_longest_run_units = declared
        .min_longest_run_units
        .as_ref()
        .map_or(fallback.min_longest_run_units(), |value| *value.value());
    if min_longest_run_units == 0 {
        return Err(config_error(
            general,
            "min_longest_run_units",
            "debe ser un entero positivo",
        ));
    }
    Ok(EffectiveRepeatAdmission {
        min_purity,
        min_longest_run_units,
    })
}

/// Mínimos **generales**, que sirven de valor de partida para los loci que no
/// declaren los suyos.
fn effective_repeat_admission(
    general: &DeclaredGeneralConfig,
) -> Result<EffectiveRepeatAdmission, Diagnostic> {
    resolve_admission(
        general,
        general.reference().repeat_admission.as_ref(),
        EffectiveRepeatAdmission::DEFAULT,
    )
    .map_err(|diagnostic| diagnostic.with("fields", "reference.repeat_admission"))
}

fn repeat_error(
    general: &DeclaredGeneralConfig,
    fields: &str,
    locus_id: &str,
    detail: impl Into<String>,
) -> Diagnostic {
    config_error(
        general,
        fields,
        format!("locus {locus_id}: {}", detail.into()),
    )
}

fn resolve_link(general_path: &Path, declared: &str) -> PathBuf {
    // Los campos de ruta que introduce Marimo son texto. Recortamos únicamente
    // espacio accidental en los extremos antes de decidir si es absoluta; los
    // espacios internos se conservan como parte legítima del nombre.
    let declared = Path::new(declared.trim());
    if declared.is_absolute() {
        declared.to_path_buf()
    } else {
        general_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(declared)
    }
}

fn input_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

fn config_error(
    general: &DeclaredGeneralConfig,
    fields: &str,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", general.source().path().display().to_string())
    .with("fields", fields)
    .with("detail", detail.into())
}

// ---------------------------------------------------------------------------
// F5.6.1 — cadena productiva de `run`
// ---------------------------------------------------------------------------

/// Lo que una ejecución F5 publicada deja al llamador.
///
/// No es un artefacto: el paquete ya está en disco y este valor solo nombra
/// dónde quedó y qué contaron sus auditores, para que la frontera CLI y las
/// baterías puedan contrastarlo sin releer el destino.
#[derive(Clone, Copy, Debug)]
pub struct F5ExecutionOutcome<'a> {
    destination: &'a Path,
    closure: crate::f5_closure::F5ClosureReport,
}

impl F5ExecutionOutcome<'_> {
    /// Destino final, ya renombrado.
    #[must_use]
    pub const fn destination(&self) -> &Path {
        self.destination
    }
    /// Recuentos de la clausura de identidad sobre ese destino.
    #[must_use]
    pub const fn closure(&self) -> crate::f5_closure::F5ClosureReport {
        self.closure
    }
}

/// Material canónico del `execution_id` de F5.
///
/// F1 a F4 lo derivan del digest del **único** documento de autoría. F5 tiene
/// dos y el general **no** declara el hash del escenario, de modo que usar
/// solo el general daría la misma identidad de ejecución a dos escenarios
/// distintos: `block_id` y `site_id` la incorporan, así que dos escenarios
/// distintos compartirían identificadores de bloque y de sitio. El material
/// es, por tanto, la concatenación de los dos SHA-256 observados, en orden de
/// rol —general y después escenario—, calculados sobre el mismo handle que
/// los leyó y sin volver a abrirlos. Los dos cubren todo lo declarado: el
/// general fija referencia y su hash, y el escenario fija perfil, máscara,
/// límites, objetivo, subpoblaciones y sustituciones con sus hashes.
///
/// El orden **no** es indiferente: intercambiarlo cambia todos los
/// `execution_id`, `block_id` y `site_id` de F5. `f5_execution_id` lo fija con
/// un vector y una rederivación independiente desde semilla y digest.
fn f5_execution_material(documents: &F5BootstrapDocuments) -> String {
    let mut material = String::with_capacity(128);
    material.push_str(documents.general().sha256());
    material.push_str(documents.scenario().sha256());
    material
}

/// Identidad de ejecución F5, derivada del mismo camino que F1 a F4.
///
/// # Errors
///
/// `GPS-CONF-017` si la semilla efectiva no es canónica y el diagnóstico de
/// `configuration_digest` o de la derivación si el material desborda.
pub fn f5_execution_id(config: &EffectiveF5Config) -> Result<ExecutionId, Diagnostic> {
    let material = f5_execution_material(config.bootstrap_documents());
    let root_seed = crate::randomness::RootSeed::parse("seed", config.seed())?;
    crate::identity::derive_execution_id(
        root_seed,
        &crate::integrity::configuration_digest(material.as_bytes())?,
    )
}

/// Ejecuta la cadena productiva F5 y publica su paquete en `destination`.
///
/// # Orden, y por qué el provisional llega antes que el layout
///
/// DTC-045 distingue el rechazo de pre-simulación —que no publica destino— del
/// fallo controlado **posterior al comienzo**, que «publica atómicamente el
/// directorio final como fallido». La simulación F5 comienza cuando se abre el
/// primer ámbito RNG derivado del `execution_id`, y eso ocurre ya en el
/// layout. Por eso el provisional se crea **antes** que él: layout, búsqueda
/// haplotípica, asignación fundadora, materialización de identidades y
/// publicación tienen las cinco un productor de paquete `failed`. Definir
/// «comenzar» como «cuando este código abre el provisional» haría la regla
/// autocumplida y dependiente del orden de las líneas.
///
/// La admisión física **no** se repite: llega ya sellada en `config`, que el
/// dispatcher materializó una sola vez.
///
/// # Obligación de proceso
///
/// Publica **exclusivamente** por [`crate::f5_closure::seal_and_publish_f5_success`]
/// y [`crate::f5_closure::seal_and_publish_f5_failure`]. Llamar a
/// `publish_f5_success` o a `publish_f5_failure` daría un destino **sin
/// testigo** y dejaría la clausura de F5.5.1 inaplicable sin ninguna señal.
///
/// # Auditores
///
/// `audit_f5_package` entra como el recuento MAN-70 que la publicación exige,
/// sobre la raíz de carga útil que `publish` renombra **sin tocar** al
/// destino; `audit_f5_identity_closure` corre después sobre el destino ya
/// renombrado, que es el que su testigo acredita. En la mitad «fallo» corren
/// los dos gemelos sobre el destino ya renombrado: `audit_f5_failure_closure`
/// acredita el testigo y `audit_f5_failure_package` relee el paquete.
///
/// # Errors
///
/// Cualquier diagnóstico de las cinco etapas, de la publicación o de los
/// auditores. Un fallo controlado publica el paquete `failed` y **aun así**
/// devuelve el diagnóstico de su informe, para que la frontera CLI conserve su
/// código de salida.
pub fn execute_validated_f5_config<'a>(
    config: &EffectiveF5Config,
    destination: &'a Path,
) -> Result<F5ExecutionOutcome<'a>, Diagnostic> {
    execute_validated_f5_config_with_recount(config, destination, |payload, declaration| {
        crate::f5_oracle::audit_f5_package(payload, declaration)
            .map(|report| report.initialization_summaries().to_vec())
            .map_err(crate::f5_oracle::F5OracleRejection::into_diagnostic)
    })
}

/// Audita el paquete `failed` ya renombrado con el oráculo independiente.
///
/// Es la mitad que F5.6.1 no pudo componer. `audit_f5_failure_package` exige
/// una [`crate::f5_oracle::F5FailureDeclaration`] que el orquestador no
/// conserva, y reconstruirla desde el informe en memoria haría la auditoría
/// autocertificante: el productor declararía lo que el productor escribió.
///
/// Tomarla entera del manifiesto **tampoco** basta, y la Revisión F5.6 lo
/// demostró: `diagnostic_document` y `manifest_document` salen los dos de
/// `common_json`, de modo que los nueve términos comunes son la misma
/// serialización emitida dos veces y el auditor los compararía consigo mismos.
/// Los dos que tienen fuente **contractual** —`scenario_id` y los siete
/// `limits`— se derivan por eso de la autoría admitida y de la admisión
/// sellada, no del manifiesto: así el auditor contrasta **contrato contra
/// bytes**. Los demás términos —código, etapa, contadores, mejor intento— son
/// medida del productor y no tienen otra fuente que comparar; lo que el
/// auditor exige de ellos es coherencia entre los dos documentos, el árbol
/// enumerado y el rehash del enlace.
///
/// # Instrumental
///
/// Es `pub` y `doc(hidden)` porque es la **fuente de verdad de la
/// composición**: una batería que construya la declaración a mano ejercita el
/// oráculo, no esto, y un refactor que devolviera `scenario_id` al manifiesto
/// no despertaría a nadie. La Revisión F5.6 lo comprobó y por eso se exporta.
///
/// # Errors
///
/// El diagnóstico del oráculo si el paquete publicado no sostiene lo que su
/// autoría admitida y su propio manifiesto declaran.
#[doc(hidden)]
pub fn audit_published_f5_failure(
    destination: &Path,
    config: &EffectiveF5Config,
    sealed: &crate::f5_preflight::SealedF5Admission,
) -> Result<(), Diagnostic> {
    use crate::f5_oracle::{F5DeclaredLink, F5FailureDeclaration, audit_f5_failure_package};
    use crate::f5_publication::F5_MANIFEST_PATH;

    let path = destination.join(F5_MANIFEST_PATH);
    let bytes = std::fs::read(&path).map_err(|error| {
        Diagnostic::f5_configuration(
            path.to_string_lossy().as_ref(),
            "manifest",
            format!("no se pudo releer el manifiesto publicado: {error}"),
        )
    })?;
    let manifest: serde_json::Value = serde_json::from_slice(&bytes).map_err(|error| {
        Diagnostic::f5_configuration(
            path.to_string_lossy().as_ref(),
            "manifest",
            format!("el manifiesto publicado no es JSON: {error}"),
        )
    })?;
    let failure = &manifest["f5"]["f5_failure"];
    let invalid = |field: &str| {
        Diagnostic::f5_configuration(
            path.to_string_lossy().as_ref(),
            field,
            "el manifiesto publicado no describe su propio fallo",
        )
    };
    let text = |value: &serde_json::Value, field: &str| {
        value
            .as_str()
            .map(str::to_owned)
            .ok_or_else(|| invalid(field))
    };
    let numbers = |object: &serde_json::Value, keys: &[&str], field: &'static str| {
        keys.iter()
            .map(|key| object[*key].as_u64().ok_or_else(|| invalid(field)))
            .collect::<Result<Vec<_>, Diagnostic>>()
    };
    // Los siete límites salen de la **autoría admitida**, con la misma regla
    // que el productor usó al publicarlos. Un manifiesto que declare otros
    // muere en el auditor, que los contrasta contra los bytes del documento.
    let contractual_limits = crate::f5_publication::limits_json(config, sealed)?;
    let limits = numbers(
        &contractual_limits,
        &[
            "max_layout_attempts",
            "max_haplotype_attempts",
            "max_swap_proposals",
            "max_haplotype_cells",
            "max_operation_units",
            "max_working_memory_bytes",
            "max_output_bytes",
        ],
        "limits",
    )?;
    let counters = numbers(
        &failure["counters"],
        &[
            "layout_attempts_used",
            "haplotype_attempts_used",
            "swap_proposals",
            "swap_accepted",
        ],
        "counters",
    )?;
    // `population_id` **sí** es acotable por el contrato, a diferencia del
    // código, la etapa y los contadores: si el fallo nombra una subpoblación,
    // tiene que ser una de las admitidas. No se deriva —cuál falló es medida
    // del productor— pero sí se exige **pertenencia** al conjunto sellado.
    let population_id = failure["population_id"].as_str().map(str::to_owned);
    if let Some(named) = population_id.as_deref()
        && !config
            .declared_scenario()
            .subpopulations
            .iter()
            .any(|subpopulation| subpopulation.population_id.value() == named)
    {
        return Err(Diagnostic::f5_configuration(
            path.to_string_lossy().as_ref(),
            "population_id",
            format!("{named} no es ninguna de las subpoblaciones admitidas"),
        ));
    }
    let link = &failure["artifacts"][0];
    let declaration = F5FailureDeclaration {
        code: text(&failure["code"], "code")?,
        stage: text(&failure["stage"], "stage")?,
        // El escenario es el de la autoría admitida, no el que el manifiesto
        // diga que es: es el término que la revisión fabricó sin que nadie lo
        // notara.
        scenario_id: config.scenario_id().to_owned(),
        population_id,
        block_id: failure["block_id"].as_str().map(str::to_owned),
        limits: limits.try_into().map_err(|_| invalid("limits"))?,
        counters: counters.try_into().map_err(|_| invalid("counters"))?,
        best_layout: failure["best_layout"].clone(),
        best_attempt: failure["best_attempt"].clone(),
        link: F5DeclaredLink {
            artifact_id: text(&link["artifact_id"], "artifacts")?,
            artifact_kind: text(&link["artifact_kind"], "artifacts")?,
            schema_version: text(&link["schema_version"], "artifacts")?,
            path: text(&link["path"], "artifacts")?,
            size_bytes: link["size_bytes"]
                .as_u64()
                .ok_or_else(|| invalid("artifacts"))?,
            sha256: text(&link["sha256"], "artifacts")?,
        },
    };
    audit_f5_failure_package(destination, &declaration)
        .map_err(crate::f5_oracle::F5OracleRejection::into_diagnostic)
}

/// Variante instrumental con el recuento MAN-70 inyectado.
///
/// Producción usa siempre `audit_f5_package`; esta forma existe para que una
/// batería pueda recorrer la rama de **publicación fallida**, que de otro modo
/// solo se alcanzaría corrompiendo el provisional desde fuera del proceso.
///
/// # Errors
///
/// Los mismos que [`execute_validated_f5_config`].
#[doc(hidden)]
pub fn execute_validated_f5_config_with_recount<'a, R>(
    config: &EffectiveF5Config,
    destination: &'a Path,
    recount: R,
) -> Result<F5ExecutionOutcome<'a>, Diagnostic>
where
    R: FnOnce(
        &Path,
        &crate::f5_oracle::F5PackageDeclaration,
    ) -> Result<Vec<crate::f5_oracle::F5InitializationSummary>, Diagnostic>,
{
    execute_validated_f5_config_with_probes(config, destination, recount, Ok)
}

/// Variante instrumental con el recuento **y** la clausura de identidad
/// inyectables.
///
/// La clausura real corre siempre; la sonda solo puede convertir su **éxito**
/// en un rechazo, que es la única forma de recorrer desde una batería la
/// degradación del destino sin corromper el paquete desde fuera del proceso.
/// Un defecto real del productor llega por el mismo camino.
///
/// # Errors
///
/// Los mismos que [`execute_validated_f5_config`].
#[doc(hidden)]
pub fn execute_validated_f5_config_with_probes<'a, R, C>(
    config: &EffectiveF5Config,
    destination: &'a Path,
    recount: R,
    closure_probe: C,
) -> Result<F5ExecutionOutcome<'a>, Diagnostic>
where
    R: FnOnce(
        &Path,
        &crate::f5_oracle::F5PackageDeclaration,
    ) -> Result<Vec<crate::f5_oracle::F5InitializationSummary>, Diagnostic>,
    C: FnOnce(
        crate::f5_closure::F5ClosureReport,
    ) -> Result<crate::f5_closure::F5ClosureReport, Diagnostic>,
{
    use crate::f5_closure::{
        F5ClosureRejection, audit_f5_failure_closure, audit_f5_identity_closure,
        seal_and_publish_f5_failure, seal_and_publish_f5_success,
    };
    use crate::f5_founders::{
        assign_f5_founders, empty_f5_founder_identities, materialize_f5_founder_identities,
    };
    use crate::f5_layout::{F5LayoutFailure, F5LayoutRequest, build_f5_physical_layout};
    use crate::f5_preflight::{F5OtherResources, recover_sealed_f5_admission};
    use crate::f5_publication::{F5FailureReport, F5FailureStage, F5PublicationRequest};
    use crate::f5_search::{F5SearchFailure, search_f5_haplotypes};
    use crate::publication::ExecutionWorkspace;

    let sealed = recover_sealed_f5_admission(config, F5OtherResources::default())?;
    let execution_id = f5_execution_id(config)?;
    let empty = empty_f5_founder_identities(config, &execution_id)?;

    // El destino que ya existe es un rechazo previo al comienzo: no se abre
    // provisional, no se consume RNG y no se publica nada.
    let mut workspace = Some(ExecutionWorkspace::create(destination)?);
    // Cada etapa entrega el provisional a la publicación del fallo, o lo
    // conserva para la siguiente. `take` hace imposible usarlo dos veces.
    let published_failure = |workspace: &mut Option<ExecutionWorkspace>,
                             report: &F5FailureReport|
     -> Diagnostic {
        let Some(provisional) = workspace.take() else {
            return report.diagnostic().clone();
        };
        match seal_and_publish_f5_failure(
            provisional,
            config,
            &sealed,
            &execution_id,
            report,
            &empty,
        ) {
            Ok(failed) => match audit_f5_failure_closure(failed.destination(), &failed, &empty) {
                Ok(_) => match audit_published_f5_failure(failed.destination(), config, &sealed) {
                    Ok(()) => report.diagnostic().clone(),
                    Err(diagnostic) => diagnostic,
                },
                Err(rejection) => rejection.into_diagnostic(),
            },
            Err(diagnostic) => diagnostic,
        }
    };

    let layout =
        match build_f5_physical_layout(F5LayoutRequest::new(config, &sealed, &execution_id)?) {
            Ok(layout) => layout,
            Err(failure) => {
                let report = F5FailureReport::from_layout(config, &failure);
                drop::<F5LayoutFailure>(failure);
                return Err(published_failure(&mut workspace, &report));
            }
        };
    let search = match search_f5_haplotypes(config, &sealed, &execution_id, &layout) {
        Ok(search) => search,
        Err(failure) => {
            let report = F5FailureReport::from_search(config, &layout, &failure);
            drop::<F5SearchFailure>(failure);
            return Err(published_failure(&mut workspace, &report));
        }
    };
    let assignment = match assign_f5_founders(config, &sealed, &search) {
        Ok(assignment) => assignment,
        Err(diagnostic) => {
            let report = F5FailureReport::from_diagnostic(
                config,
                Some(&layout),
                F5FailureStage::Founders,
                diagnostic,
            );
            return Err(published_failure(&mut workspace, &report));
        }
    };
    let mut identities = match materialize_f5_founder_identities(
        config,
        &sealed,
        &execution_id,
        &layout,
        &search,
        &assignment,
    ) {
        Ok(identities) => identities,
        Err(diagnostic) => {
            let report = F5FailureReport::from_diagnostic(
                config,
                Some(&layout),
                F5FailureStage::Founders,
                diagnostic,
            );
            return Err(published_failure(&mut workspace, &report));
        }
    };

    let provisional = workspace.take().ok_or_else(|| {
        Diagnostic::internal_invariant("f5_execution", "workspace_consumed_before_publication")
    })?;
    let request = F5PublicationRequest::new(
        config,
        &sealed,
        &execution_id,
        &layout,
        &search,
        &assignment,
    )?;
    let witness = match seal_and_publish_f5_success(provisional, &request, &mut identities, recount)
    {
        Ok(witness) => witness,
        Err(failure) => {
            let (returned, diagnostic) = failure.into_parts();
            // Sin provisional vivo no hay nada que publicar: es la
            // interrupción que DTC-045 reserva al `INCOMPLETE`.
            let Some(returned) = returned else {
                return Err(diagnostic);
            };
            workspace = Some(returned);
            let report = F5FailureReport::from_diagnostic(
                config,
                Some(&layout),
                F5FailureStage::Publication,
                diagnostic,
            );
            return Err(published_failure(&mut workspace, &report));
        }
    };
    let audited = audit_f5_identity_closure(destination, &witness, &identities)
        .map_err(F5ClosureRejection::into_diagnostic)
        .and_then(closure_probe);
    let closure = match audited {
        Ok(closure) => closure,
        Err(rejection) => return Err(degrade_rejected_destination(destination, rejection)),
    };
    Ok(F5ExecutionOutcome {
        destination,
        closure,
    })
}

/// Aparta el destino que la clausura acaba de rechazar y enriquece el motivo.
///
/// Sin esto, el paquete defectuoso se queda en el destino final con
/// `state="complete"`, sin marcador ni documento de fallo, y el único testigo
/// del rechazo es el código de salida del proceso. La mitad «éxito» de la
/// cadena ya convierte en paquete `failed` cualquier rechazo **anterior** al
/// renombrado; esta es la mitad que solo puede descubrirse después.
fn degrade_rejected_destination(destination: &Path, rejection: Diagnostic) -> Diagnostic {
    match crate::publication::demote_rejected_destination(destination) {
        Ok(demoted) => rejection.with("demoted_destination", demoted.display().to_string()),
        // No poder apartarlo es peor noticia que el rechazo: se dice entera.
        Err(failure) => rejection.with("demotion_failed", "true").with(
            "demotion_detail",
            failure
                .context()
                .get("detail")
                .cloned()
                .unwrap_or_else(|| failure.code().to_owned()),
        ),
    }
}

#[cfg(test)]
mod bounded_toml_tests {
    use std::io::{self, Read};
    use std::path::Path;
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use super::{F5_PRE_ADMISSION_TOML_LIMIT_BYTES, read_bounded_f5_toml_from_reader};
    use crate::f5_admission::InputRole;

    struct CountingReader {
        bytes: Vec<u8>,
        offset: usize,
        observed: Arc<AtomicUsize>,
    }

    impl Read for CountingReader {
        fn read(&mut self, destination: &mut [u8]) -> io::Result<usize> {
            let remaining = &self.bytes[self.offset..];
            let count = remaining.len().min(destination.len());
            destination[..count].copy_from_slice(&remaining[..count]);
            self.offset += count;
            self.observed.fetch_add(count, Ordering::Relaxed);
            Ok(count)
        }
    }

    #[test]
    fn bounded_toml_reader_observes_only_the_excess_byte_and_never_materializes_it() {
        let limit = usize::try_from(F5_PRE_ADMISSION_TOML_LIMIT_BYTES).expect("limit usize");
        let observed = Arc::new(AtomicUsize::new(0));
        let reader = CountingReader {
            bytes: vec![b'#'; limit + 4_096],
            offset: 0,
            observed: Arc::clone(&observed),
        };

        let error = read_bounded_f5_toml_from_reader(
            Path::new("unclassified.toml"),
            InputRole::GeneralToml,
            reader,
        )
        .expect_err("the byte after the inclusive limit rejects the input");

        assert_eq!(error.code(), "GPS-INPUT-001");
        assert_eq!(observed.load(Ordering::Relaxed), limit + 1);
        assert!(error.render_human().contains("byte 65537"));
        assert!(!error.context().contains_key("scenario_id"));
    }
}

#[cfg(test)]
mod profile_resource_tests {
    use super::{
        FrequencyAssignmentMode, FrequencySampling, FrequencySpectrumComponent,
        human_1000g_spectrum,
    };

    #[test]
    fn human_1000g_spectrum_matches_its_versioned_resource_literals() {
        let resource = include_str!("../resources/f2_2_6/human_1000g_global_sfs_approx_v1.tsv");
        let rows = resource
            .lines()
            .skip(1)
            .map(|line| line.split('\t').collect::<Vec<_>>())
            .collect::<Vec<_>>();
        assert_eq!(rows.len(), 3);
        let spectrum = human_1000g_spectrum("exact").unwrap();
        assert_eq!(spectrum.mode, FrequencyAssignmentMode::Exact);
        for (component, fields) in spectrum.components.iter().zip(rows) {
            assert_eq!(fields[0], "interval");
            let (min, max, weight, sampling) = match component {
                FrequencySpectrumComponent::Interval {
                    min,
                    max,
                    weight,
                    sampling,
                } => (min, max, weight, sampling),
                FrequencySpectrumComponent::Fixed { .. } => {
                    panic!("the SFS resource has only intervals")
                }
            };
            assert_eq!(min.to_bits(), fields[1].parse::<f64>().unwrap().to_bits());
            assert_eq!(max.to_bits(), fields[2].parse::<f64>().unwrap().to_bits());
            assert_eq!(
                weight.to_bits(),
                fields[3].parse::<f64>().unwrap().to_bits()
            );
            assert_eq!(
                *sampling,
                match fields[4] {
                    "log_uniform" => FrequencySampling::LogUniform,
                    "linear_uniform" => FrequencySampling::LinearUniform,
                    _ => panic!("unexpected versioned SFS sampling literal"),
                }
            );
        }
    }
}
