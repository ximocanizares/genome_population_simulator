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
    DeclaredReductionTarget, DeclaredReferenceSource, DeclaredStutter,
    DeclaredTechnicalEventTarget, ExpectedConfigKind, SequenceClass,
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
use crate::pcr_amplification::PcrSettings;
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
use crate::reference::Reference;
use crate::reference_generation::normalize_composition_probabilities;
use crate::repeat_slippage::{SlippageSettings, generate_allele_fan};
use crate::stutter::{StutterError, StutterMotifLengthRates, StutterSettings, StutterStepWeights};
use crate::tandem_repeat::{RepeatMotif, TandemRepeatLocus, TandemRepeatTract};
use crate::technical_profile::LibraryPreparation;
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
    origin: RepeatAlleleOrigin,
}

/// De dónde salió un alelo de repetición.
///
/// Se publica **por alelo**, no por locus, aunque hoy un locus tenga todos sus
/// alelos del mismo origen: un documento con loci de los dos tipos tiene que
/// poder distinguirlos alelo a alelo sin volver a la configuración, y quien
/// encadene un análisis necesita saber si una frecuencia la escribió una
/// persona o la produjo un modelo.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RepeatAlleleOrigin {
    /// Lo escribió el usuario, con su frecuencia.
    Declared,
    /// Lo produjo el generador de deslizamiento biológico partiendo del alelo
    /// de la referencia.
    Generated,
}

impl RepeatAlleleOrigin {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Declared => "declared",
            Self::Generated => "generated",
        }
    }
}

impl EffectiveRepeatAllele {
    #[must_use]
    pub const fn tract(&self) -> &TandemRepeatTract {
        &self.tract
    }

    /// Origen del alelo: declarado por el usuario o generado por deslizamiento.
    #[must_use]
    pub const fn origin(&self) -> RepeatAlleleOrigin {
        self.origin
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
/// Cota **superior** de las bases de secuencia de tracto que una ejecución
/// llegará a publicar, calculada sobre lo declarado y antes de materializar a
/// nadie.
///
/// # De dónde sale cada término, en forma cerrada
///
/// - **Verdad.** Un documento de verdad publica el tracto de una copia en
///   **tres** sitios: la vista vigente de la copia, su historial, y la lista de
///   alelos del caso genotípico. Los dos primeros son uno por copia; el tercero
///   enumera los alelos **distintos** del caso, a lo sumo la ploidía por
///   individuo, que sumada sobre los individuos no pasa del número de copias.
///   De ahí `3 * copias * alelo más largo` por locus.
/// - **VCF.** **No se estima: se mide** (DTC-097, decidido el 2 de octubre
///   de 2026). Las variantes que pisan un registro de repetición se funden en
///   él (DTC-064) y lo ensanchan y le añaden `ALT` de maneras que ninguna
///   fórmula sobre lo declarado acotaba: la que había se quedaba corta. La
///   publicación del VCF cuenta, antes de escribir nada, las bases de tracto
///   que va a escribir, y esa cifra es la que se publica en lugar de la cota;
///   el tope declarado se aplica entonces a la verdad ya **realizada** más
///   esa cifra. Aquí, antes de materializar a nadie, solo se conoce la cota de
///   la verdad, y el tope se le aplica ya **por precaución**: es una cota
///   superior, y una ejecución que la supere podría caber, porque lo
///   realizado suele quedar por debajo. A la inversa, una configuración que
///   `validate` admite puede ser rechazada por `run` al publicar el VCF.
///
/// # Por qué es superior y no exacta, y por qué eso está bien
///
/// Supone que **toda** copia saca el alelo más largo y que **todos** los alelos
/// declarados llegan a sortearse. Lo realizado será casi siempre menor: una
/// copia que saque el alelo de la referencia no publica nada, y un alelo que no
/// le toque a nadie no produce `ALT`. Esa diferencia es el comportamiento
/// correcto, no un defecto que haya que «arreglar» igualando las dos cifras.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveRepeatOutputBudget {
    per_locus: Vec<RepeatLocusOutputBound>,
    declared_cap: Option<u64>,
}

/// Lo que un locus concreto aporta a la cota de la verdad.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RepeatLocusOutputBound {
    locus_id: String,
    truth_bases: u64,
}

impl RepeatLocusOutputBound {
    #[must_use]
    pub fn locus_id(&self) -> &str {
        &self.locus_id
    }

    #[must_use]
    pub const fn truth_bases(&self) -> u64 {
        self.truth_bases
    }
}

impl EffectiveRepeatOutputBudget {
    /// Las tres copias del tracto que un documento de verdad publica por copia
    /// haplotípica: la vista vigente, el historial y el caso genotípico.
    pub const TRUTH_COPIES_PER_HAPLOTYPE_COPY: u64 = 3;

    #[must_use]
    pub fn per_locus(&self) -> &[RepeatLocusOutputBound] {
        &self.per_locus
    }

    #[must_use]
    pub const fn declared_cap(&self) -> Option<u64> {
        self.declared_cap
    }

    #[must_use]
    pub fn truth_bases(&self) -> u64 {
        self.per_locus
            .iter()
            .fold(0, |total, locus| total.saturating_add(locus.truth_bases))
    }

    /// La parte de la verdad, que es todo lo que se conoce antes de
    /// materializar a nadie: el VCF se mide al publicarlo.
    #[must_use]
    pub fn total_bases(&self) -> u64 {
        self.truth_bases()
    }

    /// El locus que más aporta, que es lo que un rechazo tiene que nombrar
    /// para ser accionable. El identificador desempata para que el mensaje no
    /// dependa del orden de recorrido; se compara prestado, sin copiarlo.
    #[must_use]
    pub fn heaviest_locus(&self) -> Option<&RepeatLocusOutputBound> {
        self.per_locus.iter().max_by(|left, right| {
            left.truth_bases
                .cmp(&right.truth_bases)
                .then_with(|| left.locus_id.cmp(&right.locus_id))
        })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveRepeatLocus {
    /// El locus del dominio, con su identificador tipado, su intervalo en el
    /// servicio común de coordenadas y su motivo validado.
    ///
    /// La frontera de configuración **no** mantiene una representación paralela
    /// con campos crudos: había dos del mismo concepto y una de ellas, la
    /// tipada, no tenía usuarios.
    declaration: TandemRepeatLocus,
    /// Unidades enteras de cada bloque del tracto que la referencia simulada
    /// lleva colocado. Solo existe en un locus compuesto sobre referencia
    /// simulada, ya validado contra el tramo.
    reference_units: Option<Vec<u64>>,
    /// Lo que el usuario escribió para situar el locus, en coordenadas del
    /// **genoma de partida**. Solo existe cuando la ejecución inserta algún
    /// tracto (DTC-098): entonces el intervalo de `declaration` está en
    /// coordenadas del genoma **final**, y esto es lo que permite publicar
    /// además lo declarado. Sin inserción no hay traducción ni nada que
    /// publicar aparte, y queda en `None`.
    declared_position: Option<DeclaredLocusPosition>,
    max_units: u64,
    admission: EffectiveRepeatAdmission,
    slippage: Option<EffectiveRepeatSlippage>,
    alleles: Vec<EffectiveRepeatAllele>,
    normalization: RepeatFrequencyNormalization,
    /// De dónde sale el locus: declarado por el usuario o detectado en la
    /// referencia final (REQ-CONF-56, REQ-MAN-98).
    origin: RepeatLocusOrigin,
    /// La puntuación del detector para un locus detectado, que publica el
    /// catálogo aparte (REQ-MAN-98); `None` en un declarado, cuyas fronteras
    /// son las del usuario.
    detector_score: Option<i32>,
}

/// De dónde sale un locus del catálogo de microsatélites.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RepeatLocusOrigin {
    /// Lo declaró el usuario; con la detección, ajusta uno detectado
    /// (REQ-REF-40) o es un motivo de 7 pb o más, que la detección no ve.
    Declared,
    /// Lo encontró la detección sobre la referencia final (REQ-REF-38).
    Detected,
}

impl RepeatLocusOrigin {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Declared => "declared",
            Self::Detected => "detected",
        }
    }
}

/// Versión del algoritmo de detección que se publica con el catálogo
/// (DTC-103): el núcleo copiado con su única desviación y los criterios del
/// apartado 11.7 del plan.
pub const REPEAT_DETECTION_ALGORITHM_VERSION: &str = "ssr-detection-v1";

/// La sección de microsatélites ya resuelta (REQ-CONF-142).
#[derive(Clone, Debug, PartialEq)]
pub struct EffectiveRepeatDetection {
    min_copies: [u64; crate::repeat_detection::DETECTION_MOTIF_LENGTHS],
    max_units: u64,
    slippage: SlippageSettings,
    /// En el orden del catálogo de la referencia final y por inicio; `None`
    /// es la referencia entera.
    regions: Option<Vec<EffectiveDetectionRegion>>,
    /// Loci detectados que un SNV suministrado o escrito a mano retiró del
    /// catálogo (REQ-CONF-143), en el orden del catálogo.
    retired: Vec<RetiredRepeatLocus>,
    /// Lo que dejó el barrido, para el catálogo aparte (REQ-MAN-98); vacío
    /// hasta que se detecta (`with_detected_repeat_loci`).
    scan: DetectionScan,
}

/// Lo que el barrido encontró y adónde fue a parar cada tramo aceptado, para
/// publicarlo en el metadato del catálogo (REQ-MAN-98).
#[derive(Clone, Debug, Default, PartialEq)]
pub struct DetectionScan {
    /// Intervalos del detector, descartes por etapa y aceptados, por longitud
    /// de motivo, sumados sobre todas las secuencias barridas.
    counts: crate::repeat_detection::DetectionCounts,
    /// Tramos aceptados que un locus declarado sustituye (REQ-REF-40), en el
    /// orden del catálogo.
    replaced: Vec<CatalogTract>,
    /// Aceptados que no entran por falta de copias haplotípicas (decisión del
    /// 30/10/2026), por longitud de motivo.
    without_haplotype_copies: [u64; crate::repeat_detection::DETECTION_MOTIF_LENGTHS],
}

impl DetectionScan {
    #[must_use]
    pub const fn counts(&self) -> &crate::repeat_detection::DetectionCounts {
        &self.counts
    }

    #[must_use]
    pub fn replaced(&self) -> &[CatalogTract] {
        &self.replaced
    }

    #[must_use]
    pub const fn without_haplotype_copies(
        &self,
    ) -> &[u64; crate::repeat_detection::DETECTION_MOTIF_LENGTHS] {
        &self.without_haplotype_copies
    }
}

/// Un tramo detectado tal y como lo publica el catálogo: base uno, cerrado,
/// en la referencia final.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CatalogTract {
    sequence_id: String,
    start: u64,
    end: u64,
    motif: String,
    score: i32,
}

impl CatalogTract {
    #[must_use]
    pub fn sequence_id(&self) -> &str {
        &self.sequence_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    #[must_use]
    pub fn motif(&self) -> &str {
        &self.motif
    }

    #[must_use]
    pub const fn score(&self) -> i32 {
        self.score
    }
}

impl EffectiveRepeatDetection {
    /// Mínimo de copias completas por longitud de motivo, de 1 a 6 pb.
    #[must_use]
    pub const fn min_copies(&self) -> &[u64; crate::repeat_detection::DETECTION_MOTIF_LENGTHS] {
        &self.min_copies
    }

    /// Tope de unidades de un locus detectado.
    #[must_use]
    pub const fn max_units(&self) -> u64 {
        self.max_units
    }

    /// El deslizamiento común.
    #[must_use]
    pub const fn slippage(&self) -> SlippageSettings {
        self.slippage
    }

    #[must_use]
    pub fn regions(&self) -> Option<&[EffectiveDetectionRegion]> {
        self.regions.as_deref()
    }

    #[must_use]
    pub fn retired(&self) -> &[RetiredRepeatLocus] {
        &self.retired
    }

    /// Lo que dejó el barrido.
    #[must_use]
    pub const fn scan(&self) -> &DetectionScan {
        &self.scan
    }

    /// Los criterios del detector: los de la sección más los mínimos
    /// **generales** de admisión, sin una segunda perilla (DP-SSR-02).
    #[must_use]
    pub fn criteria(
        &self,
        admission: EffectiveRepeatAdmission,
    ) -> crate::repeat_detection::DetectionCriteria {
        crate::repeat_detection::DetectionCriteria {
            min_copies: self.min_copies,
            min_purity: admission.min_purity(),
            min_longest_run_units: admission.min_longest_run_units(),
            max_units: self.max_units,
        }
    }

    /// Las regiones de una secuencia en base cero y semiabiertas, ordenadas;
    /// `None` si no hay lista (toda la secuencia).
    fn regions_of(&self, sequence_id: &str) -> Option<Vec<(u64, u64)>> {
        self.regions.as_ref().map(|regions| {
            regions
                .iter()
                .filter(|region| region.sequence_id == sequence_id)
                .map(|region| (region.start - 1, region.end))
                .collect()
        })
    }
}

/// Una región de detección: base uno, cerrada, en la referencia final.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveDetectionRegion {
    sequence_id: String,
    start: u64,
    end: u64,
}

impl EffectiveDetectionRegion {
    #[must_use]
    pub fn sequence_id(&self) -> &str {
        &self.sequence_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }
}

/// Un locus detectado que una variante pequeña retiró del catálogo
/// (REQ-CONF-143): queda sin variación y se publica con su causa.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RetiredRepeatLocus {
    locus_id: String,
    sequence_id: String,
    /// Base uno, cerrado, en la referencia final.
    start: u64,
    end: u64,
    motif: String,
    cause: String,
    score: i32,
    variant: RetiringVariant,
}

/// La variante que retiró un locus detectado (REQ-CONF-143): de dónde sale,
/// su identificador y su posición en base uno en la referencia final.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RetiringVariant {
    source: RetiringVariantSource,
    id: String,
    position: u64,
}

impl RetiringVariant {
    #[must_use]
    pub const fn source(&self) -> RetiringVariantSource {
        self.source
    }

    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Base uno, en la referencia final y en la secuencia del locus.
    #[must_use]
    pub const fn position(&self) -> u64 {
        self.position
    }
}

/// De dónde sale una variante que retira un locus.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum RetiringVariantSource {
    /// Un cambio de base escrito a mano (`direct_events`).
    DirectEvent,
    /// Una fila de la matriz de SNV suministrada o reutilizada.
    PopulationSnvs,
}

impl RetiringVariantSource {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::DirectEvent => "direct_event",
            Self::PopulationSnvs => "population_snvs",
        }
    }
}

impl RetiredRepeatLocus {
    #[must_use]
    pub fn locus_id(&self) -> &str {
        &self.locus_id
    }

    #[must_use]
    pub fn sequence_id(&self) -> &str {
        &self.sequence_id
    }

    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    #[must_use]
    pub fn motif(&self) -> &str {
        &self.motif
    }

    /// Qué variante cae dentro, dicho para leerlo.
    #[must_use]
    pub fn cause(&self) -> &str {
        &self.cause
    }

    /// La puntuación del detector del tramo retirado.
    #[must_use]
    pub const fn score(&self) -> i32 {
        self.score
    }

    /// La variante que lo retiró.
    #[must_use]
    pub const fn variant(&self) -> &RetiringVariant {
        &self.variant
    }
}

/// Cómo situó el usuario un locus, en coordenadas del genoma de partida y en
/// base uno.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum DeclaredLocusPosition {
    /// Tramo de un tracto que ya está en el FASTA suministrado.
    Interval { start: u64, end: u64 },
    /// El tracto se inserta justo después de esta base.
    InsertAfter(u64),
}

/// Un punto de inserción resuelto: dónde entra el tracto de un locus y qué
/// bases lleva (DTC-098).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveInsertionPoint {
    locus_id: String,
    sequence_id: String,
    /// Base uno del genoma de partida tras la que entra el tracto.
    insert_after: u64,
    bases: Arc<[u8]>,
}

impl EffectiveInsertionPoint {
    #[must_use]
    pub fn locus_id(&self) -> &str {
        &self.locus_id
    }

    #[must_use]
    pub fn sequence_id(&self) -> &str {
        &self.sequence_id
    }

    /// Base uno del genoma de partida tras la que entra el tracto.
    #[must_use]
    pub const fn insert_after(&self) -> u64 {
        self.insert_after
    }

    /// Bases del tracto: cada motivo repetido tantas veces como digan sus
    /// unidades, en orden.
    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    fn length(&self) -> u64 {
        self.bases.len() as u64
    }
}

/// **El plan de inserción de tractos de una ejecución** (REQ-REF-37,
/// DTC-098).
///
/// Es lo único que hace falta para traducir cualquier posición del genoma de
/// partida al final: los puntos y las longitudes. No necesita las bases del
/// FASTA, y por eso se resuelve con la configuración, antes de leer nada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct EffectiveTractInsertion {
    derived_reference_id: String,
    /// Catálogo del genoma final: el de partida con cada secuencia alargada
    /// por los tractos que recibe.
    sequences: Vec<EffectiveSequence>,
    /// En el orden del catálogo y, dentro de cada secuencia, del punto.
    points: Vec<EffectiveInsertionPoint>,
}

impl EffectiveTractInsertion {
    #[must_use]
    pub fn derived_reference_id(&self) -> &str {
        &self.derived_reference_id
    }

    #[must_use]
    pub fn sequences(&self) -> &[EffectiveSequence] {
        &self.sequences
    }

    #[must_use]
    pub fn points(&self) -> &[EffectiveInsertionPoint] {
        &self.points
    }

    /// Bases insertadas **antes** de la base `position` (base uno, partida)
    /// en su secuencia: las de todo tracto que entra tras una base anterior.
    fn shift_before(&self, sequence_id: &str, position: u64) -> u64 {
        self.points
            .iter()
            .filter(|point| point.sequence_id == sequence_id && point.insert_after < position)
            .map(EffectiveInsertionPoint::length)
            .sum()
    }

    /// Traduce una base del genoma de partida (base uno) a la del final.
    #[must_use]
    pub fn translate_position(&self, sequence_id: &str, position: u64) -> u64 {
        position + self.shift_before(sequence_id, position)
    }

    /// Traduce un tramo cerrado en base uno. Devuelve el punto que lo
    /// **parte** si un tracto se insertaría dentro: entre dos bases del tramo.
    fn translate_interval(
        &self,
        sequence_id: &str,
        start: u64,
        end: u64,
    ) -> Result<(u64, u64), &EffectiveInsertionPoint> {
        if let Some(inside) = self.points.iter().find(|point| {
            point.sequence_id == sequence_id
                && start <= point.insert_after
                && point.insert_after < end
        }) {
            return Err(inside);
        }
        Ok((
            self.translate_position(sequence_id, start),
            self.translate_position(sequence_id, end),
        ))
    }

    /// Primera base (base uno, genoma final) del tracto insertado tras
    /// `insert_after`: la que sigue a la base ancla ya traducida.
    fn inserted_start(&self, sequence_id: &str, insert_after: u64) -> u64 {
        self.translate_position(sequence_id, insert_after) + 1
    }
}

/// El generador de alelos por deslizamiento biológico de un locus.
///
/// Guarda **los parámetros** que el usuario declaró y, una vez generado el
/// abanico, **la contabilidad de esa generación**. Las dos cosas se publican:
/// los parámetros dicen qué se pidió y la contabilidad qué salió, incluidos los
/// eventos que **no se aplicaron** al agotar el tracto, que es la cifra que
/// distingue un servicio que los ignora de uno que los cuenta.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EffectiveRepeatSlippage {
    settings: SlippageSettings,
    generation: Option<EffectiveSlippageGeneration>,
}

impl EffectiveRepeatSlippage {
    #[must_use]
    pub const fn settings(&self) -> SlippageSettings {
        self.settings
    }

    /// La contabilidad de la generación, o `None` mientras no se haya generado.
    #[must_use]
    pub const fn generation(&self) -> Option<EffectiveSlippageGeneration> {
        self.generation
    }
}

/// Lo que la generación de un abanico consumió y produjo.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct EffectiveSlippageGeneration {
    lineages: u64,
    drawn_events: u64,
    applied_events: u64,
    unapplied_events: u64,
}

impl EffectiveSlippageGeneration {
    /// Linajes simulados: una copia haplotípica, un linaje.
    #[must_use]
    pub const fn lineages(self) -> u64 {
        self.lineages
    }

    /// Pasos **sorteados** por la Poisson, sumados sobre todos los linajes.
    #[must_use]
    pub const fn drawn_events(self) -> u64 {
        self.drawn_events
    }

    /// Pasos que llegaron a cambiar el tracto.
    #[must_use]
    pub const fn applied_events(self) -> u64 {
        self.applied_events
    }

    /// Pasos que **no se aplicaron** al agotar el tracto o al topar.
    #[must_use]
    pub const fn unapplied_events(self) -> u64 {
        self.unapplied_events
    }
}

impl EffectiveRepeatLocus {
    #[must_use]
    pub const fn declaration(&self) -> &TandemRepeatLocus {
        &self.declaration
    }

    /// Las unidades por bloque del tracto que la referencia simulada lleva
    /// colocado, en el orden de los motivos. Solo en un locus compuesto sobre
    /// referencia simulada; la suma ya cuadra con el tramo.
    #[must_use]
    pub fn reference_units(&self) -> Option<&[u64]> {
        self.reference_units.as_deref()
    }

    #[must_use]
    pub fn id(&self) -> &str {
        self.declaration.id().as_str()
    }

    /// Lo declarado en coordenadas del genoma de partida, solo cuando la
    /// ejecución inserta algún tracto (DTC-098).
    #[must_use]
    pub const fn declared_position(&self) -> Option<DeclaredLocusPosition> {
        self.declared_position
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

    /// El generador de este locus, si lo declara.
    #[must_use]
    pub const fn slippage(&self) -> Option<&EffectiveRepeatSlippage> {
        self.slippage.as_ref()
    }

    #[must_use]
    pub const fn normalization(&self) -> RepeatFrequencyNormalization {
        self.normalization
    }

    /// Declarado por el usuario o detectado.
    #[must_use]
    pub const fn origin(&self) -> RepeatLocusOrigin {
        self.origin
    }

    /// La puntuación del detector, solo en un locus detectado.
    #[must_use]
    pub const fn detector_score(&self) -> Option<i32> {
        self.detector_score
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
    repeat_output_budget: EffectiveRepeatOutputBudget,
    /// Plan de inserción de tractos, si algún locus inserta el suyo en la
    /// referencia suministrada (DTC-098). Con él, la referencia de la
    /// simulación es la final y no la suministrada.
    tract_insertion: Option<EffectiveTractInsertion>,
    /// La sección de microsatélites, si el documento la declara
    /// (REQ-CONF-142). Sin ella no hay ningún locus.
    repeat_detection: Option<EffectiveRepeatDetection>,
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
    stutter: StutterSettings,
    /// Preparación efectiva de la biblioteca. Por omisión, sin amplificar.
    preparation: LibraryPreparation,
    /// Parámetros efectivos de la amplificación. `None` con preparación sin
    /// PCR, y entonces ni el documento efectivo ni el manifiesto la mencionan.
    pcr: Option<PcrSettings>,
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

    /// Identificador del perfil incorporado que corresponde a la preparación.
    #[must_use]
    pub const fn profile(&self) -> &'static str {
        self.preparation.profile_id()
    }

    /// Preparación efectiva de la biblioteca.
    #[must_use]
    pub const fn preparation(&self) -> LibraryPreparation {
        self.preparation
    }

    /// Parámetros efectivos de la amplificación, si la hay.
    #[must_use]
    pub const fn pcr(&self) -> Option<&PcrSettings> {
        self.pcr.as_ref()
    }

    /// Plan tipado que se entrega al muestreador PE. Apagado sin amplificación.
    #[must_use]
    pub fn pcr_plan(&self) -> crate::pcr_amplification::PcrPlan {
        self.pcr.clone().map_or_else(
            crate::pcr_amplification::PcrPlan::default,
            crate::pcr_amplification::PcrPlan::new,
        )
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

    /// Parámetros efectivos del tartamudeo de base. Por omisión, apagado.
    #[must_use]
    pub const fn stutter(&self) -> &StutterSettings {
        &self.stutter
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
        self.pcr.is_some()
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

    /// Plan de inserción de tractos, si la ejecución inserta alguno.
    #[must_use]
    pub const fn tract_insertion(&self) -> Option<&EffectiveTractInsertion> {
        self.tract_insertion.as_ref()
    }

    /// La sección de microsatélites, si el documento la declara.
    #[must_use]
    pub const fn repeat_detection(&self) -> Option<&EffectiveRepeatDetection> {
        self.repeat_detection.as_ref()
    }

    /// **Identificador de la referencia sobre la que se simula.** Es el de la
    /// referencia final cuando hay inserción de tractos (DTC-098) y el de la
    /// declarada en los demás casos. Todo consumidor que siembre un sorteo o
    /// nombre la referencia de la verdad, el VCF o las lecturas usa este, y
    /// no [`Self::reference_id`], que es el del fichero suministrado.
    #[must_use]
    pub fn simulation_reference_id(&self) -> &str {
        self.tract_insertion
            .as_ref()
            .map_or(&self.reference_id, |insertion| {
                &insertion.derived_reference_id
            })
    }

    /// **Detecta los microsatélites de la referencia final** y devuelve la
    /// configuración con el catálogo entero: los detectados más los
    /// declarados, ordenados por secuencia y por inicio (REQ-REF-38 a
    /// REQ-REF-40, DTC-103).
    ///
    /// Devuelve `None` sin la sección de microsatélites: entonces no hay
    /// ningún locus y nada cambia.
    ///
    /// # Qué hace con cada declarado
    ///
    /// Tiene que **coincidir** con un detectado: su intervalo se solapa con un
    /// locus detectado de su secuencia cuyo motivo es una rotación de uno de
    /// los suyos. Entonces mandan sus fronteras y **sustituye** a todos los
    /// detectados que solapa, conservando su identificador y sus alelos. Un
    /// declarado cuyos motivos tienen todos 7 pb o más queda exento, porque
    /// la detección no los ve (REQ-CONF-60), y también sustituye a lo que
    /// solapa. Sin coincidencia se rechaza nombrando el locus, el detectado
    /// más cercano y la etapa en la que cayó el tramo de su motivo.
    ///
    /// # Sin individuos
    ///
    /// Sin copias haplotípicas no se generan alelos de los detectados
    /// (decisión del responsable del 30/10/2026; REQ-CONF-140), y un locus
    /// detectado no tiene otros, mientras que un locus publicado necesita al
    /// menos uno: **los detectados no entran en el catálogo** y la ejecución
    /// no se rechaza. Los declarados siguen sus reglas.
    ///
    /// # Por qué no puede hacerse al resolver la configuración
    ///
    /// Por la misma razón que el generador de deslizamiento: la referencia
    /// final no existe hasta que se prepara.
    pub fn with_detected_repeat_loci(
        &self,
        reference: &Reference,
    ) -> Result<Option<Self>, Diagnostic> {
        let Some(detection) = &self.repeat_detection else {
            return Ok(None);
        };
        let criteria = detection.criteria(self.repeat_admission);
        let reference_id = ReferenceId::parse(self.simulation_reference_id(), "execution")
            .map_err(|diagnostic| diagnostic.with("path", &self.source_path))?;
        let mut catalog: Vec<(usize, u64, EffectiveRepeatLocus)> =
            Vec::with_capacity(self.repeat_loci.len());
        let declared_ids = self
            .repeat_loci
            .iter()
            .map(EffectiveRepeatLocus::id)
            .collect::<BTreeSet<_>>();
        let lineages = self.haplotype_copy_count()?;
        let mut scan = DetectionScan::default();
        for (sequence_index, sequence) in reference.sequences().iter().enumerate() {
            let mut declared: Vec<&EffectiveRepeatLocus> = self
                .repeat_loci
                .iter()
                .filter(|locus| locus.declaration.interval().sequence_id() == sequence.id())
                .collect();
            declared.sort_by_key(|locus| locus.declaration.interval().start());
            let regions = detection.regions_of(sequence.id().as_str());
            if regions.as_ref().is_some_and(Vec::is_empty) && declared.is_empty() {
                continue;
            }
            let bases = sequence.sequence().as_bytes().map_err(|error| {
                Diagnostic::internal_invariant("application", "repeat_detection").with(
                    "detail",
                    format!("las bases de la referencia no están materializadas: {error:?}"),
                )
            })?;
            let found = crate::repeat_detection::detect_repeat_candidates_in_regions(
                bases,
                &criteria,
                regions.as_deref(),
            );
            scan.counts.add(found.counts());
            let candidates = found.candidates();
            // Los candidatos no se solapan y van por inicio: los que pisan
            // `[start, end)` empiezan en el primero que acaba después de
            // `start`.
            let overlapping = |start: u64, end: u64| {
                let first = candidates.partition_point(|candidate| candidate.end() <= start);
                candidates[first..]
                    .iter()
                    .take_while(move |candidate| candidate.start() < end)
            };
            for locus in &declared {
                let interval = locus.declaration.interval();
                let (start, end) = (interval.start(), interval.end());
                let exempt = locus.motifs().iter().all(|motif| {
                    motif.length() > crate::repeat_detection::DETECTION_MOTIF_LENGTHS as u64
                });
                let matches = overlapping(start, end).any(|candidate| {
                    locus.motifs().iter().any(|motif| {
                        crate::repeat_detection::is_rotation(
                            candidate.motif().bases(),
                            motif.bases(),
                        )
                    })
                });
                if !exempt && !matches {
                    return Err(self.unmatched_declared_locus(
                        locus,
                        sequence.id().as_str(),
                        bases,
                        &criteria,
                        regions.as_deref(),
                        candidates,
                    ));
                }
                catalog.push((sequence_index, start, (*locus).clone()));
            }
            for candidate in candidates {
                let first = declared.partition_point(|locus| {
                    locus.declaration.interval().end() <= candidate.start()
                });
                if declared
                    .get(first)
                    .is_some_and(|locus| locus.declaration.interval().start() < candidate.end())
                {
                    scan.replaced.push(CatalogTract {
                        sequence_id: sequence.id().as_str().to_owned(),
                        start: candidate.start() + 1,
                        end: candidate.end(),
                        motif: candidate.motif().as_text().to_owned(),
                        score: candidate.score(),
                    });
                    continue;
                }
                // **Sin copias haplotípicas no hay alelos que generar**, y un
                // locus detectado no tiene otros: no entra en el catálogo
                // (decisión del responsable del 30/10/2026), y se cuenta. Un
                // declarado sí entra, con sus reglas.
                if lineages == 0 {
                    scan.without_haplotype_copies[candidate.motif().length() as usize - 1] += 1;
                    continue;
                }
                let motif = candidate.motif().clone();
                let id = format!(
                    "SSR-{}-{}-{}-{}",
                    sequence_index + 1,
                    candidate.start() + 1,
                    candidate.end(),
                    motif.as_text()
                );
                if declared_ids.contains(id.as_str()) {
                    return Err(Diagnostic::configuration(
                        "GPS-CONF-008",
                        "La configuración contiene una combinación de campos incompatible.",
                    )
                    .with("path", &self.source_path)
                    .with("fields", "reference.repeat_loci.id")
                    .with(
                        "detail",
                        format!(
                            "identificador de locus repetido: {id} es también el de un locus \
                             detectado en otro tramo; los detectados se llaman \
                             SSR-<secuencia>-<inicio>-<fin>-<motivo>, y un declarado no puede \
                             usar ese nombre para otro tramo"
                        ),
                    ));
                }
                let locus_id = RepeatLocusId::parse(&id, "repeat_locus")
                    .map_err(|diagnostic| diagnostic.with("path", &self.source_path))?;
                let sequence_id = SequenceId::parse(sequence.id().as_str(), reference_id.as_str())
                    .map_err(|diagnostic| diagnostic.with("path", &self.source_path))?;
                let interval = GenomicInterval::new(
                    reference_id.clone(),
                    sequence_id,
                    candidate.start(),
                    candidate.end(),
                    sequence.length(),
                )
                .map_err(|error| {
                    Diagnostic::internal_invariant("application", "repeat_detection")
                        .with("detail", format!("{id}: {error:?}"))
                })?;
                let declaration =
                    TandemRepeatLocus::new(locus_id, interval, vec![motif]).map_err(|error| {
                        Diagnostic::internal_invariant("application", "repeat_detection")
                            .with("detail", format!("{id}: {}", error.detail()))
                    })?;
                catalog.push((
                    sequence_index,
                    candidate.start(),
                    EffectiveRepeatLocus {
                        declaration,
                        reference_units: None,
                        declared_position: None,
                        max_units: detection.max_units,
                        admission: self.repeat_admission,
                        slippage: Some(EffectiveRepeatSlippage {
                            settings: detection.slippage,
                            generation: None,
                        }),
                        alleles: Vec::new(),
                        normalization: RepeatFrequencyNormalization {
                            tolerance: self.normalization.tolerance(),
                            original_sum: 0.0,
                            original_delta: 0.0,
                            effective_sum: 0.0,
                            effective_delta: 0.0,
                        },
                        origin: RepeatLocusOrigin::Detected,
                        detector_score: Some(candidate.score()),
                    },
                ));
            }
        }
        catalog.sort_by_key(|(sequence_index, start, _)| (*sequence_index, *start));
        let mut resolved = self.clone();
        resolved.repeat_loci = catalog.into_iter().map(|(_, _, locus)| locus).collect();
        if let Some(detection) = resolved.repeat_detection.as_mut() {
            detection.scan = scan;
        }
        resolved.admit_detected_loci_against_stutter()?;
        Ok(Some(resolved))
    }

    /// **El último paso sobre el catálogo final**, después de detectar, de
    /// retirar lo que invade una variante y de generar los abanicos (en
    /// `validate`, en `run` y en la reconstrucción del juez físico). Solo hace
    /// algo con la sección de microsatélites:
    ///
    /// - con algún locus, la PCR exige el tartamudeo de base
    ///   (`admit_amplified_preparation_against_the_catalog`);
    /// - sin ninguno, **desconecta el tartamudeo clonal**, que se conectó al
    ///   resolver la configuración porque entonces aún no se sabía qué quedaría:
    ///   una biblioteca amplificada sin loci publica lo mismo que sin la sección
    ///   (DTC-092, catálogo de azar v5).
    ///
    /// Devuelve `None` si no cambia nada.
    pub fn with_final_repeat_catalog(&self) -> Result<Option<Self>, Diagnostic> {
        if self.repeat_detection.is_none() {
            return Ok(None);
        }
        if !self.repeat_loci.is_empty() {
            self.admit_amplified_preparation_against_the_catalog()?;
            return Ok(None);
        }
        Ok(self.without_clonal_stutter())
    }

    /// La misma configuración con el tartamudeo clonal desconectado, o `None`
    /// si no lo tenía.
    #[must_use]
    pub fn without_clonal_stutter(&self) -> Option<Self> {
        let event = self.technical_event.as_ref()?;
        event.pcr.as_ref()?.clonal_stutter()?;
        let mut resolved = self.clone();
        if let Some(event) = resolved.technical_event.as_mut() {
            event.pcr = event.pcr.take().map(PcrSettings::without_clonal_stutter);
        }
        Some(resolved)
    }

    /// **La PCR exige el tartamudeo de base si el catálogo final tiene algún
    /// locus** (REQ-CONF-138, con «loci declarados» como los del catálogo,
    /// según REQ-CONF-142). Es la misma regla que
    /// `admit_amplified_preparation_against_loci` aplica al resolver la
    /// configuración con loci declarados, ahora sobre lo que queda después de
    /// detectar y de retirar: sin ningún locus, la PCR se admite sin
    /// tartamudeo (decisión del responsable del 30/10/2026).
    fn admit_amplified_preparation_against_the_catalog(&self) -> Result<(), Diagnostic> {
        let Some(event) = self.technical_event.as_ref() else {
            return Ok(());
        };
        let Some(first) = self.repeat_loci.first() else {
            return Ok(());
        };
        if event.preparation() != LibraryPreparation::Pcr || event.stutter().enabled() {
            return Ok(());
        }
        Err(Diagnostic::configuration(
            "GPS-CONF-008",
            "La configuración contiene una combinación de campos incompatible.",
        )
        .with("path", &self.source_path)
        .with(
            "fields",
            "technical_event.library_preparation,technical_event.stutter.enabled",
        )
        .with(
            "detail",
            format!(
                "la preparación amplificada por PCR con microsatélites en el catálogo ({} loci, el \
                 primero es {}) exige technical_event.stutter.enabled = true: el tartamudeo clonal \
                 de la PCR solo despeja el incremento de la amplificación, y sin la parte de base \
                 la biblioteca subestimaría el tartamudeo. O se activa el tartamudeo de base, o \
                 se limita la detección con regiones que no contengan ninguno, o se usa \
                 library_preparation = \"{}\"",
                self.repeat_loci.len(),
                first.id(),
                LibraryPreparation::PcrFree.as_str()
            ),
        ))
    }

    /// El tartamudeo encendido alcanza a los detectados con las tasas por
    /// longitud de motivo (REQ-CONF-132): una longitud sin tasa declarada se
    /// rechaza como con un declarado, nombrando el locus.
    fn admit_detected_loci_against_stutter(&self) -> Result<(), Diagnostic> {
        let Some(event) = self.technical_event.as_ref() else {
            return Ok(());
        };
        if !event.stutter().enabled() {
            return Ok(());
        }
        for locus in &self.repeat_loci {
            if locus.origin != RepeatLocusOrigin::Detected {
                continue;
            }
            for motif in locus.motifs() {
                let length = motif.length();
                if event
                    .stutter()
                    .motif_length_rates()
                    .rate_for_motif_length(length)
                    .is_none()
                {
                    return Err(Diagnostic::configuration(
                        "GPS-CONF-008",
                        "La configuración contiene una combinación de campos incompatible.",
                    )
                    .with("path", &self.source_path)
                    .with(
                        "fields",
                        "reference.repeat_detection,technical_event.stutter.motif_length_rates",
                    )
                    .with(
                        "detail",
                        format!(
                            "el locus detectado {} tiene un motivo de {length} pb y el tartamudeo \
                             no declara tasa para esa longitud; declárela en \
                             technical_event.stutter.motif_length_rates o limite la detección \
                             con regiones",
                            locus.id()
                        ),
                    ));
                }
            }
        }
        Ok(())
    }

    /// El rechazo de un declarado que no coincide con ningún detectado
    /// (REQ-REF-40): nombra el locus, lo detectado más cercano de su
    /// secuencia y la etapa en la que cayó el tramo de su motivo que lo
    /// solapaba.
    #[allow(clippy::too_many_arguments)]
    fn unmatched_declared_locus(
        &self,
        locus: &EffectiveRepeatLocus,
        sequence_id: &str,
        bases: &[u8],
        criteria: &crate::repeat_detection::DetectionCriteria,
        regions: Option<&[(u64, u64)]>,
        candidates: &[crate::repeat_detection::RepeatCandidate],
    ) -> Diagnostic {
        let interval = locus.declaration.interval();
        let (start, end) = (interval.start(), interval.end());
        let why = crate::repeat_detection::explain_missing_tract(
            bases,
            criteria,
            regions,
            start,
            end,
            locus.motifs(),
        )
        .map_or_else(
            || "el detector no ve ningún tramo de ese motivo que lo solape".to_owned(),
            |missing| {
                format!(
                    "el tramo {sequence_id}:{}..{} de periodo {} no superó la etapa {}",
                    missing.start + 1,
                    missing.end,
                    missing.period,
                    missing.stage.name()
                )
            },
        );
        let nearest = candidates
            .iter()
            .min_by_key(|candidate| {
                if candidate.end() <= start {
                    start - candidate.end()
                } else {
                    candidate.start().saturating_sub(end)
                }
            })
            .map_or_else(
                || format!("no hay ningún locus detectado en {sequence_id}"),
                |candidate| {
                    format!(
                        "el detectado más cercano es {sequence_id}:{}..{} con motivo {}",
                        candidate.start() + 1,
                        candidate.end(),
                        candidate.motif().as_text()
                    )
                },
            );
        Diagnostic::configuration(
            "GPS-CONF-008",
            "La configuración contiene una combinación de campos incompatible.",
        )
        .with("path", &self.source_path)
        .with("fields", "reference.repeat_loci,reference.repeat_detection")
        .with(
            "detail",
            format!(
                "locus {}: no coincide con ningún microsatélite detectado en la referencia final \
                 ({sequence_id}:{}..{}, motivos {}): {why}; {nearest}. Un locus declarado solo \
                 ajusta los parámetros de uno detectado",
                locus.id(),
                start + 1,
                end,
                locus
                    .motifs()
                    .iter()
                    .map(RepeatMotif::as_text)
                    .collect::<Vec<_>>()
                    .join(", ")
            ),
        )
    }

    /// **Retira del catálogo los loci detectados que contienen una variante
    /// pequeña** suministrada o escrita a mano (REQ-CONF-143, DP-SSR-11,
    /// provisional hasta DTF-066): un cambio de base directo, o una fila de la
    /// matriz de SNV, en coordenadas de la referencia final. El locus queda
    /// sin variación y la sección publica cuál y por qué.
    ///
    /// Los **declarados** no se retiran: una variante dentro de uno se sigue
    /// rechazando (S-32) más adelante, en la admisión de los loci frente a las
    /// variantes pequeñas. Devuelve `None` sin la sección o si no retira
    /// nada.
    ///
    /// `matrix_rows` son `(secuencia, posición en base uno, identificador)`.
    pub fn with_retired_repeat_loci(
        &self,
        matrix_rows: &[(&str, u64, &str)],
    ) -> Result<Option<Self>, Diagnostic> {
        if self.repeat_detection.is_none() {
            return Ok(None);
        }
        let mut variants: std::collections::BTreeMap<(&str, u64), (String, RetiringVariant)> =
            std::collections::BTreeMap::new();
        for event in &self.direct_events {
            if let EffectiveDirectEvent::Snp {
                sequence_id,
                position_1based,
                ..
            } = event
                && *position_1based > 0
            {
                variants
                    .entry((sequence_id.as_str(), position_1based - 1))
                    .or_insert_with(|| {
                        (
                            format!(
                                "el cambio de base {} en {sequence_id}:{position_1based}",
                                event.event_id()
                            ),
                            RetiringVariant {
                                source: RetiringVariantSource::DirectEvent,
                                id: event.event_id().as_str().to_owned(),
                                position: *position_1based,
                            },
                        )
                    });
            }
        }
        for (contig, position_1based, variant_id) in matrix_rows {
            if *position_1based > 0 {
                variants
                    .entry((contig, position_1based - 1))
                    .or_insert_with(|| {
                        (
                            format!(
                                "la variante {variant_id} de la matriz en {contig}:{position_1based}"
                            ),
                            RetiringVariant {
                                source: RetiringVariantSource::PopulationSnvs,
                                id: (*variant_id).to_owned(),
                                position: *position_1based,
                            },
                        )
                    });
            }
        }
        let mut retired = Vec::new();
        let mut kept = Vec::with_capacity(self.repeat_loci.len());
        for locus in &self.repeat_loci {
            let interval = locus.declaration.interval();
            let sequence_id = interval.sequence_id().as_str();
            let inside = (locus.origin == RepeatLocusOrigin::Detected)
                .then(|| {
                    variants
                        .range((sequence_id, interval.start())..(sequence_id, interval.end()))
                        .next()
                })
                .flatten();
            match inside {
                Some((_, (cause, variant))) => retired.push(RetiredRepeatLocus {
                    locus_id: locus.id().to_owned(),
                    sequence_id: sequence_id.to_owned(),
                    start: interval.start() + 1,
                    end: interval.end(),
                    motif: locus
                        .motifs()
                        .iter()
                        .map(RepeatMotif::as_text)
                        .collect::<Vec<_>>()
                        .join(","),
                    cause: format!("{cause} cae dentro del tracto"),
                    score: locus.detector_score.unwrap_or_default(),
                    variant: variant.clone(),
                }),
                None => kept.push(locus.clone()),
            }
        }
        if retired.is_empty() {
            return Ok(None);
        }
        let mut resolved = self.clone();
        resolved.repeat_loci = kept;
        if let Some(detection) = resolved.repeat_detection.as_mut() {
            detection.retired = retired;
        }
        Ok(Some(resolved))
    }

    /// Genera el abanico de alelos de los loci que lo piden, partiendo del
    /// **alelo que el genoma de referencia lleva** en cada uno.
    ///
    /// # Por qué no puede hacerse al resolver la configuración
    ///
    /// El alelo ancestral son bases de la referencia, y la referencia no existe
    /// hasta que se prepara. Por eso el generador se aplica aquí, sobre una
    /// configuración ya resuelta, y devuelve **otra** configuración con los
    /// alelos dentro; la original no se toca.
    ///
    /// Devuelve `None` cuando **ningún locus declara generador**. No es una
    /// optimización: es lo que garantiza que un escenario que no lo usa siga
    /// publicando exactamente los mismos bytes.
    ///
    /// # La cota de bases de tracto se recomputa
    ///
    /// La cota sale de los alelos, y los alelos acaban de nacer. Dejarla como
    /// estaba publicaría una cota falsa —cero bases para un locus que sí va a
    /// escribir tracto— y el juez independiente encontraría lo realizado por
    /// encima de lo acotado. El tope declarado se vuelve a aplicar aquí, que
    /// sigue siendo antes de publicar nada.
    pub fn with_generated_repeat_alleles(
        &self,
        reference: &Reference,
    ) -> Result<Option<Self>, Diagnostic> {
        if self
            .repeat_loci
            .iter()
            .all(|locus| locus.slippage.is_none())
        {
            return Ok(None);
        }
        let lineages = self.haplotype_copy_count()?;
        let mut resolved = self.clone();
        for locus in &mut resolved.repeat_loci {
            let Some(slippage) = locus.slippage else {
                continue;
            };
            let ancestral = reference_repeat_tract(reference, locus)?;
            let fan = generate_allele_fan(
                self.reproducibility.root_seed(),
                self.simulation_reference_id(),
                locus.declaration.id().as_str(),
                &ancestral,
                locus.max_units,
                slippage.settings(),
                lineages,
            )
            .map_err(|error| {
                error.into_diagnostic(&self.source_path, locus.declaration.id().as_str())
            })?;
            locus.slippage = Some(EffectiveRepeatSlippage {
                settings: slippage.settings(),
                generation: Some(EffectiveSlippageGeneration {
                    lineages: fan.lineages(),
                    drawn_events: fan.drawn_events(),
                    applied_events: fan.applied_events(),
                    unapplied_events: fan.unapplied_events(),
                }),
            });
            // Las frecuencias generadas ya suman uno por construcción, pero la
            // auditoría se **recomputa** con la suma compensada de DTC-049 en
            // vez de darlo por hecho: es lo que el resto del proyecto hace y lo
            // que hace publicable la normalización.
            let declared: Vec<f64> = fan
                .alleles()
                .iter()
                .map(crate::repeat_slippage::SlippageAllele::frequency)
                .collect();
            let original_sum = neumaier_sum(declared.iter().copied());
            let effective: Vec<f64> = declared
                .iter()
                .map(|weight| weight / original_sum)
                .collect();
            let effective_sum = neumaier_sum(effective.iter().copied());
            locus.normalization = RepeatFrequencyNormalization {
                tolerance: self.normalization.tolerance(),
                original_sum,
                original_delta: (original_sum - 1.0).abs(),
                effective_sum,
                effective_delta: (effective_sum - 1.0).abs(),
            };
            locus.alleles = fan
                .into_alleles()
                .into_iter()
                .zip(effective)
                .map(|(allele, frequency)| {
                    let (tract, _, declared_frequency) = allele.into_parts();
                    EffectiveRepeatAllele {
                        tract,
                        declared_frequency,
                        frequency,
                        origin: RepeatAlleleOrigin::Generated,
                    }
                })
                .collect();
        }
        resolved.repeat_output_budget = recompute_repeat_output_budget(
            &resolved.repeat_loci,
            lineages,
            resolved.repeat_output_budget.declared_cap(),
        );
        if let Some(cap) = resolved.repeat_output_budget.declared_cap()
            && resolved.repeat_output_budget.total_bases() > cap
        {
            let heaviest = resolved.repeat_output_budget.heaviest_locus().map_or_else(
                || "ninguno".to_owned(),
                |locus| format!("{} con {} bases", locus.locus_id(), locus.truth_bases()),
            );
            return Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("path", &self.source_path)
            .with("fields", "reference.max_tract_output_bases")
            .with(
                "detail",
                format!(
                    "los alelos generados por deslizamiento elevan la cota de bases de tracto de la \
verdad a {} y el tope declarado es {cap}; el locus que más aporta es {heaviest}",
                    resolved.repeat_output_budget.total_bases()
                ),
            ));
        }
        Ok(Some(resolved))
    }

    /// Copias haplotípicas de todas las poblaciones iniciales juntas.
    ///
    /// Es el número de **linajes** que el generador simula: una copia, un
    /// linaje.
    pub fn haplotype_copy_count(&self) -> Result<u64, Diagnostic> {
        self.initial_populations
            .iter()
            .try_fold(0_u64, |total, population| {
                u64::try_from(population.individual_count())
                    .ok()
                    .and_then(|individuals| {
                        individuals.checked_mul(u64::from(population.ploidy().copy_count()))
                    })
                    .and_then(|copies| total.checked_add(copies))
                    .ok_or_else(|| {
                        Diagnostic::internal_invariant("application", "haplotype_copy_count")
                            .with("detail", "el recuento de copias haplotípicas desborda")
                    })
            })
    }

    #[must_use]
    pub const fn repeat_admission(&self) -> EffectiveRepeatAdmission {
        self.repeat_admission
    }

    /// Lo que el evento técnico necesita para tartamudear: los parámetros
    /// efectivos y los loci declarados, en un solo objeto.
    ///
    /// Sin evento técnico o sin loci el plan queda apagado, y un plan apagado
    /// no abre la etapa `stutter` ni consume una sola palabra aleatoria.
    #[must_use]
    pub fn stutter_plan(&self) -> crate::stutter::StutterPlan {
        let settings = self
            .technical_event
            .as_ref()
            .map_or_else(StutterSettings::default, |event| event.stutter().clone());
        let loci = self
            .repeat_loci
            .iter()
            .map(|locus| {
                crate::stutter::StutterLocus::new(
                    locus.declaration().id().clone(),
                    locus.declaration().interval().clone(),
                    locus.declaration().motifs().to_vec(),
                )
            })
            .collect();
        crate::stutter::StutterPlan::new(settings, loci)
    }

    /// Plan de amplificación del evento técnico de esta ejecución.
    ///
    /// Sin evento técnico, o con preparación sin PCR, el plan está **apagado**
    /// y el recorrido no abre la etapa `amplification` ni consume una sola
    /// palabra aleatoria.
    #[must_use]
    pub fn pcr_plan(&self) -> crate::pcr_amplification::PcrPlan {
        self.technical_event
            .as_ref()
            .map_or_else(crate::pcr_amplification::PcrPlan::default, |event| {
                event.pcr_plan()
            })
    }

    /// Cota superior de bases de tracto publicadas, calculada sobre lo
    /// declarado.
    #[must_use]
    pub const fn repeat_output_budget(&self) -> &EffectiveRepeatOutputBudget {
        &self.repeat_output_budget
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
    // **El plan de inserción va antes que los loci**: es lo que traduce los
    // intervalos de los demás al genoma final (DTC-098).
    let tract_insertion = effective_tract_insertion(general, &sequences)?;
    // Las regiones de detección viven en el genoma **final** (DP-SSR-18).
    let repeat_detection = effective_repeat_detection(
        general,
        tract_insertion
            .as_ref()
            .map_or(sequences.as_slice(), EffectiveTractInsertion::sequences),
    )?;
    let repeat_loci = effective_repeat_loci(
        general,
        &sequences,
        normalization.tolerance,
        repeat_admission,
        tract_insertion.as_ref(),
        repeat_detection
            .as_ref()
            .map(EffectiveRepeatDetection::slippage),
    )?;
    let reproducibility = ReproducibilityMetadata::new(root_seed, seed_origin);
    let initial_populations = effective_initial_populations(general, root_seed)?;
    let mut population_snvs = effective_population_snvs(general)?;
    if let Some(insertion) = &tract_insertion {
        translate_generation_regions(general, insertion, population_snvs.as_mut())?;
    }
    // Las regiones se comprueban contra el catálogo del genoma sobre el que se
    // generará la matriz: el final si hay inserción.
    let generation_sequences = tract_insertion
        .as_ref()
        .map_or(sequences.as_slice(), EffectiveTractInsertion::sequences);
    let direct_events = effective_direct_events(general, &initial_populations, &sequences)?;
    let mut technical_event =
        effective_technical_event(general, &initial_populations, &direct_events)?;
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
            let valid = generation_sequences
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
    // La cota se calcula **antes** de ceder las poblaciones al documento
    // efectivo, y su rechazo ocurre aquí: antes de materializar a nadie.
    let repeat_output_budget = repeat_output_budget(general, &repeat_loci, &initial_populations)?;
    admit_loci_against_stutter(general, technical_event.as_ref(), &repeat_loci)?;
    // **La PCR sin tartamudeo de base se rechaza aquí solo con loci
    // declarados**: con ellos el catálogo no puede quedar vacío. Los
    // detectados solo se conocen al barrer la referencia, después de esta
    // frontera, y la misma regla se aplica entonces sobre el catálogo entero
    // (`with_detected_repeat_loci`): sin ningún locus, la PCR se admite sin
    // tartamudeo (decisión del responsable del 30/10/2026). El tartamudeo
    // clonal, en cambio, se conecta aquí con la sección, porque su tabla se
    // resuelve con el documento, y se desconecta después si el catálogo final
    // se queda vacío (`with_final_repeat_catalog`).
    let microsatellites = !repeat_loci.is_empty() || repeat_detection.is_some();
    admit_amplified_preparation_against_loci(general, technical_event.as_ref(), &repeat_loci)?;
    if repeat_detection.is_some()
        && let Some(event) = technical_event.as_ref()
        && event.stutter().enabled()
        && repeat_admission.min_longest_run_units() < 2
    {
        return Err(config_error(
            general,
            "reference.repeat_admission.min_longest_run_units,technical_event.stutter",
            format!(
                "la tirada mínima general es {} unidad(es) y el tartamudeo exige al menos 2 en los \
                 loci detectados, que la usan: una tirada de una sola copia no puede contraerse \
                 nunca",
                repeat_admission.min_longest_run_units()
            ),
        ));
    }
    attach_clonal_stutter(general, technical_event.as_mut(), microsatellites)?;
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
        repeat_output_budget,
        repeat_loci,
        repeat_admission,
        tract_insertion,
        repeat_detection,
    })
}

/// **La sección de microsatélites** (REQ-CONF-142), resuelta y comprobada
/// antes de leer la referencia: los mínimos y el tope, el deslizamiento común
/// y las regiones contra el catálogo de la referencia **final**.
fn effective_repeat_detection(
    general: &DeclaredGeneralConfig,
    frame_sequences: &[EffectiveSequence],
) -> Result<Option<EffectiveRepeatDetection>, Diagnostic> {
    let Some(declared) = &general.reference().repeat_detection else {
        return Ok(None);
    };
    // **Detección y reducción no van en la misma ejecución** (DP-SSR-12): la
    // simulación correría sobre la suministrada entera y no sobre la reducida.
    if !general.reduced_references().is_empty() {
        return Err(config_error(
            general,
            "reference.repeat_detection,reduced_references",
            "la detección de microsatélites no se combina con una reducción en la misma \
             ejecución: la simulación correría sobre la referencia suministrada entera, no sobre \
             la reducida. Hágalo en dos ejecuciones: en la primera, solo la reducción, que \
             publica el FASTA reducido; en la segunda, ese FASTA como reference.fasta y la \
             sección [reference.repeat_detection]",
        ));
    }
    let mut min_copies = crate::repeat_detection::DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH;
    if let Some(declared_copies) = &declared.min_copies {
        let values = declared_copies.value();
        if values.len() != min_copies.len() {
            return Err(config_error(
                general,
                "reference.repeat_detection.min_copies",
                format!(
                    "min_copies lleva {} valores y exige exactamente {}, uno por longitud de \
                     motivo de 1 a 6 pb; por omisión es {:?}",
                    values.len(),
                    min_copies.len(),
                    crate::repeat_detection::DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH
                ),
            ));
        }
        if let Some(position) = values.iter().position(|value| *value == 0) {
            return Err(config_error(
                general,
                "reference.repeat_detection.min_copies",
                format!(
                    "el mínimo de copias de los motivos de {} pb es 0; cada mínimo es un entero \
                     positivo",
                    position + 1
                ),
            ));
        }
        min_copies.copy_from_slice(values);
    }
    let max_units = declared.max_units.as_ref().map_or(
        crate::repeat_detection::DEFAULT_DETECTION_MAX_UNITS,
        |value| *value.value(),
    );
    if max_units == 0 {
        return Err(config_error(
            general,
            "reference.repeat_detection.max_units",
            "el tope de unidades de los loci detectados debe ser un entero positivo",
        ));
    }
    let slippage = &declared.slippage;
    let number = |field: &str, value: &crate::configuration::DeclaredNumber| {
        value.literal().parse::<f64>().map_err(|_| {
            config_error(
                general,
                &format!("reference.repeat_detection.slippage.{field}"),
                format!("{field} no es un decimal válido: {}", value.literal()),
            )
        })
    };
    let settings = SlippageSettings::new(
        slippage.model.value().as_str(),
        slippage
            .step_units
            .as_ref()
            .map_or(crate::repeat_slippage::SLIPPAGE_STEP_UNITS, |value| {
                *value.value()
            }),
        number("mean_events", &slippage.mean_events)?,
        match &slippage.expansion_bias {
            None => crate::repeat_slippage::SLIPPAGE_DEFAULT_EXPANSION_BIAS,
            Some(value) => number("expansion_bias", value)?,
        },
    )
    .map_err(|error| {
        config_error(
            general,
            &error.parameter().replace(
                "reference.repeat_loci.slippage",
                "reference.repeat_detection.slippage",
            ),
            format!("deslizamiento común: {}", error.detail()),
        )
    })?;
    let regions = match &declared.regions {
        None => None,
        Some(declared_regions) => {
            const FIELD: &str = "reference.repeat_detection.regions";
            if declared_regions.is_empty() {
                return Err(config_error(
                    general,
                    FIELD,
                    "la lista de regiones está vacía: no se detectaría nada. Declare al menos una \
                     región o retire la lista para barrer la referencia entera",
                ));
            }
            let mut resolved = Vec::with_capacity(declared_regions.len());
            for region in declared_regions {
                let (start, end) = (*region.start.value(), *region.end.value());
                let (index, sequence) = frame_sequences
                    .iter()
                    .enumerate()
                    .find(|(_, sequence)| sequence.id() == region.sequence.value())
                    .ok_or_else(|| {
                        config_error(
                            general,
                            FIELD,
                            format!(
                                "la región nombra la secuencia {}, que no está en la referencia \
                                 final",
                                region.sequence.value()
                            ),
                        )
                    })?;
                if start == 0 || start > end || end > sequence.length() {
                    return Err(config_error(
                        general,
                        FIELD,
                        format!(
                            "la región {}:{start}..{end} no cabe en la secuencia ({} bases de la \
                             referencia final; base uno, cerrada)",
                            sequence.id(),
                            sequence.length()
                        ),
                    ));
                }
                resolved.push((index, start, end, sequence.id().to_owned()));
            }
            resolved.sort();
            for pair in resolved.windows(2) {
                if pair[0].0 == pair[1].0 && pair[1].1 <= pair[0].2 {
                    return Err(config_error(
                        general,
                        FIELD,
                        format!(
                            "las regiones {}:{}..{} y {}:{}..{} se solapan; declárelas como una sola",
                            pair[0].3, pair[0].1, pair[0].2, pair[1].3, pair[1].1, pair[1].2
                        ),
                    ));
                }
            }
            Some(
                resolved
                    .into_iter()
                    .map(|(_, start, end, sequence_id)| EffectiveDetectionRegion {
                        sequence_id,
                        start,
                        end,
                    })
                    .collect(),
            )
        }
    };
    Ok(Some(EffectiveRepeatDetection {
        min_copies,
        max_units,
        slippage: settings,
        regions,
        retired: Vec::new(),
        scan: DetectionScan::default(),
    }))
}

/// **El plan de inserción de tractos** (REQ-REF-37, DTC-098).
///
/// Reúne los loci que declaran `insert_after`, construye el tracto de cada uno
/// desde su motivo y sus `reference_units`, y rechaza antes de simular todo lo
/// que esta primera versión no admite junto a una inserción. Sin ningún locus
/// que inserte devuelve `None`, y entonces `inserted_reference_id` no puede
/// estar declarado.
///
/// # Lo que se rechaza, y por qué
///
/// - Dos inserciones en el mismo punto: fijar su orden sería una decisión que
///   nadie ha tomado.
/// - Un punto antes de la primera base o después de la última: toda inserción
///   configurada se ancla en una base anterior y deja un punto interno
///   (REQ-IND-40).
/// - Los cambios directos: sus inserciones y deleciones se declaran en
///   coordenadas **vigentes** de la copia, que no tienen traducción sin un
///   segundo modelo de coordenadas.
/// - La matriz de SNV suministrada o reutilizada: está atada por huella al
///   genoma de partida.
/// - La reducción: se haría sobre la suministrada, sin el tracto.
fn effective_tract_insertion(
    general: &DeclaredGeneralConfig,
    sequences: &[EffectiveSequence],
) -> Result<Option<EffectiveTractInsertion>, Diagnostic> {
    const FIELD: &str = "reference.repeat_loci.insert_after";
    let reference = general.reference();
    let inserting = reference
        .repeat_loci
        .iter()
        .filter(|locus| locus.insert_after.is_some())
        .collect::<Vec<_>>();
    let Some(first) = inserting.first() else {
        if let Some(declared) = &reference.inserted_reference_id {
            return Err(config_error(
                general,
                "reference.inserted_reference_id",
                format!(
                    "inserted_reference_id = {} se declara sin ningún locus que inserte su tracto: solo nombra la referencia final de una inserción",
                    declared.value()
                ),
            ));
        }
        return Ok(None);
    };
    if !matches!(reference.source, DeclaredReferenceSource::Supplied { .. }) {
        return Err(repeat_error(
            general,
            FIELD,
            first.id.value(),
            "insert_after solo se admite con referencia suministrada: con referencia simulada el tracto se coloca en su tramo",
        ));
    }
    let Some(declared_id) = &reference.inserted_reference_id else {
        return Err(repeat_error(
            general,
            "reference.inserted_reference_id",
            first.id.value(),
            "el locus inserta su tracto y la referencia no declara inserted_reference_id: el genoma final no es el fichero suministrado y necesita identidad propia",
        ));
    };
    ReferenceId::parse(declared_id.value(), "configuration").map_err(|diagnostic| {
        diagnostic
            .with("path", general.source().path().display().to_string())
            .with("location", declared_id.location().to_string())
    })?;
    if declared_id.value() == reference.id.value() {
        return Err(config_error(
            general,
            "reference.inserted_reference_id",
            format!(
                "inserted_reference_id repite el identificador de la referencia suministrada ({}): el genoma final no es el fichero suministrado",
                declared_id.value()
            ),
        ));
    }
    if !general.direct_events().is_empty() {
        return Err(repeat_error(
            general,
            "reference.repeat_loci.insert_after,direct_events",
            first.id.value(),
            "insertar un tracto junto a cambios directos no se admite en esta versión: sus inserciones y deleciones se declaran en coordenadas vigentes de la copia, que no tienen traducción al genoma final",
        ));
    }
    match general.population_snvs() {
        Some(DeclaredPopulationSnvs::Supplied { .. }) => {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.insert_after,population_snvs.mode",
                first.id.value(),
                "insertar un tracto junto a una matriz de SNV suministrada no se admite en esta versión: sus posiciones y su huella son las del genoma de partida; use mode = \"generate\"",
            ));
        }
        Some(DeclaredPopulationSnvs::Reuse { .. }) => {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.insert_after,population_snvs.mode",
                first.id.value(),
                "insertar un tracto junto a una matriz de SNV reutilizada no se admite en esta versión: está atada por huella al genoma de partida; use mode = \"generate\"",
            ));
        }
        Some(DeclaredPopulationSnvs::Generate { .. }) | None => {}
    }
    if !general.reduced_references().is_empty() {
        return Err(repeat_error(
            general,
            "reference.repeat_loci.insert_after,reduced_references",
            first.id.value(),
            "insertar un tracto junto a una reducción de referencia en la misma ejecución no se admite: la reducción se haría sobre la suministrada, sin el tracto. Hágalo en dos ejecuciones: en la primera, solo la reducción, que publica el FASTA reducido; en la segunda, ese FASTA como reference.fasta, donde ya puede recibir la inserción",
        ));
    }
    let mut points: Vec<(usize, EffectiveInsertionPoint)> = Vec::with_capacity(inserting.len());
    for locus in &inserting {
        let locus_id = locus.id.value();
        let insert_after = locus
            .insert_after
            .as_ref()
            .map(|value| *value.value())
            .ok_or_else(|| {
                Diagnostic::internal_invariant("application", "tract_insertion")
                    .with("detail", "un locus filtrado por insert_after no lo tiene")
            })?;
        let (sequence_index, sequence) = sequences
            .iter()
            .enumerate()
            .find(|(_, sequence)| sequence.id() == locus.sequence.value())
            .ok_or_else(|| {
                repeat_error(
                    general,
                    "reference.repeat_loci.sequence",
                    locus_id,
                    format!("secuencia desconocida: {}", locus.sequence.value()),
                )
            })?;
        if insert_after == 0 || insert_after >= sequence.length() {
            return Err(repeat_error(
                general,
                FIELD,
                locus_id,
                format!(
                    "insert_after = {insert_after} no deja una base a cada lado en {} ({} bases): va de 1 a {}",
                    sequence.id(),
                    sequence.length(),
                    sequence.length().saturating_sub(1)
                ),
            ));
        }
        let Some(units) = &locus.reference_units else {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.reference_units",
                locus_id,
                "el locus inserta su tracto y no declara reference_units: sin tramo en el genoma de partida, las unidades de cada bloque son lo único que dice qué se inserta",
            ));
        };
        if units.len() != locus.motifs.len() {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.reference_units",
                locus_id,
                format!(
                    "reference_units declara {} valores y el locus tiene {} motivos: hace falta uno por bloque",
                    units.len(),
                    locus.motifs.len()
                ),
            ));
        }
        let max_units = locus
            .max_units
            .as_ref()
            .map_or(DEFAULT_REPEAT_MAX_UNITS, |declared| *declared.value());
        let mut total_units = 0_u64;
        let mut bases = Vec::new();
        for (block, (declared_units, declared_motif)) in units.iter().zip(&locus.motifs).enumerate()
        {
            let block_units = *declared_units.value();
            if block_units == 0 {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.reference_units",
                    locus_id,
                    format!(
                        "reference_units declara cero unidades en el bloque {block}: cada bloque lleva al menos una"
                    ),
                ));
            }
            total_units = total_units.saturating_add(block_units);
            // **Antes de construir nada**, como con los alelos: un error de
            // dedo de nueve cifras no puede reservar gigabytes.
            if total_units > max_units {
                return Err(repeat_error(
                    general,
                    "reference.repeat_loci.reference_units",
                    locus_id,
                    format!(
                        "el tracto insertado mide al menos {total_units} unidades y rebasa el tope de {max_units}"
                    ),
                ));
            }
            let motif = RepeatMotif::new(declared_motif.value().as_bytes()).map_err(|error| {
                repeat_error(
                    general,
                    "reference.repeat_loci.motif",
                    locus_id,
                    error.detail(),
                )
            })?;
            for _ in 0..block_units {
                bases.extend_from_slice(motif.bases());
            }
        }
        if let Some((_, other)) = points.iter().find(|(_, other)| {
            other.sequence_id == sequence.id() && other.insert_after == insert_after
        }) {
            return Err(repeat_error(
                general,
                FIELD,
                locus_id,
                format!(
                    "inserta su tracto en el mismo punto que el locus {} ({}:{insert_after}): el orden entre los dos no está definido",
                    other.locus_id,
                    sequence.id()
                ),
            ));
        }
        points.push((
            sequence_index,
            EffectiveInsertionPoint {
                locus_id: locus_id.clone(),
                sequence_id: sequence.id().to_owned(),
                insert_after,
                bases: Arc::from(bases),
            },
        ));
    }
    points.sort_by_key(|(index, point)| (*index, point.insert_after));
    let points = points
        .into_iter()
        .map(|(_, point)| point)
        .collect::<Vec<_>>();
    let mut derived_sequences = Vec::with_capacity(sequences.len());
    for sequence in sequences {
        let added = points
            .iter()
            .filter(|point| point.sequence_id == sequence.id)
            .map(EffectiveInsertionPoint::length)
            .sum::<u64>();
        let length = sequence.length.checked_add(added).ok_or_else(|| {
            config_error(
                general,
                "reference.repeat_loci.reference_units",
                format!(
                    "los tractos insertados desbordan la longitud de {}",
                    sequence.id
                ),
            )
        })?;
        derived_sequences.push(EffectiveSequence {
            id: sequence.id.clone(),
            class: sequence.class,
            length,
        });
    }
    Ok(Some(EffectiveTractInsertion {
        derived_reference_id: declared_id.value().clone(),
        sequences: derived_sequences,
        points,
    }))
}

/// Traduce al genoma final las regiones del generador de la matriz de SNV
/// (DTC-098). Una región que un tracto partiría se rechaza: sus bases dejarían
/// de ser contiguas.
fn translate_generation_regions(
    general: &DeclaredGeneralConfig,
    insertion: &EffectiveTractInsertion,
    population_snvs: Option<&mut EffectivePopulationSnvs>,
) -> Result<(), Diagnostic> {
    let Some(EffectivePopulationSnvs::Generate { parameters }) = population_snvs else {
        return Ok(());
    };
    let GeneratedSpatialMode::WeightedRegions { regions, .. } = &mut parameters.spatial else {
        return Ok(());
    };
    for region in regions {
        let (start, end) = insertion
            .translate_interval(&region.contig, region.start_1based, region.end_1based)
            .map_err(|inside| {
                config_error(
                    general,
                    "population_snvs.generate.regions,reference.repeat_loci.insert_after",
                    format!(
                        "la región {}:{}..{} contiene el punto de inserción del locus {} (tras la base {}): un tracto insertado la partiría",
                        region.contig,
                        region.start_1based,
                        region.end_1based,
                        inside.locus_id,
                        inside.insert_after
                    ),
                )
            })?;
        region.start_1based = start;
        region.end_1based = end;
    }
    Ok(())
}

/// **Guarda de admisión que cruza el evento técnico con los loci declarados.**
///
/// # Qué exige
///
/// Un locus que vaya a tartamudear necesita una tirada mínima **efectiva** de
/// al menos dos unidades: la suya si la declara, y si no la general. Si es 1,
/// la configuración se rechaza nombrando el locus, el valor declarado y el
/// mínimo exigido.
///
/// # Por qué
///
/// El paso de unidades no permite que una contracción vacíe su tirada. La
/// consecuencia es que **un tracto cuya tirada más larga sea de una sola copia
/// no puede contraerse nunca, mientras que expandirse sí puede**: una clase
/// entera de tractos quedaría congelada frente a la contracción con la
/// expansión libre, lo que **sesga el modelo hacia el alargamiento**. Hasta que
/// cada locus pudo declarar su propia tirada mínima esa clase no llegaba a
/// entrar, porque la general es 3; desde entonces un locus con tirada mínima 1
/// la admite. Esta guarda cierra esa puerta.
///
/// Lo que **no** cierra, y está bien: con una tirada de dos copias siguen
/// prohibidos los pasos de dos y tres unidades, porque la vaciarían. Esa
/// asimetría residual es **físicamente correcta** —en la realidad tampoco se
/// pueden perder tres unidades de una tirada que solo tiene dos— y convierte un
/// sesgo artificial del modelo en una restricción física real.
///
/// # Por qué aquí y no al leer el locus
///
/// El microsatélite se declara **junto a la referencia** y el tartamudeo se
/// configura **en el evento técnico**: son dos sitios distintos del fichero. La
/// condición «si va a tartamudear, exige 2» no puede comprobarse al leer el
/// locus, porque en ese momento todavía no se sabe si habrá tartamudeo. Solo
/// puede comprobarse cuando las dos configuraciones están resueltas, que es
/// aquí.
/// **La preparación amplificada sobre microsatélites exige el tartamudeo de
/// base.**
///
/// # Qué cambió aquí, y por qué no es el mismo rechazo
///
/// Hasta la entrega 4 esta guarda rechazaba **toda** combinación de
/// preparación amplificada con loci de repetición declarados, porque el
/// tartamudeo clonal heredable por linaje no existía. Ya existe, y ese rechazo
/// temporal está levantado. Lo que queda en su sitio es una condición
/// permanente y más estrecha: la combinación se admite **si el tartamudeo de
/// base está activo**, y solo entonces.
///
/// # Por qué la condición sigue haciendo falta
///
/// Las dos columnas de la calibración son totales. El tartamudeo clonal
/// **despeja** el incremento propio de la amplificación,
/// `r_pcr = 1 − (1−r_con)/(1−r_sin)`, y lo compone con la parte de base para
/// devolver `r_con`, la fracción medida. Sin la parte de base, una biblioteca
/// amplificada se quedaría en el incremento despejado y quedaría **por debajo
/// de lo medido**: con la calibración por omisión, entre un 14 % y un 28 % menos
/// en motivos de uno a cinco pares de base (con el multiplicador de los
/// homopolímeros de DTC-099, la cifra del homopolímero cambia). No es una
/// preferencia de configuración: es que la aritmética no cierra.
///
/// Reutiliza el código de diagnóstico de configuración que ya existe, en vez
/// de abrir uno nuevo: el manifiesto del recorrido general valida sus
/// diagnósticos contra el catálogo histórico congelado, y un código nuevo se
/// vería rechazado al publicar precisamente el fallo que debía informar. La
/// limitación está registrada aparte, en DT-007.
fn admit_amplified_preparation_against_loci(
    general: &DeclaredGeneralConfig,
    technical_event: Option<&EffectiveTechnicalEvent>,
    repeat_loci: &[EffectiveRepeatLocus],
) -> Result<(), Diagnostic> {
    let Some(event) = technical_event else {
        return Ok(());
    };
    if event.preparation() != LibraryPreparation::Pcr
        || repeat_loci.is_empty()
        || event.stutter().enabled()
    {
        return Ok(());
    }
    Err(config_error(
        general,
        "technical_event.library_preparation,technical_event.stutter.enabled",
        format!(
            "la preparación amplificada por PCR con loci de repetición declarados ({} \
             declarados, el primero es {}) exige technical_event.stutter.enabled = true: el \
             tartamudeo clonal de la PCR solo despeja el incremento de la amplificación, y sin \
             la parte de base la biblioteca subestimaría el tartamudeo (con la calibración por \
             omisión, entre un 14 % y un 28 % según la longitud del motivo). O se activa el \
             tartamudeo de base, o se quitan los \
             loci, o se usa library_preparation = \"{}\"",
            repeat_loci.len(),
            repeat_loci.first().map_or("", EffectiveRepeatLocus::id),
            LibraryPreparation::PcrFree.as_str()
        ),
    ))
}

/// Conecta el tartamudeo clonal al evento cuando hay amplificación **y** loci.
///
/// # Por qué aquí y no al leer la tabla de la amplificación
///
/// Los dos hechos que lo deciden viven en sitios distintos del fichero: la
/// amplificación se declara en el evento técnico y los microsatélites junto a
/// la referencia. Al leer `[technical_event.pcr]` todavía no se sabe si habrá
/// loci, exactamente igual que pasa con la guarda de la tirada mínima.
///
/// Sin loci no se conecta nada, y una biblioteca amplificada publica entonces
/// los mismos bytes que antes de esta entrega.
///
/// # El multiplicador de los homopolímeros se admite antes que nada
///
/// Su dominio y su combinación con la sustitución directa del homopolímero
/// (DTC-099) se comprueban **aunque no haya loci**: un valor imposible es un
/// error del autor tenga o no efecto, y rechazarlo no cambia nada de lo que
/// publica una configuración válida.
fn attach_clonal_stutter(
    general: &DeclaredGeneralConfig,
    technical_event: Option<&mut EffectiveTechnicalEvent>,
    microsatellites: bool,
) -> Result<(), Diagnostic> {
    let Some(event) = technical_event else {
        return Ok(());
    };
    let homopolymer_pcr_fold = declared_homopolymer_pcr_fold(general)?;
    if event.preparation != LibraryPreparation::Pcr || !microsatellites {
        return Ok(());
    }
    let substitutions = declared_per_cycle_slips(general)?;
    let Some(settings) = event.pcr.take() else {
        return Ok(());
    };
    let cycles = settings.cycles();
    let clonal = crate::clonal_stutter::ClonalStutterSettings::new(
        cycles,
        &substitutions,
        homopolymer_pcr_fold,
    )
    .map_err(|error| clonal_config_error(general, &error))?;
    event.pcr = Some(settings.with_clonal_stutter(clonal));
    Ok(())
}

/// Multiplicador de la PCR en los homopolímeros, ya leído y admitido.
///
/// Comprueba el dominio, `[1, 37]` y finito, y que no se declare junto a una
/// sustitución directa de la probabilidad por ciclo de los homopolímeros.
/// Para lo segundo mira solo la longitud de motivo de cada fila, sin leer su
/// probabilidad: la combinación es incompatible valga lo que valga.
fn declared_homopolymer_pcr_fold(
    general: &DeclaredGeneralConfig,
) -> Result<Option<f64>, Diagnostic> {
    let Some(stutter) = general
        .technical_event()
        .and_then(|event| event.pcr.as_ref())
        .and_then(|pcr| pcr.stutter.as_ref())
    else {
        return Ok(None);
    };
    let Some(declared) = stutter.homopolymer_pcr_fold.as_ref() else {
        return Ok(None);
    };
    let value = declared.literal().parse::<f64>().map_err(|_| {
        config_error(
            general,
            "technical_event.pcr.stutter.homopolymer_pcr_fold",
            "debe ser un f64 válido",
        )
    })?;
    let value = crate::clonal_stutter::validate_homopolymer_pcr_fold(value)
        .map_err(|error| clonal_config_error(general, &error))?;
    if stutter
        .per_cycle_slip
        .iter()
        .any(|entry| *entry.motif_length_bases.value() == 1)
    {
        return Err(clonal_config_error(
            general,
            &crate::clonal_stutter::ClonalStutterError::HomopolymerFoldWithSubstitution { value },
        ));
    }
    Ok(Some(value))
}

/// Sustituciones directas de la probabilidad por ciclo, ya leídas del TOML.
///
/// Una longitud de motivo declarada **dos veces** se rechaza en vez de dejar
/// que gane la última: dos valores para la misma fila son una contradicción
/// del autor, no una preferencia.
fn declared_per_cycle_slips(
    general: &DeclaredGeneralConfig,
) -> Result<std::collections::BTreeMap<u64, f64>, Diagnostic> {
    let mut substitutions = std::collections::BTreeMap::new();
    let Some(event) = general.technical_event() else {
        return Ok(substitutions);
    };
    let Some(pcr) = event.pcr.as_ref() else {
        return Ok(substitutions);
    };
    let Some(stutter) = pcr.stutter.as_ref() else {
        return Ok(substitutions);
    };
    for entry in &stutter.per_cycle_slip {
        let motif_length = *entry.motif_length_bases.value();
        let probability = entry.probability.literal().parse::<f64>().map_err(|_| {
            config_error(
                general,
                "technical_event.pcr.stutter.per_cycle_slip.probability",
                "debe ser un f64 válido",
            )
        })?;
        if substitutions.insert(motif_length, probability).is_some() {
            return Err(config_error(
                general,
                "technical_event.pcr.stutter.per_cycle_slip.motif_length_bases",
                format!(
                    "la longitud de motivo {motif_length} se declara dos veces; cada longitud \
                     admite una sola probabilidad por ciclo"
                ),
            ));
        }
    }
    Ok(substitutions)
}

/// Traduce el rechazo del modelo clonal a un diagnóstico que nombra el
/// parámetro y su valor.
fn clonal_config_error(
    general: &DeclaredGeneralConfig,
    error: &crate::clonal_stutter::ClonalStutterError,
) -> Diagnostic {
    use crate::clonal_stutter::ClonalStutterError as Error;
    let (field, detail) = match error {
        Error::UnknownMotifLength { motif_length } => (
            "technical_event.pcr.stutter.per_cycle_slip.motif_length_bases",
            format!(
                "la longitud de motivo {motif_length} está fuera del tramo calibrado, que va de \
                 1 a {}; no se extrapola",
                crate::stutter::STUTTER_MOTIF_LENGTH_MAX
            ),
        ),
        Error::InvalidProbability {
            motif_length,
            value,
        } => (
            "technical_event.pcr.stutter.per_cycle_slip.probability",
            format!(
                "la probabilidad por ciclo {value} de la longitud de motivo {motif_length} debe \
                 ser finita y estar en [0, 1]; es la probabilidad de UN ciclo, no la fracción de \
                 extremo a extremo"
            ),
        ),
        Error::InvalidHomopolymerFold { value } => (
            "technical_event.pcr.stutter.homopolymer_pcr_fold",
            format!(
                "el multiplicador {value} de los homopolímeros debe ser finito y estar entre {} y \
                 {}, ambos incluidos; por debajo de 1 diría que la PCR quita tartamudeo, y por \
                 encima de 37 la fracción con PCR, multiplicador × {}, llegaría a 1",
                crate::clonal_stutter::HOMOPOLYMER_PCR_FOLD_MINIMUM,
                crate::clonal_stutter::HOMOPOLYMER_PCR_FOLD_MAXIMUM,
                crate::stutter::STUTTER_DEFAULT_MOTIF_LENGTH_RATES[0]
            ),
        ),
        Error::HomopolymerFoldWithSubstitution { value } => (
            "technical_event.pcr.stutter.homopolymer_pcr_fold,technical_event.pcr.stutter.per_cycle_slip.motif_length_bases",
            format!(
                "el multiplicador {value} de los homopolímeros no puede declararse junto a una \
                 sustitución directa de la probabilidad por ciclo de la longitud de motivo 1: \
                 serían dos valores para la misma probabilidad. O uno, o el otro"
            ),
        ),
    };
    config_error(general, field, detail)
}

fn admit_loci_against_stutter(
    general: &DeclaredGeneralConfig,
    technical_event: Option<&EffectiveTechnicalEvent>,
    repeat_loci: &[EffectiveRepeatLocus],
) -> Result<(), Diagnostic> {
    let Some(event) = technical_event else {
        return Ok(());
    };
    // **Una sustitución que nombre un locus inexistente se rechaza.** Sin esta
    // guarda se aceptaba y se tiraba en silencio: el usuario declaraba una tasa
    // propia, no ocurría nada y nadie se enteraba. Un identificador mal escrito
    // es exactamente la forma que toma ese error.
    let declared_locus_ids = repeat_loci
        .iter()
        .map(EffectiveRepeatLocus::id)
        .collect::<std::collections::BTreeSet<_>>();
    for locus_id in event.stutter().locus_rates().keys() {
        if !declared_locus_ids.contains(locus_id.as_str()) {
            return Err(config_error(
                general,
                "technical_event.stutter.locus_rates.locus_id",
                format!(
                    "la sustitución de tasa nombra el locus {locus_id}, que no está declarado en \
                     reference.repeat_loci; los declarados son {}",
                    declared_locus_ids
                        .iter()
                        .copied()
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ));
        }
    }
    // Lo que sigue solo tiene sentido si el tartamudeo va a ocurrir. La
    // comprobación de arriba, en cambio, vale siempre: un identificador mal
    // escrito lo está con el interruptor encendido y con él apagado.
    if !event.stutter().enabled() {
        return Ok(());
    }
    const MINIMUM_LONGEST_RUN_UNITS: u64 = 2;
    for locus in repeat_loci {
        let declared = locus.admission().min_longest_run_units();
        if declared < MINIMUM_LONGEST_RUN_UNITS {
            return Err(config_error(
                general,
                "reference.repeat_loci.min_longest_run_units,technical_event.stutter",
                format!(
                    "el locus {} declara una tirada mínima efectiva de {declared} unidad(es) y el \
                     tartamudeo exige al menos {MINIMUM_LONGEST_RUN_UNITS}: una tirada de una sola \
                     copia no puede contraerse nunca y el modelo quedaría sesgado hacia el \
                     alargamiento",
                    locus.id()
                ),
            ));
        }
        for motif in locus.motifs() {
            let length = motif.length();
            if event
                .stutter()
                .motif_length_rates()
                .rate_for_motif_length(length)
                .is_none()
            {
                return Err(config_error(
                    general,
                    "reference.repeat_loci.motif,technical_event.stutter.motif_length_rates",
                    format!(
                        "el locus {} declara un motivo de {length} pb y la calibración del \
                         tartamudeo solo cubre de 1 a {} pb; no se extrapola",
                        locus.id(),
                        crate::stutter::STUTTER_MOTIF_LENGTH_MAX
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// Calcula la cota de bases de tracto y **rechaza el escenario** si supera el
/// tope declarado.
///
/// Se hace aquí, en la frontera de configuración, antes de materializar a
/// nadie: el rechazo tiene que llegar antes de gastar una sola copia.
///
/// # Lo que esta función **no** hace
///
/// No cuenta nada de lo publicado. Su aritmética sale entera de lo declarado
/// —alelos, longitudes, individuos y ploidía— y no comparte ni una línea con
/// el contador de lo realizado. Que sean dos caminos distintos es lo único que
/// hace que compararlos signifique algo; si uno llamara al otro, la
/// comprobación sería un espejo incapaz de fallar.
fn repeat_output_budget(
    general: &DeclaredGeneralConfig,
    loci: &[EffectiveRepeatLocus],
    populations: &[EffectiveInitialPopulation],
) -> Result<EffectiveRepeatOutputBudget, Diagnostic> {
    let copies = populations.iter().try_fold(0_u64, |total, population| {
        u64::try_from(population.individual_count())
            .ok()
            .and_then(|individuals| {
                individuals.checked_mul(u64::from(population.ploidy().copy_count()))
            })
            .and_then(|copies| total.checked_add(copies))
            .ok_or_else(|| {
                config_error(
                    general,
                    "initial_populations",
                    "el número de copias haplotípicas desborda el recuento",
                )
            })
    })?;
    let declared_cap = general
        .reference()
        .max_tract_output_bases
        .as_ref()
        .map(|declared| *declared.value());
    let budget = recompute_repeat_output_budget(loci, copies, declared_cap);
    if let Some(cap) = declared_cap
        && budget.total_bases() > cap
    {
        let heaviest = budget.heaviest_locus().map_or_else(
            || "ninguno".to_owned(),
            |locus| format!("{} con {} bases", locus.locus_id(), locus.truth_bases()),
        );
        return Err(config_error(
            general,
            "reference.max_tract_output_bases",
            format!(
                "solo la verdad publicaría hasta {} bases de tracto y el tope declarado es {cap}; el locus que más aporta es {heaviest}",
                budget.total_bases()
            ),
        ));
    }
    Ok(budget)
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
    let stutter = effective_stutter(general, declared.stutter.as_ref())?;
    let (preparation, pcr) = effective_pcr(general, declared)?;
    Ok(Some(EffectiveTechnicalEvent {
        event_id,
        individual_ids,
        target,
        sequencing_errors,
        read_length,
        chimera_rate,
        index_hopping_rate,
        biological_contamination,
        stutter,
        preparation,
        pcr,
        libraries,
    }))
}

/// Resuelve `[technical_event.library_preparation]` y `[technical_event.pcr]`.
///
/// # Las tres reglas que gobiernan esta tabla
///
/// 1. La tabla **solo** se admite con preparación `pcr`. Presente con la
///    preparación sin amplificar, se rechaza: aceptarla en silencio dejaría al
///    usuario creyendo que declaró una amplificación que nadie iba a hacer.
/// 2. Con preparación `pcr` la tabla es **obligatoria**, y dentro de ella lo
///    son la fracción de duplicados, los ciclos y la eficiencia. Ninguno tiene
///    valor por omisión.
/// 3. Todo dominio se comprueba **antes** de simular, nombrando parámetro y
///    valor.
fn effective_pcr(
    general: &DeclaredGeneralConfig,
    declared: &crate::configuration::DeclaredTechnicalEvent,
) -> Result<(LibraryPreparation, Option<PcrSettings>), Diagnostic> {
    let preparation = match declared.library_preparation.as_ref() {
        None => LibraryPreparation::PcrFree,
        Some(value) => LibraryPreparation::parse(value.value()).ok_or_else(|| {
            config_error(
                general,
                "technical_event.library_preparation",
                format!(
                    "la preparación {} no existe; las admitidas son {} y {}",
                    value.value(),
                    LibraryPreparation::PcrFree.as_str(),
                    LibraryPreparation::Pcr.as_str()
                ),
            )
        })?,
    };
    let Some(table) = declared.pcr.as_ref() else {
        if preparation == LibraryPreparation::Pcr {
            return Err(config_error(
                general,
                "technical_event.pcr",
                "la preparación pcr exige la tabla technical_event.pcr con la fracción de \
                 duplicados esperada, los ciclos y la eficiencia por ciclo, que son obligatorios \
                 y no tienen valor por omisión",
            ));
        }
        return Ok((preparation, None));
    };
    if preparation != LibraryPreparation::Pcr {
        return Err(config_error(
            general,
            "technical_event.pcr",
            format!(
                "la tabla technical_event.pcr solo se admite con library_preparation = \"{}\"; \
                 la declarada es \"{}\"",
                LibraryPreparation::Pcr.as_str(),
                preparation.as_str()
            ),
        ));
    }
    let cycles = u32::try_from(*declared_cycles(table)).map_err(|_| {
        config_error(
            general,
            "technical_event.pcr.cycles",
            format!(
                "los ciclos deben ser un entero entre 1 y {}",
                crate::technical_profile::PCR_MAXIMUM_CYCLES
            ),
        )
    })?;
    let parameters = crate::pcr_amplification::DeclaredPcrParameters {
        duplicate_fraction: parse_pcr_number(
            general,
            &table.duplicate_fraction,
            "technical_event.pcr.duplicate_fraction",
        )?,
        cycles,
        efficiency: parse_pcr_number(general, &table.efficiency, "technical_event.pcr.efficiency")?,
        early_cycles: optional_pcr_count(
            general,
            table.early_cycles.as_ref(),
            "technical_event.pcr.early_cycles",
        )?,
        max_family_size: optional_pcr_count(
            general,
            table.max_family_size.as_ref(),
            "technical_event.pcr.max_family_size",
        )?,
        gc_max_attempts: optional_pcr_count(
            general,
            table.gc_max_attempts.as_ref(),
            "technical_event.pcr.gc_max_attempts",
        )?,
        gc_curve: if table.gc_relative_efficiency.is_empty() {
            None
        } else {
            Some(
                table
                    .gc_relative_efficiency
                    .iter()
                    .map(|knot| {
                        Ok((
                            parse_pcr_number(
                                general,
                                &knot.gc_fraction,
                                "technical_event.pcr.gc_relative_efficiency.gc_fraction",
                            )?,
                            parse_pcr_number(
                                general,
                                &knot.relative_efficiency,
                                "technical_event.pcr.gc_relative_efficiency.relative_efficiency",
                            )?,
                        ))
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()?,
            )
        },
        efficiency_classes: table
            .efficiency_classes
            .iter()
            .map(|class| {
                Ok((
                    parse_pcr_number(
                        general,
                        &class.efficiency,
                        "technical_event.pcr.efficiency_classes.efficiency",
                    )?,
                    parse_pcr_number(
                        general,
                        &class.weight,
                        "technical_event.pcr.efficiency_classes.weight",
                    )?,
                ))
            })
            .collect::<Result<Vec<_>, Diagnostic>>()?,
    };
    let tolerance = effective_normalization_tolerance(general)?;
    let settings = PcrSettings::new(&parameters, tolerance)
        .map_err(|error| pcr_config_error(general, &error))?;
    Ok((preparation, Some(settings)))
}

const fn declared_cycles(table: &crate::configuration::DeclaredPcr) -> &u64 {
    table.cycles.value()
}

fn parse_pcr_number(
    general: &DeclaredGeneralConfig,
    value: &crate::configuration::DeclaredNumber,
    field: &str,
) -> Result<f64, Diagnostic> {
    value
        .literal()
        .parse::<f64>()
        .map_err(|_| config_error(general, field, "debe ser un f64 válido"))
}

fn optional_pcr_count(
    general: &DeclaredGeneralConfig,
    value: Option<&crate::configuration::Located<u64>>,
    field: &str,
) -> Result<Option<u32>, Diagnostic> {
    value
        .map(|located| {
            u32::try_from(*located.value())
                .map_err(|_| config_error(general, field, "está fuera del dominio del contrato"))
        })
        .transpose()
}

/// Traduce el rechazo del modelo puro a un diagnóstico que nombra el parámetro
/// y su valor, que es lo que permite corregir la configuración sin adivinar.
fn pcr_config_error(
    general: &DeclaredGeneralConfig,
    error: &crate::pcr_amplification::PcrAmplificationError,
) -> Diagnostic {
    use crate::pcr_amplification::PcrAmplificationError as Error;
    let (field, detail) = match error {
        Error::InvalidDuplicateFraction { value } => (
            "technical_event.pcr.duplicate_fraction",
            format!(
                "la fracción de duplicados {value} debe ser finita y estar en [0, 1); es la \
                 fracción ESPERADA de una biblioteca homogénea, no la realizada, y no tiene \
                 valor por omisión"
            ),
        ),
        Error::UnreachableDuplicateFraction {
            value,
            maximum_reachable,
            maximum_mean_pairs_per_molecule,
        } => (
            "technical_event.pcr.duplicate_fraction",
            format!(
                "la fracción de duplicados {value} está por encima de la máxima alcanzable \
                 {maximum_reachable}, que es la que produce el tope de \
                 {maximum_mean_pairs_per_molecule} parejas por molécula; no se pinza en silencio"
            ),
        ),
        Error::InvalidCycles { value, maximum } => (
            "technical_event.pcr.cycles",
            format!("los ciclos {value} deben ser un entero entre 1 y {maximum}"),
        ),
        Error::InvalidEfficiency { value } => (
            "technical_event.pcr.efficiency",
            format!("la eficiencia por ciclo {value} debe estar en (0, 1]"),
        ),
        Error::InvalidEarlyCycles {
            value,
            maximum,
            cycles,
        } => (
            "technical_event.pcr.early_cycles",
            format!(
                "los ciclos explícitos {value} deben estar entre 0 y {maximum} y no superar los \
                 {cycles} ciclos declarados"
            ),
        ),
        Error::InvalidMaxFamilySize { value, maximum } => (
            "technical_event.pcr.max_family_size",
            format!("el tamaño máximo de familia {value} debe estar entre 1 y {maximum}"),
        ),
        Error::InvalidGcAttempts { value, maximum } => (
            "technical_event.pcr.gc_max_attempts",
            format!("los intentos de composición {value} deben estar entre 1 y {maximum}"),
        ),
        Error::InvalidGcCurveKnot {
            gc_fraction,
            relative_efficiency,
        } => (
            "technical_event.pcr.gc_relative_efficiency",
            format!(
                "el nudo de fracción {gc_fraction} con eficiencia relativa \
                 {relative_efficiency} está fuera de dominio: la fracción va en [0, 1] y la \
                 eficiencia relativa en (0, 1]"
            ),
        ),
        Error::NonIncreasingGcCurveKnots { previous, found } => (
            "technical_event.pcr.gc_relative_efficiency",
            format!(
                "los nudos de la curva deben ser estrictamente crecientes: tras {previous} llega \
                 {found}"
            ),
        ),
        Error::ShortGcCurve { knots } => (
            "technical_event.pcr.gc_relative_efficiency",
            format!("la curva declara {knots} nudos y necesita al menos dos"),
        ),
        Error::TooManyEfficiencyClasses { value, maximum } => (
            "technical_event.pcr.efficiency_classes",
            format!("se declaran {value} clases de eficiencia y el máximo es {maximum}"),
        ),
        Error::InvalidEfficiencyClass { efficiency, weight } => (
            "technical_event.pcr.efficiency_classes",
            format!(
                "la clase de eficiencia {efficiency} con peso {weight} está fuera de dominio: la \
                 eficiencia va en (0, 1] y el peso no puede ser negativo"
            ),
        ),
        Error::EfficiencyClassWeightsNotNormalized { sum, tolerance } => (
            "technical_event.pcr.efficiency_classes",
            format!(
                "los pesos de las clases de eficiencia suman {sum} y deben sumar 1 dentro de la \
                 tolerancia de normalización {tolerance}"
            ),
        ),
        Error::NumericFailure => (
            "technical_event.pcr",
            "la derivación de los parámetros de amplificación falló numéricamente".to_owned(),
        ),
    };
    config_error(general, field, detail)
}

/// Resuelve `[technical_event.stutter]` a los parámetros efectivos.
///
/// # Qué se configura y qué no
///
/// La **estructura** del paso no se configura: la fija el algoritmo. Lo que se
/// declara son los **pesos relativos** dentro de esa estructura, y son
/// **obligatorios y sin valor por omisión** cuando el tartamudeo está activo,
/// porque la fuente publica el método de ajuste pero no los valores
/// (pendiente P-3). El **coeficiente de la tirada** arranca en cero y se
/// publica como no calibrado (pendiente P-1).
fn effective_stutter(
    general: &DeclaredGeneralConfig,
    declared: Option<&DeclaredStutter>,
) -> Result<StutterSettings, Diagnostic> {
    let Some(declared) = declared else {
        return Ok(StutterSettings::default());
    };
    let enabled = declared
        .enabled
        .as_ref()
        .is_some_and(|value| *value.value());
    let mut rates = crate::stutter::STUTTER_DEFAULT_MOTIF_LENGTH_RATES;
    let mut seen_lengths = std::collections::BTreeSet::new();
    for entry in &declared.motif_length_rates {
        let length = *entry.motif_length_bases.value();
        let index = usize::try_from(length)
            .ok()
            .filter(|length| (1..=crate::stutter::STUTTER_MOTIF_LENGTH_MAX).contains(length))
            .ok_or_else(|| {
                config_error(
                    general,
                    "technical_event.stutter.motif_length_rates.motif_length_bases",
                    format!(
                        "la longitud de motivo {length} está fuera del tramo calibrado 1..={}",
                        crate::stutter::STUTTER_MOTIF_LENGTH_MAX
                    ),
                )
            })?;
        if !seen_lengths.insert(length) {
            return Err(config_error(
                general,
                "technical_event.stutter.motif_length_rates.motif_length_bases",
                format!("la longitud de motivo {length} se declara dos veces"),
            ));
        }
        rates[index - 1] = parse_stutter_number(
            general,
            &entry.rate,
            "technical_event.stutter.motif_length_rates.rate",
        )?;
    }
    let rates = StutterMotifLengthRates::new(rates).map_err(|error| {
        stutter_config_error(
            general,
            "technical_event.stutter.motif_length_rates",
            &error,
        )
    })?;
    let coefficient = declared.run_length_coefficient.as_ref().map_or(
        Ok(crate::stutter::STUTTER_DEFAULT_RUN_LENGTH_COEFFICIENT),
        |value| {
            parse_stutter_number(
                general,
                value,
                "technical_event.stutter.run_length_coefficient",
            )
        },
    )?;
    let weights = if declared.steps.is_empty() {
        None
    } else {
        let mut resolved = [0.0_f64; crate::stutter::STUTTER_ADMITTED_STEPS.len()];
        let mut seen = std::collections::BTreeSet::new();
        for step in &declared.steps {
            let units = *step.units.value();
            let position = crate::stutter::STUTTER_ADMITTED_STEPS
                .iter()
                .position(|admitted| *admitted == units)
                .ok_or_else(|| {
                    config_error(
                        general,
                        "technical_event.stutter.steps.units",
                        format!(
                            "el paso {units} no pertenece a la estructura del algoritmo, que es \
                             exactamente {:?} en unidades enteras",
                            crate::stutter::STUTTER_ADMITTED_STEPS
                        ),
                    )
                })?;
            if !seen.insert(units) {
                return Err(config_error(
                    general,
                    "technical_event.stutter.steps.units",
                    format!("el paso {units} se declara dos veces"),
                ));
            }
            resolved[position] = parse_stutter_number(
                general,
                &step.weight,
                "technical_event.stutter.steps.weight",
            )?;
        }
        if seen.len() != crate::stutter::STUTTER_ADMITTED_STEPS.len() {
            return Err(config_error(
                general,
                "technical_event.stutter.steps",
                format!(
                    "hay que declarar el peso de los cuatro pasos admitidos {:?}; faltan {}",
                    crate::stutter::STUTTER_ADMITTED_STEPS,
                    crate::stutter::STUTTER_ADMITTED_STEPS
                        .iter()
                        .filter(|admitted| !seen.contains(*admitted))
                        .map(i64::to_string)
                        .collect::<Vec<_>>()
                        .join(", ")
                ),
            ));
        }
        Some(StutterStepWeights::new(resolved).map_err(|error| {
            stutter_config_error(general, "technical_event.stutter.steps.weight", &error)
        })?)
    };
    let mut locus_rates = std::collections::BTreeMap::new();
    for entry in &declared.locus_rates {
        let locus_id = entry.locus_id.value().clone();
        let rate = parse_stutter_number(
            general,
            &entry.rate,
            "technical_event.stutter.locus_rates.rate",
        )?;
        if locus_rates.insert(locus_id.clone(), rate).is_some() {
            return Err(config_error(
                general,
                "technical_event.stutter.locus_rates.locus_id",
                format!("el locus {locus_id} se sustituye dos veces"),
            ));
        }
    }
    StutterSettings::new(enabled, rates, coefficient, weights, locus_rates)
        .map_err(|error| stutter_config_error(general, "technical_event.stutter", &error))
}

fn parse_stutter_number(
    general: &DeclaredGeneralConfig,
    value: &crate::configuration::DeclaredNumber,
    field: &str,
) -> Result<f64, Diagnostic> {
    value
        .literal()
        .parse::<f64>()
        .map_err(|_| config_error(general, field, "debe ser un f64 válido"))
}

/// Traduce el rechazo del modelo puro a un diagnóstico que nombra el parámetro
/// y su valor, que es lo que permite corregir la configuración sin adivinar.
fn stutter_config_error(
    general: &DeclaredGeneralConfig,
    field: &str,
    error: &StutterError,
) -> Diagnostic {
    let detail = match error {
        StutterError::InvalidRate { motif_length, rate } => format!(
            "la tasa {rate} del motivo de {motif_length} pb debe ser una fracción finita entre 0 y 1"
        ),
        StutterError::InvalidLocusRate { locus_id, rate } => format!(
            "la tasa {rate} sustituida para el locus {locus_id} debe ser una fracción finita entre 0 y 1"
        ),
        StutterError::InvalidStepWeight { step, weight } => {
            format!("el peso {weight} del paso {step} debe ser finito y no negativo")
        }
        StutterError::InvalidRunLengthCoefficient { coefficient } => {
            format!("el coeficiente de la tirada {coefficient} debe ser finito y no negativo")
        }
        StutterError::NoPositiveStepWeight => {
            "al menos uno de los cuatro pesos tiene que ser positivo".to_owned()
        }
        StutterError::MissingStepWeights => {
            // **El rechazo trae el ejemplo a propósito.** Sin él, quien tropieza
            // aquí no sabe por dónde empezar y el ejemplo ilustrativo se queda
            // enterrado en el plan. Lo que no se hace es darle valor por
            // omisión: el coeficiente puede arrancar en cero porque el cero es
            // neutro, pero **ningún reparto de pesos lo es** —cualquiera afirma
            // algo sobre la biología—, así que el usuario tiene que mirarlo.
            "technical_event.stutter.steps es obligatorio cuando el tartamudeo está activo y no \
             tiene valor por omisión: la fuente publica el método de ajuste, no los valores, y \
             ningún reparto es neutro. A modo ILUSTRATIVO, sin calibrar, un reparto compatible \
             con la estructura sería -1 con peso 100, +1 con peso 20, -2 con peso 10 y -3 con \
             peso 1"
                .to_owned()
        }
        other => format!("{other:?}"),
    };
    config_error(general, field, detail)
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

/// El intervalo de un locus, en base uno, en el **marco de la simulación**:
/// el genoma final si hay inserción y el declarado si no (DTC-098).
///
/// Devuelve además lo declarado, solo cuando hay inserción, que es cuando lo
/// declarado y el intervalo dejan de coincidir.
fn locus_interval_in_frame(
    general: &DeclaredGeneralConfig,
    locus: &crate::configuration::DeclaredRepeatLocus,
    sequence: &EffectiveSequence,
    insertion: Option<&EffectiveTractInsertion>,
) -> Result<(u64, u64, Option<DeclaredLocusPosition>), Diagnostic> {
    if let Some(insert_after) = &locus.insert_after {
        let insert_after = *insert_after.value();
        let (Some(plan), Some(point)) = (
            insertion,
            insertion.and_then(|plan| {
                plan.points()
                    .iter()
                    .find(|point| point.locus_id() == locus.id.value().as_str())
            }),
        ) else {
            return Err(
                Diagnostic::internal_invariant("application", "tract_insertion").with(
                    "detail",
                    "un locus que inserta no está en el plan de inserción",
                ),
            );
        };
        let start = plan.inserted_start(sequence.id(), insert_after);
        let end = start + point.length() - 1;
        return Ok((
            start,
            end,
            Some(DeclaredLocusPosition::InsertAfter(insert_after)),
        ));
    }
    let (Some(start), Some(end)) = (&locus.start, &locus.end) else {
        return Err(
            Diagnostic::internal_invariant("application", "repeat_locus_interval").with(
                "detail",
                "un locus sin insert_after llega sin start o sin end",
            ),
        );
    };
    let (start, end) = (*start.value(), *end.value());
    let Some(plan) = insertion else {
        return Ok((start, end, None));
    };
    // Primero contra el genoma de partida, que es donde se declaró: un tramo
    // que no cabe allí no puede traducirse a ningún sitio.
    if start == 0 || start > end || end > sequence.length() {
        return Err(repeat_error(
            general,
            "reference.repeat_loci.start,reference.repeat_loci.end",
            locus.id.value(),
            format!(
                "el intervalo {start}..{end} no cabe en {} ({} bases del genoma de partida)",
                sequence.id(),
                sequence.length()
            ),
        ));
    }
    let (final_start, final_end) = plan
        .translate_interval(sequence.id(), start, end)
        .map_err(|inside| {
            repeat_error(
                general,
                "reference.repeat_loci.start,reference.repeat_loci.end,reference.repeat_loci.insert_after",
                locus.id.value(),
                format!(
                    "el locus {} inserta su tracto tras la base {} de {}, dentro del tracto {start}..{end}, y lo partiría",
                    inside.locus_id(),
                    inside.insert_after(),
                    sequence.id()
                ),
            )
        })?;
    Ok((
        final_start,
        final_end,
        Some(DeclaredLocusPosition::Interval { start, end }),
    ))
}

/// Intervalos `[inicio, fin)` de los loci ya admitidos, por secuencia y por
/// inicio, con la posición de cada uno en el orden de declaración.
type OccupiedLoci =
    std::collections::BTreeMap<SequenceId, std::collections::BTreeMap<u64, (u64, usize)>>;

/// El locus admitido **antes**, en orden de declaración, de entre los que se
/// solapan con `interval`; `None` si ninguno lo hace.
///
/// Los intervalos admitidos no se solapan entre sí —el primero que lo hiciera
/// se habría rechazado—, de modo que, ordenados por inicio, también lo están
/// por fin. Por eso los que pisan `interval` son, como mucho, el último que
/// empieza antes que él (si acaba después de su inicio) y todos los que
/// empiezan dentro. Es el mismo conjunto que comparar con todos, sin
/// recorrerlos: con miles de loci, aquella comparación era cuadrática.
fn earliest_overlapping_locus(
    occupied: &OccupiedLoci,
    interval: &GenomicInterval,
) -> Option<usize> {
    let starts = occupied.get(interval.sequence_id())?;
    let before = starts
        .range(..interval.start())
        .next_back()
        .filter(|(_, (end, _))| *end > interval.start())
        .map(|(_, (_, index))| *index);
    let inside = starts
        .range(interval.start()..interval.end())
        .map(|(_, (_, index))| *index);
    before.into_iter().chain(inside).min()
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
    insertion: Option<&EffectiveTractInsertion>,
    common_slippage: Option<SlippageSettings>,
) -> Result<Vec<EffectiveRepeatLocus>, Diagnostic> {
    // **Con inserción, los intervalos viven en el genoma final** (DTC-098):
    // llevan su identificador y miden contra su catálogo. Sin ella, todo es
    // como siempre.
    let frame_reference_id = insertion.map_or(general.reference().id.value().as_str(), |plan| {
        plan.derived_reference_id()
    });
    let reference_id =
        ReferenceId::parse(frame_reference_id, "execution").map_err(|diagnostic| {
            diagnostic.with("path", general.source().path().display().to_string())
        })?;
    let frame_sequences = insertion.map_or(sequences, EffectiveTractInsertion::sequences);
    let mut resolved: Vec<EffectiveRepeatLocus> =
        Vec::with_capacity(general.reference().repeat_loci.len());
    // Los intervalos ya admitidos, por secuencia y por inicio, con su orden de
    // declaración: la comprobación de solapes mira solo los vecinos del
    // intervalo nuevo en vez de recorrer todos los anteriores (F6.4.a.6).
    let mut occupied = OccupiedLoci::new();
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
        let frame_length = frame_sequences
            .iter()
            .find(|frame| frame.id() == sequence.id())
            .map_or(sequence.length(), EffectiveSequence::length);
        let (start, end, declared_position) =
            locus_interval_in_frame(general, locus, sequence, insertion)?;
        // El intervalo pasa por el servicio común de coordenadas, que es quien
        // convierte de base uno cerrada a base cero semiabierta y quien
        // comprueba que cabe en la secuencia. Aquí no se resta ni se suma uno.
        let closed = ClosedInterval1::new(start, end).map_err(|error| {
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
            frame_length,
        )
        .map_err(|error| {
            repeat_error(
                general,
                "reference.repeat_loci.start,reference.repeat_loci.end",
                locus.id.value(),
                format!(
                    "el intervalo {start}..{end} no cabe en {} ({frame_length} bases): {error:?}",
                    sequence.id(),
                ),
            )
        })?;
        // Dos loci que se solapan describirían el mismo tracto dos veces y
        // producirían dos alelos sobre las mismas bases. El diagnóstico nombra
        // el **primero en orden de declaración** de los anteriores que se
        // solapan, como cuando se comparaba con todos uno a uno.
        if let Some(previous) = earliest_overlapping_locus(&occupied, &interval) {
            return Err(repeat_error(
                general,
                "reference.repeat_loci.start,reference.repeat_loci.end",
                locus.id.value(),
                format!(
                    "su intervalo se solapa con el del locus {}",
                    resolved[previous].id()
                ),
            ));
        }
        // Cada vuelta del bucle publica exactamente un locus en `resolved` o
        // termina con un error, así que su posición es `resolved.len()`.
        occupied
            .entry(interval.sequence_id().clone())
            .or_default()
            .insert(interval.start(), (interval.end(), resolved.len()));
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
        let locus_admission = resolve_admission(
            general,
            locus.admission.as_ref(),
            admission,
            "reference.repeat_loci",
            Some(locus.id.value()),
        )?;
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

        // **El generador se valida aquí, antes de simular.** Un modelo distinto
        // del escalonado o un salto de varias unidades se rechazan nombrando el
        // parámetro y el valor, sin haber gastado una sola copia.
        let slippage = match &locus.slippage {
            // **Hereda el deslizamiento común** (REQ-CONF-139): un declarado
            // sin alelos ni generador propio solo existe con la sección de
            // microsatélites, que el lector ya exigió. Como un generador
            // propio, solo admite un motivo.
            None if locus.alleles.is_empty() => {
                let Some(settings) = common_slippage else {
                    return Err(Diagnostic::internal_invariant(
                        "application",
                        "repeat_locus_common_slippage",
                    )
                    .with(
                        "detail",
                        "un locus sin alelos ni generador llegó sin la sección de microsatélites",
                    ));
                };
                if motifs.len() != 1 {
                    return Err(crate::repeat_slippage::SlippageError::CompositeLocus {
                        motifs: motifs.len(),
                    }
                    .into_diagnostic(
                        &general.source().path().display().to_string(),
                        locus.id.value(),
                    ));
                }
                Some(EffectiveRepeatSlippage {
                    settings,
                    generation: None,
                })
            }
            None => None,
            Some(declared) => {
                let settings = SlippageSettings::new(
                    declared.model.value().as_str(),
                    declared
                        .step_units
                        .as_ref()
                        .map_or(crate::repeat_slippage::SLIPPAGE_STEP_UNITS, |value| {
                            *value.value()
                        }),
                    parse_slippage_number(
                        general,
                        locus.id.value(),
                        "mean_events",
                        &declared.mean_events,
                    )?,
                    match &declared.expansion_bias {
                        None => crate::repeat_slippage::SLIPPAGE_DEFAULT_EXPANSION_BIAS,
                        Some(value) => parse_slippage_number(
                            general,
                            locus.id.value(),
                            "expansion_bias",
                            value,
                        )?,
                    },
                )
                .map_err(|error| {
                    error.into_diagnostic(
                        &general.source().path().display().to_string(),
                        locus.id.value(),
                    )
                })?;
                // Sin fronteras de bloque declaradas no hay forma de partir el
                // tracto de la referencia, de modo que el generador solo admite
                // un motivo. El rechazo llega aquí, no al generar.
                if motifs.len() != 1 {
                    return Err(crate::repeat_slippage::SlippageError::CompositeLocus {
                        motifs: motifs.len(),
                    }
                    .into_diagnostic(
                        &general.source().path().display().to_string(),
                        locus.id.value(),
                    ));
                }
                Some(EffectiveRepeatSlippage {
                    settings,
                    generation: None,
                })
            }
        };
        // **Después del generador, no antes.** Un compuesto que pide el
        // generador se rechaza porque el generador solo admite un motivo; pedirle
        // antes `reference_units` mandaba al usuario a declarar un campo que,
        // una vez declarado, no evitaba el rechazo.
        let reference_units = resolve_reference_units(general, locus, &declaration)?;
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
        if slippage.is_some() {
            // Un locus generado no trae frecuencias que normalizar: la
            // auditoría se recomputa sobre el abanico, cuando exista.
            resolved.push(EffectiveRepeatLocus {
                declaration,
                reference_units,
                declared_position,
                max_units,
                admission: locus_admission,
                slippage,
                alleles: Vec::new(),
                normalization: RepeatFrequencyNormalization {
                    tolerance,
                    original_sum: 0.0,
                    original_delta: 0.0,
                    effective_sum: 0.0,
                    effective_delta: 0.0,
                },
                origin: RepeatLocusOrigin::Declared,
                detector_score: None,
            });
            continue;
        }
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
                    origin: RepeatAlleleOrigin::Declared,
                },
            )
            .collect();
        resolved.push(EffectiveRepeatLocus {
            declaration,
            reference_units,
            declared_position,
            max_units,
            admission: locus_admission,
            slippage,
            alleles,
            normalization,
            origin: RepeatLocusOrigin::Declared,
            detector_score: None,
        });
    }
    Ok(resolved)
}

/// **Las unidades por bloque del tracto de una referencia simulada**
/// (decisión del responsable del 01/10/2026, DT-016).
///
/// Un locus **compuesto** sobre referencia **simulada** tiene que decir cuántas
/// unidades lleva cada bloque en la referencia: sin eso no hay forma de
/// construir su tracto, y el tramo sería relleno al azar. Por eso el campo es
/// obligatorio ahí, y por eso se rechaza donde no haría nada —en un locus de un
/// solo motivo, cuyo tramo ya fija el tracto, y con referencia suministrada,
/// cuyo tracto se lee del FASTA—: un campo que no hace nada y se acepta en
/// silencio hace creer al usuario que algo cambió.
///
/// Las unidades son **enteras** y de al menos una por bloque, como el número
/// de unidades por bloque de un alelo compuesto, que se materializa con
/// [`TandemRepeatTract::pure_blocks`]. La suma de unidades por longitud de cada
/// motivo tiene que ser **exactamente** la longitud del tramo: un tracto más
/// corto dejaría bases al azar dentro del locus y uno más largo invadiría el
/// flanco.
fn resolve_reference_units(
    general: &DeclaredGeneralConfig,
    locus: &crate::configuration::DeclaredRepeatLocus,
    declaration: &TandemRepeatLocus,
) -> Result<Option<Vec<u64>>, Diagnostic> {
    const FIELD: &str = "reference.repeat_loci.reference_units";
    let simulated = matches!(
        general.reference().source,
        crate::configuration::DeclaredReferenceSource::Simulated { .. }
    );
    let motifs = declaration.motifs();
    // **Al insertar, el campo es obligatorio también con un solo motivo**
    // (DTC-098): no hay tramo en el genoma de partida ni tracto que leer del
    // FASTA, y las unidades son lo único que dice qué se inserta. El plan de
    // inserción ya exigió el campo, un valor por bloque y ninguno a cero.
    if locus.insert_after.is_some() {
        return Ok(locus
            .reference_units
            .as_ref()
            .map(|declared| declared.iter().map(|value| *value.value()).collect()));
    }
    let Some(declared) = &locus.reference_units else {
        if simulated && motifs.len() > 1 {
            return Err(repeat_error(
                general,
                FIELD,
                locus.id.value(),
                format!(
                    "el locus es compuesto ({} motivos) y la referencia es simulada: declare reference_units con las unidades de cada bloque en la referencia, en el orden de los motivos",
                    motifs.len()
                ),
            ));
        }
        return Ok(None);
    };
    if motifs.len() == 1 {
        return Err(repeat_error(
            general,
            FIELD,
            locus.id.value(),
            "reference_units no se declara en un locus de un solo motivo con tramo: el tramo declarado ya fija el tracto; solo un locus compuesto, o uno que inserta su tracto con insert_after, lo declara",
        ));
    }
    if !simulated {
        return Err(repeat_error(
            general,
            FIELD,
            locus.id.value(),
            "reference_units solo se declara con referencia simulada: con referencia suministrada el tracto se lee del FASTA",
        ));
    }
    if declared.len() != motifs.len() {
        return Err(repeat_error(
            general,
            FIELD,
            locus.id.value(),
            format!(
                "reference_units declara {} valores y el locus tiene {} motivos: hace falta uno por bloque",
                declared.len(),
                motifs.len()
            ),
        ));
    }
    let units = declared
        .iter()
        .map(|value| *value.value())
        .collect::<Vec<_>>();
    if let Some(block) = units.iter().position(|value| *value == 0) {
        return Err(repeat_error(
            general,
            FIELD,
            locus.id.value(),
            format!(
                "reference_units declara cero unidades en el bloque {block}: cada bloque lleva al menos una"
            ),
        ));
    }
    // En u128 no hay desbordamiento posible: cada término cabe en 128 bits y
    // hay un término por motivo.
    let sum: u128 = units
        .iter()
        .zip(motifs)
        .map(|(units, motif)| u128::from(*units) * u128::from(motif.length()))
        .sum();
    let interval = declaration.interval();
    if sum != u128::from(interval.len()) {
        return Err(repeat_error(
            general,
            FIELD,
            locus.id.value(),
            format!(
                "reference_units suma {sum} bases (unidades por longitud de cada motivo) y el tramo declarado {}..{} mide {}",
                locus.start.as_ref().map_or(0, |value| *value.value()),
                locus.end.as_ref().map_or(0, |value| *value.value()),
                interval.len()
            ),
        ));
    }
    Ok(Some(units))
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
/// `scope` es el prefijo del campo que se denuncia —`reference.repeat_admission`
/// para los generales, `reference.repeat_loci` para los de un locus— y `subject`
/// identifica a quién pertenecen. Los recibe **aquí** en vez de envolver el
/// diagnóstico fuera porque `Diagnostic::with` **sobrescribe** la clave: un
/// envoltorio que añadiera `fields` pisaría el campo que el error interior ya
/// había nombrado, y mandaría al usuario a mirar una clave que no ha escrito.
fn resolve_admission(
    general: &DeclaredGeneralConfig,
    declared: Option<&crate::configuration::DeclaredRepeatAdmission>,
    fallback: EffectiveRepeatAdmission,
    scope: &str,
    subject: Option<&str>,
) -> Result<EffectiveRepeatAdmission, Diagnostic> {
    let describe = |field: &str, detail: &str| {
        let detail = subject.map_or_else(
            || detail.to_owned(),
            |subject| format!("locus {subject}: {detail}"),
        );
        config_error(general, &format!("{scope}.{field}"), detail)
    };
    let Some(declared) = declared else {
        return Ok(fallback);
    };
    let min_purity = match &declared.min_purity {
        None => fallback.min_purity(),
        Some(value) => {
            let parsed = value
                .literal()
                .parse::<f64>()
                .map_err(|_| describe("min_purity", "debe ser un decimal válido"))?;
            if !parsed.is_finite() || !(0.0..=1.0).contains(&parsed) {
                return Err(describe(
                    "min_purity",
                    "la pureza mínima debe ser finita y pertenecer a [0,1]",
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
        return Err(describe(
            "min_longest_run_units",
            "la tirada mínima debe ser un entero positivo",
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
        "reference.repeat_admission",
        None,
    )
}

/// La cota de bases de tracto a partir de los loci ya resueltos.
///
/// Vive aparte porque **se calcula dos veces**: al resolver la configuración,
/// con los alelos declarados, y otra vez cuando el generador de deslizamiento
/// crea los suyos. Dos copias de esta aritmética acabarían diciendo cosas
/// distintas del mismo escenario.
fn recompute_repeat_output_budget(
    loci: &[EffectiveRepeatLocus],
    copies: u64,
    declared_cap: Option<u64>,
) -> EffectiveRepeatOutputBudget {
    let mut per_locus = Vec::with_capacity(loci.len());
    for locus in loci {
        let longest = locus
            .alleles()
            .iter()
            .map(|allele| allele.tract().length())
            .max()
            .unwrap_or(0);
        per_locus.push(RepeatLocusOutputBound {
            locus_id: locus.id().to_owned(),
            truth_bases: EffectiveRepeatOutputBudget::TRUTH_COPIES_PER_HAPLOTYPE_COPY
                .saturating_mul(copies)
                .saturating_mul(longest),
        });
    }
    EffectiveRepeatOutputBudget {
        per_locus,
        declared_cap,
    }
}

/// El **alelo ancestral** de un locus: las bases que el genoma de referencia
/// lleva hoy en su intervalo, ya estructuradas como tracto de su motivo.
///
/// Es de donde parte cada linaje del generador. Que la referencia lleve un
/// alelo que después ningún individuo tenga es una de las tres limitaciones
/// declaradas del modelo, no un defecto.
fn reference_repeat_tract(
    reference: &Reference,
    locus: &EffectiveRepeatLocus,
) -> Result<TandemRepeatTract, Diagnostic> {
    let interval = locus.declaration().interval();
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == interval.sequence_id())
        .ok_or_else(|| {
            slippage_reference_error(
                locus.id(),
                format!("la referencia no contiene {}", interval.sequence_id()),
            )
        })?;
    let bases = sequence.sequence().as_bytes().map_err(|error| {
        slippage_reference_error(
            locus.id(),
            format!("las bases de la referencia no están materializadas: {error:?}"),
        )
    })?;
    let start = usize::try_from(interval.start())
        .map_err(|_| slippage_reference_error(locus.id(), "el inicio del intervalo desborda"))?;
    let end = usize::try_from(interval.end())
        .map_err(|_| slippage_reference_error(locus.id(), "el final del intervalo desborda"))?;
    let slice = bases
        .get(start..end)
        .ok_or_else(|| slippage_reference_error(locus.id(), "el intervalo rebasa la secuencia"))?;
    let motif = locus.single_motif().ok_or_else(|| {
        slippage_reference_error(locus.id(), "el generador solo admite un locus de un motivo")
    })?;
    TandemRepeatTract::new(motif.clone(), slice)
        .map_err(|error| slippage_reference_error(locus.id(), error.detail()))
}

fn slippage_reference_error(locus_id: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("fields", "reference.repeat_loci.slippage")
    .with("detail", format!("locus {locus_id}: {}", detail.into()))
}

/// Lee un decimal del bloque del generador sin perder el literal escrito.
fn parse_slippage_number(
    general: &DeclaredGeneralConfig,
    locus_id: &str,
    field: &str,
    declared: &crate::configuration::DeclaredNumber,
) -> Result<f64, Diagnostic> {
    declared.literal().parse::<f64>().map_err(|_| {
        repeat_error(
            general,
            "reference.repeat_loci.slippage",
            locus_id,
            format!("{field} no es un decimal válido: {}", declared.literal()),
        )
    })
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
