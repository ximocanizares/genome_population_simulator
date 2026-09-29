//! Contrato TOML mínimo de F1 para referencias iniciales y reducidas.

use std::fmt;
use std::ops::Range;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::diagnostics::Diagnostic;
pub use crate::reference::SequenceClass;

/// El recorrido vigente solo admite el documento general.
/// Los escenarios se incorporarán con su primer consumidor biológico real.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ConfigKind {
    General,
    F5Scenario,
}

impl ConfigKind {
    pub(crate) fn parse(value: &str) -> Option<Self> {
        match value {
            "general" => Some(Self::General),
            "f5_scenario" => Some(Self::F5Scenario),
            _ => None,
        }
    }

    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::General => "general",
            Self::F5Scenario => "f5_scenario",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ExpectedConfigKind {
    General,
    /// Frontera de CLI: conserva general v2 y despacha general v3 al lector
    /// F5 antes de interpretar el cuerpo del documento.
    GeneralDispatch,
    /// Entrada general exclusiva de F5. La versión 3 no es una extensión de
    /// general v2 y por ello se admite por una frontera distinta.
    F5General,
    F5Scenario,
}

impl ExpectedConfigKind {
    pub(crate) const fn accepts(self, actual: ConfigKind) -> bool {
        matches!(
            (self, actual),
            (
                Self::General | Self::GeneralDispatch | Self::F5General,
                ConfigKind::General
            ) | (Self::F5Scenario, ConfigKind::F5Scenario)
        )
    }

    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::General | Self::GeneralDispatch | Self::F5General => "general",
            Self::F5Scenario => "f5_scenario",
        }
    }
}

#[derive(Clone, Debug)]
pub struct AdmittedToml {
    path: PathBuf,
    bytes: Arc<[u8]>,
    kind: ConfigKind,
}

impl AdmittedToml {
    pub(crate) fn new(path: PathBuf, bytes: Arc<[u8]>, kind: ConfigKind) -> Self {
        Self { path, bytes, kind }
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }

    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    #[must_use]
    pub const fn kind(&self) -> ConfigKind {
        self.kind
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SourceLocation {
    pub byte_start: usize,
    pub byte_end: usize,
    pub line: usize,
    pub column: usize,
}

impl SourceLocation {
    #[must_use]
    pub fn from_span(bytes: &[u8], span: Range<usize>) -> Self {
        let mut location = Self::from_offset(bytes, span.start);
        location.byte_end = span.end;
        location
    }

    #[must_use]
    pub fn from_offset(bytes: &[u8], offset: usize) -> Self {
        let safe_offset = offset.min(bytes.len());
        let line_start = bytes[..safe_offset]
            .iter()
            .rposition(|byte| *byte == b'\n')
            .map_or(0, |index| index + 1);
        Self {
            byte_start: safe_offset,
            byte_end: safe_offset,
            line: bytes[..safe_offset]
                .iter()
                .filter(|byte| **byte == b'\n')
                .count()
                + 1,
            column: String::from_utf8_lossy(&bytes[line_start..safe_offset])
                .chars()
                .count()
                + 1,
        }
    }
}

impl fmt::Display for SourceLocation {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            formatter,
            "line {}, column {}, bytes {}..{}",
            self.line, self.column, self.byte_start, self.byte_end
        )
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Located<T> {
    value: T,
    location: SourceLocation,
}

impl<T> Located<T> {
    pub(crate) const fn new(value: T, location: SourceLocation) -> Self {
        Self { value, location }
    }

    #[must_use]
    pub const fn value(&self) -> &T {
        &self.value
    }

    #[must_use]
    pub const fn location(&self) -> &SourceLocation {
        &self.location
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredNumber {
    literal: String,
    location: SourceLocation,
}

impl DeclaredNumber {
    pub(crate) const fn new(literal: String, location: SourceLocation) -> Self {
        Self { literal, location }
    }

    #[must_use]
    pub fn literal(&self) -> &str {
        &self.literal
    }

    #[must_use]
    pub const fn location(&self) -> &SourceLocation {
        &self.location
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredCompositionValue {
    pub base: Located<String>,
    pub value: DeclaredNumber,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredSequence {
    pub id: Located<String>,
    pub class: Located<SequenceClass>,
    pub length: Located<u64>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredReferenceSource {
    Simulated {
        composition: Vec<DeclaredCompositionValue>,
        unknown_base_fraction: Option<DeclaredNumber>,
    },
    Supplied {
        fasta: Located<String>,
    },
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredReference {
    pub id: Located<String>,
    pub source: DeclaredReferenceSource,
    pub sequences: Vec<DeclaredSequence>,
    /// Loci de repetición en tándem declarados sobre esta referencia. Solo el
    /// documento general v4 puede traerlos; en v2 la clave es desconocida.
    pub repeat_loci: Vec<DeclaredRepeatLocus>,
    /// Mínimos de la guarda de admisión con referencia real (F6.4.2). Se
    /// declaran una sola vez para toda la referencia; si faltan, mandan los
    /// valores por omisión, que la configuración efectiva publica igualmente.
    pub repeat_admission: Option<DeclaredRepeatAdmission>,
    /// Tope declarable de bases de tracto publicadas, en la verdad y en el
    /// VCF juntas. Solo el documento general v4 puede traerlo. Quien no lo
    /// declara no nota nada: sin tope no hay rechazo.
    pub max_tract_output_bases: Option<Located<u64>>,
}

/// Un locus de repetición tal y como lo escribe el usuario.
///
/// El intervalo se declara **entero, inicio y final**, y no solo el inicio: la
/// pureza y la tirada más larga se miden **dentro de un intervalo**, y sin
/// final no hay dónde medirlas. La alternativa —extender el tracto desde el
/// inicio mientras el motivo siga apareciendo— es un algoritmo que nadie ha
/// fijado y queda fuera de este trabajo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredRepeatLocus {
    pub id: Located<String>,
    pub sequence: Located<String>,
    /// Base uno, inclusivo.
    pub start: Located<u64>,
    /// Base uno, inclusivo.
    pub end: Located<u64>,
    /// Los motivos del locus, **en orden**. Un locus corriente declara uno con
    /// la clave `motif`; un microsatélite compuesto declara varios con `motifs`,
    /// y entonces cada alelo da su número de unidades por bloque.
    pub motifs: Vec<Located<String>>,
    /// Tope propio del tamaño del tracto, **en unidades** (F6.4.3.1).
    pub max_units: Option<Located<u64>>,
    /// Mínimos de admisión **propios de este locus**, que sustituyen a los
    /// generales. El mismo número no significa lo mismo según el motivo: tres
    /// unidades de un mononucleótido son tres bases, que salen por azar a cada
    /// paso; tres de un tetranucleótido son doce, que no.
    pub admission: Option<DeclaredRepeatAdmission>,
    /// Generador de alelos por **deslizamiento biológico** (F6.4.11).
    ///
    /// Es **excluyente** con `alleles`: o el usuario escribe los alelos y sus
    /// frecuencias, o los genera el simulador partiendo del alelo que el genoma
    /// de referencia lleva en este locus. Declarar los dos sería ambiguo sobre
    /// el origen de cada alelo, que es justo lo que el manifiesto publica.
    pub slippage: Option<DeclaredRepeatSlippage>,
    pub alleles: Vec<DeclaredRepeatAllele>,
}

/// Parámetros del generador de alelos por deslizamiento, tal y como se
/// escriben.
///
/// No se validan aquí: la frontera de configuración solo los lee. Quien
/// rechaza un modelo distinto del escalonado o un salto de varias unidades es
/// `repeat_slippage::SlippageSettings::new`, que nombra el parámetro y el
/// valor.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredRepeatSlippage {
    /// El modelo. Solo se admite el escalonado.
    pub model: Located<String>,
    /// Media de la Poisson que decide cuántos pasos recibe cada linaje.
    pub mean_events: DeclaredNumber,
    /// Probabilidad de que un paso sea una expansión. Si falta, es simétrico.
    pub expansion_bias: Option<DeclaredNumber>,
    /// Unidades por paso. Si falta, una. Cualquier otro valor se rechaza.
    pub step_units: Option<Located<u64>>,
}

/// Un alelo declarado del locus, con su frecuencia.
///
/// Se declara **o** por número de unidades más sus interrupciones, **o** por su
/// secuencia exacta, que es lo que permite una referencia simulada. Nunca se
/// declara una pureza: la pureza es siempre una salida calculada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredRepeatAllele {
    /// Número de unidades **por bloque**, en el orden de los motivos del locus.
    /// Un locus de un solo motivo lleva aquí un único valor, que es lo que el
    /// usuario escribió como escalar.
    pub units: Vec<Located<u64>>,
    /// Secuencia exacta del tracto. Solo se admite con **un** motivo: unas bases
    /// sueltas no pueden decir dónde acaba un bloque y empieza el siguiente, y
    /// esa es justamente la información que la estructura en bloques aporta.
    pub sequence: Option<Located<String>>,
    pub frequency: DeclaredNumber,
    pub interruptions: Vec<DeclaredRepeatInterruption>,
}

/// Una interrupción declarada dentro de un alelo definido por unidades.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredRepeatInterruption {
    /// Desplazamiento en base cero desde el comienzo del tracto.
    pub offset: Located<u64>,
    /// Bases que sustituyen a las del motivo a partir de ese desplazamiento.
    pub bases: Located<String>,
}

/// Mínimos configurables de la guarda de admisión con referencia real.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredRepeatAdmission {
    pub min_purity: Option<DeclaredNumber>,
    pub min_longest_run_units: Option<Located<u64>>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredReductionTarget {
    Bases(Located<u64>),
    ReductionPercent(DeclaredNumber),
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredReducedReference {
    pub id: Located<String>,
    pub source_reference_id: Located<String>,
    pub region_size_bases: Located<u64>,
    pub target: DeclaredReductionTarget,
}

/// Solicitud explícita de la vista VCF de verdad. Ausente significa que la
/// ejecución conserva únicamente `truth.jsonl`, conforme a REQ-SAL-06.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredVcfOutput {
    pub mode: Located<String>,
}

/// Único evento técnico que F3 puede consumir. La tabla conserva únicamente
/// los parámetros de autoría: el perfil mínimo se resuelve después, sin abrir
/// aún el catálogo de perfiles, pares o referencias de mapeo de F4.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredTechnicalEvent {
    pub event_id: Located<String>,
    pub individual_ids: Vec<Located<String>>,
    pub target: DeclaredTechnicalEventTarget,
    /// Modo explícito de sustituciones de secuenciación F4. Ausente conserva
    /// el modelo probabilístico incorporado.
    pub sequencing_errors: Option<Located<String>>,
    pub read_length_bases: Option<Located<u64>>,
    /// Fracción Bernoulli por pareja limpia; el valor por omisión se resuelve
    /// en la configuración efectiva del perfil F4.
    pub chimera_rate: Option<DeclaredNumber>,
    /// Fracción Bernoulli por pareja posterior a contaminación biológica. La
    /// reasignación se materializa exclusivamente en F4.5.
    pub index_hopping_rate: Option<DeclaredNumber>,
    /// Límites planos de la mezcla biológica previa al demultiplexado. Si se
    /// omiten ambos, el perfil inicial se mantiene neutro (0--0 %).
    pub contamination_min_percent: Option<DeclaredNumber>,
    pub contamination_max_percent: Option<DeclaredNumber>,
    pub contamination_recipient_ids: Option<Vec<Located<String>>>,
    /// Tartamudeo de base de la polimerasa sobre los microsatélites
    /// declarados. Ausente significa **desactivado**, que es lo que hace que
    /// una ejecución anterior a este trabajo publique los mismos bytes.
    pub stutter: Option<DeclaredStutter>,
    /// Preparación de biblioteca declarada. Ausente significa `pcr_free`, que
    /// es la única que existía y la que conserva los bytes publicados.
    pub library_preparation: Option<Located<String>>,
    /// Tabla de la amplificación. Solo se admite con preparación `pcr`, y con
    /// ella es **obligatoria**.
    pub pcr: Option<DeclaredPcr>,
    pub libraries: Vec<DeclaredTechnicalLibrary>,
}

/// Amplificación declarada por el autor del escenario.
///
/// Los tres primeros parámetros son **obligatorios y sin valor por omisión**:
/// la fracción de duplicados depende de la masa de ADN de entrada, del número
/// de ciclos y del laboratorio, de modo que ningún valor por omisión sería
/// defendible. La mejor fuente para ese número son los informes de calidad del
/// propio laboratorio, no una media de genomas humanos ajenos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredPcr {
    pub duplicate_fraction: DeclaredNumber,
    pub cycles: Located<u64>,
    pub efficiency: DeclaredNumber,
    pub early_cycles: Option<Located<u64>>,
    pub max_family_size: Option<Located<u64>>,
    pub gc_max_attempts: Option<Located<u64>>,
    pub gc_relative_efficiency: Vec<DeclaredGcCurveKnot>,
    pub efficiency_classes: Vec<DeclaredEfficiencyClass>,
    /// Tartamudeo clonal heredable. Ausente significa que las probabilidades
    /// por ciclo se **derivan** de las dos columnas de la calibración y del
    /// número de ciclos, que es lo normal; la tabla solo existe para
    /// sustituirlas.
    pub stutter: Option<DeclaredPcrStutter>,
}

/// Sustituciones directas de la probabilidad de resbalón **por ciclo**.
///
/// Solo tiene efecto con loci de repetición declarados. Lo que se declara es
/// la probabilidad de un **ciclo**, no la fracción de extremo a extremo: son
/// dos magnitudes distintas y confundirlas cambiaría el resultado en un orden
/// de magnitud.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredPcrStutter {
    pub per_cycle_slip: Vec<DeclaredPerCycleSlip>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredPerCycleSlip {
    pub motif_length_bases: Located<u64>,
    pub probability: DeclaredNumber,
}

/// Un nudo de la curva de composición: fracción de GC y eficiencia relativa.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredGcCurveKnot {
    pub gc_fraction: DeclaredNumber,
    pub relative_efficiency: DeclaredNumber,
}

/// Una clase de eficiencia declarada, con su peso relativo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredEfficiencyClass {
    pub efficiency: DeclaredNumber,
    pub weight: DeclaredNumber,
}

/// Tartamudeo declarado por el autor del escenario.
///
/// La **estructura** del paso no se declara: la fija el algoritmo. Lo que se
/// declara son los **pesos relativos** dentro de esa estructura, que son
/// obligatorios cuando el tartamudeo está activo y no tienen valor por
/// omisión, y opcionalmente las tasas por longitud de motivo y el coeficiente
/// de la tirada.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredStutter {
    pub enabled: Option<Located<bool>>,
    /// Tasas por longitud de motivo, en pares de base, declaradas una a una.
    pub motif_length_rates: Vec<DeclaredStutterMotifRate>,
    pub run_length_coefficient: Option<DeclaredNumber>,
    pub steps: Vec<DeclaredStutterStep>,
    pub locus_rates: Vec<DeclaredStutterLocusRate>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredStutterMotifRate {
    pub motif_length_bases: Located<u64>,
    pub rate: DeclaredNumber,
}

/// Un paso declarado, **en unidades enteras**. Un paso expresado en bases no
/// existe en este contrato y su clave se rechaza.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredStutterStep {
    pub units: Located<i64>,
    pub weight: DeclaredNumber,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredStutterLocusRate {
    pub locus_id: Located<String>,
    pub rate: DeclaredNumber,
}

/// UDI de la única biblioteca de un individuo dentro de la tanda F4.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeclaredTechnicalLibrary {
    pub individual_id: Located<String>,
    pub i7: Located<String>,
    pub i5: Located<String>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DeclaredTechnicalEventTarget {
    CoverageWgs(DeclaredNumber),
    ReadCount(Located<u64>),
}

/// Selección declarada del recurso común de frecuencias SNV.
///
/// Solo describe la frontera de autoría. La preparación, la reutilización
/// verificable y los parámetros de generación tienen sus propios pasos F2.2.
#[derive(Clone, Debug, PartialEq)]
pub enum DeclaredPopulationSnvs {
    Supplied {
        input_path: Located<String>,
    },
    Reuse {
        metadata_path: Located<String>,
    },
    Generate {
        parameters: Box<DeclaredPopulationSnvsGeneration>,
    },
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsGeneration {
    pub count_mode: Located<String>,
    pub locus_count: Option<Located<u64>>,
    pub mean_loci: Option<DeclaredNumber>,
    pub mean_density_per_base: Option<DeclaredNumber>,
    pub spatial_mode: Located<String>,
    pub substitution_profile: Located<String>,
    pub outside_weight: Option<DeclaredNumber>,
    pub regions: Vec<DeclaredPopulationSnvsRegion>,
    pub populations: Vec<DeclaredPopulationSnvsPopulation>,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsRegion {
    pub contig: Located<String>,
    pub start: Located<u64>,
    pub end: Located<u64>,
    pub weight: DeclaredNumber,
}
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredPopulationSnvsPopulation {
    pub population_id: Located<String>,
    pub frequency_model: Located<String>,
    pub frequency_mode: Option<Located<String>>,
    pub min_alt_af: Option<DeclaredNumber>,
    pub max_alt_af: Option<DeclaredNumber>,
    pub spectrum: Vec<DeclaredFrequencySpectrumComponent>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredFrequencySpectrumComponent {
    pub kind: Located<String>,
    pub fixed_alt_af: Option<DeclaredNumber>,
    pub min_alt_af: Option<DeclaredNumber>,
    pub max_alt_af: Option<DeclaredNumber>,
    pub weight: DeclaredNumber,
    pub sampling: Option<Located<String>>,
}

/// Una cohorte de creación directa procedente de una población de la matriz.
///
/// No representa una población evolutiva: únicamente agrupa la procedencia y
/// el modelo genotípico que consumirán los pasos posteriores de F2.3.
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredInitialPopulation {
    pub population_id: Located<String>,
    pub individual_count: Located<u64>,
    pub ploidy: Located<u64>,
    pub genotype_model: Located<String>,
    pub inbreeding_coefficient: Option<DeclaredNumber>,
    pub individual_ids: Option<Vec<Located<String>>>,
    pub founder_keys: Option<Vec<Located<String>>>,
}

/// Cambio directo mínimo de F2.3.4 dirigido a una copia individual.
#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredDirectEvent {
    pub event_id: Located<String>,
    pub individual_id: Located<String>,
    pub phase: Located<u64>,
    pub kind: Located<String>,
    pub contig: Located<String>,
    pub position_1based: Option<Located<u64>>,
    pub reference: Option<Located<String>>,
    pub alternate: Option<Located<String>>,
    pub point: Option<Located<u64>>,
    pub bases: Option<Located<String>>,
    pub start: Option<Located<u64>>,
    pub end: Option<Located<u64>>,
}

#[derive(Clone, Debug)]
pub struct DeclaredGeneralConfig {
    source: AdmittedToml,
    /// Versión declarada del documento general: `"2"` o `"4"`. Decide qué
    /// versión de configuración efectiva se publica y si las repeticiones son
    /// un campo del contrato o un campo desconocido.
    document_version: String,
    seed: Option<Located<String>>,
    ploidy: Option<Located<u64>>,
    normalization_tolerance: Option<DeclaredNumber>,
    reference: DeclaredReference,
    reduced_references: Vec<DeclaredReducedReference>,
    population_snvs: Option<DeclaredPopulationSnvs>,
    vcf_output: Option<DeclaredVcfOutput>,
    initial_populations: Vec<DeclaredInitialPopulation>,
    direct_events: Vec<DeclaredDirectEvent>,
    technical_event: Option<DeclaredTechnicalEvent>,
}

impl DeclaredGeneralConfig {
    pub(crate) fn new(
        source: AdmittedToml,
        seed: Option<Located<String>>,
        ploidy: Option<Located<u64>>,
        normalization_tolerance: Option<DeclaredNumber>,
        reference: DeclaredReference,
        reduced_references: Vec<DeclaredReducedReference>,
        population_snvs: Option<DeclaredPopulationSnvs>,
    ) -> Self {
        Self {
            source,
            // El lector la fija con `with_document_version` en cuanto la ha
            // leído; la versión 2 es la histórica y sirve de valor de partida.
            document_version: "2".to_owned(),
            seed,
            ploidy,
            normalization_tolerance,
            reference,
            reduced_references,
            population_snvs,
            vcf_output: None,
            initial_populations: Vec::new(),
            direct_events: Vec::new(),
            technical_event: None,
        }
    }

    pub(crate) fn with_document_version(mut self, document_version: String) -> Self {
        self.document_version = document_version;
        self
    }

    pub(crate) fn with_initial_populations(
        mut self,
        initial_populations: Vec<DeclaredInitialPopulation>,
    ) -> Self {
        self.initial_populations = initial_populations;
        self
    }

    pub(crate) fn with_direct_events(mut self, direct_events: Vec<DeclaredDirectEvent>) -> Self {
        self.direct_events = direct_events;
        self
    }

    pub(crate) fn with_vcf_output(mut self, vcf_output: Option<DeclaredVcfOutput>) -> Self {
        self.vcf_output = vcf_output;
        self
    }

    pub(crate) fn with_technical_event(
        mut self,
        technical_event: Option<DeclaredTechnicalEvent>,
    ) -> Self {
        self.technical_event = technical_event;
        self
    }

    #[must_use]
    pub const fn source(&self) -> &AdmittedToml {
        &self.source
    }

    #[must_use]
    pub const fn seed(&self) -> Option<&Located<String>> {
        self.seed.as_ref()
    }

    /// Ploidía global solicitada; ausente conserva el valor diploide compatible con F1.
    #[must_use]
    pub const fn ploidy(&self) -> Option<&Located<u64>> {
        self.ploidy.as_ref()
    }

    /// Tolerancia general DTC-049 solicitada para las frecuencias de F2.
    #[must_use]
    pub const fn normalization_tolerance(&self) -> Option<&DeclaredNumber> {
        self.normalization_tolerance.as_ref()
    }

    /// Versión declarada del documento general.
    #[must_use]
    pub fn document_version(&self) -> &str {
        &self.document_version
    }

    #[must_use]
    pub const fn reference(&self) -> &DeclaredReference {
        &self.reference
    }

    #[must_use]
    pub fn reduced_references(&self) -> &[DeclaredReducedReference] {
        &self.reduced_references
    }

    /// Fuente de matriz SNV solicitada, si la configuración la declara.
    #[must_use]
    pub const fn population_snvs(&self) -> Option<&DeclaredPopulationSnvs> {
        self.population_snvs.as_ref()
    }

    #[must_use]
    pub const fn vcf_output(&self) -> Option<&DeclaredVcfOutput> {
        self.vcf_output.as_ref()
    }

    /// Cohortes directas declaradas; no son estructura poblacional evolutiva.
    #[must_use]
    pub fn initial_populations(&self) -> &[DeclaredInitialPopulation] {
        &self.initial_populations
    }

    #[must_use]
    pub fn direct_events(&self) -> &[DeclaredDirectEvent] {
        &self.direct_events
    }

    #[must_use]
    pub const fn technical_event(&self) -> Option<&DeclaredTechnicalEvent> {
        self.technical_event.as_ref()
    }
}

pub(crate) fn semantic_error(
    source: &AdmittedToml,
    fields: impl Into<String>,
    detail: impl Into<String>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", source.path().display().to_string())
    .with("fields", fields.into())
    .with("detail", detail.into())
}

/// Contrato de autoría F5. Se mantiene separado de `DeclaredGeneralConfig`
/// porque v3 reemplaza —no amplía— el documento general v2.
#[derive(Clone, Debug)]
pub struct DeclaredF5GeneralConfig {
    pub source: AdmittedToml,
    pub seed: Option<Located<String>>,
    pub ploidy: Located<u64>,
    pub reference_id: Located<String>,
    pub reference_path: Located<String>,
    pub reference_sha256: Located<String>,
    pub scenario_id: Located<String>,
    pub scenario_path: Located<String>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum DeclaredF5Target {
    BlockCount(Located<u64>),
    CoveredBases(Located<u64>),
    CoverageFraction(DeclaredNumber),
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredF5Subpopulation {
    pub population_id: Located<String>,
    pub individual_count: Located<u64>,
    pub sampling: Option<Located<String>>,
    pub inbreeding_f: Option<DeclaredNumber>,
    pub founder_ids: Option<Vec<Located<String>>>,
    pub founder_keys: Option<Vec<Located<String>>>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct DeclaredF5SamplingOverride {
    pub population_id: Located<String>,
    pub block_index: Located<u64>,
    pub mode: Located<String>,
    pub controlled_max_abs_error: Option<DeclaredNumber>,
}

#[derive(Clone, Debug)]
pub struct DeclaredF5Scenario {
    pub source: AdmittedToml,
    pub scenario_id: Located<String>,
    pub founder_count: Located<u64>,
    pub profile_metadata_path: Located<String>,
    pub profile_metadata_sha256: Located<String>,
    pub excluded_regions_path: Option<Located<String>>,
    pub excluded_regions_sha256: Option<Located<String>>,
    pub inbreeding_f: Option<DeclaredNumber>,
    pub max_layout_attempts: Option<Located<u64>>,
    pub max_haplotype_attempts: Option<Located<u64>>,
    pub max_swap_proposals: Option<Located<u64>>,
    pub max_haplotype_cells: Option<Located<u64>>,
    pub max_operation_units: Option<Located<u64>>,
    pub max_working_memory_bytes: Option<Located<u64>>,
    pub max_output_bytes: Option<Located<u64>>,
    pub target: DeclaredF5Target,
    pub subpopulations: Vec<DeclaredF5Subpopulation>,
    pub sampling_overrides: Vec<DeclaredF5SamplingOverride>,
}
