//! Contrato físico mínimo VCF 4.5 para F2.5.1.
//!
//! Este módulo escribe la cabecera, los DTO de registros ya preparados y el
//! envoltorio físico de salida. La proyección de alelos, `GT` y campos `INFO`
//! pertenece a los pasos F2.5.2--F2.5.4; en particular, no consulta
//! `truth.jsonl` ni valida su semántica.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use noodles_csi::{
    self as csi,
    binning_index::{Indexer, index::reference_sequence::index::BinnedIndex},
};
use noodles_tabix as tabix;
use noodles_vcf::variant::Record as _;
use noodles_vcf::{
    self as vcf,
    header::{
        FileFormat,
        record::value::{Map, map::Contig},
    },
};

use crate::identity::{ExecutionId, ReferenceId, SequenceId};
use crate::vcf_grouping::{GroupedVcfRecord, VcfRepeatSequence};

pub const VCF_FILE_FORMAT: &str = "VCFv4.5";
pub const TBI_MAXIMUM_POSITION_EXCLUSIVE: u64 = 1 << 29;
pub const VCF_REFERENCE_META_KEY: &str = "reference";
pub const VCF_EXECUTION_META_KEY: &str = "gps_execution_id";
pub const VCF_QUAL_META_KEY: &str = "gps_qual";
pub const VCF_QUAL_SEMANTICS: &str = ".";
pub const VCF_PASS_FILTER_DESCRIPTION: &str = "All filters passed";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VcfOutputMode {
    Text,
    Bgzf,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum VcfIndexKind {
    Tbi,
    Csi,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfContig {
    name: SequenceId,
    length: u64,
}

impl VcfContig {
    pub fn new(name: impl Into<String>, length: u64) -> Result<Self, VcfContractError> {
        let name = name.into();
        let name = SequenceId::parse(name.clone(), "VCF header contig")
            .map_err(|_| VcfContractError::InvalidContigIdentifier { value: name })?;
        if length == 0 {
            return Err(VcfContractError::ZeroLengthContig {
                name: name.to_string(),
            });
        }
        Ok(Self { name, length })
    }

    #[must_use]
    pub fn name(&self) -> &str {
        self.name.as_str()
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfHeaderContract {
    identity: VcfHeaderIdentity,
    contigs: Vec<VcfContig>,
    samples: Vec<String>,
}

/// Identidad de la única vista VCF directa de F2.
///
/// F2 aún no tiene conjuntos de escenario: el conjunto es exactamente la
/// ejecución que lo publica. Por ello el ID de ejecución ya existente, junto
/// al ID tipado de la referencia efectiva, identifica sin ambigüedad esta
/// vista sin introducir otra identidad de configuración.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfHeaderIdentity {
    reference_id: ReferenceId,
    execution_id: ExecutionId,
}

impl VcfHeaderIdentity {
    #[must_use]
    pub const fn new(reference_id: ReferenceId, execution_id: ExecutionId) -> Self {
        Self {
            reference_id,
            execution_id,
        }
    }

    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn execution_id(&self) -> &ExecutionId {
        &self.execution_id
    }
}

impl VcfHeaderContract {
    pub fn new(
        identity: VcfHeaderIdentity,
        contigs: Vec<VcfContig>,
        samples: Vec<String>,
    ) -> Result<Self, VcfContractError> {
        if contigs.is_empty() {
            return Err(VcfContractError::NoContigs);
        }
        let mut names = BTreeSet::<String>::new();
        for contig in &contigs {
            if !names.insert(contig.name.to_string()) {
                return Err(VcfContractError::DuplicateContig(contig.name.to_string()));
            }
        }
        for sample in &samples {
            validate_identifier("sample", sample)?;
            if !names.insert(format!("sample:{sample}")) {
                return Err(VcfContractError::DuplicateSample(sample.clone()));
            }
        }
        Ok(Self {
            identity,
            contigs,
            samples,
        })
    }

    #[must_use]
    pub const fn identity(&self) -> &VcfHeaderIdentity {
        &self.identity
    }

    #[must_use]
    pub fn contigs(&self) -> &[VcfContig] {
        &self.contigs
    }

    #[must_use]
    pub fn samples(&self) -> &[String] {
        &self.samples
    }

    #[must_use]
    pub fn selected_index(&self) -> VcfIndexKind {
        if self
            .contigs
            .iter()
            .all(|contig| contig.length < TBI_MAXIMUM_POSITION_EXCLUSIVE)
        {
            VcfIndexKind::Tbi
        } else {
            VcfIndexKind::Csi
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfOutputPaths {
    pub vcf: PathBuf,
    pub index: Option<PathBuf>,
    pub index_kind: Option<VcfIndexKind>,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfContractError {
    InvalidContigIdentifier { value: String },
    EmptyIdentifier { kind: &'static str },
    InvalidIdentifier { kind: &'static str, value: String },
    ZeroLengthContig { name: String },
    NoContigs,
    DuplicateContig(String),
    DuplicateSample(String),
    LengthTooLarge { name: String, length: u64 },
    Writer(String),
    Io(String),
    UnknownRecordContig { name: String },
    PositionOverflow { sequence_id: String, position: u64 },
    SampleCountMismatch { sequence_id: String, position: u64 },
}

impl std::fmt::Display for VcfContractError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for VcfContractError {}

/// Escribe una cabecera VCF válida y, en modo BGZF, exactamente un índice.
///
/// No escribe registros: F2.5.1 fija solamente la frontera de formato. El
/// llamador de F2.5.5 la integrará en el provisional y publicación atómica.
pub fn write_vcf_contract(
    destination: &Path,
    mode: VcfOutputMode,
    contract: &VcfHeaderContract,
) -> Result<VcfOutputPaths, VcfContractError> {
    match mode {
        VcfOutputMode::Text => {
            let vcf = destination.with_extension("vcf");
            write_vcf_records(File::create(&vcf).map_err(io_error)?, contract, &[])?;
            Ok(VcfOutputPaths {
                vcf,
                index: None,
                index_kind: None,
            })
        }
        VcfOutputMode::Bgzf => {
            let vcf = destination.with_extension("vcf.gz");
            let file = File::create(&vcf).map_err(io_error)?;
            let bgzf = noodles_bgzf::io::Writer::new(file);
            write_vcf_records(bgzf, contract, &[])?
                .finish()
                .map_err(io_error)?;

            let index_kind = contract.selected_index();
            let index = match index_kind {
                VcfIndexKind::Tbi => write_tbi(&vcf)?,
                VcfIndexKind::Csi => write_csi(&vcf, contract)?,
            };
            Ok(VcfOutputPaths {
                vcf,
                index: Some(index),
                index_kind: Some(index_kind),
            })
        }
    }
}

/// Escribe la cabecera y los registros VCF ya normalizados y agrupados.
///
/// La proyección biológica no se repite aquí: esta frontera solo serializa el
/// DTO de F2.5.3/4. `W` puede ser un fichero de texto o un escritor BGZF; en
/// ambos casos el encabezado se construye con `noodles-vcf`.
pub fn write_vcf_records<W>(
    destination: W,
    contract: &VcfHeaderContract,
    records: &[GroupedVcfRecord],
) -> Result<W, VcfContractError>
where
    W: Write,
{
    let header = build_header(contract, records_carry_tandem_repeats(records))?;
    let mut writer = vcf::io::Writer::new(destination);
    writer.write_header(&header).map_err(io_error)?;
    for record in records {
        let line = record_line(contract, record)?;
        writer
            .get_mut()
            .write_all(line.as_bytes())
            .map_err(io_error)?;
    }
    Ok(writer.into_inner())
}

/// Escribe exactamente el índice exigido para un VCF BGZF ya finalizado.
///
/// El llamador lo escribe como un artefacto separado del mismo provisional; la
/// publicación F1 los hace visibles juntos con el resto del árbol atómico.
pub fn write_vcf_index<W>(
    vcf_path: &Path,
    contract: &VcfHeaderContract,
    destination: W,
) -> Result<VcfIndexKind, VcfContractError>
where
    W: Write,
{
    match contract.selected_index() {
        VcfIndexKind::Tbi => {
            let index = vcf::fs::index(vcf_path).map_err(io_error)?;
            write_tbi_index(&index, destination)?;
            Ok(VcfIndexKind::Tbi)
        }
        VcfIndexKind::Csi => {
            let index = build_csi_index(vcf_path, contract)?;
            let mut writer = csi::io::Writer::new(destination);
            writer.write_index(&index).map_err(io_error)?;
            writer.into_inner().finish().map_err(io_error)?;
            Ok(VcfIndexKind::Csi)
        }
    }
}

/// Hay que declarar los campos de repetición si algún registro los usa.
///
/// Se decide por lo **publicado**, no por lo declarado en la configuración,
/// igual que la promoción de versión de la verdad: un escenario cuyo único
/// alelo es el de la referencia no produce ningún registro de repetición y su
/// VCF sigue siendo byte a byte el de antes.
#[must_use]
pub fn records_carry_tandem_repeats(records: &[GroupedVcfRecord]) -> bool {
    records.iter().any(|record| {
        record
            .alternate_repeats
            .iter()
            .any(|sequences| !sequences.is_empty())
    })
}

fn build_header(
    contract: &VcfHeaderContract,
    declares_tandem_repeats: bool,
) -> Result<vcf::Header, VcfContractError> {
    let mut builder = vcf::Header::builder().set_file_format(FileFormat::new(4, 5));
    for contig in contract.contigs() {
        let length =
            usize::try_from(contig.length()).map_err(|_| VcfContractError::LengthTooLarge {
                name: contig.name().to_owned(),
                length: contig.length(),
            })?;
        let declaration = Map::<Contig>::builder()
            .set_length(length)
            .build()
            .map_err(|error| VcfContractError::Writer(error.to_string()))?;
        builder = builder.add_contig(contig.name(), declaration);
    }
    for sample in contract.samples() {
        builder = builder.add_sample_name(sample);
    }
    builder = builder
        .insert(
            VCF_REFERENCE_META_KEY.parse().map_err(
                |error: noodles_vcf::header::record::key::other::ParseError| {
                    VcfContractError::Writer(error.to_string())
                },
            )?,
            contract
                .identity()
                .reference_id()
                .as_str()
                .to_owned()
                .into(),
        )
        .map_err(|error| VcfContractError::Writer(error.to_string()))?
        .insert(
            VCF_EXECUTION_META_KEY.parse().map_err(
                |error: noodles_vcf::header::record::key::other::ParseError| {
                    VcfContractError::Writer(error.to_string())
                },
            )?,
            contract
                .identity()
                .execution_id()
                .as_str()
                .to_owned()
                .into(),
        )
        .map_err(|error| VcfContractError::Writer(error.to_string()))?
        .insert(
            VCF_QUAL_META_KEY.parse().map_err(
                |error: noodles_vcf::header::record::key::other::ParseError| {
                    VcfContractError::Writer(error.to_string())
                },
            )?,
            VCF_QUAL_SEMANTICS.into(),
        )
        .map_err(|error| VcfContractError::Writer(error.to_string()))?
        .add_filter(
            "PASS",
            Map::<noodles_vcf::header::record::value::map::Filter>::new(
                VCF_PASS_FILTER_DESCRIPTION,
            ),
        )
        .add_info(
            "GPS_CASE",
            Map::<noodles_vcf::header::record::value::map::Info>::new(
                noodles_vcf::header::record::value::map::info::Number::Unknown,
                noodles_vcf::header::record::value::map::info::Type::String,
                "Canonical truth case identifiers contributing an ALT allele",
            ),
        );
    if declares_tandem_repeats {
        builder = add_tandem_repeat_info(builder);
    }
    builder = builder.add_format(
        "GT",
        Map::<noodles_vcf::header::record::value::map::Format>::from("GT"),
    );
    Ok(builder.build())
}

/// Declara los cinco campos de repetición en tándem del propio `VCFv4.5`.
///
/// Número, tipo y descripción son **literales de la especificación**, sección
/// «INFO keys used for structural variants».
///
/// # De dónde salen estos campos, dicho sin adornos
///
/// El estándar los define **para el alelo simbólico de repetición en tándem**:
/// su texto dice «Used by `<CNV:TR>` tandem repeat alleles». Aquí se emiten
/// sobre alelos **no simbólicos**, porque un simulador siempre conoce la
/// secuencia exacta y la decisión S-09 fija que la forma explícita aplica
/// siempre. El estándar no lo prohíbe y contempla expresamente el caso —la
/// regla de omisión del recuento habla de «los demás» alelos, que son los no
/// simbólicos—, pero conviene que quede escrito: es el sustrato de la única
/// pregunta abierta sobre esta proyección.
///
/// # Por qué se publica la longitud del motivo, que el estándar sugiere omitir
///
/// El estándar dice, literal, «RUL should be omitted when RUS is present (as
/// it is redundant when RS is present)». Se publica igualmente, y por tres
/// razones:
///
/// - **es una recomendación, no una obligación**, y está en las notas sobre
///   los registros simbólicos de resumen; el propio estándar fija la regla de
///   consistencia para cuando aparecen los dos —la longitud de cada motivo
///   tiene que coincidir con la de su secuencia—, que es la que aquí se cumple;
/// - **su motivo declarado es la redundancia**, y aquí la redundancia se
///   quiere: es una comprobación viva que el oráculo ejerce. Conviene no
///   confundir esta recomendación con la que la precede —«STRs should use RUS,
///   whereas VNTRs should use RUL to ensure the VCF records are not
///   excessively large»—, que es la que habla del tamaño del registro y que
///   trata de elegir entre los dos campos, no de omitir este. Dicho eso, el
///   tamaño tampoco es un problema: estos registros ya llevan el tracto entero
///   en `REF` y `ALT`, de modo que al lado de un alelo de diez mil bases una
///   lista de longitudes de motivo no pesa nada;
/// - la decisión S-09 **enumera la longitud del motivo** entre los campos que
///   se publican, y quitarla sería contradecir una decisión normativa para
///   seguir una recomendación que no lo es. Además es una comprobación viva:
///   el oráculo la compara contra la longitud de su motivo.
///
/// Dos precisiones que no se pueden perder:
///
/// - el número de copias (`RUC`) es **decimal**, porque un tracto puede
///   terminar en una unidad parcial —`(CAG)5(CA)1(CAG)4` es el ejemplo del
///   propio estándar—; declararlo entero produciría un VCF que una herramienta
///   conforme rechaza;
/// - el recuento de secuencias repetidas (`RN`) vale **uno por alelo
///   alternativo**, y el estándar dice que si se omite se asume uno para los
///   alelos simbólicos de repetición y **cero para los demás**. Omitirlo no
///   sería callar, sino afirmar cero, de modo que se emite explícito para cada
///   alelo alternativo. El alelo de referencia no se describe con estos
///   campos: `RN` es de número `A` y no alcanza a `REF`.
fn add_tandem_repeat_info(builder: vcf::header::Builder) -> vcf::header::Builder {
    use noodles_vcf::header::record::value::map::Info;
    use noodles_vcf::header::record::value::map::info::{Number, Type};

    builder
        .add_info(
            "RN",
            Map::<Info>::new(
                Number::AlternateBases,
                Type::Integer,
                "Total number of repeat sequences in this allele",
            ),
        )
        .add_info(
            "RUS",
            Map::<Info>::new(
                Number::Unknown,
                Type::String,
                "Repeat unit sequence of the corresponding repeat sequence",
            ),
        )
        .add_info(
            "RUL",
            Map::<Info>::new(
                Number::Unknown,
                Type::Integer,
                "Repeat unit length of the corresponding repeat sequence",
            ),
        )
        .add_info(
            "RUC",
            Map::<Info>::new(
                Number::Unknown,
                Type::Float,
                "Repeat unit count of corresponding repeat sequence",
            ),
        )
        .add_info(
            "RB",
            Map::<Info>::new(
                Number::Unknown,
                Type::Integer,
                "Total number of bases in the corresponding repeat sequence",
            ),
        )
}

fn record_line(
    contract: &VcfHeaderContract,
    record: &GroupedVcfRecord,
) -> Result<String, VcfContractError> {
    if !contract
        .contigs()
        .iter()
        .any(|contig| contig.name() == record.sequence_id)
    {
        return Err(VcfContractError::UnknownRecordContig {
            name: record.sequence_id.clone(),
        });
    }
    if record.samples.len() != contract.samples().len()
        || record
            .samples
            .iter()
            .map(|(sample, _)| sample)
            .ne(contract.samples().iter())
    {
        return Err(VcfContractError::SampleCountMismatch {
            sequence_id: record.sequence_id.clone(),
            position: record.position,
        });
    }
    let position =
        record
            .position
            .checked_add(1)
            .ok_or_else(|| VcfContractError::PositionOverflow {
                sequence_id: record.sequence_id.clone(),
                position: record.position,
            })?;
    let samples = record
        .samples
        .iter()
        .map(|(_, genotype)| {
            genotype
                .alleles
                .iter()
                .map(u32::to_string)
                .collect::<Vec<_>>()
                .join("|")
        })
        .collect::<Vec<_>>()
        .join("\t");
    let mut info = format!("GPS_CASE={}", record.gps_cases.join(","));
    info.push_str(&tandem_repeat_info(record));
    Ok(format!(
        "{}\t{}\t.\t{}\t{}\t.\tPASS\t{info}\tGT\t{}\n",
        record.sequence_id,
        position,
        record.reference,
        record.alternates.join(","),
        samples,
    ))
}

/// Los cinco campos de repetición del registro, o la cadena vacía.
///
/// Las cuatro listas por secuencia repetida van **aplanadas y concatenadas en
/// el orden de los alelos alternativos**, que es exactamente lo que el
/// estándar define: el recuento de cada alelo dice cuántas entradas le
/// corresponden, y la longitud total de cada lista es la suma de esos
/// recuentos. Un alelo alternativo que no procede de un tracto aporta cero
/// entradas y su recuento vale cero.
fn tandem_repeat_info(record: &GroupedVcfRecord) -> String {
    if record
        .alternate_repeats
        .iter()
        .all(|sequences| sequences.is_empty())
    {
        return String::new();
    }
    let counts = record
        .alternate_repeats
        .iter()
        .map(|sequences| sequences.len().to_string())
        .collect::<Vec<_>>()
        .join(",");
    let flattened = record
        .alternate_repeats
        .iter()
        .flatten()
        .collect::<Vec<_>>();
    let motifs = flattened
        .iter()
        .map(|sequence| sequence.motif.clone())
        .collect::<Vec<_>>()
        .join(",");
    let motif_lengths = flattened
        .iter()
        .map(|sequence| motif_length(sequence).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let unit_counts = flattened
        .iter()
        .map(|sequence| repeat_unit_count(sequence).to_string())
        .collect::<Vec<_>>()
        .join(",");
    let bases = flattened
        .iter()
        .map(|sequence| sequence.bases.to_string())
        .collect::<Vec<_>>()
        .join(",");
    format!(";RN={counts};RUS={motifs};RUL={motif_lengths};RUC={unit_counts};RB={bases}")
}

/// La longitud del motivo, que el estándar exige igual a la de su secuencia.
///
/// La capa de representación ya ha rechazado el motivo vacío, de modo que
/// nunca es cero.
fn motif_length(sequence: &VcfRepeatSequence) -> u64 {
    sequence.motif.len() as u64
}

/// El número de copias de una secuencia repetida: sus bases entre la longitud
/// de su motivo.
///
/// Es decimal a propósito y la división se hace aquí, no se guarda: guardar el
/// cociente junto a sus dos operandos permitiría que se contradijeran.
fn repeat_unit_count(sequence: &VcfRepeatSequence) -> f64 {
    sequence.bases as f64 / motif_length(sequence) as f64
}

fn write_tbi(vcf_path: &Path) -> Result<PathBuf, VcfContractError> {
    let index = vcf::fs::index(vcf_path).map_err(io_error)?;
    let path = PathBuf::from(format!("{}.tbi", vcf_path.display()));
    write_tbi_index(&index, File::create(&path).map_err(io_error)?)?;
    Ok(path)
}

fn write_tbi_index(index: &tabix::Index, destination: impl Write) -> Result<(), VcfContractError> {
    let mut writer = tabix::io::Writer::new(destination);
    writer.write_index(index).map_err(io_error)?;
    writer.into_inner().finish().map_err(io_error)?;
    Ok(())
}

fn write_csi(vcf_path: &Path, contract: &VcfHeaderContract) -> Result<PathBuf, VcfContractError> {
    let path = PathBuf::from(format!("{}.csi", vcf_path.display()));
    write_vcf_index(vcf_path, contract, File::create(&path).map_err(io_error)?)?;
    Ok(path)
}

fn build_csi_index(
    vcf_path: &Path,
    contract: &VcfHeaderContract,
) -> Result<csi::Index, VcfContractError> {
    use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;

    let mut header = csi::binning_index::index::header::Builder::vcf().build();
    for contig in contract.contigs() {
        header
            .reference_sequence_names_mut()
            .insert(contig.name().into());
    }
    let mut reader = vcf::io::Reader::new(noodles_bgzf::io::Reader::new(
        File::open(vcf_path).map_err(io_error)?,
    ));
    let vcf_header = reader.read_header().map_err(io_error)?;
    let mut indexer = Indexer::<BinnedIndex>::new(14, 6).set_header(header);
    let mut record = vcf::Record::default();
    let mut start_position = reader.get_ref().virtual_position();
    while reader.read_record(&mut record).map_err(io_error)? != 0 {
        let end_position = reader.get_ref().virtual_position();
        let name = record.reference_sequence_name();
        let reference_sequence_id = contract
            .contigs()
            .iter()
            .position(|contig| contig.name() == name)
            .ok_or_else(|| VcfContractError::UnknownRecordContig {
                name: name.to_owned(),
            })?;
        let start = record
            .variant_start()
            .transpose()
            .map_err(io_error)?
            .ok_or_else(|| VcfContractError::Writer("VCF record has no POS".to_owned()))?;
        let end = record.variant_end(&vcf_header).map_err(io_error)?;
        indexer
            .add_record(
                Some((reference_sequence_id, start, end, true)),
                Chunk::new(start_position, end_position),
            )
            .map_err(io_error)?;
        start_position = end_position;
    }
    Ok(indexer.build(contract.contigs().len()))
}

fn validate_identifier(kind: &'static str, value: &str) -> Result<(), VcfContractError> {
    if value.is_empty() {
        return Err(VcfContractError::EmptyIdentifier { kind });
    }
    if value
        .bytes()
        .any(|byte| byte.is_ascii_whitespace() || byte == b'\0')
    {
        return Err(VcfContractError::InvalidIdentifier {
            kind,
            value: value.to_owned(),
        });
    }
    Ok(())
}

fn io_error(error: io::Error) -> VcfContractError {
    VcfContractError::Io(error.to_string())
}

#[cfg(test)]
mod tests {
    use std::io::{self, Write};

    use super::*;

    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("injected finalization failure"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn tbi_finalization_error_is_propagated() {
        let index = tabix::index::Indexer::default().build();
        let error =
            write_tbi_index(&index, FailingWriter).expect_err("the BGZF EOF write must fail");
        assert!(matches!(error, VcfContractError::Io(detail) if detail.contains("injected")));
    }
}
