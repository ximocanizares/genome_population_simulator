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

use crate::identity::SequenceId;
use crate::vcf_grouping::GroupedVcfRecord;

pub const VCF_FILE_FORMAT: &str = "VCFv4.5";
pub const TBI_MAXIMUM_POSITION_EXCLUSIVE: u64 = 1 << 29;

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
    contigs: Vec<VcfContig>,
    samples: Vec<String>,
}

impl VcfHeaderContract {
    pub fn new(contigs: Vec<VcfContig>, samples: Vec<String>) -> Result<Self, VcfContractError> {
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
        Ok(Self { contigs, samples })
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
    let header = build_header(contract)?;
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

fn build_header(contract: &VcfHeaderContract) -> Result<vcf::Header, VcfContractError> {
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
        .add_info(
            "GPS_CASE",
            Map::<noodles_vcf::header::record::value::map::Info>::new(
                noodles_vcf::header::record::value::map::info::Number::Unknown,
                noodles_vcf::header::record::value::map::info::Type::String,
                "Canonical truth case identifiers contributing an ALT allele",
            ),
        )
        .add_format(
            "GT",
            Map::<noodles_vcf::header::record::value::map::Format>::from("GT"),
        );
    Ok(builder.build())
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
    Ok(format!(
        "{}\t{}\t.\t{}\t{}\t.\tPASS\tGPS_CASE={}\tGT\t{}\n",
        record.sequence_id,
        position,
        record.reference,
        record.alternates.join(","),
        record.gps_cases.join(","),
        samples,
    ))
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
