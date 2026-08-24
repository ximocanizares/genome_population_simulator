//! Contrato físico mínimo VCF 4.5 para F2.5.1.
//!
//! Este módulo solo escribe la cabecera y el envoltorio de salida. La
//! proyección de alelos, `GT`, campos `INFO` y registros pertenece a los pasos
//! F2.5.2 y posteriores. En particular, no consulta `truth.jsonl`.

use std::collections::BTreeSet;
use std::fs::File;
use std::io::{self, Write};
use std::path::{Path, PathBuf};

use noodles_csi::{
    self as csi,
    binning_index::{Indexer, index::reference_sequence::index::BinnedIndex},
};
use noodles_tabix as tabix;
use noodles_vcf::{
    self as vcf,
    header::{
        FileFormat,
        record::value::{Map, map::Contig},
    },
};

use crate::identity::SequenceId;

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
    let header = build_header(contract)?;
    match mode {
        VcfOutputMode::Text => {
            let vcf = destination.with_extension("vcf");
            let mut writer = vcf::io::Writer::new(File::create(&vcf).map_err(io_error)?);
            writer.write_header(&header).map_err(io_error)?;
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
            let mut writer = vcf::io::Writer::new(bgzf);
            writer.write_header(&header).map_err(io_error)?;
            writer.into_inner().finish().map_err(io_error)?;

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
    Ok(builder.build())
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
    let mut header = csi::binning_index::index::header::Builder::vcf().build();
    for contig in contract.contigs() {
        header
            .reference_sequence_names_mut()
            .insert(contig.name().into());
    }
    let index = Indexer::<BinnedIndex>::new(14, 6)
        .set_header(header)
        .build(contract.contigs().len());
    let path = PathBuf::from(format!("{}.csi", vcf_path.display()));
    let mut writer = csi::io::Writer::new(File::create(&path).map_err(io_error)?);
    writer.write_index(&index).map_err(io_error)?;
    writer.into_inner().finish().map_err(io_error)?;
    Ok(path)
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
