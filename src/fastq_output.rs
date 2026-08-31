//! Escritura FASTQ mínima por flujo desde orígenes WGS ya faseados.
//!
//! F3.1.4 solo convierte cada origen transitorio en su registro FASTQ final.
//! No conserva lecturas ni procedencia por lectura: al cerrar quedan el
//! catálogo de artefactos de `ExecutionWorkspace` y dos contadores agregados.

use std::fs::File;
use std::io::{BufRead, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};

use crate::application::{
    ILLUMINA_SE_MINIMAL_ACGT_PHRED, ILLUMINA_SE_MINIMAL_N_PHRED, ILLUMINA_SE_MINIMAL_READ_LENGTH,
};
use crate::demultiplexing::DemultiplexedPair;
use crate::diagnostics::Diagnostic;
use crate::error_quality::ObservedPair;
use crate::fastq_validation::validate_final_paired_fastq;
use crate::identity::{LibraryId, compose_template_id};
use crate::publication::ExecutionWorkspace;
use crate::results::ArtifactRecord;
use crate::technical_identity::TechnicalMate;
use crate::wgs_sampling::{WGS_READ_LENGTH, WgsLibrarySampler, WgsSamplingError};
use noodles_bgzf as bgzf;

pub const FASTQ_ARTIFACT_KIND: &str = "fastq";
const FASTQ_BUFFER_CAPACITY: usize = 64 * 1024;
/// Perfil físico fijo de FASTQ comprimido v1: BGZF, que es un flujo gzip
/// interoperable, con DEFLATE nivel 6. BGZF fija `MTIME = 0` y el tamaño de
/// bloque, para que los bytes sean repetibles bajo el lock y la plataforma
/// admitidos.
const FASTQ_GZIP_COMPRESSION_LEVEL: u8 = 6;
/// § 4.1.2 de la especificación BGZF (2020-12-03). `noodles-bgzf` lo emite
/// al terminar el escritor, pero su constante no forma parte de la API pública.
const BGZF_EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];

#[cfg(test)]
thread_local! {
    static FORCE_FINAL_R2_FAILURE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
    static FORCE_FINAL_PAIR_VALIDATION_FAILURE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

/// Codifica calidades numéricas estrictas como FASTQ Sanger (Phred+33).
///
/// F4.2 usa esta frontera para las vistas observadas transitorias; F3 conserva
/// su escritor propio, que ya produce el mismo alfabeto para su perfil fijo.
pub fn encode_fastq_sanger_qualities(phred: &[u8]) -> Result<Vec<u8>, Diagnostic> {
    phred
        .iter()
        .copied()
        .map(|value| {
            if value > 93 {
                Err(invariant(
                    "encode_fastq_sanger_qualities",
                    "Phred value is outside 0..=93",
                ))
            } else {
                Ok(value + 33)
            }
        })
        .collect()
}

/// Artefacto FASTQ final de la única biblioteca de un individuo en F3.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastqPublication {
    artifact: ArtifactRecord,
    read_count: u64,
    base_count: u64,
}

/// Los dos artefactos físicos y el único resumen agregado que sobreviven a la
/// publicación PE. No conserva claves PH, diferencias por base ni relaciones
/// de procedencia de las parejas que ya se escribieron.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FinalPairedFastqPublication {
    r1: ArtifactRecord,
    r2: ArtifactRecord,
    pair_count: u64,
    read_count: u64,
    base_count: u64,
    own_pair_count: u64,
    contaminant_pair_count: u64,
    own_base_count: u64,
    contaminant_base_count: u64,
    n_base_count: u64,
}

impl FinalPairedFastqPublication {
    #[must_use]
    pub const fn r1(&self) -> &ArtifactRecord {
        &self.r1
    }
    #[must_use]
    pub const fn r2(&self) -> &ArtifactRecord {
        &self.r2
    }
    #[must_use]
    pub const fn pair_count(&self) -> u64 {
        self.pair_count
    }
    #[must_use]
    pub const fn read_count(&self) -> u64 {
        self.read_count
    }
    #[must_use]
    pub const fn base_count(&self) -> u64 {
        self.base_count
    }
    #[must_use]
    pub const fn own_pair_count(&self) -> u64 {
        self.own_pair_count
    }
    #[must_use]
    pub const fn contaminant_pair_count(&self) -> u64 {
        self.contaminant_pair_count
    }
    #[must_use]
    pub const fn own_base_count(&self) -> u64 {
        self.own_base_count
    }
    #[must_use]
    pub const fn contaminant_base_count(&self) -> u64 {
        self.contaminant_base_count
    }
    #[must_use]
    pub const fn n_base_count(&self) -> u64 {
        self.n_base_count
    }
}

impl FastqPublication {
    #[must_use]
    pub const fn artifact(&self) -> &ArtifactRecord {
        &self.artifact
    }

    #[must_use]
    pub const fn read_count(&self) -> u64 {
        self.read_count
    }

    #[must_use]
    pub const fn base_count(&self) -> u64 {
        self.base_count
    }
}

/// Materializa un FASTQ single-end final para una biblioteca WGS.
///
/// Cada lectura se materializa, se escribe y se descarta dentro del mismo
/// ciclo. El temporal gzip/BGZF se cierra y se relee descomprimido antes de
/// calcular SHA-256 sobre sus bytes comprimidos y renombrarlo como artefacto
/// del provisional, comprobando el formato, el contador y cada ID y calidad
/// canónicos.
pub fn write_wgs_fastq(
    workspace: &mut ExecutionWorkspace,
    sampler: &mut WgsLibrarySampler<'_>,
) -> Result<FastqPublication, Diagnostic> {
    let event_id = sampler.event_id().as_str().to_owned();
    let individual_id = sampler.individual_id().as_str().to_owned();
    let library_id = sampler
        .library_id()
        .ok_or_else(|| invariant("write_wgs_fastq", "the sampler has no technical library"))?
        .clone();
    let expected_count = sampler.plan().read_count();
    let relative_path = fastq_relative_path(&event_id, &individual_id);
    let artifact_id = format!("fastq-{}", library_id.as_str());
    let path_for_write = relative_path.clone();
    let path_for_validation = relative_path.clone();
    let mut summary = FastqSummary::default();

    let artifact = workspace.write_streamed_artifact_with_validation(
        &artifact_id,
        FASTQ_ARTIFACT_KIND,
        &relative_path,
        |file| {
            let level = bgzf::io::writer::CompressionLevel::new(FASTQ_GZIP_COMPRESSION_LEVEL)
                .ok_or_else(|| invariant("write_wgs_fastq", "invalid fixed gzip level"))?;
            let mut compressed = bgzf::io::writer::Builder::default()
                .set_compression_level(level)
                .build_from_writer(file);
            summary = write_records(&mut compressed, sampler, expected_count, &path_for_write)?;
            compressed
                .finish()
                .map_err(|error| io_error(&path_for_write, "finish_fastq_gzip", error))?;
            Ok(())
        },
        |temporary| {
            validate_temporary_fastq(temporary, &path_for_validation, &library_id, expected_count)
        },
    )?;

    Ok(FastqPublication {
        artifact,
        read_count: summary.read_count,
        base_count: summary.base_count,
    })
}

#[must_use]
pub fn fastq_relative_path(event_id: &str, individual_id: &str) -> PathBuf {
    PathBuf::from("outputs")
        .join("fastq")
        .join(event_id)
        .join(format!("{individual_id}.fastq.gz"))
}

/// Rutas finales reservadas para una pareja PE. F4.1 solo las define junto al
/// productor de fragmentos; F4.5 escribirá las cabeceras finales después del
/// demultiplexado, sin volver a introducir la ruta SE de F3.
#[must_use]
pub fn paired_fastq_relative_paths(event_id: &str, individual_id: &str) -> (PathBuf, PathBuf) {
    let base = PathBuf::from("outputs").join("fastq").join(event_id);
    (
        base.join(format!("{individual_id}_R1.fastq.gz")),
        base.join(format!("{individual_id}_R2.fastq.gz")),
    )
}

/// Publica la única pareja FASTQ PE posterior al demultiplexado para una
/// biblioteca. Las parejas se consumen tras escribir R1/R2: nunca se genera
/// una copia limpia ni un fichero de procedencia.
pub fn write_final_paired_fastq(
    workspace: &mut ExecutionWorkspace,
    event_id: &str,
    individual_id: &str,
    library_id: &LibraryId,
    pairs: Vec<DemultiplexedPair<ObservedPair>>,
) -> Result<FinalPairedFastqPublication, Diagnostic> {
    let (r1_path, r2_path) = paired_fastq_relative_paths(event_id, individual_id);
    let summary = summarize_final_pairs(&pairs)?;
    let r1 = write_final_mate(
        workspace,
        &format!("fastq-{}-r1", library_id.as_str()),
        &r1_path,
        event_id,
        library_id,
        &pairs,
        TechnicalMate::R1,
    )?;
    let r2 = match write_final_mate(
        workspace,
        &format!("fastq-{}-r2", library_id.as_str()),
        &r2_path,
        event_id,
        library_id,
        &pairs,
        TechnicalMate::R2,
    ) {
        Ok(artifact) => artifact,
        Err(error) => {
            workspace.discard_last_artifact(&r1)?;
            return Err(error);
        }
    };
    if let Err(error) = validate_final_paired_fastq(
        &workspace.root().join(&r1_path),
        &r1,
        &workspace.root().join(&r2_path),
        &r2,
    )
    .and(maybe_inject_final_pair_validation_failure())
    {
        // R1/R2 son una única salida lógica: una falla conjunta no puede
        // dejar artefactos ni catálogo parciales en el provisional.
        workspace.discard_last_artifact(&r2)?;
        workspace.discard_last_artifact(&r1)?;
        return Err(invariant("write_final_paired_fastq", &error.to_string()));
    }
    Ok(FinalPairedFastqPublication {
        r1,
        r2,
        pair_count: summary.pair_count,
        read_count: summary.read_count,
        base_count: summary.base_count,
        own_pair_count: summary.own_pair_count,
        contaminant_pair_count: summary.contaminant_pair_count,
        own_base_count: summary.own_base_count,
        contaminant_base_count: summary.contaminant_base_count,
        n_base_count: summary.n_base_count,
    })
}

fn maybe_inject_final_pair_validation_failure()
-> Result<(), crate::fastq_validation::FastqValidationError> {
    #[cfg(test)]
    if FORCE_FINAL_PAIR_VALIDATION_FAILURE.replace(false) {
        return Err(crate::fastq_validation::FastqValidationError::Structural {
            path: "injected-final-pair".to_owned(),
            record_ordinal: 0,
            detail: "injected final-pair validation failure",
        });
    }
    Ok(())
}

#[derive(Clone, Copy, Debug, Default)]
struct FinalPairSummary {
    pair_count: u64,
    read_count: u64,
    base_count: u64,
    own_pair_count: u64,
    contaminant_pair_count: u64,
    own_base_count: u64,
    contaminant_base_count: u64,
    n_base_count: u64,
}

fn write_final_mate(
    workspace: &mut ExecutionWorkspace,
    artifact_id: &str,
    relative_path: &Path,
    event_id: &str,
    library_id: &LibraryId,
    pairs: &[DemultiplexedPair<ObservedPair>],
    mate: TechnicalMate,
) -> Result<ArtifactRecord, Diagnostic> {
    #[cfg(test)]
    if mate == TechnicalMate::R2 && FORCE_FINAL_R2_FAILURE.replace(false) {
        return Err(invariant(
            "write_final_paired_fastq",
            "injected R2 publication failure",
        ));
    }
    let path = relative_path.to_owned();
    workspace.write_streamed_artifact_with_validation(
        artifact_id,
        FASTQ_ARTIFACT_KIND,
        relative_path,
        |file| {
            let level = bgzf::io::writer::CompressionLevel::new(FASTQ_GZIP_COMPRESSION_LEVEL)
                .ok_or_else(|| invariant("write_final_paired_fastq", "invalid fixed gzip level"))?;
            let mut compressed = bgzf::io::writer::Builder::default()
                .set_compression_level(level)
                .build_from_writer(file);
            for pair in pairs {
                if pair.final_library_id() != library_id {
                    return Err(invariant(
                        "write_final_paired_fastq",
                        "pair belongs to a different final library",
                    ));
                }
                write_final_record(&mut compressed, pair, mate, event_id, library_id)?;
            }
            compressed
                .finish()
                .map_err(|error| io_error(&path, "finish_final_fastq_gzip", error))?;
            Ok(())
        },
        |temporary| {
            crate::fastq_validation::validate_final_fastq_mate(temporary, mate)
                .map_err(|error| invariant("write_final_paired_fastq", &error.to_string()))
        },
    )
}

fn write_final_record(
    output: &mut impl Write,
    pair: &DemultiplexedPair<ObservedPair>,
    mate: TechnicalMate,
    event_id: &str,
    library_id: &LibraryId,
) -> Result<(), Diagnostic> {
    let observed = match mate {
        TechnicalMate::R1 => pair.payload().r1(),
        TechnicalMate::R2 => pair.payload().r2(),
    };
    if observed.mate() != mate {
        return Err(invariant(
            "write_final_paired_fastq",
            "observed mate is inconsistent",
        ));
    }
    output
        .write_all(
            format!(
                "@GPS:{event_id}:{}:{} {}:N:0:{}\n",
                library_id.as_str(),
                pair.technical_identity().technical_pair_id().as_str(),
                mate.number(),
                pair.final_udi().literal(),
            )
            .as_bytes(),
        )
        .map_err(|error| io_error(Path::new("final-fastq"), "write_header", error))?;
    output
        .write_all(observed.sequence())
        .and_then(|()| output.write_all(b"\n+\n"))
        .map_err(|error| io_error(Path::new("final-fastq"), "write_sequence", error))?;
    let qualities = encode_fastq_sanger_qualities(observed.phred())?;
    output
        .write_all(&qualities)
        .and_then(|()| output.write_all(b"\n"))
        .map_err(|error| io_error(Path::new("final-fastq"), "write_quality", error))
}

fn summarize_final_pairs(
    pairs: &[DemultiplexedPair<ObservedPair>],
) -> Result<FinalPairSummary, Diagnostic> {
    let mut summary = FinalPairSummary::default();
    for pair in pairs {
        let r1_bases = u64::try_from(pair.payload().r1().sequence().len())
            .map_err(|_| invariant("write_final_paired_fastq", "read length overflows u64"))?;
        let r2_bases = u64::try_from(pair.payload().r2().sequence().len())
            .map_err(|_| invariant("write_final_paired_fastq", "read length overflows u64"))?;
        let bases = r1_bases
            .checked_add(r2_bases)
            .ok_or_else(|| invariant("write_final_paired_fastq", "read length overflows u64"))?;
        summary.pair_count = summary
            .pair_count
            .checked_add(1)
            .ok_or_else(|| invariant("write_final_paired_fastq", "pair count overflows u64"))?;
        summary.read_count = summary
            .read_count
            .checked_add(2)
            .ok_or_else(|| invariant("write_final_paired_fastq", "read count overflows u64"))?;
        summary.base_count = summary
            .base_count
            .checked_add(bases)
            .ok_or_else(|| invariant("write_final_paired_fastq", "base count overflows u64"))?;
        if pair.is_contaminant() {
            summary.contaminant_pair_count = summary
                .contaminant_pair_count
                .checked_add(1)
                .ok_or_else(|| {
                    invariant(
                        "write_final_paired_fastq",
                        "contaminant pair count overflows u64",
                    )
                })?;
            summary.contaminant_base_count = summary
                .contaminant_base_count
                .checked_add(bases)
                .ok_or_else(|| {
                    invariant(
                        "write_final_paired_fastq",
                        "contaminant base count overflows u64",
                    )
                })?;
        } else {
            summary.own_pair_count = summary.own_pair_count.checked_add(1).ok_or_else(|| {
                invariant("write_final_paired_fastq", "own pair count overflows u64")
            })?;
            summary.own_base_count =
                summary.own_base_count.checked_add(bases).ok_or_else(|| {
                    invariant("write_final_paired_fastq", "own base count overflows u64")
                })?;
        }
        summary.n_base_count = summary
            .n_base_count
            .checked_add(
                u64::try_from(
                    pair.payload()
                        .r1()
                        .sequence()
                        .iter()
                        .chain(pair.payload().r2().sequence())
                        .filter(|base| **base == b'N')
                        .count(),
                )
                .map_err(|_| invariant("write_final_paired_fastq", "N base count overflows u64"))?,
            )
            .ok_or_else(|| invariant("write_final_paired_fastq", "N base count overflows u64"))?;
    }
    Ok(summary)
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct FastqSummary {
    read_count: u64,
    base_count: u64,
}

fn write_records(
    file: &mut impl Write,
    sampler: &mut WgsLibrarySampler<'_>,
    expected_count: u64,
    relative_path: &Path,
) -> Result<FastqSummary, Diagnostic> {
    let mut writer = BufWriter::with_capacity(FASTQ_BUFFER_CAPACITY, file);
    let mut summary = FastqSummary::default();
    let mut qualities = [0_u8; ILLUMINA_SE_MINIMAL_READ_LENGTH as usize];

    for ordinal in 0..expected_count {
        let read = sampler
            .next()
            .ok_or_else(|| invariant("write_records", "sampler ended before its read plan"))?
            .map_err(sampling_error)?;
        if read.template().template_ordinal().get() != ordinal {
            return Err(invariant(
                "write_records",
                "sampler emitted a non-canonical template ordinal",
            ));
        }
        let sequence = read.true_sequence();
        if sequence.len() != qualities.len() {
            return Err(invariant(
                "write_records",
                "sampled sequence has a non-contractual read length",
            ));
        }
        for (index, base) in sequence.iter().copied().enumerate() {
            qualities[index] = quality_for_base(base).ok_or_else(|| {
                invariant(
                    "write_records",
                    "sampled sequence contains an invalid DNA base",
                )
            })?;
        }

        writer
            .write_all(b"@")
            .and_then(|()| writer.write_all(read.template().template_id().as_str().as_bytes()))
            .and_then(|()| writer.write_all(b"\n"))
            .and_then(|()| writer.write_all(sequence))
            .and_then(|()| writer.write_all(b"\n+\n"))
            .and_then(|()| writer.write_all(&qualities))
            .and_then(|()| writer.write_all(b"\n"))
            .map_err(|error| io_error(relative_path, "write_fastq_record", error))?;
        summary.read_count = summary
            .read_count
            .checked_add(1)
            .ok_or_else(|| invariant("write_records", "FASTQ read count exceeds u64"))?;
        summary.base_count = summary
            .base_count
            .checked_add(WGS_READ_LENGTH)
            .ok_or_else(|| invariant("write_records", "FASTQ base count exceeds u64"))?;
    }
    if sampler.next().is_some() {
        return Err(invariant(
            "write_records",
            "sampler yielded more reads than its read plan",
        ));
    }
    writer
        .flush()
        .map_err(|error| io_error(relative_path, "flush_fastq", error))?;
    Ok(summary)
}

fn validate_temporary_fastq(
    temporary: &Path,
    relative_path: &Path,
    library_id: &LibraryId,
    expected_count: u64,
) -> Result<(), Diagnostic> {
    validate_canonical_bgzf_eof(temporary, relative_path)?;
    let file = File::open(temporary).map_err(|error| io_error(temporary, "open_fastq", error))?;
    // `bgzf::io::Reader` ya conserva un bloque descomprimido; envolverlo en
    // `BufReader` está prohibido por su contrato, pero satisface `BufRead` por
    // sí mismo y permite validar el FASTQ sin materializarlo.
    let mut reader = bgzf::io::Reader::new(file);
    let mut header = Vec::with_capacity(64);
    let mut sequence = Vec::with_capacity(usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1);
    let mut plus = Vec::with_capacity(2);
    let mut qualities = Vec::with_capacity(usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1);

    for ordinal in 0..expected_count {
        read_line(
            &mut reader,
            &mut header,
            temporary,
            relative_path,
            ordinal,
            "header",
            "read_fastq_header",
        )?;
        let expected_id = compose_template_id(library_id, ordinal);
        let mut expected_header = Vec::with_capacity(expected_id.as_str().len() + 2);
        expected_header.push(b'@');
        expected_header.extend_from_slice(expected_id.as_str().as_bytes());
        expected_header.push(b'\n');
        if header != expected_header {
            return Err(format_error(
                relative_path,
                ordinal,
                "FASTQ header is not the canonical technical template ID",
            ));
        }

        read_line(
            &mut reader,
            &mut sequence,
            temporary,
            relative_path,
            ordinal,
            "sequence",
            "read_fastq_sequence",
        )?;
        read_line(
            &mut reader,
            &mut plus,
            temporary,
            relative_path,
            ordinal,
            "separator",
            "read_fastq_plus",
        )?;
        read_line(
            &mut reader,
            &mut qualities,
            temporary,
            relative_path,
            ordinal,
            "quality",
            "read_fastq_quality",
        )?;
        validate_record_lines(relative_path, ordinal, &sequence, &plus, &qualities)?;
    }

    let mut trailing = [0_u8; 1];
    if reader
        .read(&mut trailing)
        .map_err(|error| io_error(temporary, "check_fastq_eof", error))?
        != 0
    {
        return Err(format_error(
            relative_path,
            expected_count,
            "FASTQ contains records beyond the planned count",
        ));
    }
    Ok(())
}

fn validate_canonical_bgzf_eof(temporary: &Path, relative_path: &Path) -> Result<(), Diagnostic> {
    let mut file = File::open(temporary)
        .map_err(|error| io_error(temporary, "open_fastq_eof_marker", error))?;
    let size = file
        .metadata()
        .map_err(|error| io_error(temporary, "metadata_fastq_eof_marker", error))?
        .len();
    let marker_length = u64::try_from(BGZF_EOF_MARKER.len()).map_err(|_| {
        invariant(
            "validate_canonical_bgzf_eof",
            "BGZF EOF marker length overflows",
        )
    })?;
    if size < marker_length {
        return Err(format_error(
            relative_path,
            0,
            "FASTQ BGZF stream is missing its canonical EOF marker",
        ));
    }
    file.seek(SeekFrom::End(-i64::try_from(marker_length).map_err(
        |_| {
            invariant(
                "validate_canonical_bgzf_eof",
                "BGZF EOF marker offset overflows",
            )
        },
    )?))
    .map_err(|error| io_error(temporary, "seek_fastq_eof_marker", error))?;
    let mut observed = [0_u8; BGZF_EOF_MARKER.len()];
    file.read_exact(&mut observed)
        .map_err(|error| io_error(temporary, "read_fastq_eof_marker", error))?;
    if observed != BGZF_EOF_MARKER {
        return Err(format_error(
            relative_path,
            0,
            "FASTQ BGZF stream is missing its canonical EOF marker",
        ));
    }
    Ok(())
}

fn read_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    temporary: &Path,
    relative_path: &Path,
    ordinal: u64,
    field: &str,
    operation: &str,
) -> Result<(), Diagnostic> {
    line.clear();
    let length = reader
        .read_until(b'\n', line)
        .map_err(|error| io_error(temporary, operation, error))?;
    if length == 0 || line.last() != Some(&b'\n') {
        return Err(format_error(
            relative_path,
            ordinal,
            &format!("FASTQ {field} line is truncated"),
        ));
    }
    Ok(())
}

fn validate_record_lines(
    relative_path: &Path,
    ordinal: u64,
    sequence: &[u8],
    plus: &[u8],
    qualities: &[u8],
) -> Result<(), Diagnostic> {
    let expected_length = usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1;
    if sequence.len() != expected_length || sequence.last() != Some(&b'\n') {
        return Err(format_error(
            relative_path,
            ordinal,
            "FASTQ sequence length does not match the configured read length",
        ));
    }
    if plus != b"+\n" {
        return Err(format_error(
            relative_path,
            ordinal,
            "FASTQ separator line is invalid",
        ));
    }
    if qualities.len() != expected_length || qualities.last() != Some(&b'\n') {
        return Err(format_error(
            relative_path,
            ordinal,
            "FASTQ quality length does not match its sequence",
        ));
    }
    for (base, quality) in sequence[..sequence.len() - 1]
        .iter()
        .copied()
        .zip(qualities[..qualities.len() - 1].iter().copied())
    {
        if quality_for_base(base) != Some(quality) {
            return Err(format_error(
                relative_path,
                ordinal,
                "FASTQ base or Phred quality is not valid for illumina-se-minimal-v1",
            ));
        }
    }
    Ok(())
}

fn quality_for_base(base: u8) -> Option<u8> {
    match base {
        b'A' | b'C' | b'G' | b'T' => Some(ILLUMINA_SE_MINIMAL_ACGT_PHRED + b'!'),
        b'N' => Some(ILLUMINA_SE_MINIMAL_N_PHRED + b'!'),
        _ => None,
    }
}

fn sampling_error(error: WgsSamplingError) -> Diagnostic {
    Diagnostic::internal_invariant("fastq_output", "sample_wgs_read")
        .with("detail", format!("{error:?}"))
}

fn format_error(path: &Path, ordinal: u64, detail: &str) -> Diagnostic {
    // F3.1.4 solo relee el temporal que acaba de producir este mismo proceso.
    // Una discrepancia aquí no procede de una entrada de usuario: revela una
    // invariante rota del productor o de su frontera de E/S. F3.5 clasificará
    // los FASTQ externos/releídos como errores de formato de su contrato.
    Diagnostic::internal_invariant("fastq_output", "validate_temporary_fastq")
        .with("path", path.display().to_string())
        .with("record_ordinal", ordinal.to_string())
        .with("detail", detail)
}

fn io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", error.to_string())
}

fn invariant(operation: &str, detail: &str) -> Diagnostic {
    Diagnostic::internal_invariant("fastq_output", operation).with("detail", detail)
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    use noodles_bgzf as bgzf;

    use super::validate_temporary_fastq;
    use super::{
        FORCE_FINAL_PAIR_VALIDATION_FAILURE, FORCE_FINAL_R2_FAILURE,
        ILLUMINA_SE_MINIMAL_READ_LENGTH, LibraryId, compose_template_id, write_final_paired_fastq,
    };
    use crate::demultiplexing::DemultiplexedPair;
    use crate::error_quality::ObservedPair;
    use crate::publication::ExecutionWorkspace;

    static FASTQ_VALIDATION_COUNTER: AtomicU64 = AtomicU64::new(0);

    struct InvalidFastqCase {
        name: &'static str,
        bytes: Vec<u8>,
        expected_count: u64,
        ordinal: u64,
        detail: &'static str,
    }

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-fastq-validation-unit-{}-{}-{label}",
            std::process::id(),
            FASTQ_VALIDATION_COUNTER.fetch_add(1, Ordering::Relaxed)
        ))
    }

    fn valid_record(library_id: &LibraryId, ordinal: u64) -> Vec<u8> {
        let mut bytes =
            format!("@{}\n", compose_template_id(library_id, ordinal).as_str()).into_bytes();
        bytes.extend(std::iter::repeat_n(
            b'A',
            usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH),
        ));
        bytes.extend(b"\n+\n");
        bytes.extend(std::iter::repeat_n(
            b'I',
            usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH),
        ));
        bytes.push(b'\n');
        bytes
    }

    fn write_case(root: &Path, name: &str, bytes: &[u8]) -> PathBuf {
        let temporary = root.join(format!("{name}.fastq.gz"));
        let file = fs::File::create(&temporary).expect("temporary FASTQ");
        let mut writer = bgzf::io::Writer::new(file);
        use std::io::Write;
        writer.write_all(bytes).expect("compress temporary FASTQ");
        writer.finish().expect("finish temporary FASTQ gzip");
        temporary
    }

    #[test]
    fn temporary_fastq_validation_rejects_each_contract_violation() {
        let root = temporary("matrix");
        fs::create_dir(&root).expect("root");
        let library_id = LibraryId::parse("library-fastq-test", "test").expect("library ID");
        let valid = valid_record(&library_id, 0);
        let header_length = format!("@{}\n", compose_template_id(&library_id, 0).as_str()).len();
        let sequence_start = header_length;
        let plus_start = sequence_start + usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1;
        let quality_start = plus_start + 2;

        let mut invalid_header = valid.clone();
        invalid_header[0] = b'X';
        let mut invalid_id = valid.clone();
        invalid_id[header_length - 2] = b'1';
        let mut invalid_separator = valid.clone();
        invalid_separator[plus_start] = b'-';
        let mut invalid_sequence_length = valid.clone();
        invalid_sequence_length.remove(sequence_start);
        let mut invalid_base = valid.clone();
        invalid_base[sequence_start] = b'X';
        let mut invalid_quality = valid.clone();
        invalid_quality[quality_start] = b'!';
        let mut short_quality = valid.clone();
        short_quality.remove(short_quality.len() - 2);
        let mut extra_eof = valid.clone();
        extra_eof.push(b'x');
        let mut truncated_quality_in_second_record = valid.clone();
        truncated_quality_in_second_record.extend(valid_record(&library_id, 1));
        truncated_quality_in_second_record.pop();

        let cases = vec![
            InvalidFastqCase {
                name: "invalid-header",
                bytes: invalid_header,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ header is not the canonical technical template ID",
            },
            InvalidFastqCase {
                name: "invalid-id",
                bytes: invalid_id,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ header is not the canonical technical template ID",
            },
            InvalidFastqCase {
                name: "truncated-quality",
                bytes: truncated_quality_in_second_record,
                expected_count: 2,
                ordinal: 1,
                detail: "FASTQ quality line is truncated",
            },
            InvalidFastqCase {
                name: "invalid-separator",
                bytes: invalid_separator,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ separator line is invalid",
            },
            InvalidFastqCase {
                name: "invalid-sequence-length",
                bytes: invalid_sequence_length,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ sequence length does not match the configured read length",
            },
            InvalidFastqCase {
                name: "invalid-base",
                bytes: invalid_base,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ base or Phred quality is not valid for illumina-se-minimal-v1",
            },
            InvalidFastqCase {
                name: "invalid-quality",
                bytes: invalid_quality,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ base or Phred quality is not valid for illumina-se-minimal-v1",
            },
            InvalidFastqCase {
                name: "short-quality",
                bytes: short_quality,
                expected_count: 1,
                ordinal: 0,
                detail: "FASTQ quality length does not match its sequence",
            },
            InvalidFastqCase {
                name: "extra-eof",
                bytes: extra_eof,
                expected_count: 1,
                ordinal: 1,
                detail: "FASTQ contains records beyond the planned count",
            },
        ];

        for case in cases {
            let temporary = write_case(&root, case.name, &case.bytes);
            let error = validate_temporary_fastq(
                &temporary,
                Path::new("outputs/fastq/event/sample.fastq.gz"),
                &library_id,
                case.expected_count,
            )
            .expect_err(case.name);
            assert_eq!(
                error.context().get("path"),
                Some(&"outputs/fastq/event/sample.fastq.gz".to_owned()),
                "{}",
                case.name
            );
            assert_eq!(
                error.context().get("record_ordinal"),
                Some(&case.ordinal.to_string()),
                "{}",
                case.name
            );
            assert_eq!(
                error.context().get("detail"),
                Some(&case.detail.to_owned()),
                "{}",
                case.name
            );
        }
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn workspace_rejects_invalid_fastq_and_removes_its_temporary_file() {
        let root = temporary("workspace-cleanup");
        fs::create_dir(&root).expect("root");
        let destination = root.join("published");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let library_id = LibraryId::parse("library-fastq-test", "test").expect("library ID");
        let relative_path = Path::new("outputs/fastq/event/sample.fastq.gz");
        let invalid = b"@not-the-canonical-template\n";

        let error = workspace
            .write_streamed_artifact_with_validation(
                "fastq-library-fastq-test",
                "fastq",
                relative_path,
                |file| {
                    use std::io::Write;
                    let mut writer = bgzf::io::Writer::new(file);
                    writer.write_all(invalid).map_err(|source| {
                        super::io_error(relative_path, "write_test_fastq", source)
                    })?;
                    writer.finish().map_err(|source| {
                        super::io_error(relative_path, "finish_test_fastq_gzip", source)
                    })?;
                    Ok(())
                },
                |temporary| validate_temporary_fastq(temporary, relative_path, &library_id, 1),
            )
            .expect_err("invalid FASTQ must be rejected");

        assert_eq!(
            error.context().get("detail"),
            Some(&"FASTQ header is not the canonical technical template ID".to_owned())
        );
        assert!(workspace.artifacts().is_empty());
        let artifact_parent = workspace.root().join("outputs/fastq/event");
        assert!(artifact_parent.is_dir());
        assert!(
            fs::read_dir(&artifact_parent)
                .expect("temporary artifact directory")
                .next()
                .is_none(),
            "the rejected temporary FASTQ must be removed"
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn paired_publication_rolls_back_on_injected_r2_or_joint_validation_failure() {
        for (label, inject_r2) in [("r2", true), ("joint", false)] {
            let root = temporary(label);
            fs::create_dir(&root).expect("root");
            let destination = root.join("published");
            let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
            let library_id = LibraryId::parse("library-final-test", "test").expect("library ID");
            if inject_r2 {
                FORCE_FINAL_R2_FAILURE.set(true);
            } else {
                FORCE_FINAL_PAIR_VALIDATION_FAILURE.set(true);
            }
            let result = write_final_paired_fastq(
                &mut workspace,
                "event-final-test",
                "sample-final-test",
                &library_id,
                Vec::<DemultiplexedPair<ObservedPair>>::new(),
            );
            assert!(result.is_err(), "{label} failure must propagate");
            assert!(
                workspace.artifacts().is_empty(),
                "{label} must not catalog a mate"
            );
            let r1 = workspace
                .root()
                .join("outputs/fastq/event-final-test/sample-final-test_R1.fastq.gz");
            let r2 = workspace
                .root()
                .join("outputs/fastq/event-final-test/sample-final-test_R2.fastq.gz");
            assert!(
                !r1.exists() && !r2.exists(),
                "{label} must not leave a mate"
            );
            fs::remove_dir_all(root).expect("cleanup");
        }
    }
}
