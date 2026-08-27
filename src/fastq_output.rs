//! Escritura FASTQ mínima por flujo desde orígenes WGS ya faseados.
//!
//! F3.1.4 solo convierte cada origen transitorio en su registro FASTQ final.
//! No conserva lecturas ni procedencia por lectura: al cerrar quedan el
//! catálogo de artefactos de `ExecutionWorkspace` y dos contadores agregados.

use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::application::{
    ILLUMINA_SE_MINIMAL_ACGT_PHRED, ILLUMINA_SE_MINIMAL_N_PHRED, ILLUMINA_SE_MINIMAL_READ_LENGTH,
};
use crate::diagnostics::Diagnostic;
use crate::identity::{LibraryId, compose_template_id};
use crate::publication::ExecutionWorkspace;
use crate::results::ArtifactRecord;
use crate::wgs_sampling::{WGS_READ_LENGTH, WgsLibrarySampler, WgsSamplingError};

pub const FASTQ_ARTIFACT_KIND: &str = "fastq";
const FASTQ_BUFFER_CAPACITY: usize = 64 * 1024;

/// Artefacto FASTQ final de la única biblioteca de un individuo en F3.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastqPublication {
    artifact: ArtifactRecord,
    read_count: u64,
    base_count: u64,
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
/// ciclo. El temporal se relee antes de calcular SHA-256 y renombrarlo como
/// artefacto del provisional, comprobando el formato, el contador y cada ID y
/// calidad canónicos.
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
            summary = write_records(file, sampler, expected_count, &path_for_write)?;
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
        .join(format!("{individual_id}.fastq"))
}

#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
struct FastqSummary {
    read_count: u64,
    base_count: u64,
}

fn write_records(
    file: &mut File,
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
    let file = File::open(temporary).map_err(|error| io_error(temporary, "open_fastq", error))?;
    let mut reader = BufReader::with_capacity(FASTQ_BUFFER_CAPACITY, file);
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

    use super::validate_temporary_fastq;
    use super::{ILLUMINA_SE_MINIMAL_READ_LENGTH, LibraryId, compose_template_id};
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
        let temporary = root.join(format!("{name}.fastq"));
        fs::write(&temporary, bytes).expect("temporary FASTQ");
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
                Path::new("outputs/fastq/event/sample.fastq"),
                &library_id,
                case.expected_count,
            )
            .expect_err(case.name);
            assert_eq!(
                error.context().get("path"),
                Some(&"outputs/fastq/event/sample.fastq".to_owned()),
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
        let relative_path = Path::new("outputs/fastq/event/sample.fastq");
        let invalid = b"@not-the-canonical-template\n";

        let error = workspace
            .write_streamed_artifact_with_validation(
                "fastq-library-fastq-test",
                "fastq",
                relative_path,
                |file| {
                    use std::io::Write;

                    file.write_all(invalid).map_err(|source| {
                        super::io_error(relative_path, "write_test_fastq", source)
                    })
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
}
