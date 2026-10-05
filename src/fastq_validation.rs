//! Validación independiente por flujo de FASTQ WGS mínimos.
//!
//! Esta frontera no reconstruye procedencia durable ni repara artefactos. Vuelve
//! a derivar cada origen desde la subsemilla de su plantilla y compara el
//! registro FASTQ contra esa vista de solo lectura del estado haplotípico.

use crate::guard_trace::rejected;
use std::fs::File;
use std::io::{BufRead, Read, Seek, SeekFrom};
use std::path::Path;

use noodles_bgzf as bgzf;
use sha2::{Digest, Sha256};

use crate::application::{
    EffectiveTechnicalEventTarget, ILLUMINA_SE_MINIMAL_ACGT_PHRED, ILLUMINA_SE_MINIMAL_N_PHRED,
    ILLUMINA_SE_MINIMAL_READ_LENGTH,
};
use crate::fastq_output::{FASTQ_ARTIFACT_KIND, fastq_relative_path};
use crate::haplotype::{HaplotypeCopy, Individual};
use crate::identity::{EventId, LibraryId, TechnicalPairId};
use crate::materialization::{MaterializationError, materialize_copy_interval};
use crate::randomness::{RandomSamplingError, poisson_count, uniform_u64_below};
use crate::reference::Reference;
use crate::results::ArtifactRecord;
use crate::technical_identity::TechnicalIdentityPlan;
use crate::technical_identity::TemplateOrdinal;

const WGS_READ_LENGTH: u64 = ILLUMINA_SE_MINIMAL_READ_LENGTH as u64;
/// § 4.1.2 de la especificación BGZF (2020-12-03). El lector BGZF acepta un
/// final truncado como EOF, por lo que la frontera de validación debe exigir
/// físicamente este marcador canónico completo.
const BGZF_EOF_MARKER: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0x00, 0x00, 0x00, 0x00, 0x00, 0xff, 0x06, 0x00, 0x42, 0x43, 0x02, 0x00,
    0x1b, 0x00, 0x03, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00,
];
const FINAL_FASTQ_MAXIMUM_LINE_LENGTH: usize = 256;

/// Desviación observada de cobertura WGS. Es una observación estadística y no
/// una política de aceptación: una desviación nunca convierte por sí sola el
/// resultado en un error.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WgsCoverageObservation {
    target_coverage: Option<f64>,
    realized_coverage: f64,
    absolute_deviation: Option<f64>,
    relative_deviation: Option<f64>,
}

/// Las tres medidas de cobertura F4 de una biblioteca.
///
/// El objetivo solo se compara con la WGS básica previa a efectos técnicos.
/// La cobertura FASTQ final se informa separadamente y, cuando F4.4 todavía
/// no ha sustituido parejas, todas sus bases son propias y ninguna es
/// contaminante. `n_bases` cuenta las `N` observadas en el FASTQ final, por lo
/// que puede superar el denominador estructural a coberturas mayores que una.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct F4CoverageReport {
    target_x: f64,
    basic_wgs_x: f64,
    final_fastq_own_x: f64,
    final_fastq_contaminant_x: f64,
    denominator_bases: u64,
    n_bases: u64,
}

impl F4CoverageReport {
    /// Construye el resumen desde contadores ya validados por el productor o
    /// por la relectura independiente.
    pub fn from_base_counts(
        target_x: f64,
        basic_wgs_bases: u64,
        final_fastq_own_bases: u64,
        final_fastq_contaminant_bases: u64,
        denominator_bases: u64,
        n_bases: u64,
    ) -> Result<Self, FastqValidationError> {
        if !target_x.is_finite() || target_x < 0.0 {
            return Err(rejected(FastqValidationError::InvalidCoverageSummary));
        }
        let final_fastq_bases = final_fastq_own_bases
            .checked_add(final_fastq_contaminant_bases)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
        if denominator_bases == 0 || n_bases > final_fastq_bases {
            return Err(rejected(FastqValidationError::InvalidCoverageSummary));
        }
        let denominator = denominator_bases as f64;
        Ok(Self {
            target_x,
            basic_wgs_x: (basic_wgs_bases as f64) / denominator,
            final_fastq_own_x: (final_fastq_own_bases as f64) / denominator,
            final_fastq_contaminant_x: (final_fastq_contaminant_bases as f64) / denominator,
            denominator_bases,
            n_bases,
        })
    }

    #[must_use]
    pub const fn target_x(self) -> f64 {
        self.target_x
    }

    #[must_use]
    pub const fn basic_wgs_x(self) -> f64 {
        self.basic_wgs_x
    }

    #[must_use]
    pub const fn final_fastq_own_x(self) -> f64 {
        self.final_fastq_own_x
    }

    #[must_use]
    pub const fn final_fastq_contaminant_x(self) -> f64 {
        self.final_fastq_contaminant_x
    }

    #[must_use]
    pub const fn denominator_bases(self) -> u64 {
        self.denominator_bases
    }

    #[must_use]
    pub const fn n_bases(self) -> u64 {
        self.n_bases
    }
}

impl WgsCoverageObservation {
    #[must_use]
    pub const fn target_coverage(self) -> Option<f64> {
        self.target_coverage
    }

    #[must_use]
    pub const fn realized_coverage(self) -> f64 {
        self.realized_coverage
    }

    #[must_use]
    pub const fn absolute_deviation(self) -> Option<f64> {
        self.absolute_deviation
    }

    #[must_use]
    pub const fn relative_deviation(self) -> Option<f64> {
        self.relative_deviation
    }
}

/// Resultado agregado de una validación FASTQ. No retiene lecturas ni origen
/// individual alguno después de comprobar cada registro.
#[derive(Clone, Debug, PartialEq)]
pub struct FastqValidationReport {
    fastq_sha256: String,
    read_count: u64,
    base_count: u64,
    unknown_base_count: u64,
    verified_origin_count: u64,
    coverage: WgsCoverageObservation,
    f4_coverage: F4CoverageReport,
}

impl FastqValidationReport {
    #[must_use]
    pub fn fastq_sha256(&self) -> &str {
        &self.fastq_sha256
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
    pub const fn unknown_base_count(&self) -> u64 {
        self.unknown_base_count
    }

    #[must_use]
    pub const fn verified_origin_count(&self) -> u64 {
        self.verified_origin_count
    }

    #[must_use]
    pub const fn coverage(&self) -> WgsCoverageObservation {
        self.coverage
    }

    /// Medidas F4 sin equiparar la cobertura final con el objetivo.
    #[must_use]
    pub const fn f4_coverage(&self) -> F4CoverageReport {
        self.f4_coverage
    }
}

/// Error de una inconsistencia determinista del FASTQ o de la reproducción de
/// su verdad. Ninguna variante representa una desviación estocástica de
/// cobertura.
#[derive(Debug)]
pub enum FastqValidationError {
    Io {
        path: String,
        operation: &'static str,
        detail: String,
    },
    Replay(FastqReplayError),
    Artifact {
        path: String,
        field: &'static str,
        expected: String,
        observed: String,
    },
    Structural {
        path: String,
        record_ordinal: u64,
        detail: &'static str,
    },
    CounterOverflow,
    InvalidCoverageSummary,
}

/// Comprueba una sola salida PE descomprimida, sin consultar verdad biológica
/// ni claves transitorias. La comparación contra esa verdad pertenece a F4.5.3.
pub fn validate_final_fastq_mate(
    path: &Path,
    expected_mate: crate::technical_identity::TechnicalMate,
) -> Result<(), FastqValidationError> {
    let file = File::open(path).map_err(|error| io_error(path, "open_final_fastq", error))?;
    validate_canonical_bgzf_eof(
        &mut file
            .try_clone()
            .map_err(|error| io_error(path, "clone_final_fastq", error))?,
        path,
        0,
    )?;
    let mut reader = bgzf::io::Reader::new(file);
    let mut lines = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    let mut ordinal = 0_u64;
    loop {
        lines
            .iter_mut()
            .try_for_each(|line| read_final_line(&mut reader, line, path, ordinal))?;
        if lines[0].is_empty() && lines[1].is_empty() && lines[2].is_empty() && lines[3].is_empty()
        {
            return Ok(());
        }
        if lines.iter().any(Vec::is_empty) {
            return Err(rejected(FastqValidationError::Structural {
                path: path.display().to_string(),
                record_ordinal: ordinal,
                detail: "FASTQ final has a truncated record",
            }));
        }
        validate_final_record(path, ordinal, &lines, expected_mate)?;
        ordinal = ordinal
            .checked_add(1)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
    }
}

/// Comprueba que R1/R2 contienen exactamente las mismas parejas en el mismo
/// orden y que sus hashes corresponden a los bytes comprimidos publicados.
pub fn validate_final_paired_fastq(
    r1_path: &Path,
    r1_artifact: &ArtifactRecord,
    r2_path: &Path,
    r2_artifact: &ArtifactRecord,
) -> Result<(), FastqValidationError> {
    validate_final_fastq_mate(r1_path, crate::technical_identity::TechnicalMate::R1)?;
    validate_final_fastq_mate(r2_path, crate::technical_identity::TechnicalMate::R2)?;
    validate_final_artifact(r1_path, r1_artifact)?;
    validate_final_artifact(r2_path, r2_artifact)?;
    let mut r1 = bgzf::io::Reader::new(
        File::open(r1_path).map_err(|error| io_error(r1_path, "open_final_r1", error))?,
    );
    let mut r2 = bgzf::io::Reader::new(
        File::open(r2_path).map_err(|error| io_error(r2_path, "open_final_r2", error))?,
    );
    let mut ordinal = 0_u64;
    loop {
        let left = read_final_header(
            &mut r1,
            r1_path,
            ordinal,
            crate::technical_identity::TechnicalMate::R1,
        )?;
        let right = read_final_header(
            &mut r2,
            r2_path,
            ordinal,
            crate::technical_identity::TechnicalMate::R2,
        )?;
        if left.is_none() && right.is_none() {
            return Ok(());
        }
        if left.is_none() || right.is_none() || left != right {
            return Err(rejected(FastqValidationError::Structural {
                path: format!("{} / {}", r1_path.display(), r2_path.display()),
                record_ordinal: ordinal,
                detail: "FASTQ final mates are orphaned or have different pair IDs",
            }));
        }
        ordinal = ordinal
            .checked_add(1)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
    }
}

fn validate_final_artifact(
    path: &Path,
    artifact: &ArtifactRecord,
) -> Result<(), FastqValidationError> {
    let identity =
        crate::integrity::identify_file(path, Path::new(&artifact.path)).map_err(|error| {
            rejected(FastqValidationError::Artifact {
                path: path.display().to_string(),
                field: "sha256",
                expected: artifact.digest.clone(),
                observed: error.render_human(),
            })
        })?;
    if identity.size_bytes() != artifact.size_bytes {
        return Err(rejected(FastqValidationError::Artifact {
            path: path.display().to_string(),
            field: "size_bytes",
            expected: artifact.size_bytes.to_string(),
            observed: identity.size_bytes().to_string(),
        }));
    }
    if identity.digest() != artifact.digest {
        return Err(rejected(FastqValidationError::Artifact {
            path: path.display().to_string(),
            field: "sha256",
            expected: artifact.digest.clone(),
            observed: identity.digest().to_owned(),
        }));
    }
    Ok(())
}

fn read_final_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    path: &Path,
    ordinal: u64,
) -> Result<(), FastqValidationError> {
    line.clear();
    loop {
        let buffer = reader
            .fill_buf()
            .map_err(|error| io_error(path, "read_final_fastq", error))?;
        if buffer.is_empty() {
            break;
        }
        let take = buffer
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(buffer.len(), |index| index + 1);
        if line
            .len()
            .checked_add(take)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?
            > FINAL_FASTQ_MAXIMUM_LINE_LENGTH
        {
            return Err(rejected(FastqValidationError::Structural {
                path: path.display().to_string(),
                record_ordinal: ordinal,
                detail: "FASTQ final line exceeds the profile bound",
            }));
        }
        line.extend_from_slice(&buffer[..take]);
        let reached_newline = buffer[take - 1] == b'\n';
        reader.consume(take);
        if reached_newline {
            break;
        }
    }
    if !line.is_empty() && line.last() != Some(&b'\n') {
        return Err(rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: "FASTQ final line is truncated",
        }));
    }
    Ok(())
}

fn validate_final_record(
    path: &Path,
    ordinal: u64,
    lines: &[Vec<u8>; 4],
    expected_mate: crate::technical_identity::TechnicalMate,
) -> Result<(), FastqValidationError> {
    let structural = |detail| FastqValidationError::Structural {
        path: path.display().to_string(),
        record_ordinal: ordinal,
        detail,
    };
    let header = std::str::from_utf8(&lines[0])
        .map_err(|_| rejected(structural("FASTQ final header is not UTF-8")))?
        .strip_suffix('\n')
        .ok_or_else(|| rejected(structural("FASTQ final header is truncated")))?;
    let Some(header) = header.strip_prefix('@') else {
        return Err(rejected(structural("FASTQ final header lacks @")));
    };
    parse_final_header(header, expected_mate).map_err(structural)?;
    if lines[2] != b"+\n"
        || lines[1].len() < 2
        || lines[1].last() != Some(&b'\n')
        || lines[3].len() != lines[1].len()
        || lines[3].last() != Some(&b'\n')
    {
        return Err(rejected(structural(
            "FASTQ final record structure is invalid",
        )));
    }
    if !lines[1][..lines[1].len() - 1]
        .iter()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
        || !lines[3][..lines[3].len() - 1]
            .iter()
            .all(|quality| (b'!'..=b'~').contains(quality))
    {
        return Err(rejected(structural(
            "FASTQ final sequence or quality alphabet is invalid",
        )));
    }
    Ok(())
}

fn read_final_header(
    reader: &mut impl BufRead,
    path: &Path,
    ordinal: u64,
    expected_mate: crate::technical_identity::TechnicalMate,
) -> Result<Option<FinalFastqHeader>, FastqValidationError> {
    let mut lines = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];
    lines
        .iter_mut()
        .try_for_each(|line| read_final_line(reader, line, path, ordinal))?;
    if lines.iter().all(Vec::is_empty) {
        return Ok(None);
    }
    if lines.iter().any(Vec::is_empty) {
        return Err(rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: "FASTQ final has a truncated record",
        }));
    }
    let header = std::str::from_utf8(&lines[0]).map_err(|_| {
        rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: "FASTQ final header is not UTF-8",
        })
    })?;
    let header = header.strip_suffix('\n').ok_or_else(|| {
        rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: "FASTQ final header is truncated",
        })
    })?;
    let header = header.strip_prefix('@').ok_or_else(|| {
        rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: "FASTQ final header lacks @",
        })
    })?;
    parse_final_header(header, expected_mate)
        .map(Some)
        .map_err(|detail| FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail,
        })
}

/// Identidad publicable de una pareja FASTQ final. No contiene ninguna clave
/// previa al salto ni procedencia por lectura.
#[derive(Debug, Eq, PartialEq)]
struct FinalFastqHeader {
    event_id: String,
    library_id: String,
    technical_pair_id: String,
    udi: String,
}

fn parse_final_header(
    header: &str,
    expected_mate: crate::technical_identity::TechnicalMate,
) -> Result<FinalFastqHeader, &'static str> {
    let mut fields = header.split(' ');
    let identity = fields
        .next()
        .ok_or_else(|| rejected("FASTQ final header lacks identity"))?;
    let mate = fields
        .next()
        .ok_or_else(|| rejected("FASTQ final header lacks mate metadata"))?;
    if fields.next().is_some() {
        return Err(rejected("FASTQ final header has extra fields"));
    }
    let mut identity = identity.split(':');
    if identity.next() != Some("GPS") {
        return Err(rejected("FASTQ final header lacks GPS prefix"));
    }
    let event_id = identity
        .next()
        .ok_or_else(|| rejected("FASTQ final header lacks event ID"))?;
    let library_id = identity
        .next()
        .ok_or_else(|| rejected("FASTQ final header lacks library ID"))?;
    let technical_pair_id = identity
        .next()
        .ok_or_else(|| rejected("FASTQ final header lacks technical pair ID"))?;
    if identity.next().is_some() {
        return Err(rejected("FASTQ final header has invalid identity fields"));
    }
    EventId::parse(event_id, "final_fastq")
        .map_err(|_| rejected("FASTQ final event ID is invalid"))?;
    LibraryId::parse(library_id, "final_fastq")
        .map_err(|_| rejected("FASTQ final library ID is invalid"))?;
    TechnicalPairId::parse(technical_pair_id, "final_fastq")
        .map_err(|_| rejected("FASTQ final pair ID is invalid"))?;
    let expected_prefix = format!("{}:N:0:", expected_mate.number());
    let udi = mate
        .strip_prefix(&expected_prefix)
        .ok_or_else(|| rejected("FASTQ final mate metadata is invalid"))?;
    let mut udi_components = udi.split('+');
    if !udi_components.next().is_some_and(valid_udi_component)
        || !udi_components.next().is_some_and(valid_udi_component)
        || udi_components.next().is_some()
    {
        return Err(rejected("FASTQ final UDI is invalid"));
    }
    Ok(FinalFastqHeader {
        event_id: event_id.to_owned(),
        library_id: library_id.to_owned(),
        technical_pair_id: technical_pair_id.to_owned(),
        udi: udi.to_owned(),
    })
}

fn valid_udi_component(value: &str) -> bool {
    value.len() == 8
        && value
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

#[cfg(test)]
mod final_header_tests {
    use super::parse_final_header;
    use crate::technical_identity::TechnicalMate;

    const R1: &str =
        "GPS:event-main:library-main:TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 1:N:0:ACGTACGT+TGCATGCA";
    const R2: &str =
        "GPS:event-main:library-main:TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 2:N:0:ACGTACGT+TGCATGCA";

    #[test]
    fn final_headers_require_gps_identity_and_matching_pair_components() {
        let r1 = parse_final_header(R1, TechnicalMate::R1).expect("canonical R1");
        let r2 = parse_final_header(R2, TechnicalMate::R2).expect("canonical R2");
        assert_eq!(r1, r2);
        for mutated in [
            "GPS:event-other:library-main:TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 2:N:0:ACGTACGT+TGCATGCA",
            "GPS:event-main:library-other:TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 2:N:0:ACGTACGT+TGCATGCA",
            "GPS:event-main:library-main:TP_ZYXWVUTSRQPONMLKJIHGFEDCBA 2:N:0:ACGTACGT+TGCATGCA",
            "GPS:event-main:library-main:TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 2:N:0:TGCATGCA+ACGTACGT",
        ] {
            let parsed =
                parse_final_header(mutated, TechnicalMate::R2).expect("well-formed mutation");
            assert_ne!(
                r1, parsed,
                "independent identity/UDI mutation must break the pair"
            );
        }
        assert!(
            parse_final_header(
                "TP_ABCDEFGHIJKLMNOPQRSTUVWXYZ 1:N:0:ACGTACGT+TGCATGCA",
                TechnicalMate::R1,
            )
            .is_err()
        );
    }
}

/// Error al reproducir localmente las decisiones que originaron una biblioteca.
///
/// No reutiliza el iterador de producción: esta frontera deriva las mismas
/// causas estables (identidad, subsemilla, contratos RNG y materialización)
/// mediante un recorrido propio para que pueda detectar regresiones comunes.
#[derive(Debug)]
pub enum FastqReplayError {
    TechnicalPlanReferenceMismatch,
    IndividualReferenceMismatch,
    MissingTechnicalLibrary,
    NoEligibleSequence { read_length: u64 },
    LengthOverflow,
    InvalidCoverageMean,
    Random(RandomSamplingError),
    TechnicalIdentity(crate::diagnostics::Diagnostic),
    Materialization(MaterializationError),
}

impl std::fmt::Display for FastqValidationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for FastqValidationError {}

/// Relee y valida un FASTQ final frente al estado faseado y las subsemillas.
///
/// El lector usa memoria acotada: cuatro búferes de una línea FASTQ y los
/// candidatos de copia del recorrido de reproducción, sin acumular registros
/// ni procedencia.
pub fn validate_wgs_fastq(
    fastq_path: &Path,
    expected_artifact: &ArtifactRecord,
    reference: &Reference,
    individual: &Individual,
    identities: &TechnicalIdentityPlan,
    target: EffectiveTechnicalEventTarget,
) -> Result<FastqValidationReport, FastqValidationError> {
    let replay = ReplayedWgsPlan::new(reference, individual, identities, target)
        .map_err(FastqValidationError::Replay)?;
    validate_declared_artifact(expected_artifact, identities, individual)?;
    let expected_count = replay.read_count;
    let mut file = File::open(fastq_path).map_err(|error| io_error(fastq_path, "open", error))?;
    validate_canonical_bgzf_eof(&mut file, fastq_path, expected_count)?;
    file.seek(SeekFrom::Start(0))
        .map_err(|error| io_error(fastq_path, "rewind_after_eof_check", error))?;
    // BGZF es un perfil gzip interoperable. Su lector verifica cada bloque
    // comprimido y expone el FASTQ descomprimido sin acumularlo ni envolverlo
    // en otro búfer, que su contrato prohíbe.
    let mut reader = bgzf::io::Reader::new(HashingReader::new(file));
    let mut header = Vec::with_capacity(64);
    let mut sequence = Vec::with_capacity(usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1);
    let mut plus = Vec::with_capacity(2);
    let mut qualities = Vec::with_capacity(usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1);
    let mut unknown_base_count = 0_u64;

    for ordinal in 0..expected_count {
        let expected = replay
            .origin(ordinal)
            .map_err(FastqValidationError::Replay)?;
        let header_maximum = expected
            .template_id
            .as_str()
            .len()
            .checked_add(2)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
        read_line(
            &mut reader,
            &mut header,
            fastq_path,
            ordinal,
            "header",
            header_maximum,
        )?;
        read_line(
            &mut reader,
            &mut sequence,
            fastq_path,
            ordinal,
            "sequence",
            usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1,
        )?;
        read_line(&mut reader, &mut plus, fastq_path, ordinal, "separator", 2)?;
        read_line(
            &mut reader,
            &mut qualities,
            fastq_path,
            ordinal,
            "quality",
            usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1,
        )?;
        unknown_base_count = validate_record(
            fastq_path,
            ordinal,
            &header,
            &sequence,
            &plus,
            &qualities,
            &expected,
            unknown_base_count,
        )?;
    }
    let mut trailing = [0_u8; 1];
    if reader
        .read(&mut trailing)
        .map_err(|error| io_error(fastq_path, "check_eof", error))?
        != 0
    {
        return Err(rejected(FastqValidationError::Structural {
            path: fastq_path.display().to_string(),
            record_ordinal: expected_count,
            detail: "FASTQ contains records beyond the reproduced read plan",
        }));
    }

    // La identidad se actualiza desde el mismo descriptor que alimentó al
    // lector BGZF y se cierra solo después de alcanzar su EOF. Así no puede
    // asociarse la semántica de una versión con el SHA-256 de otra apertura.
    let (observed_size_bytes, observed_digest) = reader.into_inner().finish();

    let base_count = expected_count
        .checked_mul(WGS_READ_LENGTH)
        .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
    let target_coverage = match target {
        EffectiveTechnicalEventTarget::CoverageWgs(value) => Some(value),
        EffectiveTechnicalEventTarget::ReadCount(_) => None,
    };
    let structural_denominator_bases = reference.total_length();
    let f4_coverage = F4CoverageReport::from_base_counts(
        target_coverage.unwrap_or(0.0),
        base_count,
        base_count,
        0,
        structural_denominator_bases,
        unknown_base_count,
    )?;
    // Esta observación histórica conserva su API, pero su valor realizado es
    // la WGS básica: es la única medida que se compara con el objetivo.
    let realized_coverage = f4_coverage.basic_wgs_x();
    let absolute_deviation = target_coverage.map(|value| (realized_coverage - value).abs());
    let relative_deviation = target_coverage
        .and_then(|value| (value != 0.0).then(|| (realized_coverage - value).abs() / value.abs()));
    validate_observed_artifact(
        expected_artifact,
        fastq_path,
        observed_size_bytes,
        &observed_digest,
    )?;
    Ok(FastqValidationReport {
        fastq_sha256: observed_digest,
        read_count: expected_count,
        base_count,
        unknown_base_count,
        verified_origin_count: expected_count,
        coverage: WgsCoverageObservation {
            target_coverage,
            realized_coverage,
            absolute_deviation,
            relative_deviation,
        },
        f4_coverage,
    })
}

fn read_line(
    reader: &mut impl BufRead,
    line: &mut Vec<u8>,
    path: &Path,
    ordinal: u64,
    field: &'static str,
    maximum_length: usize,
) -> Result<(), FastqValidationError> {
    line.clear();
    loop {
        let buffer = reader
            .fill_buf()
            .map_err(|error| io_error(path, "read_record", error))?;
        if buffer.is_empty() {
            break;
        }
        let take = buffer
            .iter()
            .position(|byte| *byte == b'\n')
            .map_or(buffer.len(), |index| index + 1);
        let new_length = line
            .len()
            .checked_add(take)
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
        if new_length > maximum_length {
            return Err(rejected(FastqValidationError::Structural {
                path: path.display().to_string(),
                record_ordinal: ordinal,
                detail: "FASTQ line exceeds the minimal profile bound",
            }));
        }
        line.extend_from_slice(&buffer[..take]);
        let reached_newline = buffer[take - 1] == b'\n';
        reader.consume(take);
        if reached_newline {
            break;
        }
    }
    if line.is_empty() || line.last() != Some(&b'\n') {
        return Err(rejected(FastqValidationError::Structural {
            path: path.display().to_string(),
            record_ordinal: ordinal,
            detail: match field {
                "header" => "FASTQ header line is truncated",
                "sequence" => "FASTQ sequence line is truncated",
                "separator" => "FASTQ separator line is truncated",
                "quality" => "FASTQ quality line is truncated",
                _ => "FASTQ line is truncated",
            },
        }));
    }
    Ok(())
}

#[derive(Default)]
struct StreamDigest {
    hasher: Sha256,
    size_bytes: u64,
}

struct HashingReader<R> {
    inner: R,
    digest: StreamDigest,
}

impl<R> HashingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            digest: StreamDigest::default(),
        }
    }

    fn finish(self) -> (u64, String) {
        self.digest.finish()
    }
}

impl<R: Read> Read for HashingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.digest.update(&buffer[..read]).map_err(|error| {
            std::io::Error::new(std::io::ErrorKind::InvalidData, error.to_string())
        })?;
        Ok(read)
    }
}

fn validate_canonical_bgzf_eof(
    file: &mut File,
    path: &Path,
    expected_count: u64,
) -> Result<(), FastqValidationError> {
    let size = file
        .metadata()
        .map_err(|error| io_error(path, "metadata_eof_marker", error))?
        .len();
    let marker_length = u64::try_from(BGZF_EOF_MARKER.len())
        .map_err(|_| rejected(FastqValidationError::CounterOverflow))?;
    if size < marker_length {
        return Err(missing_canonical_eof(path, expected_count));
    }
    file.seek(SeekFrom::End(
        -i64::try_from(marker_length)
            .map_err(|_| rejected(FastqValidationError::CounterOverflow))?,
    ))
    .map_err(|error| io_error(path, "seek_eof_marker", error))?;
    let mut observed = [0_u8; BGZF_EOF_MARKER.len()];
    file.read_exact(&mut observed)
        .map_err(|error| io_error(path, "read_eof_marker", error))?;
    if observed != BGZF_EOF_MARKER {
        return Err(missing_canonical_eof(path, expected_count));
    }
    Ok(())
}

#[track_caller]
fn missing_canonical_eof(path: &Path, record_ordinal: u64) -> FastqValidationError {
    crate::guard_trace::record();
    FastqValidationError::Structural {
        path: path.display().to_string(),
        record_ordinal,
        detail: "FASTQ BGZF stream is missing its canonical EOF marker",
    }
}

impl StreamDigest {
    fn update(&mut self, bytes: &[u8]) -> Result<(), FastqValidationError> {
        self.size_bytes = self
            .size_bytes
            .checked_add(
                u64::try_from(bytes.len())
                    .map_err(|_| rejected(FastqValidationError::CounterOverflow))?,
            )
            .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
        self.hasher.update(bytes);
        Ok(())
    }

    fn finish(self) -> (u64, String) {
        (self.size_bytes, format!("{:x}", self.hasher.finalize()))
    }
}

#[allow(clippy::too_many_arguments)]
fn validate_record(
    path: &Path,
    ordinal: u64,
    header: &[u8],
    sequence: &[u8],
    plus: &[u8],
    qualities: &[u8],
    expected: &ReplayedWgsOrigin,
    unknown_base_count: u64,
) -> Result<u64, FastqValidationError> {
    let structural = |detail| FastqValidationError::Structural {
        path: path.display().to_string(),
        record_ordinal: ordinal,
        detail,
    };
    let mut expected_header = Vec::with_capacity(expected.template_id.as_str().len() + 2);
    expected_header.push(b'@');
    expected_header.extend_from_slice(expected.template_id.as_str().as_bytes());
    expected_header.push(b'\n');
    if header != expected_header {
        return Err(rejected(structural(
            "FASTQ template ID does not match the reproduced origin",
        )));
    }
    let line_length = usize::from(ILLUMINA_SE_MINIMAL_READ_LENGTH) + 1;
    if sequence.len() != line_length || sequence.last() != Some(&b'\n') {
        return Err(rejected(structural(
            "FASTQ sequence length does not match the minimal profile",
        )));
    }
    if plus != b"+\n" {
        return Err(rejected(structural("FASTQ separator line is invalid")));
    }
    if qualities.len() != line_length || qualities.last() != Some(&b'\n') {
        return Err(rejected(structural(
            "FASTQ quality length does not match the minimal profile",
        )));
    }
    if sequence[..sequence.len() - 1] != *expected.true_sequence {
        return Err(rejected(structural(
            "FASTQ sequence does not match the reproduced haplotype origin",
        )));
    }
    let mut unknowns = unknown_base_count;
    for (base, quality) in sequence[..sequence.len() - 1]
        .iter()
        .copied()
        .zip(qualities[..qualities.len() - 1].iter().copied())
    {
        let expected_quality = match base {
            b'A' | b'C' | b'G' | b'T' => ILLUMINA_SE_MINIMAL_ACGT_PHRED + b'!',
            b'N' => {
                unknowns = unknowns
                    .checked_add(1)
                    .ok_or_else(|| rejected(FastqValidationError::CounterOverflow))?;
                ILLUMINA_SE_MINIMAL_N_PHRED + b'!'
            }
            _ => {
                return Err(rejected(structural(
                    "FASTQ sequence contains a non-ACGTN base",
                )));
            }
        };
        if quality != expected_quality {
            return Err(rejected(structural(
                "FASTQ Phred quality does not match its literal base",
            )));
        }
    }
    Ok(unknowns)
}

struct ReplayedEligibleCopy<'a> {
    copy: &'a HaplotypeCopy,
    start_count: u64,
}

struct ReplayedWgsPlan<'a> {
    reference: &'a Reference,
    individual: &'a Individual,
    identities: &'a TechnicalIdentityPlan,
    eligible_copies: Vec<ReplayedEligibleCopy<'a>>,
    total_start_count: u64,
    read_count: u64,
}

struct ReplayedWgsOrigin {
    template_id: crate::identity::TemplateId,
    true_sequence: Vec<u8>,
}

impl<'a> ReplayedWgsPlan<'a> {
    fn new(
        reference: &'a Reference,
        individual: &'a Individual,
        identities: &'a TechnicalIdentityPlan,
        target: EffectiveTechnicalEventTarget,
    ) -> Result<Self, FastqReplayError> {
        if identities.reference_id() != reference.id() {
            return Err(rejected(FastqReplayError::TechnicalPlanReferenceMismatch));
        }
        if individual.reference_id() != reference.id() {
            return Err(rejected(FastqReplayError::IndividualReferenceMismatch));
        }
        if identities.library(individual.id()).is_none() {
            return Err(rejected(FastqReplayError::MissingTechnicalLibrary));
        }
        let haploid_eligible_length = reference
            .sequences()
            .iter()
            .filter(|sequence| sequence.length() >= WGS_READ_LENGTH)
            .try_fold(0_u64, |total, sequence| {
                total
                    .checked_add(sequence.length())
                    .ok_or_else(|| rejected(FastqReplayError::LengthOverflow))
            })?;
        let mut eligible_copies = Vec::new();
        let mut total_start_count = 0_u64;
        for genotype in individual.sequences() {
            let reference_sequence = reference
                .sequences()
                .iter()
                .find(|sequence| sequence.id() == genotype.sequence_id())
                .ok_or_else(|| rejected(FastqReplayError::IndividualReferenceMismatch))?;
            if reference_sequence.length() < WGS_READ_LENGTH {
                continue;
            }
            for copy in genotype.copies() {
                let current_length = copy
                    .segments()
                    .current_length()
                    .map_err(|_| rejected(FastqReplayError::LengthOverflow))?;
                let Some(start_count) = current_length
                    .checked_sub(WGS_READ_LENGTH)
                    .and_then(|value| value.checked_add(1))
                else {
                    continue;
                };
                total_start_count = total_start_count
                    .checked_add(start_count)
                    .ok_or_else(|| rejected(FastqReplayError::LengthOverflow))?;
                eligible_copies.push(ReplayedEligibleCopy { copy, start_count });
            }
        }
        if total_start_count == 0 {
            return Err(rejected(FastqReplayError::NoEligibleSequence {
                read_length: WGS_READ_LENGTH,
            }));
        }
        let read_count = match target {
            EffectiveTechnicalEventTarget::ReadCount(read_count) => read_count,
            EffectiveTechnicalEventTarget::CoverageWgs(coverage) => {
                let lambda = coverage * (haploid_eligible_length as f64) / (WGS_READ_LENGTH as f64);
                if !lambda.is_finite() {
                    return Err(rejected(FastqReplayError::InvalidCoverageMean));
                }
                let count_seed = identities
                    .derive_library_count_subseed(individual.id())
                    .map_err(FastqReplayError::TechnicalIdentity)?;
                poisson_count(&mut count_seed.local_rng(), lambda)
                    .map_err(FastqReplayError::Random)?
            }
        };
        Ok(Self {
            reference,
            individual,
            identities,
            eligible_copies,
            total_start_count,
            read_count,
        })
    }

    fn origin(&self, ordinal: u64) -> Result<ReplayedWgsOrigin, FastqReplayError> {
        let template = self
            .identities
            .derive_template(self.individual.id(), TemplateOrdinal::new(ordinal))
            .map_err(FastqReplayError::TechnicalIdentity)?;
        let mut rng = template.subseed().local_rng();
        let selection = uniform_u64_below(&mut rng, self.total_start_count)
            .map_err(FastqReplayError::Random)?;
        let mut cursor = 0_u64;
        let mut candidate = None;
        for replayed in &self.eligible_copies {
            cursor = cursor
                .checked_add(replayed.start_count)
                .ok_or_else(|| rejected(FastqReplayError::LengthOverflow))?;
            if selection < cursor {
                candidate = Some(replayed);
                break;
            }
        }
        let candidate = candidate.ok_or_else(|| rejected(FastqReplayError::LengthOverflow))?;
        let start =
            uniform_u64_below(&mut rng, candidate.start_count).map_err(FastqReplayError::Random)?;
        let end = start
            .checked_add(WGS_READ_LENGTH)
            .ok_or_else(|| rejected(FastqReplayError::LengthOverflow))?;
        let mut true_sequence =
            materialize_copy_interval(candidate.copy, self.reference, start, end)
                .map_err(FastqReplayError::Materialization)?;
        if uniform_u64_below(&mut rng, 2).map_err(FastqReplayError::Random)? == 1 {
            reverse_complement_in_place(&mut true_sequence);
        }
        Ok(ReplayedWgsOrigin {
            template_id: template.template_id().clone(),
            true_sequence,
        })
    }
}

fn reverse_complement_in_place(sequence: &mut [u8]) {
    sequence.reverse();
    for base in sequence {
        *base = match *base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => unreachable!("the materialization contract only returns ACGTN"),
        };
    }
}

fn validate_declared_artifact(
    artifact: &ArtifactRecord,
    identities: &TechnicalIdentityPlan,
    individual: &Individual,
) -> Result<(), FastqValidationError> {
    let library = identities.library(individual.id()).ok_or_else(|| {
        rejected(FastqValidationError::Replay(
            FastqReplayError::MissingTechnicalLibrary,
        ))
    })?;
    let expected_id = format!("fastq-{}", library.library_id().as_str());
    let expected_path =
        fastq_relative_path(identities.event_id().as_str(), individual.id().as_str())
            .to_string_lossy()
            .replace('\\', "/");
    validate_artifact_field(artifact, "id", &expected_id, &artifact.id)?;
    validate_artifact_field(artifact, "kind", FASTQ_ARTIFACT_KIND, &artifact.kind)?;
    validate_artifact_field(artifact, "path", &expected_path, &artifact.path)
}

fn validate_observed_artifact(
    artifact: &ArtifactRecord,
    path: &Path,
    observed_size_bytes: u64,
    observed_digest: &str,
) -> Result<(), FastqValidationError> {
    if artifact.size_bytes != observed_size_bytes {
        return Err(rejected(FastqValidationError::Artifact {
            path: path.display().to_string(),
            field: "size_bytes",
            expected: artifact.size_bytes.to_string(),
            observed: observed_size_bytes.to_string(),
        }));
    }
    if artifact.digest != observed_digest {
        return Err(rejected(FastqValidationError::Artifact {
            path: path.display().to_string(),
            field: "digest",
            expected: artifact.digest.clone(),
            observed: observed_digest.to_owned(),
        }));
    }
    Ok(())
}

fn validate_artifact_field(
    artifact: &ArtifactRecord,
    field: &'static str,
    expected: &str,
    observed: &str,
) -> Result<(), FastqValidationError> {
    if expected == observed {
        return Ok(());
    }
    Err(rejected(FastqValidationError::Artifact {
        path: artifact.path.clone(),
        field,
        expected: expected.to_owned(),
        observed: observed.to_owned(),
    }))
}

fn io_error(path: &Path, operation: &'static str, error: std::io::Error) -> FastqValidationError {
    FastqValidationError::Io {
        path: path.display().to_string(),
        operation,
        detail: error.to_string(),
    }
}
