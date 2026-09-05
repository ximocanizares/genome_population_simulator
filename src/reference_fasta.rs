//! Adaptador FASTA síncrono tras el contrato propio de referencia.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Cursor, Read};
use std::num::NonZero;
use std::path::Path;

use noodles_fasta as fasta;

use crate::diagnostics::Diagnostic;
use crate::identity::{ReferenceId, SequenceId};
use crate::integrity::{ArtifactIdentity, Sha256Writer};
use crate::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceError, ReferenceIndex, ReferenceIndexKind,
    ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use sha2::{Digest, Sha256};

pub const CANONICAL_FASTA_VERSION: &str = "1";
pub const CANONICAL_LINE_BASE_COUNT: usize = 80;

/// Metadatos declarados que deben corresponder a una secuencia del FASTA.
#[derive(Clone, Copy, Debug)]
pub struct SequenceExpectation<'a> {
    pub id: &'a str,
    pub class: SequenceClass,
    pub length: u64,
}

/// FASTA suministrado validado junto con la identidad exacta observada.
#[derive(Debug)]
pub struct ValidatedSuppliedReference {
    reference: Reference,
    source_identity: ArtifactIdentity,
}

impl ValidatedSuppliedReference {
    /// Reúne una referencia ya construida con la identidad observada del FASTA
    /// del que procede. Permite entregar una referencia cuyo catálogo se validó
    /// sin conservar las bases.
    #[must_use]
    pub const fn from_parts(reference: Reference, source_identity: ArtifactIdentity) -> Self {
        Self {
            reference,
            source_identity,
        }
    }

    #[must_use]
    pub const fn reference(&self) -> &Reference {
        &self.reference
    }

    pub(crate) const fn reference_mut(&mut self) -> &mut Reference {
        &mut self.reference
    }

    #[must_use]
    pub const fn source_identity(&self) -> &ArtifactIdentity {
        &self.source_identity
    }

    #[must_use]
    pub fn into_reference(self) -> Reference {
        self.reference
    }
}

/// Valida por registros sin reescribir el origen y calcula su SHA-256 sobre los
/// bytes exactos leídos. El parser externo queda encapsulado en esta frontera;
/// las reglas de alfabeto, identidad y cabecera pertenecen al dominio propio.
pub fn validate_supplied_fasta(
    path: &Path,
    reference_id: &str,
    expected: &[SequenceExpectation<'_>],
) -> Result<ValidatedSuppliedReference, Diagnostic> {
    validate_supplied_fasta_with_origin(path, path, reference_id, expected)
}

/// Valida un FASTA suministrado cuando el consumidor solo necesita su orden
/// físico y sus bases, no un catálogo de secuencias declarado por separado.
///
/// Conserva las mismas comprobaciones de bytes, cabeceras y alfabeto que la
/// frontera con expectativas explícitas. F5 usa esta variante tras haber
/// ligado la identidad completa de la referencia en su preflight.
pub fn validate_supplied_fasta_unconstrained(
    path: &Path,
    reference_id: &str,
) -> Result<ValidatedSuppliedReference, Diagnostic> {
    validate_supplied_fasta_inner(path, path, reference_id, None)
}

/// Una base REF solicitada por coordenada durante la pasada de recuperación F5.
///
/// El lector que lo produce nunca retiene la secuencia completa: conserva la
/// cabecera vigente y una fila por sitio solicitado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct StreamedFastaSiteBase {
    sequence_ordinal: u64,
    sequence_id: SequenceId,
    position_1based: u64,
    base: u8,
}

/// Resultado de una pasada secuencial de recuperación de bases REF F5.
///
/// Además de las bases solicitadas conserva el número exacto de bases FASTA
/// visitadas. El consumidor lo compara con el recuento del preflight: así una
/// referencia que cambiase entre ambas pasadas no puede quedar oculta detrás
/// de una coincidencia accidental de los sitios elegidos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub(crate) struct StreamedFastaSiteBases {
    sites: Vec<StreamedFastaSiteBase>,
    reference_bases_processed: u64,
}

impl StreamedFastaSiteBases {
    #[must_use]
    pub(crate) fn sites(&self) -> &[StreamedFastaSiteBase] {
        &self.sites
    }

    #[must_use]
    pub(crate) const fn reference_bases_processed(&self) -> u64 {
        self.reference_bases_processed
    }
}

impl StreamedFastaSiteBase {
    #[must_use]
    pub const fn sequence_ordinal(&self) -> u64 {
        self.sequence_ordinal
    }

    #[must_use]
    pub fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn position_1based(&self) -> u64 {
        self.position_1based
    }

    #[must_use]
    pub const fn base(&self) -> u8 {
        self.base
    }
}

/// Recupera solamente las bases solicitadas con una pasada secuencial del
/// FASTA y verifica al final que los bytes siguen teniendo el SHA-256 ligado
/// por el preflight. `sites` debe estar ordenado y no contener duplicados.
///
/// La rutina está deliberadamente separada del validador que materializa un
/// `Reference`: F5 ya validó la referencia completa antes de la admisión y
/// aquí solo necesita las bases REF de los SNP ya elegidos.
pub(crate) fn stream_fasta_site_bases(
    path: &Path,
    reference_id: &str,
    expected_sha256: &str,
    sites: &[(u64, u64)],
) -> Result<StreamedFastaSiteBases, Diagnostic> {
    ReferenceId::parse(reference_id, "f5_layout")?;
    if sites.iter().any(|(_, position)| *position == 0)
        || sites.windows(2).any(|pair| pair[0] >= pair[1])
    {
        return Err(internal_error(
            "F5 solicitó sitios FASTA no ordenados, repetidos o con coordenada cero",
        ));
    }
    let file = File::open(path).map_err(|error| input_error(path, error))?;
    let mut reader = BufReader::new(DigestingReader::new(file));
    if reader
        .fill_buf()
        .map_err(|error| input_error(path, error))?
        .starts_with(&[0x1f, 0x8b])
    {
        return Err(Diagnostic::format(
            "GPS-FMT-001",
            "El FASTA suministrado está comprimido y v1 solo admite texto sin comprimir.",
        )
        .with("path", path.display().to_string())
        .with("detected_format", "gzip"));
    }

    let mut result = Vec::with_capacity(sites.len());
    let mut requested_index = 0_usize;
    let mut current_ordinal = None;
    let mut current_id = None::<SequenceId>;
    let mut current_length = 0_u64;
    let mut reference_bases_processed = 0_u64;
    let mut next_ordinal = 0_u64;
    let mut at_line_start = true;
    let mut in_header = false;
    let mut header_id = Vec::new();
    let mut header_token_finished = false;
    let mut buffer = [0_u8; 65_536];

    loop {
        let read = reader
            .read(&mut buffer)
            .map_err(|error| input_error(path, error))?;
        if read == 0 {
            break;
        }
        for raw in &buffer[..read] {
            if in_header {
                if *raw == b'\n' {
                    let id = std::str::from_utf8(&header_id).map_err(|_| empty_header(path, 0))?;
                    if id.is_empty() {
                        return Err(empty_header(path, 0));
                    }
                    current_id = Some(SequenceId::parse(id, reference_id)?);
                    current_ordinal = Some(next_ordinal);
                    next_ordinal = next_ordinal.checked_add(1).ok_or_else(|| {
                        internal_error("el ordinal de secuencia FASTA excede u64")
                    })?;
                    current_length = 0;
                    in_header = false;
                    at_line_start = true;
                } else if *raw == b'\r' || raw.is_ascii_whitespace() {
                    header_token_finished = true;
                } else if !header_token_finished {
                    header_id.push(*raw);
                }
                continue;
            }

            if at_line_start && *raw == b'>' {
                ensure_streamed_sequence_complete(
                    path,
                    sites,
                    requested_index,
                    current_ordinal,
                    current_length,
                )?;
                current_ordinal = None;
                current_id = None;
                current_length = 0;
                header_id.clear();
                header_token_finished = false;
                in_header = true;
                at_line_start = false;
                continue;
            }
            if *raw == b'\n' {
                at_line_start = true;
                continue;
            }
            if *raw == b'\r' {
                continue;
            }
            at_line_start = false;
            let ordinal = current_ordinal.ok_or_else(|| empty_header(path, 0))?;
            let sequence_id = current_id.as_ref().ok_or_else(|| empty_header(path, 0))?;
            current_length = current_length
                .checked_add(1)
                .ok_or_else(|| internal_error("la posición FASTA excede u64"))?;
            reference_bases_processed = reference_bases_processed
                .checked_add(1)
                .ok_or_else(|| internal_error("el recuento de bases FASTA visitadas excede u64"))?;
            let base = raw.to_ascii_uppercase();
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(invalid_symbol(
                    path,
                    sequence_id.as_str(),
                    current_length,
                    *raw,
                ));
            }
            if let Some((requested_ordinal, requested_position)) = sites.get(requested_index) {
                if *requested_ordinal == ordinal && *requested_position < current_length {
                    return Err(internal_error(
                        "la pasada FASTA sobrepasó un sitio F5 que no pudo recuperar",
                    ));
                }
                if *requested_ordinal == ordinal && *requested_position == current_length {
                    result.push(StreamedFastaSiteBase {
                        sequence_ordinal: ordinal,
                        sequence_id: sequence_id.clone(),
                        position_1based: current_length,
                        base,
                    });
                    requested_index = requested_index
                        .checked_add(1)
                        .ok_or_else(|| internal_error("el índice de sitios FASTA excede usize"))?;
                }
            }
        }
    }
    if in_header {
        let id = std::str::from_utf8(&header_id).map_err(|_| empty_header(path, 0))?;
        if id.is_empty() {
            return Err(empty_header(path, 0));
        }
        let _ = SequenceId::parse(id, reference_id)?;
        current_ordinal = Some(next_ordinal);
        current_length = 0;
    }
    ensure_streamed_sequence_complete(
        path,
        sites,
        requested_index,
        current_ordinal,
        current_length,
    )?;
    if requested_index != sites.len() {
        return Err(internal_error(
            "la pasada FASTA terminó antes de recuperar todos los sitios F5",
        ));
    }
    let observed = reader.into_inner().finish(path);
    if observed.digest() != expected_sha256 {
        return Err(Diagnostic::input(
            "GPS-INPUT-002",
            "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
        )
        .with("path", path.display().to_string())
        .with("expected_digest", expected_sha256)
        .with("observed_digest", observed.digest()));
    }
    Ok(StreamedFastaSiteBases {
        sites: result,
        reference_bases_processed,
    })
}

fn ensure_streamed_sequence_complete(
    path: &Path,
    sites: &[(u64, u64)],
    requested_index: usize,
    ordinal: Option<u64>,
    length: u64,
) -> Result<(), Diagnostic> {
    let Some(ordinal) = ordinal else {
        return Ok(());
    };
    if length == 0 {
        return Err(empty_fasta(path));
    }
    if sites
        .get(requested_index)
        .is_some_and(|(requested_ordinal, _)| *requested_ordinal == ordinal)
    {
        return Err(internal_error(
            "un sitio F5 elegido no existe en la secuencia FASTA verificada",
        ));
    }
    Ok(())
}

/// Valida los bytes leídos desde `read_path`, conservando `origin_path` como
/// procedencia declarada. Esta separación permite que los consumidores usen
/// la copia verificada de resultados sin perder la ruta de entrada original.
pub fn validate_supplied_fasta_with_origin(
    read_path: &Path,
    origin_path: &Path,
    reference_id: &str,
    expected: &[SequenceExpectation<'_>],
) -> Result<ValidatedSuppliedReference, Diagnostic> {
    validate_supplied_fasta_inner(read_path, origin_path, reference_id, Some(expected))
}

fn validate_supplied_fasta_inner(
    read_path: &Path,
    origin_path: &Path,
    reference_id: &str,
    expected: Option<&[SequenceExpectation<'_>]>,
) -> Result<ValidatedSuppliedReference, Diagnostic> {
    let file = File::open(read_path).map_err(|error| input_error(read_path, error))?;
    let digesting = DigestingReader::new(file);
    let mut reader = BufReader::new(digesting);
    if reader
        .fill_buf()
        .map_err(|error| input_error(read_path, error))?
        .starts_with(&[0x1f, 0x8b])
    {
        return Err(Diagnostic::format(
            "GPS-FMT-001",
            "El FASTA suministrado está comprimido y v1 solo admite texto sin comprimir.",
        )
        .with("path", read_path.display().to_string())
        .with("detected_format", "gzip"));
    }

    let expectations = expected.map(|items| {
        items
            .iter()
            .map(|item| (item.id, *item))
            .collect::<BTreeMap<_, _>>()
    });
    let mut records = Vec::with_capacity(expected.map_or(0, <[SequenceExpectation<'_>]>::len));
    let mut locations = BTreeMap::<String, u64>::new();
    let mut line = Vec::new();
    let mut line_number = 0_u64;
    let mut pending: Option<(FastaHeader, Vec<u8>, u64)> = None;

    loop {
        line.clear();
        let read = reader
            .read_until(b'\n', &mut line)
            .map_err(|error| input_error(read_path, error))?;
        if read == 0 {
            break;
        }
        line_number = line_number.checked_add(1).ok_or_else(|| {
            internal_error("FASTA line number exceeds u64")
                .with("path", read_path.display().to_string())
        })?;
        strip_line_ending(&mut line);

        if let Some(raw_header) = line.strip_prefix(b">") {
            if let Some((header, bases, header_line)) = pending.take() {
                records.push(finish_record(
                    read_path,
                    reference_id,
                    expectations.as_ref(),
                    &mut locations,
                    header,
                    bases,
                    header_line,
                )?);
            }
            let header = parse_header_with_noodles(read_path, raw_header, line_number)?;
            pending = Some((header, Vec::new(), line_number));
        } else {
            let Some((header, bases, _)) = pending.as_mut() else {
                return Err(empty_header(read_path, line_number));
            };
            for raw in &line {
                let base = raw.to_ascii_uppercase();
                if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                    let position = u64::try_from(bases.len())
                        .ok()
                        .and_then(|value| value.checked_add(1))
                        .ok_or_else(|| internal_error("FASTA sequence position exceeds u64"))?;
                    return Err(invalid_symbol(read_path, header.id_text(), position, *raw));
                }
                bases.push(base);
            }
        }
    }

    if let Some((header, bases, header_line)) = pending {
        records.push(finish_record(
            read_path,
            reference_id,
            expectations.as_ref(),
            &mut locations,
            header,
            bases,
            header_line,
        )?);
    }
    if records.is_empty() {
        return Err(empty_fasta(read_path));
    }
    if let Some(expected) = expected
        && (records.len() != expected.len()
            || records
                .iter()
                .zip(expected)
                .any(|(record, expectation)| record.id().as_str() != expectation.id))
    {
        return Err(configuration_mismatch(
            read_path,
            "reference.sequences",
            "los IDs y el orden declarados no coinciden exactamente con el FASTA",
        ));
    }

    let reference = Reference::new(
        ReferenceId::parse(reference_id, "execution")?,
        ReferenceOrigin::Supplied {
            original_path: origin_path.to_path_buf(),
        },
        records,
    )
    .map_err(|error| reference_error(read_path, error))?;
    let digesting = reader.into_inner();
    let source_identity = digesting.finish(read_path);
    Ok(ValidatedSuppliedReference {
        reference,
        source_identity,
    })
}

/// Escribe el FASTA generado canónico v1: cabecera por ID, bases mayúsculas,
/// 80 bases por línea, LF y orden del catálogo de referencia.
pub fn write_generated_fasta(
    reference: &Reference,
    destination: &mut impl io::Write,
) -> Result<(), Diagnostic> {
    let line_base_count = NonZero::new(CANONICAL_LINE_BASE_COUNT)
        .ok_or_else(|| internal_error("canonical line width must be nonzero"))?;
    let mut writer = fasta::io::writer::Builder::default()
        .set_line_base_count(line_base_count)
        .build_from_writer(destination);
    for sequence in reference.sequences() {
        let definition = fasta::record::Definition::new(sequence.id().as_str(), None);
        let bases = fasta::record::Sequence::from(
            sequence
                .sequence()
                .as_bytes()
                .map_err(|error| reference_error(Path::new("reference.fa"), error))?
                .to_vec(),
        );
        writer
            .write_record(&fasta::Record::new(definition, bases))
            .map_err(|error| {
                Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                    .with("path", "generated-reference")
                    .with("operation", "write_canonical_fasta")
                    .with("detail", error.to_string())
            })?;
    }
    Ok(())
}

/// Resumen de una secuencia validada sin conservar sus bases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SequenceSummary {
    id: SequenceId,
    class: SequenceClass,
    header: FastaHeader,
    length: u64,
    unknown_base_count: u64,
}

impl SequenceSummary {
    #[must_use]
    pub const fn id(&self) -> &SequenceId {
        &self.id
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn header(&self) -> &FastaHeader {
        &self.header
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }

    #[must_use]
    pub const fn unknown_base_count(&self) -> u64 {
        self.unknown_base_count
    }
}

/// Catálogo de un FASTA validado por completo sin materializar sus bases.
#[derive(Clone, Debug)]
pub struct ValidatedFastaSummary {
    sequences: Vec<SequenceSummary>,
    source_identity: ArtifactIdentity,
    content_sha256: String,
}

impl ValidatedFastaSummary {
    #[must_use]
    pub fn sequences(&self) -> &[SequenceSummary] {
        &self.sequences
    }

    #[must_use]
    pub const fn source_identity(&self) -> &ArtifactIdentity {
        &self.source_identity
    }

    /// Huella canónica del contenido, calculada en la misma pasada que valida
    /// las bases. Reproduce el algoritmo de `reference_content_sha256`, que no
    /// puede aplicarse sobre una referencia sin bases materializadas.
    #[must_use]
    pub fn content_sha256(&self) -> &str {
        &self.content_sha256
    }
}

/// Construye una referencia a partir de un catálogo validado sin bases.
///
/// Cada secuencia conserva su ficha —identificador, clase, cabecera, longitud y
/// recuento de bases desconocidas— y declara que sus bases no se
/// materializaron. Sirve a la ejecución que solo deriva referencias reducidas:
/// planifica con el catálogo y pide después los tramos por índice.
pub fn reference_from_summary(
    summary: &ValidatedFastaSummary,
    reference_id: &str,
    origin: crate::reference::ReferenceOrigin,
) -> Result<Reference, Diagnostic> {
    let mut sequences = Vec::with_capacity(summary.sequences().len());
    for entry in summary.sequences() {
        let bases = DnaSequence::summary(entry.length(), entry.unknown_base_count())
            .map_err(|error| reference_error(Path::new("reference.fa"), error))?;
        sequences.push(
            ReferenceSequence::from_owned(
                entry.id().clone(),
                entry.header().clone(),
                entry.class(),
                bases,
            )
            .map_err(|error| reference_error(Path::new("reference.fa"), error))?,
        );
    }
    let mut reference = Reference::new(
        ReferenceId::parse(reference_id, "execution")?,
        origin,
        sequences,
    )
    .map_err(|error| reference_error(Path::new("reference.fa"), error))?;
    reference.bind_content_sha256(summary.content_sha256().to_owned());
    Ok(reference)
}

/// Valida un FASTA suministrado sin retener sus bases.
///
/// Aplica exactamente las mismas comprobaciones que [`validate_supplied_fasta`]
/// —alfabeto A/C/G/T/N, identificadores únicos, cabeceras y longitudes
/// declaradas— pero descarta cada base en cuanto la valida, de modo que la
/// memoria no depende del tamaño del genoma. Sirve al consumidor que solo
/// necesita el catálogo: la planificación de una reducción, que trabaja con
/// longitudes, y la validación física de una publicación, que compara recuentos
/// de `N`, cabeceras y alfabeto.
///
/// `reference_fasta_summary_matches_full_read` comprueba que ambos lectores
/// coinciden secuencia a secuencia.
pub fn validate_supplied_fasta_summary(
    read_path: &Path,
    reference_id: &str,
    expected: &[SequenceExpectation<'_>],
) -> Result<ValidatedFastaSummary, Diagnostic> {
    let file = File::open(read_path).map_err(|error| input_error(read_path, error))?;
    let digesting = DigestingReader::new(file);
    let mut reader = BufReader::new(digesting);
    if reader
        .fill_buf()
        .map_err(|error| input_error(read_path, error))?
        .starts_with(&[0x1f, 0x8b])
    {
        return Err(Diagnostic::format(
            "GPS-FMT-001",
            "El FASTA suministrado está comprimido y v1 solo admite texto sin comprimir.",
        )
        .with("path", read_path.display().to_string())
        .with("detected_format", "gzip"));
    }
    let expectations = expected
        .iter()
        .map(|item| (item.id, *item))
        .collect::<BTreeMap<_, _>>();
    let mut sequences = Vec::with_capacity(expected.len());
    let mut locations = BTreeMap::<String, u64>::new();
    let mut line = Vec::new();
    let mut line_number = 0_u64;
    let mut pending: Option<(FastaHeader, u64, u64, u64)> = None;
    // Reproduce el algoritmo de `reference_content_sha256` sobre las bases que
    // van pasando, para no tener que conservarlas solo por su huella.
    let mut content = <Sha256 as Digest>::new();
    content.update(b"genome_population_simulator/reference-content/v1\0");
    content.update(
        u64::try_from(reference_id.len())
            .map_err(|_| internal_error("reference identifier length exceeds u64"))?
            .to_be_bytes(),
    );
    content.update(reference_id.as_bytes());
    let mut sequence_started = false;
    let mut normalized = Vec::new();

    loop {
        line.clear();
        let read = reader
            .read_until(b'\n', &mut line)
            .map_err(|error| input_error(read_path, error))?;
        if read == 0 {
            break;
        }
        line_number = line_number
            .checked_add(1)
            .ok_or_else(|| internal_error("FASTA line number exceeds u64"))?;
        while line
            .last()
            .is_some_and(|byte| *byte == b'\n' || *byte == b'\r')
        {
            line.pop();
        }
        if line.first() == Some(&b'>') {
            let raw_header = line[1..].to_vec();
            if let Some((header, length, unknown, header_line)) = pending.take() {
                sequences.push(finish_summary(
                    read_path,
                    reference_id,
                    &expectations,
                    &mut locations,
                    header,
                    length,
                    unknown,
                    header_line,
                )?);
            }
            let header = parse_header_with_noodles(read_path, &raw_header, line_number)?;
            let expectation = expectations.get(header.id_text()).ok_or_else(|| {
                configuration_mismatch(
                    read_path,
                    "reference.sequences.id",
                    format!("el FASTA contiene el ID no declarado {}", header.id_text()),
                )
            })?;
            content.update(
                u64::try_from(expectation.id.len())
                    .map_err(|_| internal_error("sequence identifier length exceeds u64"))?
                    .to_be_bytes(),
            );
            content.update(expectation.id.as_bytes());
            content.update(expectation.length.to_be_bytes());
            sequence_started = true;
            pending = Some((header, 0, 0, line_number));
        } else {
            let Some((header, length, unknown, _)) = pending.as_mut() else {
                return Err(empty_header(read_path, line_number));
            };
            // La huella canónica se calcula sobre las bases normalizadas, no
            // sobre los bytes tal como vienen del fichero.
            normalized.clear();
            for raw in &line {
                let base = raw.to_ascii_uppercase();
                if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                    let position = length
                        .checked_add(1)
                        .ok_or_else(|| internal_error("FASTA sequence position exceeds u64"))?;
                    return Err(invalid_symbol(read_path, header.id_text(), position, *raw));
                }
                if base == b'N' {
                    *unknown = unknown
                        .checked_add(1)
                        .ok_or_else(|| internal_error("FASTA unknown base count exceeds u64"))?;
                }
                *length = length
                    .checked_add(1)
                    .ok_or_else(|| internal_error("FASTA sequence length exceeds u64"))?;
                normalized.push(base);
            }
            if sequence_started {
                content.update(&normalized);
            }
        }
    }

    if let Some((header, length, unknown, header_line)) = pending {
        sequences.push(finish_summary(
            read_path,
            reference_id,
            &expectations,
            &mut locations,
            header,
            length,
            unknown,
            header_line,
        )?);
    }
    if sequences.is_empty() {
        return Err(empty_fasta(read_path));
    }
    if sequences.len() != expected.len()
        || sequences
            .iter()
            .zip(expected)
            .any(|(summary, expectation)| summary.id.as_str() != expectation.id)
    {
        return Err(configuration_mismatch(
            read_path,
            "reference.sequences",
            "los IDs y el orden declarados no coinciden exactamente con el FASTA",
        ));
    }
    let digesting = reader.into_inner();
    let source_identity = digesting.finish(read_path);
    Ok(ValidatedFastaSummary {
        sequences,
        source_identity,
        content_sha256: format!("{:x}", content.finalize()),
    })
}

#[allow(clippy::too_many_arguments)]
fn finish_summary(
    path: &Path,
    reference_id: &str,
    expected: &BTreeMap<&str, SequenceExpectation<'_>>,
    locations: &mut BTreeMap<String, u64>,
    header: FastaHeader,
    length: u64,
    unknown_base_count: u64,
    header_line: u64,
) -> Result<SequenceSummary, Diagnostic> {
    if length == 0 {
        return Err(empty_fasta(path).with("sequence", header.id_text()));
    }
    if let Some(first_line) = locations.insert(header.id_text().to_owned(), header_line) {
        return Err(Diagnostic::format(
            "GPS-FMT-004",
            "Un identificador de secuencia FASTA está repetido.",
        )
        .with("path", path.display().to_string())
        .with("identifier", header.id_text())
        .with(
            "locations",
            format!("line {first_line}, line {header_line}"),
        ));
    }
    let expectation = expected.get(header.id_text()).ok_or_else(|| {
        configuration_mismatch(
            path,
            "reference.sequences.id",
            format!("el FASTA contiene el ID no declarado {}", header.id_text()),
        )
    })?;
    if length != expectation.length {
        return Err(configuration_mismatch(
            path,
            "reference.sequences.length",
            format!(
                "{} declara longitud {} y el FASTA contiene {}",
                expectation.id, expectation.length, length
            ),
        ));
    }
    Ok(SequenceSummary {
        id: SequenceId::parse(expectation.id, reference_id)?,
        class: expectation.class,
        header,
        length,
        unknown_base_count,
    })
}

/// Crea un FAI solo para el consumidor indexado de aceptación. El runtime F1
/// normal no llama a esta función ni publica el índice resultante.
pub fn create_fai_for_acceptance(
    fasta_path: &Path,
    fai_path: &Path,
) -> Result<ReferenceIndex, Diagnostic> {
    let file = File::open(fasta_path).map_err(|error| input_error(fasta_path, error))?;
    let mut indexer = fasta::io::Indexer::new(BufReader::new(file));
    let mut records = Vec::new();
    loop {
        match indexer.index_record() {
            Ok(Some(record)) => records.push(record),
            Ok(None) => break,
            Err(error) => {
                return Err(Diagnostic::format(
                    "GPS-FMT-005",
                    "Una secuencia contiene un símbolo de DNA no admitido.",
                )
                .with("path", fasta_path.display().to_string())
                .with("sequence", "unknown")
                .with("position", "unknown")
                .with("found", error.to_string())
                .with("accepted", "indexable FASTA validated before FAI creation"));
            }
        }
    }
    let index = fasta::fai::Index::from(records);
    fasta::fai::fs::write(fai_path, &index).map_err(|error| {
        Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
            .with("path", fai_path.display().to_string())
            .with("operation", "write_fai_acceptance")
            .with("detail", error.to_string())
    })?;
    Ok(ReferenceIndex::new(
        ReferenceIndexKind::Fai,
        fai_path.to_path_buf(),
    ))
}

/// Consulta de aceptación que demuestra el consumidor real del FAI bajo
/// demanda sin introducir la herramienta externa en el runtime.
pub fn query_fai_for_acceptance(
    fasta_path: &Path,
    fai_path: &Path,
    region: &str,
) -> Result<Vec<u8>, Diagnostic> {
    let index = fasta::fai::fs::read(fai_path).map_err(|error| input_error(fai_path, error))?;
    let mut reader = fasta::io::indexed_reader::Builder::default()
        .set_index(index)
        .build_from_path(fasta_path)
        .map_err(|error| input_error(fasta_path, error))?;
    let region = region.parse().map_err(|error| {
        Diagnostic::format(
            "GPS-FMT-005",
            "Una secuencia contiene un símbolo de DNA no admitido.",
        )
        .with("path", fasta_path.display().to_string())
        .with("sequence", region)
        .with("position", "region")
        .with("found", format!("{error}"))
        .with("accepted", "valid indexed region")
    })?;
    reader
        .query(&region)
        .map(|record| record.sequence().as_ref().to_vec())
        .map_err(|error| input_error(fasta_path, error))
}

#[allow(clippy::too_many_arguments)]
fn finish_record(
    path: &Path,
    reference_id: &str,
    expected: Option<&BTreeMap<&str, SequenceExpectation<'_>>>,
    locations: &mut BTreeMap<String, u64>,
    header: FastaHeader,
    bases: Vec<u8>,
    header_line: u64,
) -> Result<ReferenceSequence, Diagnostic> {
    if bases.is_empty() {
        return Err(empty_fasta(path).with("sequence", header.id_text()));
    }
    if let Some(first_line) = locations.insert(header.id_text().to_owned(), header_line) {
        return Err(Diagnostic::format(
            "GPS-FMT-004",
            "Un identificador de secuencia FASTA está repetido.",
        )
        .with("path", path.display().to_string())
        .with("identifier", header.id_text())
        .with(
            "locations",
            format!("line {first_line}, line {header_line}"),
        ));
    }
    let sequence =
        DnaSequence::try_from_bytes(&bases).map_err(|error| reference_error(path, error))?;
    let (sequence_id, class) = if let Some(expected) = expected {
        let expectation = expected.get(header.id_text()).ok_or_else(|| {
            configuration_mismatch(
                path,
                "reference.sequences.id",
                format!("el FASTA contiene el ID no declarado {}", header.id_text()),
            )
        })?;
        if sequence.len() != expectation.length {
            return Err(configuration_mismatch(
                path,
                "reference.sequences.length",
                format!(
                    "{} declara longitud {} y el FASTA contiene {}",
                    expectation.id,
                    expectation.length,
                    sequence.len()
                ),
            ));
        }
        (
            SequenceId::parse(expectation.id, reference_id)?,
            expectation.class,
        )
    } else {
        (
            SequenceId::parse(header.id_text(), reference_id)?,
            SequenceClass::Chromosome,
        )
    };
    ReferenceSequence::from_owned(sequence_id, header, class, sequence)
        .map_err(|error| reference_error(path, error))
}

pub(crate) fn parse_header_with_noodles(
    path: &Path,
    raw: &[u8],
    line_number: u64,
) -> Result<FastaHeader, Diagnostic> {
    let header = FastaHeader::parse(raw).map_err(|_| empty_header(path, line_number))?;
    let first = raw
        .iter()
        .position(|byte| !byte.is_ascii_whitespace())
        .ok_or_else(|| empty_header(path, line_number))?;
    let mut canonical = Vec::with_capacity(raw.len() - first + 2);
    canonical.push(b'>');
    canonical.extend_from_slice(&raw[first..]);
    canonical.push(b'\n');
    let mut reader = fasta::io::Reader::new(Cursor::new(canonical));
    let mut definition = fasta::record::Definition::default();
    reader
        .read_definition(&mut definition)
        .map_err(|_| empty_header(path, line_number))?;
    if definition.name() != header.id_text().as_bytes() {
        return Err(internal_error(
            "noodles and the domain disagreed on the FASTA ID",
        ));
    }
    Ok(header)
}

fn strip_line_ending(line: &mut Vec<u8>) {
    if line.last() == Some(&b'\n') {
        line.pop();
    }
    if line.last() == Some(&b'\r') {
        line.pop();
    }
}

struct DigestingReader<R> {
    inner: R,
    digest: Sha256Writer,
}

impl<R> DigestingReader<R> {
    fn new(inner: R) -> Self {
        Self {
            inner,
            digest: Sha256Writer::default(),
        }
    }

    fn finish(self, path: &Path) -> ArtifactIdentity {
        self.digest.finish(path)
    }
}

impl<R: Read> Read for DigestingReader<R> {
    fn read(&mut self, buffer: &mut [u8]) -> io::Result<usize> {
        let read = self.inner.read(buffer)?;
        self.digest
            .update(&buffer[..read])
            .map_err(|diagnostic| io::Error::other(diagnostic.render_human()))?;
        Ok(read)
    }
}

pub(crate) fn input_error(path: &Path, error: io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

pub(crate) fn empty_header(path: &Path, line: u64) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-003",
        "Una cabecera FASTA no contiene un identificador de secuencia.",
    )
    .with("path", path.display().to_string())
    .with("location", format!("line {line}"))
}

pub(crate) fn empty_fasta(path: &Path) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-002",
        "El FASTA no contiene ninguna secuencia no vacía.",
    )
    .with("path", path.display().to_string())
}

pub(crate) fn invalid_symbol(path: &Path, sequence: &str, position: u64, found: u8) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-005",
        "Una secuencia contiene un símbolo de DNA no admitido.",
    )
    .with("path", path.display().to_string())
    .with("sequence", sequence)
    .with("position", position.to_string())
    .with("found", format!("0x{found:02x}"))
    .with("accepted", "A,C,G,T,N (case-insensitive)")
}

fn configuration_mismatch(path: &Path, fields: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", path.display().to_string())
    .with("fields", fields)
    .with("detail", detail.into())
}

fn reference_error(path: &Path, error: ReferenceError) -> Diagnostic {
    match error {
        ReferenceError::InvalidDnaSymbol { position, symbol } => {
            invalid_symbol(path, "unknown", position, symbol)
        }
        ReferenceError::EmptyDnaSequence | ReferenceError::EmptyReference => empty_fasta(path),
        ReferenceError::EmptyHeaderId | ReferenceError::NonUtf8HeaderId => empty_header(path, 0),
        other => internal_error(format!(
            "validated FASTA domain construction failed: {other:?}"
        )),
    }
}

fn internal_error(detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "reference_fasta")
    .with("operation", "validate_supplied_fasta")
    .with("detail", detail.into())
}
