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

/// Valida los bytes leídos desde `read_path`, conservando `origin_path` como
/// procedencia declarada. Esta separación permite que los consumidores usen
/// la copia verificada de resultados sin perder la ruta de entrada original.
pub fn validate_supplied_fasta_with_origin(
    read_path: &Path,
    origin_path: &Path,
    reference_id: &str,
    expected: &[SequenceExpectation<'_>],
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

    let expectations = expected
        .iter()
        .map(|item| (item.id, *item))
        .collect::<BTreeMap<_, _>>();
    let mut records = Vec::with_capacity(expected.len());
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
                    &expectations,
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
            &expectations,
            &mut locations,
            header,
            bases,
            header_line,
        )?);
    }
    if records.is_empty() {
        return Err(empty_fasta(read_path));
    }
    if records.len() != expected.len()
        || records
            .iter()
            .zip(expected)
            .any(|(record, expectation)| record.id().as_str() != expectation.id)
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
        let bases = fasta::record::Sequence::from(sequence.sequence().as_bytes().to_vec());
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
    expected: &BTreeMap<&str, SequenceExpectation<'_>>,
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
    let expectation = expected.get(header.id_text()).ok_or_else(|| {
        configuration_mismatch(
            path,
            "reference.sequences.id",
            format!("el FASTA contiene el ID no declarado {}", header.id_text()),
        )
    })?;
    let sequence =
        DnaSequence::try_from_bytes(&bases).map_err(|error| reference_error(path, error))?;
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
    let sequence_id = SequenceId::parse(expectation.id, reference_id)?;
    ReferenceSequence::from_owned(sequence_id, header, expectation.class, sequence)
        .map_err(|error| reference_error(path, error))
}

fn parse_header_with_noodles(
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

fn input_error(path: &Path, error: io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

fn empty_header(path: &Path, line: u64) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-003",
        "Una cabecera FASTA no contiene un identificador de secuencia.",
    )
    .with("path", path.display().to_string())
    .with("location", format!("line {line}"))
}

fn empty_fasta(path: &Path) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-002",
        "El FASTA no contiene ninguna secuencia no vacía.",
    )
    .with("path", path.display().to_string())
}

fn invalid_symbol(path: &Path, sequence: &str, position: u64, found: u8) -> Diagnostic {
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
