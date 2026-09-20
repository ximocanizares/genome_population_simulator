//! Contrato de dominio para referencias y secuencias, independiente de FASTA.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use crate::identity::{ReferenceId, SequenceId};

/// Clasificación estructural obligatoria de una secuencia de referencia.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SequenceClass {
    Chromosome,
    Scaffold,
}

/// Regla de alfabeto que superó una secuencia de DNA.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AlphabetValidation {
    AcgtnCaseInsensitive,
}

/// Errores propios del contrato de referencia, anteriores a su diagnóstico público.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceError {
    EmptyDnaSequence,
    /// Se pidieron las bases de una secuencia validada sin conservarlas.
    BasesNotMaterialized,
    InvalidDnaSymbol {
        position: u64,
        symbol: u8,
    },
    EmptyHeaderId,
    NonUtf8HeaderId,
    HeaderIdMismatch {
        expected: String,
        found: String,
    },
    EmptyReference,
    DuplicateSequenceId(String),
    LengthOverflow,
    DuplicateIndexKind(ReferenceIndexKind),
    InvalidAuthoritativeResource,
    DuplicateAuthoritativeResource,
}

/// Secuencia validada, normalizada solo en la representación interna.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DnaSequence {
    /// Bases materializadas, o `None` cuando solo se validó la ficha.
    ///
    /// Una ejecución que solo deriva referencias reducidas valida el FASTA de
    /// partida sin conservarlo: le bastan la longitud y el recuento de bases
    /// desconocidas, y pide después por índice los tramos que necesita.
    /// Separar la ficha de las bases impide que un consumidor use unas bases
    /// que nunca se cargaron.
    bases: Option<Vec<u8>>,
    length: u64,
    unknown_base_count: u64,
}

impl DnaSequence {
    pub fn try_from_bytes(input: &[u8]) -> Result<Self, ReferenceError> {
        Self::try_from_owned_bytes(input.to_vec())
    }

    /// Construye una secuencia a partir de un búfer que ya pertenece al
    /// llamador. Mantiene la misma normalización y validación que
    /// [`Self::try_from_bytes`] sin una copia adicional de las bases.
    pub(crate) fn try_from_owned_bytes(mut bases: Vec<u8>) -> Result<Self, ReferenceError> {
        if bases.is_empty() {
            return Err(ReferenceError::EmptyDnaSequence);
        }

        let mut unknown_base_count = 0_u64;
        for (offset, raw) in bases.iter_mut().enumerate() {
            let original = *raw;
            let base = original.to_ascii_uppercase();
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                let position = offset
                    .checked_add(1)
                    .and_then(|value| u64::try_from(value).ok())
                    .ok_or(ReferenceError::LengthOverflow)?;
                return Err(ReferenceError::InvalidDnaSymbol {
                    position,
                    symbol: original,
                });
            }
            unknown_base_count += u64::from(base == b'N');
            *raw = base;
        }

        let length = u64::try_from(bases.len()).map_err(|_| ReferenceError::LengthOverflow)?;
        Ok(Self {
            bases: Some(bases),
            length,
            unknown_base_count,
        })
    }

    /// Ficha de una secuencia validada cuyas bases no se conservaron.
    pub(crate) const fn summary(
        length: u64,
        unknown_base_count: u64,
    ) -> Result<Self, ReferenceError> {
        if length == 0 {
            return Err(ReferenceError::EmptyDnaSequence);
        }
        Ok(Self {
            bases: None,
            length,
            unknown_base_count,
        })
    }

    /// Bases de la secuencia, cuando se materializaron.
    pub fn as_bytes(&self) -> Result<&[u8], ReferenceError> {
        self.bases
            .as_deref()
            .ok_or(ReferenceError::BasesNotMaterialized)
    }

    /// Indica si las bases están disponibles sin intentar leerlas.
    #[must_use]
    pub const fn is_materialized(&self) -> bool {
        self.bases.is_some()
    }

    #[must_use]
    pub const fn len(&self) -> u64 {
        self.length
    }

    #[must_use]
    pub const fn is_empty(&self) -> bool {
        self.length == 0
    }

    #[must_use]
    pub const fn unknown_base_count(&self) -> u64 {
        self.unknown_base_count
    }

    #[must_use]
    pub const fn validation(&self) -> AlphabetValidation {
        AlphabetValidation::AcgtnCaseInsensitive
    }
}

/// Texto conservado sin pérdida para su serialización estructurada posterior.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReversibleText {
    Utf8(String),
    Bytes(Vec<u8>),
}

impl ReversibleText {
    fn from_bytes(bytes: &[u8]) -> Self {
        match std::str::from_utf8(bytes) {
            Ok(text) => Self::Utf8(text.to_owned()),
            Err(_) => Self::Bytes(bytes.to_vec()),
        }
    }

    #[must_use]
    pub const fn encoding(&self) -> &'static str {
        match self {
            Self::Utf8(_) => "utf-8",
            Self::Bytes(_) => "bytes-hex",
        }
    }

    #[must_use]
    pub fn encoded(&self) -> String {
        match self {
            Self::Utf8(text) => text.clone(),
            Self::Bytes(bytes) => encode_hex(bytes),
        }
    }
}

/// Cabecera FASTA separada en identidad, descripción reversible y bytes crudos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct FastaHeader {
    id_text: String,
    description: Option<ReversibleText>,
    raw: Vec<u8>,
}

impl FastaHeader {
    /// Analiza los bytes posteriores a `>` y anteriores al fin de línea.
    pub fn parse(raw: &[u8]) -> Result<Self, ReferenceError> {
        let first = raw
            .iter()
            .position(|byte| !byte.is_ascii_whitespace())
            .ok_or(ReferenceError::EmptyHeaderId)?;
        let id_end = raw[first..]
            .iter()
            .position(u8::is_ascii_whitespace)
            .map_or(raw.len(), |offset| first + offset);
        let id = &raw[first..id_end];
        if id.is_empty() {
            return Err(ReferenceError::EmptyHeaderId);
        }
        let id_text = std::str::from_utf8(id)
            .map_err(|_| ReferenceError::NonUtf8HeaderId)?
            .to_owned();
        let description_bytes = raw.get(id_end.saturating_add(1)..).unwrap_or_default();
        let description =
            (!description_bytes.is_empty()).then(|| ReversibleText::from_bytes(description_bytes));

        Ok(Self {
            id_text,
            description,
            raw: raw.to_vec(),
        })
    }

    #[must_use]
    pub fn id_text(&self) -> &str {
        &self.id_text
    }

    #[must_use]
    pub fn raw_bytes(&self) -> &[u8] {
        &self.raw
    }

    #[must_use]
    pub const fn description(&self) -> Option<&ReversibleText> {
        self.description.as_ref()
    }

    #[must_use]
    pub fn description_encoding(&self) -> Option<&'static str> {
        self.description.as_ref().map(ReversibleText::encoding)
    }

    #[must_use]
    pub fn encoded_description(&self) -> Option<String> {
        self.description.as_ref().map(ReversibleText::encoded)
    }
}

/// Metadatos validados de una secuencia y su posición en la referencia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceSequence {
    id: SequenceId,
    header: FastaHeader,
    class: SequenceClass,
    sequence: DnaSequence,
}

impl ReferenceSequence {
    pub fn new(
        id: SequenceId,
        header: FastaHeader,
        class: SequenceClass,
        sequence: &DnaSequence,
    ) -> Result<Self, ReferenceError> {
        if id.as_str() != header.id_text() {
            return Err(ReferenceError::HeaderIdMismatch {
                expected: id.as_str().to_owned(),
                found: header.id_text().to_owned(),
            });
        }
        Ok(Self {
            id,
            header,
            class,
            sequence: sequence.clone(),
        })
    }

    pub(crate) fn from_owned(
        id: SequenceId,
        header: FastaHeader,
        class: SequenceClass,
        sequence: DnaSequence,
    ) -> Result<Self, ReferenceError> {
        if id.as_str() != header.id_text() {
            return Err(ReferenceError::HeaderIdMismatch {
                expected: id.as_str().to_owned(),
                found: header.id_text().to_owned(),
            });
        }
        Ok(Self {
            id,
            header,
            class,
            sequence,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &SequenceId {
        &self.id
    }

    #[must_use]
    pub const fn header(&self) -> &FastaHeader {
        &self.header
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn sequence(&self) -> &DnaSequence {
        &self.sequence
    }

    #[must_use]
    pub fn length(&self) -> u64 {
        self.sequence.len()
    }

    #[must_use]
    pub const fn unknown_base_count(&self) -> u64 {
        self.sequence.unknown_base_count()
    }

    #[must_use]
    pub const fn alphabet_validation(&self) -> AlphabetValidation {
        self.sequence.validation()
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum ReferenceOrigin {
    Generated,
    Supplied {
        original_path: PathBuf,
    },
    Reduced {
        source_reference_id: ReferenceId,
    },
    /// Referencia derivada por **inserción de tractos** sobre una real. Lleva
    /// identidad propia por la misma razón que la reducida: el genoma final no
    /// es el fichero suministrado y ningún resultado puede quedar atribuido a
    /// un genoma que no existe (decisión S-03, punto 3).
    TractInsertion {
        source_reference_id: ReferenceId,
    },
}

#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum ReferenceIndexKind {
    Fai,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReferenceIndex {
    kind: ReferenceIndexKind,
    path: PathBuf,
}

impl ReferenceIndex {
    #[must_use]
    pub fn new(kind: ReferenceIndexKind, path: PathBuf) -> Self {
        Self { kind, path }
    }

    #[must_use]
    pub const fn kind(&self) -> ReferenceIndexKind {
        self.kind
    }

    #[must_use]
    pub fn path(&self) -> &Path {
        &self.path
    }
}

/// Catálogo ordenado de una referencia y de sus índices realmente producidos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Reference {
    id: ReferenceId,
    origin: ReferenceOrigin,
    sequences: Vec<ReferenceSequence>,
    indexes: Vec<ReferenceIndex>,
    authoritative_resource: Option<PathBuf>,
    total_length: u64,
    /// Huella del contenido calculada al validar sin conservar las bases.
    content_sha256: Option<String>,
}

impl Reference {
    pub fn new(
        id: ReferenceId,
        origin: ReferenceOrigin,
        sequences: Vec<ReferenceSequence>,
    ) -> Result<Self, ReferenceError> {
        if sequences.is_empty() {
            return Err(ReferenceError::EmptyReference);
        }
        let mut ids = BTreeSet::new();
        let mut total_length = 0_u64;
        for sequence in &sequences {
            if !ids.insert(sequence.id().as_str()) {
                return Err(ReferenceError::DuplicateSequenceId(
                    sequence.id().as_str().to_owned(),
                ));
            }
            total_length = total_length
                .checked_add(sequence.length())
                .ok_or(ReferenceError::LengthOverflow)?;
        }
        Ok(Self {
            id,
            origin,
            content_sha256: None,
            sequences,
            indexes: Vec::new(),
            authoritative_resource: None,
            total_length,
        })
    }

    #[must_use]
    pub const fn id(&self) -> &ReferenceId {
        &self.id
    }

    #[must_use]
    pub const fn origin(&self) -> &ReferenceOrigin {
        &self.origin
    }

    #[must_use]
    pub fn sequences(&self) -> &[ReferenceSequence] {
        &self.sequences
    }

    #[must_use]
    pub fn indexes(&self) -> &[ReferenceIndex] {
        &self.indexes
    }

    #[must_use]
    pub const fn total_length(&self) -> u64 {
        self.total_length
    }

    #[must_use]
    /// Huella del contenido, si se calculó al validar sin conservar las bases.
    pub fn content_sha256(&self) -> Option<&str> {
        self.content_sha256.as_deref()
    }

    /// Fija la huella calculada durante una validación sin bases.
    pub fn bind_content_sha256(&mut self, digest: String) {
        self.content_sha256 = Some(digest);
    }

    /// Entrega las secuencias por valor.
    ///
    /// Permite reconstruir una referencia con otro origen sin clonar las bases:
    /// clonarlas mantendría dos copias vivas y doblaría el pico de memoria.
    #[must_use]
    pub fn into_sequences(self) -> Vec<ReferenceSequence> {
        self.sequences
    }

    pub fn authoritative_resource(&self) -> Option<&Path> {
        self.authoritative_resource.as_deref()
    }

    pub fn bind_authoritative_resource(&mut self, path: PathBuf) -> Result<(), ReferenceError> {
        if path.is_absolute()
            || path
                .components()
                .any(|component| matches!(component, std::path::Component::ParentDir))
        {
            return Err(ReferenceError::InvalidAuthoritativeResource);
        }
        if self.authoritative_resource.replace(path).is_some() {
            return Err(ReferenceError::DuplicateAuthoritativeResource);
        }
        Ok(())
    }

    pub fn add_index(&mut self, index: ReferenceIndex) -> Result<(), ReferenceError> {
        if self
            .indexes
            .iter()
            .any(|existing| existing.kind() == index.kind())
        {
            return Err(ReferenceError::DuplicateIndexKind(index.kind()));
        }
        self.indexes.push(index);
        self.indexes.sort_by_key(ReferenceIndex::kind);
        Ok(())
    }
}

fn encode_hex(bytes: &[u8]) -> String {
    const HEX: &[u8; 16] = b"0123456789abcdef";
    let mut output = String::with_capacity(bytes.len() * 2);
    for byte in bytes {
        output.push(char::from(HEX[usize::from(byte >> 4)]));
        output.push(char::from(HEX[usize::from(byte & 0x0f)]));
    }
    output
}

#[cfg(test)]
mod tests {
    use super::{DnaSequence, ReferenceError};

    #[test]
    fn owned_dna_constructor_normalizes_in_place_and_counts_unknown_bases() {
        let sequence = DnaSequence::try_from_owned_bytes(b"aCnTn".to_vec()).expect("valid DNA");
        assert_eq!(sequence.as_bytes().expect("bases"), b"ACNTN");
        assert_eq!(sequence.unknown_base_count(), 2);
    }

    #[test]
    fn owned_dna_constructor_preserves_the_public_validation_error() {
        assert_eq!(
            DnaSequence::try_from_owned_bytes(b"AC-R".to_vec()),
            Err(ReferenceError::InvalidDnaSymbol {
                position: 3,
                symbol: b'-',
            })
        );
    }
}
