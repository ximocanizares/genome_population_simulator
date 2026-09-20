//! F1.7.2: inserción de tractos en una referencia real y publicación del
//! resultado como **referencia derivada con identidad propia**, con la tabla de
//! correspondencia entre las coordenadas del genoma de partida y las del final.
//!
//! # Las dos reglas que esto sostiene (decisión S-03)
//!
//! 1. **Las posiciones que declara el usuario se refieren al genoma de
//!    partida**, y el sistema las traduce a las del genoma final y publica esa
//!    correspondencia. Sin esa regla, declarar dos inserciones seguidas es
//!    ambiguo: la segunda no sabría si habla del original o del intermedio.
//! 2. El genoma final **no es el fichero suministrado**, así que se publica con
//!    identidad propia, igual que las reducciones, para que ningún resultado
//!    quede atribuido a un genoma que no existe.
//!
//! # Por qué esto es una extensión y no maquinaria nueva
//!
//! La reducción de referencia **ya reubica** cada segmento y
//! `schemas/reference-reduction-map-v1.schema.json` **ya publica** el intervalo
//! de origen y el derivado por segmento. Este módulo extiende esa misma forma a
//! una segunda operación, la inserción, que necesita una clase de segmento más:
//! el segmento **insertado**, que tiene intervalo derivado y no tiene intervalo
//! de origen porque sus bases no vienen de ninguna parte del genoma de partida.
//! Esa es la única diferencia, y es la razón de la versión 2 del mapa.

use serde_json::{Value, json};

use crate::coordinates::GenomicInterval;
use crate::diagnostics::Diagnostic;
use crate::identity::{ReferenceId, SequenceId};
use crate::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};

pub const COORDINATE_MAP_DOCUMENT_KIND: &str = "reference_coordinate_map";
pub const COORDINATE_MAP_SCHEMA_VERSION: &str = "2";
pub const COORDINATE_MAP_OPERATION: &str = "tract_insertion";

/// Inserción declarada por el usuario **en coordenadas del genoma de partida**.
#[derive(Clone, Copy, Debug)]
pub struct DeclaredInsertion<'a> {
    pub sequence_id: &'a str,
    /// Punto de inserción en base cero: las bases entran **antes** de esta
    /// posición del genoma de partida.
    pub source_offset: u64,
    pub bases: &'a [u8],
}

/// Clase de segmento del mapa de coordenadas.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CoordinateSegmentClass {
    /// Tramo conservado del genoma de partida.
    Retained,
    /// Bases nuevas, sin intervalo de origen.
    Inserted,
}

/// Una fila de la tabla de correspondencia.
///
/// Los dos intervalos son [`GenomicInterval`], el servicio común de coordenadas
/// del proyecto: llevan dentro a qué referencia y a qué secuencia pertenecen y
/// cuánto mide esa secuencia, y son ellos quienes convierten entre la base cero
/// semiabierta interna y la base uno cerrada que se publica. Aquí no se suma ni
/// se resta uno a mano: el desplazamiento de una base es el error clásico de
/// esta frontera y no hay ninguna razón para volver a arriesgarlo.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoordinateSegment {
    materialization_order: u64,
    class: CoordinateSegmentClass,
    sequence_class: SequenceClass,
    source_interval: Option<GenomicInterval>,
    derived_interval: GenomicInterval,
}

impl CoordinateSegment {
    #[must_use]
    pub const fn materialization_order(&self) -> u64 {
        self.materialization_order
    }

    #[must_use]
    pub const fn class(&self) -> CoordinateSegmentClass {
        self.class
    }

    #[must_use]
    pub fn sequence_id(&self) -> &str {
        self.derived_interval.sequence_id().as_str()
    }

    #[must_use]
    pub const fn sequence_class(&self) -> SequenceClass {
        self.sequence_class
    }

    /// Ausente en los segmentos insertados: sus bases no vienen del origen.
    #[must_use]
    pub const fn source_interval(&self) -> Option<&GenomicInterval> {
        self.source_interval.as_ref()
    }

    #[must_use]
    pub const fn derived_interval(&self) -> &GenomicInterval {
        &self.derived_interval
    }
}

/// Tabla de correspondencia completa entre dos referencias.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CoordinateMap {
    source_reference_id: ReferenceId,
    derived_reference_id: ReferenceId,
    segments: Vec<CoordinateSegment>,
}

impl CoordinateMap {
    #[must_use]
    pub const fn source_reference_id(&self) -> &ReferenceId {
        &self.source_reference_id
    }

    #[must_use]
    pub const fn derived_reference_id(&self) -> &ReferenceId {
        &self.derived_reference_id
    }

    #[must_use]
    pub fn segments(&self) -> &[CoordinateSegment] {
        &self.segments
    }

    /// Traduce una posición **del genoma de partida** a la del genoma final.
    ///
    /// Devuelve `None` si la secuencia no existe o la posición cae fuera de
    /// ella. Es la operación que hace útil el mapa: sin ella, publicar los dos
    /// intervalos no serviría de nada a quien lea el documento.
    #[must_use]
    pub fn translate_source_position(&self, sequence_id: &str, source_offset: u64) -> Option<u64> {
        self.segments
            .iter()
            .filter(|segment| segment.sequence_id() == sequence_id)
            .find_map(|segment| {
                let source = segment.source_interval.as_ref()?;
                (source.start()..source.end())
                    .contains(&source_offset)
                    .then(|| segment.derived_interval.start() + (source_offset - source.start()))
            })
    }
}

/// Referencia derivada por inserción de tractos, con su mapa.
#[derive(Clone, Debug)]
pub struct InsertedReference {
    reference: Reference,
    map: CoordinateMap,
}

impl InsertedReference {
    #[must_use]
    pub const fn reference(&self) -> &Reference {
        &self.reference
    }

    #[must_use]
    pub const fn map(&self) -> &CoordinateMap {
        &self.map
    }

    #[must_use]
    pub fn into_parts(self) -> (Reference, CoordinateMap) {
        (self.reference, self.map)
    }
}

/// Inserta los tractos declarados y devuelve la referencia derivada y su mapa.
///
/// Todas las posiciones declaradas se leen en coordenadas del genoma de
/// partida, **también la segunda y las siguientes**: ninguna se reinterpreta
/// sobre el resultado parcial de las anteriores.
pub fn insert_tracts(
    source: &Reference,
    derived_reference_id: &str,
    insertions: &[DeclaredInsertion<'_>],
) -> Result<InsertedReference, Diagnostic> {
    let derived_id = ReferenceId::parse(derived_reference_id, "execution")?;
    if derived_id == *source.id() {
        return Err(configuration_error(
            "reference.derived_id",
            format!(
                "the derived reference must not reuse the supplied identifier {}",
                source.id()
            ),
        ));
    }
    for insertion in insertions {
        if insertion.bases.is_empty() {
            return Err(configuration_error(
                "reference.insertions.bases",
                format!("insertion on {} declares no bases", insertion.sequence_id),
            ));
        }
        for &base in insertion.bases {
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(configuration_error(
                    "reference.insertions.bases",
                    format!(
                        "insertion on {} contains 0x{base:02x}, which is not A, C, G, T or N",
                        insertion.sequence_id
                    ),
                ));
            }
        }
        if !source
            .sequences()
            .iter()
            .any(|sequence| sequence.id().as_str() == insertion.sequence_id)
        {
            return Err(configuration_error(
                "reference.insertions.sequence",
                format!(
                    "unknown sequence {} for a declared insertion",
                    insertion.sequence_id
                ),
            ));
        }
    }

    let mut segments: Vec<CoordinateSegment> = Vec::new();
    let mut derived_sequences = Vec::with_capacity(source.sequences().len());
    let mut materialization_order = 0_u64;

    for sequence in source.sequences() {
        let bases = sequence.sequence().as_bytes().map_err(|error| {
            internal_error(format!(
                "the supplied reference sequence is not materialized: {error:?}"
            ))
        })?;
        let length = sequence.length();
        // Orden estable por punto de inserción: dos inserciones declaradas en
        // el mismo punto conservan el orden en que se declararon.
        let mut here: Vec<&DeclaredInsertion<'_>> = insertions
            .iter()
            .filter(|insertion| insertion.sequence_id == sequence.id().as_str())
            .collect();
        here.sort_by_key(|insertion| insertion.source_offset);
        for insertion in &here {
            if insertion.source_offset > length {
                return Err(configuration_error(
                    "reference.insertions.start",
                    format!(
                        "insertion on {} at {} runs past the {length}-base sequence",
                        insertion.sequence_id, insertion.source_offset
                    ),
                ));
            }
        }
        // La longitud de la secuencia derivada se conoce antes de construirla,
        // y hace falta antes: `GenomicInterval` la lleva dentro y comprueba con
        // ella que ningún intervalo se sale.
        let added: u64 = here
            .iter()
            .map(|insertion| insertion.bases.len() as u64)
            .sum();
        let derived_length = length.checked_add(added).ok_or_else(|| {
            configuration_error(
                "reference.insertions",
                format!("the insertions on {} overflow the sequence", sequence.id()),
            )
        })?;
        let source_sequence_id = SequenceId::parse(sequence.id().as_str(), source.id().as_str())?;
        let derived_sequence_id = SequenceId::parse(sequence.id().as_str(), derived_id.as_str())?;

        let mut derived_bases: Vec<u8> =
            Vec::with_capacity(bases.len() + usize::try_from(added).unwrap_or(0));
        let mut source_cursor = 0_u64;
        let mut derived_cursor = 0_u64;

        for insertion in &here {
            let point = insertion.source_offset;
            if point > source_cursor {
                derived_bases.extend_from_slice(&bases[source_cursor as usize..point as usize]);
                materialization_order += 1;
                segments.push(retained_segment(
                    materialization_order,
                    sequence.class(),
                    source.id(),
                    &source_sequence_id,
                    length,
                    source_cursor,
                    point,
                    &derived_id,
                    &derived_sequence_id,
                    derived_length,
                    derived_cursor,
                )?);
                derived_cursor += point - source_cursor;
                source_cursor = point;
            }

            derived_bases.extend_from_slice(insertion.bases);
            materialization_order += 1;
            let inserted_length = insertion.bases.len() as u64;
            segments.push(CoordinateSegment {
                materialization_order,
                class: CoordinateSegmentClass::Inserted,
                sequence_class: sequence.class(),
                source_interval: None,
                derived_interval: interval(
                    &derived_id,
                    &derived_sequence_id,
                    derived_cursor,
                    derived_cursor + inserted_length,
                    derived_length,
                )?,
            });
            derived_cursor += inserted_length;
        }
        if length > source_cursor {
            derived_bases.extend_from_slice(&bases[source_cursor as usize..length as usize]);
            materialization_order += 1;
            segments.push(retained_segment(
                materialization_order,
                sequence.class(),
                source.id(),
                &source_sequence_id,
                length,
                source_cursor,
                length,
                &derived_id,
                &derived_sequence_id,
                derived_length,
                derived_cursor,
            )?);
        }

        let header =
            FastaHeader::parse(derived_sequence_id.as_str().as_bytes()).map_err(|error| {
                internal_error(format!("derived sequence header was rejected: {error:?}"))
            })?;
        let dna = DnaSequence::try_from_bytes(&derived_bases)
            .map_err(|error| internal_error(format!("derived sequence was rejected: {error:?}")))?;
        derived_sequences.push(
            ReferenceSequence::new(derived_sequence_id, header, sequence.class(), &dna).map_err(
                |error| internal_error(format!("derived sequence construction failed: {error:?}")),
            )?,
        );
    }

    let reference = Reference::new(
        derived_id.clone(),
        ReferenceOrigin::TractInsertion {
            source_reference_id: source.id().clone(),
        },
        derived_sequences,
    )
    .map_err(|error| internal_error(format!("derived reference construction failed: {error:?}")))?;

    Ok(InsertedReference {
        reference,
        map: CoordinateMap {
            source_reference_id: source.id().clone(),
            derived_reference_id: derived_id,
            segments,
        },
    })
}

#[allow(clippy::too_many_arguments)]
fn retained_segment(
    materialization_order: u64,
    sequence_class: SequenceClass,
    source_reference_id: &ReferenceId,
    source_sequence_id: &SequenceId,
    source_length: u64,
    source_start: u64,
    source_end: u64,
    derived_reference_id: &ReferenceId,
    derived_sequence_id: &SequenceId,
    derived_length: u64,
    derived_start: u64,
) -> Result<CoordinateSegment, Diagnostic> {
    Ok(CoordinateSegment {
        materialization_order,
        class: CoordinateSegmentClass::Retained,
        sequence_class,
        source_interval: Some(interval(
            source_reference_id,
            source_sequence_id,
            source_start,
            source_end,
            source_length,
        )?),
        derived_interval: interval(
            derived_reference_id,
            derived_sequence_id,
            derived_start,
            derived_start + (source_end - source_start),
            derived_length,
        )?,
    })
}

fn interval(
    reference_id: &ReferenceId,
    sequence_id: &SequenceId,
    start: u64,
    end: u64,
    sequence_length: u64,
) -> Result<GenomicInterval, Diagnostic> {
    GenomicInterval::new(
        reference_id.clone(),
        sequence_id.clone(),
        start,
        end,
        sequence_length,
    )
    .map_err(|error| {
        internal_error(format!(
            "coordinate map segment [{start}, {end}) is not a valid interval: {error:?}"
        ))
    })
}

/// Documento publicable del mapa, en la forma del mapa de reducción y con la
/// clase de segmento nueva. Los intervalos se publican **cerrados y en base
/// uno**, igual que los del mapa de reducción.
pub fn coordinate_map_document(map: &CoordinateMap) -> Result<Value, Diagnostic> {
    Ok(json!({
        "document_kind": COORDINATE_MAP_DOCUMENT_KIND,
        "schema_version": COORDINATE_MAP_SCHEMA_VERSION,
        "operation": COORDINATE_MAP_OPERATION,
        "source_reference_id": map.source_reference_id.as_str(),
        "derived_reference_id": map.derived_reference_id.as_str(),
        "segments": map.segments.iter().map(|segment| {
            let derived = published_interval(&segment.derived_interval)?;
            let mut entry = json!({
                "materialization_order": segment.materialization_order,
                "segment_class": match segment.class {
                    CoordinateSegmentClass::Retained => "retained",
                    CoordinateSegmentClass::Inserted => "inserted",
                },
                "sequence_id": segment.sequence_id(),
                "sequence_class": match segment.sequence_class {
                    SequenceClass::Chromosome => "chromosome",
                    SequenceClass::Scaffold => "scaffold",
                },
                "derived_interval": derived,
                "orientation": "forward",
            });
            if let Some(source) = &segment.source_interval {
                entry["source_interval"] = published_interval(source)?;
            }
            Ok(entry)
        }).collect::<Result<Vec<_>, Diagnostic>>()?,
    }))
}

/// Convierte a la forma publicada —cerrada y en base uno— con el servicio
/// común, no con una suma a mano.
fn published_interval(interval: &GenomicInterval) -> Result<Value, Diagnostic> {
    let closed = interval.to_closed_1().map_err(|error| {
        internal_error(format!(
            "interval cannot be published in base one: {error:?}"
        ))
    })?;
    Ok(json!({"start": closed.start(), "end": closed.end()}))
}

fn configuration_error(field: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", "effective-config")
    .with("fields", field)
    .with("detail", detail.into())
}

fn internal_error(detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal_invariant("reference_tract_insertion", "insert_tracts")
        .with("detail", detail.into())
}
