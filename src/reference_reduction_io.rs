//! Materialización canónica y mapa de una referencia reducida ya planificada.

use std::collections::{BTreeMap, BTreeSet};

use serde_json::{Value, json};

use crate::diagnostics::Diagnostic;
use crate::identity::SequenceId;
use crate::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use crate::reference_reduction::{MaterializedPlan, ReductionSegmentClass};

pub const REDUCTION_MAP_DOCUMENT_KIND: &str = "reference_reduction_map";
pub const REDUCTION_MAP_SCHEMA_VERSION: &str = "1";

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ReductionMapParameters {
    pub region_size_bases: u64,
    pub target_bases: u64,
}

pub fn build_reduced_reference(
    source: &Reference,
    plan: &MaterializedPlan,
) -> Result<Reference, Diagnostic> {
    if source.id() != plan.source_reference_id() {
        return Err(map_error(
            plan.derived_reference_id().as_str(),
            "source_reference",
            "el plan no pertenece al FASTA autoritativo",
        ));
    }
    let mut derived_sequences = Vec::new();
    let mut segments = plan.segments().iter().peekable();
    for source_sequence in source.sequences() {
        let Some(first_segment) = segments.peek() else {
            break;
        };
        if first_segment.source_interval().sequence_id() != source_sequence.id() {
            continue;
        }
        let mut bases = Vec::new();
        while let Some(segment) = segments
            .next_if(|segment| segment.source_interval().sequence_id() == source_sequence.id())
        {
            if source_sequence.class() != segment.sequence_class() {
                return Err(map_error(
                    plan.derived_reference_id().as_str(),
                    "sequence_class",
                    "la clase del segmento no coincide con la fuente",
                ));
            }
            let start = usize::try_from(segment.source_interval().start()).map_err(|_| {
                map_error(
                    plan.derived_reference_id().as_str(),
                    "source_interval",
                    "la coordenada inicial no cabe en memoria direccionable",
                )
            })?;
            let end = usize::try_from(segment.source_interval().end()).map_err(|_| {
                map_error(
                    plan.derived_reference_id().as_str(),
                    "source_interval",
                    "la coordenada final no cabe en memoria direccionable",
                )
            })?;
            let fragment = source_sequence
                .sequence()
                .as_bytes()
                .get(start..end)
                .ok_or_else(|| {
                    map_error(
                        plan.derived_reference_id().as_str(),
                        "source_interval",
                        "el intervalo fuente queda fuera de la secuencia",
                    )
                })?;
            bases.extend_from_slice(fragment);
        }
        let id = SequenceId::parse(
            source_sequence.id().as_str(),
            plan.derived_reference_id().as_str(),
        )?;
        let header = FastaHeader::parse(id.as_str().as_bytes())
            .map_err(|error| internal_error("build_reduced_header", error))?;
        let dna = DnaSequence::try_from_owned_bytes(bases)
            .map_err(|error| internal_error("build_reduced_sequence", error))?;
        derived_sequences.push(
            ReferenceSequence::from_owned(id, header, source_sequence.class(), dna)
                .map_err(|error| internal_error("build_reduced_reference_sequence", error))?,
        );
    }
    if segments.next().is_some() {
        return Err(map_error(
            plan.derived_reference_id().as_str(),
            "source_sequence",
            "el mapa referencia una secuencia inexistente o fuera de orden",
        ));
    }
    let source_reference_id = source.id().clone();
    Reference::new(
        plan.derived_reference_id().clone(),
        ReferenceOrigin::Reduced {
            source_reference_id,
        },
        derived_sequences,
    )
    .map_err(|error| internal_error("build_reduced_reference", error))
}

#[must_use]
pub fn reduction_map_document(plan: &MaterializedPlan) -> Value {
    json!({
        "document_kind": REDUCTION_MAP_DOCUMENT_KIND,
        "schema_version": REDUCTION_MAP_SCHEMA_VERSION,
        "source_reference_id": plan.source_reference_id().as_str(),
        "derived_reference_id": plan.derived_reference_id().as_str(),
        "segments": plan.segments().iter().map(|segment| {
            json!({
                "materialization_order": segment.materialization_order(),
                "selection_order": segment.selection_order(),
                "segment_class": match segment.class() {
                    ReductionSegmentClass::Full => "full",
                    ReductionSegmentClass::Residual => "residual",
                },
                "sequence_id": segment.source_interval().sequence_id().as_str(),
                "sequence_class": match segment.sequence_class() {
                    SequenceClass::Chromosome => "chromosome",
                    SequenceClass::Scaffold => "scaffold",
                },
                "source_interval": {"start": segment.source_interval().start() + 1, "end": segment.source_interval().end()},
                "derived_interval": {"start": segment.derived_interval().start() + 1, "end": segment.derived_interval().end()},
                "orientation": "forward",
            })
        }).collect::<Vec<_>>(),
    })
}

pub fn validate_reduced_reference(
    source: &Reference,
    derived: &Reference,
    plan: &MaterializedPlan,
) -> Result<(), Diagnostic> {
    let rebuilt = build_reduced_reference(source, plan)?;
    if rebuilt.id() != derived.id() || rebuilt.sequences().len() != derived.sequences().len() {
        return Err(map_error(
            derived.id().as_str(),
            "derived_catalog",
            "el catálogo derivado no coincide con la reconstrucción",
        ));
    }
    for (expected, observed) in rebuilt.sequences().iter().zip(derived.sequences()) {
        if expected.id() != observed.id()
            || expected.class() != observed.class()
            || expected.sequence().as_bytes() != observed.sequence().as_bytes()
            || expected.unknown_base_count() != observed.unknown_base_count()
        {
            return Err(map_error(
                derived.id().as_str(),
                "derived_sequence",
                "una secuencia derivada no coincide base a base con origen+mapa",
            ));
        }
    }
    Ok(())
}

/// Valida el documento serializado y reconstruye sus bases sin reutilizar el plan
/// que lo produjo. De este modo una divergencia del escritor no puede declararse
/// validada únicamente porque el plan en memoria era correcto.
pub fn validate_reduction_map_document(
    source: &Reference,
    derived: &Reference,
    document: &Value,
    parameters: ReductionMapParameters,
) -> Result<(), Diagnostic> {
    let derived_id = derived.id().as_str();
    let object = document
        .as_object()
        .ok_or_else(|| map_error(derived_id, "document", "la raíz no es un objeto JSON"))?;
    require_exact_fields(
        object,
        &[
            "document_kind",
            "schema_version",
            "source_reference_id",
            "derived_reference_id",
            "segments",
        ],
        derived_id,
        "document_fields",
    )?;
    require_text(object, "document_kind", derived_id, "document")
        .filter(|value| *value == REDUCTION_MAP_DOCUMENT_KIND)
        .ok_or_else(|| map_error(derived_id, "document_kind", "clase de documento inesperada"))?;
    require_text(object, "schema_version", derived_id, "document")
        .filter(|value| *value == REDUCTION_MAP_SCHEMA_VERSION)
        .ok_or_else(|| map_error(derived_id, "schema_version", "versión de mapa inesperada"))?;
    if require_text(object, "source_reference_id", derived_id, "document")
        != Some(source.id().as_str())
        || require_text(object, "derived_reference_id", derived_id, "document") != Some(derived_id)
    {
        return Err(map_error(
            derived_id,
            "reference_ids",
            "los identificadores del mapa no coinciden con las referencias",
        ));
    }
    let segments = object
        .get("segments")
        .and_then(Value::as_array)
        .ok_or_else(|| map_error(derived_id, "segments", "falta el catálogo de segmentos"))?;
    if segments.is_empty() {
        return Err(map_error(
            derived_id,
            "segments",
            "el catálogo de segmentos está vacío",
        ));
    }
    if parameters.region_size_bases == 0 || parameters.target_bases == 0 {
        return Err(map_error(
            derived_id,
            "effective_parameters",
            "región y objetivo efectivos deben ser mayores que cero",
        ));
    }
    let expected_full_count = parameters.target_bases / parameters.region_size_bases;
    let expected_residual = parameters.target_bases % parameters.region_size_bases;
    let expected_segment_count = expected_full_count
        .checked_add(u64::from(expected_residual > 0))
        .ok_or_else(|| map_error(derived_id, "effective_parameters", "recuento desbordado"))?;
    if u64::try_from(segments.len()).ok() != Some(expected_segment_count) {
        return Err(map_error(
            derived_id,
            "segment_count",
            "el número de segmentos no coincide con región y objetivo efectivos",
        ));
    }

    let mut rebuilt = BTreeMap::<&str, Vec<u8>>::new();
    let mut rebuilt_order = Vec::<&str>::new();
    let mut previous_source_end = BTreeMap::<&str, u64>::new();
    let mut previous_derived_end = BTreeMap::<&str, u64>::new();
    let mut selection_orders = BTreeSet::new();
    let mut residual_count = 0_usize;
    let mut previous_sequence_index = 0_usize;
    for (index, segment) in segments.iter().enumerate() {
        let segment = segment
            .as_object()
            .ok_or_else(|| map_error(derived_id, "segment", "un segmento no es un objeto JSON"))?;
        require_exact_fields(
            segment,
            &[
                "materialization_order",
                "selection_order",
                "segment_class",
                "sequence_id",
                "sequence_class",
                "source_interval",
                "derived_interval",
                "orientation",
            ],
            derived_id,
            "segment_fields",
        )?;
        let expected_order = u64::try_from(index)
            .ok()
            .and_then(|value| value.checked_add(1))
            .ok_or_else(|| map_error(derived_id, "materialization_order", "orden desbordado"))?;
        if require_u64(segment, "materialization_order", derived_id)? != expected_order {
            return Err(map_error(
                derived_id,
                "materialization_order",
                "el orden materializado no es contiguo",
            ));
        }
        let selection_order = require_u64(segment, "selection_order", derived_id)?;
        if selection_order == 0 || !selection_orders.insert(selection_order) {
            return Err(map_error(
                derived_id,
                "selection_order",
                "el orden de selección es cero o está repetido",
            ));
        }
        let declared_class = require_text(segment, "segment_class", derived_id, "segment");
        match declared_class {
            Some("full") => {}
            Some("residual") => residual_count += 1,
            _ => {
                return Err(map_error(
                    derived_id,
                    "segment_class",
                    "clase de segmento desconocida",
                ));
            }
        }
        if residual_count > 1
            || require_text(segment, "orientation", derived_id, "segment") != Some("forward")
        {
            return Err(map_error(
                derived_id,
                "segment_contract",
                "orientación o cardinalidad residual inválida",
            ));
        }
        let sequence_id = require_text(segment, "sequence_id", derived_id, "segment")
            .ok_or_else(|| map_error(derived_id, "sequence_id", "falta la secuencia"))?;
        let (sequence_index, source_sequence) = source
            .sequences()
            .iter()
            .enumerate()
            .find(|(_, sequence)| sequence.id().as_str() == sequence_id)
            .ok_or_else(|| {
                map_error(
                    derived_id,
                    "source_sequence",
                    "el mapa referencia una secuencia inexistente",
                )
            })?;
        if index > 0 && sequence_index < previous_sequence_index {
            return Err(map_error(
                derived_id,
                "source_order",
                "los segmentos no respetan el catálogo de secuencias fuente",
            ));
        }
        previous_sequence_index = sequence_index;
        let expected_class = match source_sequence.class() {
            SequenceClass::Chromosome => "chromosome",
            SequenceClass::Scaffold => "scaffold",
        };
        if require_text(segment, "sequence_class", derived_id, "segment") != Some(expected_class) {
            return Err(map_error(
                derived_id,
                "sequence_class",
                "la clase del mapa no coincide con el origen",
            ));
        }
        let (source_start, source_end) = require_interval(segment, "source_interval", derived_id)?;
        let (derived_start, derived_end) =
            require_interval(segment, "derived_interval", derived_id)?;
        let segment_length = source_end - source_start + 1;
        let (expected_class, expected_length) = if selection_order <= expected_full_count {
            ("full", parameters.region_size_bases)
        } else if expected_residual > 0 && selection_order == expected_segment_count {
            ("residual", expected_residual)
        } else {
            return Err(map_error(
                derived_id,
                "selection_order",
                "el orden de selección no corresponde a región y residual efectivos",
            ));
        };
        if declared_class != Some(expected_class) || segment_length != expected_length {
            return Err(map_error(
                derived_id,
                "segment_class",
                "la etiqueta o longitud full/residual contradice los parámetros efectivos",
            ));
        }
        if source_end > source_sequence.length()
            || previous_source_end
                .get(sequence_id)
                .is_some_and(|previous| source_start <= *previous)
        {
            return Err(map_error(
                derived_id,
                "source_interval",
                "intervalo fuente fuera de rango, solapado o desordenado",
            ));
        }
        let expected_derived_start = previous_derived_end
            .get(sequence_id)
            .copied()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or_else(|| map_error(derived_id, "derived_interval", "coordenada desbordada"))?;
        if derived_start != expected_derived_start
            || derived_end - derived_start != source_end - source_start
        {
            return Err(map_error(
                derived_id,
                "derived_interval",
                "la proyección derivada no es contigua o cambia la longitud",
            ));
        }
        let start = usize::try_from(source_start - 1)
            .map_err(|_| map_error(derived_id, "source_interval", "inicio no direccionable"))?;
        let end = usize::try_from(source_end)
            .map_err(|_| map_error(derived_id, "source_interval", "fin no direccionable"))?;
        if !rebuilt.contains_key(sequence_id) {
            rebuilt_order.push(sequence_id);
        }
        rebuilt
            .entry(sequence_id)
            .or_default()
            .extend_from_slice(&source_sequence.sequence().as_bytes()[start..end]);
        previous_source_end.insert(sequence_id, source_end);
        previous_derived_end.insert(sequence_id, derived_end);
    }
    if selection_orders.len() != segments.len()
        || selection_orders
            .iter()
            .copied()
            .ne(1..=u64::try_from(segments.len()).map_err(|_| {
                map_error(derived_id, "selection_order", "catálogo demasiado grande")
            })?)
    {
        return Err(map_error(
            derived_id,
            "selection_order",
            "el orden de selección no forma una secuencia completa",
        ));
    }
    if rebuilt.len() != derived.sequences().len() {
        return Err(map_error(
            derived_id,
            "derived_catalog",
            "el mapa omite o inventa secuencias derivadas",
        ));
    }
    for (index, sequence) in derived.sequences().iter().enumerate() {
        if rebuilt_order.get(index).copied() != Some(sequence.id().as_str()) {
            return Err(map_error(
                derived_id,
                "derived_catalog_order",
                "el orden del catálogo derivado no coincide con el mapa materializado",
            ));
        }
        let source_class = source
            .sequences()
            .iter()
            .find(|source_sequence| source_sequence.id().as_str() == sequence.id().as_str())
            .map(ReferenceSequence::class);
        if source_class != Some(sequence.class()) {
            return Err(map_error(
                derived_id,
                "derived_sequence_class",
                "la clase del catálogo derivado no coincide con la fuente y el mapa",
            ));
        }
        let Some(expected) = rebuilt.get(sequence.id().as_str()) else {
            return Err(map_error(
                derived_id,
                "derived_sequence",
                "falta una secuencia derivada en el mapa",
            ));
        };
        if expected.as_slice() != sequence.sequence().as_bytes()
            || previous_derived_end.get(sequence.id().as_str()) != Some(&sequence.length())
        {
            return Err(map_error(
                derived_id,
                "derived_sequence",
                "las bases o la cobertura derivadas no coinciden con el mapa",
            ));
        }
    }
    let rebuilt_bases = rebuilt.values().try_fold(0_u64, |total, bases| {
        total.checked_add(u64::try_from(bases.len()).ok()?)
    });
    if rebuilt_bases != Some(parameters.target_bases) {
        return Err(map_error(
            derived_id,
            "target_bases",
            "la reconstrucción no tiene el objetivo efectivo exacto",
        ));
    }
    Ok(())
}

fn require_exact_fields(
    object: &serde_json::Map<String, Value>,
    expected: &[&str],
    reference: &str,
    invariant: &str,
) -> Result<(), Diagnostic> {
    if object.len() != expected.len()
        || object
            .keys()
            .any(|field| !expected.contains(&field.as_str()))
    {
        return Err(map_error(
            reference,
            invariant,
            "el objeto contiene campos ausentes o desconocidos",
        ));
    }
    Ok(())
}

fn require_text<'a>(
    object: &'a serde_json::Map<String, Value>,
    field: &str,
    _reference: &str,
    _invariant: &str,
) -> Option<&'a str> {
    object.get(field).and_then(Value::as_str)
}

fn require_u64(
    object: &serde_json::Map<String, Value>,
    field: &str,
    reference: &str,
) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| map_error(reference, field, "falta un entero sin signo"))
}

fn require_interval(
    segment: &serde_json::Map<String, Value>,
    field: &str,
    reference: &str,
) -> Result<(u64, u64), Diagnostic> {
    let interval = segment
        .get(field)
        .and_then(Value::as_object)
        .ok_or_else(|| map_error(reference, field, "falta el intervalo"))?;
    require_exact_fields(interval, &["start", "end"], reference, field)?;
    let start = require_u64(interval, "start", reference)?;
    let end = require_u64(interval, "end", reference)?;
    if start == 0 || start > end {
        return Err(map_error(
            reference,
            field,
            "el intervalo base 1 cerrado está vacío o empieza en cero",
        ));
    }
    Ok((start, end))
}

fn map_error(reference: &str, invariant: &str, detail: &str) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-007",
        "El mapa de reducción no permite reconstruir exactamente la referencia derivada.",
    )
    .with("path", "reference-reduction-map")
    .with("reference", reference)
    .with("invariant", invariant)
    .with("detail", detail)
}

fn internal_error(operation: &str, error: impl std::fmt::Debug) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "reference_reduction_io")
    .with("operation", operation)
    .with("detail", format!("{error:?}"))
}
