//! Materialización canónica y mapa de una referencia reducida ya planificada.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::BufReader;
use std::path::{Path, PathBuf};

use noodles_fasta as fasta;

use serde_json::{Value, json};

use crate::diagnostics::Diagnostic;
use crate::identity::{ReferenceId, SequenceId};
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

/// Origen de las bases de una referencia para reconstruir una reducción.
///
/// La reconstrucción solo necesita recorrer el catálogo en orden y pedir los
/// tramos que indica el mapa, nunca la secuencia entera. Separarlo del tipo
/// `Reference` permite servir esos tramos desde memoria o desde un índice, sin
/// que la lógica de reconstrucción sepa cuál de los dos la atiende.
pub trait ReductionSource {
    fn reference_id(&self) -> &ReferenceId;

    /// Catálogo en el orden físico del FASTA.
    fn catalog(&self) -> &[SourceSequenceEntry];

    /// Bases del intervalo semiabierto `[start, end)` de una secuencia.
    fn interval(
        &self,
        sequence_id: &SequenceId,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, Diagnostic>;
}

/// Entrada de catálogo: lo único que la reconstrucción necesita saber de una
/// secuencia antes de pedir sus tramos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SourceSequenceEntry {
    id: SequenceId,
    class: SequenceClass,
    length: u64,
}

impl SourceSequenceEntry {
    #[must_use]
    pub const fn new(id: SequenceId, class: SequenceClass, length: u64) -> Self {
        Self { id, class, length }
    }

    #[must_use]
    pub const fn id(&self) -> &SequenceId {
        &self.id
    }

    #[must_use]
    pub const fn class(&self) -> SequenceClass {
        self.class
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

/// Sirve los tramos desde una referencia ya materializada en memoria.
pub struct MaterializedSource<'a> {
    reference: &'a Reference,
    catalog: Vec<SourceSequenceEntry>,
}

impl<'a> MaterializedSource<'a> {
    #[must_use]
    pub fn new(reference: &'a Reference) -> Self {
        let catalog = reference
            .sequences()
            .iter()
            .map(|sequence| {
                SourceSequenceEntry::new(sequence.id().clone(), sequence.class(), sequence.length())
            })
            .collect();
        Self { reference, catalog }
    }
}

impl ReductionSource for MaterializedSource<'_> {
    fn reference_id(&self) -> &ReferenceId {
        self.reference.id()
    }

    fn catalog(&self) -> &[SourceSequenceEntry] {
        &self.catalog
    }

    fn interval(
        &self,
        sequence_id: &SequenceId,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, Diagnostic> {
        let sequence = self
            .reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == sequence_id)
            .ok_or_else(|| {
                map_error(
                    sequence_id.as_str(),
                    "source_sequence",
                    "el mapa referencia una secuencia inexistente",
                )
            })?;
        sequence
            .sequence()
            .as_bytes()
            .map_err(|_| {
                map_error(
                    sequence_id.as_str(),
                    "source_sequence",
                    "la referencia se validó sin conservar sus bases",
                )
            })?
            .get(start..end)
            .map(<[u8]>::to_vec)
            .ok_or_else(|| {
                map_error(
                    sequence_id.as_str(),
                    "source_interval",
                    "el intervalo fuente queda fuera de la secuencia",
                )
            })
    }
}

/// Sirve los tramos desde el FASTA en disco mediante un índice FAI.
///
/// El índice se construye en un temporal del llamador y no se publica: DTC-041
/// fija que los índices solo se crean cuando se necesitan. La memoria no depende
/// del tamaño del genoma, solo de los tramos que se piden.
pub struct IndexedSource {
    reference_id: ReferenceId,
    catalog: Vec<SourceSequenceEntry>,
    fasta_path: PathBuf,
    index: fasta::fai::Index,
}

impl IndexedSource {
    /// Construye el índice del FASTA ya validado y lo deja listo para consultas.
    ///
    /// El catálogo procede de la validación previa, no del índice: el FAI
    /// describe posiciones, no la clase de cada secuencia ni su conformidad con
    /// lo declarado.
    pub fn new(
        reference_id: ReferenceId,
        catalog: Vec<SourceSequenceEntry>,
        fasta_path: &Path,
    ) -> Result<Self, Diagnostic> {
        let file = File::open(fasta_path).map_err(|error| {
            Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                .with("path", fasta_path.display().to_string())
                .with("operation", "open_for_index")
                .with("detail", error.to_string())
        })?;
        let mut indexer = fasta::io::Indexer::new(BufReader::new(file));
        let mut records = Vec::new();
        loop {
            match indexer.index_record() {
                Ok(Some(record)) => records.push(record),
                Ok(None) => break,
                Err(error) => {
                    return Err(Diagnostic::io(
                        "GPS-IO-001",
                        "Una operación de lectura o escritura falló.",
                    )
                    .with("path", fasta_path.display().to_string())
                    .with("operation", "index_reference")
                    .with("detail", error.to_string()));
                }
            }
        }
        Ok(Self {
            reference_id,
            catalog,
            fasta_path: fasta_path.to_path_buf(),
            index: fasta::fai::Index::from(records),
        })
    }
}

impl ReductionSource for IndexedSource {
    fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    fn catalog(&self) -> &[SourceSequenceEntry] {
        &self.catalog
    }

    fn interval(
        &self,
        sequence_id: &SequenceId,
        start: usize,
        end: usize,
    ) -> Result<Vec<u8>, Diagnostic> {
        if end < start {
            return Err(out_of_range(sequence_id));
        }
        // Un intervalo vacío no es representable como región del índice, pero la
        // fuente en memoria lo acepta devolviendo un fragmento vacío. Ambas
        // deben ser indistinguibles, así que se conserva ese comportamiento.
        if end == start {
            return if self.catalog.iter().any(|entry| entry.id() == sequence_id) {
                Ok(Vec::new())
            } else {
                Err(out_of_range(sequence_id))
            };
        }
        // El FAI usa coordenadas de una base y extremos inclusivos; el mapa de
        // reducción usa un intervalo semiabierto de base cero.
        let region = format!("{}:{}-{}", sequence_id.as_str(), start + 1, end);
        // Inalcanzable el rechazo de «región ilegible»: ver
        // `the_indexed_region_text_always_parses`.
        let region = region.parse().map_err(|_| out_of_range(sequence_id))?;
        let mut reader = fasta::io::indexed_reader::Builder::default()
            .set_index(self.index.clone())
            .build_from_path(&self.fasta_path)
            .map_err(|error| {
                Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                    .with("path", self.fasta_path.display().to_string())
                    .with("operation", "open_indexed_reference")
                    .with("detail", error.to_string())
            })?;
        let record = reader
            .query(&region)
            .map_err(|_| out_of_range(sequence_id))?;
        let bases = record.sequence().as_ref().to_ascii_uppercase();
        if bases.len() != end - start {
            return Err(out_of_range(sequence_id));
        }
        Ok(bases)
    }
}

/// Elige de dónde salen las bases de una referencia según estén o no cargadas.
///
/// Si la referencia se validó sin conservarlas, las sirve por índice desde el
/// FASTA indicado; si están en memoria, las sirve desde ahí. El reconstructor no
/// distingue una de otra.
pub fn source_for<'a>(
    reference: &'a Reference,
    fasta_path: Option<&Path>,
) -> Result<Box<dyn ReductionSource + 'a>, Diagnostic> {
    let materialized = reference
        .sequences()
        .first()
        .is_none_or(|sequence| sequence.sequence().is_materialized());
    if materialized {
        return Ok(Box::new(MaterializedSource::new(reference)));
    }
    let catalog = reference
        .sequences()
        .iter()
        .map(|sequence| {
            SourceSequenceEntry::new(sequence.id().clone(), sequence.class(), sequence.length())
        })
        .collect();
    // La ruta explícita manda; si no la hay, sirve la que la referencia declara
    // como origen o como recurso autoritativo.
    let path = fasta_path
        .map(Path::to_path_buf)
        .or_else(|| reference.authoritative_resource().map(Path::to_path_buf))
        .or_else(|| match reference.origin() {
            crate::reference::ReferenceOrigin::Supplied { original_path } => {
                Some(original_path.clone())
            }
            _ => None,
        })
        .ok_or_else(|| {
            map_error(
                reference.id().as_str(),
                "source_sequence",
                "no hay FASTA del que leer las bases por índice",
            )
        })?;
    IndexedSource::new(reference.id().clone(), catalog, &path).map(|source| {
        let boxed: Box<dyn ReductionSource + 'a> = Box::new(source);
        boxed
    })
}

pub fn build_reduced_reference(
    source: &dyn ReductionSource,
    plan: &MaterializedPlan,
) -> Result<Reference, Diagnostic> {
    if source.reference_id() != plan.source_reference_id() {
        return Err(map_error(
            plan.derived_reference_id().as_str(),
            "source_reference",
            "el plan no pertenece al FASTA autoritativo",
        ));
    }
    let mut derived_sequences = Vec::new();
    let mut segments = plan.segments().iter().peekable();
    for source_sequence in source.catalog() {
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
            // Inalcanzables las dos conversiones a `usize` que siguen en la plataforma
            // objetivo: `usize` y `u64` miden lo mismo en x86_64 y la conversión no falla.
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
            let fragment = source.interval(source_sequence.id(), start, end)?;
            bases.extend_from_slice(&fragment);
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
    let source_reference_id = source.reference_id().clone();
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
    source: &dyn ReductionSource,
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
    source: &dyn ReductionSource,
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
        != Some(source.reference_id().as_str())
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
        // Inalcanzable: el cociente es como mucho `u64::MAX` y solo hay resto cuando
        // `region_size_bases > 1`, de modo que sumar uno no puede desbordar.
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
            // Inalcanzable: `index` es un índice de un vector y `index + 1` no desborda `u64`.
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
            .catalog()
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
            // Inalcanzable: `previous_derived_end` solo guarda finales de segmentos que ya
            // pasaron la comprobación de longitud (la de «región o residual»), así que están
            // acotados por la suma de longitudes y sumarles uno no desborda.
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
        // Inalcanzables este y el siguiente `usize::try_from`: ver arriba (x86_64).
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
            .extend_from_slice(&source.interval(source_sequence.id(), start, end)?);
        previous_source_end.insert(sequence_id, source_end);
        previous_derived_end.insert(sequence_id, derived_end);
    }
    // Inalcanzable la comprobación que sigue, en sus dos mitades: los órdenes de
    // selección ya son únicos (`insert`) y cada uno cae en `1..=full` o es igual al
    // último (`expected_segment_count`), de modo que `n` órdenes únicos con
    // `n == expected_segment_count` forman exactamente `1..=n`. Tampoco la conversión
    // `usize` a `u64` de su interior puede fallar.
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
            .catalog()
            .iter()
            .find(|source_sequence| source_sequence.id().as_str() == sequence.id().as_str())
            .map(SourceSequenceEntry::class);
        if source_class != Some(sequence.class()) {
            return Err(map_error(
                derived_id,
                "derived_sequence_class",
                "la clase del catálogo derivado no coincide con la fuente y el mapa",
            ));
        }
        // Inalcanzable el `else`: `rebuilt_order[index]` es este mismo identificador
        // (comprobado arriba) y todo identificador de `rebuilt_order` está en `rebuilt`.
        let Some(expected) = rebuilt.get(sequence.id().as_str()) else {
            return Err(map_error(
                derived_id,
                "derived_sequence",
                "falta una secuencia derivada en el mapa",
            ));
        };
        let observed_bases = sequence.sequence().as_bytes().map_err(|_| {
            map_error(
                derived_id,
                "derived_sequence",
                "la referencia derivada se validó sin conservar sus bases",
            )
        })?;
        // La mitad de cobertura (`previous_derived_end != length`) nunca decide
        // sola: si las bases coinciden con la reconstrucción, la longitud también,
        // porque la contigüidad ya comprobada hace que el último final derivado
        // sea la suma de las longitudes de los segmentos. Quitarla es equivalente;
        // la prueba de la mitad de bases es `a_derived_reference_that_differs_from_the_map_is_rejected`.
        if expected.as_slice() != observed_bases
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
    // Inalcanzable el rechazo: las dos comprobaciones de arriba (bases iguales
    // a la reconstrucción y cobertura derivada) y el recuento y las longitudes
    // de los segmentos ya fijan que la reconstrucción suma `target_bases`:
    // cada segmento mide `region` o el residuo y los órdenes de selección
    // forman `1..=n`.
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

/// Rechazo del intervalo fuente que queda fuera de la secuencia. Es una función
/// y no una clausura para que lleve `#[track_caller]`: cada uso es su propio
/// punto de rechazo para la puerta de cobertura.
#[track_caller]
fn out_of_range(sequence_id: &SequenceId) -> Diagnostic {
    crate::guard_trace::record();
    map_error(
        sequence_id.as_str(),
        "source_interval",
        "el intervalo fuente queda fuera de la secuencia",
    )
}

#[track_caller]
fn map_error(reference: &str, invariant: &str, detail: &str) -> Diagnostic {
    crate::guard_trace::record();
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

#[cfg(test)]
mod tests {
    //! Pruebas de las guardas del mapa de reducción que ninguna prueba
    //! alcanzaba (DT-001).
    //!
    //! El documento se escribe **a mano** para un caso pequeño y conocido, de
    //! modo que cada alteración toque un solo punto y el motivo exacto del
    //! rechazo (`invariant` y `detail`) diga qué guarda ha saltado: así una
    //! guarda rota no la tapa otra anterior (DT-032). Las comprobaciones sobre
    //! el plan usan un plan real, elegido por el selector.

    use super::*;
    use crate::randomness::{RootSeed, RuntimeSubseedScope, derive_subseed};
    use crate::reference_reduction::{
        EligibleSequence, materialize_selection, select_reduced_reference,
    };
    use std::fs;
    use std::sync::atomic::{AtomicU64, Ordering};

    static COUNTER: AtomicU64 = AtomicU64::new(0);

    const PARAMETERS: ReductionMapParameters = ReductionMapParameters {
        region_size_bases: 5,
        target_bases: 13,
    };

    fn reference_id(value: &str) -> ReferenceId {
        ReferenceId::parse(value, "test").expect("reference id")
    }

    fn sequence_id(value: &str) -> SequenceId {
        SequenceId::parse(value, "test").expect("sequence id")
    }

    fn sequence(id: &str, class: SequenceClass, bases: &[u8]) -> ReferenceSequence {
        ReferenceSequence::new(
            sequence_id(id),
            FastaHeader::parse(id.as_bytes()).expect("header"),
            class,
            &DnaSequence::try_from_bytes(bases).expect("bases"),
        )
        .expect("sequence")
    }

    const CHR1: &[u8] = b"ACGTACGTACGT";
    const SCAF: &[u8] = b"TTGGCCAA";

    fn source_reference(id: &str) -> Reference {
        Reference::new(
            reference_id(id),
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, CHR1),
                sequence("scaf", SequenceClass::Scaffold, SCAF),
            ],
        )
        .expect("source")
    }

    /// Referencia derivada de la maqueta: `chr1[0..10]` y `scaf[0..3]`.
    fn derived_reference(chr1: &[u8], scaf: &[u8]) -> Reference {
        Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, chr1),
                sequence("scaf", SequenceClass::Scaffold, scaf),
            ],
        )
        .expect("derived")
    }

    fn valid_derived() -> Reference {
        derived_reference(&CHR1[..10], &SCAF[..3])
    }

    #[allow(clippy::too_many_arguments)]
    fn segment(
        materialization: u64,
        selection: u64,
        class: &str,
        id: &str,
        sequence_class: &str,
        source: (u64, u64),
        derived: (u64, u64),
    ) -> Value {
        json!({
            "materialization_order": materialization,
            "selection_order": selection,
            "segment_class": class,
            "sequence_id": id,
            "sequence_class": sequence_class,
            "source_interval": {"start": source.0, "end": source.1},
            "derived_interval": {"start": derived.0, "end": derived.1},
            "orientation": "forward",
        })
    }

    fn document(segments: Vec<Value>) -> Value {
        json!({
            "document_kind": REDUCTION_MAP_DOCUMENT_KIND,
            "schema_version": REDUCTION_MAP_SCHEMA_VERSION,
            "source_reference_id": "ref-source",
            "derived_reference_id": "ref-derived",
            "segments": segments,
        })
    }

    fn valid_document() -> Value {
        document(vec![
            segment(1, 1, "full", "chr1", "chromosome", (1, 5), (1, 5)),
            segment(2, 2, "full", "chr1", "chromosome", (6, 10), (6, 10)),
            segment(3, 3, "residual", "scaf", "scaffold", (1, 3), (1, 3)),
        ])
    }

    fn judge(
        derived: &Reference,
        document: &Value,
        parameters: ReductionMapParameters,
    ) -> Result<(), Diagnostic> {
        let source = source_reference("ref-source");
        validate_reduction_map_document(
            &MaterializedSource::new(&source),
            derived,
            document,
            parameters,
        )
    }

    /// Exige el rechazo `GPS-FMT-007` con el `invariant` y un trozo del
    /// `detail` dados.
    #[track_caller]
    fn assert_map_rejection(
        result: Result<(), Diagnostic>,
        invariant: &str,
        detail: &str,
        label: &str,
    ) {
        let error = result.expect_err(label);
        assert_eq!(
            error.code(),
            "GPS-FMT-007",
            "{label}: {}",
            error.render_human()
        );
        assert_eq!(
            error.context().get("invariant").map(String::as_str),
            Some(invariant),
            "{label}: {}",
            error.render_human()
        );
        let observed = error.context().get("detail").expect("detail");
        assert!(
            observed.contains(detail),
            "{label}: detail {observed:?} does not contain {detail:?}"
        );
    }

    fn edited(pointer: &str, replacement: Option<Value>) -> Value {
        let mut document = valid_document();
        let (parent, key) = pointer.rsplit_once('/').expect("pointer has a parent");
        let parent = if parent.is_empty() {
            &mut document
        } else {
            document.pointer_mut(parent).expect("parent exists")
        };
        match (replacement, parent) {
            (Some(new), Value::Object(map)) => {
                map.insert(key.to_owned(), new);
            }
            (Some(new), Value::Array(items)) => items[key.parse::<usize>().expect("index")] = new,
            (None, Value::Object(map)) => {
                map.remove(key).expect("key exists");
            }
            _ => panic!("unsupported edit at {pointer}"),
        }
        document
    }

    #[test]
    fn the_hand_written_map_is_valid() {
        judge(&valid_derived(), &valid_document(), PARAMETERS).expect("valid map");
    }

    #[test]
    fn every_single_point_alteration_of_the_map_is_rejected_by_its_own_guard() {
        // (caso, documento, invariante, trozo del detalle)
        let cases: Vec<(&str, Value, &str, &str)> = vec![
            (
                "root is a list",
                json!([]),
                "document",
                "la raíz no es un objeto",
            ),
            (
                "unexpected document kind",
                edited("/document_kind", Some(json!("other"))),
                "document_kind",
                "clase de documento inesperada",
            ),
            (
                "unexpected schema version",
                edited("/schema_version", Some(json!("2"))),
                "schema_version",
                "versión de mapa inesperada",
            ),
            (
                "wrong source reference id",
                edited("/source_reference_id", Some(json!("other"))),
                "reference_ids",
                "no coinciden con las referencias",
            ),
            (
                "wrong derived reference id",
                edited("/derived_reference_id", Some(json!("other"))),
                "reference_ids",
                "no coinciden con las referencias",
            ),
            (
                "segments is not a list",
                edited("/segments", Some(json!("x"))),
                "segments",
                "falta el catálogo de segmentos",
            ),
            (
                "no segments",
                edited("/segments", Some(json!([]))),
                "segments",
                "está vacío",
            ),
            (
                "a segment is not an object",
                edited("/segments/0", Some(json!(5))),
                "segment",
                "no es un objeto JSON",
            ),
            (
                "materialization order is not a number",
                edited("/segments/0/materialization_order", Some(json!("1"))),
                "materialization_order",
                "falta un entero sin signo",
            ),
            (
                "selection order is zero",
                edited("/segments/0/selection_order", Some(json!(0))),
                "selection_order",
                "cero o está repetido",
            ),
            (
                "selection order is repeated",
                edited("/segments/1/selection_order", Some(json!(1))),
                "selection_order",
                "cero o está repetido",
            ),
            (
                "unknown segment class",
                edited("/segments/0/segment_class", Some(json!("half"))),
                "segment_class",
                "clase de segmento desconocida",
            ),
            (
                "backwards orientation",
                edited("/segments/0/orientation", Some(json!("reverse"))),
                "segment_contract",
                "orientación o cardinalidad",
            ),
            (
                "sequence id is not text",
                edited("/segments/0/sequence_id", Some(json!(5))),
                "sequence_id",
                "falta la secuencia",
            ),
            (
                "sequence missing from the source",
                edited("/segments/0/sequence_id", Some(json!("nowhere"))),
                "source_sequence",
                "secuencia inexistente",
            ),
            (
                "sequence class differs from the source",
                edited("/segments/0/sequence_class", Some(json!("scaffold"))),
                "sequence_class",
                "no coincide con el origen",
            ),
            (
                "selection order outside the effective layout",
                edited("/segments/0/selection_order", Some(json!(9))),
                "selection_order",
                "no corresponde a región y residual",
            ),
            (
                "source interval overlaps the previous one",
                edited(
                    "/segments/1/source_interval",
                    Some(json!({"start": 5, "end": 9})),
                ),
                "source_interval",
                "fuera de rango, solapado o desordenado",
            ),
            (
                "source interval past the end of the sequence",
                edited(
                    "/segments/1/source_interval",
                    Some(json!({"start": 9, "end": 13})),
                ),
                "source_interval",
                "fuera de rango, solapado o desordenado",
            ),
            (
                "materialization order is not contiguous",
                edited("/segments/1/materialization_order", Some(json!(5))),
                "materialization_order",
                "no es contiguo",
            ),
            (
                "full segment labelled residual",
                edited("/segments/0/segment_class", Some(json!("residual"))),
                "segment_class",
                "etiqueta o longitud full/residual",
            ),
            (
                "full segment with the wrong length",
                {
                    let mut document = valid_document();
                    document["segments"][0]["source_interval"] = json!({"start": 1, "end": 4});
                    document["segments"][0]["derived_interval"] = json!({"start": 1, "end": 4});
                    document
                },
                "segment_class",
                "etiqueta o longitud full/residual",
            ),
            (
                // El segundo segmento sigue al primero ya ensanchado, de modo
                // que la contigüidad cuadra y solo falla la longitud.
                "derived interval changes the length",
                {
                    let mut document = valid_document();
                    document["segments"][0]["derived_interval"] = json!({"start": 1, "end": 6});
                    document["segments"][1]["derived_interval"] = json!({"start": 7, "end": 11});
                    document
                },
                "derived_interval",
                "no es contigua o cambia la longitud",
            ),
            (
                "derived projection shifted and compensated",
                // El primer segmento se proyecta en (2,6) y el segundo en
                // (6,10): la longitud de cada uno cuadra y el final derivado
                // (10) también; solo la contigüidad lo rechaza.
                {
                    let mut document = valid_document();
                    document["segments"][0]["derived_interval"] = json!({"start": 2, "end": 6});
                    document["segments"][1]["derived_interval"] = json!({"start": 6, "end": 10});
                    document
                },
                "derived_interval",
                "no es contigua o cambia la longitud",
            ),
            (
                "unknown root field",
                edited("/extra", Some(json!(1))),
                "document_fields",
                "campos ausentes o desconocidos",
            ),
            (
                "segment missing a field",
                edited("/segments/0/orientation", None),
                "segment_fields",
                "campos ausentes o desconocidos",
            ),
            (
                "segment with a renamed field",
                {
                    let mut document = valid_document();
                    let segment = document["segments"][0].as_object_mut().expect("segment");
                    let value = segment.remove("orientation").expect("orientation");
                    segment.insert("orientationx".to_owned(), value);
                    document
                },
                "segment_fields",
                "campos ausentes o desconocidos",
            ),
            (
                "interval with an unknown field",
                edited("/segments/0/source_interval/extra", Some(json!(1))),
                "source_interval",
                "campos ausentes o desconocidos",
            ),
            (
                "source interval is not an object",
                edited("/segments/0/source_interval", Some(json!(5))),
                "source_interval",
                "falta el intervalo",
            ),
            (
                "source interval starts at zero",
                edited(
                    "/segments/0/source_interval",
                    Some(json!({"start": 0, "end": 5})),
                ),
                "source_interval",
                "vacío o empieza en cero",
            ),
            (
                "source interval ends before it starts",
                edited(
                    "/segments/0/source_interval",
                    Some(json!({"start": 5, "end": 4})),
                ),
                "source_interval",
                "vacío o empieza en cero",
            ),
        ];
        for (label, document, invariant, detail) in cases {
            assert_map_rejection(
                judge(&valid_derived(), &document, PARAMETERS),
                invariant,
                detail,
                label,
            );
        }
    }

    #[test]
    fn the_effective_parameters_must_match_the_map() {
        for (label, parameters, invariant, detail) in [
            (
                "zero region size",
                ReductionMapParameters {
                    region_size_bases: 0,
                    target_bases: 13,
                },
                "effective_parameters",
                "mayores que cero",
            ),
            (
                "zero target",
                ReductionMapParameters {
                    region_size_bases: 5,
                    target_bases: 0,
                },
                "effective_parameters",
                "mayores que cero",
            ),
            (
                "a segment count the parameters do not produce",
                ReductionMapParameters {
                    region_size_bases: 5,
                    target_bases: 18,
                },
                "segment_count",
                "no coincide con región y objetivo efectivos",
            ),
        ] {
            assert_map_rejection(
                judge(&valid_derived(), &valid_document(), parameters),
                invariant,
                detail,
                label,
            );
        }
    }

    #[test]
    fn a_segment_that_goes_back_in_the_source_catalog_is_rejected() {
        // El segundo segmento cae en `scaf` y el tercero vuelve a `chr1`.
        let document = document(vec![
            segment(1, 1, "full", "chr1", "chromosome", (1, 5), (1, 5)),
            segment(2, 2, "full", "scaf", "scaffold", (1, 5), (1, 5)),
            segment(3, 3, "residual", "chr1", "chromosome", (6, 8), (6, 8)),
        ]);
        assert_map_rejection(
            judge(&valid_derived(), &document, PARAMETERS),
            "source_order",
            "no respetan el catálogo",
            "segment back in the catalog",
        );
    }

    #[test]
    fn a_second_residual_segment_is_rejected_even_after_a_valid_first_one() {
        // Residual primero (selección 3, longitud 3) y completos después: el
        // documento es válido. Marcar además el último como residual sobrepasa
        // la cardinalidad.
        let valid = document(vec![
            segment(1, 3, "residual", "chr1", "chromosome", (1, 3), (1, 3)),
            segment(2, 1, "full", "chr1", "chromosome", (6, 10), (4, 8)),
            segment(3, 2, "full", "scaf", "scaffold", (1, 5), (1, 5)),
        ]);
        let derived = {
            let mut chr1 = CHR1[..3].to_vec();
            chr1.extend_from_slice(&CHR1[5..10]);
            derived_reference(&chr1, &SCAF[..5])
        };
        judge(&derived, &valid, PARAMETERS).expect("residual-first map is valid");
        let mut doubled = valid;
        doubled["segments"][2]["segment_class"] = json!("residual");
        assert_map_rejection(
            judge(&derived, &doubled, PARAMETERS),
            "segment_contract",
            "orientación o cardinalidad",
            "two residual segments",
        );
    }

    #[test]
    fn a_derived_reference_that_differs_from_the_map_is_rejected() {
        let chr1_changed = {
            let mut bases = CHR1[..10].to_vec();
            bases[4] = b'G';
            bases
        };
        let extra = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, &CHR1[..10]),
                sequence("scaf", SequenceClass::Scaffold, &SCAF[..3]),
                sequence("other", SequenceClass::Scaffold, b"AC"),
            ],
        )
        .expect("extra");
        let missing = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![sequence("chr1", SequenceClass::Chromosome, &CHR1[..10])],
        )
        .expect("missing");
        let without_bases = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, &CHR1[..10]),
                ReferenceSequence::from_owned(
                    sequence_id("scaf"),
                    FastaHeader::parse(b"scaf").expect("header"),
                    SequenceClass::Scaffold,
                    DnaSequence::summary(3, 0).expect("summary"),
                )
                .expect("sequence"),
            ],
        )
        .expect("summary derived");
        let scaf_as_chromosome = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![
                sequence("chr1", SequenceClass::Chromosome, &CHR1[..10]),
                sequence("scaf", SequenceClass::Chromosome, &SCAF[..3]),
            ],
        )
        .expect("class-altered derived");
        let swapped_order = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![
                sequence("scaf", SequenceClass::Scaffold, &SCAF[..3]),
                sequence("chr1", SequenceClass::Chromosome, &CHR1[..10]),
            ],
        )
        .expect("reordered derived");
        for (label, derived, invariant, detail) in [
            (
                "derived catalog in another order",
                swapped_order,
                "derived_catalog_order",
                "no coincide con el mapa materializado",
            ),
            (
                "derived class differs from the source",
                scaf_as_chromosome,
                "derived_sequence_class",
                "no coincide con la fuente y el mapa",
            ),
            (
                "derived has an extra sequence",
                extra,
                "derived_catalog",
                "omite o inventa",
            ),
            (
                "derived lacks a sequence",
                missing,
                "derived_catalog",
                "omite o inventa",
            ),
            (
                "derived bases differ",
                derived_reference(&chr1_changed, &SCAF[..3]),
                "derived_sequence",
                "bases o la cobertura",
            ),
            (
                // Más corta en una base: altera las bases y la longitud a la
                // vez, así que la caza la mitad de bases; la de cobertura nunca
                // decide sola (ver el comentario junto a la guarda).
                "derived is one base shorter (caught by the bases half)",
                derived_reference(&CHR1[..9], &SCAF[..3]),
                "derived_sequence",
                "bases o la cobertura",
            ),
            (
                "derived has no bases",
                without_bases,
                "derived_sequence",
                "sin conservar sus bases",
            ),
        ] {
            assert_map_rejection(
                judge(&derived, &valid_document(), PARAMETERS),
                invariant,
                detail,
                label,
            );
        }
    }

    // ---- fuentes de bases -------------------------------------------------

    #[test]
    fn the_materialized_source_rejects_a_sequence_it_does_not_have() {
        let source = source_reference("ref-source");
        let materialized = MaterializedSource::new(&source);
        let error = materialized
            .interval(&sequence_id("nowhere"), 0, 1)
            .expect_err("unknown sequence");
        assert_eq!(error.code(), "GPS-FMT-007");
        assert_eq!(
            error.context().get("detail").map(String::as_str),
            Some("el mapa referencia una secuencia inexistente")
        );
    }

    struct Fasta {
        directory: PathBuf,
        path: PathBuf,
    }

    impl Fasta {
        fn new(bytes: &[u8]) -> Self {
            let directory = std::env::temp_dir().join(format!(
                "gps-reduction-io-unit-{}-{}",
                std::process::id(),
                COUNTER.fetch_add(1, Ordering::Relaxed)
            ));
            fs::create_dir(&directory).expect("scratch directory");
            let path = directory.join("reference.fa");
            fs::write(&path, bytes).expect("write FASTA");
            Self { directory, path }
        }
    }

    impl Drop for Fasta {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.directory);
        }
    }

    fn catalog() -> Vec<SourceSequenceEntry> {
        vec![SourceSequenceEntry::new(
            sequence_id("chr1"),
            SequenceClass::Chromosome,
            12,
        )]
    }

    #[test]
    fn the_indexed_source_rejects_each_interval_that_it_cannot_serve() {
        let fasta = Fasta::new(b">chr1\nACGTACGTACGT\n");
        let source = IndexedSource::new(reference_id("ref-source"), catalog(), &fasta.path)
            .expect("indexed source");
        let chr1 = sequence_id("chr1");
        for (label, id, start, end) in [
            ("end before start", &chr1, 5, 3),
            (
                "empty interval on an unknown sequence",
                &sequence_id("nowhere"),
                5,
                5,
            ),
            ("interval past the sequence", &chr1, 10, 13),
        ] {
            let error = source.interval(id, start, end).expect_err(label);
            assert_eq!(error.code(), "GPS-FMT-007", "{label}");
            assert_eq!(
                error.context().get("detail").map(String::as_str),
                Some("el intervalo fuente queda fuera de la secuencia"),
                "{label}"
            );
        }
        assert_eq!(source.interval(&chr1, 2, 6).expect("valid"), b"GTAC");
        assert_eq!(
            source.interval(&chr1, 5, 5).expect("empty"),
            Vec::<u8>::new()
        );
    }

    #[test]
    fn a_reference_without_bases_needs_a_fasta_to_serve_its_intervals() {
        let summary_reference = Reference::new(
            reference_id("ref-source"),
            ReferenceOrigin::Generated,
            vec![
                ReferenceSequence::from_owned(
                    sequence_id("chr1"),
                    FastaHeader::parse(b"chr1").expect("header"),
                    SequenceClass::Chromosome,
                    DnaSequence::summary(12, 0).expect("summary"),
                )
                .expect("sequence"),
            ],
        )
        .expect("reference");
        // Servir un tramo desde memoria cuando las bases no se conservaron.
        let error = MaterializedSource::new(&summary_reference)
            .interval(&sequence_id("chr1"), 0, 3)
            .expect_err("bases not kept");
        assert_eq!(error.code(), "GPS-FMT-007");
        assert_eq!(
            error.context().get("detail").map(String::as_str),
            Some("la referencia se validó sin conservar sus bases")
        );
        let error = source_for(&summary_reference, None)
            .err()
            .expect("no FASTA to read from");
        assert_eq!(error.code(), "GPS-FMT-007");
        assert_eq!(
            error.context().get("detail").map(String::as_str),
            Some("no hay FASTA del que leer las bases por índice")
        );
    }

    // ---- reconstrucción desde un plan real --------------------------------

    fn plan_over(
        sequences: &[EligibleSequence],
        target: u64,
    ) -> crate::reference_reduction::MaterializedPlan {
        let mut rng = derive_subseed(
            RootSeed::from_bytes([7; 32]),
            RuntimeSubseedScope::ReferenceReduction {
                source_reference_id: "ref-source",
                derived_reference_id: "ref-derived",
            },
        )
        .expect("subseed")
        .local_rng();
        let selection =
            select_reduced_reference(&reference_id("ref-source"), sequences, target, 5, &mut rng)
                .expect("selection");
        materialize_selection(&selection, &reference_id("ref-derived"), sequences).expect("plan")
    }

    fn eligible(id: &str, class: SequenceClass, length: u64) -> EligibleSequence {
        EligibleSequence::new(sequence_id(id), class, length).expect("eligible")
    }

    fn two_chromosome_source(id: &str, second: Option<SequenceClass>) -> Reference {
        let mut sequences = vec![sequence("chr1", SequenceClass::Chromosome, b"ACGTACGTAC")];
        if let Some(class) = second {
            sequences.push(sequence("chr2", class, b"TTGGCCAATT"));
        }
        Reference::new(reference_id(id), ReferenceOrigin::Generated, sequences).expect("source")
    }

    #[test]
    fn the_reconstruction_rejects_a_source_that_does_not_match_the_plan() {
        let eligible_pair = [
            eligible("chr1", SequenceClass::Chromosome, 10),
            eligible("chr2", SequenceClass::Chromosome, 10),
        ];
        let plan = plan_over(&eligible_pair, 15);
        let touched: BTreeSet<_> = plan
            .segments()
            .iter()
            .map(|segment| segment.source_interval().sequence_id().as_str().to_owned())
            .collect();
        assert_eq!(touched.len(), 2, "the fixture plan must use both sequences");

        let matching = two_chromosome_source("ref-source", Some(SequenceClass::Chromosome));
        build_reduced_reference(&MaterializedSource::new(&matching), &plan)
            .expect("a matching source rebuilds");

        // Fuente con otro identificador de referencia.
        let other_id = two_chromosome_source("ref-other", Some(SequenceClass::Chromosome));
        let error = build_reduced_reference(&MaterializedSource::new(&other_id), &plan)
            .expect_err("different reference");
        assert_eq!(error.context()["invariant"], "source_reference");

        // Fuente que no tiene una de las secuencias del plan: sobran segmentos.
        let short = two_chromosome_source("ref-source", None);
        let error = build_reduced_reference(&MaterializedSource::new(&short), &plan)
            .expect_err("a sequence of the plan is missing");
        assert_eq!(error.context()["invariant"], "source_sequence");
        assert!(error.context()["detail"].contains("inexistente o fuera de orden"));

        // Fuente cuya segunda secuencia es un scaffold donde el plan puso un
        // cromosoma.
        let other_class = two_chromosome_source("ref-source", Some(SequenceClass::Scaffold));
        let error = build_reduced_reference(&MaterializedSource::new(&other_class), &plan)
            .expect_err("class differs");
        assert_eq!(error.context()["invariant"], "sequence_class");
    }

    #[test]
    fn a_derived_reference_is_judged_against_the_rebuild_of_its_plan() {
        let eligible_pair = [
            eligible("chr1", SequenceClass::Chromosome, 10),
            eligible("chr2", SequenceClass::Chromosome, 10),
        ];
        let plan = plan_over(&eligible_pair, 15);
        let source = two_chromosome_source("ref-source", Some(SequenceClass::Chromosome));
        let materialized = MaterializedSource::new(&source);
        let derived = build_reduced_reference(&materialized, &plan).expect("derived");
        validate_reduced_reference(&materialized, &derived, &plan).expect("faithful derived");

        let rebuilt_with = |id: &str, class: SequenceClass, flip_first_base: bool| {
            let sequences = derived
                .sequences()
                .iter()
                .enumerate()
                .map(|(index, original)| {
                    let mut bases = original.sequence().as_bytes().expect("bases").to_vec();
                    if flip_first_base && index == 0 {
                        bases[0] = if bases[0] == b'A' { b'C' } else { b'A' };
                    }
                    let class = if index == 0 { class } else { original.class() };
                    sequence(original.id().as_str(), class, &bases)
                })
                .collect();
            Reference::new(reference_id(id), ReferenceOrigin::Generated, sequences)
                .expect("altered derived")
        };

        let other_id = rebuilt_with("ref-elsewhere", SequenceClass::Chromosome, false);
        let error = validate_reduced_reference(&materialized, &other_id, &plan)
            .expect_err("different derived id");
        assert_eq!(error.context()["invariant"], "derived_catalog");

        // Un catálogo derivado con una secuencia menos que la reconstrucción.
        let first = &derived.sequences()[0];
        let fewer = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            vec![sequence(
                first.id().as_str(),
                first.class(),
                first.sequence().as_bytes().expect("bases"),
            )],
        )
        .expect("shorter derived");
        let error =
            validate_reduced_reference(&materialized, &fewer, &plan).expect_err("fewer sequences");
        assert_eq!(error.context()["invariant"], "derived_catalog");
        // Mismo número de secuencias, la primera con otro identificador.
        let renamed = Reference::new(
            reference_id("ref-derived"),
            ReferenceOrigin::Generated,
            derived
                .sequences()
                .iter()
                .enumerate()
                .map(|(index, original)| {
                    let id = if index == 0 {
                        "chrX"
                    } else {
                        original.id().as_str()
                    };
                    sequence(
                        id,
                        original.class(),
                        original.sequence().as_bytes().expect("bases"),
                    )
                })
                .collect(),
        )
        .expect("renamed derived");
        let error =
            validate_reduced_reference(&materialized, &renamed, &plan).expect_err("another id");
        assert_eq!(error.context()["invariant"], "derived_sequence");

        for (label, altered) in [
            (
                "another base",
                rebuilt_with("ref-derived", SequenceClass::Chromosome, true),
            ),
            (
                "another class",
                rebuilt_with("ref-derived", SequenceClass::Scaffold, false),
            ),
        ] {
            let error =
                validate_reduced_reference(&materialized, &altered, &plan).expect_err(label);
            assert_eq!(error.context()["invariant"], "derived_sequence", "{label}");
        }
    }

    /// `IndexedSource::interval` no puede rechazar por «región ilegible»: la
    /// región es `{id}:{start+1}-{end}` con un identificador de la sintaxis
    /// `[A-Za-z0-9][A-Za-z0-9._-]{0,63}` y dos decimales con `start < end`, y
    /// noodles acepta cualquier nombre de esa forma, incluso con guiones y
    /// cifras que se parecen a un intervalo. La prueba lo ejecuta con los
    /// nombres más sospechosos y con las coordenadas más grandes; no
    /// discrimina esa guarda (no se alcanza).
    #[test]
    fn the_indexed_region_text_always_parses() {
        for id in [
            "chr1", "1-5", "a-1-5", "a.1-5", "x_1", "0", "a-", "9-9", "chr1-1-2",
        ] {
            for (start, end) in [(0_usize, 1_usize), (3, 9), (usize::MAX - 1, usize::MAX)] {
                let region = format!("{id}:{}-{}", start + 1, end);
                assert!(region.parse::<noodles_core::Region>().is_ok(), "{region}");
            }
        }
    }
}
