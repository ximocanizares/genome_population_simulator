//! Lector, validador de esquema y reconstructor independiente de `truth.jsonl`.
//!
//! Este módulo no llama a `truth_projection`: convierte el contrato JSONL a
//! tipos externos, valida su estructura y materializa la vista vigente desde
//! las piezas y la referencia autoritativa. Tampoco acepta ni consulta VCF.

use crate::guard_trace::rejected;
use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::io::BufRead;

use serde::de::{DeserializeSeed, MapAccess, SeqAccess, Unexpected, Visitor};
use serde::{Deserialize, Deserializer, de::Error as _};
use serde_json::Number;
use sha2::{Digest, Sha256};

use crate::population_snvs::PreparedPopulationSnvs;
use crate::reference::{Reference, ReferenceSequence};
use crate::truth_schema::{
    AlleleCarrier, CopyIdentity, CopyProvenance, GenotypeCaseKind, GenotypeRecord,
    HaplotypeCopyRecord, MatrixProvenance, ModificationProvenance, ModificationProvenanceKind,
    NamedDigest, SegmentSource, TruthAllele, TruthContractError, TruthCurrentPiece,
    TruthCurrentSnp, TruthDocument, TruthHeader, TruthHistoryEvent, TruthLocus, TruthModification,
    TruthPieceLink, TruthRecord, TruthRepeatBlock, TruthSourceKind, TruthSourceRef,
    TruthStructuralEvent, VersionRule, validate_lexical_genotype, validate_schema_version,
};

/// Resultado reconstruido exclusivamente desde `truth.jsonl` y la referencia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ReconstructedTruth {
    pub document: TruthDocument,
    pub copies: BTreeMap<CopyIdentity, Vec<u8>>,
}

/// Vista compacta de una verdad que ya ha superado la validación semántica
/// completa necesaria para derivar VCF.
///
/// La verdad mantiene piezas fuente, SNPs dispersos, inserciones y genotipos;
/// las bases completas reconstruidas viven únicamente durante la validación de
/// cada copia. El campo es privado para que los consumidores VCF no puedan
/// omitir la frontera de validación.
#[derive(Clone, Copy, Debug)]
pub(crate) struct ValidatedTruthForVcf<'truth> {
    document: &'truth TruthDocument,
    reference_id: &'truth str,
    reference_content_sha256: &'truth str,
}

impl<'truth> ValidatedTruthForVcf<'truth> {
    #[must_use]
    pub(crate) fn document(self) -> &'truth TruthDocument {
        self.document
    }

    /// Rechaza que una vista aprobada para una referencia se reutilice contra
    /// otra en los consumidores internos de VCF.
    pub(crate) fn validate_reference(
        self,
        reference: &Reference,
    ) -> Result<(), TruthValidationError> {
        if self.reference_id != reference.id().as_str() {
            return Err(rejected(TruthValidationError::ReferenceIdMismatch {
                declared: self.reference_id.to_owned(),
                actual: reference.id().as_str().to_owned(),
            }));
        }
        let observed = match reference_content_sha256(reference) {
            Ok(observed) => observed,
            Err(error) => return Err(rejected(error)),
        };
        if self.reference_content_sha256 != observed {
            return Err(rejected(TruthValidationError::ReferenceDigestMismatch));
        }
        Ok(())
    }
}

/// Recuento obtenido al validar por flujo. No conserva el documento ni las
/// secuencias de las copias ya entregadas al consumidor.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TruthStreamSummary {
    pub copy_count: u64,
    pub genotype_count: u64,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthValidationError {
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
    Io {
        line: u64,
        detail: String,
    },
    EmptyDocument,
    MissingFinalNewline,
    Schema {
        line: u64,
        detail: String,
    },
    HeaderMustBeFirst,
    DuplicateHeader,
    Contract(TruthContractError),
    ReferenceIdMismatch {
        declared: String,
        actual: String,
    },
    ReferenceDigestMismatch,
    MissingOrDuplicateReferenceDigest,
    MissingReferenceSequence {
        sequence_id: String,
    },
    SourceDoesNotNameReference {
        identity: CopyIdentity,
    },
    SourceLengthMismatch {
        sequence_id: String,
        declared: u64,
        actual: u64,
    },
    SourceIntervalOutOfBounds {
        identity: CopyIdentity,
    },
    SnpReferenceMismatch {
        event_id: String,
        position: u64,
    },
    HistoryDoesNotReplay {
        identity: CopyIdentity,
    },
    DuplicateCurrentSourcePosition {
        identity: CopyIdentity,
        position: u64,
    },
    /// El estado heredado de la copia —sus segmentos, antes del historial— no
    /// recorre su secuencia de la referencia exactamente una vez y en orden
    /// creciente (DT-048). El recorrido había llegado hasta `reached` y el
    /// segmento siguiente empieza en `next_start`; tras el último segmento,
    /// `next_start` es la longitud de la secuencia. Un hueco, un solape, una
    /// duplicación, una rotación o un desorden lo rompen en el primer segmento
    /// que no empieza donde iba el recorrido. Con `next_start == reached`, el
    /// segmento empieza donde debía pero está vacío, invertido o es de otra
    /// secuencia (DT-050).
    SegmentsNotCollinear {
        identity: CopyIdentity,
        reached: u64,
        next_start: u64,
    },
    GenotypeMismatch {
        case_id: String,
        detail: &'static str,
    },
    /// Falta un caso que la verdad tiene que publicar: un cambio de una copia
    /// con la referencia que ningún caso explica, un evento estructural del
    /// historial sin su caso o una fila de la matriz sin el caso del individuo
    /// (DT-020 y DT-021). El intervalo está en coordenadas de la referencia,
    /// salvo para unas bases insertadas sin evento, que no tienen coordenada
    /// de referencia: entonces es su tramo en la copia.
    /// Una copia declara en su procedencia (`provenance.origin`) una
    /// población que no es la de su individuo.
    CopyOriginMismatch {
        identity: CopyIdentity,
        expected: String,
    },
    MissingCase {
        individual_id: String,
        sequence_id: String,
        start: u64,
        end: u64,
        detail: &'static str,
    },
}

impl std::fmt::Display for TruthValidationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

/// Bases de una secuencia, o el error propio si no se materializaron.
///
/// **Inalcanzable** (DT-001, tanda 5): el rechazo de `BasesNotMaterialized` de esta función no lo
/// alcanza ninguna entrada. Todos los llamantes públicos contrastan primero la
/// huella de la referencia (`validate_reference_digest`), que se calcula sobre
/// las bases de **cada** secuencia: con una sola sin bases ya devuelve
/// `BasesNotMaterialized` desde allí, antes de reconstruir ninguna copia ni
/// juzgar ningún caso. Con las bases materializadas, `as_bytes` no falla.
/// Intento ejecutado: `a_reference_without_bases_is_rejected_before_any_copy_is_rebuilt`
/// (no discrimina esta guarda; muestra que el rechazo es el de la huella). Lo
/// mismo vale para los otros dos sitios que leen bases (la réplica de un cambio
/// y la materialización).
fn bases_of(sequence: &crate::reference::ReferenceSequence) -> Result<&[u8], TruthValidationError> {
    sequence
        .sequence()
        .as_bytes()
        .map_err(|_| rejected(TruthValidationError::BasesNotMaterialized))
}

impl std::error::Error for TruthValidationError {}

/// Lee y valida cada línea frente al modelo cerrado del esquema v1.
///
/// Los tipos de entrada tienen `deny_unknown_fields`, discriminadores cerrados
/// y deserialización de estructuras (que rechaza claves duplicadas), de modo
/// que esta frontera no interpreta JSON parcial antes de validar su forma.
pub fn read_truth_jsonl(reader: &mut impl BufRead) -> Result<TruthDocument, TruthValidationError> {
    let mut raw = Vec::new();
    let mut header = None;
    let mut records = Vec::new();
    let mut line_number = 0_u64;
    loop {
        if !read_jsonl_line(reader, &mut raw, line_number.saturating_add(1))? {
            break;
        }
        line_number = line_number.saturating_add(1);
        if raw.last() != Some(&b'\n') {
            return Err(rejected(TruthValidationError::MissingFinalNewline));
        }
        let content = &raw[..raw.len() - 1];
        let content = content.strip_suffix(b"\r").unwrap_or(content);
        let line = std::str::from_utf8(content).map_err(|error| {
            rejected(TruthValidationError::Schema {
                line: line_number,
                detail: error.to_string(),
            })
        })?;
        if line.is_empty() {
            return Err(rejected(TruthValidationError::Schema {
                line: line_number,
                detail: "línea JSONL vacía".to_owned(),
            }));
        }
        let value = parse_wire_text(line, line_number)?;
        match value {
            WireLine::Header(value) => {
                if line_number != 1 {
                    return Err(rejected(TruthValidationError::HeaderMustBeFirst));
                }
                // **Inalcanzable** (DT-001, tanda 5): una cabecera solo se acepta en la línea 1 (la guarda de
                // arriba) y el número de línea crece en cada vuelta, de modo que `header`
                // siempre está vacío aquí. Una segunda cabecera da `HeaderMustBeFirst`
                // (`a_header_on_line_two_is_rejected_as_not_first`).
                if header.replace(value.into()).is_some() {
                    return Err(rejected(TruthValidationError::DuplicateHeader));
                }
            }
            WireLine::HaplotypeCopy(value) => {
                if header.is_none() {
                    return Err(rejected(TruthValidationError::HeaderMustBeFirst));
                }
                records.push(TruthRecord::HaplotypeCopy(value.into()));
            }
            WireLine::Genotype(value) => {
                if header.is_none() {
                    return Err(rejected(TruthValidationError::HeaderMustBeFirst));
                }
                records.push(TruthRecord::Genotype(value.into()));
            }
        }
    }
    let document = TruthDocument {
        header: header.ok_or_else(|| rejected(TruthValidationError::EmptyDocument))?,
        records,
    };
    document
        .validate_structure()
        .map_err(TruthValidationError::Contract)?;
    Ok(document)
}

/// Lee una sola entrada JSONL. La memoria queda acotada a un único registro:
/// el contrato no impone un máximo de tamaño a un alelo final y por ello este
/// lector no rechaza una línea válida por una cota local arbitraria.
fn read_jsonl_line(
    reader: &mut impl BufRead,
    destination: &mut Vec<u8>,
    line: u64,
) -> Result<bool, TruthValidationError> {
    destination.clear();
    let count = reader.read_until(b'\n', destination).map_err(|error| {
        rejected(TruthValidationError::Io {
            line,
            detail: error.to_string(),
        })
    })?;
    if count == 0 {
        return Ok(false);
    }
    Ok(true)
}

fn parse_wire_line(raw: &[u8], line_number: u64) -> Result<WireLine, TruthValidationError> {
    if raw.last() != Some(&b'\n') {
        return Err(rejected(TruthValidationError::MissingFinalNewline));
    }
    let content = &raw[..raw.len() - 1];
    let content = content.strip_suffix(b"\r").unwrap_or(content);
    if content.is_empty() {
        return Err(rejected(TruthValidationError::Schema {
            line: line_number,
            detail: "línea JSONL vacía".to_owned(),
        }));
    }
    let line = std::str::from_utf8(content).map_err(|error| {
        rejected(TruthValidationError::Schema {
            line: line_number,
            detail: error.to_string(),
        })
    })?;
    parse_wire_text(line, line_number)
}

/// Lee una línea de texto: primero comprueba que ninguna estructura venga
/// escrita como lista (DT-032) y después la deserializa a su tipo.
fn parse_wire_text(line: &str, line_number: u64) -> Result<WireLine, TruthValidationError> {
    let schema_error = |error: serde_json::Error| {
        rejected(TruthValidationError::Schema {
            line: line_number,
            detail: error.to_string(),
        })
    };
    ObjectShape {
        list_allowed: false,
    }
    .deserialize(&mut serde_json::Deserializer::from_str(line))
    .map_err(schema_error)?;
    serde_json::from_str(line).map_err(schema_error)
}

/// SHA-256 canónico del contenido de `Reference`: etiqueta de versión, id de
/// referencia y, en orden, id, longitud y bases de cada secuencia. No depende
/// de rutas, encabezados FASTA ni de un digest suministrado por el llamador.
/// Huella del contenido de una referencia materializada.
///
/// Falla si se validó sin conservar las bases: la huella se calcula sobre
/// ellas y no puede derivarse de la ficha.
pub fn reference_content_sha256(reference: &Reference) -> Result<String, TruthValidationError> {
    let mut digest = Sha256::new();
    digest.update(b"genome_population_simulator/reference-content/v1\0");
    digest.update(
        u64::try_from(reference.id().as_str().len())
            .expect("reference identifier length fits u64")
            .to_be_bytes(),
    );
    digest.update(reference.id().as_str().as_bytes());
    for sequence in reference.sequences() {
        digest.update(
            u64::try_from(sequence.id().as_str().len())
                .expect("sequence identifier length fits u64")
                .to_be_bytes(),
        );
        digest.update(sequence.id().as_str().as_bytes());
        digest.update(sequence.length().to_be_bytes());
        digest.update(
            sequence
                .sequence()
                .as_bytes()
                .map_err(|_| TruthValidationError::BasesNotMaterialized)?,
        );
    }
    Ok(format!("{:x}", digest.finalize()))
}

fn validate_reference_digest(
    header: &TruthHeader,
    reference: &Reference,
) -> Result<(), TruthValidationError> {
    let matches = header
        .source_digests
        .iter()
        .filter(|digest| {
            digest.source.kind == TruthSourceKind::Reference
                && digest.source.source_id == header.reference_id
        })
        .collect::<Vec<_>>();
    if matches.len() != 1 {
        return Err(rejected(
            TruthValidationError::MissingOrDuplicateReferenceDigest,
        ));
    }
    let observed = match reference_content_sha256(reference) {
        Ok(observed) => observed,
        Err(error) => return Err(rejected(error)),
    };
    if header.reference_content_sha256 != observed {
        return Err(rejected(TruthValidationError::ReferenceDigestMismatch));
    }
    Ok(())
}

/// Reconstruye todas las copias y comprueba genotipos desde la verdad leída.
pub fn validate_and_reconstruct_truth_jsonl(
    reader: &mut impl BufRead,
    reference: &Reference,
) -> Result<ReconstructedTruth, TruthValidationError> {
    let document = read_truth_jsonl(reader)?;
    validate_and_reconstruct_truth(document, reference)
}

/// Valida una verdad ya deserializada sin retener secuencias reconstruidas.
///
/// La publicación VCF comprueba el contrato estructural antes de invocar esta
/// fase, que valida referencia y semántica completa —incluido cada genotipo—
/// y solo entonces entrega una vista que permite materializar ventanas
/// puntuales.
pub(crate) fn validate_truth_for_vcf<'truth>(
    document: &'truth TruthDocument,
    reference: &Reference,
) -> Result<ValidatedTruthForVcf<'truth>, TruthValidationError> {
    // Una sola implementación del juez: la de la publicación, con una
    // inspección que no hace nada, para que las dos no puedan divergir.
    // El conjunto de casos no se juzga aquí: lo juzga, con la matriz, la
    // publicación de la verdad, que corre siempre antes que cualquier VCF.
    // Estas fronteras solo reciben el documento y no conocen las filas de la
    // matriz.
    validate_truth_with_copy_inspection(document, reference, CaseSetCheck::Skip, |_, _| {
        Ok::<(), std::convert::Infallible>(())
    })
    .map_err(|error| match error {
        InspectedValidationError::CopyReconstruction { error, .. }
        | InspectedValidationError::Semantics(error) => *error,
        InspectedValidationError::Inspection(never) => match never {},
    })?;
    Ok(ValidatedTruthForVcf {
        document,
        reference_id: &document.header.reference_id,
        reference_content_sha256: &document.header.reference_content_sha256,
    })
}

/// El juez completo de la verdad sobre una verdad **ya escrita**, con las
/// filas de la matriz y la población de cada individuo: cada caso y el
/// conjunto de casos, con los todo-REF de cada fila.
///
/// Lo usa el juez de los recursos publicados, que tiene la matriz publicada
/// y la configuración efectiva, para juzgar los bytes que de verdad quedaron
/// en el paquete y no solo el documento en memoria del que salieron.
pub(crate) fn validate_published_truth(
    document: &TruthDocument,
    reference: &Reference,
    matrix: &MatrixCaseLoci,
) -> Result<(), TruthValidationError> {
    validate_truth_with_copy_inspection(
        document,
        reference,
        CaseSetCheck::Judge(Some(matrix)),
        |_, _| Ok::<(), std::convert::Infallible>(()),
    )
    .map_err(|error| match error {
        InspectedValidationError::CopyReconstruction { error, .. }
        | InspectedValidationError::Semantics(error) => *error,
        InspectedValidationError::Inspection(never) => match never {},
    })
}

/// Fallo de [`validate_truth_with_copy_inspection`], con la fase en que ocurrió.
#[derive(Debug)]
pub(crate) enum InspectedValidationError<E> {
    /// No se pudo reconstruir la copia con esa identidad.
    CopyReconstruction {
        identity: CopyIdentity,
        error: Box<TruthValidationError>,
    },
    /// Falló la referencia, o un caso genotípico no concuerda con las copias.
    Semantics(Box<TruthValidationError>),
    /// Falló la inspección que el llamador hace de cada copia reconstruida.
    Inspection(E),
}

/// El juez semántico de la verdad —[`validate_truth_for_vcf`] es este mismo con
/// una inspección vacía—. Entrega a `inspect` la identidad y la huella SHA-256
/// de cada copia reconstruida, y suelta las bases **antes** de llamarla.
///
/// Existe para que la publicación contraste las huellas de las copias y juzgue
/// los casos genotípicos con **una sola** reconstrucción por copia, en lugar de
/// reconstruir una vez para cada cosa. Las bases no se retienen: solo los
/// resúmenes compactos que el contraste de genotipos necesita.
///
/// Con [`CaseSetCheck::Judge`] juzga además el **conjunto** de casos (DT-020
/// y DT-021): que no falte ni sobre ninguno.
pub(crate) fn validate_truth_with_copy_inspection<E>(
    document: &TruthDocument,
    reference: &Reference,
    case_set: CaseSetCheck<'_>,
    mut inspect: impl FnMut(&HaplotypeCopyRecord, [u8; 32]) -> Result<(), E>,
) -> Result<(), InspectedValidationError<E>> {
    if document.header.reference_id != reference.id().as_str() {
        return Err(rejected(InspectedValidationError::Semantics(Box::new(
            TruthValidationError::ReferenceIdMismatch {
                declared: document.header.reference_id.clone(),
                actual: reference.id().as_str().to_owned(),
            },
        ))));
    }
    validate_reference_digest(&document.header, reference)
        .map_err(|error| InspectedValidationError::Semantics(Box::new(error)))?;

    let mut details = BTreeMap::new();
    for record in &document.records {
        if let TruthRecord::HaplotypeCopy(copy) = record {
            let result = reconstruct_copy(copy, reference).map_err(|error| {
                InspectedValidationError::CopyReconstruction {
                    identity: copy.identity.clone(),
                    error: Box::new(error),
                }
            })?;
            let ReconstructedCopy {
                bases,
                source_slices,
                pieces,
                event_loci,
                inherited_loss,
            } = result;
            // La huella se calcula y las bases se sueltan **antes** de la
            // inspección: esta materializa la copia del productor, y con las
            // dos secuencias vivas a la vez un cromosoma de 100 Mb pasaba de
            // unos 204 MB de memoria máxima a unos 300 MB.
            let digest: [u8; 32] = Sha256::digest(&bases).into();
            drop(bases);
            inspect(copy, digest).map_err(InspectedValidationError::Inspection)?;
            details.insert(
                copy.identity.clone(),
                CopySummary {
                    source_slices,
                    pieces,
                    repeat_history: repeat_history_by_interval(&copy.history),
                    history: Cow::Borrowed(&copy.history),
                    event_loci,
                    inherited_loss,
                    origin: Cow::Borrowed(&copy.provenance.origin),
                },
            );
        }
    }
    for record in &document.records {
        if let TruthRecord::Genotype(genotype) = record {
            validate_genotype_semantics(genotype, reference, &details)
                .map_err(|error| InspectedValidationError::Semantics(Box::new(error)))?;
        }
    }
    if let CaseSetCheck::Judge(matrix) = case_set {
        validate_case_set(document_genotypes(document), &details, matrix)
            .map_err(|error| InspectedValidationError::Semantics(Box::new(error)))?;
    }
    Ok(())
}

/// Si el juez comprueba también el **conjunto** de casos —que no falte ni
/// sobre ninguno— además de cada caso.
#[derive(Clone, Copy, Debug)]
pub(crate) enum CaseSetCheck<'a> {
    /// Solo cada caso. Es lo que hacen las fronteras del VCF, que solo
    /// reciben el documento y corren después de que la publicación de la
    /// verdad haya juzgado ya el conjunto.
    Skip,
    /// También el conjunto. Con los loci de la matriz, exige además el caso
    /// todo-REF de cada fila y ata a su fila todo caso SNP que no sea cambio
    /// de una copia; sin ellos, ese caso se acepta sin poder comprobarlo.
    Judge(Option<&'a MatrixCaseLoci>),
}

fn document_genotypes(document: &TruthDocument) -> impl Iterator<Item = &GenotypeRecord> {
    document.records.iter().filter_map(|record| match record {
        TruthRecord::Genotype(genotype) => Some(genotype),
        TruthRecord::HaplotypeCopy(_) => None,
    })
}

/// Valida y reconstruye por flujo, reteniendo únicamente metadatos compactos
/// por copia para contrastar los genotipos que aparecen después en el orden
/// canónico. Cada secuencia materializada se entrega al callback y se libera.
pub fn validate_truth_jsonl_streaming(
    reader: &mut impl BufRead,
    reference: &Reference,
    mut on_copy: impl FnMut(&CopyIdentity, &[u8]),
) -> Result<TruthStreamSummary, TruthValidationError> {
    let mut raw = Vec::new();
    let mut line_number = 0_u64;
    let mut header = None;
    let mut carries_repeats = false;
    let mut copy_count = 0_u64;
    let mut genotype_count = 0_u64;
    let mut previous_copy = None;
    let mut previous_genotype = None;
    let mut summaries = BTreeMap::<CopyIdentity, CopySummary>::new();
    let mut copy_identities = BTreeSet::new();
    let mut seen_genotypes = BTreeSet::new();
    // El conjunto de casos se juzga individuo a individuo: el orden canónico
    // trae juntos los casos de cada uno, y solo se retienen los del actual.
    let mut individual_cases = Vec::<GenotypeRecord>::new();
    let mut judged_individuals = BTreeSet::<String>::new();
    let judge_individual =
        |cases: &[GenotypeRecord],
         individual_id: &str,
         summaries: &BTreeMap<CopyIdentity, CopySummary<'_>>| {
            // Las copias de un individuo son contiguas en el orden de las
            // identidades: se toman por rango, sin recorrer las de todos.
            let first = CopyIdentity {
                reference_id: reference.id().as_str().to_owned(),
                individual_id: individual_id.to_owned(),
                sequence_id: String::new(),
                copy_index: 0,
            };
            let copies = summaries
                .range(first..)
                .take_while(|(identity, _)| identity.individual_id == individual_id)
                .collect::<Vec<_>>();
            validate_individual_case_set(
                individual_id,
                &cases.iter().collect::<Vec<_>>(),
                &copies,
                None,
            )
        };
    while read_jsonl_line(reader, &mut raw, line_number.saturating_add(1))? {
        line_number = line_number.saturating_add(1);
        let value = parse_wire_line(&raw, line_number)?;
        match value {
            WireLine::Header(value) => {
                if line_number != 1 {
                    return Err(rejected(TruthValidationError::HeaderMustBeFirst));
                }
                let parsed: TruthHeader = value.into();
                let mut isolated = parsed.clone();
                isolated.copy_count = 0;
                isolated.genotype_count = 0;
                TruthDocument {
                    header: isolated,
                    records: Vec::new(),
                }
                .validate_structure_with_version(VersionRule::AnyKnown)
                .map_err(TruthValidationError::Contract)?;
                if parsed.reference_id != reference.id().as_str() {
                    return Err(rejected(TruthValidationError::ReferenceIdMismatch {
                        declared: parsed.reference_id,
                        actual: reference.id().as_str().to_owned(),
                    }));
                }
                validate_reference_digest(&parsed, reference)?;
                // **Inalcanzable** (DT-001, tanda 5): igual que en el lector completo; la cabecera solo se
                // acepta en la línea 1 y `header` siempre está vacío aquí.
                if header.replace(parsed).is_some() {
                    return Err(rejected(TruthValidationError::DuplicateHeader));
                }
            }
            WireLine::HaplotypeCopy(value) => {
                let header = header
                    .as_ref()
                    .ok_or_else(|| rejected(TruthValidationError::HeaderMustBeFirst))?;
                let copy: HaplotypeCopyRecord = value.into();
                // **Equivalente** (DT-001, tanda 5) solo `summaries.contains_key`:
                // `previous >= identity` ya rechaza toda identidad repetida (la anterior es
                // siempre la mayor de las vistas, así que toda identidad ya vista es menor o
                // igual que ella). `previous >= identity` **no** es redundante: anularla
                // sola muere por `two_copies_in_decreasing_order_are_out_of_order`, y la
                // repetida la cubre `a_repeated_copy_is_out_of_order_for_the_stream_and_duplicate_for_the_whole`.
                if previous_genotype.is_some()
                    || previous_copy
                        .as_ref()
                        .is_some_and(|previous| previous >= &copy.identity)
                    || summaries.contains_key(&copy.identity)
                {
                    return Err(rejected(TruthValidationError::Contract(
                        TruthContractError::RecordsOutOfOrder,
                    )));
                }
                let mut isolated_header = header.clone();
                isolated_header.copy_count = 1;
                isolated_header.genotype_count = 0;
                let isolated = TruthDocument {
                    header: isolated_header,
                    records: vec![TruthRecord::HaplotypeCopy(copy.clone())],
                };
                isolated
                    .validate_structure_with_version(VersionRule::AnyKnown)
                    .map_err(TruthValidationError::Contract)?;
                // Un alelo de repetición del historial cuenta aunque una deleción
                // posterior lo borre de las piezas vigentes y sin que se exija caso
                // de repetición (`records_carry_repeat_alleles` lo ve igual).
                carries_repeats |= isolated.carries_repeat_alleles();
                let result = reconstruct_copy(&copy, reference)?;
                on_copy(&copy.identity, &result.bases);
                let copy_identity = copy.identity.clone();
                summaries.insert(
                    copy_identity.clone(),
                    CopySummary::from_reconstructed(
                        result,
                        Cow::Owned(copy.history),
                        Cow::Owned(copy.provenance.origin),
                    ),
                );
                copy_identities.insert(copy_identity.clone());
                previous_copy = Some(copy_identity);
                // **Inalcanzable** (DT-001, tanda 5): el contador de copias es un `u64` y desbordarlo exige 2^64
                // registros, que ningún fichero ni memoria alcanza.
                copy_count = copy_count.checked_add(1).ok_or_else(|| {
                    rejected(TruthValidationError::Contract(
                        TruthContractError::CopyCountMismatch {
                            declared: 0,
                            actual: u64::MAX,
                        },
                    ))
                })?;
            }
            WireLine::Genotype(value) => {
                let header = header
                    .as_ref()
                    .ok_or_else(|| rejected(TruthValidationError::HeaderMustBeFirst))?;
                let genotype: GenotypeRecord = value.into();
                let key = (
                    genotype.reference_id.clone(),
                    genotype.individual_id.clone(),
                    genotype.locus.clone(),
                    genotype.case_id.clone(),
                );
                // **Equivalente** (DT-001, tanda 5) solo `seen_genotypes.insert`, por la
                // misma razón que `contains_key` en las copias: `previous >= key` ya rechaza
                // todo caso repetido. `previous >= key` **no** es redundante: anularla sola
                // muere por `two_genotypes_in_decreasing_order_are_out_of_order`.
                if previous_genotype
                    .as_ref()
                    .is_some_and(|previous| previous >= &key)
                    || !seen_genotypes.insert(key.clone())
                {
                    return Err(rejected(TruthValidationError::Contract(
                        TruthContractError::RecordsOutOfOrder,
                    )));
                }
                // El flujo debe aplicar el mismo contrato estructural que el
                // documento completo antes de comparar semántica: incluidos
                // portadores duplicados, índices y dosificaciones.
                validate_lexical_genotype(&genotype).map_err(TruthValidationError::Contract)?;
                // Con un documento sin copia con repetición, este acumulador solo
                // cambia el error (`UnsupportedVersion` en vez de lo que diga el
                // juez del caso), no la aceptación.
                carries_repeats |= matches!(genotype.case_kind, GenotypeCaseKind::Repeat);
                crate::truth_schema::validate_genotype(&genotype, header, &copy_identities)
                    .map_err(TruthValidationError::Contract)?;
                validate_genotype_semantics(&genotype, reference, &summaries)?;
                if let Some(previous) = individual_cases.first()
                    && previous.individual_id != genotype.individual_id
                {
                    let individual_id = previous.individual_id.clone();
                    judge_individual(&individual_cases, &individual_id, &summaries)?;
                    judged_individuals.insert(individual_id);
                    individual_cases.clear();
                }
                individual_cases.push(genotype);
                previous_genotype = Some(key);
                // **Inalcanzable** (DT-001, tanda 5): igual que el de copias, exige 2^64 registros.
                genotype_count = genotype_count.checked_add(1).ok_or_else(|| {
                    rejected(TruthValidationError::Contract(
                        TruthContractError::GenotypeCountMismatch {
                            declared: 0,
                            actual: u64::MAX,
                        },
                    ))
                })?;
            }
        }
    }
    let header = header.ok_or_else(|| rejected(TruthValidationError::EmptyDocument))?;
    // La versión desconocida ya se rechazó en la cabecera. El desajuste entre
    // la 4 y la 5 depende de todos los registros y solo puede decidirse aquí,
    // así que con un segundo defecto los errores de los dos lectores pueden
    // diferir; la aceptación coincide.
    validate_schema_version(&header, carries_repeats).map_err(TruthValidationError::Contract)?;
    if let Some(last) = individual_cases.first() {
        let individual_id = last.individual_id.clone();
        judge_individual(&individual_cases, &individual_id, &summaries)?;
        judged_individuals.insert(individual_id);
    }
    let without_cases = summaries
        .keys()
        .map(|identity| identity.individual_id.as_str())
        .filter(|individual_id| !judged_individuals.contains(*individual_id))
        .collect::<BTreeSet<_>>();
    for individual_id in without_cases {
        judge_individual(&[], individual_id, &summaries)?;
    }
    if header.copy_count != copy_count {
        return Err(rejected(TruthValidationError::Contract(
            TruthContractError::CopyCountMismatch {
                declared: header.copy_count,
                actual: copy_count,
            },
        )));
    }
    if header.genotype_count != genotype_count {
        return Err(rejected(TruthValidationError::Contract(
            TruthContractError::GenotypeCountMismatch {
                declared: header.genotype_count,
                actual: genotype_count,
            },
        )));
    }
    Ok(TruthStreamSummary {
        copy_count,
        genotype_count,
    })
}

/// Variante para quien ya haya leído el JSONL por flujo.
pub fn validate_and_reconstruct_truth(
    document: TruthDocument,
    reference: &Reference,
) -> Result<ReconstructedTruth, TruthValidationError> {
    if document.header.reference_id != reference.id().as_str() {
        return Err(rejected(TruthValidationError::ReferenceIdMismatch {
            declared: document.header.reference_id.clone(),
            actual: reference.id().as_str().to_owned(),
        }));
    }
    validate_reference_digest(&document.header, reference)?;
    let mut reconstructed = BTreeMap::new();
    let mut details = BTreeMap::new();
    for record in &document.records {
        if let TruthRecord::HaplotypeCopy(copy) = record {
            let result = reconstruct_copy(copy, reference)?;
            let ReconstructedCopy {
                bases,
                source_slices,
                pieces,
                event_loci,
                inherited_loss,
            } = result;
            reconstructed.insert(copy.identity.clone(), bases);
            details.insert(
                copy.identity.clone(),
                CopySummary {
                    source_slices,
                    pieces,
                    repeat_history: repeat_history_by_interval(&copy.history),
                    history: Cow::Borrowed(&copy.history),
                    event_loci,
                    inherited_loss,
                    origin: Cow::Borrowed(&copy.provenance.origin),
                },
            );
        }
    }
    for record in &document.records {
        if let TruthRecord::Genotype(genotype) = record {
            validate_genotype_semantics(genotype, reference, &details)?;
        }
    }
    validate_case_set(document_genotypes(&document), &details, None)?;
    Ok(ReconstructedTruth {
        document,
        copies: reconstructed,
    })
}

#[derive(Clone, Debug)]
struct ReconstructedCopy {
    bases: Vec<u8>,
    /// Cortes de referencia y SNPs dispersos. Nunca se conserva un mapa por
    /// base para el contraste posterior de genotipos en modo streaming.
    source_slices: Vec<SourceSlice>,
    pieces: Vec<PieceLayout>,
    event_loci: EventLoci,
    /// Primer tramo de la referencia que falta en el estado **heredado** de
    /// la copia —el de sus segmentos, antes del historial—: una deleción de
    /// segmento (un tramo que ningún segmento cubre ya lo rechazó
    /// `require_collinear_segments`). Ningún caso puede explicarlo, porque un
    /// caso *indel* exige su evento en el historial.
    inherited_loss: Option<(u64, u64)>,
}

/// Locus en coordenadas de la referencia de cada inserción y deleción del
/// historial, por identificador de evento, tal como lo deja la réplica.
type EventLoci = BTreeMap<String, (u64, u64)>;

#[derive(Clone, Debug)]
struct CopySummary<'truth> {
    source_slices: Vec<SourceSlice>,
    pieces: Vec<PieceLayout>,
    history: Cow<'truth, [TruthHistoryEvent]>,
    /// Posiciones en `history` de los alelos de repetición, por su intervalo y
    /// en el orden del historial (ver [`repeat_history_by_interval`]).
    repeat_history: BTreeMap<(u64, u64), Vec<usize>>,
    event_loci: EventLoci,
    inherited_loss: Option<(u64, u64)>,
    /// Procedencia declarada de la copia (`provenance.origin`).
    origin: Cow<'truth, str>,
}

impl<'truth> CopySummary<'truth> {
    fn from_reconstructed(
        value: ReconstructedCopy,
        history: Cow<'truth, [TruthHistoryEvent]>,
        origin: Cow<'truth, str>,
    ) -> Self {
        Self {
            source_slices: value.source_slices,
            pieces: value.pieces,
            repeat_history: repeat_history_by_interval(&history),
            history,
            event_loci: value.event_loci,
            inherited_loss: value.inherited_loss,
            origin,
        }
    }
}

/// Dónde están en el historial los alelos de repetición de una copia, por su
/// intervalo y en el orden del historial.
///
/// El juez de cada caso de repetición buscaba el alelo de la copia recorriendo
/// el historial entero, que tiene una entrada por locus aplicado: casos por
/// loci (paso F6.4.a.8j del plan de detección de microsatélites). Con la tabla
/// mira solo los alelos de ese intervalo, en el mismo orden, de modo que sigue
/// encontrando **el primero** que cumple.
fn repeat_history_by_interval(history: &[TruthHistoryEvent]) -> BTreeMap<(u64, u64), Vec<usize>> {
    let mut table = BTreeMap::<(u64, u64), Vec<usize>>::new();
    for (position, event) in history.iter().enumerate() {
        if let TruthHistoryEvent::RepeatAllele { start, end, .. } = event {
            table.entry((*start, *end)).or_default().push(position);
        }
    }
    table
}

#[derive(Clone, Debug)]
struct SourceSlice {
    start: u64,
    end: u64,
    /// Ordenado por posición; los cambios son escasos y se consultan por
    /// búsqueda binaria sobre el intervalo que sigue vigente.
    snps: Vec<TruthCurrentSnp>,
    /// Cierto cuando el tramo lo ocupa un **tracto realizado** y no las bases
    /// de la referencia.
    ///
    /// Sigue registrándose como tramo ocupado para que dos fragmentos no
    /// puedan reclamar la misma coordenada, pero un caso SNP que cayera dentro
    /// **no** puede leer de él la base: la copia no tiene ahí una base de
    /// referencia que consultar, igual que no la tiene una copia que perdió la
    /// coordenada en una deleción.
    ///
    /// # Es defensiva, y ninguna prueba la ejercita. Dicho sin adornos
    ///
    /// Desde que la decisión S-32 rechaza **antes de simular** cualquier
    /// cambio pequeño dentro del intervalo de un locus, este programa no puede
    /// producir un documento en el que un caso SNP caiga dentro de un tracto
    /// sustituido. La marca no sobra por eso: **este módulo no produce la
    /// verdad, la verifica**, y su trabajo es no fiarse de quien la escribió,
    /// venga de esta versión, de otra o de otro productor. Sin ella, un
    /// documento así leería en ese tramo la base de la referencia —que no es la
    /// que la copia tiene— y el contraste de genotipos daría por bueno un dato
    /// falso.
    ///
    /// Queda anotada como una de las comprobaciones de verificador que ninguna
    /// prueba dispara, que es la primera deuda del registro de deudas técnicas.
    replaced_by_repeat: bool,
}

#[derive(Clone, Debug, Eq, PartialEq)]
enum ReplayPiece {
    Source {
        source: SegmentSource,
        snps: Vec<TruthCurrentSnp>,
    },
    Inserted {
        bases: String,
        provenance: ModificationProvenance,
    },
    RepeatAllele {
        locus_id: String,
        source: SegmentSource,
        bases: String,
        blocks: Vec<TruthRepeatBlock>,
        provenance: ModificationProvenance,
    },
}

impl ReplayPiece {
    fn len(&self) -> Option<u64> {
        match self {
            Self::Source { source, .. } => source.end.checked_sub(source.start),
            Self::Inserted { bases, .. } | Self::RepeatAllele { bases, .. } => {
                u64::try_from(bases.len()).ok()
            }
        }
    }

    fn current(&self) -> TruthCurrentPiece {
        match self {
            Self::Source { source, snps } => TruthCurrentPiece::Source {
                source: source.clone(),
                snps: snps.clone(),
            },
            Self::Inserted { bases, provenance } => TruthCurrentPiece::Inserted {
                bases: bases.clone(),
                provenance: provenance.clone(),
            },
            Self::RepeatAllele {
                locus_id,
                source,
                bases,
                blocks,
                provenance,
            } => TruthCurrentPiece::RepeatAllele {
                locus_id: locus_id.clone(),
                source: source.clone(),
                bases: bases.clone(),
                blocks: blocks.clone(),
                provenance: provenance.clone(),
            },
        }
    }
}

/// Aplica un evento del historial a los fragmentos replicados.
///
/// Para una inserción o una deleción devuelve además su locus en coordenadas
/// de la referencia, calculado aquí con una regla propia y no leído del caso:
/// la inserción se ancla donde termina la última base de referencia que la
/// precede; la deleción cubre de la primera a la última base de referencia
/// que elimina, o es el ancla vacía de su inicio si solo elimina bases
/// insertadas.
///
/// Los dos rechazos por no haber ninguna base de referencia delante son
/// defensivos: el productor nunca inserta ni elimina en el punto cero, y la
/// réplica ya exige `point > 0` y `start > 0`, así que delante siempre queda al
/// menos el primer fragmento de la copia. Solo los alcanza un documento cuya
/// copia empieza por bases insertadas (un segmento con una inserción en su
/// primer punto), que ejercitan
/// `an_insertion_with_no_reference_base_before_it_does_not_replay` y
/// `a_deletion_of_inserted_bases_with_no_reference_base_before_it_does_not_replay`.
fn replay_event(
    copy: &HaplotypeCopyRecord,
    reference: &Reference,
    pieces: &mut Vec<ReplayPiece>,
    event: &TruthHistoryEvent,
) -> Result<Option<(u64, u64)>, TruthValidationError> {
    match event {
        TruthHistoryEvent::Snp {
            event_id,
            position,
            provenance,
        } => {
            // **Equivalentes** (DT-001, tanda 5) las igualdades de posición y de
            // procedencia de este filtro: un cambio con el mismo identificador y otra
            // posición o procedencia ya lo ha rechazado `validate_structural_event_ids`;
            // y la de identificador, porque una modificación con otro identificador pero
            // la misma posición y procedencia ya está como cambio heredado en el fragmento
            // y el evento chocaría con él en la comprobación de «un cambio ya en esa
            // posición», con el mismo error.
            let mut modifications = copy
                .segments
                .iter()
                .flat_map(|segment| segment.modifications.iter())
                .filter(|modification| {
                    matches!(modification, TruthModification::Snp { event_id: id, position: candidate, provenance: origin, .. } if id == event_id && candidate == position && origin == provenance)
                });
            let modification = modifications.next();
            let Some(TruthModification::Snp {
                reference: expected,
                alternate,
                ..
            }) = modification
            else {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            };
            // **Inalcanzable** (DT-001, tanda 5): dos modificaciones de segmento para el mismo evento exigen
            // el mismo identificador repetido, y `validate_structural_event_ids` (que
            // `reconstruct_copy` ejecuta antes que ninguna réplica) rechaza todo
            // identificador heredado repetido. Intento ejecutado:
            // `unreachable_guards_are_shadowed_by_other_checks` (da el mismo error, pero de
            // aquella comprobación).
            if modifications.next().is_some() {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            let mut matching = pieces
                .iter_mut()
                .filter_map(|piece| match piece {
                    ReplayPiece::Source { source, snps }
                        if source.start <= *position && *position < source.end =>
                    {
                        Some((source, snps))
                    }
                    _ => None,
                })
                .collect::<Vec<_>>();
            let Some((source, snps)) = matching.pop() else {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            };
            // **Redundante para la seguridad desde DT-048**: con segmentos colineales los
            // fragmentos de fuente no se solapan y ninguna posición cae en dos, así que
            // solo llega aquí un documento que además incumple la colinealidad, y
            // `require_collinear_segments` lo rechazaría igual después. Quitarla no
            // dejaría pasar nada; se mantiene por el diagnóstico, más preciso
            // (`a_base_change_covered_by_two_fragments_does_not_replay`).
            if !matching.is_empty() {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            let sequence = validated_source(copy, source, reference)?;
            // **Inalcanzable** (DT-001, tanda 5) en el objetivo x86_64 (`usize` y `u64` tienen el mismo
            // ancho): la conversión no falla. En un objetivo de 32 bits sería alcanzable.
            let index = usize::try_from(*position).map_err(|_| {
                rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })
            })?;
            // **Inalcanzable** (DT-001, tanda 5): ver `bases_of`; la huella de la referencia falla antes.
            let bases = sequence
                .sequence()
                .as_bytes()
                .map_err(|_| rejected(TruthValidationError::BasesNotMaterialized))?;
            if bases.get(index) != Some(&expected.as_bytes()[0]) {
                return Err(rejected(TruthValidationError::SnpReferenceMismatch {
                    event_id: event_id.clone(),
                    position: *position,
                }));
            }
            // Guarda defensiva. Un documento con dos cambios en la misma posición lo
            // rechaza el contrato (`ModificationsNotOrdered` dentro de un segmento;
            // los segmentos solapados que lo evitarían con una deleción que salta el
            // cambio, `OverlappingModifications`), así que solo se llega aquí con un
            // documento en memoria que no pasó el contrato:
            // `two_base_changes_at_one_position_do_not_replay` lo ejercita y su
            // mutante muere. La rama «idéntico» (`Ok(None)`) es
            // inalcanzable: exigiría el mismo identificador dos veces, que
            // `validate_structural_event_ids` rechaza antes.
            if let Some(existing) = snps.iter().find(|snp| snp.position == *position) {
                if existing.event_id == *event_id
                    && existing.reference == *expected
                    && existing.alternate == *alternate
                    && existing.provenance == *provenance
                {
                    return Ok(None);
                }
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            // La verdad publica los cambios de base de cada tramo ordenados
            // por posición, no por orden de llegada al historial: se insertan
            // en su sitio, sin ordenar después ni relajar la comparación
            // final. La posición no está repetida (se rechazó arriba).
            let insert_at = snps.partition_point(|snp| snp.position < *position);
            snps.insert(
                insert_at,
                TruthCurrentSnp {
                    event_id: event_id.clone(),
                    position: *position,
                    reference: expected.clone(),
                    alternate: alternate.clone(),
                    provenance: provenance.clone(),
                },
            );
        }
        TruthHistoryEvent::Insertion {
            point,
            bases,
            provenance,
            ..
        } => {
            let length = replay_length(copy, pieces)?;
            // **Equivalente** (DT-001, tanda 5) la mitad `*point == 0`: sin ella, el punto cero deja
            // vacío lo de delante y rechaza el ancla de más abajo con el mismo error.
            if *point == 0 || *point >= length {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            split_replay_at(copy, pieces, *point)?;
            let index = piece_boundary_index(copy, pieces, *point)?;
            let Some(anchor) = replay_reference_anchor(&pieces[..index]) else {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            };
            pieces.insert(
                index,
                ReplayPiece::Inserted {
                    bases: bases.clone(),
                    provenance: provenance.clone(),
                },
            );
            return Ok(Some((anchor, anchor)));
        }
        TruthHistoryEvent::Deletion { start, end, .. } => {
            let length = replay_length(copy, pieces)?;
            // **Equivalente** (DT-001, tanda 5) la mitad `start >= end`: con ella ignorada no se quita
            // ningún fragmento (tras partir en los dos extremos no hay fragmento que
            // empiece antes del final y acabe después del inicio) y el `match` de
            // abajo rechaza con el mismo error. `*start == 0` y `*end >= length`
            // sí se distinguen (`a_deletion_may_not_reach_the_end_of_the_copy` y
            // `a_deletion_from_the_first_base_does_not_replay`).
            if *start == 0 || start >= end || *end >= length {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            split_replay_at(copy, pieces, *start)?;
            split_replay_at(copy, pieces, *end)?;
            let mut cursor = 0_u64;
            let mut retained = Vec::with_capacity(pieces.len());
            let mut anchor = None;
            let mut removed: Option<(u64, u64)> = None;
            // **Inalcanzable** (DT-001, tanda 5) los dos rechazos de `piece.len()` y de la suma: ningún
            // fragmento de fuente tiene el final por detrás del inicio (los iniciales se
            // forman con `push_initial_source`, que ignora los vacíos y rechaza los
            // invertidos, y los que parte la réplica lo hacen dentro de su intervalo) y
            // `replay_length`, en esta misma rama, ya sumó las longitudes de todos sin
            // desbordar.
            for piece in pieces.drain(..) {
                let next = cursor
                    .checked_add(piece.len().ok_or_else(|| {
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?)
                    .ok_or_else(|| {
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?;
                if next <= *start || cursor >= *end {
                    retained.push(piece);
                } else {
                    // `retained` contiene aquí justo lo que precede al
                    // primer fragmento eliminado.
                    if anchor.is_none() {
                        anchor = Some(replay_reference_anchor(&retained));
                    }
                    if let Some((first, last)) = replay_reference_interval(&piece) {
                        removed = Some(removed.map_or((first, last), |(low, high)| {
                            (low.min(first), high.max(last))
                        }));
                    }
                }
                cursor = next;
            }
            *pieces = retained;
            return match (removed, anchor.flatten()) {
                (Some(span), _) => Ok(Some(span)),
                (None, Some(anchor)) => Ok(Some((anchor, anchor))),
                (None, None) => Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })),
            };
        }
        // **El alelo de repetición se replica en coordenadas de referencia**,
        // no en las vigentes: sustituye el intervalo declarado del locus, que
        // tiene que llegar entero dentro de un solo fragmento de fuente y sin
        // ningún cambio de base dentro. Es la misma regla que impone el estado
        // haplotípico, escrita aquí por segunda vez y a propósito: si las dos
        // no coincidieran, la comparación de fragmentos vigentes lo diría.
        TruthHistoryEvent::RepeatAllele {
            locus_id,
            start,
            end,
            bases,
            blocks,
            provenance,
            ..
        } => {
            if start >= end || bases.is_empty() {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            }
            let index = pieces
                .iter()
                .position(|piece| match piece {
                    ReplayPiece::Source { source, snps } => {
                        source.start <= *start
                            && *end <= source.end
                            && !snps
                                .iter()
                                .any(|snp| (*start..*end).contains(&snp.position))
                    }
                    ReplayPiece::Inserted { .. } | ReplayPiece::RepeatAllele { .. } => false,
                })
                .ok_or_else(|| {
                    rejected(TruthValidationError::HistoryDoesNotReplay {
                        identity: copy.identity.clone(),
                    })
                })?;
            // **Inalcanzable** (DT-001, tanda 5): `index` salió de buscar un fragmento de fuente, de modo
            // que este fragmento lo es siempre.
            let ReplayPiece::Source { source, snps } = pieces[index].clone() else {
                return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                }));
            };
            let mut replacement = Vec::with_capacity(3);
            if source.start < *start {
                let mut left = source.clone();
                left.end = *start;
                replacement.push(ReplayPiece::Source {
                    source: left,
                    snps: snps
                        .iter()
                        .filter(|snp| snp.position < *start)
                        .cloned()
                        .collect(),
                });
            }
            let mut replaced = source.clone();
            replaced.start = *start;
            replaced.end = *end;
            replacement.push(ReplayPiece::RepeatAllele {
                locus_id: locus_id.clone(),
                source: replaced,
                bases: bases.clone(),
                blocks: blocks.clone(),
                provenance: provenance.clone(),
            });
            if *end < source.end {
                let mut right = source.clone();
                right.start = *end;
                replacement.push(ReplayPiece::Source {
                    source: right,
                    snps: snps
                        .iter()
                        .filter(|snp| snp.position >= *end)
                        .cloned()
                        .collect(),
                });
            }
            pieces.splice(index..=index, replacement);
        }
    }
    Ok(None)
}

/// Fin de la última base de referencia entre `before`, o nada si no hay
/// ninguna.
fn replay_reference_anchor(before: &[ReplayPiece]) -> Option<u64> {
    before
        .iter()
        .rev()
        .find_map(|piece| replay_reference_interval(piece).map(|(_, end)| end))
}

fn replay_reference_interval(piece: &ReplayPiece) -> Option<(u64, u64)> {
    match piece {
        ReplayPiece::Source { source, .. } | ReplayPiece::RepeatAllele { source, .. } => {
            Some((source.start, source.end))
        }
        ReplayPiece::Inserted { .. } => None,
    }
}

/// Longitud actual de la copia. De los dos rechazos, el de la suma que desborda
/// es alcanzable (dos segmentos de casi 2^64 bases; lo ejercita
/// `lengths_that_overflow_do_not_replay`); el de `piece.len()` (**Inalcanzable** (DT-001, tanda 5))
/// no, porque ningún fragmento de fuente tiene el final por detrás del inicio.
fn replay_length(
    copy: &HaplotypeCopyRecord,
    pieces: &[ReplayPiece],
) -> Result<u64, TruthValidationError> {
    pieces.iter().try_fold(0_u64, |length, piece| {
        length
            .checked_add(piece.len().ok_or_else(|| {
                rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })
            })?)
            .ok_or_else(|| {
                rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })
            })
    })
}

/// Índice del fragmento que empieza en `point`.
///
/// **Inalcanzable** (DT-001, tanda 5) los tres rechazos: `piece.len()` y la suma ya los validó
/// `replay_length` en la misma rama sobre los mismos fragmentos, y el rechazo
/// final exigiría que no hubiera frontera en `point`, cuando
/// `split_replay_at(point)` acaba de crearla (el punto es menor que la
/// longitud: lo exigen las dos ramas que llaman, antes de partir).
fn piece_boundary_index(
    copy: &HaplotypeCopyRecord,
    pieces: &[ReplayPiece],
    point: u64,
) -> Result<usize, TruthValidationError> {
    let mut cursor = 0_u64;
    for (index, piece) in pieces.iter().enumerate() {
        if cursor == point {
            return Ok(index);
        }
        cursor = cursor
            .checked_add(piece.len().ok_or_else(|| {
                rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })
            })?)
            .ok_or_else(|| {
                rejected(TruthValidationError::HistoryDoesNotReplay {
                    identity: copy.identity.clone(),
                })
            })?;
    }
    if cursor == point {
        Ok(pieces.len())
    } else {
        // **Inalcanzable** (DT-001, tanda 5): la frontera existe siempre; ver la nota de la función.
        Err(rejected(TruthValidationError::HistoryDoesNotReplay {
            identity: copy.identity.clone(),
        }))
    }
}

/// Parte la réplica en `point`.
///
/// **Inalcanzable** (DT-001, tanda 5) `piece.len()` y la suma de las longitudes (ya validadas por
/// `replay_length` en la rama que llama), la suma de `source.start + split`
/// (menor que el final del fragmento, que cabe), la conversión a `usize`
/// (x86_64) y el rechazo de después del bucle (los dos llamantes rechazan
/// antes un punto que no sea menor que la longitud, y todo punto menor cae en
/// algún fragmento). Las mitades `point == cursor` y `point == end` de la
/// primera comprobación son **equivalentes** entre sí: la frontera se ve desde
/// los dos lados.
fn split_replay_at(
    copy: &HaplotypeCopyRecord,
    pieces: &mut Vec<ReplayPiece>,
    point: u64,
) -> Result<(), TruthValidationError> {
    let mut cursor = 0_u64;
    for index in 0..pieces.len() {
        let length = pieces[index].len().ok_or_else(|| {
            rejected(TruthValidationError::HistoryDoesNotReplay {
                identity: copy.identity.clone(),
            })
        })?;
        let end = cursor.checked_add(length).ok_or_else(|| {
            rejected(TruthValidationError::HistoryDoesNotReplay {
                identity: copy.identity.clone(),
            })
        })?;
        if point == cursor || point == end {
            return Ok(());
        }
        if cursor < point && point < end {
            let right = match &pieces[index] {
                ReplayPiece::Source { source, snps } => {
                    let split = point - cursor;
                    let source_point = source.start.checked_add(split).ok_or_else(|| {
                        // **Inalcanzable** (DT-001, tanda 5): `source.start + split` no desborda; ver la nota de la función.
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?;
                    let mut left_source = source.clone();
                    left_source.end = source_point;
                    let mut right_source = source.clone();
                    right_source.start = source_point;
                    let left_snps = snps
                        .iter()
                        .filter(|snp| snp.position < source_point)
                        .cloned()
                        .collect();
                    let right_snps = snps
                        .iter()
                        .filter(|snp| snp.position >= source_point)
                        .cloned()
                        .collect();
                    pieces[index] = ReplayPiece::Source {
                        source: left_source,
                        snps: left_snps,
                    };
                    ReplayPiece::Source {
                        source: right_source,
                        snps: right_snps,
                    }
                }
                ReplayPiece::RepeatAllele { .. } => {
                    return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                        identity: copy.identity.clone(),
                    }));
                }
                ReplayPiece::Inserted { bases, provenance } => {
                    let split = usize::try_from(point - cursor).map_err(|_| {
                        // **Inalcanzable** (DT-001, tanda 5): `usize` y `u64` del mismo ancho; ver la nota de la función.
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?;
                    let (left, right) = bases.split_at(split);
                    let right = ReplayPiece::Inserted {
                        bases: right.to_owned(),
                        provenance: provenance.clone(),
                    };
                    pieces[index] = ReplayPiece::Inserted {
                        bases: left.to_owned(),
                        provenance: provenance.clone(),
                    };
                    right
                }
            };
            pieces.insert(index + 1, right);
            return Ok(());
        }
        cursor = end;
    }
    // **Inalcanzable** (DT-001, tanda 5): todo punto menor que la longitud cae en un fragmento; ver la nota de la función.
    Err(rejected(TruthValidationError::HistoryDoesNotReplay {
        identity: copy.identity.clone(),
    }))
}

#[derive(Clone, Debug)]
struct PieceLayout {
    start: u64,
    end: u64,
    inserted_provenance: Option<ModificationProvenance>,
    /// El fragmento es un tracto de repetición, no bases insertadas, aunque
    /// también lleve procedencia.
    tract: bool,
}

fn reconstruct_copy(
    copy: &HaplotypeCopyRecord,
    reference: &Reference,
) -> Result<ReconstructedCopy, TruthValidationError> {
    let mut pieces = initial_replay_from_segments(copy)?;
    let mut event_loci = EventLoci::new();
    for event in &copy.history {
        if let Some(locus) = replay_event(copy, reference, &mut pieces, event)? {
            // Con identificadores repetidos, el contraste de genotipos lee el
            // primer evento del historial; el locus es el de ese mismo.
            event_loci
                .entry(history_event_id(event).to_owned())
                .or_insert(locus);
        }
    }
    let replayed = pieces.iter().map(ReplayPiece::current).collect::<Vec<_>>();
    if replayed != copy.current_pieces {
        return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
            identity: copy.identity.clone(),
        }));
    }
    let mut result = materialize_replayed_copy(copy, reference, &pieces)?;
    require_collinear_segments(copy, reference)?;
    result.event_loci = event_loci;
    result.inherited_loss = inherited_loss(copy);
    Ok(result)
}

/// El estado heredado de la copia —sus segmentos, antes del historial— recorre
/// su secuencia de la referencia **una sola vez y en orden creciente** (DT-048):
/// cada segmento es de la secuencia de la copia y no está vacío ni invertido, el
/// primero empieza en 0, cada uno empieza donde acaba el anterior y el último
/// acaba en la longitud de la secuencia. Un hueco, un solape, una duplicación,
/// una rotación o un desorden rompen el recorrido en el primer segmento que no
/// empieza donde iba; un segmento vacío, invertido o de otra secuencia (DT-050),
/// en ese mismo segmento aunque empiece donde iba.
///
/// La regla no depende de ninguna otra guarda: el contrato ya rechaza por bytes
/// un segmento vacío, invertido o de otra secuencia (`InvalidSegmentInterval`,
/// `SegmentSourceMismatch`), pero un documento en memoria no pasa por él, y un
/// segmento invertido acompañado de deleciones de segmento que reponen el cursor
/// supera la réplica.
///
/// Las modificaciones de un segmento no cuentan aquí: un cambio de base y una
/// inserción heredados no consumen referencia, y una deleción de segmento la
/// rechaza el juez del conjunto de casos (ver `inherited_loss`).
///
/// Se exige **después** de replicar y materializar la copia, y no al formar la
/// réplica inicial, a propósito: un documento que además incumple otra regla de
/// la réplica —dos fragmentos que se solapan, longitudes que desbordan, un
/// fragmento más allá del final, un cambio de base en la última posición
/// representable— lo sigue parando la guarda de esa regla, con su motivo, y
/// esas guardas siguen siendo alcanzables. Esta para lo que ninguna otra ve: el
/// hueco, la rotación, el desorden y el tramo duplicado que una deleción del
/// historial vuelve a quitar.
fn require_collinear_segments(
    copy: &HaplotypeCopyRecord,
    reference: &Reference,
) -> Result<(), TruthValidationError> {
    // Si la copia conserva algún fragmento de fuente o tracto, la materialización
    // ya ha rechazado una secuencia que la referencia no tiene (el contrato ata la
    // secuencia de cada fragmento a la de la copia). Aquí solo llega sin ella un
    // documento en memoria que no pasó el contrato, como una copia sin segmentos.
    let length = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == copy.identity.sequence_id)
        .map(ReferenceSequence::length)
        .ok_or_else(|| {
            rejected(TruthValidationError::MissingReferenceSequence {
                sequence_id: copy.identity.sequence_id.clone(),
            })
        })?;
    let mut reached = 0;
    for segment in &copy.segments {
        let SegmentSource {
            sequence_id,
            start,
            end,
            ..
        } = &segment.source;
        if *start != reached || *end <= *start || *sequence_id != copy.identity.sequence_id {
            return Err(rejected(TruthValidationError::SegmentsNotCollinear {
                identity: copy.identity.clone(),
                reached,
                next_start: *start,
            }));
        }
        // **Equivalente** (DT-048) `reached = reached.max(*end)`: la condición de
        // arriba acaba de exigir `start == reached` y `end > start`, así que `end`
        // es siempre mayor.
        reached = *end;
    }
    if reached != length {
        return Err(rejected(TruthValidationError::SegmentsNotCollinear {
            identity: copy.identity.clone(),
            reached,
            next_start: length,
        }));
    }
    Ok(())
}

/// Primera deleción de segmento de la copia: un tramo de la referencia que el
/// estado heredado no conserva.
///
/// Con los segmentos recorriendo la referencia una vez y en orden (lo exige
/// antes `require_collinear_segments`) y sin deleciones de segmento, cada base
/// de referencia que falta en la copia final la quitó una deleción del
/// historial, y esa deleción tiene su caso. Por eso basta con mirar los
/// segmentos, y no los loci de las deleciones: el locus de una deleción va de
/// la primera a la última base que elimina y puede abarcar un hueco anterior,
/// de modo que una pérdida heredada escondida dentro de él pasaría por
/// explicada.
fn inherited_loss(copy: &HaplotypeCopyRecord) -> Option<(u64, u64)> {
    copy.segments
        .iter()
        .flat_map(|segment| segment.modifications.iter())
        .find_map(|modification| match modification {
            TruthModification::Deletion { start, end, .. } => Some((*start, *end)),
            TruthModification::Snp { .. } | TruthModification::Insertion { .. } => None,
        })
}

fn history_event_id(event: &TruthHistoryEvent) -> &str {
    match event {
        TruthHistoryEvent::Snp { event_id, .. }
        | TruthHistoryEvent::Insertion { event_id, .. }
        | TruthHistoryEvent::Deletion { event_id, .. }
        | TruthHistoryEvent::RepeatAllele { event_id, .. } => event_id,
    }
}

/// Materializa el estado heredado de los segmentos antes de los eventos en
/// coordenadas vigentes. Inserciones y deleciones del segmento son, por tanto,
/// parte del punto de partida y no exigen un duplicado en `history`.
fn initial_replay_from_segments(
    copy: &HaplotypeCopyRecord,
) -> Result<Vec<ReplayPiece>, TruthValidationError> {
    validate_structural_event_ids(copy)?;
    let mut pieces = Vec::new();
    for segment in &copy.segments {
        let mut cursor = segment.source.start;
        for modification in &segment.modifications {
            match modification {
                TruthModification::Snp {
                    event_id,
                    position,
                    reference,
                    alternate,
                    provenance,
                } => {
                    // `with_snp` conserva el descriptor de segmento y añade
                    // el mismo SNP al historial activo. En ese caso el
                    // historial es quien lo aplica; un SNP heredado sin tal
                    // evento sí forma parte del estado inicial.
                    // **Equivalente** (DT-001, tanda 5) `historical_position == position`:
                    // un evento del historial con el mismo identificador y otra posición lo
                    // rechaza antes `validate_structural_event_ids`.
                    if copy.history.iter().any(|event| {
                        matches!(event, TruthHistoryEvent::Snp { event_id: historical_id, position: historical_position, provenance: historical_provenance }
                            if historical_id == event_id
                                && historical_position == position
                                && historical_provenance == provenance)
                    }) {
                        continue;
                    }
                    push_initial_source(
                        copy,
                        &mut pieces,
                        &segment.source,
                        cursor,
                        *position,
                        Vec::new(),
                    )?;
                    let end = position.checked_add(1).ok_or_else(|| {
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?;
                    push_initial_source(
                        copy,
                        &mut pieces,
                        &segment.source,
                        *position,
                        end,
                        vec![TruthCurrentSnp {
                            event_id: event_id.clone(),
                            position: *position,
                            reference: reference.clone(),
                            alternate: alternate.clone(),
                            provenance: provenance.clone(),
                        }],
                    )?;
                    cursor = end;
                }
                TruthModification::Insertion {
                    point,
                    bases,
                    provenance,
                    ..
                } => {
                    push_initial_source(
                        copy,
                        &mut pieces,
                        &segment.source,
                        cursor,
                        *point,
                        Vec::new(),
                    )?;
                    pieces.push(ReplayPiece::Inserted {
                        bases: bases.clone(),
                        provenance: provenance.clone(),
                    });
                    cursor = *point;
                }
                TruthModification::Deletion { start, end, .. } => {
                    push_initial_source(
                        copy,
                        &mut pieces,
                        &segment.source,
                        cursor,
                        *start,
                        Vec::new(),
                    )?;
                    cursor = *end;
                }
            }
        }
        push_initial_source(
            copy,
            &mut pieces,
            &segment.source,
            cursor,
            segment.source.end,
            Vec::new(),
        )?;
    }
    Ok(pieces)
}

/// Los identificadores heredados son globales por copia y conservan su clase.
/// Repetirlos entre segmentos o reutilizarlos desde `history` haría ambiguo el
/// origen, aun cuando las piezas finales parezcan reproducibles. La única
/// excepción es el mismo SNP de segmento e historial que `with_snp` conserva
/// deliberadamente; nunca habilita una colisión entre clases.
fn validate_structural_event_ids(copy: &HaplotypeCopyRecord) -> Result<(), TruthValidationError> {
    #[derive(Clone, Copy, Eq, PartialEq)]
    enum SegmentEventClass {
        Snp,
        Structural,
    }

    let mut inherited = BTreeMap::new();
    for modification in copy
        .segments
        .iter()
        .flat_map(|segment| segment.modifications.iter())
    {
        let (event_id, class) = match modification {
            TruthModification::Snp { event_id, .. } => (event_id, SegmentEventClass::Snp),
            TruthModification::Insertion { event_id, .. }
            | TruthModification::Deletion { event_id, .. } => {
                (event_id, SegmentEventClass::Structural)
            }
        };
        if inherited.insert(event_id, class).is_some() {
            return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                identity: copy.identity.clone(),
            }));
        }
    }
    for event in &copy.history {
        let event_id = match event {
            TruthHistoryEvent::Snp { event_id, .. }
            | TruthHistoryEvent::Insertion { event_id, .. }
            | TruthHistoryEvent::Deletion { event_id, .. }
            | TruthHistoryEvent::RepeatAllele { event_id, .. } => event_id,
        };
        let Some(class) = inherited.get(event_id) else {
            continue;
        };
        // **Equivalentes** (DT-001, tanda 5): de las tres condiciones solo decide la
        // tercera. `segment_snp_matches_history` exige una modificación de
        // **cambio de base** con ese identificador, posición y procedencia, de modo
        // que ya es falsa si la clase heredada es estructural o si el evento no es
        // un cambio; y las tres igualdades que contiene son las mismas que usa el
        // filtro de la réplica de un cambio de base, que rechazaría igual con el
        // mismo error. Quitar cualquier condición sola sobrevive; no hay entrada
        // que las distinga (`a_history_event_sharing_an_inherited_identifier_must_be_the_same_base_change`
        // prueba las cuatro formas y todas dan el mismo rechazo).
        let permitted_snp_duplicate = *class == SegmentEventClass::Snp
            && matches!(event, TruthHistoryEvent::Snp { .. })
            && segment_snp_matches_history(copy, event);
        if !permitted_snp_duplicate {
            return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
                identity: copy.identity.clone(),
            }));
        }
    }
    Ok(())
}

fn segment_snp_matches_history(copy: &HaplotypeCopyRecord, event: &TruthHistoryEvent) -> bool {
    let TruthHistoryEvent::Snp {
        event_id,
        position,
        provenance,
    } = event
    else {
        return false;
    };
    copy.segments
        .iter()
        .flat_map(|segment| segment.modifications.iter())
        .filter(|modification| {
            matches!(modification, TruthModification::Snp { event_id: candidate_id, position: candidate_position, provenance: candidate_provenance, .. }
                if candidate_id == event_id
                    && candidate_position == position
                    && candidate_provenance == provenance)
        })
        .take(2)
        .count()
        == 1
}

fn push_initial_source(
    copy: &HaplotypeCopyRecord,
    pieces: &mut Vec<ReplayPiece>,
    source: &SegmentSource,
    start: u64,
    end: u64,
    snps: Vec<TruthCurrentSnp>,
) -> Result<(), TruthValidationError> {
    if start > end {
        return Err(rejected(TruthValidationError::HistoryDoesNotReplay {
            identity: copy.identity.clone(),
        }));
    }
    if start < end {
        let mut source = source.clone();
        source.start = start;
        source.end = end;
        pieces.push(ReplayPiece::Source { source, snps });
    }
    Ok(())
}

fn materialize_replayed_copy(
    copy: &HaplotypeCopyRecord,
    reference: &Reference,
    pieces: &[ReplayPiece],
) -> Result<ReconstructedCopy, TruthValidationError> {
    // **Inalcanzable** (DT-001, tanda 5) las conversiones `usize` a `u64` y de `u64` a `usize`
    // de esta función: x86_64 tiene los dos del mismo ancho. El rechazo de la suma
    // `snp.position + 1` tampoco se alcanza: la posición de todo cambio es menor que
    // el final de su fragmento (la de un cambio heredado en `u64::MAX` ya la rechazó
    // `initial_replay_from_segments`) y el final cabe. Y `BasesNotMaterialized`,
    // como en `bases_of`.
    let mut bases = Vec::new();
    let mut source_slices = Vec::<SourceSlice>::new();
    let mut layouts = Vec::with_capacity(pieces.len());
    for piece in pieces {
        let start = u64::try_from(bases.len()).map_err(|_| {
            rejected(TruthValidationError::SourceIntervalOutOfBounds {
                identity: copy.identity.clone(),
            })
        })?;
        let inserted_provenance = match piece {
            ReplayPiece::Source { source, snps } => {
                let sequence = validated_source(copy, source, reference)?;
                for snp in snps {
                    let bases = sequence
                        .sequence()
                        .as_bytes()
                        .map_err(|_| rejected(TruthValidationError::BasesNotMaterialized))?;
                    let found = bases[usize::try_from(snp.position).map_err(|_| {
                        rejected(TruthValidationError::SourceIntervalOutOfBounds {
                            identity: copy.identity.clone(),
                        })
                    })?];
                    if found != snp.reference.as_bytes()[0] {
                        return Err(rejected(TruthValidationError::SnpReferenceMismatch {
                            event_id: snp.event_id.clone(),
                            position: snp.position,
                        }));
                    }
                }
                // **Redundante para la seguridad desde DT-048**: dos fragmentos que se
                // solapan exigen segmentos que no recorren la referencia una vez y en
                // orden, y `require_collinear_segments` los rechazaría igual después.
                // Quitarla no dejaría pasar nada; se mantiene por el diagnóstico, que
                // nombra la posición duplicada (`two_overlapping_source_fragments_are_rejected`).
                for prior in &source_slices {
                    if prior.start < source.end && source.start < prior.end {
                        return Err(rejected(
                            TruthValidationError::DuplicateCurrentSourcePosition {
                                identity: copy.identity.clone(),
                                position: prior.start.max(source.start),
                            },
                        ));
                    }
                }
                let mut cursor = source.start;
                for snp in snps {
                    let position = usize::try_from(snp.position).map_err(|_| {
                        // **Inalcanzable** (DT-001, tanda 5): conversión de ancho; ver `materialize_replayed_copy`.
                        rejected(TruthValidationError::SourceIntervalOutOfBounds {
                            identity: copy.identity.clone(),
                        })
                    })?;
                    let start_index = usize::try_from(cursor).map_err(|_| {
                        // **Inalcanzable** (DT-001, tanda 5): conversión de ancho; ver `materialize_replayed_copy`.
                        rejected(TruthValidationError::SourceIntervalOutOfBounds {
                            identity: copy.identity.clone(),
                        })
                    })?;
                    bases.extend_from_slice(&bases_of(sequence)?[start_index..position]);
                    bases.push(snp.alternate.as_bytes()[0]);
                    cursor = snp.position.checked_add(1).ok_or_else(|| {
                        // **Inalcanzable** (DT-001, tanda 5): `snp.position + 1` no desborda; ver `materialize_replayed_copy`.
                        rejected(TruthValidationError::HistoryDoesNotReplay {
                            identity: copy.identity.clone(),
                        })
                    })?;
                }
                let start_index = usize::try_from(cursor).map_err(|_| {
                    // **Inalcanzable** (DT-001, tanda 5): conversión de ancho; ver `materialize_replayed_copy`.
                    rejected(TruthValidationError::SourceIntervalOutOfBounds {
                        identity: copy.identity.clone(),
                    })
                })?;
                let end_index = usize::try_from(source.end).map_err(|_| {
                    // **Inalcanzable** (DT-001, tanda 5): conversión de ancho; ver `materialize_replayed_copy`.
                    rejected(TruthValidationError::SourceIntervalOutOfBounds {
                        identity: copy.identity.clone(),
                    })
                })?;
                bases.extend_from_slice(&bases_of(sequence)?[start_index..end_index]);
                source_slices.push(SourceSlice {
                    start: source.start,
                    end: source.end,
                    snps: snps.clone(),
                    replaced_by_repeat: false,
                });
                None
            }
            ReplayPiece::Inserted {
                bases: inserted,
                provenance,
            } => {
                bases.extend_from_slice(inserted.as_bytes());
                Some(provenance.clone())
            }
            ReplayPiece::RepeatAllele {
                source,
                bases: tract,
                provenance,
                ..
            } => {
                // El tracto sustituye al intervalo, de modo que ese tramo de la
                // referencia queda ocupado y ningún otro fragmento puede
                // reclamarlo.
                // **Redundante para la seguridad desde DT-048**: dos fragmentos que se
                // solapan exigen segmentos que no recorren la referencia una vez y en
                // orden, y `require_collinear_segments` los rechazaría igual después.
                // Quitarla no dejaría pasar nada; se mantiene por el diagnóstico, que
                // nombra la posición duplicada (`a_repeat_tract_over_an_occupied_stretch_is_rejected`).
                for prior in &source_slices {
                    if prior.start < source.end && source.start < prior.end {
                        return Err(rejected(
                            TruthValidationError::DuplicateCurrentSourcePosition {
                                identity: copy.identity.clone(),
                                position: prior.start.max(source.start),
                            },
                        ));
                    }
                }
                validated_source(copy, source, reference)?;
                bases.extend_from_slice(tract.as_bytes());
                source_slices.push(SourceSlice {
                    start: source.start,
                    end: source.end,
                    snps: Vec::new(),
                    replaced_by_repeat: true,
                });
                Some(provenance.clone())
            }
        };
        let end = u64::try_from(bases.len()).map_err(|_| {
            // **Inalcanzable** (DT-001, tanda 5): conversión de ancho; ver `materialize_replayed_copy`.
            rejected(TruthValidationError::SourceIntervalOutOfBounds {
                identity: copy.identity.clone(),
            })
        })?;
        layouts.push(PieceLayout {
            start,
            end,
            inserted_provenance,
            tract: matches!(piece, ReplayPiece::RepeatAllele { .. }),
        });
    }
    Ok(ReconstructedCopy {
        bases,
        source_slices,
        pieces: layouts,
        event_loci: EventLoci::new(),
        inherited_loss: None,
    })
}

fn validated_source<'a>(
    copy: &HaplotypeCopyRecord,
    source: &SegmentSource,
    reference: &'a Reference,
) -> Result<&'a ReferenceSequence, TruthValidationError> {
    if source.reference_id != reference.id().as_str()
        || source.source.kind != TruthSourceKind::Reference
        || source.source.source_id != reference.id().as_str()
    {
        return Err(rejected(TruthValidationError::SourceDoesNotNameReference {
            identity: copy.identity.clone(),
        }));
    }
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == source.sequence_id)
        .ok_or_else(|| {
            rejected(TruthValidationError::MissingReferenceSequence {
                sequence_id: source.sequence_id.clone(),
            })
        })?;
    if source.sequence_length != sequence.length() {
        return Err(rejected(TruthValidationError::SourceLengthMismatch {
            sequence_id: source.sequence_id.clone(),
            declared: source.sequence_length,
            actual: sequence.length(),
        }));
    }
    // **Equivalente** (DT-001, tanda 5) la mitad `source.start >= source.end`: todo fragmento que llega
    // aquí tiene el inicio por delante del final (los iniciales los forma
    // `push_initial_source` y los de un tracto los exige la guarda de
    // `replay_event`). `source.end > sequence.length()` solo se alcanza con un
    // documento en memoria (el contrato exige el final menor o igual que la
    // longitud declarada, y la comprobación de arriba que esta sea la real).
    if source.start >= source.end || source.end > sequence.length() {
        return Err(rejected(TruthValidationError::SourceIntervalOutOfBounds {
            identity: copy.identity.clone(),
        }));
    }
    Ok(sequence)
}

fn validate_genotype_semantics(
    genotype: &GenotypeRecord,
    reference: &Reference,
    copies: &BTreeMap<CopyIdentity, CopySummary<'_>>,
) -> Result<(), TruthValidationError> {
    match genotype.case_kind {
        GenotypeCaseKind::Snp => validate_snp_genotype(genotype, reference, copies),
        GenotypeCaseKind::Indel => validate_indel_genotype(genotype, copies),
        GenotypeCaseKind::Repeat => validate_repeat_genotype(genotype, reference, copies),
    }
}

fn validate_snp_genotype(
    genotype: &GenotypeRecord,
    reference: &Reference,
    copies: &BTreeMap<CopyIdentity, CopySummary<'_>>,
) -> Result<(), TruthValidationError> {
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == genotype.locus.sequence_id)
        .ok_or_else(|| mismatch(genotype, "la secuencia del locus no existe"))?;
    let position = usize::try_from(genotype.locus.start)
        .ok()
        .filter(|position| *position < sequence.length() as usize)
        .ok_or_else(|| mismatch(genotype, "el locus SNP sale de la referencia"))?;
    let reference_base = bases_of(sequence)?[position];
    let mut expected = BTreeMap::new();
    for (identity, copy) in copies.iter().filter(|(identity, _)| {
        identity.reference_id == genotype.reference_id
            && identity.individual_id == genotype.individual_id
            && identity.sequence_id == genotype.locus.sequence_id
    }) {
        let Some(slice) = copy
            .source_slices
            .iter()
            .find(|slice| slice.start <= genotype.locus.start && genotype.locus.start < slice.end)
            .filter(|slice| !slice.replaced_by_repeat)
        else {
            continue;
        };
        let base = match slice
            .snps
            .binary_search_by_key(&genotype.locus.start, |snp| snp.position)
        {
            Ok(index) => slice.snps[index].alternate.as_bytes()[0],
            Err(_) => reference_base,
        };
        let allele_index = genotype
            .alleles
            .iter()
            .position(|allele| match allele {
                TruthAllele::Reference => base == reference_base,
                TruthAllele::Snp { base: allele } => allele.as_bytes() == [base],
                _ => false,
            })
            .ok_or_else(|| {
                mismatch(
                    genotype,
                    "un alelo SNP reconstruido no aparece en el genotipo",
                )
            })?;
        expected.insert(identity.clone(), allele_index);
    }
    carriers_match(genotype, expected)
}

/// Contrasta un caso de repetición con las copias reconstruidas.
///
/// El alelo de índice cero es el de la **referencia** y no lleva secuencia: en
/// un caso de repetición significa «esta copia conserva el tracto que hay en la
/// referencia». Los alternativos llevan su secuencia y su estructura en
/// bloques, y se reconocen por ellas y no por ningún contador.
fn validate_repeat_genotype(
    genotype: &GenotypeRecord,
    reference: &Reference,
    copies: &BTreeMap<CopyIdentity, CopySummary<'_>>,
) -> Result<(), TruthValidationError> {
    let locus_id = genotype
        .alleles
        .iter()
        .find_map(|allele| match allele {
            TruthAllele::Repeat { locus_id, .. } => Some(locus_id),
            // Enumerada y no con comodín: si mañana se añade una clase de
            // alelo, el compilador obliga a decidir si pertenece a un caso de
            // repetición en vez de absorberla en silencio.
            TruthAllele::Reference
            | TruthAllele::Snp { .. }
            | TruthAllele::Insertion { .. }
            | TruthAllele::Deletion { .. } => None,
        })
        .ok_or_else(|| {
            mismatch(
                genotype,
                "un caso de repetición publica al menos un alelo alternativo",
            )
        })?;
    if !genotype.alleles.iter().all(|allele| {
        matches!(allele, TruthAllele::Reference)
            || matches!(allele, TruthAllele::Repeat { locus_id: other, .. } if other == locus_id)
    }) {
        return Err(mismatch(
            genotype,
            "todos los alelos de un caso de repetición pertenecen a su locus",
        ));
    }
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == genotype.locus.sequence_id)
        .ok_or_else(|| mismatch(genotype, "la secuencia del locus no existe"))?;
    // **Inalcanzable** (DT-001, tanda 5) en x86_64: `usize` y `u64` tienen el mismo ancho.
    let start = usize::try_from(genotype.locus.start)
        .map_err(|_| mismatch(genotype, "el locus de repetición sale de la referencia"))?;
    let end = usize::try_from(genotype.locus.end)
        .ok()
        .filter(|end| *end <= sequence.length() as usize && start < *end)
        .ok_or_else(|| mismatch(genotype, "el locus de repetición sale de la referencia"))?;
    // Un alternativo que repitiera las bases de la referencia no sería un
    // alelo: sería el mismo tracto publicado dos veces.
    let reference_tract = &bases_of(sequence)?[start..end];
    if genotype.alleles.iter().any(|allele| {
        matches!(allele, TruthAllele::Repeat { bases, .. } if bases.as_bytes() == reference_tract)
    }) {
        return Err(mismatch(
            genotype,
            "un alelo alternativo repite el tracto de la referencia",
        ));
    }
    let mut expected = BTreeMap::new();
    // Las copias de este individuo y esta secuencia son un tramo contiguo del
    // mapa, que ordena por referencia, individuo, secuencia y copia: se toma
    // ese tramo en vez de filtrar todas las copias de todos los individuos.
    let copy_of = |copy_index| CopyIdentity {
        reference_id: genotype.reference_id.clone(),
        individual_id: genotype.individual_id.clone(),
        sequence_id: genotype.locus.sequence_id.clone(),
        copy_index,
    };
    for (identity, copy) in copies.range(copy_of(u8::MIN)..=copy_of(u8::MAX)) {
        let interval = (genotype.locus.start, genotype.locus.end);
        let applied = copy
            .repeat_history
            .get(&interval)
            .into_iter()
            .flatten()
            .map(|&position| &copy.history[position])
            .find_map(|event| match event {
                TruthHistoryEvent::RepeatAllele {
                    locus_id: applied_locus,
                    start,
                    end,
                    bases,
                    blocks,
                    ..
                } if applied_locus == locus_id
                    && *start == genotype.locus.start
                    && *end == genotype.locus.end =>
                {
                    Some((bases, blocks))
                }
                TruthHistoryEvent::RepeatAllele { .. }
                | TruthHistoryEvent::Snp { .. }
                | TruthHistoryEvent::Insertion { .. }
                | TruthHistoryEvent::Deletion { .. } => None,
            });
        let index = match applied {
            None => 0,
            Some((bases, blocks)) => genotype
                .alleles
                .iter()
                .position(|allele| {
                    matches!(
                        allele,
                        TruthAllele::Repeat { bases: catalogued, blocks: catalogued_blocks, .. }
                            if catalogued == bases && catalogued_blocks == blocks
                    )
                })
                .ok_or_else(|| {
                    mismatch(
                        genotype,
                        "la copia lleva un alelo de repetición que el caso no cataloga",
                    )
                })?,
        };
        expected.insert(identity.clone(), index);
    }
    carriers_match(genotype, expected)
}

fn validate_indel_genotype(
    genotype: &GenotypeRecord,
    copies: &BTreeMap<CopyIdentity, CopySummary<'_>>,
) -> Result<(), TruthValidationError> {
    let alternates = genotype
        .alleles
        .iter()
        .enumerate()
        .filter(|(_, allele)| !matches!(allele, TruthAllele::Reference))
        .collect::<Vec<_>>();
    if alternates.len() != 1 || !matches!(genotype.alleles.first(), Some(TruthAllele::Reference)) {
        return Err(mismatch(
            genotype,
            "F2 requiere REF y un único alelo indel alternativo",
        ));
    }
    let (alternate_index, alternate) = alternates[0];
    let structural_event = match alternate {
        TruthAllele::Insertion { event_id, .. } | TruthAllele::Deletion { event_id, .. } => {
            Some(event_id.as_str())
        }
        _ => None,
    };
    let mut expected = BTreeMap::new();
    for (identity, copy) in copies.iter().filter(|(identity, _)| {
        identity.reference_id == genotype.reference_id
            && identity.individual_id == genotype.individual_id
            && identity.sequence_id == genotype.locus.sequence_id
    }) {
        let active = match alternate {
            TruthAllele::Insertion { event_id, .. } => {
                let history = copy.history.iter().find(|event| matches!(event, TruthHistoryEvent::Insertion { event_id: historical, .. } if historical == event_id));
                match history {
                    Some(history @ TruthHistoryEvent::Insertion { .. }) => {
                        active_insertion_matches(copy, alternate, history)
                    }
                    _ => false,
                }
            }
            TruthAllele::Deletion {
                event_id,
                start,
                end,
                source_event,
                piece_links,
            } => {
                let history = copy.history.iter().any(|event| matches!(event, TruthHistoryEvent::Deletion { event_id: historical, start: historical_start, end: historical_end, .. } if historical == event_id && historical_start == start && historical_end == end));
                history && deletion_matches(source_event, piece_links)
            }
            _ => false,
        };
        // El locus del caso está en coordenadas de la referencia: tiene que
        // ser el que la réplica de esta copia asigna a su evento, no las
        // coordenadas vigentes con que se declaró.
        if active
            && let Some(event_id) = structural_event
            && copy.event_loci.get(event_id) != Some(&(genotype.locus.start, genotype.locus.end))
        {
            return Err(mismatch(
                genotype,
                "el locus del caso no es el que su evento ocupa en la referencia",
            ));
        }
        expected.insert(identity.clone(), if active { alternate_index } else { 0 });
    }
    carriers_match(genotype, expected)
}

fn active_insertion_matches(
    copy: &CopySummary<'_>,
    allele: &TruthAllele,
    history: &TruthHistoryEvent,
) -> bool {
    let TruthAllele::Insertion {
        event_id,
        bases,
        source_event,
        piece_links,
    } = allele
    else {
        return false;
    };
    let TruthHistoryEvent::Insertion {
        point,
        bases: historical_bases,
        provenance,
        ..
    } = history
    else {
        return false;
    };
    if historical_bases != bases
        || !matches!(source_event, TruthStructuralEvent::Insertion { event_id: source_id, point: source_point, bases: source_bases } if source_id == event_id && *source_point == *point && source_bases == bases)
    {
        return false;
    }
    let expected = copy
        .pieces
        .iter()
        .enumerate()
        .filter_map(|(index, piece)| {
            (piece.inserted_provenance.as_ref() == Some(provenance)).then_some(TruthPieceLink {
                piece_index: u64::try_from(index).expect("piece index fits u64"),
                start: piece.start,
                end: piece.end,
            })
        })
        .collect::<Vec<_>>();
    expected.as_slice() == piece_links.as_slice()
}

/// Una deleción no deja fragmentos propios en la copia reconstruida, de modo
/// que no puede enlazar ninguno.
fn deletion_matches(source_event: &TruthStructuralEvent, piece_links: &[TruthPieceLink]) -> bool {
    matches!(source_event, TruthStructuralEvent::Deletion { .. }) && piece_links.is_empty()
}

fn carriers_match(
    genotype: &GenotypeRecord,
    expected: BTreeMap<CopyIdentity, usize>,
) -> Result<(), TruthValidationError> {
    let observed = genotype
        .carriers
        .iter()
        .map(|carrier| (carrier.copy.clone(), carrier.allele_index))
        .collect::<BTreeMap<_, _>>();
    if observed != expected {
        return Err(mismatch(
            genotype,
            "las copias portadoras no coinciden con la reconstrucción",
        ));
    }
    Ok(())
}

#[track_caller]
fn mismatch(genotype: &GenotypeRecord, detail: &'static str) -> TruthValidationError {
    crate::guard_trace::record();
    TruthValidationError::GenotypeMismatch {
        case_id: genotype.case_id.clone(),
        detail,
    }
}

/// Loci de la matriz de SNV en que **cada** individuo de la población publica
/// un caso SNP, aunque ninguna de sus copias lleve el alternativo (DTC-061: la
/// matriz «conserva también loci todo-REF»).
///
/// `truth.jsonl` no lleva las filas de la matriz, solo su huella: las filas
/// se publican en su propio artefacto (`population_snvs.tsv` y sus metadatos,
/// ligados por SHA-256 en el manifiesto y en `group_digests`), y la población
/// de cada individuo está en la configuración efectiva. Por eso este dato no
/// se lee del documento: se lo dan al juez la publicación de la verdad, desde
/// la matriz preparada en memoria, y el juez de los recursos publicados, desde
/// la matriz del paquete. Con él, el juez exige el caso de cada fila, rechaza
/// cualquier caso SNP que no sea ni fila ni cambio de una copia, exige los
/// alelos exactos de cada caso SNP y comprueba la población de las copias y
/// el identificador de cada caso. Sin él —al derivar el VCF o al leer solo
/// `truth.jsonl`— nada de eso puede comprobarse, y un caso SNP sin cambio en
/// ninguna copia se acepta.
#[derive(Clone, Debug, Default)]
pub(crate) struct MatrixCaseLoci {
    /// Alternativos de cada fila, por población y por secuencia y posición en
    /// base 0. Un locus multialélico ocupa varias filas.
    loci_by_population: BTreeMap<String, BTreeMap<(String, u64), BTreeSet<String>>>,
    population_by_individual: BTreeMap<String, String>,
}

impl MatrixCaseLoci {
    /// Las filas de la matriz preparada y la población de cada individuo.
    ///
    /// Es lo que necesitan los dos jueces que conocen la matriz —la
    /// publicación de la verdad y el de los recursos publicados—, y los dos
    /// la construyen aquí, con una sola regla para pasar la posición de la
    /// fila, en base 1, a la del caso, en base 0.
    pub(crate) fn from_matrix<'a>(
        matrix: &PreparedPopulationSnvs,
        individuals: impl IntoIterator<Item = (&'a str, &'a str)>,
    ) -> Self {
        let mut loci = Self::default();
        for row in matrix.rows() {
            loci.loci_by_population
                .entry(row.population_id().as_str().to_owned())
                .or_default()
                .entry((
                    row.contig().as_str().to_owned(),
                    row.position_1based().saturating_sub(1),
                ))
                .or_default()
                .insert(char::from(row.alternate()).to_string());
        }
        loci.add_populations(individuals);
        loci
    }

    /// Solo la población de cada individuo, para una simulación sin matriz.
    pub(crate) fn from_populations<'a>(
        individuals: impl IntoIterator<Item = (&'a str, &'a str)>,
    ) -> Self {
        let mut loci = Self::default();
        loci.add_populations(individuals);
        loci
    }

    fn add_populations<'a>(&mut self, individuals: impl IntoIterator<Item = (&'a str, &'a str)>) {
        for (individual_id, population_id) in individuals {
            self.population_by_individual
                .insert(individual_id.to_owned(), population_id.to_owned());
        }
    }
}

/// Lo que ata un caso a su fuente, dentro de un individuo: el locus de un
/// SNP, el evento de un *indel* o el intervalo de una repetición. Dos casos
/// del mismo individuo con la misma clave serían el mismo caso publicado dos
/// veces.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
enum CaseAnchor<'a> {
    Snp {
        sequence: &'a str,
        position: u64,
    },
    Event {
        sequence: &'a str,
        event_id: &'a str,
    },
    Repeat {
        sequence: &'a str,
        start: u64,
        end: u64,
    },
}

fn case_anchor(genotype: &GenotypeRecord) -> Result<CaseAnchor<'_>, TruthValidationError> {
    let sequence = genotype.locus.sequence_id.as_str();
    match genotype.case_kind {
        GenotypeCaseKind::Snp => Ok(CaseAnchor::Snp {
            sequence,
            position: genotype.locus.start,
        }),
        GenotypeCaseKind::Repeat => Ok(CaseAnchor::Repeat {
            sequence,
            start: genotype.locus.start,
            end: genotype.locus.end,
        }),
        GenotypeCaseKind::Indel => genotype
            .alleles
            .iter()
            .find_map(|allele| match allele {
                TruthAllele::Insertion { event_id, .. }
                | TruthAllele::Deletion { event_id, .. } => Some(CaseAnchor::Event {
                    sequence,
                    event_id: event_id.as_str(),
                }),
                TruthAllele::Reference | TruthAllele::Snp { .. } | TruthAllele::Repeat { .. } => {
                    None
                }
            })
            // No se ejercita: el contrato estructural rechaza antes un caso
            // indel sin alelo de inserción ni de deleción. Es defensiva.
            .ok_or_else(|| mismatch(genotype, "un caso indel nombra su evento estructural")),
    }
}

/// **El conjunto de casos**, no cada caso: la verdad publica todos los casos
/// que debe y ninguno más (DT-020 y DT-021).
///
/// El contraste caso a caso (`validate_genotype_semantics`) comprueba que cada
/// caso **que está** concuerda con las copias. Esto comprueba lo que aquel no
/// puede ver, individuo a individuo:
///
/// - **del cambio al caso**: cada diferencia de una copia con la referencia
///   —cambio de base, bases insertadas, bases de referencia que faltan, tracto
///   de repetición—, leída de los fragmentos que el juez replica por su
///   cuenta, tiene su caso; también cada deleción del historial, aunque ya no
///   cambie la copia (una que borra entera una inserción anterior), y cada
///   fila de la matriz cuando se conoce;
/// - **del caso a su fuente**: cada caso es uno de esos —un cambio de alguna
///   copia del individuo, un evento de su historial o una fila de la matriz—,
///   y ninguno se publica dos veces.
///
/// Un caso SNP que no es cambio de ninguna copia solo puede ser una fila de la
/// matriz; si `matrix` no se conoce, se acepta sin poder comprobarlo.
///
/// Recibe los casos y las copias de **un** individuo; quien la llama agrupa.
fn validate_individual_case_set<'a>(
    individual_id: &str,
    cases: &[&'a GenotypeRecord],
    copies: &[(&'a CopyIdentity, &'a CopySummary<'a>)],
    matrix: Option<&'a MatrixCaseLoci>,
) -> Result<(), TruthValidationError> {
    let mut published = BTreeMap::<CaseAnchor<'a>, (&'a GenotypeRecord, bool)>::new();
    for genotype in cases {
        if published
            .insert(case_anchor(genotype)?, (*genotype, false))
            .is_some()
        {
            return Err(mismatch(
                genotype,
                "otro caso del individuo publica ya ese locus o ese evento",
            ));
        }
    }
    let missing = |sequence: &str, start: u64, end: u64, detail: &'static str| {
        TruthValidationError::MissingCase {
            individual_id: individual_id.to_owned(),
            sequence_id: sequence.to_owned(),
            start,
            end,
            detail,
        }
    };
    let mut require =
        |anchor: CaseAnchor<'a>, start: u64, end: u64, detail: &'static str| match published
            .get_mut(&anchor)
        {
            Some((_, required)) => {
                *required = true;
                Ok(())
            }
            None => {
                let sequence = match anchor {
                    CaseAnchor::Snp { sequence, .. }
                    | CaseAnchor::Event { sequence, .. }
                    | CaseAnchor::Repeat { sequence, .. } => sequence,
                };
                Err(rejected(missing(sequence, start, end, detail)))
            }
        };
    for (identity, copy) in copies {
        let sequence = identity.sequence_id.as_str();
        for slice in &copy.source_slices {
            if slice.replaced_by_repeat {
                require(
                    CaseAnchor::Repeat {
                        sequence,
                        start: slice.start,
                        end: slice.end,
                    },
                    slice.start,
                    slice.end,
                    "un tracto de repetición de la copia no tiene caso",
                )?;
            } else {
                for snp in &slice.snps {
                    require(
                        CaseAnchor::Snp {
                            sequence,
                            position: snp.position,
                        },
                        snp.position,
                        snp.position.saturating_add(1),
                        "un cambio de base de la copia no tiene caso",
                    )?;
                }
            }
        }
        // Inserciones del historial por la fuente de su procedencia, para no
        // recorrer el historial entero por cada fragmento insertado.
        let mut insertions = BTreeMap::<&TruthSourceRef, Vec<_>>::new();
        for event in copy.history.iter() {
            if let TruthHistoryEvent::Insertion {
                event_id,
                provenance,
                ..
            } = event
            {
                insertions
                    .entry(&provenance.source)
                    .or_default()
                    .push((provenance, event_id.as_str()));
            }
        }
        for piece in copy.pieces.iter().filter(|piece| !piece.tract) {
            let Some(provenance) = &piece.inserted_provenance else {
                continue;
            };
            let event_id = insertions.get(&provenance.source).and_then(|candidates| {
                candidates
                    .iter()
                    .find(|(historical, _)| *historical == provenance)
                    .map(|(_, event_id)| *event_id)
            });
            // Bases insertadas sin evento en el historial solo pueden venir de
            // un segmento heredado, y ningún caso puede explicarlas: un caso
            // indel exige su evento en el historial.
            let event_id = event_id.ok_or_else(|| {
                rejected(missing(
                    sequence,
                    piece.start,
                    piece.end,
                    "unas bases insertadas de la copia no tienen evento ni caso",
                ))
            })?;
            let (start, end) = copy.event_loci.get(event_id).copied().unwrap_or_default();
            require(
                CaseAnchor::Event { sequence, event_id },
                start,
                end,
                "unas bases insertadas de la copia no tienen caso",
            )?;
        }
        for event in copy.history.iter() {
            if let TruthHistoryEvent::Deletion { event_id, .. } = event {
                let (start, end) = copy.event_loci.get(event_id).copied().unwrap_or_default();
                require(
                    CaseAnchor::Event {
                        sequence,
                        event_id: event_id.as_str(),
                    },
                    start,
                    end,
                    "una deleción del historial no tiene caso",
                )?;
            }
        }
        if let Some((start, end)) = copy.inherited_loss {
            return Err(rejected(missing(
                sequence,
                start,
                end,
                "faltan bases de referencia en el estado heredado de los segmentos, sin evento ni caso",
            )));
        }
    }
    let mut population = None;
    if let Some(matrix) = matrix {
        // Un individuo del documento que no es de la simulación. Al proyectar
        // no se alcanza —el contraste de las copias con lo simulado lo
        // rechaza antes—, pero en el juez de los recursos publicados sí: un
        // paquete con un individuo renombrado y sus identificadores de caso
        // rehechos llega aquí, y lo ejercita
        // `published_truth_is_judged_with_the_published_matrix`.
        let known = matrix
            .population_by_individual
            .get(individual_id)
            .ok_or_else(|| {
                rejected(missing(
                    "",
                    0,
                    0,
                    "el individuo no tiene población conocida",
                ))
            })?;
        population = Some(known.as_str());
        let rows = matrix.loci_by_population.get(known.as_str());
        for (sequence, position) in rows.into_iter().flat_map(BTreeMap::keys) {
            require(
                CaseAnchor::Snp {
                    sequence,
                    position: *position,
                },
                *position,
                position.saturating_add(1),
                "una fila de la matriz no tiene el caso del individuo",
            )?;
        }
    }
    for (anchor, (genotype, required)) in &published {
        if *required || (matrix.is_none() && matches!(anchor, CaseAnchor::Snp { .. })) {
            continue;
        }
        return Err(mismatch(
            genotype,
            "el caso no corresponde a ningún cambio de las copias, evento del historial ni fila de la matriz",
        ));
    }
    if let Some((matrix, population)) = matrix.zip(population) {
        let rows = matrix.loci_by_population.get(population);
        // La población que declaran las copias es la del individuo.
        let origin = format!("{COPY_ORIGIN_PREFIX}{population}");
        if let Some((identity, _)) = copies.iter().find(|(_, copy)| copy.origin != origin) {
            return Err(rejected(TruthValidationError::CopyOriginMismatch {
                identity: (*identity).clone(),
                expected: origin,
            }));
        }
        for genotype in cases {
            if genotype.case_kind == GenotypeCaseKind::Snp {
                let rows = rows.and_then(|rows| {
                    rows.get(&(genotype.locus.sequence_id.clone(), genotype.locus.start))
                });
                if genotype.alleles != expected_snp_alleles(genotype, copies, rows) {
                    return Err(mismatch(
                        genotype,
                        "los alelos del caso SNP no son los de la fila de la matriz y los observados en las copias",
                    ));
                }
            }
            if genotype.case_kind == GenotypeCaseKind::Repeat
                && genotype.alleles != expected_repeat_alleles(genotype, copies)
            {
                return Err(mismatch(
                    genotype,
                    "los alelos del caso de repetición no son los de sus copias en su orden",
                ));
            }
            if genotype.case_id != expected_case_id(population, genotype) {
                return Err(mismatch(
                    genotype,
                    "el identificador del caso no es el que dan su población, individuo, locus y alelos",
                ));
            }
        }
    }
    Ok(())
}

/// Prefijo de `provenance.origin` en las copias de una población inicial; le
/// sigue el identificador de la población (`schemas/truth-v4.md`).
const COPY_ORIGIN_PREFIX: &str = "initial_population:";

/// Catálogo exacto de un caso SNP: la referencia; después, si el locus es una
/// fila de la matriz, sus alternativos en orden de base, y a continuación los
/// observados en las copias que la fila no tiene, en el orden en que aparecen
/// de la primera copia a la última; si no es fila, solo los observados, en
/// orden de base.
///
/// La rama de «observados que la fila no tiene» **no se alcanza** desde una
/// configuración: un cambio de base directo en una fila de la matriz con un
/// alternativo que la fila no tiene se rechaza antes de simular
/// (`ExplicitSnpAlleleOutsideMatrix`, comprobado ejecutándolo). Se conserva
/// porque es la regla del productor, y ninguna prueba distingue si esos
/// observados se añaden o no.
fn expected_snp_alleles(
    genotype: &GenotypeRecord,
    copies: &[(&CopyIdentity, &CopySummary<'_>)],
    rows: Option<&BTreeSet<String>>,
) -> Vec<TruthAllele> {
    let mut observed = Vec::<String>::new();
    for (identity, copy) in copies {
        if identity.sequence_id != genotype.locus.sequence_id {
            continue;
        }
        let alternate = copy
            .source_slices
            .iter()
            .find(|slice| {
                !slice.replaced_by_repeat
                    && slice.start <= genotype.locus.start
                    && genotype.locus.start < slice.end
            })
            .and_then(|slice| {
                slice
                    .snps
                    .binary_search_by_key(&genotype.locus.start, |snp| snp.position)
                    .ok()
                    .map(|index| slice.snps[index].alternate.clone())
            });
        if let Some(alternate) = alternate
            && !observed.contains(&alternate)
        {
            observed.push(alternate);
        }
    }
    let alternates = match rows {
        Some(rows) => {
            let mut alternates = rows.iter().cloned().collect::<Vec<_>>();
            alternates.extend(observed.into_iter().filter(|base| !rows.contains(base)));
            alternates
        }
        None => {
            observed.sort_unstable();
            observed
        }
    };
    std::iter::once(TruthAllele::Reference)
        .chain(alternates.into_iter().map(|base| TruthAllele::Snp { base }))
        .collect()
}

/// Catálogo exacto de un caso de repetición: la referencia y después cada
/// tracto distinto que llevan las copias del individuo en ese intervalo, en el
/// orden en que aparece de la primera copia a la última.
fn expected_repeat_alleles(
    genotype: &GenotypeRecord,
    copies: &[(&CopyIdentity, &CopySummary<'_>)],
) -> Vec<TruthAllele> {
    let mut alleles = vec![TruthAllele::Reference];
    for (identity, copy) in copies {
        if identity.sequence_id != genotype.locus.sequence_id {
            continue;
        }
        let tract = copy.history.iter().find_map(|event| match event {
            TruthHistoryEvent::RepeatAllele {
                locus_id,
                start,
                end,
                bases,
                blocks,
                ..
            } if *start == genotype.locus.start && *end == genotype.locus.end => {
                Some(TruthAllele::Repeat {
                    locus_id: locus_id.clone(),
                    bases: bases.clone(),
                    blocks: blocks.clone(),
                })
            }
            TruthHistoryEvent::RepeatAllele { .. }
            | TruthHistoryEvent::Snp { .. }
            | TruthHistoryEvent::Insertion { .. }
            | TruthHistoryEvent::Deletion { .. } => None,
        });
        if let Some(tract) = tract
            && !alleles.contains(&tract)
        {
            alleles.push(tract);
        }
    }
    alleles
}

/// El `case_id` que corresponde a un caso: `CASE_` seguido de los 59
/// primeros dígitos hexadecimales del SHA-256 de población, individuo,
/// secuencia, posición en base 1 del inicio del locus y la lista de claves de
/// sus alelos separadas por comas, cada campo terminado en un byte nulo
/// (`schemas/truth-v4.md`). Implementación propia del juez.
fn expected_case_id(population: &str, genotype: &GenotypeRecord) -> String {
    let alleles = genotype
        .alleles
        .iter()
        .map(|allele| match allele {
            TruthAllele::Reference => "REF".to_owned(),
            TruthAllele::Snp { base } => format!("SNP:{base}"),
            TruthAllele::Insertion {
                event_id, bases, ..
            } => format!("INS:{event_id}:{bases}"),
            TruthAllele::Deletion {
                event_id,
                start,
                end,
                ..
            } => format!("DEL:{event_id}:{start}:{end}"),
            TruthAllele::Repeat {
                locus_id, bases, ..
            } => format!("STR:{locus_id}:{bases}"),
        })
        .collect::<Vec<_>>()
        .join(",");
    let position = genotype.locus.start.saturating_add(1).to_string();
    let mut hasher = Sha256::new();
    for field in [
        population,
        genotype.individual_id.as_str(),
        genotype.locus.sequence_id.as_str(),
        position.as_str(),
        alleles.as_str(),
    ] {
        hasher.update(field.as_bytes());
        hasher.update([0]);
    }
    let digest = format!("{:x}", hasher.finalize());
    format!("CASE_{}", &digest[..59])
}

/// Agrupa por individuo los casos y las copias de un documento completo y
/// juzga el conjunto de casos de cada uno.
fn validate_case_set<'truth>(
    genotypes: impl Iterator<Item = &'truth GenotypeRecord>,
    copies: &BTreeMap<CopyIdentity, CopySummary<'_>>,
    matrix: Option<&MatrixCaseLoci>,
) -> Result<(), TruthValidationError> {
    let mut cases = BTreeMap::<&str, Vec<&GenotypeRecord>>::new();
    for genotype in genotypes {
        cases
            .entry(genotype.individual_id.as_str())
            .or_default()
            .push(genotype);
    }
    let mut by_individual = BTreeMap::<&str, Vec<(&CopyIdentity, &CopySummary<'_>)>>::new();
    for (identity, copy) in copies {
        by_individual
            .entry(identity.individual_id.as_str())
            .or_default()
            .push((identity, copy));
    }
    let individuals = cases
        .keys()
        .chain(by_individual.keys())
        .copied()
        .collect::<BTreeSet<_>>();
    for individual_id in individuals {
        validate_individual_case_set(
            individual_id,
            cases.get(individual_id).map_or(&[], Vec::as_slice),
            by_individual.get(individual_id).map_or(&[], Vec::as_slice),
            matrix,
        )?;
    }
    Ok(())
}

#[derive(Deserialize)]
#[serde(tag = "record_type", rename_all = "snake_case", deny_unknown_fields)]
enum WireLine {
    Header(WireHeader),
    HaplotypeCopy(WireCopy),
    Genotype(WireGenotype),
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireHeader {
    document_kind: String,
    schema_version: String,
    execution_id: String,
    reference_id: String,
    reference_content_sha256: String,
    coordinate_convention: String,
    counts: WireCounts,
    source_digests: Vec<WireDigest>,
    group_digests: Vec<WireDigest>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireCounts {
    #[serde(deserialize_with = "deserialize_u64")]
    copies: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    genotypes: u64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireDigest {
    source: WireSource,
    sha256: String,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct WireSource {
    kind: WireSourceKind,
    source_id: String,
}
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum WireSourceKind {
    Reference,
    Founder,
    Explicit,
    Random,
    Matrix,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireCopy {
    identity: WireIdentity,
    #[serde(deserialize_with = "deserialize_u8")]
    phase: u8,
    provenance: WireCopyProvenance,
    segments: Vec<WireSegment>,
    current_pieces: Vec<WireCurrentPiece>,
    history: Vec<WireHistoryEvent>,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct WireIdentity {
    reference_id: String,
    individual_id: String,
    sequence_id: String,
    #[serde(deserialize_with = "deserialize_u8")]
    copy_index: u8,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireCopyProvenance {
    origin: String,
    sources: Vec<WireSource>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireSegment {
    source: WireSegmentSource,
    modifications: Vec<WireModification>,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct WireSegmentSource {
    source: WireSource,
    reference_id: String,
    sequence_id: String,
    #[serde(deserialize_with = "deserialize_u64")]
    start: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    end: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    sequence_length: u64,
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum WireCurrentPiece {
    Source {
        source: WireSegmentSource,
        snps: Vec<WireCurrentSnp>,
    },
    Inserted {
        bases: String,
        provenance: WireProvenance,
    },
    RepeatAllele {
        locus_id: String,
        source: WireSegmentSource,
        bases: String,
        blocks: Vec<WireRepeatBlock>,
        provenance: WireProvenance,
    },
}

/// Bloque leído del documento.
///
/// El documento publica cuatro campos y **dos de ellos son derivados**: la
/// longitud del motivo es la de su secuencia y el número de unidades es las
/// bases entre esa longitud. Aquí se **comprueban** en vez de ignorarse, y solo
/// se conservan los dos que no se derivan de nada. Un documento cuyas listas
/// estuvieran desalineadas —el error silencioso más probable de estos
/// campos— se rechaza al leerlo.
struct WireRepeatBlock {
    motif: String,
    bases: u64,
}

impl<'de> Deserialize<'de> for WireRepeatBlock {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        #[derive(Deserialize)]
        #[serde(deny_unknown_fields)]
        struct Raw {
            motif: String,
            #[serde(deserialize_with = "deserialize_u64")]
            motif_length: u64,
            #[serde(deserialize_with = "deserialize_u64")]
            bases: u64,
            // `Number` y no `f64`: serde almacena el contenido de una
            // enumeración etiquetada por dentro antes de repartirlo, y desde
            // ese almacén un `f64` no se puede leer —falla con «invalid type:
            // map»—. `serde_json::Number` sí, y es lo que el resto del lector
            // usa para los enteros por la misma razón.
            units: Number,
        }
        let raw = Raw::deserialize(deserializer)?;
        let units = raw
            .units
            .as_f64()
            .ok_or_else(|| rejected(D::Error::custom("el número de unidades no es un decimal")))?;
        // **Inalcanzable** (DT-001, tanda 5): `usize` a `u64` no falla en x86_64.
        let observed = u64::try_from(raw.motif.len())
            .map_err(|_| rejected(D::Error::custom("el motivo del bloque es inabarcable")))?;
        if observed == 0 || raw.motif_length != observed {
            return Err(rejected(D::Error::custom(
                "la longitud declarada del motivo no es la de su secuencia",
            )));
        }
        if raw.bases == 0 {
            return Err(rejected(D::Error::custom(
                "un bloque sin bases no es un bloque",
            )));
        }
        let derived = raw.bases as f64 / observed as f64;
        if (units - derived).abs() > 1e-9 * derived.max(1.0) {
            return Err(rejected(D::Error::custom(
                "el número de unidades no es las bases entre la longitud del motivo",
            )));
        }
        Ok(Self {
            motif: raw.motif,
            bases: raw.bases,
        })
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireCurrentSnp {
    event_id: String,
    #[serde(deserialize_with = "deserialize_u64")]
    position: u64,
    reference: String,
    alternate: String,
    provenance: WireProvenance,
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum WireHistoryEvent {
    Snp {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        position: u64,
        provenance: WireProvenance,
    },
    RepeatAllele {
        event_id: String,
        locus_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        start: u64,
        #[serde(deserialize_with = "deserialize_u64")]
        end: u64,
        bases: String,
        blocks: Vec<WireRepeatBlock>,
        provenance: WireProvenance,
    },
    Insertion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        point: u64,
        bases: String,
        provenance: WireProvenance,
    },
    Deletion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        start: u64,
        #[serde(deserialize_with = "deserialize_u64")]
        end: u64,
        provenance: WireProvenance,
    },
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum WireModification {
    Snp {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        position: u64,
        reference: String,
        alternate: String,
        provenance: WireProvenance,
    },
    Insertion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        point: u64,
        bases: String,
        provenance: WireProvenance,
    },
    Deletion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        start: u64,
        #[serde(deserialize_with = "deserialize_u64")]
        end: u64,
        provenance: WireProvenance,
    },
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct WireProvenance {
    kind: WireProvenanceKind,
    version: String,
    source: WireSource,
    matrix: Option<WireMatrixProvenance>,
}
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum WireProvenanceKind {
    Explicit,
    Random,
    Matrix,
}
#[derive(Clone, Deserialize)]
#[serde(deny_unknown_fields)]
struct WireMatrixProvenance {
    matrix_source_id: String,
    population_id: String,
    variant_id: String,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireGenotype {
    case_id: String,
    reference_id: String,
    individual_id: String,
    case_kind: WireGenotypeKind,
    locus: WireLocus,
    alleles: Vec<WireAllele>,
    #[serde(deserialize_with = "deserialize_vec_u8")]
    allele_dosages: Vec<u8>,
    carriers: Vec<WireCarrier>,
}
#[derive(Clone, Copy, Deserialize)]
#[serde(rename_all = "snake_case")]
enum WireGenotypeKind {
    Snp,
    Indel,
    Repeat,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireLocus {
    sequence_id: String,
    #[serde(deserialize_with = "deserialize_u64")]
    start: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    end: u64,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WireCarrier {
    #[serde(deserialize_with = "deserialize_usize")]
    allele_index: usize,
    copy: WireIdentity,
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum WireAllele {
    // Variante con llaves vacías a propósito (DT-028): una variante sin campos
    // de un enum etiquetado por dentro no aplica `deny_unknown_fields`, y el
    // esquema de la verdad cierra también el alelo de referencia.
    Reference {},
    Snp {
        base: String,
    },
    Repeat {
        locus_id: String,
        bases: String,
        blocks: Vec<WireRepeatBlock>,
    },
    Insertion {
        event_id: String,
        bases: String,
        source_event: WireStructuralEvent,
        piece_links: Vec<WirePieceLink>,
    },
    Deletion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        start: u64,
        #[serde(deserialize_with = "deserialize_u64")]
        end: u64,
        source_event: WireStructuralEvent,
        piece_links: Vec<WirePieceLink>,
    },
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum WireStructuralEvent {
    Insertion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        point: u64,
        bases: String,
    },
    Deletion {
        event_id: String,
        #[serde(deserialize_with = "deserialize_u64")]
        start: u64,
        #[serde(deserialize_with = "deserialize_u64")]
        end: u64,
    },
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct WirePieceLink {
    #[serde(deserialize_with = "deserialize_u64")]
    piece_index: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    start: u64,
    #[serde(deserialize_with = "deserialize_u64")]
    end: u64,
}

/// Nombres de campo cuyo valor es una lista en el esquema de la verdad
/// (`schemas/truth-v4.schema.json` y `truth-v5.schema.json`). **Cualquier
/// otro valor que sea una lista se rechaza** (DT-032).
///
/// `serde` acepta, para una estructura, tanto un objeto como una secuencia con
/// los valores en el orden de los campos (`"counts":[6,8]`), y
/// `deny_unknown_fields` no lo impide. En vez de cerrar cada tipo, se
/// recorre la línea una vez sin construir nada y se exige que una lista solo
/// aparezca bajo uno de estos nombres y que sus elementos sean objetos o
/// escalares, nunca otra lista. Un tipo nuevo queda cubierto sin tocar nada; un
/// campo de lista nuevo hay que añadirlo aquí, y la prueba
/// `list_field_names_are_exactly_the_schema_array_properties` lo exige.
const TRUTH_LIST_FIELDS: [&str; 13] = [
    "allele_dosages",
    "alleles",
    "blocks",
    "carriers",
    "current_pieces",
    "group_digests",
    "history",
    "modifications",
    "piece_links",
    "segments",
    "snps",
    "source_digests",
    "sources",
];

/// Recorre un valor JSON sin conservarlo y rechaza toda lista que no esté en
/// posición de lista. Ver [`TRUTH_LIST_FIELDS`].
struct ObjectShape {
    list_allowed: bool,
}

impl<'de> DeserializeSeed<'de> for ObjectShape {
    type Value = ();
    fn deserialize<D: Deserializer<'de>>(self, deserializer: D) -> Result<(), D::Error> {
        deserializer.deserialize_any(self)
    }
}

impl<'de> Visitor<'de> for ObjectShape {
    type Value = ();

    fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str("a JSON value")
    }

    fn visit_bool<E>(self, _: bool) -> Result<(), E> {
        Ok(())
    }
    fn visit_i64<E>(self, _: i64) -> Result<(), E> {
        Ok(())
    }
    fn visit_u64<E>(self, _: u64) -> Result<(), E> {
        Ok(())
    }
    fn visit_f64<E>(self, _: f64) -> Result<(), E> {
        Ok(())
    }
    fn visit_str<E>(self, _: &str) -> Result<(), E> {
        Ok(())
    }
    fn visit_unit<E>(self) -> Result<(), E> {
        Ok(())
    }

    fn visit_map<A: MapAccess<'de>>(self, mut map: A) -> Result<(), A::Error> {
        while let Some(IsListField(list_allowed)) = map.next_key()? {
            map.next_value_seed(ObjectShape { list_allowed })?;
        }
        Ok(())
    }

    fn visit_seq<A: SeqAccess<'de>>(self, mut seq: A) -> Result<(), A::Error> {
        if !self.list_allowed {
            return Err(rejected(A::Error::invalid_type(
                Unexpected::Seq,
                &"an object, not a list of its values",
            )));
        }
        while seq
            .next_element_seed(ObjectShape {
                list_allowed: false,
            })?
            .is_some()
        {}
        Ok(())
    }
}

/// Clave de objeto leída solo para saber si nombra un campo de lista.
struct IsListField(bool);

impl<'de> Deserialize<'de> for IsListField {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        struct KeyVisitor;
        impl Visitor<'_> for KeyVisitor {
            type Value = IsListField;
            fn expecting(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                formatter.write_str("an object key")
            }
            fn visit_str<E>(self, key: &str) -> Result<IsListField, E> {
                Ok(IsListField(TRUTH_LIST_FIELDS.contains(&key)))
            }
        }
        deserializer.deserialize_str(KeyVisitor)
    }
}

fn deserialize_u64<'de, D>(deserializer: D) -> Result<u64, D::Error>
where
    D: Deserializer<'de>,
{
    let number = Number::deserialize(deserializer)?;
    parse_json_integer(&number.to_string()).ok_or_else(|| {
        rejected(D::Error::custom(
            "se esperaba un entero JSON no negativo exacto en u64",
        ))
    })
}

fn deserialize_u8<'de, D>(deserializer: D) -> Result<u8, D::Error>
where
    D: Deserializer<'de>,
{
    let value = deserialize_u64(deserializer)?;
    u8::try_from(value).map_err(|_| rejected(D::Error::custom("entero fuera de u8")))
}

fn deserialize_usize<'de, D>(deserializer: D) -> Result<usize, D::Error>
where
    D: Deserializer<'de>,
{
    let value = deserialize_u64(deserializer)?;
    // **Inalcanzable** (DT-001, tanda 5) en x86_64: `usize` y `u64` tienen el mismo ancho.
    usize::try_from(value).map_err(|_| rejected(D::Error::custom("entero fuera de usize")))
}

fn deserialize_vec_u8<'de, D>(deserializer: D) -> Result<Vec<u8>, D::Error>
where
    D: Deserializer<'de>,
{
    let values = Vec::<Number>::deserialize(deserializer)?;
    values
        .into_iter()
        .map(|number| {
            let value = parse_json_integer(&number.to_string())
                .ok_or_else(|| rejected(D::Error::custom("dosis no es un entero exacto")))?;
            u8::try_from(value).map_err(|_| rejected(D::Error::custom("dosis fuera de u8")))
        })
        .collect()
}

/// JSON Schema considera `2.0` y `2e0` enteros. Aceptamos esas grafías solo
/// si su valor decimal es exacto y cabe en `u64`, sin convertir por `f64`.
fn parse_json_integer(text: &str) -> Option<u64> {
    let (negative, unsigned) = text
        .strip_prefix('-')
        .map_or((false, text), |value| (true, value));
    let (mantissa, exponent_text) = match unsigned.find(['e', 'E']) {
        Some(index) => (&unsigned[..index], Some(&unsigned[index + 1..])),
        None => (unsigned, None),
    };
    let (whole, fractional) = match mantissa.split_once('.') {
        Some(parts) => parts,
        None => (mantissa, ""),
    };
    let digits = format!("{whole}{fractional}");
    if digits.is_empty() || !digits.bytes().all(|byte| byte.is_ascii_digit()) {
        return None;
    }
    if let Some(exponent) = exponent_text
        && !valid_json_exponent(exponent)
    {
        return None;
    }
    // JSON Schema considera entero cualquier mantisa nula, incluida -0.0 o
    // 0e999...; no hace falta ni es seguro convertir su exponente a un tipo
    // finito para conocer el resultado exacto.
    if digits.bytes().all(|byte| byte == b'0') {
        return Some(0);
    }
    if negative {
        return None;
    }
    let exponent = match exponent_text {
        Some(text) => text.parse::<i64>().ok()?,
        None => 0,
    };
    let shift = exponent.checked_sub(i64::try_from(fractional.len()).ok()?)?;
    let significant = digits.trim_start_matches('0');
    if shift >= 0 {
        let shift = usize::try_from(shift).ok()?;
        // Un número no nulo con más de veinte cifras decimales no cabe en
        // u64. Esta comprobación impide recorrer exponentes hostiles.
        if significant.len().checked_add(shift)? > u64::MAX.to_string().len() {
            return None;
        }
        let mut value = significant.parse::<u64>().ok()?;
        for _ in 0..shift {
            value = value.checked_mul(10)?;
        }
        return Some(value);
    }
    let places = usize::try_from(shift.unsigned_abs()).ok()?;
    if places >= digits.len()
        || !digits.as_bytes()[digits.len() - places..]
            .iter()
            .all(|byte| *byte == b'0')
    {
        return None;
    }
    digits[..digits.len() - places]
        .trim_start_matches('0')
        .parse::<u64>()
        .ok()
}

fn valid_json_exponent(text: &str) -> bool {
    let digits = text
        .strip_prefix('+')
        .or_else(|| text.strip_prefix('-'))
        .unwrap_or(text);
    !digits.is_empty() && digits.bytes().all(|byte| byte.is_ascii_digit())
}

impl From<WireHeader> for TruthHeader {
    fn from(value: WireHeader) -> Self {
        Self {
            document_kind: value.document_kind,
            schema_version: value.schema_version,
            execution_id: value.execution_id,
            reference_id: value.reference_id,
            reference_content_sha256: value.reference_content_sha256,
            coordinate_convention: value.coordinate_convention,
            copy_count: value.counts.copies,
            genotype_count: value.counts.genotypes,
            source_digests: value.source_digests.into_iter().map(Into::into).collect(),
            group_digests: value.group_digests.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<WireDigest> for NamedDigest {
    fn from(value: WireDigest) -> Self {
        Self {
            source: value.source.into(),
            sha256: value.sha256,
        }
    }
}
impl From<WireSource> for TruthSourceRef {
    fn from(value: WireSource) -> Self {
        Self {
            kind: value.kind.into(),
            source_id: value.source_id,
        }
    }
}
impl From<WireSourceKind> for TruthSourceKind {
    fn from(value: WireSourceKind) -> Self {
        match value {
            WireSourceKind::Reference => Self::Reference,
            WireSourceKind::Founder => Self::Founder,
            WireSourceKind::Explicit => Self::Explicit,
            WireSourceKind::Random => Self::Random,
            WireSourceKind::Matrix => Self::Matrix,
        }
    }
}
impl From<WireIdentity> for CopyIdentity {
    fn from(value: WireIdentity) -> Self {
        Self {
            reference_id: value.reference_id,
            individual_id: value.individual_id,
            sequence_id: value.sequence_id,
            copy_index: value.copy_index,
        }
    }
}
impl From<WireCopy> for HaplotypeCopyRecord {
    fn from(value: WireCopy) -> Self {
        Self {
            identity: value.identity.into(),
            phase: value.phase,
            provenance: CopyProvenance {
                origin: value.provenance.origin,
                sources: value
                    .provenance
                    .sources
                    .into_iter()
                    .map(Into::into)
                    .collect(),
            },
            segments: value.segments.into_iter().map(Into::into).collect(),
            current_pieces: value.current_pieces.into_iter().map(Into::into).collect(),
            history: value.history.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<WireSegment> for crate::truth_schema::TruthSegment {
    fn from(value: WireSegment) -> Self {
        Self {
            source: value.source.into(),
            modifications: value.modifications.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<WireSegmentSource> for SegmentSource {
    fn from(value: WireSegmentSource) -> Self {
        Self {
            source: value.source.into(),
            reference_id: value.reference_id,
            sequence_id: value.sequence_id,
            start: value.start,
            end: value.end,
            sequence_length: value.sequence_length,
        }
    }
}
impl From<WireCurrentPiece> for TruthCurrentPiece {
    fn from(value: WireCurrentPiece) -> Self {
        match value {
            WireCurrentPiece::Source { source, snps } => Self::Source {
                source: source.into(),
                snps: snps.into_iter().map(Into::into).collect(),
            },
            WireCurrentPiece::Inserted { bases, provenance } => Self::Inserted {
                bases,
                provenance: provenance.into(),
            },
            WireCurrentPiece::RepeatAllele {
                locus_id,
                source,
                bases,
                blocks,
                provenance,
            } => Self::RepeatAllele {
                locus_id,
                source: source.into(),
                bases,
                blocks: blocks.into_iter().map(Into::into).collect(),
                provenance: provenance.into(),
            },
        }
    }
}
impl From<WireRepeatBlock> for TruthRepeatBlock {
    fn from(value: WireRepeatBlock) -> Self {
        Self {
            motif: value.motif,
            bases: value.bases,
        }
    }
}
impl From<WireCurrentSnp> for TruthCurrentSnp {
    fn from(value: WireCurrentSnp) -> Self {
        Self {
            event_id: value.event_id,
            position: value.position,
            reference: value.reference,
            alternate: value.alternate,
            provenance: value.provenance.into(),
        }
    }
}
impl From<WireHistoryEvent> for TruthHistoryEvent {
    fn from(value: WireHistoryEvent) -> Self {
        match value {
            WireHistoryEvent::Snp {
                event_id,
                position,
                provenance,
            } => Self::Snp {
                event_id,
                position,
                provenance: provenance.into(),
            },
            WireHistoryEvent::Insertion {
                event_id,
                point,
                bases,
                provenance,
            } => Self::Insertion {
                event_id,
                point,
                bases,
                provenance: provenance.into(),
            },
            WireHistoryEvent::Deletion {
                event_id,
                start,
                end,
                provenance,
            } => Self::Deletion {
                event_id,
                start,
                end,
                provenance: provenance.into(),
            },
            WireHistoryEvent::RepeatAllele {
                event_id,
                locus_id,
                start,
                end,
                bases,
                blocks,
                provenance,
            } => Self::RepeatAllele {
                event_id,
                locus_id,
                start,
                end,
                bases,
                blocks: blocks.into_iter().map(Into::into).collect(),
                provenance: provenance.into(),
            },
        }
    }
}
impl From<WireModification> for TruthModification {
    fn from(value: WireModification) -> Self {
        match value {
            WireModification::Snp {
                event_id,
                position,
                reference,
                alternate,
                provenance,
            } => Self::Snp {
                event_id,
                position,
                reference,
                alternate,
                provenance: provenance.into(),
            },
            WireModification::Insertion {
                event_id,
                point,
                bases,
                provenance,
            } => Self::Insertion {
                event_id,
                point,
                bases,
                provenance: provenance.into(),
            },
            WireModification::Deletion {
                event_id,
                start,
                end,
                provenance,
            } => Self::Deletion {
                event_id,
                start,
                end,
                provenance: provenance.into(),
            },
        }
    }
}
impl From<WireProvenance> for ModificationProvenance {
    fn from(value: WireProvenance) -> Self {
        Self {
            kind: value.kind.into(),
            version: value.version,
            source: value.source.into(),
            matrix: value.matrix.map(Into::into),
        }
    }
}
impl From<WireProvenanceKind> for ModificationProvenanceKind {
    fn from(value: WireProvenanceKind) -> Self {
        match value {
            WireProvenanceKind::Explicit => Self::Explicit,
            WireProvenanceKind::Random => Self::Random,
            WireProvenanceKind::Matrix => Self::Matrix,
        }
    }
}
impl From<WireMatrixProvenance> for MatrixProvenance {
    fn from(value: WireMatrixProvenance) -> Self {
        Self {
            matrix_source_id: value.matrix_source_id,
            population_id: value.population_id,
            variant_id: value.variant_id,
        }
    }
}
impl From<WireGenotype> for GenotypeRecord {
    fn from(value: WireGenotype) -> Self {
        Self {
            case_id: value.case_id,
            reference_id: value.reference_id,
            individual_id: value.individual_id,
            case_kind: value.case_kind.into(),
            locus: value.locus.into(),
            alleles: value.alleles.into_iter().map(Into::into).collect(),
            allele_dosages: value.allele_dosages,
            carriers: value.carriers.into_iter().map(Into::into).collect(),
        }
    }
}
impl From<WireGenotypeKind> for GenotypeCaseKind {
    fn from(value: WireGenotypeKind) -> Self {
        match value {
            WireGenotypeKind::Snp => Self::Snp,
            WireGenotypeKind::Indel => Self::Indel,
            WireGenotypeKind::Repeat => Self::Repeat,
        }
    }
}
impl From<WireLocus> for TruthLocus {
    fn from(value: WireLocus) -> Self {
        Self {
            sequence_id: value.sequence_id,
            start: value.start,
            end: value.end,
        }
    }
}
impl From<WireCarrier> for AlleleCarrier {
    fn from(value: WireCarrier) -> Self {
        Self {
            allele_index: value.allele_index,
            copy: value.copy.into(),
        }
    }
}
impl From<WireAllele> for TruthAllele {
    fn from(value: WireAllele) -> Self {
        match value {
            WireAllele::Reference {} => Self::Reference,
            WireAllele::Snp { base } => Self::Snp { base },
            WireAllele::Repeat {
                locus_id,
                bases,
                blocks,
            } => Self::Repeat {
                locus_id,
                bases,
                blocks: blocks.into_iter().map(Into::into).collect(),
            },
            WireAllele::Insertion {
                event_id,
                bases,
                source_event,
                piece_links,
            } => Self::Insertion {
                event_id,
                bases,
                source_event: source_event.into(),
                piece_links: piece_links.into_iter().map(Into::into).collect(),
            },
            WireAllele::Deletion {
                event_id,
                start,
                end,
                source_event,
                piece_links,
            } => Self::Deletion {
                event_id,
                start,
                end,
                source_event: source_event.into(),
                piece_links: piece_links.into_iter().map(Into::into).collect(),
            },
        }
    }
}
impl From<WireStructuralEvent> for TruthStructuralEvent {
    fn from(value: WireStructuralEvent) -> Self {
        match value {
            WireStructuralEvent::Insertion {
                event_id,
                point,
                bases,
            } => Self::Insertion {
                event_id,
                point,
                bases,
            },
            WireStructuralEvent::Deletion {
                event_id,
                start,
                end,
            } => Self::Deletion {
                event_id,
                start,
                end,
            },
        }
    }
}
impl From<WirePieceLink> for TruthPieceLink {
    fn from(value: WirePieceLink) -> Self {
        Self {
            piece_index: value.piece_index,
            start: value.start,
            end: value.end,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::TRUTH_LIST_FIELDS;
    use serde_json::Value;
    use std::collections::BTreeSet;

    use super::{
        CaseSetCheck, InspectedValidationError, MatrixCaseLoci, TruthValidationError,
        read_truth_jsonl, reference_content_sha256, validate_published_truth,
        validate_truth_for_vcf, validate_truth_with_copy_inspection,
    };
    use crate::identity::{ReferenceId, SequenceId};
    use crate::reference::{
        DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
    };
    use crate::truth_schema::TruthDocument;

    /// Referencia de ocho bases con el identificador dado, con las bases
    /// materializadas o solo su ficha.
    fn reference_named(id: &str, bases: &[u8], materialized: bool) -> Reference {
        let dna = if materialized {
            DnaSequence::try_from_bytes(bases).expect("DNA")
        } else {
            DnaSequence::summary(u64::try_from(bases.len()).expect("longitud"), 0)
                .expect("ficha sin bases")
        };
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "truth guards").expect("sequence ID"),
            FastaHeader::parse(b"chr1").expect("FASTA header"),
            SequenceClass::Chromosome,
            &dna,
        )
        .expect("reference sequence");
        Reference::new(
            ReferenceId::parse(id, "truth guards").expect("reference ID"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference")
    }

    /// La verdad mínima de `tests/fixtures/truth-v4`, con la huella de la
    /// referencia `REF_main` de ocho `A`.
    fn minimal_truth() -> (TruthDocument, Reference) {
        let reference = reference_named("REF_main", b"AAAAAAAA", true);
        let text = include_str!("../tests/fixtures/truth-v4/valid.jsonl").replacen(
            "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb",
            &reference_content_sha256(&reference).expect("huella"),
            1,
        );
        (
            read_truth_jsonl(&mut std::io::Cursor::new(text)).expect("verdad mínima"),
            reference,
        )
    }

    /// Una vista validada contra una referencia no se reutiliza contra otra:
    /// con otro identificador, o con el mismo identificador y una referencia
    /// cuyas bases no se conservaron (no se puede calcular su huella).
    ///
    /// Los consumidores de la vista (la proyección y el agrupador del VCF)
    /// reciben siempre la misma referencia con que se validó, de modo que
    /// ninguna ejecución llega aquí: se ejercita la guarda directamente.
    #[test]
    fn a_validated_view_is_not_reused_against_another_reference() {
        let (document, reference) = minimal_truth();
        let view = validate_truth_for_vcf(&document, &reference).expect("verdad válida");
        view.validate_reference(&reference)
            .expect("la misma referencia");
        assert_eq!(
            view.validate_reference(&reference_named("REF_other", b"AAAAAAAA", true))
                .err(),
            Some(TruthValidationError::ReferenceIdMismatch {
                declared: "REF_main".to_owned(),
                actual: "REF_other".to_owned(),
            })
        );
        assert_eq!(
            view.validate_reference(&reference_named("REF_main", b"AAAAAAAA", false))
                .err(),
            Some(TruthValidationError::BasesNotMaterialized)
        );
        assert_eq!(
            view.validate_reference(&reference_named("REF_main", b"AAAAAAAC", true))
                .err(),
            Some(TruthValidationError::ReferenceDigestMismatch)
        );
    }

    /// El juez completo de una verdad ya escrita, y el que inspecciona cada
    /// copia, rechazan una referencia de otro identificador. En el paquete
    /// publicado, la referencia sale de la configuración efectiva y la verdad
    /// de un fichero, de modo que pueden diferir.
    #[test]
    fn the_judges_reject_a_reference_of_another_identifier() {
        let (document, _) = minimal_truth();
        let other = reference_named("REF_other", b"AAAAAAAA", true);
        let expected = TruthValidationError::ReferenceIdMismatch {
            declared: "REF_main".to_owned(),
            actual: "REF_other".to_owned(),
        };
        assert_eq!(
            validate_published_truth(&document, &other, &MatrixCaseLoci::default()).err(),
            Some(expected.clone())
        );
        let inspected =
            validate_truth_with_copy_inspection(&document, &other, CaseSetCheck::Skip, |_, _| {
                Ok::<(), std::convert::Infallible>(())
            });
        assert!(
            matches!(
                inspected,
                Err(InspectedValidationError::Semantics(ref error)) if **error == expected
            ),
            "se esperaba la semántica con el identificador, y se obtuvo {inspected:?}"
        );
    }

    /// Si el valor de una propiedad es una lista, siguiendo `#/$defs/..`,
    /// `type` como cadena o como lista de tipos y las ramas de
    /// `oneOf`/`anyOf`/`allOf` (basta que una lo sea).
    fn is_array(node: &Value, root: &Value) -> bool {
        if let Some(reference) = node.get("$ref").and_then(Value::as_str) {
            let name = reference.strip_prefix("#/$defs/").expect("local reference");
            return is_array(&root["$defs"][name], root);
        }
        let by_type = match node.get("type") {
            Some(Value::String(kind)) => kind == "array",
            Some(Value::Array(kinds)) => kinds.iter().any(|kind| kind == "array"),
            _ => false,
        };
        by_type
            || ["oneOf", "anyOf", "allOf"].iter().any(|key| {
                node.get(key)
                    .and_then(Value::as_array)
                    .is_some_and(|branches| branches.iter().any(|branch| is_array(branch, root)))
            })
    }

    /// Nombres de propiedad del esquema, separados en los que son lista y los
    /// que no, en cualquier objeto `properties`.
    fn properties_by_kind(schema: &Value) -> (BTreeSet<String>, BTreeSet<String>) {
        fn walk(
            node: &Value,
            root: &Value,
            lists: &mut BTreeSet<String>,
            others: &mut BTreeSet<String>,
        ) {
            match node {
                Value::Object(map) => {
                    if let Some(Value::Object(properties)) = map.get("properties") {
                        for (name, property) in properties {
                            if is_array(property, root) {
                                lists.insert(name.clone());
                            } else {
                                others.insert(name.clone());
                            }
                        }
                    }
                    for value in map.values() {
                        walk(value, root, lists, others);
                    }
                }
                Value::Array(items) => items
                    .iter()
                    .for_each(|item| walk(item, root, lists, others)),
                _ => {}
            }
        }
        let (mut lists, mut others) = (BTreeSet::new(), BTreeSet::new());
        walk(schema, schema, &mut lists, &mut others);
        (lists, others)
    }

    /// Si el esquema añade una propiedad de lista que esta lista no conoce, el
    /// lector rechazaría una verdad válida; si esta lista conoce un nombre que
    /// el esquema no tiene como lista, el lector aceptaría una estructura
    /// escrita como lista. En los dos sentidos falla aquí. No ejercita nada del
    /// lector: solo fija la lista contra el esquema, de ahí que el barrido
    /// de `tests/truth_validation.rs` sea quien prueba el rechazo. Reconoce una
    /// lista por `type` (cadena o lista de tipos), `$ref` a `#/$defs/..` y
    /// ramas de `oneOf`/`anyOf`/`allOf`; otra forma de declararla no se vería.
    #[test]
    fn list_field_names_are_exactly_the_schema_array_properties() {
        let expected = TRUTH_LIST_FIELDS
            .iter()
            .map(|name| (*name).to_owned())
            .collect::<BTreeSet<_>>();
        let v5: Value = serde_json::from_str(include_str!("../schemas/truth-v5.schema.json"))
            .expect("schema v5 JSON");
        let (lists, others) = properties_by_kind(&v5);
        assert_eq!(lists, expected, "schema v5");
        // Un nombre de la lista que en otro sitio del esquema es un objeto o
        // un escalar dejaría pasar esa estructura escrita como lista.
        assert_eq!(
            others.intersection(&expected).collect::<Vec<_>>(),
            Vec::<&String>::new(),
            "a list name is also a non-list property"
        );
        // La v4 no tiene repeticiones: sus listas son las de la v5 menos `blocks`.
        let v4: Value = serde_json::from_str(include_str!("../schemas/truth-v4.schema.json"))
            .expect("schema v4 JSON");
        let (lists, others) = properties_by_kind(&v4);
        assert!(lists.is_subset(&expected), "schema v4");
        assert!(
            others.is_disjoint(&expected),
            "schema v4: a list name is not a list"
        );
    }
}
