//! Proyección pura de alelos pequeños normalizados para la vista VCF F2.5.2.
//!
//! No escribe VCF ni decide registros multialélicos, `GT` o procedencia INFO.
//! Esas responsabilidades empiezan en F2.5.3. Esta capa recibe la verdad
//! `truth.jsonl` ya estructurada, reconstruye los alelos finales desde sus
//! piezas vigentes y devuelve una entrada por copia y alelo alternativo.

use std::collections::BTreeMap;

use crate::reference::{Reference, ReferenceSequence};
use crate::truth_schema::{
    CopyIdentity, TruthCurrentPiece, TruthDocument, TruthRecord, TruthRepeatBlock,
};
use crate::truth_validation::{TruthValidationError, ValidatedTruthForVcf, validate_truth_for_vcf};

/// Versión estable del algoritmo interno de DTC-064.
pub const VCF_NORMALIZATION_VERSION: &str = "vcf-small-allele-normalization-v1";

/// Alelo alternativo de una copia, listo para que F2.5.3 lo agrupe.
///
/// `position` usa coordenadas internas zero-based. La escritura VCF posterior
/// la convertirá una sola vez a `POS` one-based.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NormalizedAllele {
    pub sequence_id: String,
    pub position: u64,
    pub reference: String,
    pub alternate: String,
    /// Intervalo de referencia semiaabierto que compuso este alelo final antes
    /// de recortar y alinear. F2.5.3 lo usa para enlazar todos los casos fuente
    /// que contribuyeron a un ALT compuesto, aunque su representación VCF se
    /// haya desplazado a la izquierda.
    pub source_start: u64,
    pub source_end: u64,
    pub copy: CopyIdentity,
    /// Secuencias repetidas que este alelo contiene, en el orden del tracto.
    ///
    /// Está vacío cuando el alelo no procede de ningún tracto declarado. Un
    /// alelo que sustituye un tracto lleva **los bloques tal y como la verdad
    /// los publica**: la estructura en bloques se declara y no puede derivarse
    /// de las bases (decisión S-30), de modo que aquí se transporta, no se
    /// recalcula.
    pub repeat_blocks: Vec<TruthRepeatBlock>,
}

/// Alelo normalizado de uso interno en la ruta de publicación.
///
/// La API pública conserva `NormalizedAllele` propietario. La publicación VCF
/// usa esta forma hasta que el agrupador conoce cuántos ALT distintos habrá:
/// una inserción directa grande puede entonces prestarse de `truth` para cada
/// copia portadora y materializarse una sola vez por valor de ALT final.
#[derive(Clone, Debug)]
pub(crate) struct ProjectedAllele<'truth> {
    pub(crate) sequence_id: &'truth str,
    pub(crate) position: u64,
    pub(crate) reference: String,
    alternate: ProjectedAlternate<'truth>,
    pub(crate) source_start: u64,
    pub(crate) source_end: u64,
    pub(crate) copy: &'truth CopyIdentity,
    /// Igual que en [`NormalizedAllele::repeat_blocks`], pero prestado.
    pub(crate) repeat_blocks: Vec<&'truth TruthRepeatBlock>,
}

#[derive(Clone, Debug)]
enum ProjectedAlternate<'truth> {
    Owned(String),
    /// La secuencia inicial era `anchor + bases`; `rotations` describe los
    /// desplazamientos a la izquierda como rotaciones a la derecha de esa
    /// secuencia. No se crea un buffer proporcional a `bases` por copia.
    DirectInsertion {
        /// Ancla alcanzada por la alineación izquierda; sustituye el primer
        /// byte de la ALT compacta exactamente como `left_align`.
        anchor: u8,
        /// Ancla original, necesaria para reconstruir el núcleo de la ALT sin
        /// materializar una cadena por copia.
        initial_anchor: u8,
        bases: &'truth str,
        /// Desplazamientos de ese núcleo, reducidos módulo `bases.len()`.
        rotations: usize,
        left_shifted: bool,
    },
}

/// Descriptor comparable por bytes sin concatenar los ALT de cada portador.
/// El sufijo corresponde a la expansión al REF común del registro.
#[derive(Clone, Copy, Debug)]
pub(crate) enum ProjectedAlternateValue<'truth> {
    Owned {
        alternate: &'truth str,
        suffix: &'truth str,
    },
    DirectInsertion {
        anchor: u8,
        initial_anchor: u8,
        bases: &'truth str,
        rotations: usize,
        left_shifted: bool,
        suffix: &'truth str,
    },
}

impl<'truth> ProjectedAllele<'truth> {
    pub(crate) fn alternate_with_suffix<'value>(
        &'value self,
        suffix: &'value str,
    ) -> ProjectedAlternateValue<'value> {
        match &self.alternate {
            ProjectedAlternate::Owned(alternate) => {
                ProjectedAlternateValue::Owned { alternate, suffix }
            }
            ProjectedAlternate::DirectInsertion {
                anchor,
                initial_anchor,
                bases,
                rotations,
                left_shifted,
            } => ProjectedAlternateValue::DirectInsertion {
                anchor: *anchor,
                initial_anchor: *initial_anchor,
                bases,
                rotations: *rotations,
                left_shifted: *left_shifted,
                suffix,
            },
        }
    }

    fn into_public(self) -> NormalizedAllele {
        let alternate = self.alternate_with_suffix("").materialize();
        NormalizedAllele {
            sequence_id: self.sequence_id.to_owned(),
            position: self.position,
            reference: self.reference,
            alternate,
            source_start: self.source_start,
            source_end: self.source_end,
            copy: self.copy.to_owned(),
            repeat_blocks: self.repeat_blocks.into_iter().cloned().collect(),
        }
    }
}

impl ProjectedAlternateValue<'_> {
    pub(crate) fn materialize(self) -> String {
        let mut result = String::with_capacity(self.len());
        for index in 0..self.len() {
            result.push(char::from(self.byte_at(index)));
        }
        result
    }

    pub(crate) fn len(self) -> usize {
        match self {
            Self::Owned { alternate, suffix } => alternate.len() + suffix.len(),
            Self::DirectInsertion { bases, suffix, .. } => bases.len() + 1 + suffix.len(),
        }
    }

    fn byte_at(self, index: usize) -> u8 {
        match self {
            Self::Owned { alternate, suffix } => alternate
                .as_bytes()
                .get(index)
                .copied()
                .unwrap_or_else(|| suffix.as_bytes()[index - alternate.len()]),
            Self::DirectInsertion {
                anchor,
                initial_anchor,
                bases,
                rotations,
                left_shifted,
                suffix,
            } => {
                let inserted_length = bases.len() + 1;
                if index >= inserted_length {
                    return suffix.as_bytes()[index - inserted_length];
                }
                if index == 0 {
                    return anchor;
                }
                if !left_shifted {
                    return bases.as_bytes()[index - 1];
                }
                let core_length = bases.len();
                let core_index = (core_length - rotations + index) % core_length;
                if core_index == 0 {
                    initial_anchor
                } else {
                    bases.as_bytes()[core_index - 1]
                }
            }
        }
    }
}

impl PartialEq for ProjectedAlternateValue<'_> {
    fn eq(&self, other: &Self) -> bool {
        self.len() == other.len()
            && (0..self.len()).all(|index| self.byte_at(index) == other.byte_at(index))
    }
}

impl Eq for ProjectedAlternateValue<'_> {}

impl PartialOrd for ProjectedAlternateValue<'_> {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Ord for ProjectedAlternateValue<'_> {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let common = self.len().min(other.len());
        for index in 0..common {
            match self.byte_at(index).cmp(&other.byte_at(index)) {
                std::cmp::Ordering::Equal => {}
                order => return order,
            }
        }
        self.len().cmp(&other.len())
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfProjectionError {
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
    Truth(TruthValidationError),
    MissingReferenceSequence {
        sequence_id: String,
    },
    SourceOutsideReference {
        copy: CopyIdentity,
    },
    NonCollinearSourcePieces {
        copy: CopyIdentity,
    },
    CannotAnchor {
        sequence_id: String,
        position: u64,
    },
    CoordinateOverflow,
}

impl std::fmt::Display for VcfProjectionError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for VcfProjectionError {}

/// Proyecta diferencias pequeñas de la verdad válida contra la referencia.
///
/// La validación semántica se repite al inicio para que esta frontera no pueda
/// aceptar una estructura Rust que parezca un `truth.jsonl` válido pero no
/// reconstruya las copias. No hay E/S, herramientas externas ni reparación.
/// Primero se forma el intervalo de referencia común de todo cambio solapado;
/// después se recorta y alinea cada alelo individualmente.
pub fn project_normalized_alleles(
    truth: &TruthDocument,
    reference: &Reference,
) -> Result<Vec<NormalizedAllele>, VcfProjectionError> {
    let validated = validate_truth_for_vcf(truth, reference).map_err(VcfProjectionError::Truth)?;
    project_projected_alleles_from_validated_truth(validated, reference).map(|alleles| {
        alleles
            .into_iter()
            .map(ProjectedAllele::into_public)
            .collect()
    })
}

/// Proyecta una verdad cuya estructura y semántica ya se validaron completa y
/// secuencialmente. Solo la publicación VCF puede encadenar esta variante con
/// el agrupador para no reconstruir todas las copias varias veces.
pub(crate) fn project_projected_alleles_from_validated_truth<'truth>(
    truth: ValidatedTruthForVcf<'truth>,
    reference: &Reference,
) -> Result<Vec<ProjectedAllele<'truth>>, VcfProjectionError> {
    truth
        .validate_reference(reference)
        .map_err(VcfProjectionError::Truth)?;
    let truth = truth.document();

    let catalog = reference
        .sequences()
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.id().as_str(), (index, sequence)))
        .collect::<BTreeMap<_, _>>();
    let mut paths = Vec::new();
    for record in &truth.records {
        let TruthRecord::HaplotypeCopy(copy) = record else {
            continue;
        };
        let (_, sequence) = catalog
            .get(copy.identity.sequence_id.as_str())
            .ok_or_else(|| VcfProjectionError::MissingReferenceSequence {
                sequence_id: copy.identity.sequence_id.clone(),
            })?;
        paths.push(CopyPath::from_truth(
            &copy.identity,
            &copy.current_pieces,
            sequence,
        )?);
    }

    let mut changes = BTreeMap::<String, Vec<ReferenceSpan>>::new();
    for path in &paths {
        changes
            .entry(path.sequence_id.to_owned())
            .or_default()
            .extend(path.changes());
    }

    let mut projected = Vec::new();
    for (sequence_id, spans) in changes {
        let (_, sequence) = catalog.get(sequence_id.as_str()).ok_or_else(|| {
            VcfProjectionError::MissingReferenceSequence {
                sequence_id: sequence_id.clone(),
            }
        })?;
        let components = merge_overlapping(spans);
        for path in paths.iter().filter(|path| path.sequence_id == sequence_id) {
            projected.extend(normalize_path_components(path, sequence, &components)?);
        }
    }

    projected.sort_by(|left, right| {
        let left_catalog = catalog
            .get(left.sequence_id)
            .expect("projected sequence belongs to the reference")
            .0;
        let right_catalog = catalog
            .get(right.sequence_id)
            .expect("projected sequence belongs to the reference")
            .0;
        (left_catalog, left.position, &left.reference)
            .cmp(&(right_catalog, right.position, &right.reference))
            .then_with(|| {
                left.alternate_with_suffix("")
                    .cmp(&right.alternate_with_suffix(""))
            })
            .then_with(|| left.copy.cmp(right.copy))
    });
    Ok(projected)
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
struct ReferenceSpan {
    start: u64,
    end: u64,
}

#[derive(Clone, Debug)]
struct CopyPath<'truth> {
    identity: &'truth CopyIdentity,
    sequence_id: &'truth str,
    pieces: Vec<PathPiece<'truth>>,
}

#[derive(Clone, Debug)]
enum PathPiece<'truth> {
    Source {
        start: u64,
        end: u64,
        snps: &'truth [crate::truth_schema::TruthCurrentSnp],
    },
    Inserted {
        anchor: u64,
        bases: &'truth str,
    },
    /// Tracto de repetición realizado: sustituye el intervalo
    /// `[start, end)` de la referencia por sus propias bases.
    Replaced {
        start: u64,
        end: u64,
        bases: &'truth str,
        /// Estructura en bloques declarada del tracto, que es lo que el VCF
        /// publica como secuencias repetidas del alelo.
        blocks: &'truth [TruthRepeatBlock],
    },
}

impl<'truth> CopyPath<'truth> {
    fn from_truth(
        identity: &'truth CopyIdentity,
        current_pieces: &'truth [TruthCurrentPiece],
        sequence: &ReferenceSequence,
    ) -> Result<Self, VcfProjectionError> {
        let mut pieces = Vec::with_capacity(current_pieces.len());
        let mut prior_source_end = None;
        for current in current_pieces {
            match current {
                TruthCurrentPiece::Source { source, snps } => {
                    if source.sequence_id != sequence.id().as_str()
                        || source.start >= source.end
                        || source.end > sequence.length()
                        || prior_source_end.is_some_and(|end| source.start < end)
                    {
                        return Err(VcfProjectionError::SourceOutsideReference {
                            copy: identity.clone(),
                        });
                    }
                    for snp in snps {
                        if snp.position < source.start || snp.position >= source.end {
                            return Err(VcfProjectionError::SourceOutsideReference {
                                copy: identity.clone(),
                            });
                        }
                    }
                    pieces.push(PathPiece::Source {
                        start: source.start,
                        end: source.end,
                        snps,
                    });
                    prior_source_end = Some(source.end);
                }
                TruthCurrentPiece::Inserted { bases, .. } => {
                    let anchor = prior_source_end.ok_or_else(|| {
                        VcfProjectionError::NonCollinearSourcePieces {
                            copy: identity.clone(),
                        }
                    })?;
                    pieces.push(PathPiece::Inserted { anchor, bases });
                }
                TruthCurrentPiece::RepeatAllele {
                    source,
                    bases,
                    blocks,
                    ..
                } => {
                    if source.sequence_id != sequence.id().as_str()
                        || source.start >= source.end
                        || source.end > sequence.length()
                        || prior_source_end.is_some_and(|end| source.start < end)
                    {
                        return Err(VcfProjectionError::SourceOutsideReference {
                            copy: identity.clone(),
                        });
                    }
                    pieces.push(PathPiece::Replaced {
                        start: source.start,
                        end: source.end,
                        bases,
                        blocks,
                    });
                    prior_source_end = Some(source.end);
                }
            }
        }
        Ok(Self {
            sequence_id: &identity.sequence_id,
            identity,
            pieces,
        })
    }

    fn changes(&self) -> Vec<ReferenceSpan> {
        let mut changes = Vec::new();
        let mut prior_source_end = None;
        for piece in &self.pieces {
            match piece {
                PathPiece::Source {
                    start, end, snps, ..
                } => {
                    if let Some(previous) = prior_source_end
                        && previous < *start
                    {
                        changes.push(ReferenceSpan {
                            start: previous,
                            end: *start,
                        });
                    }
                    for position in snps.iter().map(|snp| snp.position) {
                        changes.push(ReferenceSpan {
                            start: position,
                            end: position.saturating_add(1),
                        });
                    }
                    prior_source_end = Some(*end);
                }
                PathPiece::Inserted { anchor, .. } => changes.push(ReferenceSpan {
                    start: *anchor,
                    end: *anchor,
                }),
                PathPiece::Replaced { start, end, .. } => {
                    if let Some(previous) = prior_source_end
                        && previous < *start
                    {
                        changes.push(ReferenceSpan {
                            start: previous,
                            end: *start,
                        });
                    }
                    changes.push(ReferenceSpan {
                        start: *start,
                        end: *end,
                    });
                    prior_source_end = Some(*end);
                }
            }
        }
        changes
    }

    fn materialize_window(
        &self,
        sequence: &ReferenceSequence,
        start: u64,
        end: u64,
    ) -> Result<Vec<u8>, VcfProjectionError> {
        let mut output = Vec::new();
        for piece in &self.pieces {
            match piece {
                PathPiece::Source {
                    start: piece_start,
                    end: piece_end,
                    snps,
                    ..
                } if *piece_start < end && start < *piece_end => {
                    let window_start = (*piece_start).max(start);
                    let window_end = (*piece_end).min(end);
                    let local_start = usize::try_from(window_start)
                        .map_err(|_| VcfProjectionError::CoordinateOverflow)?;
                    let local_end = usize::try_from(window_end)
                        .map_err(|_| VcfProjectionError::CoordinateOverflow)?;
                    let offset = output.len();
                    output.extend_from_slice(
                        &sequence
                            .sequence()
                            .as_bytes()
                            .map_err(|_| VcfProjectionError::BasesNotMaterialized)?
                            [local_start..local_end],
                    );
                    for snp in *snps {
                        if window_start <= snp.position && snp.position < window_end {
                            let relative = usize::try_from(snp.position - window_start)
                                .map_err(|_| VcfProjectionError::CoordinateOverflow)?;
                            output[offset + relative] = snp.alternate.as_bytes()[0];
                        }
                    }
                }
                PathPiece::Inserted { anchor, bases } if start <= *anchor && *anchor <= end => {
                    output.extend_from_slice(bases.as_bytes());
                }
                // El tracto se copia entero o no se copia: la ventana procede
                // de fundir los intervalos que cambiaron, de modo que o
                // contiene el tracto completo o no lo toca. Una ventana que lo
                // cortara por la mitad produciría un ALT recortado y silencioso.
                PathPiece::Replaced {
                    start: piece_start,
                    end: piece_end,
                    bases,
                    ..
                } if *piece_start < end && start < *piece_end => {
                    if *piece_start < start || end < *piece_end {
                        return Err(VcfProjectionError::CoordinateOverflow);
                    }
                    output.extend_from_slice(bases.as_bytes());
                }
                _ => {}
            }
        }
        Ok(output)
    }

    /// Bloques de los tractos que caen dentro de la ventana `[start, end)`.
    ///
    /// El orden es el de las piezas de la copia, que es el de la secuencia, de
    /// modo que las listas paralelas que el VCF publica quedan en el mismo
    /// orden que las bases del alelo. La condición de intersección es la misma
    /// que usa [`Self::materialize_window`], que además rechaza una ventana que
    /// partiera un tracto por la mitad: aquí no puede colarse media anotación.
    fn repeat_blocks_in(&self, start: u64, end: u64) -> Vec<&'truth TruthRepeatBlock> {
        self.pieces
            .iter()
            .filter_map(|piece| match piece {
                PathPiece::Replaced {
                    start: piece_start,
                    end: piece_end,
                    blocks,
                    ..
                } if *piece_start < end && start < *piece_end => Some(*blocks),
                _ => None,
            })
            .flatten()
            .collect()
    }

    fn direct_insertion_at(&self, anchor: u64) -> Option<&'truth str> {
        let mut inserted = self.pieces.iter().filter_map(|piece| match piece {
            PathPiece::Inserted {
                anchor: piece_anchor,
                bases,
            } if *piece_anchor == anchor => Some(*bases),
            _ => None,
        });
        let bases = inserted.next()?;
        inserted.next().is_none().then_some(bases)
    }
}

fn merge_overlapping(mut spans: Vec<ReferenceSpan>) -> Vec<ReferenceSpan> {
    spans.sort_by_key(|span| (span.start, span.end));
    let mut merged = Vec::<ReferenceSpan>::new();
    let mut right_boundary_bridge = None;
    for span in spans {
        if let Some(current) = merged.last_mut() {
            let previous_end = current.end;
            let joins_through_insertion =
                right_boundary_bridge == Some(previous_end) && span.start == previous_end;
            if spans_overlap(*current, span) || joins_through_insertion {
                current.end = current.end.max(span.end);
                right_boundary_bridge = if span.start == span.end && span.start == previous_end {
                    Some(previous_end)
                } else if current.end > previous_end {
                    None
                } else {
                    right_boundary_bridge
                };
            } else {
                merged.push(span);
                right_boundary_bridge = None;
            }
        } else {
            merged.push(span);
        }
    }
    merged
}

/// Los intervalos de referencia son semiabiertos: dos cambios de una base que
/// sólo son adyacentes no son un componente. Una inserción es una frontera de
/// longitud cero; `merge_overlapping` conserva además ese punto como puente
/// entre los dos cambios que toca, para materializarlo una sola vez.
fn spans_overlap(left: ReferenceSpan, right: ReferenceSpan) -> bool {
    match (left.start == left.end, right.start == right.end) {
        (false, false) => left.start < right.end && right.start < left.end,
        (true, true) => left.start == right.start,
        (true, false) => right.start <= left.start && left.start <= right.end,
        (false, true) => left.start <= right.start && right.start <= left.end,
    }
}

/// Recompone hasta punto fijo los componentes de una copia que, después de la
/// alineación izquierda, ocupan el mismo `CHROM+POS`. La unión vuelve a leer la
/// secuencia haplotípica final, por lo que conserva inserciones intermedias y
/// nunca representa dos veces una misma copia en ese locus.
fn normalize_path_components<'truth>(
    path: &CopyPath<'truth>,
    sequence: &ReferenceSequence,
    initial: &[ReferenceSpan],
) -> Result<Vec<ProjectedAllele<'truth>>, VcfProjectionError> {
    let mut components = initial.to_vec();
    loop {
        let normalized = components
            .iter()
            .copied()
            .map(|span| normalize_path_window(path, sequence, span).map(|allele| (span, allele)))
            .collect::<Result<Vec<_>, _>>()?;
        let mut collisions = normalized_ref_collisions(&normalized)?;
        collisions.extend(normalization_traversal_closures(&normalized, &components)?);
        if collisions.is_empty() {
            return Ok(normalized
                .into_iter()
                .filter_map(|(span, allele)| {
                    allele.map(|mut allele| {
                        allele.source_start = span.start;
                        allele.source_end = span.end;
                        allele
                    })
                })
                .collect());
        }
        components = merge_normalized_collisions(components, &collisions);
    }
}

/// Devuelve los cierres que un alelo debe absorber porque quedaron atravesados
/// al desplazarlo a la izquierda. El intervalo conserva tanto el span fuente
/// como el REF normalizado, de modo que un SNP intermedio no puede sobrevivir
/// como un segundo registro incompatible con el alelo recompuesto.
fn normalization_traversal_closures(
    normalized: &[(ReferenceSpan, Option<ProjectedAllele<'_>>)],
    components: &[ReferenceSpan],
) -> Result<Vec<Vec<ReferenceSpan>>, VcfProjectionError> {
    let mut closures = Vec::new();
    for (span, allele) in normalized {
        let Some(allele) = allele else {
            continue;
        };
        let length = u64::try_from(allele.reference.len())
            .map_err(|_| VcfProjectionError::CoordinateOverflow)?;
        let normalized_end = allele
            .position
            .checked_add(length)
            .ok_or(VcfProjectionError::CoordinateOverflow)?;
        let traversed = ReferenceSpan {
            start: span.start.min(allele.position),
            end: span.end.max(normalized_end),
        };
        let contained = components
            .iter()
            .copied()
            .filter(|component| span_is_within(*component, traversed))
            .collect::<Vec<_>>();
        if contained.len() > 1 {
            closures.push(contained);
        }
    }
    Ok(closures)
}

/// Detecta componentes cuyos intervalos REF VCF ya normalizados se solapan de
/// verdad. La condición es semiabierta: el mero contacto de extremos no basta.
fn normalized_ref_collisions(
    normalized: &[(ReferenceSpan, Option<ProjectedAllele<'_>>)],
) -> Result<Vec<Vec<ReferenceSpan>>, VcfProjectionError> {
    let mut intervals = normalized
        .iter()
        .filter_map(|(span, allele)| {
            allele.as_ref().map(|allele| {
                let length = u64::try_from(allele.reference.len())
                    .map_err(|_| VcfProjectionError::CoordinateOverflow)?;
                let end = allele
                    .position
                    .checked_add(length)
                    .ok_or(VcfProjectionError::CoordinateOverflow)?;
                Ok((allele.position, end, *span))
            })
        })
        .collect::<Result<Vec<_>, VcfProjectionError>>()?;
    intervals.sort_by_key(|(start, end, span)| (*start, *end, span.start, span.end));
    let mut collisions = Vec::new();
    let mut current = Vec::<ReferenceSpan>::new();
    let mut current_end = 0_u64;
    for (start, end, span) in intervals {
        if current.is_empty() {
            current.push(span);
            current_end = end;
        } else if start < current_end {
            current.push(span);
            current_end = current_end.max(end);
        } else {
            if current.len() > 1 {
                collisions.push(current);
            }
            current = vec![span];
            current_end = end;
        }
    }
    if current.len() > 1 {
        collisions.push(current);
    }
    Ok(collisions)
}

fn merge_normalized_collisions(
    components: Vec<ReferenceSpan>,
    collisions: &[Vec<ReferenceSpan>],
) -> Vec<ReferenceSpan> {
    let indices = components
        .iter()
        .enumerate()
        .map(|(index, span)| ((span.start, span.end), index))
        .collect::<BTreeMap<_, _>>();
    let mut disjoint = DisjointComponents::new(components.len());
    for collision in collisions {
        let Some((first, remainder)) = collision.split_first() else {
            continue;
        };
        let Some(&first_index) = indices.get(&(first.start, first.end)) else {
            continue;
        };
        for span in remainder {
            if let Some(&index) = indices.get(&(span.start, span.end)) {
                disjoint.union(first_index, index);
            }
        }
    }
    // El intervalo recién recompuesto puede envolver un tercer cambio que no
    // colisionaba con ninguno de los dos VCF minimizados iniciales. Se absorbe
    // antes de reproyectar, y se repite hasta cerrar esa relación.
    loop {
        let envelopes = disjoint.envelopes(&components);
        let mut changed = false;
        for (root, envelope, member_count) in &envelopes {
            if *member_count < 2 {
                continue;
            }
            for (index, span) in components.iter().copied().enumerate() {
                if disjoint.find(index) != *root && span_is_within(span, *envelope) {
                    disjoint.union(*root, index);
                    changed = true;
                }
            }
        }
        if !changed {
            break;
        }
    }
    let envelopes = disjoint.envelopes(&components);
    let replacement = envelopes
        .into_iter()
        .map(|(root, envelope, _)| (root, envelope))
        .collect::<BTreeMap<_, _>>();
    let mut merged = components
        .iter()
        .enumerate()
        .map(|(index, span)| {
            replacement
                .get(&disjoint.find(index))
                .copied()
                .unwrap_or(*span)
        })
        .collect::<Vec<_>>();
    merged.sort_by_key(|span| (span.start, span.end));
    merged.dedup();
    merged
}

fn span_is_within(candidate: ReferenceSpan, envelope: ReferenceSpan) -> bool {
    envelope.start <= candidate.start && candidate.end <= envelope.end
}

struct DisjointComponents {
    parents: Vec<usize>,
}

impl DisjointComponents {
    fn new(count: usize) -> Self {
        Self {
            parents: (0..count).collect(),
        }
    }

    fn find(&mut self, index: usize) -> usize {
        let parent = self.parents[index];
        if parent != index {
            let root = self.find(parent);
            self.parents[index] = root;
        }
        self.parents[index]
    }

    fn union(&mut self, left: usize, right: usize) {
        let left = self.find(left);
        let right = self.find(right);
        if left != right {
            self.parents[right] = left;
        }
    }

    fn envelopes(&mut self, components: &[ReferenceSpan]) -> Vec<(usize, ReferenceSpan, usize)> {
        let mut grouped = BTreeMap::<usize, (ReferenceSpan, usize)>::new();
        for (index, span) in components.iter().copied().enumerate() {
            let root = self.find(index);
            grouped
                .entry(root)
                .and_modify(|(envelope, count)| {
                    envelope.start = envelope.start.min(span.start);
                    envelope.end = envelope.end.max(span.end);
                    *count += 1;
                })
                .or_insert((span, 1));
        }
        grouped
            .into_iter()
            .map(|(root, (envelope, count))| (root, envelope, count))
            .collect()
    }
}

fn normalize_path_window<'truth>(
    path: &CopyPath<'truth>,
    sequence: &ReferenceSequence,
    span: ReferenceSpan,
) -> Result<Option<ProjectedAllele<'truth>>, VcfProjectionError> {
    if span.start == span.end
        && let Some(bases) = path.direct_insertion_at(span.start)
    {
        return normalize_direct_insertion(path, sequence, span, bases).map(Some);
    }
    let reference_bases = sequence
        .sequence()
        .as_bytes()
        .map_err(|_| VcfProjectionError::BasesNotMaterialized)?;
    let mut start = span.start;
    let end = span.end;
    let mut reference = reference_bases[usize::try_from(start)
        .map_err(|_| VcfProjectionError::CoordinateOverflow)?
        ..usize::try_from(end).map_err(|_| VcfProjectionError::CoordinateOverflow)?]
        .to_vec();
    let mut alternate = path.materialize_window(sequence, start, end)?;
    if reference == alternate {
        return Ok(None);
    }
    if reference.is_empty() || alternate.is_empty() {
        start = start
            .checked_sub(1)
            .ok_or_else(|| VcfProjectionError::CannotAnchor {
                sequence_id: path.sequence_id.to_owned(),
                position: span.start,
            })?;
        reference = reference_bases[usize::try_from(start)
            .map_err(|_| VcfProjectionError::CoordinateOverflow)?
            ..usize::try_from(end).map_err(|_| VcfProjectionError::CoordinateOverflow)?]
            .to_vec();
        alternate = path.materialize_window(sequence, start, end)?;
    }

    // **Los registros de repetición no se recortan.**
    //
    // DTC-064 manda recortar el prefijo y el sufijo comunes de todo alelo
    // pequeño, y la decisión S-09 manda publicar los de repetición «con la
    // secuencia completa en `REF` y `ALT`». Las dos no pueden cumplirse a la
    // vez sobre el mismo registro, y el responsable resolvió a favor de S-09
    // el 23 de septiembre de 2026: **para los registros de repetición manda
    // S-09; para todos los demás, DTC-064 sigue intacta**.
    //
    // No es un capricho de formato: los campos de repetición que el registro
    // publica dicen cuántas bases tiene el tracto, y recortar dejaba un `ALT`
    // de una base anotado como si contuviera noventa. Con el tracto entero la
    // anotación es verdadera, y es además una de las dos formas que el propio
    // estándar ofrece —la de «reportar la secuencia completa de todos los
    // alelos del locus en un único registro»—.
    let repeat_blocks = path.repeat_blocks_in(span.start, span.end);
    if repeat_blocks.is_empty() {
        trim_common(&mut start, &mut reference, &mut alternate);
        if reference.len() != alternate.len() {
            left_align(reference_bases, &mut start, &mut reference, &mut alternate);
        }
    }
    Ok(Some(ProjectedAllele {
        sequence_id: path.sequence_id,
        position: start,
        reference: String::from_utf8(reference).expect("reference alphabet is UTF-8 DNA"),
        alternate: ProjectedAlternate::Owned(
            String::from_utf8(alternate).expect("truth alphabet is UTF-8 DNA"),
        ),
        // Se reemplazan en `normalize_path_components` por el componente final
        // después de cerrar colisiones de la normalización.
        source_start: span.start,
        source_end: span.end,
        copy: path.identity,
        repeat_blocks,
    }))
}

fn normalize_direct_insertion<'truth>(
    path: &CopyPath<'truth>,
    sequence: &ReferenceSequence,
    span: ReferenceSpan,
    bases: &'truth str,
) -> Result<ProjectedAllele<'truth>, VcfProjectionError> {
    let mut position =
        span.start
            .checked_sub(1)
            .ok_or_else(|| VcfProjectionError::CannotAnchor {
                sequence_id: path.sequence_id.to_owned(),
                position: span.start,
            })?;
    let reference_bases = sequence
        .sequence()
        .as_bytes()
        .map_err(|_| VcfProjectionError::BasesNotMaterialized)?;
    let mut anchor = reference_bases
        [usize::try_from(position).map_err(|_| VcfProjectionError::CoordinateOverflow)?];
    let initial_anchor = anchor;
    let core_length = bases.len();
    debug_assert!(core_length > 0, "validated direct insertion has bases");
    let mut rotations = 0_usize;
    let mut left_shifted = false;
    let mut last = *bases
        .as_bytes()
        .last()
        .expect("validated direct insertion has bases");
    while position > 0 && anchor == last {
        position -= 1;
        anchor = reference_bases
            [usize::try_from(position).map_err(|_| VcfProjectionError::CoordinateOverflow)?];
        rotations = (rotations + 1) % core_length;
        let core_last_index = (core_length - rotations) % core_length;
        last = if core_last_index == 0 {
            initial_anchor
        } else {
            bases.as_bytes()[core_last_index - 1]
        };
        left_shifted = true;
    }
    Ok(ProjectedAllele {
        sequence_id: path.sequence_id,
        position,
        reference: char::from(anchor).to_string(),
        alternate: ProjectedAlternate::DirectInsertion {
            anchor,
            initial_anchor,
            bases,
            rotations,
            left_shifted,
        },
        source_start: span.start,
        source_end: span.end,
        copy: path.identity,
        // Una inserción directa es una frontera de anchura cero: no puede
        // contener un tracto, que siempre sustituye un intervalo no vacío.
        repeat_blocks: Vec::new(),
    })
}

fn trim_common(position: &mut u64, reference: &mut Vec<u8>, alternate: &mut Vec<u8>) {
    while reference.len() > 1 && alternate.len() > 1 && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }
    while reference.len() > 1 && alternate.len() > 1 && reference.first() == alternate.first() {
        reference.remove(0);
        alternate.remove(0);
        *position += 1;
    }
}

fn left_align(
    reference_sequence: &[u8],
    position: &mut u64,
    reference: &mut Vec<u8>,
    alternate: &mut Vec<u8>,
) {
    while *position > 0 && reference.last() == alternate.last() {
        *position -= 1;
        let base = reference_sequence[usize::try_from(*position).expect("reference fits usize")];
        reference.insert(0, base);
        alternate.insert(0, base);
        reference.pop();
        alternate.pop();
    }
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::process::Command;
    use std::sync::atomic::{AtomicU64, Ordering};

    use super::*;
    use crate::identity::{ReferenceId, SequenceId};
    use crate::reference::{DnaSequence, FastaHeader, ReferenceOrigin, SequenceClass};
    use crate::truth_schema::{
        CopyProvenance, HaplotypeCopyRecord, NamedDigest, SegmentSource,
        TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TRUTH_SCHEMA_VERSION, TruthCurrentSnp,
        TruthHeader, TruthModification, TruthSourceKind, TruthSourceRef,
    };
    use crate::truth_validation::reference_content_sha256;

    static NONCE: AtomicU64 = AtomicU64::new(0);

    fn reference(bases: &[u8]) -> Reference {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "VCF projection test").expect("sequence"),
            FastaHeader::parse(b"chr1").expect("header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(bases).expect("DNA"),
        )
        .expect("sequence");
        Reference::new(
            ReferenceId::parse("REF_main", "VCF projection test").expect("reference"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference")
    }

    fn source(start: u64, end: u64, length: u64) -> SegmentSource {
        SegmentSource {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            reference_id: "REF_main".into(),
            sequence_id: "chr1".into(),
            start,
            end,
            sequence_length: length,
        }
    }

    fn copy(
        index: u8,
        pieces: Vec<TruthCurrentPiece>,
        modifications: Vec<TruthModification>,
    ) -> HaplotypeCopyRecord {
        copy_with_length(index, 5, pieces, modifications)
    }

    fn copy_with_length(
        index: u8,
        length: u64,
        pieces: Vec<TruthCurrentPiece>,
        modifications: Vec<TruthModification>,
    ) -> HaplotypeCopyRecord {
        HaplotypeCopyRecord {
            identity: CopyIdentity {
                reference_id: "REF_main".into(),
                individual_id: format!("IND_{index:02}"),
                sequence_id: "chr1".into(),
                copy_index: index,
            },
            phase: index,
            provenance: CopyProvenance {
                origin: "test".into(),
                sources: vec![TruthSourceRef {
                    kind: TruthSourceKind::Reference,
                    source_id: "REF_main".into(),
                }],
            },
            segments: vec![crate::truth_schema::TruthSegment {
                source: source(0, length, length),
                modifications,
            }],
            current_pieces: pieces,
            history: vec![],
        }
    }

    fn document(reference: &Reference, copies: Vec<HaplotypeCopyRecord>) -> TruthDocument {
        TruthDocument {
            header: TruthHeader {
                document_kind: TRUTH_DOCUMENT_KIND.into(),
                schema_version: TRUTH_SCHEMA_VERSION.into(),
                execution_id: "RUN_01".into(),
                reference_id: "REF_main".into(),
                reference_content_sha256: reference_content_sha256(reference).expect("huella"),
                coordinate_convention: TRUTH_COORDINATE_CONVENTION.into(),
                copy_count: u64::try_from(copies.len()).expect("count"),
                genotype_count: 0,
                source_digests: vec![NamedDigest {
                    source: TruthSourceRef {
                        kind: TruthSourceKind::Reference,
                        source_id: "REF_main".into(),
                    },
                    sha256: "a".repeat(64),
                }],
                group_digests: vec![],
            },
            records: copies.into_iter().map(TruthRecord::HaplotypeCopy).collect(),
        }
    }

    fn direct_insertion_document(reference: &Reference, point: u64, bases: &str) -> TruthDocument {
        let length = reference.sequences()[0].sequence().len();
        let insertion = TruthModification::Insertion {
            event_id: "INS_01".into(),
            point,
            bases: bases.into(),
            provenance: explicit(),
        };
        let insertion_copy = copy_with_length(
            0,
            length,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, point, length),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: bases.into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(point, length, length),
                    snps: vec![],
                },
            ],
            vec![insertion],
        );
        document(reference, vec![insertion_copy])
    }

    fn explicit() -> crate::truth_schema::ModificationProvenance {
        crate::truth_schema::ModificationProvenance {
            kind: crate::truth_schema::ModificationProvenanceKind::Explicit,
            version: "1".into(),
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".into(),
            },
            matrix: None,
        }
    }

    fn apply_normalized(reference: &[u8], allele: &NormalizedAllele) -> Vec<u8> {
        let start = usize::try_from(allele.position).expect("test position");
        let end = start + allele.reference.len();
        assert_eq!(&reference[start..end], allele.reference.as_bytes());
        let mut reconstructed = reference[..start].to_vec();
        reconstructed.extend_from_slice(allele.alternate.as_bytes());
        reconstructed.extend_from_slice(&reference[end..]);
        reconstructed
    }

    #[test]
    fn projects_a_snp_from_final_truth_pieces() {
        let reference = reference(b"AACGT");
        let alteration = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 2,
            reference: "C".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let mut altered = copy(
            0,
            vec![TruthCurrentPiece::Source {
                source: source(0, 5, 5),
                snps: vec![TruthCurrentSnp {
                    event_id: "SNP_01".into(),
                    position: 2,
                    reference: "C".into(),
                    alternate: "T".into(),
                    provenance: explicit(),
                }],
            }],
            vec![alteration],
        );
        altered.history = vec![crate::truth_schema::TruthHistoryEvent::Snp {
            event_id: "SNP_01".into(),
            position: 2,
            provenance: explicit(),
        }];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (2, "C", "T")
        );
    }

    #[test]
    fn source_paths_retain_sparse_metadata_and_materialize_only_the_requested_window() {
        let reference = reference(&vec![b'A'; 1_000_000]);
        let sequence = &reference.sequences()[0];
        let mut copy = copy_with_length(0, 1_000_000, Vec::new(), Vec::new());
        copy.current_pieces = vec![TruthCurrentPiece::Source {
            source: source(0, 1_000_000, 1_000_000),
            snps: vec![TruthCurrentSnp {
                event_id: "SNP_WINDOW".to_owned(),
                position: 500_000,
                reference: "A".to_owned(),
                alternate: "G".to_owned(),
                provenance: explicit(),
            }],
        }];
        let path = CopyPath::from_truth(&copy.identity, &copy.current_pieces, sequence)
            .expect("compact source path");

        assert!(matches!(
            path.pieces.as_slice(),
            [PathPiece::Source { start: 0, end: 1_000_000, snps }] if snps.len() == 1
        ));
        assert_eq!(
            path.materialize_window(sequence, 499_999, 500_002)
                .expect("window"),
            b"AGA"
        );
    }

    #[test]
    fn validated_compact_view_rejects_a_mixed_reference_at_its_effective_boundary() {
        let validated_reference = reference(b"AAAAA");
        let truth = document(
            &validated_reference,
            vec![copy_with_length(
                0,
                5,
                vec![TruthCurrentPiece::Source {
                    source: source(0, 5, 5),
                    snps: Vec::new(),
                }],
                Vec::new(),
            )],
        );
        let validated = validate_truth_for_vcf(&truth, &validated_reference)
            .expect("truth validated for its original reference");
        let mixed_reference = reference(b"AAAAT");

        assert!(matches!(
            project_projected_alleles_from_validated_truth(validated, &mixed_reference),
            Err(VcfProjectionError::Truth(
                TruthValidationError::ReferenceDigestMismatch
            ))
        ));
    }

    #[test]
    fn large_direct_insertions_are_borrowed_for_ten_and_twenty_copies() {
        let reference = reference(b"AAAAA");
        let sequence = &reference.sequences()[0];
        let requested_copies = std::env::var("GPS_VCF_INSERTION_COPIES")
            .ok()
            .and_then(|value| value.parse::<u8>().ok());
        let copy_counts: &[u8] = requested_copies
            .as_ref()
            .map_or(&[10, 20], std::slice::from_ref);
        for &copies in copy_counts {
            let identities = (0..copies)
                .map(|copy_index| CopyIdentity {
                    reference_id: "REF_main".to_owned(),
                    individual_id: format!("IND_{copy_index:02}"),
                    sequence_id: "chr1".to_owned(),
                    copy_index,
                })
                .collect::<Vec<_>>();
            let pieces = (0..copies)
                .map(|_| {
                    vec![
                        TruthCurrentPiece::Source {
                            source: source(0, 5, 5),
                            snps: Vec::new(),
                        },
                        TruthCurrentPiece::Inserted {
                            // Es un caso grande de regresión, no una cota del
                            // dominio biológico: la ruta acepta cualquier
                            // longitud validable de inserción directa.
                            bases: "G".repeat(1_000_000),
                            provenance: explicit(),
                        },
                    ]
                })
                .collect::<Vec<_>>();
            let paths = identities
                .iter()
                .zip(&pieces)
                .map(|(identity, pieces)| CopyPath::from_truth(identity, pieces, sequence))
                .collect::<Result<Vec<_>, _>>()
                .expect("compact paths");

            for (path, source_pieces) in paths.iter().zip(&pieces) {
                let TruthCurrentPiece::Inserted {
                    bases: expected, ..
                } = &source_pieces[1]
                else {
                    panic!("inserted source piece");
                };
                let PathPiece::Inserted {
                    bases: observed, ..
                } = &path.pieces[1]
                else {
                    panic!("inserted path piece");
                };
                assert_eq!(observed.len(), expected.len());
                assert_eq!(observed.as_ptr(), expected.as_ptr(), "{copies} copies");
            }
            std::hint::black_box(&paths);
            std::hint::black_box(&pieces);
        }
    }

    #[test]
    fn left_aligns_insertions_and_deletions() {
        let reference = reference(b"CAAAA");
        let insertion = TruthModification::Insertion {
            event_id: "INS_01".into(),
            point: 4,
            bases: "A".into(),
            provenance: explicit(),
        };
        let deletion = TruthModification::Deletion {
            event_id: "DEL_01".into(),
            start: 3,
            end: 4,
            provenance: explicit(),
        };
        let insertion_copy = copy(
            0,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 4, 5),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(4, 5, 5),
                    snps: vec![],
                },
            ],
            vec![insertion],
        );
        let deletion_copy = copy(
            1,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 3, 5),
                    snps: vec![],
                },
                TruthCurrentPiece::Source {
                    source: source(4, 5, 5),
                    snps: vec![],
                },
            ],
            vec![deletion],
        );
        let truth = document(&reference, vec![insertion_copy, deletion_copy]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 2);
        assert!(alleles.iter().any(|allele| (
            allele.position,
            allele.reference.as_str(),
            allele.alternate.as_str()
        ) == (0, "C", "CA")));
        assert!(alleles.iter().any(|allele| (
            allele.position,
            allele.reference.as_str(),
            allele.alternate.as_str()
        ) == (0, "CA", "C")));
    }

    #[test]
    fn direct_insertion_left_alignment_replaces_the_rotated_anchor() {
        let reference = reference(b"CAAAA");
        let truth = direct_insertion_document(&reference, 4, "A");
        let validated = validate_truth_for_vcf(&truth, &reference).expect("validated truth");
        let projected = project_projected_alleles_from_validated_truth(validated, &reference)
            .expect("projection");

        assert_eq!(projected.len(), 1);
        assert!(matches!(
            projected[0].alternate,
            ProjectedAlternate::DirectInsertion { .. }
        ));
        assert_eq!(projected[0].position, 0);
        assert_eq!(projected[0].reference, "C");
        assert_eq!(projected[0].alternate_with_suffix("").materialize(), "CA");

        let alleles = project_normalized_alleles(&truth, &reference).expect("public projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(alleles[0].position, 0);
        assert_eq!(alleles[0].reference, "C");
        assert_eq!(alleles[0].alternate, "CA");
    }

    #[test]
    fn direct_insertions_match_historical_left_alignment_in_homopolymers() {
        struct Case {
            reference_bases: &'static [u8],
            point: u64,
            bases: &'static str,
            position: u64,
            reference: &'static str,
            alternate: &'static str,
        }

        let cases = [
            Case {
                reference_bases: b"GAAAAAA",
                point: 6,
                bases: "AA",
                position: 0,
                reference: "G",
                alternate: "GAA",
            },
            Case {
                reference_bases: b"CTTTT",
                point: 4,
                bases: "TT",
                position: 0,
                reference: "C",
                alternate: "CTT",
            },
            Case {
                reference_bases: b"CATT",
                point: 3,
                bases: "AT",
                position: 0,
                reference: "C",
                alternate: "CAT",
            },
        ];

        for case in cases {
            let reference = reference(case.reference_bases);
            let truth = direct_insertion_document(&reference, case.point, case.bases);
            let validated =
                validate_truth_for_vcf(&truth, &reference).expect("validated direct insertion");
            let projected = project_projected_alleles_from_validated_truth(validated, &reference)
                .expect("projection");

            assert_eq!(projected.len(), 1, "{:?}", case.reference_bases);
            assert!(matches!(
                projected[0].alternate,
                ProjectedAlternate::DirectInsertion { .. }
            ));
            assert_eq!(
                projected[0].position, case.position,
                "{:?}",
                case.reference_bases
            );
            assert_eq!(
                projected[0].reference, case.reference,
                "{:?}",
                case.reference_bases
            );
            assert_eq!(
                projected[0].alternate_with_suffix("").materialize(),
                case.alternate,
                "{:?}",
                case.reference_bases
            );
        }
    }

    #[test]
    fn direct_insertions_match_historical_normalization_for_short_repeats() {
        for encoded_reference in 0_u8..32 {
            let reference_bases = (0..5)
                .map(|index| {
                    if encoded_reference & (1 << index) == 0 {
                        b'A'
                    } else {
                        b'C'
                    }
                })
                .collect::<Vec<_>>();
            let reference = reference(&reference_bases);
            let sequence = &reference.sequences()[0];

            for point in 1..5 {
                for bases in ["A", "C", "AA", "AC", "CA", "CC"] {
                    let copy = copy_with_length(
                        0,
                        5,
                        vec![
                            TruthCurrentPiece::Source {
                                source: source(0, point, 5),
                                snps: vec![],
                            },
                            TruthCurrentPiece::Inserted {
                                bases: bases.into(),
                                provenance: explicit(),
                            },
                            TruthCurrentPiece::Source {
                                source: source(point, 5, 5),
                                snps: vec![],
                            },
                        ],
                        Vec::new(),
                    );
                    let path = CopyPath::from_truth(&copy.identity, &copy.current_pieces, sequence)
                        .expect("direct path");
                    let actual = normalize_direct_insertion(
                        &path,
                        sequence,
                        ReferenceSpan {
                            start: point,
                            end: point,
                        },
                        bases,
                    )
                    .expect("compact normalization");

                    let mut expected_position = point - 1;
                    let mut expected_reference = vec![reference_bases[expected_position as usize]];
                    let mut expected_alternate = expected_reference.clone();
                    expected_alternate.extend_from_slice(bases.as_bytes());
                    left_align(
                        &reference_bases,
                        &mut expected_position,
                        &mut expected_reference,
                        &mut expected_alternate,
                    );

                    assert_eq!(
                        actual.position, expected_position,
                        "{reference_bases:?}, {point}, {bases}"
                    );
                    assert_eq!(
                        actual.reference.as_bytes(),
                        expected_reference,
                        "{reference_bases:?}, {point}, {bases}"
                    );
                    assert_eq!(
                        actual.alternate_with_suffix("").materialize().as_bytes(),
                        expected_alternate,
                        "{reference_bases:?}, {point}, {bases}"
                    );
                }
            }
        }
    }

    #[test]
    fn expands_overlapping_changes_to_one_common_reference_interval_before_trimming() {
        let reference = reference(b"AACCG");
        let mut altered = copy(
            0,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 2, 5),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: "T".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(4, 5, 5),
                    snps: vec![],
                },
            ],
            vec![],
        );
        // La inserción se produjo después de la deleción, en coordenadas
        // vigentes: ambas transformaciones permanecen solapadas respecto a la
        // referencia y la secuencia final es AATG.
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Deletion {
                event_id: "DEL_01".into(),
                start: 2,
                end: 4,
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_01".into(),
                point: 2,
                bases: "T".into(),
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (2, "CC", "T")
        );
    }

    #[test]
    fn orders_projected_alleles_by_reference_catalog_and_coordinate() {
        let reference = reference(b"AACGT");
        let later = TruthModification::Snp {
            event_id: "SNP_LATER".into(),
            position: 4,
            reference: "T".into(),
            alternate: "A".into(),
            provenance: explicit(),
        };
        let earlier = TruthModification::Snp {
            event_id: "SNP_EARLIER".into(),
            position: 1,
            reference: "A".into(),
            alternate: "C".into(),
            provenance: explicit(),
        };
        let mut first = copy(
            0,
            vec![TruthCurrentPiece::Source {
                source: source(0, 5, 5),
                snps: vec![TruthCurrentSnp {
                    event_id: "SNP_LATER".into(),
                    position: 4,
                    reference: "T".into(),
                    alternate: "A".into(),
                    provenance: explicit(),
                }],
            }],
            vec![later],
        );
        first.history = vec![crate::truth_schema::TruthHistoryEvent::Snp {
            event_id: "SNP_LATER".into(),
            position: 4,
            provenance: explicit(),
        }];
        let mut second = copy(
            1,
            vec![TruthCurrentPiece::Source {
                source: source(0, 5, 5),
                snps: vec![TruthCurrentSnp {
                    event_id: "SNP_EARLIER".into(),
                    position: 1,
                    reference: "A".into(),
                    alternate: "C".into(),
                    provenance: explicit(),
                }],
            }],
            vec![earlier],
        );
        second.history = vec![crate::truth_schema::TruthHistoryEvent::Snp {
            event_id: "SNP_EARLIER".into(),
            position: 1,
            provenance: explicit(),
        }];
        let truth = document(&reference, vec![first, second]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(
            alleles
                .iter()
                .map(|allele| allele.position)
                .collect::<Vec<_>>(),
            vec![1, 4]
        );
    }

    #[test]
    fn recomposes_nonoverlapping_insertions_that_collide_after_left_alignment() {
        let reference = reference(b"CAAAA");
        let mut altered = copy(
            0,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 3, 5),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(3, 4, 5),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(4, 5, 5),
                    snps: vec![],
                },
            ],
            vec![],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_01".into(),
                point: 3,
                bases: "A".into(),
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_02".into(),
                point: 5,
                bases: "A".into(),
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (0, "C", "CAA")
        );
        let reconstructed = apply_normalized(b"CAAAA", &alleles[0]);
        assert_eq!(reconstructed, b"CAAAAAA");
        assert_eq!(reconstructed.len(), 7);
    }

    #[test]
    fn recomposition_absorbs_an_interior_snp_after_distant_insertions_collide() {
        let reference = reference(b"CAAAAAAAA");
        let snp = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 5,
            reference: "A".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let mut altered = copy_with_length(
            0,
            9,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 3, 9),
                    snps: vec![],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(3, 7, 9),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_01".into(),
                        position: 5,
                        reference: "A".into(),
                        alternate: "T".into(),
                        provenance: explicit(),
                    }],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(7, 9, 9),
                    snps: vec![],
                },
            ],
            vec![snp],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_01".into(),
                point: 3,
                bases: "A".into(),
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_02".into(),
                point: 8,
                bases: "A".into(),
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_01".into(),
                position: 5,
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (4, "A", "AAT")
        );
        assert_eq!(apply_normalized(b"CAAAAAAAA", &alleles[0]), b"CAAAAATAAAA");
    }

    #[test]
    fn recomposes_a_snp_and_left_aligned_deletion_when_normalized_ref_intervals_overlap() {
        let reference = reference(b"CAAAA");
        let snp = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 1,
            reference: "A".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let mut altered = copy(
            0,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 3, 5),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_01".into(),
                        position: 1,
                        reference: "A".into(),
                        alternate: "T".into(),
                        provenance: explicit(),
                    }],
                },
                TruthCurrentPiece::Source {
                    source: source(4, 5, 5),
                    snps: vec![],
                },
            ],
            vec![snp],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_01".into(),
                position: 1,
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Deletion {
                event_id: "DEL_01".into(),
                start: 3,
                end: 4,
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (1, "AA", "T")
        );
        assert_eq!(apply_normalized(b"CAAAA", &alleles[0]), b"CTAA");
    }

    #[test]
    fn recomposition_absorbs_a_snp_crossed_while_left_aligning_a_deletion() {
        let reference = reference(b"CAAAAAAA");
        let snp = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 3,
            reference: "A".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let mut altered = copy_with_length(
            0,
            8,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 6, 8),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_01".into(),
                        position: 3,
                        reference: "A".into(),
                        alternate: "T".into(),
                        provenance: explicit(),
                    }],
                },
                TruthCurrentPiece::Source {
                    source: source(7, 8, 8),
                    snps: vec![],
                },
            ],
            vec![snp],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_01".into(),
                position: 3,
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Deletion {
                event_id: "DEL_01".into(),
                start: 6,
                end: 7,
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str()
            ),
            (3, "AA", "T")
        );
        assert_eq!(
            (alleles[0].source_start, alleles[0].source_end),
            (3, 7),
            "el DTO conserva el componente fuente compuesto para F2.5.3"
        );
        assert_eq!(apply_normalized(b"CAAAAAAA", &alleles[0]), b"CAATAAA");
    }

    #[test]
    fn insertion_at_a_boundary_bridges_two_snps_into_one_final_allele_once() {
        let reference = reference(b"ACG");
        let first = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 1,
            reference: "C".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let second = TruthModification::Snp {
            event_id: "SNP_02".into(),
            position: 2,
            reference: "G".into(),
            alternate: "A".into(),
            provenance: explicit(),
        };
        let mut altered = copy_with_length(
            0,
            3,
            vec![
                TruthCurrentPiece::Source {
                    source: source(0, 2, 3),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_01".into(),
                        position: 1,
                        reference: "C".into(),
                        alternate: "T".into(),
                        provenance: explicit(),
                    }],
                },
                TruthCurrentPiece::Inserted {
                    bases: "A".into(),
                    provenance: explicit(),
                },
                TruthCurrentPiece::Source {
                    source: source(2, 3, 3),
                    snps: vec![TruthCurrentSnp {
                        event_id: "SNP_02".into(),
                        position: 2,
                        reference: "G".into(),
                        alternate: "A".into(),
                        provenance: explicit(),
                    }],
                },
            ],
            vec![first, second],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_01".into(),
                position: 1,
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Insertion {
                event_id: "INS_01".into(),
                point: 2,
                bases: "A".into(),
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_02".into(),
                position: 2,
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 1);
        assert_eq!(
            (
                alleles[0].position,
                alleles[0].reference.as_str(),
                alleles[0].alternate.as_str(),
                alleles[0].source_start,
                alleles[0].source_end,
            ),
            (1, "CG", "TAA", 1, 3)
        );
        let reconstructed = apply_normalized(b"ACG", &alleles[0]);
        assert_eq!(reconstructed, b"ATAA");
        assert_eq!(reconstructed.len(), 4);
    }

    #[test]
    fn keeps_adjacent_snps_as_two_ordered_components() {
        let reference = reference(b"AACGT");
        let first = TruthModification::Snp {
            event_id: "SNP_01".into(),
            position: 1,
            reference: "A".into(),
            alternate: "C".into(),
            provenance: explicit(),
        };
        let second = TruthModification::Snp {
            event_id: "SNP_02".into(),
            position: 2,
            reference: "C".into(),
            alternate: "T".into(),
            provenance: explicit(),
        };
        let mut altered = copy(
            0,
            vec![TruthCurrentPiece::Source {
                source: source(0, 5, 5),
                snps: vec![
                    TruthCurrentSnp {
                        event_id: "SNP_01".into(),
                        position: 1,
                        reference: "A".into(),
                        alternate: "C".into(),
                        provenance: explicit(),
                    },
                    TruthCurrentSnp {
                        event_id: "SNP_02".into(),
                        position: 2,
                        reference: "C".into(),
                        alternate: "T".into(),
                        provenance: explicit(),
                    },
                ],
            }],
            vec![first, second],
        );
        altered.history = vec![
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_01".into(),
                position: 1,
                provenance: explicit(),
            },
            crate::truth_schema::TruthHistoryEvent::Snp {
                event_id: "SNP_02".into(),
                position: 2,
                provenance: explicit(),
            },
        ];
        let truth = document(&reference, vec![altered]);
        let alleles = project_normalized_alleles(&truth, &reference).expect("projection");
        assert_eq!(alleles.len(), 2);
        assert_eq!(
            alleles
                .iter()
                .map(|allele| (
                    allele.position,
                    allele.reference.as_str(),
                    allele.alternate.as_str()
                ))
                .collect::<Vec<_>>(),
            vec![(1, "A", "C"), (2, "C", "T")]
        );
    }

    /// Oráculo externo DTC-066, deliberadamente separado del runtime y de las
    /// puertas verdes locales: requiere bcftools 1.24, que el entorno actual
    /// todavía no suministra. Comprueba que `norm` no repara los vectores.
    #[test]
    fn external_bcftools_norm_oracle_matches_frozen_small_allele_vectors_without_repair() {
        let observed = Command::new("bcftools")
            .arg("--version")
            .output()
            .expect("bcftools installed");
        assert_eq!(
            String::from_utf8_lossy(&observed.stdout).lines().next(),
            Some("bcftools 1.24")
        );
        let root = std::env::temp_dir().join(format!(
            "gps-f2-5-2-external-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&root).expect("temporary directory");
        let fasta = root.join("reference.fa");
        let vcf = root.join("normalized.vcf");
        fs::write(&fasta, b">chr1\nCAAAA\n").expect("reference");
        fs::write(
            &vcf,
            concat!(
                "##fileformat=VCFv4.5\n",
                "##contig=<ID=chr1,length=5>\n",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
                "chr1\t1\t.\tC\tCA\t.\t.\t.\n",
                "chr1\t1\t.\tC\tCAA\t.\t.\t.\n",
                "chr1\t1\t.\tCA\tC\t.\t.\t.\n"
            ),
        )
        .expect("VCF");
        let normalized = Command::new("bcftools")
            .args([
                "norm",
                "--no-version",
                "--check-ref",
                "e",
                "--fasta-ref",
                fasta.to_str().expect("UTF-8 path"),
                vcf.to_str().expect("UTF-8 path"),
            ])
            .output()
            .expect("bcftools norm");
        fs::remove_dir_all(&root).expect("temporary cleanup");
        assert!(
            normalized.status.success(),
            "bcftools norm failed: {}",
            String::from_utf8_lossy(&normalized.stderr)
        );
        let normalized_text = String::from_utf8(normalized.stdout).expect("bcftools output UTF-8");
        let records = normalized_text
            .lines()
            .filter(|line| !line.starts_with('#'))
            .collect::<Vec<_>>();
        assert_eq!(
            records,
            [
                "chr1\t1\t.\tC\tCA\t.\t.\t.",
                "chr1\t1\t.\tC\tCAA\t.\t.\t.",
                "chr1\t1\t.\tCA\tC\t.\t.\t.",
            ]
        );
    }
}
