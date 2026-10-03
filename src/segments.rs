//! Segmentos haplotípicos y modificaciones pequeñas sin materializar referencias.
//!
//! Una [`SegmentChain`] describe la secuencia lógica de una copia haplotípica.
//! Sus segmentos ordinarios conservan solamente coordenadas de una fuente; las
//! únicas bases que posee el estado son las de una inserción nueva.

use std::sync::Arc;

use crate::coordinates::{BasePosition, GenomicInterval, InsertionAnchor, InternalDeletion};
use crate::identity::{EventId, PopulationId, RepeatLocusId, VariantId};
use crate::tandem_repeat::TandemRepeatTract;

/// Base concreta admitida como alelo de un SNP.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CanonicalBase {
    A,
    C,
    G,
    T,
}

impl CanonicalBase {
    pub fn from_byte(value: u8) -> Result<Self, SegmentError> {
        match value.to_ascii_uppercase() {
            b'A' => Ok(Self::A),
            b'C' => Ok(Self::C),
            b'G' => Ok(Self::G),
            b'T' => Ok(Self::T),
            _ => Err(SegmentError::NonCanonicalSnpBase { value }),
        }
    }

    #[must_use]
    pub const fn as_byte(self) -> u8 {
        match self {
            Self::A => b'A',
            Self::C => b'C',
            Self::G => b'G',
            Self::T => b'T',
        }
    }
}

/// Procedencia explícita de una modificación pequeña.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ModificationProvenance {
    event_id: EventId,
    origin: ModificationOrigin,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ModificationOrigin {
    Explicit,
    Random,
}

impl ModificationProvenance {
    #[must_use]
    pub const fn new(event_id: EventId) -> Self {
        Self {
            event_id,
            origin: ModificationOrigin::Explicit,
        }
    }

    #[must_use]
    pub const fn random(event_id: EventId) -> Self {
        Self {
            event_id,
            origin: ModificationOrigin::Random,
        }
    }

    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub const fn origin(&self) -> ModificationOrigin {
        self.origin
    }
}

/// Sustitución de una base de fuente por otra base canónica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SnpModification {
    position: BasePosition,
    reference: CanonicalBase,
    alternate: CanonicalBase,
    provenance: ModificationProvenance,
    matrix_variant_id: Option<VariantId>,
    matrix_population_id: Option<PopulationId>,
    matrix_source_id: Option<String>,
}

impl SnpModification {
    pub fn new(
        position: BasePosition,
        reference: CanonicalBase,
        alternate: CanonicalBase,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        if reference == alternate {
            return Err(SegmentError::SnpDoesNotChangeBase);
        }
        Ok(Self {
            position,
            reference,
            alternate,
            provenance,
            matrix_variant_id: None,
            matrix_population_id: None,
            matrix_source_id: None,
        })
    }

    /// Construye una sustitución procedente de una fila de matriz canónica.
    pub(crate) fn from_population_matrix(
        position: BasePosition,
        reference: CanonicalBase,
        alternate: CanonicalBase,
        provenance: ModificationProvenance,
        matrix_variant_id: VariantId,
        matrix_population_id: PopulationId,
        matrix_source_id: String,
    ) -> Result<Self, SegmentError> {
        let mut modification = Self::new(position, reference, alternate, provenance)?;
        modification.matrix_variant_id = Some(matrix_variant_id);
        modification.matrix_population_id = Some(matrix_population_id);
        modification.matrix_source_id = Some(matrix_source_id);
        Ok(modification)
    }

    #[must_use]
    pub const fn position(&self) -> &BasePosition {
        &self.position
    }

    #[must_use]
    pub const fn reference(&self) -> CanonicalBase {
        self.reference
    }

    #[must_use]
    pub const fn alternate(&self) -> CanonicalBase {
        self.alternate
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }

    /// Identidad de la fila de matriz que originó el SNP, si corresponde.
    #[must_use]
    pub const fn matrix_variant_id(&self) -> Option<&VariantId> {
        self.matrix_variant_id.as_ref()
    }

    /// Población de origen de la frecuencia de matriz, si corresponde.
    #[must_use]
    pub const fn matrix_population_id(&self) -> Option<&PopulationId> {
        self.matrix_population_id.as_ref()
    }

    /// Identificador de procedencia de la fila canónica, incluso si es vacío.
    #[must_use]
    pub fn matrix_source_id(&self) -> Option<&str> {
        self.matrix_source_id.as_deref()
    }
}

/// Inserción con bases propias, pues no procede de un intervalo de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InsertionModification {
    anchor: InsertionAnchor,
    bases: Arc<[u8]>,
    provenance: ModificationProvenance,
}

impl InsertionModification {
    pub fn new(
        anchor: InsertionAnchor,
        bases: impl Into<Arc<[u8]>>,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        let bases = bases.into();
        if bases.is_empty() {
            return Err(SegmentError::EmptyInsertion);
        }
        for &base in bases.iter() {
            if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
                return Err(SegmentError::InvalidInsertedBase { base });
            }
        }
        Ok(Self {
            anchor,
            bases,
            provenance,
        })
    }

    #[must_use]
    pub const fn anchor(&self) -> &InsertionAnchor {
        &self.anchor
    }

    /// Las bases devueltas son el único contenido secuencial poseído por el estado.
    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Deleción pequeña de un intervalo interno de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DeletionModification {
    deletion: InternalDeletion,
    provenance: ModificationProvenance,
}

impl DeletionModification {
    #[must_use]
    pub const fn new(deletion: InternalDeletion, provenance: ModificationProvenance) -> Self {
        Self {
            deletion,
            provenance,
        }
    }

    #[must_use]
    pub const fn deletion(&self) -> &InternalDeletion {
        &self.deletion
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

/// Una modificación pequeña tipada sobre un segmento de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SmallModification {
    Snp(SnpModification),
    Insertion(InsertionModification),
    Deletion(DeletionModification),
}

/// Segmento ordinario que puede reconstruirse desde su intervalo de fuente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct Segment {
    source: GenomicInterval,
    modifications: Vec<SmallModification>,
}

impl Segment {
    pub fn new(
        source: GenomicInterval,
        modifications: Vec<SmallModification>,
    ) -> Result<Self, SegmentError> {
        let mut previous: Option<ModificationExtent> = None;
        for modification in &modifications {
            let extent = ModificationExtent::from_modification(modification);
            if !modification_matches_source(modification, &source) {
                return Err(SegmentError::ModificationSourceMismatch);
            }
            validate_extent_in_source(&extent, &source)?;
            if let Some(previous) = previous {
                if previous.sort_key() >= extent.sort_key() {
                    return Err(SegmentError::ModificationsNotStrictlyOrdered);
                }
                if previous.overlaps(extent) {
                    return Err(SegmentError::OverlappingModifications);
                }
            }
            previous = Some(extent);
        }
        Ok(Self {
            source,
            modifications,
        })
    }

    #[must_use]
    pub const fn source(&self) -> &GenomicInterval {
        &self.source
    }

    #[must_use]
    pub fn modifications(&self) -> &[SmallModification] {
        &self.modifications
    }
}

/// Cadena ordenada, persistente y compartible de segmentos de una copia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct SegmentChain {
    segments: Arc<[Segment]>,
    /// Vista canónica de la secuencia vigente. `segments` conserva la
    /// procedencia estructural de compatibilidad, pero no es autoridad una vez
    /// que existen modificaciones pequeñas compuestas.
    indel_state: Arc<CurrentIndelState>,
}

/// Edición pequeña aceptada en las coordenadas de la secuencia vigente.
///
/// `point`, `start` y `end` son coordenadas **vigentes**, las del momento en
/// que se aplicó el evento, y así las replica el historial. No hay una
/// traducción general a la referencia: después de una inserción o una
/// deleción una base vigente puede no tener base de referencia. Lo que sí
/// existe siempre es el **locus del caso** en coordenadas de la referencia,
/// que se fija al aplicar el evento y se guarda aparte (`reference_*`):
///
/// - una inserción se ancla donde termina la última base de referencia que la
///   precede, la misma regla con la que el VCF ancla un fragmento insertado;
/// - una deleción cubre desde la primera hasta la última base de referencia
///   que elimina; si solo elimina bases insertadas, es el ancla de su inicio y
///   queda vacía.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CurrentIndelEvent {
    Snp {
        snp: SnpModification,
    },
    Insertion {
        event_id: EventId,
        point: u64,
        bases: Arc<[u8]>,
        provenance: ModificationProvenance,
        /// Punto de inserción en coordenadas de la referencia.
        reference_point: u64,
    },
    Deletion {
        event_id: EventId,
        start: u64,
        end: u64,
        provenance: ModificationProvenance,
        /// Intervalo de referencia que elimina, semiabierto.
        reference_start: u64,
        reference_end: u64,
    },
    RepeatAllele {
        allele: CurrentRepeatAllelePiece,
    },
}

impl CurrentIndelEvent {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        match self {
            Self::Snp { snp } => snp.provenance().event_id(),
            Self::RepeatAllele { allele } => allele.provenance().event_id(),
            Self::Insertion { event_id, .. } | Self::Deletion { event_id, .. } => event_id,
        }
    }
}

/// Fragmento canónico de una copia después de todas sus modificaciones
/// pequeñas vigentes. Las ausencias entre fragmentos de fuente representan
/// deleciones; el historial de eventos se consulta por separado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum CurrentIndelPiece {
    Source(CurrentSourcePiece),
    Inserted(CurrentInsertedPiece),
    RepeatAllele(CurrentRepeatAllelePiece),
}

/// **Alelo de una repetición en tándem realizado en una copia haplotípica.**
///
/// Es una clase propia de variante pequeña, hermana del cambio de base, la
/// inserción y la deleción: sustituye el **intervalo de referencia** del locus
/// declarado por la **secuencia realizada** del tracto. No guarda número de
/// unidades, ni pureza, ni interrupciones: los tres se derivan de las bases
/// cada vez que se piden, de modo que estructuralmente no pueden
/// contradecirlas (decisión S-02).
///
/// # Por qué es una clase y no una inserción con una deleción al lado
///
/// Dos alelos con **el mismo número de unidades** y la interrupción en sitio
/// distinto son alelos distintos; un contador no puede expresarlo y un par
/// inserción/deleción perdería el locus y el motivo, que es lo que la verdad y
/// el VCF necesitan anotar.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CurrentRepeatAllelePiece {
    locus_id: RepeatLocusId,
    interval: GenomicInterval,
    tract: TandemRepeatTract,
    provenance: ModificationProvenance,
}

impl CurrentRepeatAllelePiece {
    #[must_use]
    pub const fn new(
        locus_id: RepeatLocusId,
        interval: GenomicInterval,
        tract: TandemRepeatTract,
        provenance: ModificationProvenance,
    ) -> Self {
        Self {
            locus_id,
            interval,
            tract,
            provenance,
        }
    }

    #[must_use]
    pub const fn locus_id(&self) -> &RepeatLocusId {
        &self.locus_id
    }

    /// Intervalo de la **referencia** que el alelo sustituye.
    #[must_use]
    pub const fn interval(&self) -> &GenomicInterval {
        &self.interval
    }

    #[must_use]
    pub const fn tract(&self) -> &TandemRepeatTract {
        &self.tract
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }

    fn length(&self) -> u64 {
        self.tract.length()
    }
}

/// Fragmento de referencia en el espacio de coordenadas vigente.
///
/// Los SNP conservan coordenadas de fuente incluso cuando una inserción o una
/// deleción posterior divide el fragmento. Esto permite que la eliminación de
/// un SNP sea una consecuencia del evento posterior, sin reconstruir ni
/// aplanar la copia.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CurrentSourcePiece {
    source: GenomicInterval,
    snps: Arc<[SnpModification]>,
}

impl CurrentSourcePiece {
    fn new(source: GenomicInterval, snps: impl Into<Arc<[SnpModification]>>) -> Self {
        Self {
            source,
            snps: snps.into(),
        }
    }

    #[must_use]
    pub const fn source(&self) -> &GenomicInterval {
        &self.source
    }

    #[must_use]
    pub fn snps(&self) -> &[SnpModification] {
        &self.snps
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CurrentInsertedPiece {
    bases: Arc<[u8]>,
    provenance: ModificationProvenance,
}

impl CurrentInsertedPiece {
    #[must_use]
    pub const fn bases(&self) -> &Arc<[u8]> {
        &self.bases
    }

    #[must_use]
    pub const fn provenance(&self) -> &ModificationProvenance {
        &self.provenance
    }
}

impl CurrentIndelPiece {
    fn len(&self) -> u64 {
        match self {
            Self::Source(source) => source.source.len(),
            Self::Inserted(inserted) => u64::try_from(inserted.bases.len()).unwrap_or(u64::MAX),
            Self::RepeatAllele(allele) => allele.length(),
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
struct CurrentIndelState {
    pieces: Arc<[CurrentIndelPiece]>,
    events: Arc<[CurrentIndelEvent]>,
}

impl SegmentChain {
    /// Crea la cadena intacta de una copia desde un intervalo ya validado.
    #[must_use]
    pub fn intact(source: GenomicInterval) -> Self {
        Self {
            segments: Arc::from(vec![Segment {
                source: source.clone(),
                modifications: Vec::new(),
            }]),
            indel_state: Arc::new(CurrentIndelState {
                pieces: Arc::from(vec![CurrentIndelPiece::Source(CurrentSourcePiece::new(
                    source,
                    Arc::from(Vec::<SnpModification>::new()),
                ))]),
                events: Arc::from(Vec::<CurrentIndelEvent>::new()),
            }),
        }
    }

    pub fn new(segments: Vec<Segment>) -> Result<Self, SegmentError> {
        if segments.is_empty() {
            return Err(SegmentError::EmptyChain);
        }
        let first = segments[0].source();
        if segments.iter().skip(1).any(|segment| {
            let source = segment.source();
            source.reference_id() != first.reference_id()
                || source.sequence_id() != first.sequence_id()
                || source.sequence_length() != first.sequence_length()
        }) {
            return Err(SegmentError::MixedSourceSpace);
        }
        let indel_state = current_state_from_segments(&segments)?;
        Ok(Self {
            segments: Arc::from(segments),
            indel_state: Arc::new(indel_state),
        })
    }

    /// Devuelve una nueva cadena con un SNP insertado en el único segmento de
    /// fuente que contiene su posición. La cadena original permanece intacta.
    pub fn with_snp(&self, snp: SnpModification) -> Result<Self, SegmentError> {
        let current = self.with_current_snp(snp.clone())?;
        let position = snp.position().offset();
        let matching = self
            .segments
            .iter()
            .enumerate()
            .filter(|(_, segment)| {
                let source = segment.source();
                snp.position().reference_id() == source.reference_id()
                    && snp.position().sequence_id() == source.sequence_id()
                    && snp.position().sequence_length() == source.sequence_length()
                    && (source.start()..source.end()).contains(&position)
            })
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let index = match matching.as_slice() {
            [] => return Err(SegmentError::SnpPositionNotCovered),
            [index] => *index,
            _ => return Err(SegmentError::AmbiguousSnpPosition),
        };

        let mut segments = self.segments.to_vec();
        let segment = &segments[index];
        if segment.modifications().iter().any(|modification| {
            matches!(modification, SmallModification::Snp(existing) if existing.position() == snp.position())
        }) {
            return Err(SegmentError::SnpPositionAlreadyModified);
        }
        if segment.modifications().iter().any(|modification| {
            matches!(modification, SmallModification::Deletion(deletion)
                if (deletion.deletion().interval().start()..deletion.deletion().interval().end())
                    .contains(&position))
        }) {
            return Err(SegmentError::SnpPositionCoveredByDeletion);
        }
        let new_extent =
            ModificationExtent::from_modification(&SmallModification::Snp(snp.clone()));
        let insert_at = segment
            .modifications()
            .iter()
            .position(|modification| {
                new_extent.sort_key()
                    < ModificationExtent::from_modification(modification).sort_key()
            })
            .unwrap_or(segment.modifications().len());
        let mut modifications = segment.modifications().to_vec();
        modifications.insert(insert_at, SmallModification::Snp(snp));
        segments[index] = Segment::new(segment.source().clone(), modifications)?;
        Ok(Self {
            segments: Arc::from(segments),
            indel_state: current.indel_state,
        })
    }

    /// Añade un SNP a un fragmento de fuente que sigue presente tras indels.
    fn with_current_snp(&self, snp: SnpModification) -> Result<Self, SegmentError> {
        let state = self.current_indel_state()?;
        let position = snp.position().offset();
        let matching = state
            .pieces
            .iter()
            .enumerate()
            .filter(|(_, piece)| {
                matches!(piece, CurrentIndelPiece::Source(source)
                    if snp.position().reference_id() == source.source.reference_id()
                        && snp.position().sequence_id() == source.source.sequence_id()
                        && snp.position().sequence_length() == source.source.sequence_length()
                        && (source.source.start()..source.source.end()).contains(&position))
            })
            .map(|(index, _)| index)
            .collect::<Vec<_>>();
        let index = match matching.as_slice() {
            [] => {
                let deleted_by_legacy_segment = self.segments.iter().any(|segment| {
                    segment.modifications().iter().any(|modification| {
                        matches!(modification, SmallModification::Deletion(deletion)
                            if deletion.deletion().interval().reference_id()
                                == snp.position().reference_id()
                                && deletion.deletion().interval().sequence_id()
                                    == snp.position().sequence_id()
                                && (deletion.deletion().interval().start()
                                    ..deletion.deletion().interval().end())
                                    .contains(&position))
                    })
                });
                return Err(
                    if deleted_by_legacy_segment
                        || state
                            .events
                            .iter()
                            .any(|event| matches!(event, CurrentIndelEvent::Deletion { .. }))
                    {
                        SegmentError::SnpPositionCoveredByDeletion
                    } else {
                        SegmentError::SnpPositionNotCovered
                    },
                );
            }
            [index] => *index,
            _ => return Err(SegmentError::AmbiguousSnpPosition),
        };

        let mut pieces = state.pieces.to_vec();
        let CurrentIndelPiece::Source(source) = &pieces[index] else {
            return Err(SegmentError::CurrentCoordinateMismatch);
        };
        if source
            .snps
            .iter()
            .any(|existing| existing.position() == snp.position())
        {
            return Err(SegmentError::SnpPositionAlreadyModified);
        }
        let mut snps = source.snps.to_vec();
        let insert_at = snps
            .iter()
            .position(|existing| existing.position().offset() > position)
            .unwrap_or(snps.len());
        let event = snp.clone();
        snps.insert(insert_at, snp);
        pieces[index] = CurrentIndelPiece::Source(CurrentSourcePiece::new(
            source.source.clone(),
            Arc::from(snps),
        ));
        Ok(Self {
            segments: self.segments.clone(),
            indel_state: Arc::new(CurrentIndelState {
                pieces: Arc::from(pieces),
                events: Arc::from({
                    let mut events = state.events.to_vec();
                    events.push(CurrentIndelEvent::Snp { snp: event });
                    events
                }),
            }),
        })
    }

    /// Inserta bases canónicas en un punto interno de la secuencia vigente.
    ///
    /// La cadena mantiene fragmentos de fuente y de inserción, por lo que una
    /// inserción posterior puede anclarse dentro de bases que no pertenecen a
    /// la referencia. F2.2.10 conserva además los SNP de los fragmentos de
    /// fuente para que eventos posteriores puedan eliminarlos o coexistir.
    pub(crate) fn with_current_insertion(
        &self,
        point: u64,
        bases: Arc<[u8]>,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        if bases.is_empty() {
            return Err(SegmentError::EmptyInsertion);
        }
        if bases
            .iter()
            .any(|base| !matches!(*base, b'A' | b'C' | b'G' | b'T'))
        {
            return Err(SegmentError::InvalidInsertedBase {
                base: *bases
                    .iter()
                    .find(|base| !matches!(**base, b'A' | b'C' | b'G' | b'T'))
                    .expect("a non-canonical base was found"),
            });
        }
        let state = self.current_indel_state()?;
        let length = piece_length(&state.pieces)?;
        if point == 0 || point >= length {
            return Err(SegmentError::CurrentInsertionMustBeInternal { point, length });
        }
        let mut pieces = split_pieces_at(&state.pieces, point)?;
        let insertion_index = pieces
            .iter()
            .scan(0_u64, |cursor, piece| {
                let here = *cursor;
                *cursor = cursor.checked_add(piece.len())?;
                Some((here, *cursor))
            })
            .position(|(_, end)| end == point)
            .ok_or(SegmentError::CurrentCoordinateMismatch)?
            .checked_add(1)
            .ok_or(SegmentError::ArithmeticOverflow)?;
        let reference_point = reference_anchor(&pieces[..insertion_index])?;
        pieces.insert(
            insertion_index,
            CurrentIndelPiece::Inserted(CurrentInsertedPiece {
                bases: bases.clone(),
                provenance: provenance.clone(),
            }),
        );
        let mut events = state.events.to_vec();
        events.push(CurrentIndelEvent::Insertion {
            event_id: provenance.event_id().clone(),
            point,
            bases,
            provenance,
            reference_point,
        });
        Ok(Self {
            segments: self.segments.clone(),
            indel_state: Arc::new(CurrentIndelState {
                pieces: Arc::from(pieces),
                events: Arc::from(events),
            }),
        })
    }

    /// Elimina un intervalo interno de la secuencia vigente sin aplanarla.
    pub(crate) fn with_current_deletion(
        &self,
        start: u64,
        end: u64,
        provenance: ModificationProvenance,
    ) -> Result<Self, SegmentError> {
        let state = self.current_indel_state()?;
        let length = piece_length(&state.pieces)?;
        if start == 0 || start >= end || end >= length {
            return Err(SegmentError::CurrentDeletionMustBeInternal { start, end, length });
        }
        let split = split_pieces_at(&state.pieces, start)?;
        let split = split_pieces_at(&split, end)?;
        let mut cursor = 0_u64;
        let mut pieces = Vec::with_capacity(split.len());
        // El ancla se toma al llegar al primer fragmento eliminado, cuando
        // `pieces` contiene exactamente los que lo preceden.
        let mut anchor = None;
        let mut removed_reference: Option<(u64, u64)> = None;
        for piece in split {
            let next = cursor
                .checked_add(piece.len())
                .ok_or(SegmentError::ArithmeticOverflow)?;
            if next <= start || cursor >= end {
                pieces.push(piece);
            } else {
                if anchor.is_none() {
                    anchor = Some(reference_anchor(&pieces)?);
                }
                if let Some((piece_start, piece_end)) = reference_interval(&piece) {
                    removed_reference = Some(match removed_reference {
                        None => (piece_start, piece_end),
                        Some((first, last)) => (first.min(piece_start), last.max(piece_end)),
                    });
                }
            }
            cursor = next;
        }
        let (reference_start, reference_end) = match (removed_reference, anchor) {
            (Some(span), _) => span,
            (None, Some(anchor)) => (anchor, anchor),
            (None, None) => return Err(SegmentError::CurrentCoordinateMismatch),
        };
        let mut events = state.events.to_vec();
        events.push(CurrentIndelEvent::Deletion {
            event_id: provenance.event_id().clone(),
            start,
            end,
            provenance,
            reference_start,
            reference_end,
        });
        Ok(Self {
            segments: self.segments.clone(),
            indel_state: Arc::new(CurrentIndelState {
                pieces: Arc::from(pieces),
                events: Arc::from(events),
            }),
        })
    }

    /// Sustituye el intervalo de referencia de un locus de repetición por la
    /// secuencia realizada de su alelo.
    ///
    /// # El intervalo tiene que llegar **intacto**, y por eso se aplica el último
    ///
    /// Un alelo de repetición describe el tracto **de la referencia**: su
    /// número de unidades, su pureza y sus interrupciones se leen respecto de
    /// esas bases. Si un cambio de base o un indel anterior ya tocó el tramo,
    /// el alelo dejaría de describir lo que sustituye y publicaría una verdad
    /// falsa. Por eso la operación exige que el intervalo entero viva dentro de
    /// **un solo** fragmento de fuente y **sin ningún SNP dentro**, y rechaza
    /// en otro caso identificando el conflicto.
    ///
    /// Es la razón de que en el recorrido los alelos de repetición se apliquen
    /// **después** de la matriz y de los cambios directos: así lo que la guarda
    /// ve es el estado final, y no un estado intermedio que todavía puede
    /// cambiar.
    pub fn with_current_repeat_allele(
        &self,
        allele: CurrentRepeatAllelePiece,
    ) -> Result<Self, SegmentError> {
        let state = self.current_indel_state()?;
        let interval = allele.interval();
        let index = self.intact_repeat_interval(interval)?;
        let CurrentIndelPiece::Source(source) = &state.pieces[index] else {
            return Err(SegmentError::RepeatIntervalNotIntact);
        };
        let mut pieces = Vec::with_capacity(state.pieces.len() + 2);
        pieces.extend_from_slice(&state.pieces[..index]);
        push_source_piece(
            &mut pieces,
            &source.source,
            source.source.start(),
            interval.start(),
            &source
                .snps
                .iter()
                .filter(|snp| snp.position().offset() < interval.start())
                .cloned()
                .collect::<Vec<_>>(),
        )?;
        pieces.push(CurrentIndelPiece::RepeatAllele(allele.clone()));
        push_source_piece(
            &mut pieces,
            &source.source,
            interval.end(),
            source.source.end(),
            &source
                .snps
                .iter()
                .filter(|snp| snp.position().offset() >= interval.end())
                .cloned()
                .collect::<Vec<_>>(),
        )?;
        pieces.extend_from_slice(&state.pieces[index + 1..]);
        let mut events = state.events.to_vec();
        events.push(CurrentIndelEvent::RepeatAllele { allele });
        Ok(Self {
            segments: self.segments.clone(),
            indel_state: Arc::new(CurrentIndelState {
                pieces: Arc::from(pieces),
                events: Arc::from(events),
            }),
        })
    }

    /// Comprueba que el intervalo de un locus llega **intacto** y devuelve el
    /// fragmento que lo contiene.
    ///
    /// # Por qué es pública y se consulta aunque no haya nada que aplicar
    ///
    /// Cuando la copia sortea el alelo que coincide con el tracto de la
    /// referencia no hay modificación que escribir, pero **la condición sigue
    /// teniendo que cumplirse**: la verdad va a publicar que esa copia
    /// «conserva el tracto de la referencia», y eso es falso si algo lo tocó.
    /// Si la comprobación viviera solo dentro de la operación que aplica el
    /// alelo, el rechazo dependería del sorteo y, con él, de la semilla.
    pub fn intact_repeat_interval(
        &self,
        interval: &GenomicInterval,
    ) -> Result<usize, SegmentError> {
        let (index, source) = self
            .indel_state
            .pieces
            .iter()
            .enumerate()
            .find_map(|(index, piece)| match piece {
                CurrentIndelPiece::Source(source)
                    if source.source.reference_id() == interval.reference_id()
                        && source.source.sequence_id() == interval.sequence_id()
                        && source.source.sequence_length() == interval.sequence_length()
                        && source.source.start() <= interval.start()
                        && interval.end() <= source.source.end() =>
                {
                    Some((index, source))
                }
                CurrentIndelPiece::Source(_)
                | CurrentIndelPiece::Inserted(_)
                | CurrentIndelPiece::RepeatAllele(_) => None,
            })
            .ok_or(SegmentError::RepeatIntervalNotIntact)?;
        if source
            .snps
            .iter()
            .any(|snp| (interval.start()..interval.end()).contains(&snp.position().offset()))
        {
            return Err(SegmentError::RepeatIntervalAlreadyModified);
        }
        Ok(index)
    }

    fn current_indel_state(&self) -> Result<CurrentIndelState, SegmentError> {
        Ok((*self.indel_state).clone())
    }

    #[must_use]
    pub fn current_indel_events(&self) -> &[CurrentIndelEvent] {
        &self.indel_state.events
    }

    /// Historial declarado de eventos pequeños, separado de la vista vigente.
    #[must_use]
    pub fn small_event_history(&self) -> &[CurrentIndelEvent] {
        &self.indel_state.events
    }

    /// Devuelve la única vista pública canónica de las modificaciones vigentes.
    ///
    /// La lista está en orden de secuencia final y es la autoridad para
    /// materialización y consultas genotípicas. Los segmentos heredados solo
    /// describen fuentes y no deben usarse para inferir una variante vigente.
    #[must_use]
    pub fn current_modifications(&self) -> &[CurrentIndelPiece] {
        &self.indel_state.pieces
    }

    /// Longitud de la secuencia vigente después de todas las modificaciones
    /// pequeñas, sin reconstruir sus bases.
    pub fn current_length(&self) -> Result<u64, SegmentError> {
        piece_length(&self.indel_state.pieces)
    }

    #[must_use]
    /// Descriptores de fuente conservados por compatibilidad estructural.
    ///
    /// No representan necesariamente el estado pequeño vigente después de una
    /// composición; para ello debe usarse [`Self::current_modifications`].
    pub fn segments(&self) -> &[Segment] {
        &self.segments
    }

    /// Indica si dos cadenas comparten realmente la misma asignación persistente.
    #[must_use]
    pub fn is_shared_with(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.segments, &other.segments)
    }
}

#[derive(Clone, Copy)]
struct ModificationExtent {
    start: u64,
    end: u64,
    is_insertion: bool,
    is_deletion: bool,
}

impl ModificationExtent {
    fn from_modification(modification: &SmallModification) -> Self {
        match modification {
            SmallModification::Snp(snp) => Self {
                start: snp.position.offset(),
                end: snp.position.offset() + 1,
                is_insertion: false,
                is_deletion: false,
            },
            SmallModification::Insertion(insertion) => Self {
                start: insertion.anchor.point().offset(),
                end: insertion.anchor.point().offset(),
                is_insertion: true,
                is_deletion: false,
            },
            SmallModification::Deletion(deletion) => Self {
                start: deletion.deletion.interval().start(),
                end: deletion.deletion.interval().end(),
                is_insertion: false,
                is_deletion: true,
            },
        }
    }

    const fn sort_key(self) -> (u64, u64, bool) {
        (self.start, self.end, self.is_insertion)
    }

    const fn overlaps(self, other: Self) -> bool {
        if self.is_insertion && other.is_insertion {
            self.start == other.start
        } else if self.is_insertion && other.is_deletion {
            other.start < self.start && self.start <= other.end
        } else if other.is_insertion && self.is_deletion {
            self.start < other.start && other.start <= self.end
        } else if self.is_insertion || other.is_insertion {
            false
        } else {
            self.start < other.end && other.start < self.end
        }
    }
}

fn validate_extent_in_source(
    extent: &ModificationExtent,
    source: &GenomicInterval,
) -> Result<(), SegmentError> {
    if extent.start < source.start() || extent.end > source.end() {
        return Err(SegmentError::ModificationOutsideSource);
    }
    Ok(())
}

fn modification_matches_source(modification: &SmallModification, source: &GenomicInterval) -> bool {
    let (reference_id, sequence_id, sequence_length) = match modification {
        SmallModification::Snp(value) => (
            value.position().reference_id(),
            value.position().sequence_id(),
            value.position().sequence_length(),
        ),
        SmallModification::Insertion(value) => (
            value.anchor().point().reference_id(),
            value.anchor().point().sequence_id(),
            value.anchor().point().sequence_length(),
        ),
        SmallModification::Deletion(value) => (
            value.deletion().interval().reference_id(),
            value.deletion().interval().sequence_id(),
            value.deletion().interval().sequence_length(),
        ),
    };
    reference_id == source.reference_id()
        && sequence_id == source.sequence_id()
        && sequence_length == source.sequence_length()
}

/// Invariantes que evitan una cadena o una composición ambigua.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SegmentError {
    NonCanonicalSnpBase {
        value: u8,
    },
    SnpDoesNotChangeBase,
    EmptyInsertion,
    InvalidInsertedBase {
        base: u8,
    },
    EmptyChain,
    MixedSourceSpace,
    ModificationSourceMismatch,
    ModificationOutsideSource,
    ModificationsNotStrictlyOrdered,
    OverlappingModifications,
    SnpPositionNotCovered,
    AmbiguousSnpPosition,
    SnpPositionAlreadyModified,
    SnpPositionCoveredByDeletion,
    SmallIndelCompositionDeferred,
    CurrentInsertionMustBeInternal {
        point: u64,
        length: u64,
    },
    CurrentDeletionMustBeInternal {
        start: u64,
        end: u64,
        length: u64,
    },
    CurrentCoordinateMismatch,
    ArithmeticOverflow,
    /// El intervalo del locus de repetición no llega entero dentro de un solo
    /// fragmento de fuente: un indel anterior lo partió, lo borró o ya lo
    /// sustituyó otro alelo.
    RepeatIntervalNotIntact,
    /// Un cambio de base anterior cae dentro del intervalo del locus, de modo
    /// que el tracto que el alelo dice sustituir ya no es el de la referencia.
    RepeatIntervalAlreadyModified,
    /// Una inserción o una deleción posterior caería **dentro** del tracto
    /// realizado. No se admite: partir el tracto dejaría dos mitades que ya no
    /// son un alelo y la anotación no podría derivarse de ellas.
    RepeatAlleleIsNotDivisible,
}

fn piece_length(pieces: &[CurrentIndelPiece]) -> Result<u64, SegmentError> {
    pieces.iter().try_fold(0_u64, |total, piece| {
        total
            .checked_add(piece.len())
            .ok_or(SegmentError::ArithmeticOverflow)
    })
}

fn current_state_from_segments(segments: &[Segment]) -> Result<CurrentIndelState, SegmentError> {
    let mut pieces = Vec::new();
    for segment in segments {
        let source = segment.source();
        let mut cursor = source.start();
        for modification in segment.modifications() {
            match modification {
                SmallModification::Snp(snp) => {
                    push_source_piece(&mut pieces, source, cursor, snp.position().offset(), &[])?;
                    let end = snp
                        .position()
                        .offset()
                        .checked_add(1)
                        .ok_or(SegmentError::ArithmeticOverflow)?;
                    push_source_piece(
                        &mut pieces,
                        source,
                        snp.position().offset(),
                        end,
                        std::slice::from_ref(snp),
                    )?;
                    cursor = end;
                }
                SmallModification::Insertion(insertion) => {
                    let point = insertion.anchor().point().offset();
                    push_source_piece(&mut pieces, source, cursor, point, &[])?;
                    pieces.push(CurrentIndelPiece::Inserted(CurrentInsertedPiece {
                        bases: Arc::from(insertion.bases()),
                        provenance: insertion.provenance().clone(),
                    }));
                    cursor = point;
                }
                SmallModification::Deletion(deletion) => {
                    let interval = deletion.deletion().interval();
                    push_source_piece(&mut pieces, source, cursor, interval.start(), &[])?;
                    cursor = interval.end();
                }
            }
        }
        push_source_piece(&mut pieces, source, cursor, source.end(), &[])?;
    }
    Ok(CurrentIndelState {
        pieces: Arc::from(pieces),
        events: Arc::from(Vec::<CurrentIndelEvent>::new()),
    })
}

fn push_source_piece(
    pieces: &mut Vec<CurrentIndelPiece>,
    source: &GenomicInterval,
    start: u64,
    end: u64,
    snps: &[SnpModification],
) -> Result<(), SegmentError> {
    if start == end {
        return Ok(());
    }
    let interval = GenomicInterval::new(
        source.reference_id().clone(),
        source.sequence_id().clone(),
        start,
        end,
        source.sequence_length(),
    )
    .map_err(|_| SegmentError::CurrentCoordinateMismatch)?;
    pieces.push(CurrentIndelPiece::Source(CurrentSourcePiece::new(
        interval,
        Arc::from(snps),
    )));
    Ok(())
}

/// Coordenada de referencia que sigue a la última base de referencia de
/// `before`, los fragmentos vigentes que preceden a un límite.
///
/// Es el ancla con la que el VCF coloca un fragmento insertado. La primera
/// base de una copia es siempre de referencia —no se inserta ni se elimina en
/// el punto cero—, así que no encontrarla es una coordenada imposible.
fn reference_anchor(before: &[CurrentIndelPiece]) -> Result<u64, SegmentError> {
    before
        .iter()
        .rev()
        .find_map(|piece| reference_interval(piece).map(|(_, end)| end))
        .ok_or(SegmentError::CurrentCoordinateMismatch)
}

/// Intervalo de referencia que ocupa un fragmento vigente, si lo tiene.
fn reference_interval(piece: &CurrentIndelPiece) -> Option<(u64, u64)> {
    match piece {
        CurrentIndelPiece::Source(source) => Some((source.source.start(), source.source.end())),
        CurrentIndelPiece::RepeatAllele(allele) => {
            Some((allele.interval().start(), allele.interval().end()))
        }
        CurrentIndelPiece::Inserted(_) => None,
    }
}

fn split_pieces_at(
    pieces: &[CurrentIndelPiece],
    point: u64,
) -> Result<Vec<CurrentIndelPiece>, SegmentError> {
    let length = piece_length(pieces)?;
    if point > length {
        return Err(SegmentError::CurrentCoordinateMismatch);
    }
    let mut result = Vec::with_capacity(pieces.len() + 1);
    let mut cursor = 0_u64;
    for piece in pieces {
        let next = cursor
            .checked_add(piece.len())
            .ok_or(SegmentError::ArithmeticOverflow)?;
        if cursor < point && point < next {
            let left_len = point - cursor;
            let right_len = next - point;
            match piece {
                CurrentIndelPiece::Source(source) => {
                    result.push(CurrentIndelPiece::Source(CurrentSourcePiece::new(
                        GenomicInterval::new(
                            source.source.reference_id().clone(),
                            source.source.sequence_id().clone(),
                            source.source.start(),
                            source
                                .source
                                .start()
                                .checked_add(left_len)
                                .ok_or(SegmentError::ArithmeticOverflow)?,
                            source.source.sequence_length(),
                        )
                        .map_err(|_| SegmentError::CurrentCoordinateMismatch)?,
                        Arc::from(
                            source
                                .snps
                                .iter()
                                .filter(|snp| {
                                    snp.position().offset() < source.source.start() + left_len
                                })
                                .cloned()
                                .collect::<Vec<_>>(),
                        ),
                    )));
                    result.push(CurrentIndelPiece::Source(CurrentSourcePiece::new(
                        GenomicInterval::new(
                            source.source.reference_id().clone(),
                            source.source.sequence_id().clone(),
                            source
                                .source
                                .end()
                                .checked_sub(right_len)
                                .ok_or(SegmentError::ArithmeticOverflow)?,
                            source.source.end(),
                            source.source.sequence_length(),
                        )
                        .map_err(|_| SegmentError::CurrentCoordinateMismatch)?,
                        Arc::from(
                            source
                                .snps
                                .iter()
                                .filter(|snp| {
                                    snp.position().offset() >= source.source.end() - right_len
                                })
                                .cloned()
                                .collect::<Vec<_>>(),
                        ),
                    )));
                }
                CurrentIndelPiece::Inserted(inserted) => {
                    let left =
                        usize::try_from(left_len).map_err(|_| SegmentError::ArithmeticOverflow)?;
                    result.push(CurrentIndelPiece::Inserted(CurrentInsertedPiece {
                        bases: Arc::from(&inserted.bases[..left]),
                        provenance: inserted.provenance.clone(),
                    }));
                    result.push(CurrentIndelPiece::Inserted(CurrentInsertedPiece {
                        bases: Arc::from(&inserted.bases[left..]),
                        provenance: inserted.provenance.clone(),
                    }));
                }
                CurrentIndelPiece::RepeatAllele(_) => {
                    return Err(SegmentError::RepeatAlleleIsNotDivisible);
                }
            }
        } else {
            result.push(piece.clone());
        }
        cursor = next;
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::coordinates::{
        GenomicInterval, InsertionAnchor, InterbasePosition, InternalDeletion,
    };
    use crate::identity::{EventId, ReferenceId, SequenceId};

    fn interval(start: u64, end: u64) -> GenomicInterval {
        GenomicInterval::new(
            ReferenceId::parse("ref", "test").expect("reference"),
            SequenceId::parse("chr1", "test").expect("sequence"),
            start,
            end,
            10,
        )
        .expect("interval")
    }

    fn provenance() -> ModificationProvenance {
        ModificationProvenance::new(EventId::parse("event-1", "test").expect("event"))
    }

    #[test]
    fn retains_an_intact_source_segment_without_reference_bases() {
        let segment = Segment::new(interval(0, 10), vec![]).expect("intact segment");
        assert_eq!(segment.source().len(), 10);
        assert!(segment.modifications().is_empty());
    }

    #[test]
    fn keeps_each_small_edit_typed_and_its_provenance() {
        let snp = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::G,
                CanonicalBase::T,
                provenance(),
            )
            .expect("SNP"),
        );
        let insertion = SmallModification::Insertion(
            InsertionModification::new(
                InsertionAnchor::new(
                    InterbasePosition::new(
                        ReferenceId::parse("ref", "test").expect("reference"),
                        SequenceId::parse("chr1", "test").expect("sequence"),
                        4,
                        10,
                    )
                    .expect("point"),
                )
                .expect("anchor"),
                Arc::<[u8]>::from(&b"NN"[..]),
                provenance(),
            )
            .expect("insertion"),
        );
        let deletion = SmallModification::Deletion(DeletionModification::new(
            InternalDeletion::new(interval(6, 8)).expect("deletion"),
            provenance(),
        ));
        let segment =
            Segment::new(interval(0, 10), vec![snp, insertion, deletion]).expect("ordered edits");
        assert_eq!(segment.modifications().len(), 3);
        assert!(matches!(
            segment.modifications()[0],
            SmallModification::Snp(_)
        ));
        assert!(matches!(
            segment.modifications()[1],
            SmallModification::Insertion(_)
        ));
        assert!(matches!(
            segment.modifications()[2],
            SmallModification::Deletion(_)
        ));
        for modification in segment.modifications() {
            let observed = match modification {
                SmallModification::Snp(value) => value.provenance(),
                SmallModification::Insertion(value) => value.provenance(),
                SmallModification::Deletion(value) => value.provenance(),
            };
            assert_eq!(observed.event_id().as_str(), "event-1");
        }
    }

    #[test]
    fn cloned_chains_share_segments_instead_of_reference_or_segment_copies() {
        let chain = SegmentChain::new(vec![
            Segment::new(interval(0, 10), vec![]).expect("segment"),
        ])
        .expect("chain");
        let shared = chain.clone();
        assert!(chain.is_shared_with(&shared));
        assert_eq!(chain.segments()[0].source().reference_id().as_str(), "ref");
    }

    #[test]
    fn preserves_the_declared_order_of_a_short_known_composition() {
        let chain = SegmentChain::new(vec![
            Segment::new(interval(0, 3), vec![]).expect("first segment"),
            Segment::new(interval(7, 10), vec![]).expect("second segment"),
        ])
        .expect("composition");
        let observed: Vec<_> = chain
            .segments()
            .iter()
            .map(|segment| (segment.source().start(), segment.source().end()))
            .collect();
        assert_eq!(observed, vec![(0, 3), (7, 10)]);
    }

    #[test]
    fn rejects_ambiguous_or_overlapping_modifications() {
        let first = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        let repeated = first.clone();
        assert_eq!(
            Segment::new(interval(0, 10), vec![first, repeated]),
            Err(SegmentError::ModificationsNotStrictlyOrdered)
        );

        let outside = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    7,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        assert_eq!(
            Segment::new(interval(0, 4), vec![outside]),
            Err(SegmentError::ModificationOutsideSource)
        );

        let deletion = SmallModification::Deletion(DeletionModification::new(
            InternalDeletion::new(interval(2, 5)).expect("deletion"),
            provenance(),
        ));
        let inside_deleted_bases = SmallModification::Insertion(
            InsertionModification::new(
                InsertionAnchor::new(
                    InterbasePosition::new(
                        ReferenceId::parse("ref", "test").expect("reference"),
                        SequenceId::parse("chr1", "test").expect("sequence"),
                        4,
                        10,
                    )
                    .expect("point"),
                )
                .expect("anchor"),
                Arc::<[u8]>::from(&b"A"[..]),
                provenance(),
            )
            .expect("insertion"),
        );
        assert_eq!(
            Segment::new(interval(0, 10), vec![deletion, inside_deleted_bases]),
            Err(SegmentError::OverlappingModifications)
        );

        let other_space = SmallModification::Snp(
            SnpModification::new(
                BasePosition::new(
                    ReferenceId::parse("other-ref", "test").expect("reference"),
                    SequenceId::parse("chr1", "test").expect("sequence"),
                    2,
                    10,
                )
                .expect("position"),
                CanonicalBase::A,
                CanonicalBase::C,
                provenance(),
            )
            .expect("SNP"),
        );
        assert_eq!(
            Segment::new(interval(0, 10), vec![other_space]),
            Err(SegmentError::ModificationSourceMismatch)
        );
    }

    /// **Un tracto de repetición cuenta como base de referencia al anclar.**
    ///
    /// En el recorrido los alelos de repetición se aplican los últimos, de
    /// modo que ningún *indel* llega a ver un tracto sustituido y esta regla
    /// no cambia nada de lo que el programa publica hoy. La función es del
    /// dominio y no del recorrido, sin embargo, y su contrato no debe
    /// depender de ese orden: una inserción detrás de un tracto se ancla al
    /// final del tracto y una deleción que lo borra entero cubre su
    /// intervalo, no el ancla de antes del tracto.
    #[test]
    fn a_realized_repeat_tract_is_reference_for_the_reference_locus() {
        use crate::identity::RepeatLocusId;
        use crate::tandem_repeat::{RepeatMotif, TandemRepeatTract};
        let repeat = SegmentChain::intact(interval(0, 10))
            .with_current_repeat_allele(CurrentRepeatAllelePiece::new(
                RepeatLocusId::parse("STR-1", "test").expect("locus id"),
                interval(2, 6),
                TandemRepeatTract::new(RepeatMotif::new(b"CA").expect("motif"), b"CACACA")
                    .expect("tract"),
                provenance(),
            ))
            .expect("repeat allele");
        // Vigentes: ref[0,2) en [0, 2), el tracto de 6 bases en [2, 8) y
        // ref[6,10) en [8, 12).
        let inserted = repeat
            .with_current_insertion(8, Arc::from(&b"T"[..]), provenance())
            .expect("insertion after the tract");
        assert!(matches!(
            inserted.current_indel_events().last(),
            Some(CurrentIndelEvent::Insertion {
                reference_point: 6,
                ..
            })
        ));
        let deleted = repeat
            .with_current_deletion(2, 8, provenance())
            .expect("deletion of the whole tract");
        assert!(matches!(
            deleted.current_indel_events().last(),
            Some(CurrentIndelEvent::Deletion {
                reference_start: 2,
                reference_end: 6,
                ..
            })
        ));
    }
}
