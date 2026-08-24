//! Construcción y aplicación transaccional de SNP explícitos o de matriz.
//!
//! F2.2.8 no decide genotipos ni fase aleatoria. Recibe el haplotipo destino
//! ya elegido por un consumidor posterior, valida todo el lote contra la
//! referencia autoritativa y solo entonces devuelve el nuevo estado lógico.

use std::collections::BTreeSet;
use std::fmt;

use crate::coordinates::BasePosition;
use crate::haplotype::{HaplotypePhase, Individual};
use crate::identity::{EventId, PopulationId, VariantId, derive_matrix_variant_id};
use crate::materialization::SnpLocus;
use crate::population_snvs::PreparedPopulationSnv;
use crate::reference::Reference;
use crate::segments::{CanonicalBase, ModificationProvenance, SegmentError, SnpModification};

/// Origen inmutable de un SNP aplicado a una copia haplotípica.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SnpOrigin {
    Explicit,
    PopulationMatrix {
        variant_id: VariantId,
        population_id: PopulationId,
        source_id: String,
    },
}

/// SNP ya expresado en el espacio interno de una referencia y listo para aplicar.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ApplicableSnp {
    locus: SnpLocus,
    event_id: EventId,
    origin: SnpOrigin,
}

impl ApplicableSnp {
    #[must_use]
    pub const fn explicit(locus: SnpLocus, event_id: EventId) -> Self {
        Self {
            locus,
            event_id,
            origin: SnpOrigin::Explicit,
        }
    }

    /// Convierte una fila preparada de la matriz canónica sin volver a sortear
    /// locus, alelo o frecuencia.
    pub fn from_population_matrix(
        reference: &Reference,
        row: &PreparedPopulationSnv,
        event_id: EventId,
    ) -> Result<Self, SnpApplicationError> {
        let expected_id = derive_matrix_variant_id(
            reference.id(),
            row.contig(),
            row.position_1based(),
            row.reference(),
            row.alternate(),
        )
        .map_err(|error| SnpApplicationError::InvalidMatrixVariantId {
            detail: error.render_human(),
        })?;
        if &expected_id != row.variant_id() {
            return Err(SnpApplicationError::InvalidMatrixVariantId {
                detail: format!(
                    "la fila declara {} y la identidad canónica es {}",
                    row.variant_id(),
                    expected_id
                ),
            });
        }
        let sequence = reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == row.contig())
            .ok_or_else(|| SnpApplicationError::MissingReferenceSequence {
                sequence_id: row.contig().clone(),
            })?;
        let offset = row.position_1based().checked_sub(1).ok_or(
            SnpApplicationError::InvalidMatrixPosition {
                position_1based: row.position_1based(),
            },
        )?;
        let position = BasePosition::new(
            reference.id().clone(),
            row.contig().clone(),
            offset,
            sequence.sequence().len(),
        )
        .map_err(|_| SnpApplicationError::InvalidMatrixPosition {
            position_1based: row.position_1based(),
        })?;
        let reference_base = CanonicalBase::from_byte(row.reference()).map_err(|_| {
            SnpApplicationError::NonCanonicalMatrixAllele {
                allele: row.reference(),
            }
        })?;
        let alternate = CanonicalBase::from_byte(row.alternate()).map_err(|_| {
            SnpApplicationError::NonCanonicalMatrixAllele {
                allele: row.alternate(),
            }
        })?;
        let locus = SnpLocus::new(position, reference_base, alternate)
            .map_err(|_| SnpApplicationError::EqualAlleles)?;
        Ok(Self {
            locus,
            event_id,
            origin: SnpOrigin::PopulationMatrix {
                variant_id: row.variant_id().clone(),
                population_id: row.population_id().clone(),
                source_id: matrix_source_id(row.source_id()).to_owned(),
            },
        })
    }

    #[must_use]
    pub const fn locus(&self) -> &SnpLocus {
        &self.locus
    }

    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub const fn origin(&self) -> &SnpOrigin {
        &self.origin
    }

    fn modification(&self) -> Result<SnpModification, SnpApplicationError> {
        let provenance = ModificationProvenance::new(self.event_id.clone());
        match &self.origin {
            SnpOrigin::Explicit => SnpModification::new(
                self.locus.position().clone(),
                self.locus.reference(),
                self.locus.alternate(),
                provenance,
            )
            .map_err(SnpApplicationError::Segment),
            SnpOrigin::PopulationMatrix {
                variant_id,
                population_id,
                source_id,
            } => SnpModification::from_population_matrix(
                self.locus.position().clone(),
                self.locus.reference(),
                self.locus.alternate(),
                provenance,
                variant_id.clone(),
                population_id.clone(),
                source_id.clone(),
            )
            .map_err(SnpApplicationError::Segment),
        }
    }
}

/// La columna TSV `source_id` puede quedar vacía. El estado de dominio no
/// propaga ese valor no identificable a `truth.jsonl`: la matriz canónica es
/// entonces la fuente publicada común de esas filas.
fn matrix_source_id(source_id: &str) -> &str {
    if source_id.is_empty() {
        "MATRIX_DEFAULT"
    } else {
        source_id
    }
}

/// Solicitud de aplicación sobre una fase explícita, sin RNG implícito.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TargetedSnp {
    phase: HaplotypePhase,
    snp: ApplicableSnp,
}

impl TargetedSnp {
    #[must_use]
    pub const fn new(phase: HaplotypePhase, snp: ApplicableSnp) -> Self {
        Self { phase, snp }
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn snp(&self) -> &ApplicableSnp {
        &self.snp
    }
}

/// Rechazo de construcción, preflight o aplicación de SNP.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum SnpApplicationError {
    InvalidMatrixVariantId {
        detail: String,
    },
    InvalidMatrixPosition {
        position_1based: u64,
    },
    NonCanonicalMatrixAllele {
        allele: u8,
    },
    EqualAlleles,
    IndividualReferenceMismatch,
    LocusReferenceMismatch,
    MissingReferenceSequence {
        sequence_id: crate::identity::SequenceId,
    },
    ReferenceBaseUnknown {
        position: u64,
    },
    ReferenceBaseMismatch {
        position: u64,
        expected: CanonicalBase,
        found: u8,
    },
    MissingTargetCopy {
        phase: HaplotypePhase,
    },
    DuplicateTarget {
        phase: HaplotypePhase,
        position: u64,
    },
    Segment(SegmentError),
}

impl fmt::Display for SnpApplicationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for SnpApplicationError {}

/// Valida el lote completo sin cambiar el individuo. Es obligatorio llamarlo
/// antes de publicar cualquier estado derivado; [`apply_snps`] lo invoca de
/// nuevo para conservar la operación transaccional.
pub fn preflight_snps(
    individual: &Individual,
    reference: &Reference,
    applications: &[TargetedSnp],
) -> Result<(), SnpApplicationError> {
    if individual.reference_id() != reference.id() {
        return Err(SnpApplicationError::IndividualReferenceMismatch);
    }
    let mut targets = BTreeSet::new();
    for application in applications {
        let locus = application.snp().locus();
        if locus.position().reference_id() != reference.id() {
            return Err(SnpApplicationError::LocusReferenceMismatch);
        }
        let target = (
            locus.position().sequence_id().clone(),
            application.phase(),
            locus.position().offset(),
        );
        if !targets.insert(target) {
            return Err(SnpApplicationError::DuplicateTarget {
                phase: application.phase(),
                position: locus.position().offset(),
            });
        }
        validate_reference_base(reference, locus)?;
        let copy = individual
            .copies_for(locus.position().sequence_id())
            .and_then(|copies| {
                copies
                    .iter()
                    .find(|copy| copy.phase() == application.phase())
            })
            .ok_or(SnpApplicationError::MissingTargetCopy {
                phase: application.phase(),
            })?;
        // Construir una cadena candidata detecta SNP previos, deleciones y
        // segmentos no cubiertos sin alterar la cadena compartida.
        copy.segments()
            .with_snp(application.snp().modification()?)
            .map_err(SnpApplicationError::Segment)?;
    }
    Ok(())
}

/// Aplica un lote validado de SNP y devuelve un individuo nuevo. No consume
/// aleatoriedad, no modifica la matriz y no deja estados parciales ante error.
pub fn apply_snps(
    individual: &Individual,
    reference: &Reference,
    applications: &[TargetedSnp],
) -> Result<Individual, SnpApplicationError> {
    preflight_snps(individual, reference, applications)?;
    let mut updated = individual.clone();
    for application in applications {
        let locus = application.snp().locus();
        let copy = updated
            .copies_for(locus.position().sequence_id())
            .and_then(|copies| {
                copies
                    .iter()
                    .find(|copy| copy.phase() == application.phase())
            })
            .ok_or(SnpApplicationError::MissingTargetCopy {
                phase: application.phase(),
            })?;
        let segments = copy
            .segments()
            .with_snp(application.snp().modification()?)
            .map_err(SnpApplicationError::Segment)?;
        updated = updated
            .with_copy_segments(
                locus.position().sequence_id(),
                application.phase(),
                segments,
            )
            .map_err(|_| SnpApplicationError::MissingTargetCopy {
                phase: application.phase(),
            })?;
    }
    Ok(updated)
}

fn validate_reference_base(
    reference: &Reference,
    locus: &SnpLocus,
) -> Result<(), SnpApplicationError> {
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id() == locus.position().sequence_id())
        .ok_or_else(|| SnpApplicationError::MissingReferenceSequence {
            sequence_id: locus.position().sequence_id().clone(),
        })?;
    let offset = usize::try_from(locus.position().offset()).map_err(|_| {
        SnpApplicationError::ReferenceBaseMismatch {
            position: locus.position().offset(),
            expected: locus.reference(),
            found: b'?',
        }
    })?;
    let found = sequence.sequence().as_bytes().get(offset).copied().ok_or(
        SnpApplicationError::ReferenceBaseMismatch {
            position: locus.position().offset(),
            expected: locus.reference(),
            found: b'?',
        },
    )?;
    if found == b'N' {
        return Err(SnpApplicationError::ReferenceBaseUnknown {
            position: locus.position().offset(),
        });
    }
    if found != locus.reference().as_byte() {
        return Err(SnpApplicationError::ReferenceBaseMismatch {
            position: locus.position().offset(),
            expected: locus.reference(),
            found,
        });
    }
    Ok(())
}
