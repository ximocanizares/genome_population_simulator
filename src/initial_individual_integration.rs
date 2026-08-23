//! Integración transaccional de genotipos faseados y cambios directos F2.3.4.
//!
//! Esta frontera recibe una copia lógica ya creada y un único conjunto de loci
//! faseados para su identidad. Proyecta solamente los ALT sorteados de la
//! matriz sobre sus fases, conserva los cambios directos e invoca el compositor
//! ya existente para los indels. No materializa bases ni calcula estadísticas.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::application::EffectiveInitialGenotypeModel;
use crate::haplotype::{HaplotypePhase, Individual};
use crate::identity::{AutomaticIdRequest, EventId, IdClass, IdentityRegistry, VariantId};
use crate::initial_genotypes::{InitialPhasedGenotypes, InitialPhasedLocus};
use crate::population_snvs::{
    PreparedPopulationSnv, PreparedPopulationSnvs, PreparedReferenceIdentity,
};
use crate::randomness::RootSeed;
use crate::reference::Reference;
use crate::small_event_composition::{
    DeclaredSmallEvent, RandomSmallEventBlock, SmallEventCompositionError,
    apply_declared_small_events, preflight_declared_small_events,
};
use crate::snp_application::{ApplicableSnp, SnpApplicationError, SnpOrigin, TargetedSnp};

/// Confirmación específica para conservar una restricción SNP incompatible con
/// el genotipo sorteado. No equivale a una opción general de forzado.
#[derive(Clone, Copy, Debug, Default, Eq, PartialEq)]
pub struct IncompatibleGenotypeAuthorization {
    allowed: bool,
}

impl IncompatibleGenotypeAuthorization {
    #[must_use]
    pub const fn denied() -> Self {
        Self { allowed: false }
    }

    #[must_use]
    pub const fn allowed() -> Self {
        Self { allowed: true }
    }

    #[must_use]
    pub const fn is_allowed(self) -> bool {
        self.allowed
    }
}

/// Resultado observable de una restricción SNP explícita que coincide con un
/// locus de la matriz. Conserva que el estado final procede de la restricción,
/// sin cambiar la matriz ni corregir el genotipo sorteado.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExplicitSnpConstraint {
    event_id: EventId,
    phase: HaplotypePhase,
    matrix_variant_id: VariantId,
    sampled_alternate: Option<u8>,
    explicit_alternate: u8,
    incompatible: bool,
}

impl ExplicitSnpConstraint {
    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn matrix_variant_id(&self) -> &VariantId {
        &self.matrix_variant_id
    }

    #[must_use]
    pub const fn sampled_alternate(&self) -> Option<u8> {
        self.sampled_alternate
    }

    #[must_use]
    pub const fn explicit_alternate(&self) -> u8 {
        self.explicit_alternate
    }

    #[must_use]
    pub const fn is_incompatible(&self) -> bool {
        self.incompatible
    }
}

/// Estado integrado de una copia individual, todavía compacto y no
/// materializado, junto con las restricciones explícitas que reemplazaron la
/// proyección de matriz en su mismo locus/fase.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialIndividualIntegration {
    individual: Individual,
    explicit_snp_constraints: Vec<ExplicitSnpConstraint>,
}

/// Entradas compartidas para integrar una sola copia individual. Los cambios
/// siguen prestados para evitar duplicarlos al crear la frontera transaccional.
#[derive(Clone, Copy, Debug)]
pub struct InitialIndividualIntegrationRequest<'a> {
    matrix: &'a PreparedPopulationSnvs,
    phased: &'a InitialPhasedGenotypes,
    root_seed: RootSeed,
    explicit_events: &'a [DeclaredSmallEvent],
    random_blocks: &'a [RandomSmallEventBlock],
    authorization: IncompatibleGenotypeAuthorization,
    genotype_model: EffectiveInitialGenotypeModel,
}

impl<'a> InitialIndividualIntegrationRequest<'a> {
    #[must_use]
    pub const fn new(
        matrix: &'a PreparedPopulationSnvs,
        phased: &'a InitialPhasedGenotypes,
        root_seed: RootSeed,
        explicit_events: &'a [DeclaredSmallEvent],
        random_blocks: &'a [RandomSmallEventBlock],
        authorization: IncompatibleGenotypeAuthorization,
        genotype_model: EffectiveInitialGenotypeModel,
    ) -> Self {
        Self {
            matrix,
            phased,
            root_seed,
            explicit_events,
            random_blocks,
            authorization,
            genotype_model,
        }
    }
}

impl InitialIndividualIntegration {
    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub fn explicit_snp_constraints(&self) -> &[ExplicitSnpConstraint] {
        &self.explicit_snp_constraints
    }
}

/// Rechazo previo a cualquier cambio de la copia individual.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum InitialIndividualIntegrationError {
    IndividualReferenceMismatch,
    MatrixReferenceMismatch {
        matrix_reference: PreparedReferenceIdentity,
        supplied_reference: PreparedReferenceIdentity,
    },
    PhasedMatrixMismatch {
        phased_reference: PreparedReferenceIdentity,
        phased_fingerprint: String,
        supplied_fingerprint: String,
    },
    PhasedLocusReferenceMismatch {
        reference_id: String,
    },
    PhasedLocusIndividualMismatch {
        individual_id: String,
    },
    InvalidPhasedLocus {
        contig: String,
        position_1based: u64,
    },
    MissingMatrixAllele {
        variant_id: VariantId,
        population_id: String,
    },
    InconsistentMatrixAllele {
        variant_id: VariantId,
    },
    IncompatibleExplicitSnp {
        event_id: EventId,
        phase: HaplotypePhase,
        matrix_variant_id: VariantId,
        sampled_alternate: Option<u8>,
        explicit_alternate: u8,
    },
    ExplicitSnpAlleleOutsideMatrix {
        event_id: EventId,
        phase: HaplotypePhase,
        explicit_alternate: u8,
    },
    NonExplicitDirectSnp {
        event_id: EventId,
    },
    MatrixSnp(SnpApplicationError),
    EventComposition(SmallEventCompositionError),
    Identity(crate::diagnostics::Diagnostic),
}

impl fmt::Display for InitialIndividualIntegrationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for InitialIndividualIntegrationError {}

/// Integra los SNP ALT de la matriz en un individuo ya faseado y después los
/// cambios directos declarados. Todos los conflictos se validan sobre estados
/// candidatos antes de devolver una copia derivada.
pub fn integrate_initial_individual(
    individual: &Individual,
    reference: &Reference,
    request: InitialIndividualIntegrationRequest<'_>,
) -> Result<InitialIndividualIntegration, InitialIndividualIntegrationError> {
    preflight_references(individual, reference, request.matrix, request.phased)?;
    let matrix_projection = project_matrix_snps(
        individual,
        reference,
        request.matrix,
        request.phased,
        request.root_seed,
    )?;
    let (matrix_snps, constraints) = resolve_explicit_snp_constraints(
        matrix_projection,
        request.explicit_events,
        request.authorization,
        request.genotype_model,
    )?;

    let mut events = Vec::with_capacity(matrix_snps.len() + request.explicit_events.len());
    events.extend(matrix_snps.into_iter().map(DeclaredSmallEvent::Snp));
    events.extend(request.explicit_events.iter().cloned());
    preflight_declared_small_events(
        individual,
        reference,
        request.root_seed,
        &events,
        request.random_blocks,
    )
    .map_err(InitialIndividualIntegrationError::EventComposition)?;
    let result = apply_declared_small_events(
        individual,
        reference,
        request.root_seed,
        &events,
        request.random_blocks,
    )
    .map_err(InitialIndividualIntegrationError::EventComposition)?;
    Ok(InitialIndividualIntegration {
        individual: result.individual().clone(),
        explicit_snp_constraints: constraints,
    })
}

fn preflight_references(
    individual: &Individual,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
    phased: &InitialPhasedGenotypes,
) -> Result<(), InitialIndividualIntegrationError> {
    if individual.reference_id() != reference.id() {
        return Err(InitialIndividualIntegrationError::IndividualReferenceMismatch);
    }
    let supplied_reference = PreparedReferenceIdentity::from_reference(reference);
    if matrix.prepared_reference() != &supplied_reference {
        return Err(InitialIndividualIntegrationError::MatrixReferenceMismatch {
            matrix_reference: matrix.prepared_reference().clone(),
            supplied_reference,
        });
    }
    if phased.matrix_reference() != matrix.prepared_reference()
        || phased.matrix_fingerprint() != matrix.canonical_fingerprint()
    {
        return Err(InitialIndividualIntegrationError::PhasedMatrixMismatch {
            phased_reference: phased.matrix_reference().clone(),
            phased_fingerprint: phased.matrix_fingerprint().to_owned(),
            supplied_fingerprint: matrix.canonical_fingerprint(),
        });
    }
    Ok(())
}

#[derive(Clone)]
struct MatrixProjection {
    snps: Vec<TargetedSnp>,
    expected_by_target: BTreeMap<SnpTarget, ExpectedMatrixAllele>,
}

#[derive(Clone)]
struct ExpectedMatrixAllele {
    variant_ids_by_alternate: BTreeMap<u8, VariantId>,
    frequencies_by_alternate: BTreeMap<u8, f64>,
    reference_frequency: f64,
    sampled_alternate: Option<u8>,
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct SnpTarget {
    sequence_id: crate::identity::SequenceId,
    phase: HaplotypePhase,
    offset: u64,
}

fn project_matrix_snps(
    individual: &Individual,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
    phased: &InitialPhasedGenotypes,
    root_seed: RootSeed,
) -> Result<MatrixProjection, InitialIndividualIntegrationError> {
    let mut snps = Vec::new();
    let mut expected_by_target = BTreeMap::new();
    for locus in phased.loci() {
        validate_phased_locus(individual, reference, locus)?;
        let rows = rows_for_locus(matrix, locus)?;
        for allele in locus.alleles_by_phase() {
            let target = SnpTarget {
                sequence_id: locus.composition().contig().clone(),
                phase: allele.phase(),
                offset: locus.composition().position_1based() - 1,
            };
            let sampled_alternate = if allele.allele_index() == 0 {
                None
            } else {
                let alternate = locus
                    .composition()
                    .alternate_dosages()
                    .get(usize::from(allele.allele_index() - 1))
                    .ok_or_else(|| invalid_locus(locus))?;
                Some(alternate.alternate())
            };
            let variant_ids_by_alternate = rows
                .iter()
                .map(|row| (row.alternate(), row.variant_id().clone()))
                .collect();
            let frequencies_by_alternate = rows
                .iter()
                .map(|row| (row.alternate(), row.alt_af()))
                .collect::<BTreeMap<_, _>>();
            let reference_frequency = 1.0 - frequencies_by_alternate.values().sum::<f64>();
            if expected_by_target
                .insert(
                    target.clone(),
                    ExpectedMatrixAllele {
                        variant_ids_by_alternate,
                        frequencies_by_alternate,
                        reference_frequency,
                        sampled_alternate,
                    },
                )
                .is_some()
            {
                return Err(invalid_locus(locus));
            }
            if let Some(alternate) = sampled_alternate {
                let row = rows
                    .iter()
                    .find(|row| row.alternate() == alternate)
                    .ok_or_else(|| invalid_locus(locus))?;
                let event_id = derive_matrix_event_id(root_seed, individual, row, allele.phase())?;
                snps.push(TargetedSnp::new(
                    allele.phase(),
                    ApplicableSnp::from_population_matrix(reference, row, event_id)
                        .map_err(InitialIndividualIntegrationError::MatrixSnp)?,
                ));
            }
        }
    }
    snps.sort_by_key(snp_target);
    Ok(MatrixProjection {
        snps,
        expected_by_target,
    })
}

fn validate_phased_locus(
    individual: &Individual,
    reference: &Reference,
    locus: &InitialPhasedLocus,
) -> Result<(), InitialIndividualIntegrationError> {
    let composition = locus.composition();
    if composition.reference_id() != reference.id().as_str() {
        return Err(
            InitialIndividualIntegrationError::PhasedLocusReferenceMismatch {
                reference_id: composition.reference_id().to_owned(),
            },
        );
    }
    if composition.individual_id() != individual.id() {
        return Err(
            InitialIndividualIntegrationError::PhasedLocusIndividualMismatch {
                individual_id: composition.individual_id().as_str().to_owned(),
            },
        );
    }
    let expected_phases = (0..individual.ploidy().copy_count())
        .map(HaplotypePhase::from_copy_index)
        .collect::<Result<BTreeSet<_>, _>>()
        .map_err(|_| invalid_locus(locus))?;
    let observed_phases = locus
        .alleles_by_phase()
        .iter()
        .map(|allele| allele.phase())
        .collect::<BTreeSet<_>>();
    if observed_phases != expected_phases
        || locus.alleles_by_phase().len() != expected_phases.len()
        || locus
            .alleles_by_phase()
            .iter()
            .any(|allele| allele.allele_index() as usize > composition.alternate_dosages().len())
    {
        return Err(invalid_locus(locus));
    }
    Ok(())
}

fn rows_for_locus<'a>(
    matrix: &'a PreparedPopulationSnvs,
    locus: &InitialPhasedLocus,
) -> Result<Vec<&'a PreparedPopulationSnv>, InitialIndividualIntegrationError> {
    let composition = locus.composition();
    let mut rows = Vec::with_capacity(composition.alternate_dosages().len());
    for dosage in composition.alternate_dosages() {
        let row = matrix
            .rows()
            .iter()
            .find(|row| {
                row.variant_id() == dosage.variant_id()
                    && row.population_id() == composition.population_id()
            })
            .ok_or_else(|| InitialIndividualIntegrationError::MissingMatrixAllele {
                variant_id: dosage.variant_id().clone(),
                population_id: composition.population_id().as_str().to_owned(),
            })?;
        if row.contig() != composition.contig()
            || row.position_1based() != composition.position_1based()
            || row.reference() != composition.reference()
            || row.alternate() != dosage.alternate()
        {
            return Err(
                InitialIndividualIntegrationError::InconsistentMatrixAllele {
                    variant_id: dosage.variant_id().clone(),
                },
            );
        }
        rows.push(row);
    }
    Ok(rows)
}

fn resolve_explicit_snp_constraints(
    projection: MatrixProjection,
    explicit_events: &[DeclaredSmallEvent],
    authorization: IncompatibleGenotypeAuthorization,
    genotype_model: EffectiveInitialGenotypeModel,
) -> Result<(Vec<TargetedSnp>, Vec<ExplicitSnpConstraint>), InitialIndividualIntegrationError> {
    let mut matrix_snps = projection.snps;
    let mut suppressed_targets = BTreeSet::new();
    let mut constraints = Vec::new();
    let explicit_by_target = explicit_events
        .iter()
        .filter_map(|event| match event {
            DeclaredSmallEvent::Snp(snp) => Some((snp_target(snp), snp)),
            _ => None,
        })
        .collect::<BTreeMap<_, _>>();
    for event in explicit_events {
        let DeclaredSmallEvent::Snp(explicit) = event else {
            continue;
        };
        if !matches!(explicit.snp().origin(), SnpOrigin::Explicit) {
            return Err(InitialIndividualIntegrationError::NonExplicitDirectSnp {
                event_id: explicit.snp().event_id().clone(),
            });
        }
        let target = snp_target(explicit);
        let Some(expected) = projection.expected_by_target.get(&target) else {
            continue;
        };
        let explicit_alternate = explicit.snp().locus().alternate().as_byte();
        let Some(matrix_variant_id) = expected
            .variant_ids_by_alternate
            .get(&explicit_alternate)
            .cloned()
        else {
            // La autorización solo conserva una realización forzada de un
            // alelo ya declarado por la matriz; nunca relaja el catálogo de
            // alelos de su locus.
            return Err(
                InitialIndividualIntegrationError::ExplicitSnpAlleleOutsideMatrix {
                    event_id: explicit.snp().event_id().clone(),
                    phase: explicit.phase(),
                    explicit_alternate,
                },
            );
        };
        let incompatible = !genotype_has_positive_probability(
            &projection.expected_by_target,
            &target,
            &explicit_by_target,
            genotype_model,
        );
        if incompatible && !authorization.is_allowed() {
            return Err(InitialIndividualIntegrationError::IncompatibleExplicitSnp {
                event_id: explicit.snp().event_id().clone(),
                phase: explicit.phase(),
                matrix_variant_id,
                sampled_alternate: expected.sampled_alternate,
                explicit_alternate,
            });
        }
        suppressed_targets.insert(target);
        constraints.push(ExplicitSnpConstraint {
            event_id: explicit.snp().event_id().clone(),
            phase: explicit.phase(),
            matrix_variant_id,
            sampled_alternate: expected.sampled_alternate,
            explicit_alternate,
            incompatible,
        });
    }
    matrix_snps.retain(|snp| !suppressed_targets.contains(&snp_target(snp)));
    Ok((matrix_snps, constraints))
}

fn genotype_has_positive_probability(
    expected_by_target: &BTreeMap<SnpTarget, ExpectedMatrixAllele>,
    target: &SnpTarget,
    explicit_by_target: &BTreeMap<SnpTarget, &TargetedSnp>,
    model: EffectiveInitialGenotypeModel,
) -> bool {
    let Some(expected) = expected_by_target.get(target) else {
        return false;
    };
    let alleles = expected_by_target
        .iter()
        .filter(|(candidate, _)| {
            candidate.sequence_id == target.sequence_id && candidate.offset == target.offset
        })
        .map(|(candidate, sampled)| {
            explicit_by_target
                .get(candidate)
                .map(|snp| snp.snp().locus().alternate().as_byte())
                .or(sampled.sampled_alternate)
        })
        .collect::<Vec<_>>();
    let positive = |allele: Option<u8>| match allele {
        Some(value) => expected
            .frequencies_by_alternate
            .get(&value)
            .is_some_and(|frequency| *frequency > 0.0),
        None => expected.reference_frequency > 0.0,
    };
    if !alleles.iter().copied().all(positive) {
        return false;
    }
    match model {
        EffectiveInitialGenotypeModel::HardyWeinberg => true,
        EffectiveInitialGenotypeModel::FDiploid {
            inbreeding_coefficient,
        } => inbreeding_coefficient < 1.0 || alleles.windows(2).all(|pair| pair[0] == pair[1]),
    }
}

fn snp_target(snp: &TargetedSnp) -> SnpTarget {
    SnpTarget {
        sequence_id: snp.snp().locus().position().sequence_id().clone(),
        phase: snp.phase(),
        offset: snp.snp().locus().position().offset(),
    }
}

fn derive_matrix_event_id(
    root_seed: RootSeed,
    individual: &Individual,
    row: &PreparedPopulationSnv,
    phase: HaplotypePhase,
) -> Result<EventId, InitialIndividualIntegrationError> {
    let event_key = format!(
        "reference={};individual={};population={};variant={};phase={}",
        individual.reference_id().as_str(),
        individual.id().as_str(),
        row.population_id().as_str(),
        row.variant_id().as_str(),
        phase.copy_index(),
    );
    let request = AutomaticIdRequest::new(
        IdClass::Event,
        "f2.3.4-matrix-snp",
        &[("event_key", event_key.as_str())],
    )
    .map_err(InitialIndividualIntegrationError::Identity)?;
    IdentityRegistry::default()
        .derive_event(root_seed, &request)
        .map_err(InitialIndividualIntegrationError::Identity)
}

fn invalid_locus(locus: &InitialPhasedLocus) -> InitialIndividualIntegrationError {
    InitialIndividualIntegrationError::InvalidPhasedLocus {
        contig: locus.composition().contig().as_str().to_owned(),
        position_1based: locus.composition().position_1based(),
    }
}
