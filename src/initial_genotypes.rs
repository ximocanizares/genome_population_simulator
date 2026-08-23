//! Composición y fase genotípicas iniciales, antes de materializar copias.
//!
//! F2.3.2 consume una matriz ya preparada, por lo que no reinterpreta `alt_af`,
//! no elige loci y no modifica la matriz. Antes de derivar RNG exige que la
//! referencia declarada coincida con la identidad y huella canónica retenidas
//! durante la preparación. F2.3.3 asigna después las dosis ALT a copias con un
//! ámbito independiente por locus, sin aplicar todavía cambios a las copias.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::application::{
    EffectiveConfig, EffectiveInitialGenotypeModel, EffectiveInitialPopulation,
};
use crate::haplotype::HaplotypePhase;
use crate::identity::{IndividualId, PopulationId, SequenceId, VariantId};
use crate::population_snvs::{
    PreparedPopulationSnv, PreparedPopulationSnvs, PreparedReferenceIdentity,
};
use crate::randomness::{
    F_DIPLOID_GENOTYPE_ALGORITHM, HWE_COMPOSITION_ALGORITHM, PHASE_PERMUTATION_ALGORITHM,
    RandomSamplingError, RootSeed, RuntimeSubseedScope, derive_subseed, f_diploid_genotype,
    hwe_composition, phase_permutation,
};
use crate::reference::Reference;

/// La composición HWE conserva una dosis por alelo, sin fase.
pub const INITIAL_GENOTYPE_HWE_ALGORITHM: &str = HWE_COMPOSITION_ALGORITHM;
/// La composición F conserva una dosis por alelo, sin fase.
pub const INITIAL_GENOTYPE_F_DIPLOID_ALGORITHM: &str = F_DIPLOID_GENOTYPE_ALGORITHM;
/// Ámbito de subsemilla para la composición genotípica inicial.
pub const INITIAL_GENOTYPE_SCOPE: &str = "initial_genotype";
/// Ámbito de subsemilla para asignar fases de un locus tras fijar su composición.
pub const INITIAL_PHASE_SCOPE: &str = "initial_phase_locus";
/// Algoritmo que permuta copias sin alterar la dosis ya muestreada.
pub const INITIAL_PHASE_ALGORITHM: &str = PHASE_PERMUTATION_ALGORITHM;

/// Alelo ALT y su dosis no faseada en una composición inicial.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialAlternateDosage {
    variant_id: VariantId,
    alternate: u8,
    dosage: u8,
}

impl InitialAlternateDosage {
    #[must_use]
    pub const fn variant_id(&self) -> &VariantId {
        &self.variant_id
    }

    #[must_use]
    pub const fn alternate(&self) -> u8 {
        self.alternate
    }

    #[must_use]
    pub const fn dosage(&self) -> u8 {
        self.dosage
    }
}

/// Composición no faseada de un locus en un individuo.
///
/// El orden de ALT es ASCII ascendente A/C/G/T y forma parte de la clave RNG
/// multialélica. La dosis REF más las ALT suma siempre la ploidía de la cohorte.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialLocusComposition {
    reference_id: String,
    population_id: PopulationId,
    individual_id: IndividualId,
    contig: SequenceId,
    position_1based: u64,
    reference: u8,
    reference_dosage: u8,
    alternate_dosages: Vec<InitialAlternateDosage>,
}

impl InitialLocusComposition {
    #[must_use]
    pub fn reference_id(&self) -> &str {
        &self.reference_id
    }

    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn contig(&self) -> &SequenceId {
        &self.contig
    }

    #[must_use]
    pub const fn position_1based(&self) -> u64 {
        self.position_1based
    }

    #[must_use]
    pub const fn reference(&self) -> u8 {
        self.reference
    }

    #[must_use]
    pub const fn reference_dosage(&self) -> u8 {
        self.reference_dosage
    }

    #[must_use]
    pub fn alternate_dosages(&self) -> &[InitialAlternateDosage] {
        &self.alternate_dosages
    }
}

/// Resultado intermedio de F2.3.2, ordenado e independiente de la fase.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialGenotypeCompositions {
    loci: Vec<InitialLocusComposition>,
    matrix_reference: PreparedReferenceIdentity,
    matrix_fingerprint: String,
}

impl InitialGenotypeCompositions {
    #[must_use]
    pub fn loci(&self) -> &[InitialLocusComposition] {
        &self.loci
    }

    #[must_use]
    pub const fn matrix_reference(&self) -> &PreparedReferenceIdentity {
        &self.matrix_reference
    }

    #[must_use]
    pub fn matrix_fingerprint(&self) -> &str {
        &self.matrix_fingerprint
    }
}

/// Alelo canónico asignado a una copia faseada de un locus.
///
/// `allele_index=0` representa REF; `1..` indexan `alternate_dosages()` de la
/// composición retenida. Cada fase aparece exactamente una vez.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialPhasedAllele {
    phase: HaplotypePhase,
    allele_index: u8,
}

impl InitialPhasedAllele {
    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub const fn allele_index(&self) -> u8 {
        self.allele_index
    }
}

/// Composición de un locus junto con todos sus alelos asignados a fases.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialPhasedLocus {
    composition: InitialLocusComposition,
    alleles_by_phase: Vec<InitialPhasedAllele>,
}

impl InitialPhasedLocus {
    #[must_use]
    pub const fn composition(&self) -> &InitialLocusComposition {
        &self.composition
    }

    #[must_use]
    pub fn alleles_by_phase(&self) -> &[InitialPhasedAllele] {
        &self.alleles_by_phase
    }
}

/// Resultado F2.3.3, aún separado del estado de copias haplotípicas.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct InitialPhasedGenotypes {
    loci: Vec<InitialPhasedLocus>,
    matrix_reference: PreparedReferenceIdentity,
    matrix_fingerprint: String,
}

impl InitialPhasedGenotypes {
    #[must_use]
    pub fn loci(&self) -> &[InitialPhasedLocus] {
        &self.loci
    }

    /// Vista aislada de una identidad; nunca mezcla loci de otro individuo.
    #[must_use]
    pub fn for_individual(&self, individual_id: &IndividualId) -> Self {
        Self {
            loci: self
                .loci
                .iter()
                .filter(|locus| locus.composition().individual_id() == individual_id)
                .cloned()
                .collect(),
            matrix_reference: self.matrix_reference.clone(),
            matrix_fingerprint: self.matrix_fingerprint.clone(),
        }
    }

    #[must_use]
    pub const fn matrix_reference(&self) -> &PreparedReferenceIdentity {
        &self.matrix_reference
    }

    #[must_use]
    pub fn matrix_fingerprint(&self) -> &str {
        &self.matrix_fingerprint
    }
}

/// Error anterior al sorteo o procedente de la utilidad RNG normativa.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum InitialGenotypeSamplingError {
    ReferenceMismatch {
        configured_reference_id: String,
        sampled_reference: PreparedReferenceIdentity,
        prepared_reference: PreparedReferenceIdentity,
    },
    MissingMatrixPopulation {
        population_id: String,
    },
    InvalidLocusFrequencies {
        contig: String,
        position_1based: u64,
    },
    InconsistentLocus {
        contig: String,
        position_1based: u64,
    },
    Random(RandomSamplingError),
    Derivation(crate::diagnostics::Diagnostic),
}

impl fmt::Display for InitialGenotypeSamplingError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ReferenceMismatch {
                configured_reference_id,
                sampled_reference,
                prepared_reference,
            } => write!(
                formatter,
                "la configuración usa la referencia {configured_reference_id}, pero el muestreo ({}, {}) no coincide con la matriz preparada ({}, {})",
                sampled_reference.reference_id(),
                sampled_reference.canonical_fingerprint().as_str(),
                prepared_reference.reference_id(),
                prepared_reference.canonical_fingerprint().as_str(),
            ),
            Self::MissingMatrixPopulation { population_id } => {
                write!(
                    formatter,
                    "la matriz no contiene la población inicial {population_id}"
                )
            }
            Self::InvalidLocusFrequencies {
                contig,
                position_1based,
            } => write!(
                formatter,
                "frecuencias inválidas en {contig}:{position_1based}"
            ),
            Self::InconsistentLocus {
                contig,
                position_1based,
            } => write!(
                formatter,
                "locus inconsistente en {contig}:{position_1based}"
            ),
            Self::Random(error) => write!(formatter, "error RNG al muestrear genotipo: {error:?}"),
            Self::Derivation(error) => {
                write!(
                    formatter,
                    "no se pudo derivar la subsemilla: {}",
                    error.render_human()
                )
            }
        }
    }
}

impl std::error::Error for InitialGenotypeSamplingError {}

/// Error de F2.3.3 al derivar o aplicar la fase a una composición ya válida.
#[derive(Clone, Debug, Eq, PartialEq)]
pub enum InitialPhaseAssignmentError {
    InvalidComposition {
        contig: String,
        position_1based: u64,
    },
    Random(RandomSamplingError),
    Derivation(crate::diagnostics::Diagnostic),
}

impl fmt::Display for InitialPhaseAssignmentError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidComposition {
                contig,
                position_1based,
            } => write!(
                formatter,
                "la composición no tiene una ploidía válida en {contig}:{position_1based}"
            ),
            Self::Random(error) => write!(formatter, "error RNG al asignar fase: {error:?}"),
            Self::Derivation(error) => write!(
                formatter,
                "no se pudo derivar la subsemilla de fase: {}",
                error.render_human()
            ),
        }
    }
}

impl std::error::Error for InitialPhaseAssignmentError {}

/// Muestrea las composiciones declaradas en una configuración efectiva.
///
/// Esta frontera no publica estadísticas ni asigna fase; ambos pertenecen a
/// F2.3.5 y F2.3.3, respectivamente.
pub fn sample_initial_genotypes(
    config: &EffectiveConfig,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
) -> Result<InitialGenotypeCompositions, InitialGenotypeSamplingError> {
    sample_initial_genotype_compositions(
        config.reference_id(),
        reference,
        config.reproducibility().root_seed(),
        config.initial_populations(),
        matrix,
    )
}

/// Asigna la composición completa de cada locus a copias faseadas independientes.
///
/// Esta frontera no materializa ni modifica copias y, por tanto, no adelanta
/// las incompatibilidades ni la autorización exclusiva de F2.3.4.
pub fn assign_initial_phases(
    root_seed: RootSeed,
    compositions: &InitialGenotypeCompositions,
) -> Result<InitialPhasedGenotypes, InitialPhaseAssignmentError> {
    let mut loci = Vec::with_capacity(compositions.loci.len());
    for composition in &compositions.loci {
        let copy_count = composition
            .alternate_dosages
            .iter()
            .map(InitialAlternateDosage::dosage)
            .try_fold(composition.reference_dosage, |total, dosage| {
                total.checked_add(dosage)
            })
            .ok_or_else(|| invalid_composition(composition))?;
        if !(1..=10).contains(&copy_count) {
            return Err(invalid_composition(composition));
        }

        let mut allele_dosages = Vec::with_capacity(composition.alternate_dosages.len() + 1);
        allele_dosages.push(composition.reference_dosage);
        allele_dosages.extend(
            composition
                .alternate_dosages
                .iter()
                .map(InitialAlternateDosage::dosage),
        );
        let alternate_alleles = canonical_alternate_alleles(composition);
        let reference_bytes = [composition.reference];
        let reference_base =
            std::str::from_utf8(&reference_bytes).map_err(|_| invalid_composition(composition))?;
        let mut rng = derive_subseed(
            root_seed,
            RuntimeSubseedScope::InitialPhaseLocus {
                scenario_id: "direct",
                reference_id: composition.reference_id(),
                population_id: composition.population_id().as_str(),
                individual_id: composition.individual_id().as_str(),
                sequence_id: composition.contig().as_str(),
                position_1based: composition.position_1based(),
                reference_base,
                alternate_alleles: &alternate_alleles,
            },
        )
        .map_err(InitialPhaseAssignmentError::Derivation)?
        .local_rng();
        let alleles_by_phase = phase_permutation(&mut rng, &allele_dosages)
            .map_err(InitialPhaseAssignmentError::Random)?
            .into_iter()
            .enumerate()
            .map(|(copy_index, allele_index)| {
                let copy_index =
                    u8::try_from(copy_index).map_err(|_| invalid_composition(composition))?;
                Ok(InitialPhasedAllele {
                    phase: HaplotypePhase::from_copy_index(copy_index)
                        .map_err(|_| invalid_composition(composition))?,
                    allele_index,
                })
            })
            .collect::<Result<Vec<_>, InitialPhaseAssignmentError>>()?;
        loci.push(InitialPhasedLocus {
            composition: composition.clone(),
            alleles_by_phase,
        });
    }
    Ok(InitialPhasedGenotypes {
        loci,
        matrix_reference: compositions.matrix_reference.clone(),
        matrix_fingerprint: compositions.matrix_fingerprint.clone(),
    })
}

fn canonical_alternate_alleles(composition: &InitialLocusComposition) -> String {
    composition
        .alternate_dosages
        .iter()
        .map(|alternate| char::from(alternate.alternate()).to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn invalid_composition(composition: &InitialLocusComposition) -> InitialPhaseAssignmentError {
    InitialPhaseAssignmentError::InvalidComposition {
        contig: composition.contig().as_str().to_owned(),
        position_1based: composition.position_1based(),
    }
}

/// Muestrea composición por locus e individuo con subsemillas aisladas.
pub fn sample_initial_genotype_compositions(
    reference_id: &str,
    reference: &Reference,
    root_seed: RootSeed,
    initial_populations: &[EffectiveInitialPopulation],
    matrix: &PreparedPopulationSnvs,
) -> Result<InitialGenotypeCompositions, InitialGenotypeSamplingError> {
    preflight_reference(reference_id, reference, matrix)?;
    preflight_populations(initial_populations, matrix)?;
    let loci = canonical_loci(matrix)?;
    let mut cohorts = initial_populations.iter().collect::<Vec<_>>();
    cohorts.sort_by(|left, right| left.population_id().cmp(right.population_id()));

    let mut sampled = Vec::new();
    for cohort in cohorts {
        let mut individuals = cohort.individual_ids().iter().collect::<Vec<_>>();
        individuals.sort();
        for individual_id in individuals {
            for locus in &loci {
                let rows = locus
                    .rows_by_population
                    .get(cohort.population_id())
                    .ok_or_else(|| InitialGenotypeSamplingError::MissingMatrixPopulation {
                        population_id: cohort.population_id().as_str().to_owned(),
                    })?;
                let frequencies = matrix.allele_frequencies_for_rows(rows).ok_or_else(|| {
                    InitialGenotypeSamplingError::InvalidLocusFrequencies {
                        contig: locus.contig.as_str().to_owned(),
                        position_1based: locus.position_1based,
                    }
                })?;
                let mut rng = derive_subseed(
                    root_seed,
                    RuntimeSubseedScope::InitialGenotype {
                        scenario_id: "direct",
                        reference_id,
                        population_id: cohort.population_id().as_str(),
                        individual_id: individual_id.as_str(),
                        sequence_id: locus.contig.as_str(),
                        position_1based: locus.position_1based,
                        reference_base: std::str::from_utf8(&[locus.reference])
                            .expect("prepared SNV reference is ASCII"),
                        alternate_alleles: &locus.alternate_alleles,
                    },
                )
                .map_err(InitialGenotypeSamplingError::Derivation)?
                .local_rng();
                let dosage = match cohort.genotype_model() {
                    EffectiveInitialGenotypeModel::HardyWeinberg => {
                        hwe_composition(&mut rng, &frequencies, cohort.ploidy().copy_count())
                            .map_err(InitialGenotypeSamplingError::Random)?
                    }
                    EffectiveInitialGenotypeModel::FDiploid {
                        inbreeding_coefficient,
                    } => {
                        let genotype =
                            f_diploid_genotype(&mut rng, &frequencies, inbreeding_coefficient)
                                .map_err(InitialGenotypeSamplingError::Random)?;
                        let mut dosage = vec![0_u8; frequencies.len()];
                        dosage[usize::from(genotype[0])] += 1;
                        dosage[usize::from(genotype[1])] += 1;
                        dosage
                    }
                };
                sampled.push(InitialLocusComposition {
                    reference_id: reference_id.to_owned(),
                    population_id: cohort.population_id().clone(),
                    individual_id: individual_id.clone(),
                    contig: locus.contig.clone(),
                    position_1based: locus.position_1based,
                    reference: locus.reference,
                    reference_dosage: dosage[0],
                    alternate_dosages: rows
                        .iter()
                        .zip(dosage.iter().copied().skip(1))
                        .map(|(row, dosage)| InitialAlternateDosage {
                            variant_id: row.variant_id().clone(),
                            alternate: row.alternate(),
                            dosage,
                        })
                        .collect(),
                });
            }
        }
    }
    Ok(InitialGenotypeCompositions {
        loci: sampled,
        matrix_reference: matrix.prepared_reference().clone(),
        matrix_fingerprint: matrix.canonical_fingerprint(),
    })
}

fn preflight_reference(
    reference_id: &str,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
) -> Result<(), InitialGenotypeSamplingError> {
    let sampled_reference = PreparedReferenceIdentity::from_reference(reference);
    if sampled_reference.reference_id().as_str() != reference_id
        || matrix.prepared_reference() != &sampled_reference
    {
        return Err(InitialGenotypeSamplingError::ReferenceMismatch {
            configured_reference_id: reference_id.to_owned(),
            sampled_reference,
            prepared_reference: matrix.prepared_reference().clone(),
        });
    }
    Ok(())
}

fn preflight_populations(
    initial_populations: &[EffectiveInitialPopulation],
    matrix: &PreparedPopulationSnvs,
) -> Result<(), InitialGenotypeSamplingError> {
    if matrix.rows().is_empty() {
        return Ok(());
    }
    let matrix_populations = matrix.population_ids().iter().collect::<BTreeSet<_>>();
    for cohort in initial_populations {
        if !matrix_populations.contains(cohort.population_id()) {
            return Err(InitialGenotypeSamplingError::MissingMatrixPopulation {
                population_id: cohort.population_id().as_str().to_owned(),
            });
        }
    }
    Ok(())
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct LocusIdentity {
    contig: SequenceId,
    position_1based: u64,
    reference: u8,
}

struct CanonicalLocus<'a> {
    contig: SequenceId,
    position_1based: u64,
    reference: u8,
    alternate_alleles: String,
    rows_by_population: BTreeMap<PopulationId, Vec<&'a PreparedPopulationSnv>>,
}

fn canonical_loci(
    matrix: &PreparedPopulationSnvs,
) -> Result<Vec<CanonicalLocus<'_>>, InitialGenotypeSamplingError> {
    // PreparedPopulationSnvs ya usa el orden del FASTA. Lo conservamos en el
    // resultado en lugar de ordenar contigs lexicográficamente, que alteraría
    // el orden observable de una referencia con IDs no lexicográficos.
    let mut locus_indices = BTreeMap::<LocusIdentity, usize>::new();
    let mut grouped = Vec::<(
        LocusIdentity,
        BTreeMap<PopulationId, Vec<&PreparedPopulationSnv>>,
    )>::new();
    for row in matrix.rows() {
        let identity = LocusIdentity {
            contig: row.contig().clone(),
            position_1based: row.position_1based(),
            reference: row.reference(),
        };
        let index = match locus_indices.get(&identity) {
            Some(index) => *index,
            None => {
                let index = grouped.len();
                locus_indices.insert(identity.clone(), index);
                grouped.push((identity, BTreeMap::new()));
                index
            }
        };
        grouped[index]
            .1
            .entry(row.population_id().clone())
            .or_default()
            .push(row);
    }
    grouped
        .into_iter()
        .map(|(identity, mut rows_by_population)| {
            for rows in rows_by_population.values_mut() {
                rows.sort_by(|left, right| {
                    left.alternate()
                        .cmp(&right.alternate())
                        .then(left.variant_id().cmp(right.variant_id()))
                });
            }
            let Some(reference_rows) = rows_by_population.values().next() else {
                return Err(InitialGenotypeSamplingError::InconsistentLocus {
                    contig: identity.contig.as_str().to_owned(),
                    position_1based: identity.position_1based,
                });
            };
            let alternate_alleles = reference_rows
                .iter()
                .map(|row| char::from(row.alternate()).to_string())
                .collect::<Vec<_>>()
                .join(",");
            if alternate_alleles.is_empty()
                || rows_by_population.values().any(|rows| {
                    rows.len() != reference_rows.len()
                        || rows.iter().zip(reference_rows).any(|(left, right)| {
                            left.alternate() != right.alternate()
                                || left.variant_id() != right.variant_id()
                        })
                })
            {
                return Err(InitialGenotypeSamplingError::InconsistentLocus {
                    contig: identity.contig.as_str().to_owned(),
                    position_1based: identity.position_1based,
                });
            }
            Ok(CanonicalLocus {
                contig: identity.contig,
                position_1based: identity.position_1based,
                reference: identity.reference,
                alternate_alleles,
                rows_by_population,
            })
        })
        .collect()
}
