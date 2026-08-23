//! Composición genotípica inicial por locus, antes de la fase y materialización.
//!
//! F2.3.2 consume una matriz ya preparada, por lo que no reinterpreta `alt_af`,
//! no elige loci y no modifica la matriz. Antes de derivar RNG exige que la
//! referencia declarada coincida con la identidad y huella canónica retenidas
//! durante la preparación. La fase queda deliberadamente fuera de este módulo: F2.3.3
//! consumirá estas dosis con su propio ámbito RNG.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::application::{
    EffectiveConfig, EffectiveInitialGenotypeModel, EffectiveInitialPopulation,
};
use crate::identity::{IndividualId, PopulationId, SequenceId, VariantId};
use crate::population_snvs::{
    PreparedPopulationSnv, PreparedPopulationSnvs, PreparedReferenceIdentity,
};
use crate::randomness::{
    F_DIPLOID_GENOTYPE_ALGORITHM, HWE_COMPOSITION_ALGORITHM, RandomSamplingError, RootSeed,
    RuntimeSubseedScope, derive_subseed, f_diploid_genotype, hwe_composition,
};
use crate::reference::Reference;

/// La composición HWE conserva una dosis por alelo, sin fase.
pub const INITIAL_GENOTYPE_HWE_ALGORITHM: &str = HWE_COMPOSITION_ALGORITHM;
/// La composición F conserva una dosis por alelo, sin fase.
pub const INITIAL_GENOTYPE_F_DIPLOID_ALGORITHM: &str = F_DIPLOID_GENOTYPE_ALGORITHM;
/// Ámbito de subsemilla para la composición genotípica inicial.
pub const INITIAL_GENOTYPE_SCOPE: &str = "initial_genotype";

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
}

impl InitialGenotypeCompositions {
    #[must_use]
    pub fn loci(&self) -> &[InitialLocusComposition] {
        &self.loci
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
                let frequencies = allele_frequencies(rows, locus)?;
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
    Ok(InitialGenotypeCompositions { loci: sampled })
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

fn allele_frequencies(
    rows: &[&PreparedPopulationSnv],
    locus: &CanonicalLocus<'_>,
) -> Result<Vec<f64>, InitialGenotypeSamplingError> {
    let alternate_sum = neumaier_sum(rows.iter().map(|row| row.alt_af()));
    let reference = 1.0 - alternate_sum;
    if !alternate_sum.is_finite()
        || alternate_sum < 0.0
        || reference < -1e-12
        || rows
            .iter()
            .any(|row| !row.alt_af().is_finite() || row.alt_af() < 0.0)
    {
        return Err(InitialGenotypeSamplingError::InvalidLocusFrequencies {
            contig: locus.contig.as_str().to_owned(),
            position_1based: locus.position_1based,
        });
    }
    let mut frequencies = Vec::with_capacity(rows.len() + 1);
    frequencies.push(reference.max(0.0));
    frequencies.extend(rows.iter().map(|row| row.alt_af()));
    Ok(frequencies)
}

fn neumaier_sum(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut correction = 0.0;
    for value in values {
        let next = sum + value;
        if sum.abs() >= value.abs() {
            correction += (sum - next) + value;
        } else {
            correction += (value - next) + sum;
        }
        sum = next;
    }
    sum + correction
}
