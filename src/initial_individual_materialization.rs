//! Materialización final y resumen de frecuencias de fundadores F2.3.5.
//!
//! El estado canónico continúa siendo la cadena de segmentos de cada copia. La
//! materialización de bases se conserva aquí únicamente como la salida exacta
//! solicitada por F2.3.5; nunca sustituye a segmentos, matriz o fase como fuente
//! de verdad. Las frecuencias realizadas se calculan después de integrar los
//! cambios directos, sin modificar las copias sorteadas.

use std::collections::BTreeMap;
use std::fmt;

use serde_json::{Value, json};

use crate::application::{EffectiveConfig, EffectiveInitialGenotypeModel};
use crate::coordinates::BasePosition;
use crate::haplotype::{HaplotypePhase, Individual};
use crate::identity::{IndividualId, PopulationId, SequenceId};
use crate::initial_genotypes::{assign_initial_phases, sample_initial_genotypes};
use crate::initial_individual_integration::{
    ExplicitSnpConstraint, IncompatibleGenotypeAuthorization, InitialIndividualIntegrationError,
    InitialIndividualIntegrationRequest, integrate_initial_individual,
};
use crate::materialization::{
    GenotypeQueryError, MaterializationError, materialize_copy, observed_snp_allele,
};
use crate::population_snvs::{PreparedPopulationSnv, PreparedPopulationSnvs};
use crate::reference::Reference;
use crate::segments::CanonicalBase;
use crate::small_event_composition::DeclaredSmallEvent;

/// Secuencia final de una copia. Es una vista materializada de salida; el
/// individuo retiene de forma canónica únicamente sus segmentos compartidos.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaterializedInitialCopy {
    sequence_id: SequenceId,
    phase: HaplotypePhase,
    bases: Vec<u8>,
}

impl MaterializedInitialCopy {
    #[must_use]
    pub const fn sequence_id(&self) -> &SequenceId {
        &self.sequence_id
    }

    #[must_use]
    pub const fn phase(&self) -> HaplotypePhase {
        self.phase
    }

    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }
}

/// Fundador final, con sus copias compactas y sus bytes materializados.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct MaterializedInitialIndividual {
    population_id: PopulationId,
    individual: Individual,
    copies: Vec<MaterializedInitialCopy>,
    explicit_snp_constraints: Vec<ExplicitSnpConstraint>,
}

impl MaterializedInitialIndividual {
    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub const fn individual(&self) -> &Individual {
        &self.individual
    }

    #[must_use]
    pub fn copies(&self) -> &[MaterializedInitialCopy] {
        &self.copies
    }

    #[must_use]
    pub fn explicit_snp_constraints(&self) -> &[ExplicitSnpConstraint] {
        &self.explicit_snp_constraints
    }
}

/// Resultado completo de F2.3.5. El resumen conserva matriz, modelo, objetivos
/// y realizaciones; las frecuencias nunca retroalimentan el muestreo.
#[derive(Clone, Debug, PartialEq)]
pub struct InitialIndividualsMaterialization {
    individuals: Vec<MaterializedInitialIndividual>,
    statistics: Value,
}

impl InitialIndividualsMaterialization {
    #[must_use]
    pub fn individuals(&self) -> &[MaterializedInitialIndividual] {
        &self.individuals
    }

    /// Documento estable que entra en configuración efectiva y manifiesto.
    #[must_use]
    pub const fn statistics(&self) -> &Value {
        &self.statistics
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum InitialIndividualsMaterializationError {
    Integration(InitialIndividualIntegrationError),
    Materialization(MaterializationError),
    Genotype(GenotypeQueryError),
    MissingSequence {
        individual_id: IndividualId,
        sequence_id: SequenceId,
    },
    MissingReferenceSequence {
        sequence_id: SequenceId,
    },
    InvalidBase {
        contig: SequenceId,
        position_1based: u64,
    },
    InvalidStatistics,
    Sampling(String),
}

impl fmt::Display for InitialIndividualsMaterializationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for InitialIndividualsMaterializationError {}

/// Integra todos los fundadores, materializa cada copia y calcula sus
/// estadísticas desde el estado final. El orden es población e individuo por
/// identidad, por lo que añadir o reordenar entidades no cambia las existentes.
pub fn materialize_initial_individuals(
    config: &EffectiveConfig,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
    authorization: IncompatibleGenotypeAuthorization,
    direct_events: impl Fn(
        &IndividualId,
    )
        -> Result<Vec<DeclaredSmallEvent>, InitialIndividualsMaterializationError>,
) -> Result<InitialIndividualsMaterialization, InitialIndividualsMaterializationError> {
    let compositions = sample_initial_genotypes(config, reference, matrix)
        .map_err(|error| InitialIndividualsMaterializationError::Sampling(error.to_string()))?;
    let phased = assign_initial_phases(config.reproducibility().root_seed(), &compositions)
        .map_err(|error| InitialIndividualsMaterializationError::Sampling(error.to_string()))?;

    let mut populations = config.initial_populations().iter().collect::<Vec<_>>();
    populations.sort_by(|left, right| left.population_id().cmp(right.population_id()));
    let mut individuals = Vec::new();
    for population in populations {
        let mut ids = population.individual_ids().iter().collect::<Vec<_>>();
        ids.sort();
        for individual_id in ids {
            let declared_events = direct_events(individual_id)?;
            let original =
                Individual::from_reference(individual_id.clone(), reference, population.ploidy());
            let integrated = integrate_initial_individual(
                &original,
                reference,
                InitialIndividualIntegrationRequest::new(
                    matrix,
                    &phased.for_individual(individual_id),
                    config.reproducibility().root_seed(),
                    &declared_events,
                    &[],
                    authorization,
                    population.genotype_model(),
                ),
            )
            .map_err(InitialIndividualsMaterializationError::Integration)?;
            let copies = materialize_individual(integrated.individual(), reference)?;
            individuals.push(MaterializedInitialIndividual {
                population_id: population.population_id().clone(),
                individual: integrated.individual().clone(),
                copies,
                explicit_snp_constraints: integrated.explicit_snp_constraints().to_vec(),
            });
        }
    }
    let statistics = statistics_document(config, reference, matrix, &individuals)?;
    Ok(InitialIndividualsMaterialization {
        individuals,
        statistics,
    })
}

fn materialize_individual(
    individual: &Individual,
    reference: &Reference,
) -> Result<Vec<MaterializedInitialCopy>, InitialIndividualsMaterializationError> {
    let mut copies = Vec::new();
    for sequence in individual.sequences() {
        for copy in sequence.copies() {
            copies.push(MaterializedInitialCopy {
                sequence_id: sequence.sequence_id().clone(),
                phase: copy.phase(),
                bases: materialize_copy(copy, reference)
                    .map_err(InitialIndividualsMaterializationError::Materialization)?,
            });
        }
    }
    copies.sort_by(|left, right| {
        left.sequence_id
            .cmp(&right.sequence_id)
            .then(left.phase.cmp(&right.phase))
    });
    Ok(copies)
}

#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
struct LocusKey {
    population_id: PopulationId,
    contig: SequenceId,
    position_1based: u64,
    reference: u8,
}

fn statistics_document(
    config: &EffectiveConfig,
    reference: &Reference,
    matrix: &PreparedPopulationSnvs,
    individuals: &[MaterializedInitialIndividual],
) -> Result<Value, InitialIndividualsMaterializationError> {
    let mut rows_by_locus = BTreeMap::<LocusKey, Vec<&PreparedPopulationSnv>>::new();
    for row in matrix.rows() {
        rows_by_locus
            .entry(LocusKey {
                population_id: row.population_id().clone(),
                contig: row.contig().clone(),
                position_1based: row.position_1based(),
                reference: row.reference(),
            })
            .or_default()
            .push(row);
    }
    for rows in rows_by_locus.values_mut() {
        rows.sort_by(|left, right| {
            left.alternate()
                .cmp(&right.alternate())
                .then(left.variant_id().cmp(right.variant_id()))
        });
    }
    let populations = config
        .initial_populations()
        .iter()
        .map(|population| {
            let loci = rows_by_locus
                .iter()
                .filter(|(key, _)| key.population_id == *population.population_id())
                .map(|(key, rows)| {
                    locus_statistics(
                        reference,
                        population.genotype_model(),
                        population.ploidy().copy_count(),
                        key,
                        rows,
                        individuals,
                    )
                })
                .collect::<Result<Vec<_>, _>>()?;
            Ok(json!({
                "population_id": population.population_id().as_str(),
                "genotype_model": genotype_model_json(population.genotype_model()),
                "individual_count": population.individual_count(),
                "ploidy": population.ploidy().copy_count(),
                "loci": loci,
            }))
        })
        .collect::<Result<Vec<_>, InitialIndividualsMaterializationError>>()?;
    let explicit_snp_constraints = individuals
        .iter()
        .flat_map(|individual| {
            individual.explicit_snp_constraints.iter().map(move |constraint| {
                json!({
                    "population_id": individual.population_id.as_str(),
                    "individual_id": individual.individual.id().as_str(),
                    "event_id": constraint.event_id().as_str(),
                    "phase": constraint.phase().copy_index(),
                    "matrix_variant_id": constraint.matrix_variant_id().as_str(),
                    "sampled_alternate": constraint.sampled_alternate().map(|base| char::from(base).to_string()),
                    "explicit_alternate": char::from(constraint.explicit_alternate()).to_string(),
                    "incompatible": constraint.is_incompatible(),
                })
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "kind": "initial_population_statistics",
        "version": "1",
        "provenance": {
            "matrix_reference_id": matrix.reference_id().as_str(),
            "matrix_reference_fingerprint": matrix.reference_fingerprint(),
            "matrix_canonical_tsv_sha256": matrix.canonical_fingerprint(),
            "genotype_algorithms": [
                crate::initial_genotypes::INITIAL_GENOTYPE_HWE_ALGORITHM,
                crate::initial_genotypes::INITIAL_GENOTYPE_F_DIPLOID_ALGORITHM,
                crate::initial_genotypes::INITIAL_PHASE_ALGORITHM,
            ],
            "missing_locus_policy": "missing_by_deletion_excluded_from_allele_frequency_denominator_v1",
            "explicit_snp_constraints": explicit_snp_constraints,
        },
        "populations": populations,
    }))
}

fn locus_statistics(
    reference: &Reference,
    model: EffectiveInitialGenotypeModel,
    ploidy: u8,
    key: &LocusKey,
    rows: &[&PreparedPopulationSnv],
    individuals: &[MaterializedInitialIndividual],
) -> Result<Value, InitialIndividualsMaterializationError> {
    let reference_base = CanonicalBase::from_byte(key.reference).map_err(|_| {
        InitialIndividualsMaterializationError::InvalidBase {
            contig: key.contig.clone(),
            position_1based: key.position_1based,
        }
    })?;
    let position = BasePosition::from_base_1(
        reference.id().clone(),
        key.contig.clone(),
        key.position_1based,
        reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == &key.contig)
            .ok_or_else(
                || InitialIndividualsMaterializationError::MissingReferenceSequence {
                    sequence_id: key.contig.clone(),
                },
            )?
            .length(),
    )
    .map_err(|_| InitialIndividualsMaterializationError::InvalidBase {
        contig: key.contig.clone(),
        position_1based: key.position_1based,
    })?;
    let mut frequencies = Vec::with_capacity(rows.len() + 1);
    let alt_sum = neumaier_sum(rows.iter().map(|row| row.alt_af()));
    let reference_frequency = 1.0 - alt_sum;
    if !reference_frequency.is_finite() || reference_frequency < 0.0 {
        return Err(InitialIndividualsMaterializationError::InvalidStatistics);
    }
    frequencies.push(reference_frequency);
    frequencies.extend(rows.iter().map(|row| row.alt_af()));
    if frequencies
        .iter()
        .any(|frequency| !frequency.is_finite() || *frequency < 0.0)
    {
        return Err(InitialIndividualsMaterializationError::InvalidStatistics);
    }
    let mut counts = vec![0_u64; frequencies.len()];
    let mut observed_copies = 0_u64;
    let mut missing_by_deletion = 0_u64;
    let mut realized_genotypes = BTreeMap::<Vec<u8>, u64>::new();
    for individual in individuals
        .iter()
        .filter(|individual| individual.population_id == key.population_id)
    {
        let copies = individual
            .individual
            .copies_for(&key.contig)
            .ok_or_else(|| InitialIndividualsMaterializationError::MissingSequence {
                individual_id: individual.individual.id().clone(),
                sequence_id: key.contig.clone(),
            })?;
        let mut genotype = Vec::with_capacity(copies.len());
        let mut complete = true;
        for copy in copies {
            match observed_snp_allele(copy, &position, reference_base) {
                Ok(allele) => {
                    let index = if allele == reference_base {
                        0
                    } else {
                        rows.iter()
                            .position(|row| row.alternate() == allele.as_byte())
                            .map(|index| index + 1)
                            .ok_or(InitialIndividualsMaterializationError::InvalidStatistics)?
                    };
                    counts[index] += 1;
                    observed_copies += 1;
                    genotype.push(
                        u8::try_from(index).map_err(|_| {
                            InitialIndividualsMaterializationError::InvalidStatistics
                        })?,
                    );
                }
                Err(GenotypeQueryError::LocusDeletedFromCopy { .. }) => {
                    missing_by_deletion += 1;
                    complete = false;
                }
                Err(error) => return Err(InitialIndividualsMaterializationError::Genotype(error)),
            }
        }
        if complete {
            genotype.sort_unstable();
            *realized_genotypes.entry(genotype).or_default() += 1;
        }
    }
    let total_copies = u64::try_from(
        individuals
            .iter()
            .filter(|individual| individual.population_id == key.population_id)
            .count(),
    )
    .map_err(|_| InitialIndividualsMaterializationError::InvalidStatistics)?
    .checked_mul(u64::from(ploidy))
    .ok_or(InitialIndividualsMaterializationError::InvalidStatistics)?;
    let population_individual_count = individuals
        .iter()
        .filter(|individual| individual.population_id == key.population_id)
        .count();
    let population_individual_count = u64::try_from(population_individual_count)
        .map_err(|_| InitialIndividualsMaterializationError::InvalidStatistics)?;
    let allele_rows = std::iter::once(json!({
        "allele": char::from(key.reference).to_string(),
        "kind": "ref",
        "requested_frequency": reference_frequency,
        "realized_count": counts[0],
        "realized_frequency": frequency(counts[0], observed_copies)?,
    }))
    .chain(rows.iter().enumerate().map(|(index, row)| {
        json!({
            "allele": char::from(row.alternate()).to_string(),
            "kind": "alt",
            "variant_id": row.variant_id().as_str(),
            "source_id": row.source_id(),
            "requested_frequency": row.alt_af(),
            "realized_count": counts[index + 1],
            "realized_frequency": frequency(counts[index + 1], observed_copies).unwrap_or(0.0),
        })
    }))
    .collect::<Vec<_>>();
    let expected = expected_genotypes(&frequencies, ploidy, model)?;
    let genotypes = expected
        .into_iter()
        .map(|(dosages, probability)| {
            let key = expanded_genotype(&dosages)?;
            Ok(json!({
                "allele_dosages": dosages,
                "expected_probability": probability,
                "expected_count": probability * population_individual_count as f64,
                "realized_count": realized_genotypes.get(&key).copied().unwrap_or(0),
            }))
        })
        .collect::<Result<Vec<_>, InitialIndividualsMaterializationError>>()?;
    Ok(json!({
        "contig": key.contig.as_str(),
        "position_1based": key.position_1based,
        "ref": char::from(key.reference).to_string(),
        "total_copies": total_copies,
        "observed_copies": observed_copies,
        "missing_by_deletion": missing_by_deletion,
        "alleles": allele_rows,
        "complete_genotype_observations": realized_genotypes.values().sum::<u64>(),
        "genotypes": genotypes,
    }))
}

fn frequency(count: u64, denominator: u64) -> Result<f64, InitialIndividualsMaterializationError> {
    if denominator == 0 {
        return Ok(0.0);
    }
    let value = count as f64 / denominator as f64;
    if value.is_finite() {
        Ok(value)
    } else {
        Err(InitialIndividualsMaterializationError::InvalidStatistics)
    }
}

fn genotype_model_json(model: EffectiveInitialGenotypeModel) -> Value {
    match model {
        EffectiveInitialGenotypeModel::HardyWeinberg => json!({"kind": "hwe"}),
        EffectiveInitialGenotypeModel::FDiploid {
            inbreeding_coefficient,
        } => json!({"kind": "f_diploid", "inbreeding_coefficient": inbreeding_coefficient}),
    }
}

fn expected_genotypes(
    frequencies: &[f64],
    ploidy: u8,
    model: EffectiveInitialGenotypeModel,
) -> Result<Vec<(Vec<u8>, f64)>, InitialIndividualsMaterializationError> {
    match model {
        EffectiveInitialGenotypeModel::HardyWeinberg => {
            let mut values = Vec::new();
            let mut current = vec![0_u8; frequencies.len()];
            enumerate_hwe(frequencies, ploidy, 0, &mut current, &mut values)?;
            Ok(values)
        }
        EffectiveInitialGenotypeModel::FDiploid {
            inbreeding_coefficient,
        } => {
            let mut values = Vec::new();
            for (index, frequency) in frequencies.iter().copied().enumerate() {
                let mut dosage = vec![0_u8; frequencies.len()];
                dosage[index] = 2;
                values.push((
                    dosage,
                    frequency * frequency + inbreeding_coefficient * frequency * (1.0 - frequency),
                ));
            }
            for left in 0..frequencies.len() {
                for right in (left + 1)..frequencies.len() {
                    let mut dosage = vec![0_u8; frequencies.len()];
                    dosage[left] = 1;
                    dosage[right] = 1;
                    values.push((
                        dosage,
                        2.0 * (1.0 - inbreeding_coefficient)
                            * frequencies[left]
                            * frequencies[right],
                    ));
                }
            }
            Ok(values)
        }
    }
}

fn enumerate_hwe(
    frequencies: &[f64],
    remaining: u8,
    index: usize,
    current: &mut [u8],
    values: &mut Vec<(Vec<u8>, f64)>,
) -> Result<(), InitialIndividualsMaterializationError> {
    if index + 1 == frequencies.len() {
        current[index] = remaining;
        let mut probability = factorial(
            u8::try_from(current.iter().map(|value| u16::from(*value)).sum::<u16>())
                .map_err(|_| InitialIndividualsMaterializationError::InvalidStatistics)?,
        );
        for (dosage, frequency) in current.iter().zip(frequencies) {
            probability /= factorial(*dosage);
            probability *= frequency.powi(i32::from(*dosage));
        }
        if !probability.is_finite() {
            return Err(InitialIndividualsMaterializationError::InvalidStatistics);
        }
        values.push((current.to_vec(), probability));
        return Ok(());
    }
    for dosage in 0..=remaining {
        current[index] = dosage;
        enumerate_hwe(frequencies, remaining - dosage, index + 1, current, values)?;
    }
    Ok(())
}

fn factorial(value: u8) -> f64 {
    (1..=u64::from(value)).fold(1.0, |result, item| result * item as f64)
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

fn expanded_genotype(dosages: &[u8]) -> Result<Vec<u8>, InitialIndividualsMaterializationError> {
    let mut result = Vec::new();
    for (index, dosage) in dosages.iter().copied().enumerate() {
        let allele = u8::try_from(index)
            .map_err(|_| InitialIndividualsMaterializationError::InvalidStatistics)?;
        result.extend(std::iter::repeat_n(allele, usize::from(dosage)));
    }
    Ok(result)
}
