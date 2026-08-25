//! Decisión final de representación lógica VCF para F2.5.4.
//!
//! Esta capa recibe el resultado ya normalizado, agrupado y faseado de F2.5.3.
//! No escribe un VCF, no crea ficheros complementarios y no modifica la verdad.
//! Su cometido es conservar `N` literalmente cuando el alelo explícito VCF lo
//! admite sin pérdida, y rechazar el único SNP prohibido con `REF=N` antes de
//! que un escritor pueda emitirlo. F2 no produce casos con pérdida real: por
//! eso el complemento condicional es `None` para toda salida F2 válida.

use std::collections::BTreeSet;

use crate::vcf_grouping::{GroupedVcfRecord, PhasedGenotype};

/// Versión reservada para el fichero complementario que un consumidor futuro
/// podrá escribir únicamente cuando exista una pérdida real.
///
/// F2.5.4 no materializa ese fichero ni inventa entradas de complemento.
pub const VCF_COMPLEMENT_SCHEMA_VERSION: &str = "vcf-complement-v1";

/// Salida lógica que F2.5.5 podrá escribir sin volver a decidir `N`.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConditionalVcfRecords {
    /// Registros VCF completos y representables sin pérdida.
    pub records: Vec<GroupedVcfRecord>,
    /// Ausente mientras ningún caso haya perdido detalle. En F2 siempre es
    /// `None`: la primera entrada positiva pertenece a F5.5.
    complement: Option<ConditionalComplement>,
}

impl ConditionalVcfRecords {
    /// Devuelve el complemento únicamente si una fase posterior produjo una
    /// pérdida real. F2 siempre devuelve `None`.
    #[must_use]
    pub const fn complement(&self) -> Option<&ConditionalComplement> {
        self.complement.as_ref()
    }
}

/// Marcador estructurado del complemento condicional del conjunto.
///
/// La cardinalidad, los códigos de pérdida y las entradas positivas se fijan
/// con el primer caso real (F5.5); anticiparlos aquí permitiría fabricar una
/// pérdida F2 inexistente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ConditionalComplement {
    schema_version: &'static str,
}

impl ConditionalComplement {
    #[must_use]
    pub const fn schema_version(&self) -> &'static str {
        self.schema_version
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfRepresentationError {
    EmptyReference {
        sequence_id: String,
        position: u64,
    },
    EmptyAlternate {
        sequence_id: String,
        position: u64,
    },
    NoAlternates {
        sequence_id: String,
        position: u64,
    },
    InvalidAlleleBase {
        sequence_id: String,
        position: u64,
        allele: String,
    },
    ReferenceEqualsAlternate {
        sequence_id: String,
        position: u64,
        allele: String,
    },
    /// DTC-032/REQ-VER-47: un SNP no puede declarar una referencia desconocida.
    NReferenceSnp {
        sequence_id: String,
        position: u64,
    },
    NonCanonicalCases {
        sequence_id: String,
        position: u64,
    },
    AlternateCasesMismatch {
        sequence_id: String,
        position: u64,
    },
    InvalidGenotype {
        sequence_id: String,
        position: u64,
        sample: String,
    },
}

impl std::fmt::Display for VcfRepresentationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for VcfRepresentationError {}

/// Verifica qué registros agrupados se pueden expresar en VCF sin pérdida.
///
/// El resultado conserva exactamente `REF`, `ALT`, `GT`, dosis y `GPS_CASE`.
/// `N` se deja intacta en alelos explícitos de indels (o inserciones), que VCF
/// 4.5 puede expresar. Un SNP con `REF=N` no se omite silenciosamente ni se
/// convierte en A/C/G/T: es un estado que F2 ya debe haber excluido y se
/// devuelve como error estructurado.
pub fn prepare_conditional_vcf_records(
    records: Vec<GroupedVcfRecord>,
) -> Result<ConditionalVcfRecords, VcfRepresentationError> {
    for record in &records {
        validate_record(record)?;
    }

    Ok(ConditionalVcfRecords {
        records,
        complement: None,
    })
}

fn validate_record(record: &GroupedVcfRecord) -> Result<(), VcfRepresentationError> {
    let position = record.position;
    let sequence_id = &record.sequence_id;
    if record.reference.is_empty() {
        return Err(VcfRepresentationError::EmptyReference {
            sequence_id: sequence_id.clone(),
            position,
        });
    }
    validate_allele(sequence_id, position, &record.reference)?;
    if record.alternates.is_empty() {
        return Err(VcfRepresentationError::NoAlternates {
            sequence_id: sequence_id.clone(),
            position,
        });
    }
    for alternate in &record.alternates {
        if alternate.is_empty() {
            return Err(VcfRepresentationError::EmptyAlternate {
                sequence_id: sequence_id.clone(),
                position,
            });
        }
        validate_allele(sequence_id, position, alternate)?;
        if alternate == &record.reference {
            return Err(VcfRepresentationError::ReferenceEqualsAlternate {
                sequence_id: sequence_id.clone(),
                position,
                allele: alternate.clone(),
            });
        }
    }
    if record.reference == "N"
        && record
            .alternates
            .iter()
            .any(|alternate| alternate.len() == 1)
    {
        return Err(VcfRepresentationError::NReferenceSnp {
            sequence_id: sequence_id.clone(),
            position,
        });
    }
    validate_cases(record)?;
    for (sample, genotype) in &record.samples {
        validate_genotype(record, sample, genotype)?;
    }
    Ok(())
}

fn validate_allele(
    sequence_id: &str,
    position: u64,
    allele: &str,
) -> Result<(), VcfRepresentationError> {
    if allele
        .bytes()
        .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
    {
        Ok(())
    } else {
        Err(VcfRepresentationError::InvalidAlleleBase {
            sequence_id: sequence_id.to_owned(),
            position,
            allele: allele.to_owned(),
        })
    }
}

fn validate_cases(record: &GroupedVcfRecord) -> Result<(), VcfRepresentationError> {
    let invalid = record.gps_cases.is_empty()
        || record.gps_cases.iter().any(String::is_empty)
        || !is_strictly_sorted(&record.gps_cases)
        || record.alternate_cases.len() != record.alternates.len()
        || record.alternate_cases.iter().any(|cases| {
            cases.is_empty() || cases.iter().any(String::is_empty) || !is_strictly_sorted(cases)
        });
    if invalid {
        return Err(VcfRepresentationError::NonCanonicalCases {
            sequence_id: record.sequence_id.clone(),
            position: record.position,
        });
    }
    let from_alternates = record
        .alternate_cases
        .iter()
        .flatten()
        .cloned()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect::<Vec<_>>();
    if from_alternates != record.gps_cases {
        return Err(VcfRepresentationError::AlternateCasesMismatch {
            sequence_id: record.sequence_id.clone(),
            position: record.position,
        });
    }
    Ok(())
}

fn is_strictly_sorted(values: &[String]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn validate_genotype(
    record: &GroupedVcfRecord,
    sample: &str,
    genotype: &PhasedGenotype,
) -> Result<(), VcfRepresentationError> {
    let Some(allele_count) = record.alternates.len().checked_add(1) else {
        return Err(VcfRepresentationError::InvalidGenotype {
            sequence_id: record.sequence_id.clone(),
            position: record.position,
            sample: sample.to_owned(),
        });
    };
    let indices_valid = genotype
        .alleles
        .iter()
        .all(|allele| usize::try_from(*allele).is_ok_and(|index| index < allele_count));
    let dosage_matches = genotype.allele_dosages.len() == allele_count
        && genotype
            .allele_dosages
            .iter()
            .map(|dosage| usize::from(*dosage))
            .sum::<usize>()
            == genotype.alleles.len()
        && genotype
            .allele_dosages
            .iter()
            .enumerate()
            .all(|(index, dosage)| {
                usize::from(*dosage)
                    == genotype
                        .alleles
                        .iter()
                        .filter(|allele| usize::try_from(**allele).ok() == Some(index))
                        .count()
            });
    if indices_valid && dosage_matches {
        Ok(())
    } else {
        Err(VcfRepresentationError::InvalidGenotype {
            sequence_id: record.sequence_id.clone(),
            position: record.position,
            sample: sample.to_owned(),
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn record(reference: &str, alternates: &[&str]) -> GroupedVcfRecord {
        GroupedVcfRecord {
            sequence_id: "chr1".into(),
            position: 4,
            reference: reference.into(),
            alternates: alternates
                .iter()
                .map(|alternate| (*alternate).into())
                .collect(),
            samples: vec![(
                "IND_A".into(),
                PhasedGenotype {
                    alleles: vec![1, 0],
                    allele_dosages: vec![1, 1],
                },
            )],
            gps_cases: vec!["CASE_A".into()],
            alternate_cases: vec![vec!["CASE_A".into()]],
        }
    }

    #[test]
    fn keeps_n_in_explicit_indel_ref_and_alt_without_a_complement() {
        let deletion = record("AN", &["A"]);
        let insertion = record("A", &["ANT"]);
        let output = prepare_conditional_vcf_records(vec![deletion.clone(), insertion.clone()])
            .expect("N in a lossless explicit indel allele is valid");

        assert_eq!(output.records, [deletion, insertion]);
        assert_eq!(output.complement(), None);
    }

    #[test]
    fn rejects_a_snp_with_unknown_reference_without_resolving_n() {
        let error = prepare_conditional_vcf_records(vec![record("N", &["A"])])
            .expect_err("DTC-032 forbids REF=N SNPs");
        assert_eq!(
            error,
            VcfRepresentationError::NReferenceSnp {
                sequence_id: "chr1".into(),
                position: 4,
            }
        );
    }

    #[test]
    fn rejects_a_multiallelic_n_record_when_any_alt_is_a_snp() {
        let mut grouped = record("N", &["A", "NA"]);
        grouped.samples[0].1 = PhasedGenotype {
            alleles: vec![1, 2],
            allele_dosages: vec![0, 1, 1],
        };
        grouped.gps_cases = vec!["CASE_A".into(), "CASE_NA".into()];
        grouped.alternate_cases = vec![vec!["CASE_A".into()], vec!["CASE_NA".into()]];

        let error = prepare_conditional_vcf_records(vec![grouped])
            .expect_err("a REF=N SNP ALT cannot share a record with a valid indel ALT");
        assert_eq!(
            error,
            VcfRepresentationError::NReferenceSnp {
                sequence_id: "chr1".into(),
                position: 4,
            }
        );
    }

    #[test]
    fn representable_records_preserve_alt_gt_and_canonical_gps_case() {
        let mut grouped = record("C", &["A", "T"]);
        grouped.samples[0].1 = PhasedGenotype {
            alleles: vec![2, 1],
            allele_dosages: vec![0, 1, 1],
        };
        grouped.gps_cases = vec!["CASE_A".into(), "CASE_B".into()];
        grouped.alternate_cases = vec![vec!["CASE_A".into()], vec!["CASE_B".into()]];

        let output = prepare_conditional_vcf_records(vec![grouped.clone()]).expect("valid VCF");
        assert_eq!(output.records, [grouped]);
        assert!(output.complement().is_none());
    }

    #[test]
    fn rejects_noncanonical_case_relationships_before_writing() {
        let mut grouped = record("C", &["A"]);
        grouped.gps_cases = vec!["CASE_B".into(), "CASE_A".into()];
        let error = prepare_conditional_vcf_records(vec![grouped])
            .expect_err("GPS_CASE must remain canonical");
        assert!(matches!(
            error,
            VcfRepresentationError::NonCanonicalCases { .. }
        ));
    }
}
