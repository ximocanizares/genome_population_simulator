//! Generación reproducible de la referencia inicial mínima.

use crate::diagnostics::Diagnostic;
use crate::identity::{ReferenceId, SequenceId};
use crate::randomness::{RootSeed, RuntimeSubseedScope, derive_subseed};
use crate::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};

pub(crate) const COMPOSITION_TOLERANCE: f64 = 1e-9;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaseFrequency {
    pub base: u8,
    pub probability: f64,
}

#[derive(Clone, Copy, Debug)]
pub struct GeneratedSequenceSpec<'a> {
    pub id: &'a str,
    pub class: SequenceClass,
    pub length: u64,
}

pub fn generate_reference(
    reference_id: &str,
    root_seed: RootSeed,
    specs: &[GeneratedSequenceSpec<'_>],
    composition: &[BaseFrequency],
) -> Result<Reference, Diagnostic> {
    generate_reference_with_unknown_fraction(reference_id, root_seed, specs, composition, 0.0)
}

/// Genera una referencia con una fracción independiente de bases ambiguas `N`.
/// Las frecuencias A/C/G/T se aplican solo a las posiciones que no son `N`.
pub fn generate_reference_with_unknown_fraction(
    reference_id: &str,
    root_seed: RootSeed,
    specs: &[GeneratedSequenceSpec<'_>],
    composition: &[BaseFrequency],
    unknown_base_fraction: f64,
) -> Result<Reference, Diagnostic> {
    let cumulative = validate_composition(composition)?;
    let unknown_base_fraction = validate_unknown_base_fraction(unknown_base_fraction)?;
    let mut rng = derive_subseed(
        root_seed,
        RuntimeSubseedScope::ReferenceGeneration { reference_id },
    )?
    .local_rng();
    let mut sequences = Vec::with_capacity(specs.len());
    for spec in specs {
        let length = usize::try_from(spec.length).map_err(|_| resource_error(spec.length))?;
        if length == 0 {
            return Err(configuration_error(
                "reference.sequences.length",
                format!("{} has zero length", spec.id),
            ));
        }
        let mut bases = Vec::with_capacity(length);
        for _ in 0..length {
            let unit = u64_to_unit_interval(rng.next_u64());
            bases.push(if unit < unknown_base_fraction {
                b'N'
            } else {
                select_base(
                    (unit - unknown_base_fraction) / (1.0 - unknown_base_fraction),
                    &cumulative,
                )
            });
        }
        // Se consume el búfer recién construido en lugar de copiarlo:
        // `try_from_bytes` duplicaría la secuencia entera y el pico pasaría al
        // doble del genoma generado.
        let sequence = DnaSequence::try_from_owned_bytes(bases).map_err(|error| {
            internal_error(format!(
                "generated A/C/G/T sequence was rejected: {error:?}"
            ))
        })?;
        let id = SequenceId::parse(spec.id, reference_id)?;
        let header = FastaHeader::parse(spec.id.as_bytes()).map_err(|error| {
            internal_error(format!(
                "validated sequence ID was rejected as header: {error:?}"
            ))
        })?;
        sequences.push(
            ReferenceSequence::from_owned(id, header, spec.class, sequence).map_err(|error| {
                internal_error(format!("generated sequence construction failed: {error:?}"))
            })?,
        );
    }
    Reference::new(
        ReferenceId::parse(reference_id, "execution")?,
        ReferenceOrigin::Generated,
        sequences,
    )
    .map_err(|error| {
        internal_error(format!(
            "generated reference construction failed: {error:?}"
        ))
    })
}

fn validate_unknown_base_fraction(value: f64) -> Result<f64, Diagnostic> {
    if !value.is_finite() || !(0.0..=1.0).contains(&value) {
        return Err(configuration_error(
            "reference.composition.N",
            "N probability must be finite and in [0,1]",
        ));
    }
    Ok(value)
}

fn validate_composition(composition: &[BaseFrequency]) -> Result<[(u8, f64); 4], Diagnostic> {
    if composition.len() != 4 {
        return Err(configuration_error(
            "reference.composition",
            "the effective composition must contain A, C, G and T exactly once",
        ));
    }
    let mut ordered = [(b'A', None), (b'C', None), (b'G', None), (b'T', None)];
    for item in composition {
        let base = item.base;
        let Some((_, slot)) = ordered.iter_mut().find(|(candidate, _)| *candidate == base) else {
            return Err(configuration_error(
                "reference.composition",
                format!("unsupported generated base 0x{:02x}", item.base),
            ));
        };
        if slot.replace(item.probability).is_some() {
            return Err(configuration_error(
                "reference.composition",
                format!("duplicate generated base {}", char::from(base)),
            ));
        }
    }
    let mut probabilities = [0.0_f64; 4];
    for (index, (base, probability)) in ordered.into_iter().enumerate() {
        let probability = probability.ok_or_else(|| {
            configuration_error(
                "reference.composition",
                format!("missing generated base {}", char::from(base)),
            )
        })?;
        probabilities[index] = probability;
    }
    let probabilities = normalize_composition_probabilities(probabilities)
        .map_err(|detail| configuration_error("reference.composition", detail))?;
    let mut cumulative = 0.0;
    let mut result = [(b'A', 0.0), (b'C', 0.0), (b'G', 0.0), (b'T', 1.0)];
    for (index, (probability, (_, threshold))) in
        probabilities.into_iter().zip(result.iter_mut()).enumerate()
    {
        cumulative += probability;
        *threshold = if index == 3 { 1.0 } else { cumulative };
    }
    Ok(result)
}

pub(crate) fn normalize_composition_probabilities(
    mut values: [f64; 4],
) -> Result<[f64; 4], String> {
    if values
        .iter()
        .any(|value| !value.is_finite() || !(0.0..=1.0).contains(value))
    {
        return Err("each A/C/G/T probability must be finite and in [0,1]".to_owned());
    }
    let sum = neumaier(values);
    let delta = (sum - 1.0).abs();
    if !sum.is_finite() || sum <= 0.0 || delta > COMPOSITION_TOLERANCE {
        return Err(format!(
            "A/C/G/T sum {sum} differs from 1 by more than {COMPOSITION_TOLERANCE}"
        ));
    }
    for value in &mut values {
        *value /= sum;
    }
    let normalized_sum = neumaier(values);
    if !normalized_sum.is_finite() || (normalized_sum - 1.0).abs() > COMPOSITION_TOLERANCE {
        return Err("A/C/G/T composition cannot be normalized finitely".to_owned());
    }
    Ok(values)
}

fn neumaier(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0_f64;
    let mut correction = 0.0_f64;
    for value in values {
        let next = sum + value;
        correction += if sum.abs() >= value.abs() {
            (sum - next) + value
        } else {
            (value - next) + sum
        };
        sum = next;
    }
    sum + correction
}

fn u64_to_unit_interval(value: u64) -> f64 {
    const DENOMINATOR: f64 = (1_u64 << 53) as f64;
    ((value >> 11) as f64) / DENOMINATOR
}

fn select_base(value: f64, cumulative: &[(u8, f64); 4]) -> u8 {
    cumulative
        .iter()
        .find(|(_, threshold)| value < *threshold)
        .map_or(b'T', |(base, _)| *base)
}

fn configuration_error(field: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-008",
        "La configuración contiene una combinación de campos incompatible.",
    )
    .with("path", "effective-config")
    .with("fields", field)
    .with("detail", detail.into())
}

fn resource_error(length: u64) -> Diagnostic {
    Diagnostic::resource(
        "GPS-RES-002",
        "Los recursos disponibles no permiten completar el plan sin degradar resultados obligatorios.",
    )
    .with("resource", "address_space")
    .with("required", length.to_string())
    .with("available", usize::MAX.to_string())
}

fn internal_error(detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "reference_generation")
    .with("operation", "generate_reference")
    .with("detail", detail.into())
}
