//! Validación y normalización numérica reproducible de DTC-049.

use crate::configuration::{DeclaredFrequencyValue, DeclaredNumber};
use crate::diagnostics::Diagnostic;

pub const DEFAULT_NORMALIZATION_TOLERANCE: f64 = 1e-9;
pub const MAX_NORMALIZATION_TOLERANCE: f64 = 1e-6;
pub const NORMALIZATION_ALGORITHM: &str = "neumaier";
pub const NORMALIZATION_VERSION: &str = "1";
pub const CANONICAL_ORDER: &str = "id,category";

#[derive(Clone, Debug, PartialEq)]
pub struct NormalizedValue {
    pub id: String,
    pub category: String,
    pub before: f64,
    pub after: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct NormalizedCollection {
    pub parameter: String,
    pub tolerance: f64,
    pub canonical_order: &'static str,
    pub algorithm: &'static str,
    pub algorithm_version: &'static str,
    pub sum_before: f64,
    pub delta_before: f64,
    pub values: Vec<NormalizedValue>,
    pub sum_after: f64,
    pub delta_after: f64,
}

pub fn normalize_declared(
    parameter: &str,
    declared: &[DeclaredFrequencyValue],
    tolerance: Option<&DeclaredNumber>,
) -> Result<NormalizedCollection, Diagnostic> {
    let tolerance = match tolerance {
        Some(value) => parse_finite(value, parameter, "tolerance")?,
        None => DEFAULT_NORMALIZATION_TOLERANCE,
    };
    if !(0.0..=MAX_NORMALIZATION_TOLERANCE).contains(&tolerance) {
        return Err(numeric_value_error(
            parameter,
            "tolerance",
            tolerance,
            "finite [0,1e-6]",
        ));
    }

    let mut values = declared
        .iter()
        .map(|value| {
            let parsed = parse_finite(&value.value, parameter, value.id.value())?;
            if !(0.0..=1.0).contains(&parsed) {
                return Err(numeric_value_error(
                    parameter,
                    value.id.value(),
                    parsed,
                    "finite [0,1]",
                ));
            }
            Ok(NormalizedValue {
                id: value.id.value().clone(),
                category: value.category.value().clone(),
                before: parsed,
                after: parsed,
            })
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    values
        .sort_by(|first, second| (&first.id, &first.category).cmp(&(&second.id, &second.category)));

    let sum_before = neumaier(values.iter().map(|value| value.before));
    let delta_before = (sum_before - 1.0).abs();
    if !sum_before.is_finite()
        || sum_before <= 0.0
        || !delta_before.is_finite()
        || delta_before > tolerance
    {
        return Err(normalization_error(
            parameter,
            tolerance,
            sum_before,
            delta_before,
        ));
    }

    for value in &mut values {
        value.after = value.before / sum_before;
        if !value.after.is_finite() || !(0.0..=1.0).contains(&value.after) {
            return Err(normalization_error(
                parameter,
                tolerance,
                sum_before,
                delta_before,
            ));
        }
    }
    let sum_after = neumaier(values.iter().map(|value| value.after));
    let delta_after = (sum_after - 1.0).abs();
    if !sum_after.is_finite() || !delta_after.is_finite() || delta_after > tolerance {
        return Err(normalization_error(
            parameter,
            tolerance,
            sum_after,
            delta_after,
        ));
    }

    Ok(NormalizedCollection {
        parameter: parameter.to_owned(),
        tolerance,
        canonical_order: CANONICAL_ORDER,
        algorithm: NORMALIZATION_ALGORITHM,
        algorithm_version: NORMALIZATION_VERSION,
        sum_before,
        delta_before,
        values,
        sum_after,
        delta_after,
    })
}

#[must_use]
pub fn neumaier(values: impl IntoIterator<Item = f64>) -> f64 {
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

/// Representación decimal mínima cuya lectura recupera exactamente el mismo
/// `f64`, garantía proporcionada por la implementación estándar de `Display`.
#[must_use]
pub fn roundtrip_decimal(value: f64) -> String {
    value.to_string()
}

fn parse_finite(
    declared: &DeclaredNumber,
    parameter: &str,
    identity: &str,
) -> Result<f64, Diagnostic> {
    let parsed = declared.literal().parse::<f64>().map_err(|_| {
        numeric_value_error(
            parameter,
            identity,
            declared.literal(),
            "finite decimal f64",
        )
    })?;
    if parsed.is_finite() {
        Ok(parsed)
    } else {
        Err(numeric_value_error(
            parameter,
            identity,
            declared.literal(),
            "finite decimal f64",
        ))
    }
}

fn numeric_value_error(
    parameter: &str,
    identity: &str,
    found: impl ToString,
    accepted: &str,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-015",
        "Una probabilidad, frecuencia o tolerancia no es finita o está fuera de su intervalo.",
    )
    .with("parameter", parameter)
    .with("identity", identity)
    .with("found", found.to_string())
    .with("accepted", accepted)
}

fn normalization_error(parameter: &str, tolerance: f64, sum: f64, delta: f64) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-016",
        "Una colección de probabilidades o frecuencias no puede normalizarse con la tolerancia efectiva.",
    )
    .with("parameter", parameter)
    .with("tolerance", roundtrip_decimal(tolerance))
    .with("algorithm", NORMALIZATION_ALGORITHM)
    .with("sum", roundtrip_decimal(sum))
    .with("delta", roundtrip_decimal(delta))
}
