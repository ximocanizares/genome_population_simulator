//! Aritmética decimal exacta para fronteras que no admiten redondeo binario.

use std::cmp::Ordering;

use crate::diagnostics::Diagnostic;

/// Decimal finito no negativo conservado como coeficiente decimal y potencia de diez.
///
/// Esta representación no impone el límite artificial de `u128`: el tamaño admitido
/// queda limitado únicamente por el propio literal TOML y por los resultados `u64`
/// del dominio.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct ExactDecimal {
    coefficient: String,
    exponent10: i64,
}

impl ExactDecimal {
    /// Analiza la sintaxis decimal que TOML admite para enteros y flotantes.
    pub fn parse(field: &str, literal: &str) -> Result<Self, Diagnostic> {
        let ParsedTomlDecimal {
            negative,
            integer,
            fraction,
            exponent,
        } = parse_toml_decimal(literal).ok_or_else(|| invalid_decimal(field, literal))?;
        if negative {
            return Err(invalid_decimal(field, literal));
        }
        let coefficient = format!("{integer}{fraction}");
        let coefficient = coefficient.trim_start_matches('0');
        let coefficient = if coefficient.is_empty() {
            "0".to_owned()
        } else {
            coefficient.to_owned()
        };
        let fraction_digits =
            i64::try_from(fraction.len()).map_err(|_| invalid_decimal(field, literal))?;
        let exponent10 = exponent
            .checked_sub(fraction_digits)
            .ok_or_else(|| invalid_decimal(field, literal))?;
        Ok(Self {
            coefficient,
            exponent10,
        })
    }

    /// Calcula `floor(S * (100 - p) / 100)` sin pasar por `f64`.
    pub fn retained_after_percent(&self, eligible_bases: u64) -> Result<u64, Diagnostic> {
        if self.coefficient == "0" || self.compare_with_hundred() != Ordering::Less {
            return Err(percent_range(self));
        }

        // floor(S * (1 - p/100)) = S - ceil(S * p/100).
        let product = multiply_decimal_by_u64(&self.coefficient, eligible_bases);
        let scale = self
            .exponent10
            .checked_sub(2)
            .ok_or_else(|| invalid_decimal("reduction_percent", "exponent overflow"))?;
        let removed = if scale >= 0 {
            let zero_count = usize::try_from(scale)
                .map_err(|_| invalid_decimal("reduction_percent", "exponent overflow"))?;
            let mut integer = product;
            integer
                .try_reserve(zero_count)
                .map_err(|_| invalid_decimal("reduction_percent", "allocation failed"))?;
            integer.extend(std::iter::repeat_n('0', zero_count));
            parse_domain_integer(&integer)?
        } else {
            let decimal_places = usize::try_from(scale.unsigned_abs())
                .map_err(|_| invalid_decimal("reduction_percent", "exponent overflow"))?;
            ceil_divide_by_power_of_ten(&product, decimal_places)?
        };
        eligible_bases
            .checked_sub(removed)
            .ok_or_else(|| percent_range(self))
    }

    /// Calcula `floor(multiplicand * self)` sin un producto entero intermedio
    /// de ancho fijo. El llamador conserva la autoridad sobre el intervalo
    /// admitido; esta operación solo aplica la representación decimal exacta.
    pub fn floor_product(&self, multiplicand: u64) -> Result<u64, Diagnostic> {
        let product = multiply_decimal_by_u64(&self.coefficient, multiplicand);
        if self.exponent10 >= 0 {
            let zero_count = usize::try_from(self.exponent10)
                .map_err(|_| invalid_decimal("exact_decimal", "exponent overflow"))?;
            let mut integer = product;
            integer
                .try_reserve(zero_count)
                .map_err(|_| invalid_decimal("exact_decimal", "allocation failed"))?;
            integer.extend(std::iter::repeat_n('0', zero_count));
            parse_domain_integer(&integer)
        } else {
            let decimal_places = usize::try_from(self.exponent10.unsigned_abs())
                .map_err(|_| invalid_decimal("exact_decimal", "exponent overflow"))?;
            floor_divide_by_power_of_ten(&product, decimal_places)
        }
    }

    fn compare_with_hundred(&self) -> Ordering {
        if self.coefficient == "0" {
            return Ordering::Less;
        }
        let Ok(coefficient_digits) = i64::try_from(self.coefficient.len()) else {
            return Ordering::Greater;
        };
        let Some(integer_digits) = coefficient_digits.checked_add(self.exponent10) else {
            return if self.exponent10.is_positive() {
                Ordering::Greater
            } else {
                Ordering::Less
            };
        };
        match integer_digits.cmp(&3) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => {
                let mut first_three = self.coefficient.bytes().take(3).collect::<Vec<_>>();
                first_three.resize(3, b'0');
                match first_three.as_slice().cmp(b"100") {
                    Ordering::Equal => {
                        if self.coefficient.bytes().skip(3).any(|digit| digit != b'0') {
                            Ordering::Greater
                        } else {
                            Ordering::Equal
                        }
                    }
                    ordering => ordering,
                }
            }
        }
    }
}

#[derive(Debug, Eq, PartialEq)]
struct ParsedTomlDecimal {
    negative: bool,
    integer: String,
    fraction: String,
    exponent: i64,
}

/// Reconoce primero la gramática numérica decimal de TOML. Los separadores se
/// eliminan únicamente después de comprobar que están entre dos dígitos.
fn parse_toml_decimal(literal: &str) -> Option<ParsedTomlDecimal> {
    let bytes = literal.as_bytes();
    let mut cursor = 0_usize;
    let negative = match bytes.first() {
        Some(b'+') => {
            cursor = 1;
            false
        }
        Some(b'-') => {
            cursor = 1;
            true
        }
        _ => false,
    };
    let integer_end = bytes[cursor..]
        .iter()
        .position(|byte| matches!(byte, b'.' | b'e' | b'E'))
        .map_or(bytes.len(), |offset| cursor + offset);
    let integer = toml_digit_group(&bytes[cursor..integer_end])?;
    if integer.len() > 1 && integer.starts_with('0') {
        return None;
    }
    cursor = integer_end;

    let mut fraction = String::new();
    if bytes.get(cursor) == Some(&b'.') {
        cursor += 1;
        let fraction_end = bytes[cursor..]
            .iter()
            .position(|byte| matches!(byte, b'e' | b'E'))
            .map_or(bytes.len(), |offset| cursor + offset);
        fraction = toml_digit_group(&bytes[cursor..fraction_end])?;
        cursor = fraction_end;
    }

    let exponent = if matches!(bytes.get(cursor), Some(b'e' | b'E')) {
        cursor += 1;
        let exponent_negative = match bytes.get(cursor) {
            Some(b'+') => {
                cursor += 1;
                false
            }
            Some(b'-') => {
                cursor += 1;
                true
            }
            _ => false,
        };
        let digits = toml_digit_group(&bytes[cursor..])?;
        cursor = bytes.len();
        let magnitude = digits.parse::<i64>().ok()?;
        if exponent_negative {
            magnitude.checked_neg()?
        } else {
            magnitude
        }
    } else {
        0
    };
    (cursor == bytes.len()).then_some(ParsedTomlDecimal {
        negative,
        integer,
        fraction,
        exponent,
    })
}

fn toml_digit_group(bytes: &[u8]) -> Option<String> {
    if bytes.is_empty() || !bytes.first()?.is_ascii_digit() || !bytes.last()?.is_ascii_digit() {
        return None;
    }
    let mut previous_was_digit = false;
    let mut digits = String::with_capacity(bytes.len());
    for (index, byte) in bytes.iter().copied().enumerate() {
        if byte.is_ascii_digit() {
            digits.push(char::from(byte));
            previous_was_digit = true;
        } else if byte == b'_'
            && previous_was_digit
            && bytes.get(index + 1).is_some_and(u8::is_ascii_digit)
        {
            previous_was_digit = false;
        } else {
            return None;
        }
    }
    Some(digits)
}

fn multiply_decimal_by_u64(digits: &str, multiplier: u64) -> String {
    if multiplier == 0 || digits == "0" {
        return "0".to_owned();
    }
    let mut reversed = Vec::with_capacity(digits.len() + 20);
    let mut carry = 0_u128;
    for digit in digits.bytes().rev() {
        let value = u128::from(digit - b'0') * u128::from(multiplier) + carry;
        reversed.push(u8::try_from(value % 10).expect("decimal digit") + b'0');
        carry = value / 10;
    }
    while carry > 0 {
        reversed.push(u8::try_from(carry % 10).expect("decimal digit") + b'0');
        carry /= 10;
    }
    reversed.reverse();
    String::from_utf8(reversed).expect("only ASCII decimal digits")
}

fn ceil_divide_by_power_of_ten(digits: &str, decimal_places: usize) -> Result<u64, Diagnostic> {
    if decimal_places >= digits.len() {
        return Ok(u64::from(digits.bytes().any(|digit| digit != b'0')));
    }
    let split = digits.len() - decimal_places;
    let quotient = parse_domain_integer(&digits[..split])?;
    let has_remainder = digits[split..].bytes().any(|digit| digit != b'0');
    quotient
        .checked_add(u64::from(has_remainder))
        .ok_or_else(|| invalid_decimal("reduction_percent", "result overflow"))
}

fn floor_divide_by_power_of_ten(digits: &str, decimal_places: usize) -> Result<u64, Diagnostic> {
    if decimal_places >= digits.len() {
        return Ok(0);
    }
    parse_domain_integer(&digits[..digits.len() - decimal_places])
}

fn parse_domain_integer(digits: &str) -> Result<u64, Diagnostic> {
    digits
        .parse::<u64>()
        .map_err(|_| invalid_decimal("reduction_percent", "result overflow"))
}

fn invalid_decimal(field: &str, literal: &str) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-019",
        "El porcentaje de reducción no es un decimal finito representable exactamente.",
    )
    .with("field", field)
    .with("found", literal)
    .with("accepted", "finite exact decimal in (0,100)")
}

fn percent_range(value: &ExactDecimal) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-020",
        "El porcentaje de reducción está fuera del intervalo abierto admitido.",
    )
    .with("field", "reduction_percent")
    .with(
        "found",
        format!("{}e{}", value.coefficient, value.exponent10),
    )
    .with("accepted", "0 < reduction_percent < 100")
}

#[cfg(test)]
mod tests {
    use super::ExactDecimal;

    #[test]
    fn equivalent_literals_produce_the_same_exact_target() {
        for literal in ["25", "25.0", "2.5e1", "2500e-2"] {
            let value = ExactDecimal::parse("reduction_percent", literal).unwrap();
            assert_eq!(value.retained_after_percent(101).unwrap(), 75);
        }
    }

    #[test]
    fn floors_once_after_the_exact_product() {
        let value = ExactDecimal::parse("reduction_percent", "33.333333333333333333").unwrap();
        assert_eq!(value.retained_after_percent(3).unwrap(), 2);
    }

    #[test]
    fn exact_boundary_does_not_round_through_f64() {
        let value = ExactDecimal::parse("reduction_percent", "10.0000000000000000001").unwrap();
        assert_eq!(value.retained_after_percent(10).unwrap(), 8);
        let rounded = "10.0000000000000000001".parse::<f64>().unwrap();
        assert_eq!((10.0 * (100.0 - rounded) / 100.0).floor() as u64, 9);
    }

    #[test]
    fn supports_precision_beyond_fixed_width_integers() {
        let tiny = ExactDecimal::parse("reduction_percent", "1e-1000").unwrap();
        assert_eq!(tiny.retained_after_percent(100).unwrap(), 99);
        let long = ExactDecimal::parse(
            "reduction_percent",
            "99.999999999999999999999999999999999999999999999999999999999999",
        )
        .unwrap();
        assert_eq!(long.retained_after_percent(u64::MAX).unwrap(), 0);
    }

    #[test]
    fn exact_floor_product_does_not_overflow_before_division() {
        let half = ExactDecimal::parse("coverage_fraction", "0.500000000000000000").unwrap();
        assert_eq!(half.floor_product(100).unwrap(), 50);
        let almost_one = ExactDecimal::parse("coverage_fraction", "0.999999999999999999").unwrap();
        assert_eq!(almost_one.floor_product(u64::MAX).unwrap(), u64::MAX - 19);
    }

    #[test]
    fn rejects_endpoints_negative_and_values_above_hundred() {
        for literal in ["0", "100", "100.0", "-1", "1e1000"] {
            let result = ExactDecimal::parse("reduction_percent", literal)
                .and_then(|value| value.retained_after_percent(100));
            assert!(result.is_err(), "literal unexpectedly accepted: {literal}");
        }
    }

    #[test]
    fn enforces_toml_decimal_grammar_before_removing_separators() {
        for literal in [
            "1_0",
            "+1_0.0_0",
            "1.2_5e+1",
            "9_9.999_999e-1",
            "30.0000000000000000001",
        ] {
            ExactDecimal::parse("reduction_percent", literal)
                .unwrap_or_else(|_| panic!("valid TOML decimal rejected: {literal}"));
        }
        for literal in [
            "_10", "10_", "1__0", "1_.0", "1._0", ".5", "5.", "00", "01", "+01", "1e", "1e+",
            "1e_2", "1e2_", "1e2__0", "++1", " 1", "1 ",
        ] {
            assert!(
                ExactDecimal::parse("reduction_percent", literal).is_err(),
                "invalid TOML decimal accepted: {literal}"
            );
        }
    }
}
