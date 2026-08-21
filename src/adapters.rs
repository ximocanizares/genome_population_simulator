//! Adaptador TOML cerrado para la configuración general mínima de F1.

use std::path::Path;
use std::sync::Arc;

use toml::de::{DeTable, DeValue};

use crate::configuration::{
    AdmittedToml, ConfigKind, DeclaredCompositionValue, DeclaredGeneralConfig, DeclaredNumber,
    DeclaredPopulationSnvs, DeclaredReducedReference, DeclaredReductionTarget, DeclaredReference,
    DeclaredReferenceSource, DeclaredSequence, ExpectedConfigKind, Located, SequenceClass,
    SourceLocation, semantic_error,
};
use crate::diagnostics::Diagnostic;

pub fn admit_toml_header(
    path: &Path,
    bytes: Vec<u8>,
    expected: ExpectedConfigKind,
) -> Result<AdmittedToml, Diagnostic> {
    let source = std::str::from_utf8(&bytes).map_err(|error| {
        Diagnostic::configuration(
            "GPS-CONF-001",
            "El documento de configuración no tiene sintaxis TOML válida.",
        )
        .with("path", path.display().to_string())
        .with(
            "location",
            SourceLocation::from_offset(&bytes, error.valid_up_to()).to_string(),
        )
        .with("detail", error.to_string())
    })?;
    let parsed = DeTable::parse(source).map_err(|error| {
        let span = error.span().unwrap_or(0..0);
        let (code, message) = if error.message().contains("duplicate key") {
            (
                "GPS-CONF-007",
                "La configuración repite una clave o un campo.",
            )
        } else {
            (
                "GPS-CONF-001",
                "El documento de configuración no tiene sintaxis TOML válida.",
            )
        };
        let mut diagnostic = Diagnostic::configuration(code, message)
            .with("path", path.display().to_string())
            .with(
                "location",
                SourceLocation::from_span(&bytes, span.clone()).to_string(),
            )
            .with("detail", error.to_string());
        if code == "GPS-CONF-007" {
            let field = source.get(span.clone()).unwrap_or("<unknown>").trim();
            diagnostic = diagnostic.with("field", field).with(
                "locations",
                format!(
                    "{}; {}",
                    SourceLocation::from_offset(&bytes, source.find(field).unwrap_or(span.start)),
                    SourceLocation::from_span(&bytes, span)
                ),
            );
        }
        diagnostic
    })?;
    let table = parsed.get_ref();
    let kind_value = table.get("config_kind").ok_or_else(|| {
        Diagnostic::configuration(
            "GPS-CONF-002",
            "Falta el discriminador de clase del documento.",
        )
        .with("path", path.display().to_string())
        .with("field", "config_kind")
    })?;
    let kind_text = kind_value.get_ref().as_str().ok_or_else(|| {
        invalid_kind(
            path,
            kind_value.get_ref(),
            expected,
            &bytes,
            kind_value.span(),
        )
    })?;
    let kind = ConfigKind::parse(kind_text)
        .filter(|kind| expected.accepts(*kind))
        .ok_or_else(|| {
            invalid_kind(
                path,
                kind_value.get_ref(),
                expected,
                &bytes,
                kind_value.span(),
            )
        })?;
    let version_value = table.get("schema_version").ok_or_else(|| {
        Diagnostic::configuration("GPS-CONF-004", "Falta la versión de esquema del documento.")
            .with("path", path.display().to_string())
            .with("field", "schema_version")
    })?;
    if version_value.get_ref().as_str() != Some("1") {
        return Err(Diagnostic::configuration(
            "GPS-CONF-005",
            "La versión de esquema declarada no está admitida.",
        )
        .with("path", path.display().to_string())
        .with("found", render_toml_value(version_value.get_ref()))
        .with("accepted", "1")
        .with(
            "location",
            SourceLocation::from_span(&bytes, version_value.span()).to_string(),
        ));
    }
    Ok(AdmittedToml::new(
        path.to_path_buf(),
        Arc::from(bytes),
        kind,
    ))
}

fn invalid_kind(
    path: &Path,
    found: &DeValue<'_>,
    expected: ExpectedConfigKind,
    bytes: &[u8],
    span: std::ops::Range<usize>,
) -> Diagnostic {
    Diagnostic::configuration(
        "GPS-CONF-003",
        "La clase declarada del documento no está admitida en esta frontera.",
    )
    .with("path", path.display().to_string())
    .with("field", "config_kind")
    .with("found", render_toml_value(found))
    .with("accepted", expected.as_str())
    .with(
        "location",
        SourceLocation::from_span(bytes, span).to_string(),
    )
}

fn render_toml_value(value: &DeValue<'_>) -> String {
    value
        .as_str()
        .map_or_else(|| value.type_str().to_owned(), ToOwned::to_owned)
}

pub fn parse_declared_config(admitted: AdmittedToml) -> Result<DeclaredGeneralConfig, Diagnostic> {
    let source = std::str::from_utf8(admitted.bytes()).map_err(|error| {
        semantic_error(&admitted, "document", format!("UTF-8 inválido: {error}"))
    })?;
    let parsed = DeTable::parse(source)
        .map_err(|error| semantic_error(&admitted, "document", error.to_string()))?;
    let table = parsed.get_ref();
    check_unknown(
        &admitted,
        table,
        &[
            "config_kind",
            "schema_version",
            "seed",
            "ploidy",
            "normalization_tolerance",
            "reference",
            "reduced_references",
            "population_snvs",
        ],
    )?;
    let seed = optional_located_string(&admitted, table, "seed")?;
    let ploidy = optional_located_u64(&admitted, table, "ploidy")?;
    let normalization_tolerance = optional_number(&admitted, table, "normalization_tolerance")?;
    let reference = parse_reference(&admitted, required_table(&admitted, table, "reference")?)?;
    let reduced = parse_reduced_references(&admitted, table)?;
    let population_snvs = parse_population_snvs(&admitted, table)?;
    Ok(DeclaredGeneralConfig::new(
        admitted,
        seed,
        ploidy,
        normalization_tolerance,
        reference,
        reduced,
        population_snvs,
    ))
}

fn parse_population_snvs(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Option<DeclaredPopulationSnvs>, Diagnostic> {
    let Some(value) = table.get("population_snvs") else {
        return Ok(None);
    };
    let entry = value
        .get_ref()
        .as_table()
        .ok_or_else(|| wrong_type(admitted, "population_snvs", value, "table"))?;
    check_unknown(admitted, entry, &["mode", "input_path", "metadata_path"])?;
    let mode = required_string(admitted, entry, "mode")?;
    match mode.value().as_str() {
        "supplied" => {
            if entry.contains_key("metadata_path") {
                return Err(semantic_error(
                    admitted,
                    "population_snvs.mode,population_snvs.metadata_path",
                    "el modo supplied no admite metadata_path",
                ));
            }
            Ok(Some(DeclaredPopulationSnvs::Supplied {
                input_path: required_non_empty_string(
                    admitted,
                    entry,
                    "input_path",
                    "population_snvs.input_path",
                )?,
            }))
        }
        "reuse" => {
            if entry.contains_key("input_path") {
                return Err(semantic_error(
                    admitted,
                    "population_snvs.mode,population_snvs.input_path",
                    "el modo reuse no admite input_path",
                ));
            }
            Ok(Some(DeclaredPopulationSnvs::Reuse {
                metadata_path: required_non_empty_string(
                    admitted,
                    entry,
                    "metadata_path",
                    "population_snvs.metadata_path",
                )?,
            }))
        }
        "generate" => {
            if entry.contains_key("input_path") || entry.contains_key("metadata_path") {
                return Err(semantic_error(
                    admitted,
                    "population_snvs",
                    "el modo generate no admite rutas en F2.2.1",
                ));
            }
            Ok(Some(DeclaredPopulationSnvs::Generate))
        }
        other => Err(semantic_error(
            admitted,
            "population_snvs.mode",
            format!("modo de matriz SNV desconocido: {other}"),
        )),
    }
}

fn required_non_empty_string(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
    diagnostic_field: &str,
) -> Result<Located<String>, Diagnostic> {
    let value = required_string(admitted, table, field)?;
    if value.value().is_empty() {
        return Err(
            semantic_error(admitted, diagnostic_field, "la ruta no puede estar vacía")
                .with("location", value.location().to_string()),
        );
    }
    Ok(value)
}

fn parse_reference(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<DeclaredReference, Diagnostic> {
    check_unknown(
        admitted,
        table,
        &["kind", "id", "fasta", "sequences", "composition"],
    )?;
    let id = required_string(admitted, table, "id")?;
    let kind = required_string(admitted, table, "kind")?;
    let source = match kind.value().as_str() {
        "simulated" => {
            if table.contains_key("fasta") {
                return Err(semantic_error(
                    admitted,
                    "reference.kind,reference.fasta",
                    "una referencia simulada no admite fasta",
                ));
            }
            let composition_table = required_table(admitted, table, "composition")?;
            check_unknown(admitted, composition_table, &["A", "C", "G", "T"])?;
            let composition = ["A", "C", "G", "T"]
                .into_iter()
                .map(|base| {
                    let value = required_number(admitted, composition_table, base)?;
                    Ok(DeclaredCompositionValue {
                        base: Located::new(base.to_owned(), value.location().clone()),
                        value,
                    })
                })
                .collect::<Result<Vec<_>, Diagnostic>>()?;
            DeclaredReferenceSource::Simulated { composition }
        }
        "supplied" => {
            if table.contains_key("composition") {
                return Err(semantic_error(
                    admitted,
                    "reference.kind,reference.composition",
                    "una referencia suministrada no admite composition",
                ));
            }
            DeclaredReferenceSource::Supplied {
                fasta: required_string(admitted, table, "fasta")?,
            }
        }
        other => {
            return Err(semantic_error(
                admitted,
                "reference.kind",
                format!("origen de referencia desconocido: {other}"),
            ));
        }
    };
    let values = required_array(admitted, table, "sequences")?;
    if values.is_empty() {
        return Err(semantic_error(
            admitted,
            "reference.sequences",
            "se exige al menos una secuencia",
        ));
    }
    let mut sequences = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let sequence = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "reference.sequences", value, "table"))?;
        check_unknown(admitted, sequence, &["id", "class", "length"])?;
        let sequence_id = required_string(admitted, sequence, "id")?;
        if !ids.insert(sequence_id.value().clone()) {
            return Err(semantic_error(
                admitted,
                "reference.sequences.id",
                format!("identificador repetido: {}", sequence_id.value()),
            ));
        }
        let class_text = required_string(admitted, sequence, "class")?;
        let class = match class_text.value().as_str() {
            "chromosome" => SequenceClass::Chromosome,
            "scaffold" => SequenceClass::Scaffold,
            other => {
                return Err(semantic_error(
                    admitted,
                    "reference.sequences.class",
                    format!("clase de secuencia desconocida: {other}"),
                ));
            }
        };
        sequences.push(DeclaredSequence {
            id: sequence_id,
            class: Located::new(class, class_text.location().clone()),
            length: required_u64(admitted, sequence, "length")?,
        });
    }
    Ok(DeclaredReference {
        id,
        source,
        sequences,
    })
}

fn parse_reduced_references(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Vec<DeclaredReducedReference>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "reduced_references")? else {
        return Ok(Vec::new());
    };
    let mut reductions = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "reduced_references", value, "table"))?;
        check_unknown(
            admitted,
            entry,
            &[
                "id",
                "source_reference_id",
                "region_size_bases",
                "target_bases",
                "reduction_percent",
            ],
        )?;
        let id = required_string(admitted, entry, "id")?;
        if !ids.insert(id.value().clone()) {
            return Err(semantic_error(
                admitted,
                "reduced_references.id",
                format!("identificador repetido: {}", id.value()),
            ));
        }
        let target = match (
            entry.contains_key("target_bases"),
            entry.contains_key("reduction_percent"),
        ) {
            (true, false) => {
                DeclaredReductionTarget::Bases(required_u64(admitted, entry, "target_bases")?)
            }
            (false, true) => DeclaredReductionTarget::ReductionPercent(required_number(
                admitted,
                entry,
                "reduction_percent",
            )?),
            _ => {
                return Err(semantic_error(
                    admitted,
                    "reduced_references.target_bases,reduced_references.reduction_percent",
                    "cada referencia reducida exige exactamente un objetivo",
                ));
            }
        };
        reductions.push(DeclaredReducedReference {
            id,
            source_reference_id: required_string(admitted, entry, "source_reference_id")?,
            region_size_bases: required_u64(admitted, entry, "region_size_bases")?,
            target,
        });
    }
    Ok(reductions)
}

fn check_unknown(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    accepted: &[&str],
) -> Result<(), Diagnostic> {
    if let Some(key) = table
        .keys()
        .find(|key| !accepted.contains(&key.get_ref().as_ref()))
    {
        return Err(Diagnostic::configuration(
            "GPS-CONF-006",
            "La configuración contiene un campo desconocido.",
        )
        .with("path", admitted.path().display().to_string())
        .with("field", key.get_ref().as_ref())
        .with("location", location(admitted, key.span()).to_string()));
    }
    Ok(())
}

fn required_value<'a>(
    admitted: &AdmittedToml,
    table: &'a DeTable<'a>,
    field: &str,
) -> Result<&'a toml::Spanned<DeValue<'a>>, Diagnostic> {
    table
        .get(field)
        .ok_or_else(|| semantic_error(admitted, field, "falta un campo obligatorio"))
}

fn required_table<'a>(
    admitted: &AdmittedToml,
    table: &'a DeTable<'a>,
    field: &str,
) -> Result<&'a DeTable<'a>, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    value
        .get_ref()
        .as_table()
        .ok_or_else(|| wrong_type(admitted, field, value, "table"))
}

fn required_array<'a>(
    admitted: &AdmittedToml,
    table: &'a DeTable<'a>,
    field: &str,
) -> Result<&'a [toml::Spanned<DeValue<'a>>], Diagnostic> {
    let value = required_value(admitted, table, field)?;
    value
        .get_ref()
        .as_array()
        .map(AsRef::as_ref)
        .ok_or_else(|| wrong_type(admitted, field, value, "array"))
}

fn optional_array<'a>(
    admitted: &AdmittedToml,
    table: &'a DeTable<'a>,
    field: &str,
) -> Result<Option<&'a [toml::Spanned<DeValue<'a>>]>, Diagnostic> {
    let Some(value) = table.get(field) else {
        return Ok(None);
    };
    value
        .get_ref()
        .as_array()
        .map(AsRef::as_ref)
        .map(Some)
        .ok_or_else(|| wrong_type(admitted, field, value, "array"))
}

fn required_string(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Located<String>, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    value
        .get_ref()
        .as_str()
        .map(|text| Located::new(text.to_owned(), location(admitted, value.span())))
        .ok_or_else(|| wrong_type(admitted, field, value, "string"))
}

fn optional_located_string(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<Located<String>>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_string(admitted, table, field))
        .transpose()
}

fn optional_located_u64(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<Located<u64>>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_u64(admitted, table, field))
        .transpose()
}

fn optional_number(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<DeclaredNumber>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_number(admitted, table, field))
        .transpose()
}

fn required_u64(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Located<u64>, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    let integer = value.get_ref().as_integer();
    let number =
        integer.and_then(|integer| u64::from_str_radix(integer.as_str(), integer.radix()).ok());
    number
        .map(|number| Located::new(number, location(admitted, value.span())))
        .ok_or_else(|| wrong_type(admitted, field, value, "non-negative integer"))
}

fn required_number(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<DeclaredNumber, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    let literal = if let Some(float) = value.get_ref().as_float() {
        float.as_str()
    } else if let Some(integer) = value.get_ref().as_integer()
        && integer.radix() == 10
    {
        integer.as_str()
    } else {
        return Err(wrong_type(admitted, field, value, "decimal number"));
    };
    Ok(DeclaredNumber::new(
        literal.to_owned(),
        location(admitted, value.span()),
    ))
}

fn wrong_type(
    admitted: &AdmittedToml,
    field: &str,
    value: &toml::Spanned<DeValue<'_>>,
    expected: &str,
) -> Diagnostic {
    semantic_error(
        admitted,
        field,
        format!(
            "se esperaba {expected}, se encontró {} en {}",
            value.get_ref().type_str(),
            location(admitted, value.span())
        ),
    )
}

fn location(admitted: &AdmittedToml, span: std::ops::Range<usize>) -> SourceLocation {
    SourceLocation::from_span(admitted.bytes(), span)
}
