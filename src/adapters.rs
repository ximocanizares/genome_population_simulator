//! Adaptadores de configuración, formatos, persistencia y procesos externos.

use std::path::Path;
use std::sync::Arc;

use toml::de::{DeTable, DeValue};

use crate::configuration::{
    AdmittedToml, ConfigKind, DeclaredConfig, DeclaredEntity, DeclaredFrequencyCount,
    DeclaredFrequencySet, DeclaredFrequencyValue, DeclaredGeneralConfig, DeclaredGeneralParts,
    DeclaredNumber, DeclaredReference, DeclaredReferenceSource, DeclaredScenarioConfig,
    DeclaredSequence, ExpectedConfigKind, FrequencyMode, Located, Route, ScenarioLink,
    SequenceClass, SourceLocation, semantic_error,
};
use crate::diagnostics::Diagnostic;
use crate::rules::{
    ParameterKey, PopulationSelector, RegionSelector, Rule, RuleProvenance, RuleValue, SelectorSet,
    SpatialSelector, TemporalSelector,
};

/// Admite la sintaxis y la cabecera de un TOML antes de interpretar su
/// contenido semántico.
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
        let location = SourceLocation::from_span(&bytes, span.clone());
        let detail = error.to_string();
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
            .with("location", location.to_string())
            .with("detail", detail);
        if code == "GPS-CONF-007" {
            let field = source.get(span.clone()).unwrap_or("<unknown>").trim();
            let first = source.find(field).unwrap_or(span.start);
            diagnostic = diagnostic.with("field", field).with(
                "locations",
                format!(
                    "{}; {}",
                    SourceLocation::from_offset(&bytes, first),
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
    let version = version_value.get_ref().as_str();
    if version != Some("1") {
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

/// Construye tipos declarados propios después de que la frontera de clase y
/// versión haya terminado satisfactoriamente.
pub fn parse_declared_config(admitted: AdmittedToml) -> Result<DeclaredConfig, Diagnostic> {
    let source = std::str::from_utf8(admitted.bytes()).map_err(|error| {
        semantic_error(&admitted, "document", format!("UTF-8 inválido: {error}"))
    })?;
    let parsed = DeTable::parse(source)
        .map_err(|error| semantic_error(&admitted, "document", error.to_string()))?;
    match admitted.kind() {
        ConfigKind::General => parse_general(admitted.clone(), parsed.get_ref()),
        ConfigKind::Scenario => parse_scenario(admitted.clone(), parsed.get_ref()),
    }
}

fn parse_general(
    admitted: AdmittedToml,
    table: &DeTable<'_>,
) -> Result<DeclaredConfig, Diagnostic> {
    check_unknown(
        &admitted,
        table,
        &[
            "config_kind",
            "schema_version",
            "route",
            "seed",
            "normalization_tolerance",
            "reference",
            "scenarios",
            "events",
            "individuals",
            "populations",
            "subpopulations",
            "phases",
            "groups",
            "blocks",
            "rules",
            "frequency_sets",
        ],
    )?;
    let route = required_route(&admitted, table, "route")?;
    let seed = optional_located_string(&admitted, table, "seed")?;
    let normalization_tolerance = optional_number(&admitted, table, "normalization_tolerance")?;
    let reference = parse_reference(&admitted, required_table(&admitted, table, "reference")?)?;
    let scenarios = parse_scenario_links(&admitted, table)?;
    let events = parse_id_list(&admitted, table, "events")?;
    let individuals = parse_id_list(&admitted, table, "individuals")?;
    let populations = parse_id_list(&admitted, table, "populations")?;
    let subpopulations = parse_id_list(&admitted, table, "subpopulations")?;
    let phases = parse_id_list(&admitted, table, "phases")?;
    let groups = parse_id_list(&admitted, table, "groups")?;
    let blocks = parse_id_list(&admitted, table, "blocks")?;
    let rules = parse_rules(&admitted, table)?;
    let frequency_sets = parse_frequency_sets(&admitted, table)?;

    Ok(DeclaredConfig::General(Box::new(
        DeclaredGeneralConfig::new(DeclaredGeneralParts {
            source: admitted,
            route,
            seed,
            normalization_tolerance,
            reference,
            scenarios,
            events,
            individuals,
            populations,
            subpopulations,
            phases,
            groups,
            blocks,
            rules,
            frequency_sets,
        }),
    )))
}

fn parse_scenario(
    admitted: AdmittedToml,
    table: &DeTable<'_>,
) -> Result<DeclaredConfig, Diagnostic> {
    check_unknown(
        &admitted,
        table,
        &[
            "config_kind",
            "schema_version",
            "id",
            "route",
            "phases",
            "events",
            "crosses",
        ],
    )?;
    let id = required_string(&admitted, table, "id")?;
    let route = required_route(&admitted, table, "route")?;
    let phases = parse_id_list(&admitted, table, "phases")?;
    let events = parse_id_list(&admitted, table, "events")?;
    let crosses = parse_id_list(&admitted, table, "crosses")?;
    Ok(DeclaredConfig::Scenario(Box::new(
        DeclaredScenarioConfig::new(admitted, id, route, phases, events, crosses),
    )))
}

fn parse_reference(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<DeclaredReference, Diagnostic> {
    check_unknown(
        admitted,
        table,
        &[
            "kind",
            "id",
            "default_ploidy",
            "fasta",
            "sequences",
            "composition",
        ],
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
            let ploidy = required_u8(admitted, table, "default_ploidy")?;
            let composition = parse_frequency_values(
                admitted,
                required_array(admitted, table, "composition")?,
                "reference.composition",
            )?;
            DeclaredReferenceSource::Simulated {
                default_ploidy: ploidy,
                composition,
            }
        }
        "supplied" => {
            if table.contains_key("default_ploidy") || table.contains_key("composition") {
                return Err(semantic_error(
                    admitted,
                    "reference.kind,reference.default_ploidy,reference.composition",
                    "una referencia suministrada no admite campos de simulación",
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
    let sequence_values = required_array(admitted, table, "sequences")?;
    let mut sequences = Vec::with_capacity(sequence_values.len());
    for value in sequence_values {
        let sequence = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "reference.sequences", value, "table"))?;
        check_unknown(admitted, sequence, &["id", "class", "length", "ploidy"])?;
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
            id: required_string(admitted, sequence, "id")?,
            class: Located::new(class, class_text.location().clone()),
            length: required_u64(admitted, sequence, "length")?,
            ploidy: optional_u8(admitted, sequence, "ploidy")?,
        });
    }
    Ok(DeclaredReference {
        id,
        source,
        sequences,
    })
}

fn parse_rules(admitted: &AdmittedToml, table: &DeTable<'_>) -> Result<Vec<Rule>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "rules")? else {
        return Ok(Vec::new());
    };
    let mut result = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "rules", value, "table"))?;
        check_unknown(
            admitted,
            entry,
            &[
                "id",
                "parameter",
                "value",
                "individual",
                "population",
                "subpopulation",
                "phase",
                "generation",
                "reference",
                "sequence",
                "region_start",
                "region_end",
                "locus",
                "group",
                "block",
                "technical_event",
            ],
        )?;
        let id = required_string(admitted, entry, "id")?;
        if !ids.insert(id.value().clone()) {
            return Err(semantic_error(
                admitted,
                "rules.id",
                format!("identificador repetido: {}", id.value()),
            ));
        }
        let parameter_text = required_string(admitted, entry, "parameter")?;
        let parameter = ParameterKey::parse(parameter_text.value()).ok_or_else(|| {
            semantic_error(
                admitted,
                "rules.parameter",
                format!("parámetro desconocido: {}", parameter_text.value()),
            )
        })?;
        let selectors = parse_selectors(admitted, entry)?;
        let rule_value = required_rule_value(admitted, entry, "value", parameter)?;
        result.push(Rule::new(
            id.value(),
            parameter,
            rule_value,
            selectors,
            RuleProvenance {
                path: admitted.path().display().to_string(),
                location: location(admitted, value.span()).to_string(),
            },
        )?);
    }
    Ok(result)
}

fn parse_selectors(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<SelectorSet, Diagnostic> {
    let individual = optional_string(admitted, table, "individual")?;
    let population = optional_string(admitted, table, "population")?;
    let subpopulation = optional_string(admitted, table, "subpopulation")?;
    let population = match (population, subpopulation) {
        (Some(_), Some(_)) => {
            return Err(semantic_error(
                admitted,
                "rules.population,rules.subpopulation",
                "los niveles population y subpopulation son excluyentes",
            ));
        }
        (Some(id), None) => Some(PopulationSelector::Population(id)),
        (None, Some(id)) => Some(PopulationSelector::Subpopulation(id)),
        (None, None) => None,
    };
    let phase = optional_string(admitted, table, "phase")?;
    let generation = optional_u32(admitted, table, "generation")?;
    let temporal = match (phase, generation) {
        (Some(_), Some(_)) => {
            return Err(semantic_error(
                admitted,
                "rules.phase,rules.generation",
                "los niveles phase y generation son excluyentes",
            ));
        }
        (Some(id), None) => Some(TemporalSelector::Phase(id)),
        (None, Some(value)) => Some(TemporalSelector::Generation(value)),
        (None, None) => None,
    };

    let reference = optional_string(admitted, table, "reference")?;
    let sequence = optional_string(admitted, table, "sequence")?;
    let region_start = optional_u64(admitted, table, "region_start")?;
    let region_end = optional_u64(admitted, table, "region_end")?;
    let locus = optional_string(admitted, table, "locus")?;
    let spatial = match (reference, sequence, region_start, region_end, locus) {
        (None, None, None, None, None) => None,
        (Some(reference), Some(sequence), None, None, None) => Some(SpatialSelector::Sequence {
            reference,
            sequence,
        }),
        (Some(reference), Some(sequence), Some(start), Some(end), None) => {
            Some(SpatialSelector::Region(RegionSelector {
                reference,
                sequence,
                start,
                end,
            }))
        }
        (Some(reference), Some(sequence), None, None, Some(locus)) => {
            Some(SpatialSelector::Locus {
                reference,
                sequence,
                locus,
            })
        }
        _ => {
            return Err(semantic_error(
                admitted,
                "rules.reference,rules.sequence,rules.region_start,rules.region_end,rules.locus",
                "el selector espacial debe declarar referencia y secuencia, y una región completa o un locus",
            ));
        }
    };

    Ok(SelectorSet {
        individual,
        population,
        temporal,
        spatial,
        group: optional_string(admitted, table, "group")?,
        block: optional_string(admitted, table, "block")?,
        technical_event: optional_string(admitted, table, "technical_event")?,
    })
}

fn required_rule_value(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
    parameter: ParameterKey,
) -> Result<RuleValue, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    if let Some(text) = value.get_ref().as_str() {
        return Ok(RuleValue::Text(text.to_owned()));
    }
    if let Some(boolean) = value.get_ref().as_bool() {
        return Ok(RuleValue::Boolean(boolean));
    }
    if let Some(integer) = value.get_ref().as_integer() {
        if integer.radix() == 10
            && !matches!(
                parameter,
                ParameterKey::Ploidy | ParameterKey::OutputEnabled
            )
        {
            return integer
                .as_str()
                .parse::<f64>()
                .map(RuleValue::Float)
                .map_err(|_| wrong_type(admitted, field, value, "finite f64"));
        }
        return i64::from_str_radix(integer.as_str(), integer.radix())
            .map(RuleValue::Integer)
            .map_err(|_| wrong_type(admitted, field, value, "signed 64-bit integer"));
    }
    if let Some(float) = value.get_ref().as_float() {
        return float
            .as_str()
            .parse::<f64>()
            .map(RuleValue::Float)
            .map_err(|_| wrong_type(admitted, field, value, "finite f64"));
    }
    Err(wrong_type(
        admitted,
        field,
        value,
        "string, boolean, integer or float",
    ))
}

fn parse_frequency_sets(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Vec<DeclaredFrequencySet>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "frequency_sets")? else {
        return Ok(Vec::new());
    };
    let mut result = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "frequency_sets", value, "table"))?;
        check_unknown(
            admitted,
            entry,
            &["id", "mode", "tolerance", "values", "counts"],
        )?;
        let id = required_string(admitted, entry, "id")?;
        if !ids.insert(id.value().clone()) {
            return Err(semantic_error(
                admitted,
                "frequency_sets.id",
                format!("identificador repetido: {}", id.value()),
            ));
        }
        let mode_text = required_string(admitted, entry, "mode")?;
        let mode = match mode_text.value().as_str() {
            "probabilistic" => FrequencyMode::Probabilistic,
            "controlled" => FrequencyMode::Controlled,
            other => {
                return Err(semantic_error(
                    admitted,
                    "frequency_sets.mode",
                    format!("modo desconocido: {other}"),
                ));
            }
        };
        let tolerance = optional_number(admitted, entry, "tolerance")?;
        let declared_values = optional_array(admitted, entry, "values")?;
        let declared_counts = optional_array(admitted, entry, "counts")?;
        let (frequency_values, counts) = match (
            mode,
            declared_values,
            declared_counts,
            tolerance.as_ref(),
        ) {
            (FrequencyMode::Probabilistic, Some(values), None, _) => (
                parse_frequency_values(admitted, values, "frequency_sets.values")?,
                Vec::new(),
            ),
            (FrequencyMode::Controlled, None, Some(counts), None) => (
                Vec::new(),
                parse_frequency_counts(admitted, counts, "frequency_sets.counts")?,
            ),
            _ => {
                return Err(semantic_error(
                    admitted,
                    "frequency_sets.mode,frequency_sets.values,frequency_sets.counts,frequency_sets.tolerance",
                    "los modos probabilístico y controlado son excluyentes y exigen su forma propia",
                ));
            }
        };
        result.push(DeclaredFrequencySet {
            id,
            mode: Located::new(mode, mode_text.location().clone()),
            tolerance,
            values: frequency_values,
            counts,
        });
    }
    Ok(result)
}

fn parse_frequency_values(
    admitted: &AdmittedToml,
    values: &[toml::Spanned<DeValue<'_>>],
    field: &str,
) -> Result<Vec<DeclaredFrequencyValue>, Diagnostic> {
    let mut result = Vec::with_capacity(values.len());
    let mut identities = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, field, value, "table"))?;
        check_unknown(admitted, entry, &["id", "category", "value"])?;
        let id = required_string(admitted, entry, "id")?;
        let category = required_string(admitted, entry, "category")?;
        if !identities.insert((category.value().clone(), id.value().clone())) {
            return Err(semantic_error(
                admitted,
                field,
                format!(
                    "identidad y categoría repetidas: {}:{}",
                    category.value(),
                    id.value()
                ),
            ));
        }
        result.push(DeclaredFrequencyValue {
            id,
            category,
            value: required_number(admitted, entry, "value")?,
        });
    }
    Ok(result)
}

fn parse_frequency_counts(
    admitted: &AdmittedToml,
    values: &[toml::Spanned<DeValue<'_>>],
    field: &str,
) -> Result<Vec<DeclaredFrequencyCount>, Diagnostic> {
    let mut result = Vec::with_capacity(values.len());
    let mut identities = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, field, value, "table"))?;
        check_unknown(admitted, entry, &["id", "category", "count"])?;
        let id = required_string(admitted, entry, "id")?;
        let category = required_string(admitted, entry, "category")?;
        if !identities.insert((category.value().clone(), id.value().clone())) {
            return Err(semantic_error(
                admitted,
                field,
                format!(
                    "identidad y categoría repetidas: {}:{}",
                    category.value(),
                    id.value()
                ),
            ));
        }
        result.push(DeclaredFrequencyCount {
            id,
            category,
            count: required_u64(admitted, entry, "count")?,
        });
    }
    Ok(result)
}

fn parse_scenario_links(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Vec<ScenarioLink>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "scenarios")? else {
        return Ok(Vec::new());
    };
    let mut result = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "scenarios", value, "table"))?;
        check_unknown(admitted, entry, &["id", "path", "route"])?;
        let id = required_string(admitted, entry, "id")?;
        if !ids.insert(id.value().clone()) {
            return Err(semantic_error(
                admitted,
                "scenarios.id",
                format!("identificador repetido: {}", id.value()),
            ));
        }
        result.push(ScenarioLink {
            id,
            path: required_string(admitted, entry, "path")?,
            route: required_route(admitted, entry, "route")?,
        });
    }
    Ok(result)
}

fn parse_id_list(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Vec<DeclaredEntity>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, field)? else {
        return Ok(Vec::new());
    };
    let mut result = Vec::with_capacity(values.len());
    let mut ids = std::collections::BTreeSet::new();
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, field, value, "table"))?;
        check_unknown(admitted, entry, &["id"])?;
        let id = required_string(admitted, entry, "id")?;
        if !ids.insert(id.value().clone()) {
            return Err(semantic_error(
                admitted,
                format!("{field}.id"),
                format!("identificador repetido: {}", id.value()),
            ));
        }
        result.push(DeclaredEntity::new(id));
    }
    Ok(result)
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

fn optional_string(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<String>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_string(admitted, table, field).map(|value| value.value().clone()))
        .transpose()
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

fn required_route(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Located<Route>, Diagnostic> {
    let declared = required_string(admitted, table, field)?;
    Route::parse(declared.value())
        .map(|route| Located::new(route, declared.location().clone()))
        .ok_or_else(|| {
            semantic_error(
                admitted,
                field,
                format!("vía desconocida: {}", declared.value()),
            )
        })
}

fn required_u64(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Located<u64>, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    parse_u64(value.get_ref()).map_or_else(
        || Err(wrong_type(admitted, field, value, "non-negative integer")),
        |number| Ok(Located::new(number, location(admitted, value.span()))),
    )
}

fn optional_u64(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<u64>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_u64(admitted, table, field).map(|value| *value.value()))
        .transpose()
}

fn optional_u32(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<u32>, Diagnostic> {
    optional_u64(admitted, table, field)?
        .map(|value| {
            u32::try_from(value)
                .map_err(|_| semantic_error(admitted, field, "el entero no cabe en u32"))
        })
        .transpose()
}

fn required_u8(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Located<u8>, Diagnostic> {
    let value = required_u64(admitted, table, field)?;
    u8::try_from(*value.value())
        .map(|number| Located::new(number, value.location().clone()))
        .map_err(|_| semantic_error(admitted, field, "el entero no cabe en u8"))
}

fn optional_u8(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<Located<u8>>, Diagnostic> {
    table
        .contains_key(field)
        .then(|| required_u8(admitted, table, field))
        .transpose()
}

fn required_number(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<DeclaredNumber, Diagnostic> {
    let value = required_value(admitted, table, field)?;
    number_from_value(admitted, field, value)
}

fn optional_number(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<DeclaredNumber>, Diagnostic> {
    table
        .get(field)
        .map(|value| number_from_value(admitted, field, value))
        .transpose()
}

fn number_from_value(
    admitted: &AdmittedToml,
    field: &str,
    value: &toml::Spanned<DeValue<'_>>,
) -> Result<DeclaredNumber, Diagnostic> {
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

fn parse_u64(value: &DeValue<'_>) -> Option<u64> {
    let integer = value.as_integer()?;
    u64::from_str_radix(integer.as_str(), integer.radix()).ok()
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
