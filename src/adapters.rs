//! Adaptador TOML cerrado para la configuración general mínima de F1.

use std::path::Path;
use std::sync::Arc;

use toml::de::{DeTable, DeValue};

use crate::configuration::{
    AdmittedToml, ConfigKind, DeclaredCompositionValue, DeclaredDirectEvent,
    DeclaredF5GeneralConfig, DeclaredF5SamplingOverride, DeclaredF5Scenario,
    DeclaredF5Subpopulation, DeclaredF5Target, DeclaredGeneralConfig, DeclaredInitialPopulation,
    DeclaredNumber, DeclaredPopulationSnvs, DeclaredPopulationSnvsGeneration,
    DeclaredPopulationSnvsPopulation, DeclaredPopulationSnvsRegion, DeclaredReducedReference,
    DeclaredReductionTarget, DeclaredReference, DeclaredReferenceSource, DeclaredSequence,
    DeclaredTechnicalEvent, DeclaredTechnicalEventTarget, DeclaredTechnicalLibrary,
    DeclaredVcfOutput, ExpectedConfigKind, Located, SequenceClass, SourceLocation, semantic_error,
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
    // DTC-082 activa atómicamente la frontera de autoría F4: el runtime no
    // conserva lectores de TOML general v1 cuando ya puede producir PE real.
    let accepted_version = match expected {
        ExpectedConfigKind::General => "2",
        ExpectedConfigKind::GeneralDispatch => "2,3",
        ExpectedConfigKind::F5General => "3",
        ExpectedConfigKind::F5Scenario => "1",
    };
    let accepted = match expected {
        ExpectedConfigKind::GeneralDispatch => {
            matches!(version_value.get_ref().as_str(), Some("2" | "3"))
        }
        _ => version_value.get_ref().as_str() == Some(accepted_version),
    };
    if !accepted {
        return Err(Diagnostic::configuration(
            "GPS-CONF-005",
            "La versión de esquema declarada no está admitida.",
        )
        .with("path", path.display().to_string())
        .with("found", render_toml_value(version_value.get_ref()))
        .with("accepted", accepted_version)
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

/// La cabecera ya fue admitida por la frontera de CLI; por ello este lector
/// pequeño no vuelve a decidir errores de sintaxis ni de discriminador.
pub fn admitted_is_f5_general_v3(admitted: &AdmittedToml) -> bool {
    let source = std::str::from_utf8(admitted.bytes()).expect("admitted TOML is UTF-8");
    let parsed = DeTable::parse(source).expect("admitted TOML remains parseable");
    parsed
        .get_ref()
        .get("schema_version")
        .and_then(|value| value.get_ref().as_str())
        .is_some_and(|version| version == "3")
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
            "initial_populations",
            "direct_events",
            "vcf_output",
            "technical_event",
        ],
    )?;
    let seed = optional_located_string(&admitted, table, "seed")?;
    let ploidy = optional_located_u64(&admitted, table, "ploidy")?;
    let normalization_tolerance = optional_number(&admitted, table, "normalization_tolerance")?;
    let reference = parse_reference(&admitted, required_table(&admitted, table, "reference")?)?;
    let reduced = parse_reduced_references(&admitted, table)?;
    let population_snvs = parse_population_snvs(&admitted, table)?;
    let initial_populations = parse_initial_populations(&admitted, table)?;
    let direct_events = parse_direct_events(&admitted, table)?;
    let vcf_output = parse_vcf_output(&admitted, table)?;
    let technical_event = parse_technical_event(&admitted, table)?;
    Ok(DeclaredGeneralConfig::new(
        admitted,
        seed,
        ploidy,
        normalization_tolerance,
        reference,
        reduced,
        population_snvs,
    )
    .with_initial_populations(initial_populations)
    .with_direct_events(direct_events)
    .with_vcf_output(vcf_output)
    .with_technical_event(technical_event))
}

/// Lee el general v3 cerrado de F5. No comparte el lector v2: admitir v3 por
/// la ruta F1--F4 permitiría mezclar contratos de autoría incompatibles.
pub fn parse_declared_f5_general(
    admitted: AdmittedToml,
) -> Result<DeclaredF5GeneralConfig, Diagnostic> {
    if admitted.kind() != ConfigKind::General {
        return Err(semantic_error(
            &admitted,
            "config_kind",
            "F5 exige un documento general",
        ));
    }
    let source = std::str::from_utf8(admitted.bytes())
        .map_err(|error| semantic_error(&admitted, "document", error.to_string()))?;
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
            "reference_id",
            "reference_path",
            "reference_sha256",
            "f5_scenario",
        ],
    )?;
    let seed = optional_located_string(&admitted, table, "seed")?;
    if let Some(seed) = &seed {
        require_sha256(&admitted, "seed", seed)?;
    }
    let ploidy = required_u64(&admitted, table, "ploidy")?;
    if !(1..=10).contains(ploidy.value()) {
        return Err(semantic_error(
            &admitted,
            "ploidy",
            "debe pertenecer a 1..=10",
        ));
    }
    let reference_id = required_string(&admitted, table, "reference_id")?;
    require_id(&admitted, "reference_id", &reference_id)?;
    let reference_path = required_string(&admitted, table, "reference_path")?;
    require_relative_path(&admitted, "reference_path", &reference_path)?;
    let reference_sha256 = required_string(&admitted, table, "reference_sha256")?;
    require_sha256(&admitted, "reference_sha256", &reference_sha256)?;
    let link = required_table(&admitted, table, "f5_scenario")?;
    check_unknown(&admitted, link, &["scenario_id", "scenario_path"])?;
    let scenario_id = required_string(&admitted, link, "scenario_id")?;
    require_id(&admitted, "f5_scenario.scenario_id", &scenario_id)?;
    let scenario_path = required_string(&admitted, link, "scenario_path")?;
    require_relative_path(&admitted, "f5_scenario.scenario_path", &scenario_path)?;
    Ok(DeclaredF5GeneralConfig {
        source: admitted,
        seed,
        ploidy,
        reference_id,
        reference_path,
        reference_sha256,
        scenario_id,
        scenario_path,
    })
}

/// Lee el documento de política F5 v1, sin abrir perfiles ni máscaras. Esos
/// recursos pertenecen al lector de perfiles de F5.0.2.
pub fn parse_declared_f5_scenario(
    admitted: AdmittedToml,
) -> Result<DeclaredF5Scenario, Diagnostic> {
    if admitted.kind() != ConfigKind::F5Scenario {
        return Err(semantic_error(
            &admitted,
            "config_kind",
            "se esperaba f5_scenario",
        ));
    }
    let source = std::str::from_utf8(admitted.bytes())
        .map_err(|error| semantic_error(&admitted, "document", error.to_string()))?;
    let parsed = DeTable::parse(source)
        .map_err(|error| semantic_error(&admitted, "document", error.to_string()))?;
    let table = parsed.get_ref();
    check_unknown(
        &admitted,
        table,
        &[
            "config_kind",
            "schema_version",
            "scenario_id",
            "founder_count",
            "f5",
            "subpopulations",
            "sampling_overrides",
        ],
    )?;
    let scenario_id = required_string(&admitted, table, "scenario_id")?;
    require_id(&admitted, "scenario_id", &scenario_id)?;
    let founder_count = required_u64(&admitted, table, "founder_count")?;
    if *founder_count.value() == 0 {
        return Err(semantic_error(
            &admitted,
            "founder_count",
            "debe ser positivo",
        ));
    }
    let f5 = required_table(&admitted, table, "f5")?;
    check_unknown(
        &admitted,
        f5,
        &[
            "profile_metadata_path",
            "profile_metadata_sha256",
            "excluded_regions_path",
            "excluded_regions_sha256",
            "inbreeding_f",
            "max_layout_attempts",
            "max_haplotype_attempts",
            "max_swap_proposals",
            "max_haplotype_cells",
            "max_operation_units",
            "max_working_memory_bytes",
            "max_output_bytes",
            "block_count",
            "covered_bases",
            "coverage_fraction",
        ],
    )?;
    let profile_metadata_path = required_string(&admitted, f5, "profile_metadata_path")?;
    require_relative_path(
        &admitted,
        "f5.profile_metadata_path",
        &profile_metadata_path,
    )?;
    let profile_metadata_sha256 = required_string(&admitted, f5, "profile_metadata_sha256")?;
    require_sha256(
        &admitted,
        "f5.profile_metadata_sha256",
        &profile_metadata_sha256,
    )?;
    let excluded_regions_path = optional_located_string(&admitted, f5, "excluded_regions_path")?;
    let excluded_regions_sha256 =
        optional_located_string(&admitted, f5, "excluded_regions_sha256")?;
    if excluded_regions_path.is_some() != excluded_regions_sha256.is_some() {
        return Err(semantic_error(
            &admitted,
            "f5.excluded_regions_path,f5.excluded_regions_sha256",
            "máscara y hash deben aparecer juntos",
        ));
    }
    if let Some(path) = &excluded_regions_path {
        require_relative_path(&admitted, "f5.excluded_regions_path", path)?;
    }
    if let Some(hash) = &excluded_regions_sha256 {
        require_sha256(&admitted, "f5.excluded_regions_sha256", hash)?;
    }
    let inbreeding_f = optional_number(&admitted, f5, "inbreeding_f")?;
    if let Some(value) = &inbreeding_f {
        require_probability(&admitted, "f5.inbreeding_f", value)?;
    }
    let max_layout_attempts = optional_located_u64(&admitted, f5, "max_layout_attempts")?;
    let max_haplotype_attempts = optional_located_u64(&admitted, f5, "max_haplotype_attempts")?;
    let max_swap_proposals = optional_located_u64(&admitted, f5, "max_swap_proposals")?;
    let max_haplotype_cells = optional_located_u64(&admitted, f5, "max_haplotype_cells")?;
    let max_operation_units = optional_located_u64(&admitted, f5, "max_operation_units")?;
    let max_working_memory_bytes = optional_located_u64(&admitted, f5, "max_working_memory_bytes")?;
    let max_output_bytes = optional_located_u64(&admitted, f5, "max_output_bytes")?;
    require_limit(
        &admitted,
        "f5.max_layout_attempts",
        max_layout_attempts.as_ref(),
        100_000,
    )?;
    require_limit(
        &admitted,
        "f5.max_haplotype_attempts",
        max_haplotype_attempts.as_ref(),
        10_000,
    )?;
    require_limit(
        &admitted,
        "f5.max_swap_proposals",
        max_swap_proposals.as_ref(),
        10_000_000,
    )?;
    require_limit(
        &admitted,
        "f5.max_haplotype_cells",
        max_haplotype_cells.as_ref(),
        67_108_864,
    )?;
    require_limit(
        &admitted,
        "f5.max_operation_units",
        max_operation_units.as_ref(),
        268_435_456,
    )?;
    require_limit(
        &admitted,
        "f5.max_working_memory_bytes",
        max_working_memory_bytes.as_ref(),
        1_073_741_824,
    )?;
    require_limit(
        &admitted,
        "f5.max_output_bytes",
        max_output_bytes.as_ref(),
        2_147_483_648,
    )?;
    let block_count = optional_located_u64(&admitted, f5, "block_count")?;
    let covered_bases = optional_located_u64(&admitted, f5, "covered_bases")?;
    let coverage_fraction = optional_number(&admitted, f5, "coverage_fraction")?;
    if usize::from(block_count.is_some())
        + usize::from(covered_bases.is_some())
        + usize::from(coverage_fraction.is_some())
        != 1
    {
        return Err(semantic_error(
            &admitted,
            "f5.block_count,f5.covered_bases,f5.coverage_fraction",
            "se exige exactamente un objetivo",
        ));
    }
    let target = if let Some(value) = block_count {
        if *value.value() == 0 || *value.value() > u64::from(u32::MAX) {
            return Err(semantic_error(
                &admitted,
                "f5.block_count",
                "debe ser u32 positivo",
            ));
        }
        DeclaredF5Target::BlockCount(value)
    } else if let Some(value) = covered_bases {
        if *value.value() == 0 {
            return Err(semantic_error(
                &admitted,
                "f5.covered_bases",
                "debe ser positivo",
            ));
        }
        DeclaredF5Target::CoveredBases(value)
    } else {
        let value = coverage_fraction.expect("one-of checked");
        require_coverage_fraction(&admitted, &value)?;
        DeclaredF5Target::CoverageFraction(value)
    };
    let subpopulations = required_array(&admitted, table, "subpopulations")?
        .iter()
        .map(|entry| parse_f5_subpopulation(&admitted, entry))
        .collect::<Result<Vec<_>, _>>()?;
    let sampling_overrides = optional_array(&admitted, table, "sampling_overrides")?
        .unwrap_or_default()
        .iter()
        .map(|entry| parse_f5_override(&admitted, entry))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(DeclaredF5Scenario {
        source: admitted,
        scenario_id,
        founder_count,
        profile_metadata_path,
        profile_metadata_sha256,
        excluded_regions_path,
        excluded_regions_sha256,
        inbreeding_f,
        max_layout_attempts,
        max_haplotype_attempts,
        max_swap_proposals,
        max_haplotype_cells,
        max_operation_units,
        max_working_memory_bytes,
        max_output_bytes,
        target,
        subpopulations,
        sampling_overrides,
    })
}

/// Comprueba las invariantes que requieren ambos documentos. Las rutas ya se
/// han validado sintácticamente; su resolución física queda para el cargador.
pub fn validate_declared_f5_pair(
    general: &DeclaredF5GeneralConfig,
    scenario: &DeclaredF5Scenario,
) -> Result<(), Diagnostic> {
    if general.scenario_id.value() != scenario.scenario_id.value() {
        return Err(semantic_error(
            &general.source,
            "f5_scenario.scenario_id,scenario_id",
            "los IDs de escenario no coinciden",
        ));
    }
    let mut populations = std::collections::BTreeSet::new();
    let mut ids = std::collections::BTreeSet::new();
    let mut total = 0_u64;
    for population in &scenario.subpopulations {
        require_id(
            &scenario.source,
            "subpopulations.population_id",
            &population.population_id,
        )?;
        if !populations.insert(population.population_id.value()) {
            return Err(semantic_error(
                &scenario.source,
                "subpopulations.population_id",
                "población repetida",
            ));
        }
        if *population.individual_count.value() == 0 {
            return Err(semantic_error(
                &scenario.source,
                "subpopulations.individual_count",
                "debe ser positivo",
            ));
        }
        total = total
            .checked_add(*population.individual_count.value())
            .ok_or_else(|| {
                semantic_error(
                    &scenario.source,
                    "subpopulations.individual_count",
                    "la suma desborda u64",
                )
            })?;
        if population.founder_ids.is_some() == population.founder_keys.is_some() {
            return Err(semantic_error(
                &scenario.source,
                "subpopulations.founder_ids,subpopulations.founder_keys",
                "se exige exactamente uno",
            ));
        }
        let founders = population
            .founder_ids
            .as_ref()
            .or(population.founder_keys.as_ref())
            .expect("one-of checked");
        if founders.len()
            != usize::try_from(*population.individual_count.value()).unwrap_or(usize::MAX)
        {
            return Err(semantic_error(
                &scenario.source,
                "subpopulations.founder_ids,subpopulations.founder_keys",
                "la cardinalidad no coincide con individual_count",
            ));
        }
        let mut local = std::collections::BTreeSet::new();
        for founder in founders {
            require_id(&scenario.source, "subpopulations.founders", founder)?;
            if !local.insert(founder.value()) {
                return Err(semantic_error(
                    &scenario.source,
                    "subpopulations.founders",
                    "fundador repetido",
                ));
            }
        }
        if let Some(founder_ids) = &population.founder_ids {
            for founder in founder_ids {
                if !ids.insert(founder.value()) {
                    return Err(semantic_error(
                        &scenario.source,
                        "subpopulations.founder_ids",
                        "individual_id repetido globalmente",
                    ));
                }
            }
        }
        if let Some(value) = &population.inbreeding_f {
            require_probability(&scenario.source, "subpopulations.inbreeding_f", value)?;
        }
        let effective_f = population
            .inbreeding_f
            .as_ref()
            .or(scenario.inbreeding_f.as_ref())
            .map_or(0.0, |value| {
                value.literal().parse::<f64>().unwrap_or(f64::NAN)
            });
        if population
            .sampling
            .as_ref()
            .is_some_and(|mode| mode.value() == "controlled")
            && effective_f != 0.0
        {
            return Err(semantic_error(
                &scenario.source,
                "subpopulations.sampling,inbreeding_f",
                "controlled exige F efectivo igual a cero",
            ));
        }
    }
    if total != *scenario.founder_count.value() {
        return Err(semantic_error(
            &scenario.source,
            "founder_count,subpopulations.individual_count",
            "founder_count no coincide con la suma",
        ));
    }
    if *general.ploidy.value() != 2
        && (scenario.inbreeding_f.is_some()
            || scenario
                .subpopulations
                .iter()
                .any(|p| p.inbreeding_f.is_some()))
    {
        return Err(semantic_error(
            &scenario.source,
            "inbreeding_f,ploidy",
            "F solo se admite con ploidía 2",
        ));
    }
    let mut override_keys = std::collections::BTreeSet::new();
    for override_ in &scenario.sampling_overrides {
        if !populations.contains(override_.population_id.value()) {
            return Err(semantic_error(
                &scenario.source,
                "sampling_overrides.population_id",
                "población desconocida",
            ));
        }
        if !override_keys.insert((
            override_.population_id.value(),
            *override_.block_index.value(),
        )) {
            return Err(semantic_error(
                &scenario.source,
                "sampling_overrides",
                "sustitución duplicada",
            ));
        }
        if override_.mode.value() == "probabilistic" && override_.controlled_max_abs_error.is_some()
        {
            return Err(semantic_error(
                &scenario.source,
                "sampling_overrides.controlled_max_abs_error",
                "solo se admite en controlled",
            ));
        }
        if override_.mode.value() == "controlled" {
            if let Some(value) = &override_.controlled_max_abs_error {
                require_probability(
                    &scenario.source,
                    "sampling_overrides.controlled_max_abs_error",
                    value,
                )?;
            }
            let population = scenario
                .subpopulations
                .iter()
                .find(|population| {
                    population.population_id.value() == override_.population_id.value()
                })
                .expect("membership was checked");
            let effective_f = population
                .inbreeding_f
                .as_ref()
                .or(scenario.inbreeding_f.as_ref())
                .map_or(0.0, |value| {
                    value.literal().parse::<f64>().unwrap_or(f64::NAN)
                });
            if effective_f != 0.0 {
                return Err(semantic_error(
                    &scenario.source,
                    "sampling_overrides.mode,inbreeding_f",
                    "controlled exige F efectivo igual a cero",
                ));
            }
            if let Some(value) = &override_.controlled_max_abs_error {
                let copies = population
                    .individual_count
                    .value()
                    .checked_mul(*general.ploidy.value())
                    .ok_or_else(|| {
                        semantic_error(
                            &scenario.source,
                            "individual_count,ploidy",
                            "el número de copias desborda u64",
                        )
                    })?;
                let tolerance = value.literal().parse::<f64>().unwrap_or(f64::NAN);
                if tolerance > 0.5 / copies as f64 {
                    return Err(semantic_error(
                        &scenario.source,
                        "sampling_overrides.controlled_max_abs_error",
                        "excede 0.5/C",
                    ));
                }
            }
        }
    }
    Ok(())
}

fn parse_f5_subpopulation(
    admitted: &AdmittedToml,
    value: &toml::Spanned<DeValue<'_>>,
) -> Result<DeclaredF5Subpopulation, Diagnostic> {
    let table = value
        .get_ref()
        .as_table()
        .ok_or_else(|| wrong_type(admitted, "subpopulations", value, "table"))?;
    check_unknown(
        admitted,
        table,
        &[
            "population_id",
            "individual_count",
            "sampling",
            "inbreeding_f",
            "founder_ids",
            "founder_keys",
        ],
    )?;
    let sampling = optional_located_string(admitted, table, "sampling")?;
    if let Some(mode) = &sampling
        && !matches!(mode.value().as_str(), "probabilistic" | "controlled")
    {
        return Err(semantic_error(
            admitted,
            "subpopulations.sampling",
            "debe ser probabilistic o controlled",
        ));
    }
    Ok(DeclaredF5Subpopulation {
        population_id: required_string(admitted, table, "population_id")?,
        individual_count: required_u64(admitted, table, "individual_count")?,
        sampling,
        inbreeding_f: optional_number(admitted, table, "inbreeding_f")?,
        founder_ids: optional_f5_string_array(admitted, table, "founder_ids")?,
        founder_keys: optional_f5_string_array(admitted, table, "founder_keys")?,
    })
}

fn parse_f5_override(
    admitted: &AdmittedToml,
    value: &toml::Spanned<DeValue<'_>>,
) -> Result<DeclaredF5SamplingOverride, Diagnostic> {
    let table = value
        .get_ref()
        .as_table()
        .ok_or_else(|| wrong_type(admitted, "sampling_overrides", value, "table"))?;
    check_unknown(
        admitted,
        table,
        &[
            "population_id",
            "block_index",
            "mode",
            "controlled_max_abs_error",
        ],
    )?;
    let mode = required_string(admitted, table, "mode")?;
    if !matches!(mode.value().as_str(), "probabilistic" | "controlled") {
        return Err(semantic_error(
            admitted,
            "sampling_overrides.mode",
            "debe ser probabilistic o controlled",
        ));
    }
    let block_index = required_u64(admitted, table, "block_index")?;
    if *block_index.value() > u64::from(u32::MAX) {
        return Err(semantic_error(
            admitted,
            "sampling_overrides.block_index",
            "debe pertenecer a u32",
        ));
    }
    Ok(DeclaredF5SamplingOverride {
        population_id: required_string(admitted, table, "population_id")?,
        block_index,
        mode,
        controlled_max_abs_error: optional_number(admitted, table, "controlled_max_abs_error")?,
    })
}

fn optional_f5_string_array(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
    field: &str,
) -> Result<Option<Vec<Located<String>>>, Diagnostic> {
    optional_array(admitted, table, field)?
        .map(|values| {
            values
                .iter()
                .map(|value| {
                    value
                        .get_ref()
                        .as_str()
                        .map(|text| Located::new(text.to_owned(), location(admitted, value.span())))
                        .ok_or_else(|| wrong_type(admitted, field, value, "string"))
                })
                .collect()
        })
        .transpose()
}

fn require_id(
    admitted: &AdmittedToml,
    field: &str,
    value: &Located<String>,
) -> Result<(), Diagnostic> {
    let text = value.value();
    let valid = !text.is_empty()
        && text.len() <= 64
        && text.as_bytes()[0].is_ascii_alphanumeric()
        && text
            .bytes()
            .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'));
    valid
        .then_some(())
        .ok_or_else(|| semantic_error(admitted, field, "ID inválido conforme a DTC-039"))
}

fn require_sha256(
    admitted: &AdmittedToml,
    field: &str,
    value: &Located<String>,
) -> Result<(), Diagnostic> {
    let valid = value.value().len() == 64
        && value
            .value()
            .bytes()
            .all(|byte| byte.is_ascii_digit() || matches!(byte, b'a'..=b'f'));
    valid.then_some(()).ok_or_else(|| {
        semantic_error(
            admitted,
            field,
            "debe ser SHA-256 hexadecimal en minúsculas",
        )
    })
}

fn require_relative_path(
    admitted: &AdmittedToml,
    field: &str,
    value: &Located<String>,
) -> Result<(), Diagnostic> {
    // `Path::components` no sirve como frontera léxica: normaliza `//` y los
    // segmentos `.` antes de que podamos distinguir una grafía no canónica.
    // Esta comprobación replica el alfabeto y la segmentación de los esquemas
    // F5 de autoría antes de construir cualquier `Path`.
    let valid = !value.value().is_empty()
        && !value.value().starts_with('/')
        && !value.value().ends_with('/')
        && value.value().split('/').all(|segment| {
            !segment.is_empty()
                && !matches!(segment, "." | "..")
                && segment
                    .bytes()
                    .all(|byte| byte.is_ascii_alphanumeric() || matches!(byte, b'.' | b'_' | b'-'))
        });
    valid.then_some(()).ok_or_else(|| {
        semantic_error(
            admitted,
            field,
            "debe ser una ruta relativa léxicamente normalizada",
        )
    })
}

fn require_probability(
    admitted: &AdmittedToml,
    field: &str,
    value: &DeclaredNumber,
) -> Result<(), Diagnostic> {
    value
        .literal()
        .parse::<f64>()
        .ok()
        .filter(|number| number.is_finite() && (0.0..=1.0).contains(number))
        .map(|_| ())
        .ok_or_else(|| semantic_error(admitted, field, "debe ser finito y pertenecer a [0,1]"))
}

fn require_coverage_fraction(
    admitted: &AdmittedToml,
    value: &DeclaredNumber,
) -> Result<(), Diagnostic> {
    let literal = value.literal();
    let valid_decimal = literal.split_once('.').is_some_and(|(whole, fraction)| {
        (1..=18).contains(&fraction.len())
            && fraction.bytes().all(|byte| byte.is_ascii_digit())
            && match whole {
                "0" => fraction.bytes().any(|byte| byte != b'0'),
                "1" => fraction.bytes().all(|byte| byte == b'0'),
                _ => false,
            }
    });
    valid_decimal.then_some(()).ok_or_else(|| {
        semantic_error(
            admitted,
            "f5.coverage_fraction",
            "debe ser decimal literal en (0,1] con 1..=18 decimales",
        )
    })
}

fn require_limit(
    admitted: &AdmittedToml,
    field: &str,
    value: Option<&Located<u64>>,
    maximum: u64,
) -> Result<(), Diagnostic> {
    value.map_or(Ok(()), |value| {
        (*value.value() > 0 && *value.value() <= maximum)
            .then_some(())
            .ok_or_else(|| semantic_error(admitted, field, format!("debe estar en 1..={maximum}")))
    })
}

fn parse_technical_event(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Option<DeclaredTechnicalEvent>, Diagnostic> {
    let Some(event) = optional_table(admitted, table, "technical_event")? else {
        return Ok(None);
    };
    check_unknown(
        admitted,
        event,
        &[
            "event_id",
            "individual_ids",
            "target_coverage_wgs",
            "sequencing_errors",
            "read_length_bases",
            "chimera_rate",
            "index_hopping_rate",
            "contamination_min_percent",
            "contamination_max_percent",
            "contamination_recipient_ids",
            "libraries",
        ],
    )?;
    let target = optional_number(admitted, event, "target_coverage_wgs")?
        .map(DeclaredTechnicalEventTarget::CoverageWgs)
        .ok_or_else(|| {
            semantic_error(
                admitted,
                "technical_event.target_coverage_wgs",
                "F4 exige una cobertura objetivo; read_count y pair_count no están admitidos",
            )
        })?;
    let individuals = required_array(admitted, event, "individual_ids")?;
    if individuals.is_empty() {
        return Err(semantic_error(
            admitted,
            "technical_event.individual_ids",
            "debe contener al menos un individuo",
        ));
    }
    let individual_ids = individuals
        .iter()
        .map(|value| {
            value
                .get_ref()
                .as_str()
                .map(|text| Located::new(text.to_owned(), location(admitted, value.span())))
                .ok_or_else(|| {
                    wrong_type(admitted, "technical_event.individual_ids", value, "string")
                })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let libraries = required_array(admitted, event, "libraries")?
        .iter()
        .map(|value| {
            let library = value
                .get_ref()
                .as_table()
                .ok_or_else(|| wrong_type(admitted, "technical_event.libraries", value, "table"))?;
            check_unknown(admitted, library, &["individual_id", "i7", "i5"])?;
            Ok(DeclaredTechnicalLibrary {
                individual_id: required_string(admitted, library, "individual_id")?,
                i7: required_string(admitted, library, "i7")?,
                i5: required_string(admitted, library, "i5")?,
            })
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    if libraries.is_empty() {
        return Err(semantic_error(
            admitted,
            "technical_event.libraries",
            "debe declarar exactamente una biblioteca por individuo",
        ));
    }
    let contamination_recipient_ids =
        optional_array(admitted, event, "contamination_recipient_ids")?
            .map(|values| {
                values
                    .iter()
                    .map(|value| {
                        value
                            .get_ref()
                            .as_str()
                            .map(|text| {
                                Located::new(text.to_owned(), location(admitted, value.span()))
                            })
                            .ok_or_else(|| {
                                wrong_type(
                                    admitted,
                                    "technical_event.contamination_recipient_ids",
                                    value,
                                    "string",
                                )
                            })
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()
            })
            .transpose()?;
    Ok(Some(DeclaredTechnicalEvent {
        event_id: required_string(admitted, event, "event_id")?,
        individual_ids,
        target,
        sequencing_errors: optional_located_string(admitted, event, "sequencing_errors")?,
        read_length_bases: optional_located_u64(admitted, event, "read_length_bases")?,
        chimera_rate: optional_number(admitted, event, "chimera_rate")?,
        index_hopping_rate: optional_number(admitted, event, "index_hopping_rate")?,
        contamination_min_percent: optional_number(admitted, event, "contamination_min_percent")?,
        contamination_max_percent: optional_number(admitted, event, "contamination_max_percent")?,
        contamination_recipient_ids,
        libraries,
    }))
}

fn parse_vcf_output(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Option<DeclaredVcfOutput>, Diagnostic> {
    let Some(output) = optional_table(admitted, table, "vcf_output")? else {
        return Ok(None);
    };
    check_unknown(admitted, output, &["mode"])?;
    let mode = required_string(admitted, output, "mode")?;
    if !matches!(mode.value().as_str(), "text" | "bgzf") {
        return Err(
            semantic_error(admitted, "vcf_output.mode", "debe ser text o bgzf")
                .with("location", mode.location().to_string()),
        );
    }
    Ok(Some(DeclaredVcfOutput { mode }))
}

fn parse_direct_events(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Vec<DeclaredDirectEvent>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "direct_events")? else {
        return Ok(Vec::new());
    };
    values
        .iter()
        .map(|value| {
            let entry = value
                .get_ref()
                .as_table()
                .ok_or_else(|| wrong_type(admitted, "direct_events", value, "table"))?;
            check_unknown(
                admitted,
                entry,
                &[
                    "event_id",
                    "individual_id",
                    "phase",
                    "kind",
                    "contig",
                    "position_1based",
                    "ref",
                    "alt",
                    "point",
                    "bases",
                    "start",
                    "end",
                ],
            )?;
            let kind = required_string(admitted, entry, "kind")?;
            if !matches!(kind.value().as_str(), "snp" | "insertion" | "deletion") {
                return Err(semantic_error(
                    admitted,
                    "direct_events.kind",
                    "debe ser snp, insertion o deletion",
                ));
            }
            let position_1based = optional_located_u64(admitted, entry, "position_1based")?;
            let reference = optional_located_string(admitted, entry, "ref")?;
            let alternate = optional_located_string(admitted, entry, "alt")?;
            let point = optional_located_u64(admitted, entry, "point")?;
            let bases = optional_located_string(admitted, entry, "bases")?;
            let start = optional_located_u64(admitted, entry, "start")?;
            let end = optional_located_u64(admitted, entry, "end")?;
            let valid = match kind.value().as_str() {
                "snp" => {
                    position_1based.is_some()
                        && reference.is_some()
                        && alternate.is_some()
                        && point.is_none()
                        && bases.is_none()
                        && start.is_none()
                        && end.is_none()
                }
                "insertion" => {
                    point.is_some()
                        && bases.is_some()
                        && position_1based.is_none()
                        && reference.is_none()
                        && alternate.is_none()
                        && start.is_none()
                        && end.is_none()
                }
                "deletion" => {
                    start.is_some()
                        && end.is_some()
                        && position_1based.is_none()
                        && reference.is_none()
                        && alternate.is_none()
                        && point.is_none()
                        && bases.is_none()
                }
                _ => false,
            };
            if !valid {
                return Err(semantic_error(
                    admitted,
                    "direct_events",
                    "los campos no corresponden a kind",
                ));
            }
            Ok(DeclaredDirectEvent {
                event_id: required_string(admitted, entry, "event_id")?,
                individual_id: required_string(admitted, entry, "individual_id")?,
                phase: required_u64(admitted, entry, "phase")?,
                kind,
                contig: required_string(admitted, entry, "contig")?,
                position_1based,
                reference,
                alternate,
                point,
                bases,
                start,
                end,
            })
        })
        .collect()
}

fn parse_initial_populations(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<Vec<DeclaredInitialPopulation>, Diagnostic> {
    let Some(values) = optional_array(admitted, table, "initial_populations")? else {
        return Ok(Vec::new());
    };
    let mut populations = Vec::with_capacity(values.len());
    for value in values {
        let entry = value
            .get_ref()
            .as_table()
            .ok_or_else(|| wrong_type(admitted, "initial_populations", value, "table"))?;
        check_unknown(
            admitted,
            entry,
            &[
                "population_id",
                "individual_count",
                "ploidy",
                "genotype_model",
                "inbreeding_coefficient",
                "individual_ids",
                "founder_keys",
            ],
        )?;
        let genotype_model = required_string(admitted, entry, "genotype_model")?;
        if !matches!(genotype_model.value().as_str(), "hwe" | "f_diploid") {
            return Err(semantic_error(
                admitted,
                "initial_populations.genotype_model",
                "debe ser hwe o f_diploid",
            ));
        }
        let inbreeding_coefficient = optional_number(admitted, entry, "inbreeding_coefficient")?;
        match (
            genotype_model.value().as_str(),
            inbreeding_coefficient.is_some(),
        ) {
            ("hwe", true) => {
                return Err(semantic_error(
                    admitted,
                    "initial_populations.inbreeding_coefficient",
                    "hwe no admite inbreeding_coefficient",
                ));
            }
            ("f_diploid", false) => {
                return Err(semantic_error(
                    admitted,
                    "initial_populations.inbreeding_coefficient",
                    "f_diploid exige inbreeding_coefficient",
                ));
            }
            _ => {}
        }
        let individual_ids = optional_array(admitted, entry, "individual_ids")?
            .map(|values| {
                values
                    .iter()
                    .map(|value| {
                        value
                            .get_ref()
                            .as_str()
                            .map(|text| {
                                Located::new(text.to_owned(), location(admitted, value.span()))
                            })
                            .ok_or_else(|| {
                                wrong_type(
                                    admitted,
                                    "initial_populations.individual_ids",
                                    value,
                                    "string",
                                )
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        let founder_keys = optional_array(admitted, entry, "founder_keys")?
            .map(|values| {
                values
                    .iter()
                    .map(|value| {
                        value
                            .get_ref()
                            .as_str()
                            .map(|text| {
                                Located::new(text.to_owned(), location(admitted, value.span()))
                            })
                            .ok_or_else(|| {
                                wrong_type(
                                    admitted,
                                    "initial_populations.founder_keys",
                                    value,
                                    "string",
                                )
                            })
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .transpose()?;
        populations.push(DeclaredInitialPopulation {
            population_id: required_string(admitted, entry, "population_id")?,
            individual_count: required_u64(admitted, entry, "individual_count")?,
            ploidy: required_u64(admitted, entry, "ploidy")?,
            genotype_model,
            inbreeding_coefficient,
            individual_ids,
            founder_keys,
        });
    }
    Ok(populations)
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
    check_unknown(
        admitted,
        entry,
        &["mode", "input_path", "metadata_path", "generate"],
    )?;
    let mode = required_string(admitted, entry, "mode")?;
    match mode.value().as_str() {
        "supplied" => {
            if entry.contains_key("metadata_path") || entry.contains_key("generate") {
                return Err(semantic_error(
                    admitted,
                    "population_snvs.mode,population_snvs.metadata_path",
                    "el modo supplied no admite metadata_path ni generate",
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
            if entry.contains_key("input_path") || entry.contains_key("generate") {
                return Err(semantic_error(
                    admitted,
                    "population_snvs.mode,population_snvs.input_path",
                    "el modo reuse no admite input_path ni generate",
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
            Ok(Some(DeclaredPopulationSnvs::Generate {
                parameters: Box::new(parse_population_snvs_generation(
                    admitted,
                    required_table(admitted, entry, "generate")?,
                )?),
            }))
        }
        other => Err(semantic_error(
            admitted,
            "population_snvs.mode",
            format!("modo de matriz SNV desconocido: {other}"),
        )),
    }
}

fn parse_population_snvs_generation(
    admitted: &AdmittedToml,
    table: &DeTable<'_>,
) -> Result<DeclaredPopulationSnvsGeneration, Diagnostic> {
    check_unknown(
        admitted,
        table,
        &[
            "count_mode",
            "locus_count",
            "mean_loci",
            "mean_density_per_base",
            "spatial_mode",
            "substitution_profile",
            "outside_weight",
            "regions",
            "populations",
        ],
    )?;
    let count_mode = required_string(admitted, table, "count_mode")?;
    let locus_count = optional_located_u64(admitted, table, "locus_count")?;
    let mean_loci = optional_number(admitted, table, "mean_loci")?;
    let mean_density_per_base = optional_number(admitted, table, "mean_density_per_base")?;
    let count_fields = usize::from(locus_count.is_some())
        + usize::from(mean_loci.is_some())
        + usize::from(mean_density_per_base.is_some());
    if count_fields != 1
        || !matches!(
            (
                count_mode.value().as_str(),
                locus_count.is_some(),
                mean_loci.is_some(),
                mean_density_per_base.is_some()
            ),
            ("exact", true, false, false)
                | ("poisson", false, true, false)
                | ("poisson", false, false, true)
        )
    {
        return Err(semantic_error(
            admitted,
            "population_snvs.generate",
            "count_mode exige exactamente el parámetro de recuento compatible",
        ));
    }
    let spatial_mode = required_string(admitted, table, "spatial_mode")?;
    let substitution_profile = required_string(admitted, table, "substitution_profile")?;
    if !matches!(
        substitution_profile.value().as_str(),
        "uniform_snv12_v1" | "human_harris2017_snv12_approx_v1"
    ) {
        return Err(semantic_error(
            admitted,
            "population_snvs.generate.substitution_profile",
            "substitution_profile debe ser uniform_snv12_v1 o human_harris2017_snv12_approx_v1",
        ));
    }
    let outside_weight = optional_number(admitted, table, "outside_weight")?;
    let regions = optional_array(admitted, table, "regions")?
        .unwrap_or_default()
        .iter()
        .map(|value| {
            let region = value.get_ref().as_table().ok_or_else(|| {
                wrong_type(admitted, "population_snvs.generate.regions", value, "table")
            })?;
            check_unknown(admitted, region, &["contig", "start", "end", "weight"])?;
            Ok(DeclaredPopulationSnvsRegion {
                contig: required_string(admitted, region, "contig")?,
                start: required_u64(admitted, region, "start")?,
                end: required_u64(admitted, region, "end")?,
                weight: required_number(admitted, region, "weight")?,
            })
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    if (spatial_mode.value() == "weighted_regions") != outside_weight.is_some()
        || !matches!(
            spatial_mode.value().as_str(),
            "uniform" | "weighted_regions"
        )
    {
        return Err(semantic_error(
            admitted,
            "population_snvs.generate.spatial_mode",
            "weighted_regions exige outside_weight; uniform no lo admite",
        ));
    }
    if spatial_mode.value() == "uniform" && !regions.is_empty() {
        return Err(semantic_error(
            admitted,
            "population_snvs.generate.regions",
            "uniform no admite regiones",
        ));
    }
    let populations = required_array(admitted, table, "populations")?
        .iter()
        .map(|value| {
            let population = value.get_ref().as_table().ok_or_else(|| {
                wrong_type(
                    admitted,
                    "population_snvs.generate.populations",
                    value,
                    "table",
                )
            })?;
            check_unknown(
                admitted,
                population,
                &[
                    "population_id",
                    "frequency_model",
                    "frequency_mode",
                    "min_alt_af",
                    "max_alt_af",
                    "spectrum",
                ],
            )?;
            let frequency_model = required_string(admitted, population, "frequency_model")?;
            let frequency_mode = optional_located_string(admitted, population, "frequency_mode")?;
            let min_alt_af = optional_number(admitted, population, "min_alt_af")?;
            let max_alt_af = optional_number(admitted, population, "max_alt_af")?;
            let spectrum = optional_array(admitted, population, "spectrum")?
                .unwrap_or_default()
                .iter()
                .map(|value| {
                    let component = value.get_ref().as_table().ok_or_else(|| {
                        wrong_type(admitted, "population_snvs.generate.populations.spectrum", value, "table")
                    })?;
                    check_unknown(
                        admitted,
                        component,
                        &["kind", "fixed_alt_af", "min_alt_af", "max_alt_af", "weight", "sampling"],
                    )?;
                    let kind = required_string(admitted, component, "kind")?;
                    let fixed_alt_af = optional_number(admitted, component, "fixed_alt_af")?;
                    let min_alt_af = optional_number(admitted, component, "min_alt_af")?;
                    let max_alt_af = optional_number(admitted, component, "max_alt_af")?;
                    let sampling = optional_located_string(admitted, component, "sampling")?;
                    let valid = match kind.value().as_str() {
                        "fixed" => fixed_alt_af.is_some()
                            && min_alt_af.is_none()
                            && max_alt_af.is_none()
                            && sampling.is_none(),
                        "interval" => fixed_alt_af.is_none()
                            && min_alt_af.is_some()
                            && max_alt_af.is_some()
                            && sampling.is_some(),
                        _ => false,
                    };
                    if !valid {
                        return Err(semantic_error(
                            admitted,
                            "population_snvs.generate.populations.spectrum",
                            "fixed exige solo fixed_alt_af; interval exige min_alt_af, max_alt_af y sampling",
                        ));
                    }
                    Ok(crate::configuration::DeclaredFrequencySpectrumComponent {
                        kind,
                        fixed_alt_af,
                        min_alt_af,
                        max_alt_af,
                        weight: required_number(admitted, component, "weight")?,
                        sampling,
                    })
                })
                .collect::<Result<Vec<_>, Diagnostic>>()?;
            let simple_uniform = frequency_model.value() == "uniform"
                && frequency_mode.is_none()
                && min_alt_af.is_some()
                && max_alt_af.is_some()
                && spectrum.is_empty();
            let supplied_spectrum = frequency_model.value() == "spectrum"
                && min_alt_af.is_none()
                && max_alt_af.is_none()
                && !spectrum.is_empty();
            let builtin = matches!(
                frequency_model.value().as_str(),
                "human_1000g_global_sfs_approx_v1"
            ) && min_alt_af.is_none()
                && max_alt_af.is_none()
                && spectrum.is_empty();
            if !(simple_uniform || supplied_spectrum || builtin) {
                return Err(semantic_error(
                    admitted,
                    "population_snvs.generate.populations",
                    "uniform exige límites y omite frequency_mode; spectrum exige componentes y human_1000g_global_sfs_approx_v1 no admite límites",
                ));
            }
            Ok(DeclaredPopulationSnvsPopulation {
                population_id: required_string(admitted, population, "population_id")?,
                frequency_model,
                frequency_mode,
                min_alt_af,
                max_alt_af,
                spectrum,
            })
        })
        .collect::<Result<Vec<_>, Diagnostic>>()?;
    if populations.is_empty() {
        return Err(semantic_error(
            admitted,
            "population_snvs.generate.populations",
            "se exige al menos una población",
        ));
    }
    for population in &populations {
        if let Some(mode) = &population.frequency_mode
            && !matches!(mode.value().as_str(), "stochastic" | "exact")
        {
            return Err(semantic_error(
                admitted,
                "population_snvs.generate.populations",
                "frequency_mode debe ser stochastic o exact",
            ));
        }
    }
    Ok(DeclaredPopulationSnvsGeneration {
        count_mode,
        locus_count,
        mean_loci,
        mean_density_per_base,
        spatial_mode,
        substitution_profile,
        outside_weight,
        regions,
        populations,
    })
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
            check_unknown(admitted, composition_table, &["A", "C", "G", "T", "N"])?;
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
            DeclaredReferenceSource::Simulated {
                composition,
                unknown_base_fraction: optional_number(admitted, composition_table, "N")?,
            }
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
        let field = key.get_ref().as_ref();
        if let Some(deferred_to) = deferred_f4_alignment_destination(field) {
            return Err(Diagnostic::configuration(
                "GPS-CONF-023",
                "La configuración solicita una categoría PE observable o un alineamiento directo aplazado.",
            )
            .with("path", admitted.path().display().to_string())
            .with("field", field)
            .with("location", location(admitted, key.span()).to_string())
            .with("deferred_to", deferred_to));
        }
        return Err(Diagnostic::configuration(
            "GPS-CONF-006",
            "La configuración contiene un campo desconocido.",
        )
        .with("path", admitted.path().display().to_string())
        .with("field", field)
        .with("location", location(admitted, key.span()).to_string()));
    }
    Ok(())
}

/// Clasifica nombres reservados que F4 rechaza antes de crear cualquier
/// configuración efectiva, estadística o salida de alineamiento. No introduce
/// una representación de estas categorías: solo reserva sus solicitudes para
/// los destinos de diseño aprobados.
fn deferred_f4_alignment_destination(field: &str) -> Option<&'static str> {
    let family = match field {
        "paired_end_observable"
        | "paired_end_observables"
        | "orientation_anomaly"
        | "orientation_anomalies"
        | "orientation_anomaly_rate"
        | "orientation_anomaly_rates"
        | "observed_insert_size"
        | "observed_insert_sizes"
        | "observed_insert_size_distribution"
        | "observed_insert_size_distributions"
        | "observed_insert_length"
        | "observed_insert_lengths"
        | "observed_insert_length_distribution"
        | "observed_insert_length_distributions"
        | "unaligned_mate"
        | "unaligned_mates"
        | "unaligned_mate_rate"
        | "unaligned_mate_rates" => DeferredF4AlignmentFamily::ObservablePairedEnd,
        "direct_alignment"
        | "direct_alignments"
        | "direct_alignment_mode"
        | "direct_alignment_modes"
        | "direct_alignment_category"
        | "direct_alignment_categories"
        | "synthetic_alignment"
        | "synthetic_alignments"
        | "synthetic_alignment_category"
        | "synthetic_alignment_categories"
        | "alignment_category"
        | "alignment_categories" => DeferredF4AlignmentFamily::DirectAlignment,
        _ => return None,
    };
    Some(family.destination())
}

#[derive(Clone, Copy)]
enum DeferredF4AlignmentFamily {
    ObservablePairedEnd,
    DirectAlignment,
}

impl DeferredF4AlignmentFamily {
    const fn destination(self) -> &'static str {
        match self {
            Self::ObservablePairedEnd => "DTF-055",
            Self::DirectAlignment => "DTF-051",
        }
    }
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

fn optional_table<'a>(
    admitted: &AdmittedToml,
    table: &'a DeTable<'a>,
    field: &str,
) -> Result<Option<&'a DeTable<'a>>, Diagnostic> {
    let Some(value) = table.get(field) else {
        return Ok(None);
    };
    value
        .get_ref()
        .as_table()
        .map(Some)
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
