//! Lector estricto del paquete de perfiles F5 v1.
//!
//! Esta frontera solo admite los cinco recursos declarados por el metadata. No
//! decide geometría, marginales ni haplotipos: entrega tablas ya verificadas a
//! los pasos posteriores.

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

use serde_json::{Map, Value};
use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;

/// Tolerancia general DTC-049 para pesos de un perfil F5.
pub const F5_PROFILE_NORMALIZATION_TOLERANCE: f64 = 1e-9;
pub const F5_PROFILE_NORMALIZATION_ALGORITHM: &str = "neumaier-f64-v1";

#[derive(Clone, Debug, PartialEq)]
pub struct F5Profile {
    profile_id: String,
    block_lengths: Vec<WeightedInterval<u64>>,
    snp_gaps: Vec<WeightedInterval<u64>>,
    snp_spectrum: Vec<WeightedInterval<f64>>,
    haplotype_weights: Vec<f64>,
    ld_bins: BTreeMap<String, Vec<LdBin>>,
    normalizations: Vec<F5NormalizationAudit>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum F5NormalizationSet {
    BlockLengths,
    SnpGaps,
    SnpSpectrum,
    HaplotypeWeights,
}

impl F5NormalizationSet {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::BlockLengths => "block_lengths",
            Self::SnpGaps => "snp_gaps",
            Self::SnpSpectrum => "snp_spectrum",
            Self::HaplotypeWeights => "haplotype_weights",
        }
    }

    #[must_use]
    pub const fn canonical_order(self) -> &'static str {
        match self {
            Self::BlockLengths => "length_min_bp",
            Self::SnpGaps => "gap_min_bp",
            Self::SnpSpectrum => "alt_af_min",
            Self::HaplotypeWeights => "haplotype_rank",
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct F5NamedProbability {
    key: String,
    value: f64,
}

impl F5NamedProbability {
    #[must_use]
    pub fn key(&self) -> &str {
        &self.key
    }

    #[must_use]
    pub const fn value(&self) -> f64 {
        self.value
    }
}

/// Auditoría DTC-049 tipada y completa de uno de los cuatro conjuntos del perfil.
#[derive(Clone, Debug, PartialEq)]
pub struct F5NormalizationAudit {
    set: F5NormalizationSet,
    input_values: Vec<F5NamedProbability>,
    input_sum: f64,
    input_delta: f64,
    effective_values: Vec<F5NamedProbability>,
    output_sum: f64,
    output_delta: f64,
}

impl F5NormalizationAudit {
    #[must_use]
    pub const fn set(&self) -> F5NormalizationSet {
        self.set
    }

    #[must_use]
    pub const fn algorithm(&self) -> &'static str {
        F5_PROFILE_NORMALIZATION_ALGORITHM
    }

    #[must_use]
    pub const fn tolerance(&self) -> f64 {
        F5_PROFILE_NORMALIZATION_TOLERANCE
    }

    #[must_use]
    pub fn input_values(&self) -> &[F5NamedProbability] {
        &self.input_values
    }

    #[must_use]
    pub const fn input_sum(&self) -> f64 {
        self.input_sum
    }

    #[must_use]
    pub const fn input_delta(&self) -> f64 {
        self.input_delta
    }

    #[must_use]
    pub fn effective_values(&self) -> &[F5NamedProbability] {
        &self.effective_values
    }

    #[must_use]
    pub const fn output_sum(&self) -> f64 {
        self.output_sum
    }

    #[must_use]
    pub const fn output_delta(&self) -> f64 {
        self.output_delta
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct WeightedInterval<T> {
    minimum: T,
    maximum: T,
    weight: f64,
}

impl<T: Copy> WeightedInterval<T> {
    #[must_use]
    pub const fn minimum(&self) -> T {
        self.minimum
    }
    #[must_use]
    pub const fn maximum(&self) -> T {
        self.maximum
    }
    #[must_use]
    pub const fn weight(&self) -> f64 {
        self.weight
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct LdBin {
    distance_start_bp: u64,
    distance_end_bp: u64,
    target_mean_r2: f64,
}

impl LdBin {
    #[must_use]
    pub const fn distance_start_bp(&self) -> u64 {
        self.distance_start_bp
    }
    #[must_use]
    pub const fn distance_end_bp(&self) -> u64 {
        self.distance_end_bp
    }
    #[must_use]
    pub const fn target_mean_r2(&self) -> f64 {
        self.target_mean_r2
    }
}

impl F5Profile {
    #[must_use]
    pub fn profile_id(&self) -> &str {
        &self.profile_id
    }
    #[must_use]
    pub fn block_lengths(&self) -> &[WeightedInterval<u64>] {
        &self.block_lengths
    }
    #[must_use]
    pub fn snp_gaps(&self) -> &[WeightedInterval<u64>] {
        &self.snp_gaps
    }
    #[must_use]
    pub fn snp_spectrum(&self) -> &[WeightedInterval<f64>] {
        &self.snp_spectrum
    }
    #[must_use]
    pub fn haplotype_weights(&self) -> &[f64] {
        &self.haplotype_weights
    }
    #[must_use]
    pub fn ld_bins(&self, population_id: &str) -> Option<&[LdBin]> {
        self.ld_bins.get(population_id).map(Vec::as_slice)
    }
    #[must_use]
    pub fn normalizations(&self) -> &[F5NormalizationAudit] {
        &self.normalizations
    }
}

/// Lee una vez el metadata y exactamente los cinco TSV enlazados. `population_ids`
/// procede del escenario ya admitido y evita cualquier remapeo implícito.
pub fn load_f5_profile(
    metadata_path: &Path,
    expected_metadata_sha256: &str,
    population_ids: &BTreeSet<String>,
) -> Result<F5Profile, Diagnostic> {
    let bytes = read(metadata_path)?;
    verify_identity(metadata_path, &bytes, None, expected_metadata_sha256)?;
    let metadata: Value = serde_json::from_slice(&bytes)
        .map_err(|error| format_error(metadata_path, "profile.metadata.json", error.to_string()))?;
    let root = object(&metadata, metadata_path, "metadata")?;
    exact_keys(
        root,
        &[
            "profile_kind",
            "schema_version",
            "profile_id",
            "provenance_kind",
            "source_citation",
            "source_license",
            "source_sha256",
            "files",
        ],
        metadata_path,
        "metadata",
    )?;
    required_string(root, "profile_kind", metadata_path)?
        .eq("f5_generation_profile")
        .then_some(())
        .ok_or_else(|| format_error(metadata_path, "profile_kind", "literal inválido"))?;
    required_string(root, "schema_version", metadata_path)?
        .eq("1")
        .then_some(())
        .ok_or_else(|| format_error(metadata_path, "schema_version", "versión inválida"))?;
    let profile_id = required_id(root, "profile_id", metadata_path)?.to_owned();
    let provenance = required_string(root, "provenance_kind", metadata_path)?;
    if !matches!(provenance, "built_in" | "user_supplied") {
        return Err(format_error(
            metadata_path,
            "provenance_kind",
            "valor inválido",
        ));
    }
    if required_string(root, "source_citation", metadata_path)?.is_empty() {
        return Err(format_error(
            metadata_path,
            "source_citation",
            "debe ser no vacío",
        ));
    }
    nullable_string(root, "source_license", metadata_path)?;
    nullable_sha(root, "source_sha256", metadata_path)?;
    let files = root
        .get("files")
        .and_then(Value::as_array)
        .ok_or_else(|| format_error(metadata_path, "files", "debe ser una lista"))?;
    if files.len() != 5 {
        return Err(format_error(
            metadata_path,
            "files",
            "debe contener exactamente cinco enlaces",
        ));
    }
    let expected = [
        ("block_lengths.tsv", "f5_block_length_profile"),
        ("snp_gaps.tsv", "f5_snp_gap_profile"),
        ("snp_spectrum.tsv", "f5_snp_spectrum"),
        ("haplotype_weights.tsv", "f5_haplotype_weights"),
        ("ld_profile.tsv", "f5_ld_profile"),
    ];
    let mut contents = Vec::with_capacity(5);
    for (entry, (name, kind)) in files.iter().zip(expected) {
        let entry = object(entry, metadata_path, "files[]")?;
        exact_keys(
            entry,
            &[
                "artifact_kind",
                "schema_version",
                "path",
                "size_bytes",
                "sha256",
            ],
            metadata_path,
            "files[]",
        )?;
        if required_string(entry, "artifact_kind", metadata_path)? != kind
            || required_string(entry, "schema_version", metadata_path)? != "1"
            || required_string(entry, "path", metadata_path)? != name
        {
            return Err(format_error(
                metadata_path,
                "files",
                "orden, clase, versión o ruta no canónicos",
            ));
        }
        let size = required_u64(entry, "size_bytes", metadata_path)?;
        let digest = required_sha(entry, "sha256", metadata_path)?;
        let path = metadata_path
            .parent()
            .unwrap_or_else(|| Path::new("."))
            .join(name);
        let bytes = read(&path)?;
        verify_identity(&path, &bytes, Some(size), digest)?;
        contents.push((path, bytes));
    }
    let (block_lengths, block_lengths_audit) = parse_integer_intervals(
        &contents[0].0,
        &contents[0].1,
        "length_min_bp\tlength_max_bp\tweight",
        2,
        F5NormalizationSet::BlockLengths,
    )?;
    validate_continuous(&block_lengths, &contents[0].0)?;
    let l_max = block_lengths
        .last()
        .map(|interval| interval.maximum)
        .ok_or_else(|| domain(&contents[0].0, "el perfil de longitudes está vacío"))?;
    let (snp_gaps, snp_gaps_audit) = parse_integer_intervals(
        &contents[1].0,
        &contents[1].1,
        "gap_min_bp\tgap_max_bp\tweight",
        1,
        F5NormalizationSet::SnpGaps,
    )?;
    let (snp_spectrum, snp_spectrum_audit) = parse_float_intervals(&contents[2].0, &contents[2].1)?;
    let (haplotype_weights, haplotype_weights_audit) =
        parse_haplotype_weights(&contents[3].0, &contents[3].1)?;
    let ld_bins = parse_ld_profile(&contents[4].0, &contents[4].1, population_ids, l_max)?;
    Ok(F5Profile {
        profile_id,
        block_lengths,
        snp_gaps,
        snp_spectrum,
        haplotype_weights,
        ld_bins,
        normalizations: vec![
            block_lengths_audit,
            snp_gaps_audit,
            snp_spectrum_audit,
            haplotype_weights_audit,
        ],
    })
}

fn parse_integer_intervals(
    path: &Path,
    bytes: &[u8],
    header: &str,
    minimum: u64,
    set: F5NormalizationSet,
) -> Result<(Vec<WeightedInterval<u64>>, F5NormalizationAudit), Diagnostic> {
    let rows = tsv_rows(path, bytes, header)?;
    let mut result = Vec::with_capacity(rows.len());
    let mut keys = Vec::with_capacity(rows.len());
    for row in rows {
        let lo = canonical_u64(path, row[0])?;
        let hi = canonical_u64(path, row[1])?;
        let weight = canonical_f64(path, row[2])?;
        if lo < minimum || lo > hi || !is_positive_probability(weight) {
            return Err(domain(path, "intervalo o peso inválido"));
        }
        if result
            .last()
            .is_some_and(|previous: &WeightedInterval<u64>| lo <= previous.maximum)
        {
            return Err(domain(path, "filas duplicadas, solapadas o desordenadas"));
        }
        result.push(WeightedInterval {
            minimum: lo,
            maximum: hi,
            weight,
        });
        keys.push(row[0].to_owned());
    }
    let audit = normalize_weights(path, set, keys, &mut result)?;
    Ok((result, audit))
}

fn validate_continuous(intervals: &[WeightedInterval<u64>], path: &Path) -> Result<(), Diagnostic> {
    for pair in intervals.windows(2) {
        if pair[0].maximum.checked_add(1) != Some(pair[1].minimum) {
            return Err(domain(
                path,
                "los intervalos de longitudes no forman un soporte continuo",
            ));
        }
    }
    Ok(())
}

fn parse_float_intervals(
    path: &Path,
    bytes: &[u8],
) -> Result<(Vec<WeightedInterval<f64>>, F5NormalizationAudit), Diagnostic> {
    let rows = tsv_rows(path, bytes, "alt_af_min\talt_af_max\tweight")?;
    let mut result = Vec::with_capacity(rows.len());
    let mut keys = Vec::with_capacity(rows.len());
    for row in rows {
        let lo = canonical_f64(path, row[0])?;
        let hi = canonical_f64(path, row[1])?;
        let weight = canonical_f64(path, row[2])?;
        if !(lo.is_finite()
            && hi.is_finite()
            && is_probability(lo)
            && is_probability(hi)
            && lo <= hi
            && is_positive_probability(weight))
        {
            return Err(domain(path, "intervalo o peso inválido"));
        }
        if result
            .last()
            .is_some_and(|previous: &WeightedInterval<f64>| lo <= previous.maximum)
        {
            return Err(domain(path, "filas duplicadas, solapadas o desordenadas"));
        }
        result.push(WeightedInterval {
            minimum: lo,
            maximum: hi,
            weight,
        });
        keys.push(row[0].to_owned());
    }
    let audit = normalize_weights(path, F5NormalizationSet::SnpSpectrum, keys, &mut result)?;
    Ok((result, audit))
}

fn parse_haplotype_weights(
    path: &Path,
    bytes: &[u8],
) -> Result<(Vec<f64>, F5NormalizationAudit), Diagnostic> {
    let rows = tsv_rows(path, bytes, "haplotype_rank\tweight")?;
    if !(2..=4096).contains(&rows.len()) {
        return Err(domain(path, "H debe pertenecer a 2..=4096"));
    }
    let mut weights = Vec::with_capacity(rows.len());
    let mut keys = Vec::with_capacity(rows.len());
    for (rank, row) in rows.iter().enumerate() {
        if canonical_u64(path, row[0])? != rank as u64 {
            return Err(domain(path, "haplotype_rank no es consecutivo"));
        }
        let weight = canonical_f64(path, row[1])?;
        if !is_positive_probability(weight) {
            return Err(domain(path, "peso inválido"));
        }
        weights.push(weight);
        keys.push(row[0].to_owned());
    }
    let audit = normalize_slice(
        path,
        F5NormalizationSet::HaplotypeWeights,
        keys,
        &mut weights,
    )?;
    Ok((weights, audit))
}

fn parse_ld_profile(
    path: &Path,
    bytes: &[u8],
    populations: &BTreeSet<String>,
    l_max: u64,
) -> Result<BTreeMap<String, Vec<LdBin>>, Diagnostic> {
    let rows = tsv_rows(
        path,
        bytes,
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error",
    )?;
    let mut output = BTreeMap::<String, Vec<LdBin>>::new();
    let mut previous: Option<(String, u64)> = None;
    for row in rows {
        let population = valid_id(path, row[0])?.to_owned();
        let start = canonical_u64(path, row[1])?;
        let end = canonical_u64(path, row[2])?;
        let target = canonical_f64(path, row[3])?;
        let error = canonical_f64(path, row[4])?;
        if start == 0 || start > end || !is_probability(target) || error != 0.05 {
            return Err(domain(path, "bin LD inválido"));
        }
        if previous
            .as_ref()
            .is_some_and(|key| key >= &(population.clone(), start))
        {
            return Err(domain(path, "filas LD duplicadas o desordenadas"));
        }
        let bins = output.entry(population.clone()).or_default();
        if let Some(last) = bins.last() {
            if last.distance_end_bp.checked_add(1) != Some(start) {
                return Err(domain(path, "bins LD no contiguos"));
            }
        } else if start != 1 {
            return Err(domain(path, "bins LD deben empezar en 1"));
        }
        bins.push(LdBin {
            distance_start_bp: start,
            distance_end_bp: end,
            target_mean_r2: target,
        });
        previous = Some((population, start));
    }
    if output.len() != populations.len() || output.keys().any(|id| !populations.contains(id)) {
        return Err(domain(
            path,
            "population_profile_id no coincide exactamente con las subpoblaciones",
        ));
    }
    let expected_end = l_max
        .checked_sub(1)
        .ok_or_else(|| domain(path, "Lmax no permite distancias LD"))?;
    if output.values().any(|bins| {
        bins.last()
            .is_none_or(|last| last.distance_end_bp != expected_end)
    }) {
        return Err(domain(
            path,
            "cada población LD debe terminar exactamente en Lmax-1",
        ));
    }
    Ok(output)
}

fn tsv_rows<'a>(
    path: &Path,
    bytes: &'a [u8],
    expected_header: &str,
) -> Result<Vec<Vec<&'a str>>, Diagnostic> {
    if bytes.starts_with(&[0xef, 0xbb, 0xbf]) || !bytes.ends_with(b"\n") || bytes.contains(&b'\r') {
        return Err(format_error(
            path,
            "encoding",
            "se exige UTF-8 sin BOM, LF final y sin CR",
        ));
    }
    let text =
        std::str::from_utf8(bytes).map_err(|_| format_error(path, "encoding", "no es UTF-8"))?;
    let mut lines = text.strip_suffix('\n').unwrap_or(text).split('\n');
    if lines.next() != Some(expected_header) {
        return Err(format_error(path, "header", "cabecera TSV no canónica"));
    }
    let mut rows = Vec::new();
    for line in lines {
        if line.is_empty() || line.contains('"') || line.starts_with('#') {
            return Err(format_error(
                path,
                "row",
                "fila vacía o sintaxis TSV prohibida",
            ));
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        let expected_fields = expected_header.split('\t').count();
        if fields.len() != expected_fields || fields.iter().any(|field| field.is_empty()) {
            return Err(format_error(
                path,
                "row",
                "número de columnas o campo vacío inválido",
            ));
        }
        rows.push(fields);
    }
    if rows.is_empty() {
        return Err(format_error(path, "row", "el TSV debe contener filas"));
    }
    Ok(rows)
}

fn normalize_weights<T>(
    path: &Path,
    set: F5NormalizationSet,
    keys: Vec<String>,
    values: &mut [WeightedInterval<T>],
) -> Result<F5NormalizationAudit, Diagnostic> {
    let mut w = values.iter().map(|v| v.weight).collect::<Vec<_>>();
    let audit = normalize_slice(path, set, keys, &mut w)?;
    for (value, weight) in values.iter_mut().zip(w) {
        value.weight = weight;
    }
    Ok(audit)
}
fn normalize_slice(
    path: &Path,
    set: F5NormalizationSet,
    keys: Vec<String>,
    values: &mut [f64],
) -> Result<F5NormalizationAudit, Diagnostic> {
    let input_values = named_values(&keys, values);
    let input_sum = neumaier(values.iter().copied());
    let input_delta = (input_sum - 1.0).abs();
    if !(input_sum.is_finite()
        && input_sum > 0.0
        && input_delta <= F5_PROFILE_NORMALIZATION_TOLERANCE)
    {
        return Err(domain(path, "los pesos no se normalizan bajo DTC-049"));
    }
    for value in values.iter_mut() {
        *value /= input_sum;
    }
    let output_sum = neumaier(values.iter().copied());
    let output_delta = (output_sum - 1.0).abs();
    if !values.iter().copied().all(is_probability)
        || !output_sum.is_finite()
        || output_delta > F5_PROFILE_NORMALIZATION_TOLERANCE
    {
        return Err(domain(
            path,
            "la normalización DTC-049 produjo valores inválidos",
        ));
    }
    Ok(F5NormalizationAudit {
        set,
        input_values,
        input_sum,
        input_delta,
        effective_values: named_values(&keys, values),
        output_sum,
        output_delta,
    })
}
fn named_values(keys: &[String], values: &[f64]) -> Vec<F5NamedProbability> {
    keys.iter()
        .cloned()
        .zip(values.iter().copied())
        .map(|(key, value)| F5NamedProbability { key, value })
        .collect()
}
fn is_probability(value: f64) -> bool {
    value.is_finite() && (0.0..=1.0).contains(&value)
}
fn is_positive_probability(value: f64) -> bool {
    is_probability(value) && value > 0.0
}
fn neumaier(values: impl IntoIterator<Item = f64>) -> f64 {
    let mut sum = 0.0;
    let mut correction = 0.0;
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
fn canonical_u64(path: &Path, value: &str) -> Result<u64, Diagnostic> {
    if value.is_empty()
        || (value.len() > 1 && value.starts_with('0'))
        || !value.bytes().all(|b| b.is_ascii_digit())
    {
        return Err(format_error(path, "integer", "entero decimal no canónico"));
    }
    value
        .parse()
        .map_err(|_| domain(path, "entero fuera de u64"))
}
fn canonical_f64(path: &Path, value: &str) -> Result<f64, Diagnostic> {
    let parsed = value
        .parse::<f64>()
        .map_err(|_| format_error(path, "float", "f64 inválido"))?;
    if !parsed.is_finite() || parsed.to_string() != value {
        return Err(format_error(
            path,
            "float",
            "f64 no usa el decimal mínimo de ida y vuelta",
        ));
    }
    Ok(parsed)
}
fn read(path: &Path) -> Result<Vec<u8>, Diagnostic> {
    fs::read(path).map_err(|error| {
        Diagnostic::input(
            "GPS-INPUT-001",
            "No se puede abrir o leer una entrada declarada.",
        )
        .with("path", path.display().to_string())
        .with("detail", error.to_string())
    })
}
fn verify_identity(
    path: &Path,
    bytes: &[u8],
    expected_size: Option<u64>,
    expected_sha: &str,
) -> Result<(), Diagnostic> {
    let size = u64::try_from(bytes.len())
        .map_err(|_| format_error(path, "size_bytes", "tamaño no representable"))?;
    let digest = format!("{:x}", Sha256::digest(bytes));
    if expected_size.is_some_and(|expected| expected != size) || digest != expected_sha {
        return Err(format_error(
            path,
            "identity",
            format!("size={size}; sha256={digest}"),
        ));
    }
    Ok(())
}
fn object<'a>(
    value: &'a Value,
    path: &Path,
    field: &str,
) -> Result<&'a Map<String, Value>, Diagnostic> {
    value
        .as_object()
        .ok_or_else(|| format_error(path, field, "debe ser un objeto"))
}
fn exact_keys(
    object: &Map<String, Value>,
    expected: &[&str],
    path: &Path,
    field: &str,
) -> Result<(), Diagnostic> {
    let actual = object.keys().map(String::as_str).collect::<BTreeSet<_>>();
    let expected = expected.iter().copied().collect::<BTreeSet<_>>();
    (actual == expected)
        .then_some(())
        .ok_or_else(|| format_error(path, field, "faltan o sobran campos"))
}
fn required_string<'a>(
    object: &'a Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_str)
        .ok_or_else(|| format_error(path, field, "debe ser una cadena"))
}
fn nullable_string(
    object: &Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<(), Diagnostic> {
    if object
        .get(field)
        .is_some_and(|value| value.is_null() || value.as_str().is_some())
    {
        Ok(())
    } else {
        Err(format_error(path, field, "debe ser cadena o null"))
    }
}
fn nullable_sha(object: &Map<String, Value>, field: &str, path: &Path) -> Result<(), Diagnostic> {
    if object
        .get(field)
        .is_some_and(|value| value.is_null() || value.as_str().is_some_and(is_sha))
    {
        Ok(())
    } else {
        Err(format_error(path, field, "debe ser SHA-256 o null"))
    }
}
fn required_sha<'a>(
    object: &'a Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    let value = required_string(object, field, path)?;
    is_sha(value)
        .then_some(value)
        .ok_or_else(|| format_error(path, field, "debe ser SHA-256 hexadecimal minúsculo"))
}
fn required_u64(object: &Map<String, Value>, field: &str, path: &Path) -> Result<u64, Diagnostic> {
    object
        .get(field)
        .and_then(Value::as_u64)
        .ok_or_else(|| format_error(path, field, "debe ser u64"))
}
fn required_id<'a>(
    object: &'a Map<String, Value>,
    field: &str,
    path: &Path,
) -> Result<&'a str, Diagnostic> {
    valid_id(path, required_string(object, field, path)?)
}
fn valid_id<'a>(path: &Path, value: &'a str) -> Result<&'a str, Diagnostic> {
    let valid = (1..=64).contains(&value.len())
        && value
            .as_bytes()
            .first()
            .is_some_and(u8::is_ascii_alphanumeric)
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'.' | b'_' | b'-'));
    valid
        .then_some(value)
        .ok_or_else(|| domain(path, "ID inválido conforme a DTC-039"))
}
fn is_sha(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'a'..=b'f'))
}
fn format_error(path: &Path, field: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::format(
        "GPS-FMT-101",
        "El paquete de perfiles F5 no satisface su formato, versión, orden o identidad declarada.",
    )
    .with("path", path.display().to_string())
    .with("field", field)
    .with("detail", detail.into())
}
fn domain(path: &Path, detail: impl Into<String>) -> Diagnostic {
    format_error(path, "domain", detail)
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeSet;
    use std::fs;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};

    use sha2::{Digest, Sha256};

    use super::{F5NormalizationSet, load_f5_profile};

    static SERIAL: AtomicU64 = AtomicU64::new(0);

    fn fixture() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "gps-f5-profile-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&path).expect("fixture directory");
        write(
            &path,
            "block_lengths.tsv",
            "length_min_bp\tlength_max_bp\tweight\n2\t5\t0.5\n6\t10\t0.5\n",
        );
        write(
            &path,
            "snp_gaps.tsv",
            "gap_min_bp\tgap_max_bp\tweight\n1\t2\t0.5\n3\t4\t0.5\n",
        );
        write(
            &path,
            "snp_spectrum.tsv",
            "alt_af_min\talt_af_max\tweight\n0\t0.5\t0.5\n0.6\t1\t0.5\n",
        );
        write(
            &path,
            "haplotype_weights.tsv",
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        );
        write(
            &path,
            "ld_profile.tsv",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t9\t0.5\t0.05\n",
        );
        refresh_metadata(&path, "1", None, None);
        path
    }

    fn write(root: &Path, name: &str, content: &str) {
        fs::write(root.join(name), content).expect("fixture file");
    }
    fn digest(bytes: &[u8]) -> String {
        format!("{:x}", Sha256::digest(bytes))
    }
    fn refresh_metadata(
        root: &Path,
        version: &str,
        size_override: Option<u64>,
        sha_override: Option<&str>,
    ) {
        let names = [
            "block_lengths.tsv",
            "snp_gaps.tsv",
            "snp_spectrum.tsv",
            "haplotype_weights.tsv",
            "ld_profile.tsv",
        ];
        let kinds = [
            "f5_block_length_profile",
            "f5_snp_gap_profile",
            "f5_snp_spectrum",
            "f5_haplotype_weights",
            "f5_ld_profile",
        ];
        let files = names.iter().zip(kinds).enumerate().map(|(index, (name, kind))| {
            let bytes = fs::read(root.join(name)).expect("fixture bytes");
            let size = if index == 0 { size_override.unwrap_or(bytes.len() as u64) } else { bytes.len() as u64 };
            let sha = if index == 0 { sha_override.map_or_else(|| digest(&bytes), ToOwned::to_owned) } else { digest(&bytes) };
            format!(r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{size},"sha256":"{sha}"}}"#)
        }).collect::<Vec<_>>().join(",");
        fs::write(root.join("profile.metadata.json"), format!(r#"{{"profile_kind":"f5_generation_profile","schema_version":"{version}","profile_id":"test_profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{files}]}}"#)).expect("metadata");
    }
    fn load(root: &Path) -> Result<(), String> {
        let bytes = fs::read(root.join("profile.metadata.json")).expect("metadata bytes");
        let populations = BTreeSet::from(["P1".to_owned()]);
        load_f5_profile(
            &root.join("profile.metadata.json"),
            &digest(&bytes),
            &populations,
        )
        .map(|_| ())
        .map_err(|error| error.render_human())
    }

    #[test]
    fn reads_the_single_valid_canonical_profile_corpus() {
        let root = fixture();
        let profile = load_f5_profile(
            &root.join("profile.metadata.json"),
            &digest(&fs::read(root.join("profile.metadata.json")).expect("metadata")),
            &BTreeSet::from(["P1".to_owned()]),
        )
        .expect("valid profile");
        assert_eq!(profile.profile_id(), "test_profile");
        assert_eq!(profile.haplotype_weights(), &[0.5, 0.5]);
        assert_eq!(profile.ld_bins("P1").expect("P1")[0].distance_end_bp(), 9);
        assert_eq!(profile.normalizations().len(), 4);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn preserves_the_typed_dtc_049_audit_before_and_after_normalization() {
        let root = fixture();
        write(
            &root,
            "block_lengths.tsv",
            "length_min_bp\tlength_max_bp\tweight\n2\t5\t0.5000000001\n6\t10\t0.5\n",
        );
        refresh_metadata(&root, "1", None, None);

        let bytes = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let profile = load_f5_profile(
            &root.join("profile.metadata.json"),
            &digest(&bytes),
            &BTreeSet::from(["P1".to_owned()]),
        )
        .expect("normalizable profile");
        let audit = &profile.normalizations()[0];
        assert_eq!(audit.set(), F5NormalizationSet::BlockLengths);
        assert_eq!(audit.set().canonical_order(), "length_min_bp");
        assert_eq!(audit.algorithm(), "neumaier-f64-v1");
        assert_eq!(
            audit
                .input_values()
                .iter()
                .map(|value| value.key())
                .collect::<Vec<_>>(),
            ["2", "6"]
        );
        assert_eq!(audit.input_values()[0].value(), 0.500_000_000_1);
        assert_eq!(audit.input_sum(), 1.000_000_000_1);
        assert_eq!(audit.input_delta(), (audit.input_sum() - 1.0).abs());
        assert_ne!(
            audit.input_values()[0].value(),
            audit.effective_values()[0].value()
        );
        assert_eq!(audit.output_sum(), 1.0);
        assert_eq!(audit.output_delta(), 0.0);
        assert_eq!(
            profile.block_lengths()[0].weight(),
            audit.effective_values()[0].value()
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn rejects_every_non_probability_endpoint_and_weight_before_normalization() {
        type Mutation = (&'static str, fn(&Path));
        let cases: [Mutation; 8] = [
            ("block-weight", |root| {
                write(
                    root,
                    "block_lengths.tsv",
                    "length_min_bp\tlength_max_bp\tweight\n2\t10\t1.0000000005\n",
                );
            }),
            ("gap-weight", |root| {
                write(
                    root,
                    "snp_gaps.tsv",
                    "gap_min_bp\tgap_max_bp\tweight\n1\t4\t1.0000000005\n",
                );
            }),
            ("spectrum-min-nan", |root| {
                write(
                    root,
                    "snp_spectrum.tsv",
                    "alt_af_min\talt_af_max\tweight\nNaN\t1\t1\n",
                );
            }),
            ("spectrum-max-infinite", |root| {
                write(
                    root,
                    "snp_spectrum.tsv",
                    "alt_af_min\talt_af_max\tweight\n0\tinf\t1\n",
                );
            }),
            ("spectrum-weight", |root| {
                write(
                    root,
                    "snp_spectrum.tsv",
                    "alt_af_min\talt_af_max\tweight\n0\t1\t1.0000000005\n",
                );
            }),
            ("haplotype-weight", |root| {
                write(
                    root,
                    "haplotype_weights.tsv",
                    "haplotype_rank\tweight\n0\t1.0000000005\n1\t0.0000000001\n",
                );
            }),
            ("ld-target-nan", |root| {
                write(
                    root,
                    "ld_profile.tsv",
                    "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t9\tNaN\t0.05\n",
                );
            }),
            ("ld-target-over-one", |root| {
                write(
                    root,
                    "ld_profile.tsv",
                    "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t9\t1.1\t0.05\n",
                );
            }),
        ];
        for (name, mutate) in cases {
            let root = fixture();
            mutate(&root);
            refresh_metadata(&root, "1", None, None);
            assert!(load(&root).is_err(), "{name} must be rejected");
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    #[test]
    fn each_population_ld_range_must_end_at_lmax_minus_one() {
        let root = fixture();
        write(
            &root,
            "ld_profile.tsv",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t8\t0.5\t0.05\n",
        );
        refresh_metadata(&root, "1", None, None);
        let error = load(&root).expect_err("truncated LD range");
        assert!(error.contains("Lmax-1"), "{error}");
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn rejects_each_isolated_profile_mutation() {
        type Mutation = (&'static str, fn(&Path));
        let cases: [Mutation; 10] = [
            ("bom", |root| {
                let path = root.join("block_lengths.tsv");
                let mut bytes = vec![0xef, 0xbb, 0xbf];
                bytes.extend(fs::read(&path).expect("file"));
                fs::write(path, bytes).expect("write");
                refresh_metadata(root, "1", None, None);
            }),
            ("lf", |root| {
                let path = root.join("block_lengths.tsv");
                let mut bytes = fs::read(&path).expect("file");
                bytes.pop();
                fs::write(path, bytes).expect("write");
                refresh_metadata(root, "1", None, None);
            }),
            ("header", |root| {
                write(root, "block_lengths.tsv", "bad\tbad\tbad\n2\t5\t1\n");
                refresh_metadata(root, "1", None, None);
            }),
            ("row", |root| {
                write(
                    root,
                    "block_lengths.tsv",
                    "length_min_bp\tlength_max_bp\tweight\n2\t5\t1\textra\n",
                );
                refresh_metadata(root, "1", None, None);
            }),
            ("interval", |root| {
                write(
                    root,
                    "block_lengths.tsv",
                    "length_min_bp\tlength_max_bp\tweight\n2\t5\t0.5\n5\t8\t0.5\n",
                );
                refresh_metadata(root, "1", None, None);
            }),
            ("weight", |root| {
                write(
                    root,
                    "haplotype_weights.tsv",
                    "haplotype_rank\tweight\n0\t0.8\n1\t0.8\n",
                );
                refresh_metadata(root, "1", None, None);
            }),
            ("population", |root| {
                write(
                    root,
                    "ld_profile.tsv",
                    "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP2\t1\t9\t0.5\t0.05\n",
                );
                refresh_metadata(root, "1", None, None);
            }),
            ("version", |root| refresh_metadata(root, "2", None, None)),
            ("size", |root| refresh_metadata(root, "1", Some(1), None)),
            ("sha", |root| {
                refresh_metadata(root, "1", None, Some(&"0".repeat(64)))
            }),
        ];
        for (name, mutate) in cases {
            let root = fixture();
            mutate(&root);
            let error = load(&root).expect_err(name);
            assert!(error.contains("GPS-FMT-101"), "{name}: {error}");
            fs::remove_dir_all(root).expect("cleanup");
        }
    }
}
