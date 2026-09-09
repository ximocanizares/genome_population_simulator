//! Lector estricto del paquete de perfiles F5 v1.
//!
//! Esta frontera solo admite los cinco recursos declarados por el metadata. No
//! decide geometría, marginales ni haplotipos: entrega tablas ya verificadas a
//! los pasos posteriores.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::path::{Path, PathBuf};

use crate::diagnostics::Diagnostic;
use crate::f5_admission::{
    AdmissionPhase, AdmittedInputIdentity, BoundedInput, F5ResourceLedger, InputRole, LedgerCharge,
    OperationCost, PhysicalInputIdentity, ResourceDimension, WorkingMemoryPeak,
    input_identity_mismatch_diagnostic,
};
use crate::reference_fasta::input_error;
use serde_json::{Map, Value};

const F5_PROFILE_RETAINED_ITEM_BYTES: u64 = 256;

/// Tolerancia general DTC-049 para pesos de un perfil F5.
pub const F5_PROFILE_NORMALIZATION_TOLERANCE: f64 = 1e-9;
pub const F5_PROFILE_NORMALIZATION_ALGORITHM: &str = "neumaier-f64-v1";
/// Los cinco enlaces de un perfil F5 son un conjunto cerrado y ordenado. El
/// preflight usa estos nombres para obtener tamaños físicos antes de abrir el
/// metadata o materializar los TSV.
pub const F5_PROFILE_TSV_FILENAMES: [&str; 5] = [
    "block_lengths.tsv",
    "snp_gaps.tsv",
    "snp_spectrum.tsv",
    "haplotype_weights.tsv",
    "ld_profile.tsv",
];

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

/// Resultado de la fase de perfil: el contenido tipado y las seis identidades
/// físicas, ya con su identidad declarada, que el receipt sella.
#[derive(Clone, Debug)]
pub(crate) struct LoadedF5Profile {
    profile: F5Profile,
    identities: Vec<AdmittedInputIdentity>,
}

impl LoadedF5Profile {
    #[cfg(test)]
    #[must_use]
    pub(crate) fn identities(&self) -> &[AdmittedInputIdentity] {
        &self.identities
    }

    #[must_use]
    pub(crate) fn into_parts(self) -> (F5Profile, Vec<AdmittedInputIdentity>) {
        (self.profile, self.identities)
    }
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
#[cfg(test)]
fn load_f5_profile(
    metadata_path: &Path,
    expected_metadata_sha256: &str,
    population_ids: &BTreeSet<String>,
) -> Result<F5Profile, Diagnostic> {
    let mut ledger = F5ResourceLedger::new(crate::f5_admission::LedgerLimits::new(
        OperationCost::new(u64::MAX),
        WorkingMemoryPeak::new(u64::MAX),
        crate::f5_admission::OutputReservation::new(u64::MAX),
    ));
    ledger.transition_to(AdmissionPhase::ProfileMetadata)?;
    load_f5_profile_under_budget(
        metadata_path,
        expected_metadata_sha256,
        population_ids,
        &mut ledger,
    )
    .map(|loaded| loaded.profile)
}

/// Lee el paquete canónico bajo la fase actual del ledger. No consulta
/// atributos del sistema: el límite del metadata se deriva antes de abrirlo y
/// cada tabla usa el tamaño declarado más un byte para clasificar una
/// discrepancia de identidad sin materializar una entrada arbitraria.
pub(crate) fn load_f5_profile_under_budget(
    metadata_path: &Path,
    expected_metadata_sha256: &str,
    population_ids: &BTreeSet<String>,
    ledger: &mut F5ResourceLedger,
) -> Result<LoadedF5Profile, Diagnostic> {
    let metadata_limit = ledger.remaining_input_bytes(InputRole::ProfileMetadata)?;
    let (profile_id, links, metadata_identity) = read_profile_metadata(
        metadata_path,
        expected_metadata_sha256,
        metadata_limit,
        ledger,
    )?;

    // Los cinco costes de lectura ya son conocidos al completar el metadata.
    // Se reservan juntos, antes de abrir la primera tabla; la memoria es el
    // máximo de una tabla porque los búferes no coexisten.
    let mut total_table_bytes = 0_u64;
    let mut maximum_table_transient_bytes = 0_u64;
    for link in &links {
        total_table_bytes = total_table_bytes.checked_add(link.size).ok_or_else(|| {
            ledger.reject_resource(
                InputRole::ProfileMetadata,
                ResourceDimension::OperationUnits,
                "overflow al sumar tamaños declarados de TSV",
            )
        })?;
        let transient_bytes = link.size.checked_mul(2).ok_or_else(|| {
            ledger.reject_resource(
                InputRole::ProfileMetadata,
                ResourceDimension::WorkingMemoryBytes,
                "overflow al reservar búfer y parser del TSV declarado",
            )
        })?;
        maximum_table_transient_bytes = maximum_table_transient_bytes.max(transient_bytes);
    }
    let table_operation_cost = total_table_bytes.checked_mul(2).ok_or_else(|| {
        ledger.reject_resource(
            InputRole::ProfileMetadata,
            ResourceDimension::OperationUnits,
            "overflow al reservar lectura/hash y parser de los TSV",
        )
    })?;
    ledger.charge(
        InputRole::ProfileMetadata,
        LedgerCharge::Operation(OperationCost::new(table_operation_cost)),
    )?;
    ledger.reserve_transient_memory(InputRole::ProfileMetadata, maximum_table_transient_bytes)?;
    let mut retained_rows = 0_u64;

    let mut identities = Vec::with_capacity(6);
    identities.push(metadata_identity);
    let (block_path, block_bytes) = read_profile_table(
        &links[0],
        InputRole::ProfileBlockLengths,
        AdmissionPhase::ProfileTableBlockLengths,
        ledger,
        &mut identities,
        retained_rows,
    )?;
    let (block_lengths, block_lengths_audit) = parse_integer_intervals(
        &block_path,
        &block_bytes,
        "length_min_bp\tlength_max_bp\tweight",
        2,
        F5NormalizationSet::BlockLengths,
        ProfileParseContext::new(ledger, InputRole::ProfileBlockLengths, &mut retained_rows),
    )?;
    validate_continuous(&block_lengths, &block_path)?;
    let l_max = block_lengths
        .last()
        .map(|interval| interval.maximum)
        .ok_or_else(|| domain(&block_path, "el perfil de longitudes está vacío"))?;
    drop(block_bytes);
    let (gap_path, gap_bytes) = read_profile_table(
        &links[1],
        InputRole::ProfileSnpGaps,
        AdmissionPhase::ProfileTableSnpGaps,
        ledger,
        &mut identities,
        retained_rows,
    )?;
    let (snp_gaps, snp_gaps_audit) = parse_integer_intervals(
        &gap_path,
        &gap_bytes,
        "gap_min_bp\tgap_max_bp\tweight",
        1,
        F5NormalizationSet::SnpGaps,
        ProfileParseContext::new(ledger, InputRole::ProfileSnpGaps, &mut retained_rows),
    )?;
    drop(gap_bytes);
    let (spectrum_path, spectrum_bytes) = read_profile_table(
        &links[2],
        InputRole::ProfileSnpSpectrum,
        AdmissionPhase::ProfileTableSnpSpectrum,
        ledger,
        &mut identities,
        retained_rows,
    )?;
    let (snp_spectrum, snp_spectrum_audit) = parse_float_intervals(
        &spectrum_path,
        &spectrum_bytes,
        ledger,
        InputRole::ProfileSnpSpectrum,
        &mut retained_rows,
    )?;
    drop(spectrum_bytes);
    let (weights_path, weights_bytes) = read_profile_table(
        &links[3],
        InputRole::ProfileHaplotypeWeights,
        AdmissionPhase::ProfileTableHaplotypeWeights,
        ledger,
        &mut identities,
        retained_rows,
    )?;
    let (haplotype_weights, haplotype_weights_audit) = parse_haplotype_weights(
        &weights_path,
        &weights_bytes,
        ledger,
        InputRole::ProfileHaplotypeWeights,
        &mut retained_rows,
    )?;
    drop(weights_bytes);
    let (ld_path, ld_bytes) = read_profile_table(
        &links[4],
        InputRole::ProfileLd,
        AdmissionPhase::ProfileTableLd,
        ledger,
        &mut identities,
        retained_rows,
    )?;
    let ld_bins = parse_ld_profile(
        &ld_path,
        &ld_bytes,
        population_ids,
        l_max,
        ledger,
        InputRole::ProfileLd,
        &mut retained_rows,
    )?;
    drop(ld_bytes);
    Ok(LoadedF5Profile {
        profile: F5Profile {
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
        },
        identities,
    })
}

#[derive(Clone, Debug)]
struct ProfileLink {
    path: PathBuf,
    size: u64,
    digest: String,
}

fn read_bounded(
    path: &Path,
    role: InputRole,
    phase: AdmissionPhase,
    limit: u64,
    ledger: &mut F5ResourceLedger,
) -> Result<crate::f5_admission::BoundedInputRead, Diagnostic> {
    let file = File::open(path).map_err(|error| input_error(path, error))?;
    let input = BoundedInput::new(role, phase, path.to_path_buf(), limit, file);
    if role == InputRole::ProfileMetadata {
        input.read_to_end_with_ledger_counting_operations(ledger)
    } else {
        input.read_to_end_with_ledger(Some(ledger))
    }
}

/// Construye los datos persistentes mínimos del metadata y deja caer tanto el
/// búfer JSON como `serde_json::Value` antes de que se abra el primer TSV.
fn read_profile_metadata(
    metadata_path: &Path,
    expected_metadata_sha256: &str,
    metadata_limit: u64,
    ledger: &mut F5ResourceLedger,
) -> Result<(String, Vec<ProfileLink>, AdmittedInputIdentity), Diagnostic> {
    let metadata = read_bounded(
        metadata_path,
        InputRole::ProfileMetadata,
        AdmissionPhase::ProfileMetadata,
        metadata_limit,
        ledger,
    )?;
    let (bytes, identity) = metadata.into_parts();
    verify_identity(metadata_path, &identity, None, expected_metadata_sha256)?;
    ledger.charge(
        InputRole::ProfileMetadata,
        LedgerCharge::WorkingMemoryPeak(WorkingMemoryPeak::new(
            identity
                .observed_size_bytes()
                .checked_mul(2)
                .ok_or_else(|| {
                    format_error(
                        metadata_path,
                        "size_bytes",
                        "overflow al reservar lectura/hash y parser del metadata",
                    )
                })?,
        )),
    )?;
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
    if !matches!(
        required_string(root, "provenance_kind", metadata_path)?,
        "built_in" | "user_supplied"
    ) {
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
        (F5_PROFILE_TSV_FILENAMES[0], "f5_block_length_profile"),
        (F5_PROFILE_TSV_FILENAMES[1], "f5_snp_gap_profile"),
        (F5_PROFILE_TSV_FILENAMES[2], "f5_snp_spectrum"),
        (F5_PROFILE_TSV_FILENAMES[3], "f5_haplotype_weights"),
        (F5_PROFILE_TSV_FILENAMES[4], "f5_ld_profile"),
    ];
    let mut links = Vec::with_capacity(5);
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
        links.push(ProfileLink {
            path: metadata_path
                .parent()
                .unwrap_or_else(|| Path::new("."))
                .join(name),
            size: required_u64(entry, "size_bytes", metadata_path)?,
            digest: required_sha(entry, "sha256", metadata_path)?.to_owned(),
        });
    }
    // El escenario declara el SHA del metadata, pero no su tamaño.
    let identity =
        AdmittedInputIdentity::new(identity, None, Some(expected_metadata_sha256.to_owned()))?;
    Ok((profile_id, links, identity))
}

fn read_profile_table(
    link: &ProfileLink,
    role: InputRole,
    phase: AdmissionPhase,
    ledger: &mut F5ResourceLedger,
    identities: &mut Vec<AdmittedInputIdentity>,
    retained_rows: u64,
) -> Result<(PathBuf, Vec<u8>), Diagnostic> {
    ledger.transition_to(phase)?;
    charge_profile_phase(ledger, role, retained_rows, link.size)?;
    let file = File::open(&link.path).map_err(|error| input_error(&link.path, error))?;
    let bounded = BoundedInput::new(role, phase, link.path.clone(), link.size, file)
        .read_to_end_with_ledger_and_expected_identity(ledger, link.size, &link.digest)?;
    let (bytes, identity) = bounded.into_parts();
    verify_identity(&link.path, &identity, Some(link.size), &link.digest)?;
    identities.push(AdmittedInputIdentity::new(
        identity,
        Some(link.size),
        Some(link.digest.clone()),
    )?);
    Ok((link.path.clone(), bytes))
}

fn parse_integer_intervals(
    path: &Path,
    bytes: &[u8],
    header: &str,
    minimum: u64,
    set: F5NormalizationSet,
    mut context: ProfileParseContext<'_>,
) -> Result<(Vec<WeightedInterval<u64>>, F5NormalizationAudit), Diagnostic> {
    let mut result = Vec::new();
    let mut keys = Vec::new();
    tsv_for_each_row(path, bytes, header, |row| {
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
        context.retain_row(bytes.len())?;
        result.push(WeightedInterval {
            minimum: lo,
            maximum: hi,
            weight,
        });
        keys.push(row[0].to_owned());
        Ok(())
    })?;
    // La auditoría normalizada crea vectores propios. Su retención se admite
    // antes de construirlos, no después de una asignación ya consumada.
    context.retain_rows(result.len(), bytes.len())?;
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
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &mut u64,
) -> Result<(Vec<WeightedInterval<f64>>, F5NormalizationAudit), Diagnostic> {
    let mut result = Vec::new();
    let mut keys = Vec::new();
    tsv_for_each_row(path, bytes, "alt_af_min\talt_af_max\tweight", |row| {
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
        retain_profile_row(ledger, role, retained_rows, bytes.len())?;
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
        Ok(())
    })?;
    retain_profile_rows(ledger, role, retained_rows, result.len(), bytes.len())?;
    let audit = normalize_weights(path, F5NormalizationSet::SnpSpectrum, keys, &mut result)?;
    Ok((result, audit))
}

fn parse_haplotype_weights(
    path: &Path,
    bytes: &[u8],
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &mut u64,
) -> Result<(Vec<f64>, F5NormalizationAudit), Diagnostic> {
    let mut weights = Vec::new();
    let mut keys = Vec::new();
    let mut rank = 0_usize;
    tsv_for_each_row(path, bytes, "haplotype_rank\tweight", |row| {
        if canonical_u64(path, row[0])? != rank as u64 {
            return Err(domain(path, "haplotype_rank no es consecutivo"));
        }
        retain_profile_row(ledger, role, retained_rows, bytes.len())?;
        let weight = canonical_f64(path, row[1])?;
        if !is_positive_probability(weight) {
            return Err(domain(path, "peso inválido"));
        }
        weights.push(weight);
        keys.push(row[0].to_owned());
        rank = rank.checked_add(1).ok_or_else(|| {
            Diagnostic::f5_resource("<unresolved>", "overflow de filas de haplotipos")
        })?;
        Ok(())
    })?;
    if !(2..=4096).contains(&weights.len()) {
        return Err(domain(path, "H debe pertenecer a 2..=4096"));
    }
    retain_profile_rows(ledger, role, retained_rows, weights.len(), bytes.len())?;
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
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &mut u64,
) -> Result<BTreeMap<String, Vec<LdBin>>, Diagnostic> {
    let mut output = BTreeMap::<String, Vec<LdBin>>::new();
    let mut previous: Option<(String, u64)> = None;
    tsv_for_each_row(
        path,
        bytes,
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error",
        |row| {
            let population = valid_id(path, row[0])?.to_owned();
            let start = canonical_u64(path, row[1])?;
            let end = canonical_u64(path, row[2])?;
            let target = canonical_f64(path, row[3])?;
            let error = canonical_f64(path, row[4])?;
            if start == 0 || start > end || !is_probability(target) || error != 0.05 {
                return Err(domain(path, "bin LD inválido"));
            }
            retain_profile_row(ledger, role, retained_rows, bytes.len())?;
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
            Ok(())
        },
    )?;
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

fn charge_profile_phase(
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: u64,
    table_bytes: u64,
) -> Result<(), Diagnostic> {
    let retained = retained_rows
        .checked_mul(F5_PROFILE_RETAINED_ITEM_BYTES)
        .ok_or_else(|| {
            ledger.reject_resource(
                role,
                ResourceDimension::WorkingMemoryBytes,
                "overflow de retención del perfil",
            )
        })?;
    let transient = table_bytes.checked_mul(2).ok_or_else(|| {
        ledger.reject_resource(
            role,
            ResourceDimension::WorkingMemoryBytes,
            "overflow de lectura y parser TSV",
        )
    })?;
    ledger.charge_memory_phase(role, retained, transient, retained)
}

struct ProfileParseContext<'a> {
    ledger: &'a mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &'a mut u64,
}

impl<'a> ProfileParseContext<'a> {
    fn new(ledger: &'a mut F5ResourceLedger, role: InputRole, retained_rows: &'a mut u64) -> Self {
        Self {
            ledger,
            role,
            retained_rows,
        }
    }

    fn retain_row(&mut self, table_bytes: usize) -> Result<(), Diagnostic> {
        retain_profile_row(self.ledger, self.role, self.retained_rows, table_bytes)
    }

    fn retain_rows(&mut self, rows: usize, table_bytes: usize) -> Result<(), Diagnostic> {
        retain_profile_rows(
            self.ledger,
            self.role,
            self.retained_rows,
            rows,
            table_bytes,
        )
    }
}

fn retain_profile_row(
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &mut u64,
    table_bytes: usize,
) -> Result<(), Diagnostic> {
    retain_profile_rows(ledger, role, retained_rows, 1, table_bytes)
}

fn retain_profile_rows(
    ledger: &mut F5ResourceLedger,
    role: InputRole,
    retained_rows: &mut u64,
    rows: usize,
    table_bytes: usize,
) -> Result<(), Diagnostic> {
    let rows = u64::try_from(rows)
        .map_err(|_| Diagnostic::f5_resource("<unresolved>", "overflow de filas de perfil"))?;
    let after = retained_rows.checked_add(rows).ok_or_else(|| {
        ledger.reject_resource(
            role,
            ResourceDimension::WorkingMemoryBytes,
            "overflow de filas retenidas del perfil",
        )
    })?;
    let before_bytes = retained_rows
        .checked_mul(F5_PROFILE_RETAINED_ITEM_BYTES)
        .ok_or_else(|| {
            ledger.reject_resource(
                role,
                ResourceDimension::WorkingMemoryBytes,
                "overflow de retención del perfil",
            )
        })?;
    let after_bytes = after
        .checked_mul(F5_PROFILE_RETAINED_ITEM_BYTES)
        .ok_or_else(|| {
            ledger.reject_resource(
                role,
                ResourceDimension::WorkingMemoryBytes,
                "overflow de retención del perfil",
            )
        })?;
    let transient = u64::try_from(table_bytes)
        .map_err(|_| Diagnostic::f5_resource("<unresolved>", "overflow de bytes TSV"))?
        .checked_mul(2)
        .ok_or_else(|| {
            ledger.reject_resource(
                role,
                ResourceDimension::WorkingMemoryBytes,
                "overflow de lectura y parser TSV",
            )
        })?;
    ledger.charge_memory_phase(role, before_bytes, transient, after_bytes)?;
    *retained_rows = after;
    Ok(())
}

fn tsv_for_each_row(
    path: &Path,
    bytes: &[u8],
    expected_header: &str,
    mut visit: impl FnMut(&[&str]) -> Result<(), Diagnostic>,
) -> Result<(), Diagnostic> {
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
    let expected_fields = expected_header.split('\t').count();
    const MAX_CONTRACTUAL_FIELDS: usize = 5;
    if expected_fields > MAX_CONTRACTUAL_FIELDS {
        return Err(Diagnostic::internal_invariant(
            "f5_profile",
            "la cabecera contractual excede cinco campos TSV",
        ));
    }
    let mut count = 0_u64;
    for line in lines {
        if line.is_empty() || line.contains('"') || line.starts_with('#') {
            return Err(format_error(
                path,
                "row",
                "fila vacía o sintaxis TSV prohibida",
            ));
        }
        let mut fields = [""; MAX_CONTRACTUAL_FIELDS];
        let mut field_count = 0_usize;
        for field in line.split('\t') {
            if field_count == expected_fields || field.is_empty() {
                return Err(format_error(
                    path,
                    "row",
                    "número de columnas o campo vacío inválido",
                ));
            }
            fields[field_count] = field;
            field_count += 1;
        }
        if field_count != expected_fields {
            return Err(format_error(
                path,
                "row",
                "número de columnas o campo vacío inválido",
            ));
        }
        visit(&fields[..expected_fields])?;
        count = count
            .checked_add(1)
            .ok_or_else(|| Diagnostic::f5_resource("<unresolved>", "overflow de filas TSV"))?;
    }
    if count == 0 {
        return Err(format_error(path, "row", "el TSV debe contener filas"));
    }
    Ok(())
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
fn verify_identity(
    path: &Path,
    identity: &PhysicalInputIdentity,
    expected_size: Option<u64>,
    expected_sha: &str,
) -> Result<(), Diagnostic> {
    let size = identity.observed_size_bytes();
    let digest = identity.observed_sha256();
    if expected_size.is_some_and(|expected| expected != size) || digest != expected_sha {
        let mut diagnostic = Diagnostic::input(
            "GPS-INPUT-002",
            "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
        )
        .with("path", path.display().to_string())
        .with("expected_digest", expected_sha)
        .with("observed_digest", digest);
        if let Some(expected_size) = expected_size {
            diagnostic =
                input_identity_mismatch_diagnostic(path, expected_size, expected_sha, size, digest);
        }
        return Err(diagnostic);
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
    use std::time::{SystemTime, UNIX_EPOCH};

    use sha2::{Digest, Sha256};

    use crate::f5_admission::{
        AdmissionPhase, F5ResourceLedger, LedgerLimits, OperationCost, OutputReservation,
        WorkingMemoryPeak,
    };

    use super::{
        F5_PROFILE_TSV_FILENAMES, F5NormalizationSet, load_f5_profile, load_f5_profile_under_budget,
    };

    static SERIAL: AtomicU64 = AtomicU64::new(0);

    fn fixture() -> PathBuf {
        let path = std::env::temp_dir().join(format!(
            "gps-f5-profile-{}-{}-{}",
            std::process::id(),
            SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("clock after epoch")
                .as_nanos(),
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

    fn ledger(operation_limit: u64, memory_limit: u64) -> F5ResourceLedger {
        let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(operation_limit),
            WorkingMemoryPeak::new(memory_limit),
            OutputReservation::new(u64::MAX),
        ));
        ledger
            .transition_to(AdmissionPhase::ProfileMetadata)
            .expect("profile follows bootstrap");
        ledger
    }

    fn profile_bytes(root: &Path) -> (u64, u64, u64) {
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let tables = [
            "block_lengths.tsv",
            "snp_gaps.tsv",
            "snp_spectrum.tsv",
            "haplotype_weights.tsv",
            "ld_profile.tsv",
        ];
        let (total, peak) = tables.iter().fold((0_u64, 0_u64), |(total, peak), name| {
            let size = fs::metadata(root.join(name)).expect("table").len();
            (total + size, peak.max(size))
        });
        (u64::try_from(metadata.len()).expect("size"), total, peak)
    }

    fn profile_memory_oracle(root: &Path) -> (u64, u64) {
        let metadata_bytes = fs::metadata(root.join("profile.metadata.json"))
            .expect("metadata")
            .len();
        let mut peak = 2 * metadata_bytes;
        let mut retained_rows = 0_u64;
        for (index, name) in F5_PROFILE_TSV_FILENAMES.iter().enumerate() {
            let bytes = fs::read(root.join(name)).expect("table oracle bytes");
            let row_count = std::str::from_utf8(&bytes)
                .expect("fixture UTF-8")
                .lines()
                .count()
                .checked_sub(1)
                .expect("header") as u64;
            retained_rows += row_count;
            if index < 4 {
                retained_rows += row_count;
            }
            peak = peak.max(2 * bytes.len() as u64 + 256 * retained_rows);
        }
        (peak, retained_rows * 256)
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
    fn profile_budget_is_exact_sequential_and_minus_one_stops_before_tables() {
        let root = fixture();
        let (metadata_bytes, table_bytes, _) = profile_bytes(&root);
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let populations = BTreeSet::from(["P1".to_owned()]);

        let operation_budget = 2 * (metadata_bytes + table_bytes);
        let (memory_budget, retained_bytes) = profile_memory_oracle(&root);
        let mut exact = ledger(operation_budget, memory_budget);
        let loaded = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &populations,
            &mut exact,
        )
        .expect("the inclusive budget admits all five declared tables");
        assert_eq!(loaded.identities().len(), 6);
        assert_eq!(exact.phase(), AdmissionPhase::ProfileTableLd);
        assert_eq!(exact.retained_memory_bytes(), retained_bytes);
        assert_eq!(exact.totals().operation_units().value(), operation_budget);
        assert_eq!(
            exact.totals().working_memory_peak().value(),
            memory_budget,
            "el pico incluye los dominios y auditorías ya retenidos"
        );

        let mut minus_one = ledger(operation_budget - 1, memory_budget);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &populations,
            &mut minus_one,
        )
        .expect_err("all known table costs are reserved before opening the first table");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(minus_one.phase(), AdmissionPhase::Rejected);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn metadata_and_table_oversize_are_rejected_by_the_derived_bounds() {
        let root = fixture();
        let (metadata_bytes, table_bytes, table_peak) = profile_bytes(&root);
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let populations = BTreeSet::from(["P1".to_owned()]);
        let mut metadata_too_large = ledger(u64::MAX, 2 * metadata_bytes - 1);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &populations,
            &mut metadata_too_large,
        )
        .expect_err("metadata cannot reserve a byte over the remaining budget");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(metadata_too_large.phase(), AdmissionPhase::Rejected);

        refresh_metadata(&root, "1", Some(1), None);
        let undersized = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let mut table_too_large = ledger(
            2 * (metadata_bytes + table_bytes),
            2 * table_peak.max(metadata_bytes),
        );
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&undersized),
            &populations,
            &mut table_too_large,
        )
        .expect_err("the declared size plus one is the physical table bound");
        assert_eq!(error.code(), "GPS-INPUT-002");
        assert_eq!(table_too_large.phase(), AdmissionPhase::Rejected);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn tsv_dominant_memory_limit_is_exact_and_minus_one_with_live_retention() {
        let root = fixture();
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let populations = BTreeSet::from(["P1".to_owned()]);
        let (required, retained_bytes) = profile_memory_oracle(&root);
        assert!(
            required > 2 * metadata.len() as u64,
            "a TSV phase dominates"
        );
        let mut exact = ledger(u64::MAX, required);
        load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &populations,
            &mut exact,
        )
        .expect("the inclusive profile phase peak is admitted");
        assert_eq!(exact.retained_memory_bytes(), retained_bytes);

        let mut minus_one = ledger(u64::MAX, required - 1);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &populations,
            &mut minus_one,
        )
        .expect_err("the retained profile plus parser do not fit");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(minus_one.phase(), AdmissionPhase::Rejected);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn known_later_tsv_transient_is_reserved_before_opening_the_first_table() {
        let root = fixture();
        let metadata_path = root.join("profile.metadata.json");
        let mut metadata: serde_json::Value =
            serde_json::from_slice(&fs::read(&metadata_path).expect("metadata"))
                .expect("metadata JSON");
        let oversized_later_table = 10_000_u64;
        metadata["files"][1]["size_bytes"] = oversized_later_table.into();
        fs::write(
            &metadata_path,
            serde_json::to_vec(&metadata).expect("serialized metadata"),
        )
        .expect("metadata rewrite");
        fs::remove_file(root.join("block_lengths.tsv")).expect("first-table opening spy");
        let metadata = fs::read(&metadata_path).expect("rewritten metadata");
        let populations = BTreeSet::from(["P1".to_owned()]);
        let required_transient = 2 * oversized_later_table;

        let mut minus_one = ledger(u64::MAX, required_transient - 1);
        let error = load_f5_profile_under_budget(
            &metadata_path,
            &digest(&metadata),
            &populations,
            &mut minus_one,
        )
        .expect_err("the known later transient rejects before the absent first table opens");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        assert_eq!(error.context()["phase"], "profile_metadata");
        assert_eq!(error.context()["role"], "profile_metadata");
        assert_eq!(minus_one.phase(), AdmissionPhase::Rejected);

        let mut exact = ledger(u64::MAX, required_transient);
        let error = load_f5_profile_under_budget(
            &metadata_path,
            &digest(&metadata),
            &populations,
            &mut exact,
        )
        .expect_err("the exact reservation reaches the deliberately absent first table");
        assert_eq!(error.code(), "GPS-INPUT-001");
        assert_eq!(
            exact.totals().working_memory_peak().value(),
            required_transient,
            "the exact prospective peak is reserved without retained rows"
        );
        assert_eq!(exact.retained_memory_bytes(), 0);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn excessive_tsv_columns_use_a_fixed_field_array_under_a_reduced_limit() {
        let root = fixture();
        let malicious = format!("{}1\n", "1\t".repeat(200_000));
        write(
            &root,
            "block_lengths.tsv",
            &format!("length_min_bp\tlength_max_bp\tweight\n{malicious}"),
        );
        refresh_metadata(&root, "1", None, None);
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let table_bytes = fs::metadata(root.join("block_lengths.tsv"))
            .expect("table")
            .len();
        let memory_limit = 2 * table_bytes;
        let mut bounded = ledger(u64::MAX, memory_limit);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &BTreeSet::from(["P1".to_owned()]),
            &mut bounded,
        )
        .expect_err("the fourth field rejects a three-column table");
        assert_eq!(error.code(), "GPS-FMT-101");
        assert_eq!(bounded.totals().working_memory_peak().value(), memory_limit);
        assert_eq!(bounded.retained_memory_bytes(), 0);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn normalization_audit_retention_is_charged_before_allocation() {
        let root = fixture();
        let mut table = String::from("length_min_bp\tlength_max_bp\tweight\n");
        for value in 2..=101 {
            table.push_str(&format!("{value}\t{value}\t0.001\n"));
        }
        write(&root, "block_lengths.tsv", &table);
        refresh_metadata(&root, "1", None, None);
        fs::remove_file(root.join("snp_gaps.tsv")).expect("later-open spy");
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let table_bytes = fs::metadata(root.join("block_lengths.tsv"))
            .expect("table")
            .len();
        let rows_only_limit = 2 * table_bytes + 100 * 256;
        let mut bounded = ledger(u64::MAX, rows_only_limit);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &BTreeSet::from(["P1".to_owned()]),
            &mut bounded,
        )
        .expect_err("audit vectors require a prospective reservation");
        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(bounded.phase(), AdmissionPhase::Rejected);
        assert_eq!(
            error.context().get("dimension"),
            Some(&"working_memory_bytes".to_owned())
        );
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn first_table_rejection_does_not_open_a_later_table() {
        let root = fixture();
        write(
            &root,
            "block_lengths.tsv",
            "length_min_bp\tlength_max_bp\tweight\n2\t10\t0\n",
        );
        refresh_metadata(&root, "1", None, None);
        fs::remove_file(root.join("snp_gaps.tsv")).expect("later table is an opening spy");
        let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
        let mut ledger = ledger(u64::MAX, u64::MAX);
        let error = load_f5_profile_under_budget(
            &root.join("profile.metadata.json"),
            &digest(&metadata),
            &BTreeSet::from(["P1".to_owned()]),
            &mut ledger,
        )
        .expect_err("the first parser rejection stops before the absent second table");
        assert_eq!(error.code(), "GPS-FMT-101");
        assert!(error.render_human().contains("block_lengths.tsv"));
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn profile_identity_contract_uses_the_six_ordered_observations() {
        let root = fixture();
        let metadata_path = root.join("profile.metadata.json");
        let metadata = fs::read(&metadata_path).expect("metadata");
        let mut ledger = ledger(u64::MAX, u64::MAX);
        let loaded = load_f5_profile_under_budget(
            &metadata_path,
            &digest(&metadata),
            &BTreeSet::from(["P1".to_owned()]),
            &mut ledger,
        )
        .expect("valid profile");
        let roles = [
            crate::f5_admission::InputRole::ProfileMetadata,
            crate::f5_admission::InputRole::ProfileBlockLengths,
            crate::f5_admission::InputRole::ProfileSnpGaps,
            crate::f5_admission::InputRole::ProfileSnpSpectrum,
            crate::f5_admission::InputRole::ProfileHaplotypeWeights,
            crate::f5_admission::InputRole::ProfileLd,
        ];
        let names = [
            "profile.metadata.json",
            "block_lengths.tsv",
            "snp_gaps.tsv",
            "snp_spectrum.tsv",
            "haplotype_weights.tsv",
            "ld_profile.tsv",
        ];
        for ((identity, role), name) in loaded.identities().iter().zip(roles).zip(names) {
            let bytes = fs::read(root.join(name)).expect("literal identity oracle");
            assert_eq!(identity.role(), role);
            assert_eq!(identity.normalized_path(), root.join(name));
            assert_eq!(identity.observed_size_bytes(), bytes.len() as u64);
            assert_eq!(identity.observed_sha256(), digest(&bytes));
        }
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn each_tsv_excess_is_an_identity_rejection_before_the_next_open() {
        for name in F5_PROFILE_TSV_FILENAMES {
            let root = fixture();
            let mut bytes = fs::read(root.join(name)).expect("table");
            let declared_size = bytes.len().to_string();
            bytes.extend(std::iter::repeat_n(b'x', 8_193));
            fs::write(root.join(name), bytes).expect("oversized table");
            let metadata = fs::read(root.join("profile.metadata.json")).expect("metadata");
            let mut ledger = ledger(u64::MAX, u64::MAX);
            let error = load_f5_profile_under_budget(
                &root.join("profile.metadata.json"),
                &digest(&metadata),
                &BTreeSet::from(["P1".to_owned()]),
                &mut ledger,
            )
            .expect_err("a table larger than its declared identity is rejected");
            assert_eq!(error.code(), "GPS-INPUT-002", "{name}");
            assert_eq!(
                error.context().get("expected_size_bytes"),
                Some(&declared_size),
                "{name}"
            );
            assert!(error.render_human().contains(name));
            fs::remove_dir_all(root).expect("cleanup");
        }
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
            if matches!(name, "size" | "sha") {
                assert!(error.contains("GPS-INPUT-002"), "{name}: {error}");
            } else {
                assert!(error.contains("GPS-FMT-101"), "{name}: {error}");
            }
            fs::remove_dir_all(root).expect("cleanup");
        }
    }
}
