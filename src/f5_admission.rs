//! Frontera física y contabilidad incremental de la admisión F5 v7.
//!
//! Este módulo no conoce parsers de perfil, referencia ni layout. Su contrato
//! es deliberadamente pequeño: un mismo handle cuenta, hashea y acota una
//! entrada, y el ledger rechaza prospectivamente costes acumulados.

use std::fmt;
use std::io::Read;
use std::path::{Path, PathBuf};

use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;

pub const F5_RESOURCE_ACCOUNTING_ALGORITHM: &str = "f5-resource-accounting-v7";

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum InputRole {
    GeneralToml,
    ScenarioToml,
    ProfileMetadata,
    ProfileBlockLengths,
    ProfileSnpGaps,
    ProfileSnpSpectrum,
    ProfileHaplotypeWeights,
    ProfileLd,
    ReferenceFasta,
    Mask,
}

impl InputRole {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::GeneralToml => "general_toml",
            Self::ScenarioToml => "scenario_toml",
            Self::ProfileMetadata => "profile_metadata",
            Self::ProfileBlockLengths => "profile_block_lengths",
            Self::ProfileSnpGaps => "profile_snp_gaps",
            Self::ProfileSnpSpectrum => "profile_snp_spectrum",
            Self::ProfileHaplotypeWeights => "profile_haplotype_weights",
            Self::ProfileLd => "profile_ld",
            Self::ReferenceFasta => "reference_fasta",
            Self::Mask => "mask",
        }
    }
}

impl fmt::Display for InputRole {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum AdmissionPhase {
    Bootstrap,
    ProfileMetadata,
    ProfileTableBlockLengths,
    ProfileTableSnpGaps,
    ProfileTableSnpSpectrum,
    ProfileTableHaplotypeWeights,
    ProfileTableLd,
    ReferenceAndMaskScan,
    FinalFormula,
    Sealed,
    Rejected,
}

impl AdmissionPhase {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Bootstrap => "bootstrap",
            Self::ProfileMetadata => "profile_metadata",
            Self::ProfileTableBlockLengths => "profile_table_block_lengths",
            Self::ProfileTableSnpGaps => "profile_table_snp_gaps",
            Self::ProfileTableSnpSpectrum => "profile_table_snp_spectrum",
            Self::ProfileTableHaplotypeWeights => "profile_table_haplotype_weights",
            Self::ProfileTableLd => "profile_table_ld",
            Self::ReferenceAndMaskScan => "reference_and_mask_scan",
            Self::FinalFormula => "final_formula",
            Self::Sealed => "sealed",
            Self::Rejected => "rejected",
        }
    }
}

impl fmt::Display for AdmissionPhase {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ResourceDimension {
    OperationUnits,
    WorkingMemoryBytes,
    OutputBytes,
    PhysicalInputBytes,
}

impl ResourceDimension {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::OperationUnits => "operation_units",
            Self::WorkingMemoryBytes => "working_memory_bytes",
            Self::OutputBytes => "output_bytes",
            Self::PhysicalInputBytes => "physical_input_bytes",
        }
    }
}

impl fmt::Display for ResourceDimension {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

macro_rules! resource_amount {
    ($name:ident) => {
        #[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
        pub struct $name(u64);

        impl $name {
            #[must_use]
            pub const fn new(value: u64) -> Self {
                Self(value)
            }

            #[must_use]
            pub const fn value(self) -> u64 {
                self.0
            }
        }
    };
}

resource_amount!(OperationCost);
resource_amount!(WorkingMemoryPeak);
resource_amount!(OutputReservation);

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LedgerLimits {
    operation_units: OperationCost,
    working_memory_peak: WorkingMemoryPeak,
    output_bytes: OutputReservation,
}

impl LedgerLimits {
    #[must_use]
    pub const fn new(
        operation_units: OperationCost,
        working_memory_peak: WorkingMemoryPeak,
        output_bytes: OutputReservation,
    ) -> Self {
        Self {
            operation_units,
            working_memory_peak,
            output_bytes,
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct LedgerTotals {
    operation_units: OperationCost,
    working_memory_peak: WorkingMemoryPeak,
    output_bytes: OutputReservation,
}

impl LedgerTotals {
    #[must_use]
    pub const fn operation_units(self) -> OperationCost {
        self.operation_units
    }

    #[must_use]
    pub const fn working_memory_peak(self) -> WorkingMemoryPeak {
        self.working_memory_peak
    }

    #[must_use]
    pub const fn output_bytes(self) -> OutputReservation {
        self.output_bytes
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LedgerCharge {
    Operation(OperationCost),
    WorkingMemoryPeak(WorkingMemoryPeak),
    Output(OutputReservation),
}

/// Ledger acumulativo de la máquina de estados v7.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct F5ResourceLedger {
    limits: LedgerLimits,
    totals: LedgerTotals,
    phase: AdmissionPhase,
}

impl F5ResourceLedger {
    #[must_use]
    pub const fn new(limits: LedgerLimits) -> Self {
        Self {
            limits,
            totals: LedgerTotals {
                operation_units: OperationCost::new(0),
                working_memory_peak: WorkingMemoryPeak::new(0),
                output_bytes: OutputReservation::new(0),
            },
            phase: AdmissionPhase::Bootstrap,
        }
    }

    #[must_use]
    pub const fn totals(&self) -> LedgerTotals {
        self.totals
    }

    #[must_use]
    pub const fn phase(&self) -> AdmissionPhase {
        self.phase
    }

    /// Cambia a la siguiente fase de la máquina de admisión; no acepta saltos.
    pub fn transition_to(&mut self, next: AdmissionPhase) -> Result<(), Diagnostic> {
        if self.phase == AdmissionPhase::Rejected {
            return Err(terminal_rejection_diagnostic());
        }
        if next == self.phase {
            return Ok(());
        }
        if next != next_phase(self.phase) {
            return Err(resource_diagnostic(
                self.phase,
                InputRole::GeneralToml,
                ResourceDimension::OperationUnits,
                0,
                0,
                "transición de admisión fuera de orden",
            ));
        }
        self.phase = next;
        Ok(())
    }

    /// Comprueba primero el nuevo total y solo después muta el ledger.
    pub fn charge(&mut self, role: InputRole, charge: LedgerCharge) -> Result<(), Diagnostic> {
        if self.phase == AdmissionPhase::Rejected {
            return Err(terminal_rejection_diagnostic());
        }

        let prospective = match charge {
            LedgerCharge::Operation(cost) => {
                let prospective = self
                    .totals
                    .operation_units
                    .value()
                    .checked_add(cost.value())
                    .ok_or_else(|| {
                        overflow_diagnostic(
                            self.phase,
                            role,
                            ResourceDimension::OperationUnits,
                            self.limits.operation_units.value(),
                        )
                    });
                prospective.map(|prospective| {
                    (
                        ResourceDimension::OperationUnits,
                        prospective,
                        self.limits.operation_units.value(),
                        LedgerTotals {
                            operation_units: OperationCost::new(prospective),
                            ..self.totals
                        },
                    )
                })
            }
            LedgerCharge::WorkingMemoryPeak(peak) => {
                let prospective = self.totals.working_memory_peak.value().max(peak.value());
                Ok((
                    ResourceDimension::WorkingMemoryBytes,
                    prospective,
                    self.limits.working_memory_peak.value(),
                    LedgerTotals {
                        working_memory_peak: WorkingMemoryPeak::new(prospective),
                        ..self.totals
                    },
                ))
            }
            LedgerCharge::Output(reservation) => {
                let prospective = self
                    .totals
                    .output_bytes
                    .value()
                    .checked_add(reservation.value())
                    .ok_or_else(|| {
                        overflow_diagnostic(
                            self.phase,
                            role,
                            ResourceDimension::OutputBytes,
                            self.limits.output_bytes.value(),
                        )
                    });
                prospective.map(|prospective| {
                    (
                        ResourceDimension::OutputBytes,
                        prospective,
                        self.limits.output_bytes.value(),
                        LedgerTotals {
                            output_bytes: OutputReservation::new(prospective),
                            ..self.totals
                        },
                    )
                })
            }
        };
        let (dimension, accumulated, limit, updated) = match prospective {
            Ok(prospective) => prospective,
            Err(diagnostic) => return Err(self.reject(diagnostic)),
        };
        if accumulated > limit {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                role,
                dimension,
                accumulated,
                limit,
                "reserva prospectiva excede el límite inclusivo",
            )));
        }
        self.totals = updated;
        Ok(())
    }

    fn reject(&mut self, diagnostic: Diagnostic) -> Diagnostic {
        self.phase = AdmissionPhase::Rejected;
        diagnostic
    }

    fn restore_totals(&mut self, totals: LedgerTotals) {
        self.totals = totals;
    }
}

const fn next_phase(phase: AdmissionPhase) -> AdmissionPhase {
    match phase {
        AdmissionPhase::Bootstrap => AdmissionPhase::ProfileMetadata,
        AdmissionPhase::ProfileMetadata => AdmissionPhase::ProfileTableBlockLengths,
        AdmissionPhase::ProfileTableBlockLengths => AdmissionPhase::ProfileTableSnpGaps,
        AdmissionPhase::ProfileTableSnpGaps => AdmissionPhase::ProfileTableSnpSpectrum,
        AdmissionPhase::ProfileTableSnpSpectrum => AdmissionPhase::ProfileTableHaplotypeWeights,
        AdmissionPhase::ProfileTableHaplotypeWeights => AdmissionPhase::ProfileTableLd,
        AdmissionPhase::ProfileTableLd => AdmissionPhase::ReferenceAndMaskScan,
        AdmissionPhase::ReferenceAndMaskScan => AdmissionPhase::FinalFormula,
        AdmissionPhase::FinalFormula => AdmissionPhase::Sealed,
        AdmissionPhase::Sealed => AdmissionPhase::Sealed,
        AdmissionPhase::Rejected => AdmissionPhase::Rejected,
    }
}

fn terminal_rejection_diagnostic() -> Diagnostic {
    resource_diagnostic(
        AdmissionPhase::Rejected,
        InputRole::GeneralToml,
        ResourceDimension::OperationUnits,
        0,
        0,
        "el ledger ya está rechazado terminalmente",
    )
}

const fn is_bootstrap_role(role: InputRole) -> bool {
    matches!(role, InputRole::GeneralToml | InputRole::ScenarioToml)
}

const fn role_matches_phase(role: InputRole, phase: AdmissionPhase) -> bool {
    matches!(
        (role, phase),
        (
            InputRole::GeneralToml | InputRole::ScenarioToml,
            AdmissionPhase::Bootstrap
        ) | (InputRole::ProfileMetadata, AdmissionPhase::ProfileMetadata)
            | (
                InputRole::ProfileBlockLengths,
                AdmissionPhase::ProfileTableBlockLengths
            )
            | (
                InputRole::ProfileSnpGaps,
                AdmissionPhase::ProfileTableSnpGaps
            )
            | (
                InputRole::ProfileSnpSpectrum,
                AdmissionPhase::ProfileTableSnpSpectrum
            )
            | (
                InputRole::ProfileHaplotypeWeights,
                AdmissionPhase::ProfileTableHaplotypeWeights
            )
            | (InputRole::ProfileLd, AdmissionPhase::ProfileTableLd)
            | (
                InputRole::ReferenceFasta | InputRole::Mask,
                AdmissionPhase::ReferenceAndMaskScan
            )
    )
}

fn resource_diagnostic(
    phase: AdmissionPhase,
    role: InputRole,
    dimension: ResourceDimension,
    accumulated: u64,
    limit: u64,
    detail: &'static str,
) -> Diagnostic {
    Diagnostic::f5_resource("<unresolved>", detail)
        .with("algorithm", F5_RESOURCE_ACCOUNTING_ALGORITHM)
        .with("phase", phase.as_str())
        .with("role", role.as_str())
        .with("dimension", dimension.as_str())
        .with("accumulated", accumulated.to_string())
        .with("limit", limit.to_string())
}

fn overflow_diagnostic(
    phase: AdmissionPhase,
    role: InputRole,
    dimension: ResourceDimension,
    limit: u64,
) -> Diagnostic {
    resource_diagnostic(
        phase,
        role,
        dimension,
        u64::MAX,
        limit,
        "overflow aritmético prospectivo",
    )
}

/// Identidad que une una entrada a los bytes físicamente observados.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct PhysicalInputIdentity {
    role: InputRole,
    normalized_path: PathBuf,
    observed_size_bytes: u64,
    observed_sha256: String,
}

impl PhysicalInputIdentity {
    #[must_use]
    pub const fn role(&self) -> InputRole {
        self.role
    }

    #[must_use]
    pub fn normalized_path(&self) -> &Path {
        &self.normalized_path
    }

    #[must_use]
    pub const fn observed_size_bytes(&self) -> u64 {
        self.observed_size_bytes
    }

    #[must_use]
    pub fn observed_sha256(&self) -> &str {
        &self.observed_sha256
    }
}

/// Resultado EOF de una lectura física acotada.
#[derive(Debug)]
pub struct BoundedInputRead {
    bytes: Vec<u8>,
    identity: PhysicalInputIdentity,
}

impl BoundedInputRead {
    #[must_use]
    pub fn bytes(&self) -> &[u8] {
        &self.bytes
    }

    #[must_use]
    pub fn into_bytes(self) -> Vec<u8> {
        self.bytes
    }

    #[must_use]
    pub const fn identity(&self) -> &PhysicalInputIdentity {
        &self.identity
    }

    #[must_use]
    pub fn into_parts(self) -> (Vec<u8>, PhysicalInputIdentity) {
        (self.bytes, self.identity)
    }
}

/// Handle que mide, hashea y limita exactamente los bytes que entrega.
pub struct BoundedInput<R> {
    role: InputRole,
    phase: AdmissionPhase,
    normalized_path: PathBuf,
    limit: u64,
    reader: R,
}

impl<R: Read> BoundedInput<R> {
    #[must_use]
    pub fn new(
        role: InputRole,
        phase: AdmissionPhase,
        normalized_path: PathBuf,
        limit: u64,
        reader: R,
    ) -> Self {
        Self {
            role,
            phase,
            normalized_path,
            limit,
            reader,
        }
    }

    /// Lee hasta EOF. Al alcanzar la cota, observa solamente un byte adicional.
    pub fn read_to_end(self) -> Result<BoundedInputRead, Diagnostic> {
        self.read_to_end_with_ledger(None)
    }

    /// Variante que reserva prospectivamente la memoria lógica de cada aumento
    /// del búfer, antes de que `Vec` solicite nueva capacidad.
    pub fn read_to_end_with_ledger(
        self,
        ledger: Option<&mut F5ResourceLedger>,
    ) -> Result<BoundedInputRead, Diagnostic> {
        self.read_to_end_with_reserver(ledger, |bytes, additional| {
            bytes.try_reserve(additional).is_ok()
        })
    }

    fn read_to_end_with_reserver(
        self,
        ledger: Option<&mut F5ResourceLedger>,
        mut reserve: impl FnMut(&mut Vec<u8>, usize) -> bool,
    ) -> Result<BoundedInputRead, Diagnostic> {
        self.validate_read_admission(ledger.as_deref())?;
        let Self {
            role,
            phase,
            normalized_path,
            limit,
            mut reader,
        } = self;
        let mut bytes = Vec::new();
        let mut observed = 0_u64;
        let mut digest = Sha256::new();
        let mut chunk = [0_u8; 8_192];
        let mut ledger = ledger;

        loop {
            let remaining = match limit.checked_sub(observed) {
                Some(remaining) => remaining,
                None => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::PhysicalInputBytes,
                            limit,
                        ),
                    ));
                }
            };
            if remaining == 0 {
                let mut excess = [0_u8; 1];
                let read = reader
                    .read(&mut excess)
                    .map_err(|error| input_read_error(&normalized_path, error))?;
                if read != 0 {
                    let accumulated = match observed.checked_add(1) {
                        Some(accumulated) => accumulated,
                        None => {
                            return Err(reject_resource_error(
                                &mut ledger,
                                overflow_diagnostic(
                                    phase,
                                    role,
                                    ResourceDimension::PhysicalInputBytes,
                                    limit,
                                ),
                            ));
                        }
                    };
                    return Err(reject_resource_error(
                        &mut ledger,
                        resource_diagnostic(
                            phase,
                            role,
                            ResourceDimension::PhysicalInputBytes,
                            accumulated,
                            limit,
                            "el byte posterior al límite confirma el exceso físico",
                        ),
                    ));
                }
                break;
            }
            let permitted = match usize::try_from(remaining.min(chunk.len() as u64)) {
                Ok(permitted) => permitted,
                Err(_) => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::PhysicalInputBytes,
                            limit,
                        ),
                    ));
                }
            };
            let read = reader
                .read(&mut chunk[..permitted])
                .map_err(|error| input_read_error(&normalized_path, error))?;
            if read == 0 {
                break;
            }
            let new_len = match bytes.len().checked_add(read) {
                Some(new_len) => new_len,
                None => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::WorkingMemoryBytes,
                            limit,
                        ),
                    ));
                }
            };
            let new_len_u64 = match u64::try_from(new_len) {
                Ok(new_len) => new_len,
                Err(_) => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::WorkingMemoryBytes,
                            limit,
                        ),
                    ));
                }
            };
            let totals_before_reservation = ledger.as_deref().map(F5ResourceLedger::totals);
            if let Some(ledger) = ledger.as_deref_mut() {
                ledger.charge(
                    role,
                    LedgerCharge::WorkingMemoryPeak(WorkingMemoryPeak::new(new_len_u64)),
                )?;
            }
            if !reserve(&mut bytes, read) {
                if let (Some(ledger), Some(totals)) =
                    (ledger.as_deref_mut(), totals_before_reservation)
                {
                    ledger.restore_totals(totals);
                }
                return Err(resource_diagnostic(
                    phase,
                    role,
                    ResourceDimension::WorkingMemoryBytes,
                    new_len_u64,
                    new_len_u64,
                    "no se pudo reservar el crecimiento prospectivo del búfer",
                ));
            }
            digest.update(&chunk[..read]);
            bytes.extend_from_slice(&chunk[..read]);
            let read_u64 = match u64::try_from(read) {
                Ok(read) => read,
                Err(_) => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::PhysicalInputBytes,
                            limit,
                        ),
                    ));
                }
            };
            observed = match observed.checked_add(read_u64) {
                Some(observed) => observed,
                None => {
                    return Err(reject_resource_error(
                        &mut ledger,
                        overflow_diagnostic(
                            phase,
                            role,
                            ResourceDimension::PhysicalInputBytes,
                            limit,
                        ),
                    ));
                }
            };
        }

        Ok(BoundedInputRead {
            bytes,
            identity: PhysicalInputIdentity {
                role,
                normalized_path,
                observed_size_bytes: observed,
                observed_sha256: format!("{:x}", digest.finalize()),
            },
        })
    }

    fn validate_read_admission(&self, ledger: Option<&F5ResourceLedger>) -> Result<(), Diagnostic> {
        match ledger {
            None if self.phase == AdmissionPhase::Bootstrap && is_bootstrap_role(self.role) => {
                Ok(())
            }
            None => Err(resource_diagnostic(
                self.phase,
                self.role,
                ResourceDimension::OperationUnits,
                0,
                0,
                "una lectura fuera del bootstrap exige un ledger compatible",
            )),
            Some(ledger) if ledger.phase() == AdmissionPhase::Rejected => {
                Err(terminal_rejection_diagnostic())
            }
            Some(ledger)
                if self.phase == ledger.phase() && role_matches_phase(self.role, self.phase) =>
            {
                Ok(())
            }
            Some(_) => Err(resource_diagnostic(
                self.phase,
                self.role,
                ResourceDimension::OperationUnits,
                0,
                0,
                "el rol o la fase de lectura no coincide con el ledger",
            )),
        }
    }
}

fn reject_resource_error(
    ledger: &mut Option<&mut F5ResourceLedger>,
    diagnostic: Diagnostic,
) -> Diagnostic {
    match ledger.as_deref_mut() {
        Some(ledger) => ledger.reject(diagnostic),
        None => diagnostic,
    }
}

fn input_read_error(path: &Path, error: std::io::Error) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-001",
        "No se puede abrir o leer una entrada declarada.",
    )
    .with("path", path.display().to_string())
    .with("detail", error.to_string())
}

#[cfg(test)]
mod tests {
    use std::io::Cursor;
    use std::path::PathBuf;

    use super::{
        AdmissionPhase, BoundedInput, F5ResourceLedger, InputRole, LedgerLimits, OperationCost,
        OutputReservation, WorkingMemoryPeak,
    };

    #[test]
    fn injected_buffer_reservation_failure_restores_ledger_totals() {
        let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(10),
            WorkingMemoryPeak::new(10),
            OutputReservation::new(10),
        ));
        ledger
            .transition_to(AdmissionPhase::ProfileMetadata)
            .expect("profile metadata follows bootstrap");
        let before = ledger.totals();

        let error = BoundedInput::new(
            InputRole::ProfileMetadata,
            AdmissionPhase::ProfileMetadata,
            PathBuf::from("/authoring/profiles/profile.metadata.json"),
            1,
            Cursor::new(b"x".to_vec()),
        )
        .read_to_end_with_reserver(Some(&mut ledger), |_bytes, _additional| false)
        .expect_err("injected allocation failure");

        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(ledger.totals(), before);
        assert_eq!(ledger.phase(), AdmissionPhase::ProfileMetadata);
    }
}
