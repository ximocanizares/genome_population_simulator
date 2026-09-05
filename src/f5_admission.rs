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
    retained_memory_bytes: u64,
    read_operation_escrow_units: u64,
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
            retained_memory_bytes: 0,
            read_operation_escrow_units: 0,
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

    /// Memoria persistente que sigue viva al terminar la última transición
    /// registrada. No es el pico histórico y, por tanto, es la base correcta
    /// para derivar el headroom de la entrada siguiente.
    #[must_use]
    pub const fn retained_memory_bytes(&self) -> u64 {
        self.retained_memory_bytes
    }

    #[must_use]
    pub(crate) const fn remaining_operation_units(&self) -> u64 {
        self.limits
            .operation_units
            .value()
            .saturating_sub(self.totals.operation_units.value())
    }

    /// Capacidad prospectiva disponible para una lectura: headroom nuevo más
    /// una fracción de trabajo futuro ya cargada que puede respaldar la sonda
    /// EOF. El escrow sigue incluido en `totals`; no amplía el presupuesto.
    #[must_use]
    pub(crate) const fn available_read_operation_units(&self) -> u64 {
        self.remaining_operation_units()
            .saturating_add(self.read_operation_escrow_units)
    }

    /// Identifica trabajo futuro ya admitido que puede reclasificarse
    /// temporalmente como lectura. Si la lectura entrega cero bytes, no cambia
    /// ninguna magnitud; si entrega bytes, consume escrow y la fórmula final
    /// deberá admitir el coste adicional antes de continuar.
    pub(crate) fn earmark_read_operation_escrow(
        &mut self,
        role: InputRole,
        units: u64,
    ) -> Result<(), Diagnostic> {
        let prospective = self
            .read_operation_escrow_units
            .checked_add(units)
            .ok_or_else(|| {
                self.reject_resource(
                    role,
                    ResourceDimension::OperationUnits,
                    "overflow del escrow de lectura",
                )
            })?;
        if prospective > self.totals.operation_units.value() {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                role,
                ResourceDimension::OperationUnits,
                prospective,
                self.totals.operation_units.value(),
                "el escrow de lectura no está respaldado por trabajo ya cargado",
            )));
        }
        self.read_operation_escrow_units = prospective;
        Ok(())
    }

    /// Reserva el coste máximo de una lectura física antes de invocar
    /// `Read::read`. El llamador debe liquidarla incluso ante error de E/S.
    pub(crate) fn reserve_read(
        &mut self,
        role: InputRole,
        prospective_bytes: usize,
        units_per_byte: u64,
    ) -> Result<ReadReservation, Diagnostic> {
        let prospective_bytes = u64::try_from(prospective_bytes).map_err(|_| {
            self.reject_resource(
                role,
                ResourceDimension::OperationUnits,
                "el bloque físico no cabe en u64",
            )
        })?;
        let reserved_units = prospective_bytes
            .checked_mul(units_per_byte)
            .ok_or_else(|| {
                self.reject_resource(
                    role,
                    ResourceDimension::OperationUnits,
                    "overflow al reservar el bloque físico",
                )
            })?;
        let before = self.totals.operation_units.value();
        let normal_units = reserved_units.min(self.remaining_operation_units());
        let escrow_units = reserved_units - normal_units;
        if escrow_units > self.read_operation_escrow_units {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                role,
                ResourceDimension::OperationUnits,
                before.saturating_add(reserved_units),
                self.limits.operation_units.value(),
                "reserva prospectiva de lectura excede el límite inclusivo",
            )));
        }
        self.charge(
            role,
            LedgerCharge::Operation(OperationCost::new(normal_units)),
        )?;
        Ok(ReadReservation {
            role,
            before,
            prospective_bytes,
            units_per_byte,
            normal_units,
            escrow_units,
        })
    }

    /// Conserva solo el coste correspondiente a bytes realmente entregados.
    /// La capacidad no usada —incluida la sonda EOF— vuelve al ledger.
    pub(crate) fn settle_read(
        &mut self,
        reservation: ReadReservation,
        observed_bytes: usize,
    ) -> Result<(), Diagnostic> {
        let observed_bytes = u64::try_from(observed_bytes).map_err(|_| {
            self.reject_resource(
                reservation.role,
                ResourceDimension::OperationUnits,
                "los bytes observados no caben en u64",
            )
        })?;
        if observed_bytes > reservation.prospective_bytes {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                reservation.role,
                ResourceDimension::OperationUnits,
                observed_bytes,
                reservation.prospective_bytes,
                "el lector entregó más bytes que el bloque reservado",
            )));
        }
        let observed_units = observed_bytes
            .checked_mul(reservation.units_per_byte)
            .ok_or_else(|| {
                self.reject_resource(
                    reservation.role,
                    ResourceDimension::OperationUnits,
                    "overflow al liquidar el bloque físico",
                )
            })?;
        let normal_units = observed_units.min(reservation.normal_units);
        let escrow_units = observed_units - normal_units;
        if escrow_units > reservation.escrow_units
            || escrow_units > self.read_operation_escrow_units
        {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                reservation.role,
                ResourceDimension::OperationUnits,
                observed_units,
                reservation
                    .normal_units
                    .saturating_add(reservation.escrow_units),
                "la lectura entregó operaciones sin reserva prospectiva",
            )));
        }
        self.read_operation_escrow_units -= escrow_units;
        self.totals.operation_units = OperationCost::new(
            reservation
                .before
                .checked_add(normal_units)
                .ok_or_else(|| {
                    self.reject_resource(
                        reservation.role,
                        ResourceDimension::OperationUnits,
                        "overflow al liquidar el bloque físico",
                    )
                })?,
        );
        Ok(())
    }

    /// Convierte memoria nueva realmente retenida en estado vivo del ledger.
    pub(crate) fn retain_memory(
        &mut self,
        role: InputRole,
        additional_bytes: u64,
        coexisting_transient_bytes: u64,
    ) -> Result<(), Diagnostic> {
        let after = self
            .retained_memory_bytes
            .checked_add(additional_bytes)
            .ok_or_else(|| {
                self.reject_resource(
                    role,
                    ResourceDimension::WorkingMemoryBytes,
                    "overflow de memoria retenida",
                )
            })?;
        self.charge_memory_phase(
            role,
            self.retained_memory_bytes,
            coexisting_transient_bytes,
            after,
        )
    }

    /// Cota física disponible para la siguiente entrada materializada. Se
    /// deriva del presupuesto que queda en las dos dimensiones que una lectura
    /// puede consumir (trabajo y pico de búfer), nunca de `metadata`. Cada
    /// byte de documento necesita dos unidades: lectura/hash y parser, y las
    /// dos representaciones coexisten transitoriamente.
    pub fn remaining_input_bytes(&mut self, role: InputRole) -> Result<u64, Diagnostic> {
        if self.phase == AdmissionPhase::Rejected {
            return Err(terminal_rejection_diagnostic());
        }
        if !role_matches_phase(role, self.phase) {
            return Err(resource_diagnostic(
                self.phase,
                role,
                ResourceDimension::OperationUnits,
                0,
                0,
                "el rol no coincide con la fase para derivar su límite físico",
            ));
        }
        let operation_remaining = self
            .limits
            .operation_units
            .value()
            .saturating_sub(self.totals.operation_units.value());
        let memory_headroom = self
            .limits
            .working_memory_peak
            .value()
            .saturating_sub(self.retained_memory_bytes);
        let operation_input_capacity = operation_remaining / 2;
        let memory_input_capacity = memory_headroom / 2;
        let limit = operation_input_capacity.min(memory_input_capacity);
        if limit == 0 {
            let (dimension, accumulated, maximum) =
                if operation_input_capacity <= memory_input_capacity {
                    (
                        ResourceDimension::OperationUnits,
                        self.totals.operation_units.value(),
                        self.limits.operation_units.value(),
                    )
                } else {
                    (
                        ResourceDimension::WorkingMemoryBytes,
                        self.retained_memory_bytes,
                        self.limits.working_memory_peak.value(),
                    )
                };
            return Err(self.reject(resource_diagnostic(
                self.phase,
                role,
                dimension,
                accumulated,
                maximum,
                "no queda presupuesto para abrir la siguiente entrada",
            )));
        }
        Ok(limit)
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

    /// Registra el pico de una fase sin convertir picos anteriores en una
    /// reserva acumulativa.  Las estructuras persistentes de la fase anterior
    /// y posterior pueden coexistir con su transitorio, pero dos fases ya
    /// terminadas no lo hacen.
    pub fn charge_memory_phase(
        &mut self,
        role: InputRole,
        retained_before_bytes: u64,
        transient_bytes: u64,
        retained_after_bytes: u64,
    ) -> Result<(), Diagnostic> {
        if self.phase == AdmissionPhase::Rejected {
            return Err(terminal_rejection_diagnostic());
        }
        if retained_before_bytes != self.retained_memory_bytes {
            return Err(self.reject(resource_diagnostic(
                self.phase,
                role,
                ResourceDimension::WorkingMemoryBytes,
                retained_before_bytes,
                self.retained_memory_bytes,
                "la retención declarada no coincide con la memoria viva del ledger",
            )));
        }
        let retained = retained_before_bytes.max(retained_after_bytes);
        let peak = retained.checked_add(transient_bytes).ok_or_else(|| {
            self.reject(overflow_diagnostic(
                self.phase,
                role,
                ResourceDimension::WorkingMemoryBytes,
                self.limits.working_memory_peak.value(),
            ))
        })?;
        self.charge(
            role,
            LedgerCharge::WorkingMemoryPeak(WorkingMemoryPeak::new(peak)),
        )?;
        self.retained_memory_bytes = retained_after_bytes;
        Ok(())
    }

    /// Comprueba un transitorio futuro ya conocido contra la retención viva
    /// actual. Actualiza el máximo histórico admitido, pero no convierte la
    /// reserva en retención ni altera la fase.
    pub fn reserve_transient_memory(
        &mut self,
        role: InputRole,
        transient_bytes: u64,
    ) -> Result<(), Diagnostic> {
        self.charge_memory_phase(
            role,
            self.retained_memory_bytes,
            transient_bytes,
            self.retained_memory_bytes,
        )
    }

    fn reject(&mut self, diagnostic: Diagnostic) -> Diagnostic {
        self.phase = AdmissionPhase::Rejected;
        diagnostic
    }

    /// Rechaza explícitamente una aritmética del consumidor que no puede
    /// expresarse como una carga simple, preservando el estado terminal.
    pub fn reject_resource(
        &mut self,
        role: InputRole,
        dimension: ResourceDimension,
        detail: &'static str,
    ) -> Diagnostic {
        let limit = match dimension {
            ResourceDimension::OperationUnits => self.limits.operation_units.value(),
            ResourceDimension::WorkingMemoryBytes | ResourceDimension::PhysicalInputBytes => {
                self.limits.working_memory_peak.value()
            }
            ResourceDimension::OutputBytes => self.limits.output_bytes.value(),
        };
        self.reject(resource_diagnostic(
            self.phase,
            role,
            dimension,
            u64::MAX,
            limit,
            detail,
        ))
    }

    fn restore_totals(&mut self, totals: LedgerTotals) {
        self.totals = totals;
    }
}

#[derive(Clone, Copy, Debug)]
pub(crate) struct ReadReservation {
    role: InputRole,
    before: u64,
    prospective_bytes: u64,
    units_per_byte: u64,
    normal_units: u64,
    escrow_units: u64,
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
    pub(crate) fn from_observed(
        role: InputRole,
        normalized_path: PathBuf,
        observed_size_bytes: u64,
        observed_sha256: String,
    ) -> Self {
        Self {
            role,
            normalized_path,
            observed_size_bytes,
            observed_sha256,
        }
    }
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
        self.read_to_end_with_reserver(ledger, false, None, |bytes, additional| {
            bytes.try_reserve(additional).is_ok()
        })
    }

    /// Igual que la lectura bajo ledger, pero clasifica el primer byte que
    /// excede un tamaño declarado como discrepancia de identidad. Solo se
    /// incorpora al hash observado; nunca se materializa en el búfer.
    pub fn read_to_end_with_ledger_and_expected_identity(
        self,
        ledger: &mut F5ResourceLedger,
        expected_size: u64,
        expected_sha256: &str,
    ) -> Result<BoundedInputRead, Diagnostic> {
        self.read_to_end_with_reserver(
            Some(ledger),
            false,
            Some((expected_size, expected_sha256)),
            |bytes, additional| bytes.try_reserve(additional).is_ok(),
        )
    }

    /// Variante para entradas cuyo coste de operación aún no puede conocerse
    /// antes de EOF (el metadata). Reserva una unidad antes de cada byte; una
    /// sonda EOF se revierte y nunca se contabiliza como contenido.
    pub fn read_to_end_with_ledger_counting_operations(
        self,
        ledger: &mut F5ResourceLedger,
    ) -> Result<BoundedInputRead, Diagnostic> {
        self.read_to_end_with_reserver(Some(ledger), true, None, |bytes, additional| {
            bytes.try_reserve(additional).is_ok()
        })
    }

    fn read_to_end_with_reserver(
        self,
        ledger: Option<&mut F5ResourceLedger>,
        count_operations: bool,
        expected_identity: Option<(u64, &str)>,
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
                    if let Some((expected_size, expected_sha256)) = expected_identity {
                        digest.update(excess);
                        return Err(reject_resource_error(
                            &mut ledger,
                            input_identity_mismatch_diagnostic(
                                &normalized_path,
                                expected_size,
                                expected_sha256,
                                accumulated,
                                &format!("{:x}", digest.finalize()),
                            ),
                        ));
                    }
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
            let permitted = match usize::try_from(if count_operations {
                remaining.min(1)
            } else {
                remaining.min(chunk.len() as u64)
            }) {
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
            let totals_before_operation = ledger.as_deref().map(F5ResourceLedger::totals);
            if count_operations && let Some(ledger) = ledger.as_deref_mut() {
                ledger.charge(role, LedgerCharge::Operation(OperationCost::new(2)))?;
            }
            let read = reader
                .read(&mut chunk[..permitted])
                .map_err(|error| input_read_error(&normalized_path, error))?;
            if read == 0 {
                if count_operations
                    && let (Some(ledger), Some(totals)) =
                        (ledger.as_deref_mut(), totals_before_operation)
                {
                    ledger.restore_totals(totals);
                }
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
                return Err(reject_resource_error(
                    &mut ledger,
                    resource_diagnostic(
                        phase,
                        role,
                        ResourceDimension::WorkingMemoryBytes,
                        new_len_u64,
                        new_len_u64,
                        "no se pudo reservar el crecimiento prospectivo del búfer",
                    ),
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

pub(crate) fn input_identity_mismatch_diagnostic(
    path: &Path,
    expected_size: u64,
    expected_sha256: &str,
    observed_size: u64,
    observed_sha256: &str,
) -> Diagnostic {
    Diagnostic::input(
        "GPS-INPUT-002",
        "Una entrada cambió durante su copia o no coincide con su identidad declarada.",
    )
    .with("path", path.display().to_string())
    .with("expected_size_bytes", expected_size.to_string())
    .with("observed_size_bytes", observed_size.to_string())
    .with("expected_digest", expected_sha256)
    .with("observed_digest", observed_sha256)
}

#[cfg(test)]
mod tests {
    use std::io::{Cursor, Read};
    use std::path::PathBuf;
    use std::sync::{Arc, Mutex};

    use super::{
        AdmissionPhase, BoundedInput, F5ResourceLedger, InputRole, LedgerCharge, LedgerLimits,
        OperationCost, OutputReservation, WorkingMemoryPeak,
    };

    #[test]
    fn input_headroom_uses_live_retention_not_the_historical_peak() {
        let mut ledger = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(1_000),
            WorkingMemoryPeak::new(100),
            OutputReservation::new(1_000),
        ));
        ledger
            .charge(
                InputRole::GeneralToml,
                LedgerCharge::WorkingMemoryPeak(WorkingMemoryPeak::new(100)),
            )
            .expect("bootstrap peak at the inclusive limit");
        ledger
            .transition_to(AdmissionPhase::ProfileMetadata)
            .expect("profile after bootstrap");
        assert_eq!(ledger.retained_memory_bytes(), 0);
        assert_eq!(
            ledger
                .remaining_input_bytes(InputRole::ProfileMetadata)
                .expect("released bootstrap buffers leave headroom"),
            50
        );

        ledger
            .charge_memory_phase(InputRole::ProfileMetadata, 0, 20, 30)
            .expect("metadata leaves typed retained values");
        ledger
            .transition_to(AdmissionPhase::ProfileTableBlockLengths)
            .expect("first table after metadata");
        assert_eq!(ledger.retained_memory_bytes(), 30);
        assert_eq!(
            ledger
                .remaining_input_bytes(InputRole::ProfileBlockLengths)
                .expect("headroom subtracts the live 30 bytes"),
            35
        );
    }

    #[test]
    fn zero_input_limit_reports_the_limiting_divided_dimension_and_is_terminal() {
        let mut operations = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(3),
            WorkingMemoryPeak::new(100),
            OutputReservation::new(100),
        ));
        operations
            .charge(
                InputRole::GeneralToml,
                LedgerCharge::Operation(OperationCost::new(2)),
            )
            .expect("bootstrap operation charge");
        operations
            .transition_to(AdmissionPhase::ProfileMetadata)
            .expect("metadata phase");
        let error = operations
            .remaining_input_bytes(InputRole::ProfileMetadata)
            .expect_err("one remaining operation cannot admit one input byte");
        assert_eq!(error.context()["dimension"], "operation_units");
        assert_eq!(operations.phase(), AdmissionPhase::Rejected);

        let mut memory = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(100),
            WorkingMemoryPeak::new(3),
            OutputReservation::new(100),
        ));
        memory
            .charge_memory_phase(InputRole::GeneralToml, 0, 0, 2)
            .expect("bootstrap retention");
        memory
            .transition_to(AdmissionPhase::ProfileMetadata)
            .expect("metadata phase");
        let error = memory
            .remaining_input_bytes(InputRole::ProfileMetadata)
            .expect_err("one remaining memory byte cannot admit one input byte");
        assert_eq!(error.context()["dimension"], "working_memory_bytes");
        assert_eq!(memory.phase(), AdmissionPhase::Rejected);

        let mut tie = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(3),
            WorkingMemoryPeak::new(3),
            OutputReservation::new(100),
        ));
        tie.charge(
            InputRole::GeneralToml,
            LedgerCharge::Operation(OperationCost::new(2)),
        )
        .expect("bootstrap operation charge");
        tie.charge_memory_phase(InputRole::GeneralToml, 0, 0, 2)
            .expect("bootstrap retention");
        tie.transition_to(AdmissionPhase::ProfileMetadata)
            .expect("metadata phase");
        let error = tie
            .remaining_input_bytes(InputRole::ProfileMetadata)
            .expect_err("a tie is rejected deterministically");
        assert_eq!(error.context()["dimension"], "operation_units");
        assert_eq!(tie.phase(), AdmissionPhase::Rejected);
    }

    #[test]
    fn injected_buffer_reservation_failure_is_atomic_and_terminal() {
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
        .read_to_end_with_reserver(Some(&mut ledger), false, None, |_bytes, _additional| false)
        .expect_err("injected allocation failure");

        assert_eq!(error.code(), "GPS-RES-101");
        assert_eq!(ledger.totals(), before);
        assert_eq!(ledger.phase(), AdmissionPhase::Rejected);
        assert!(
            ledger
                .transition_to(AdmissionPhase::ProfileTableBlockLengths)
                .is_err(),
            "el rechazo de reserva impide transiciones posteriores"
        );

        #[derive(Clone)]
        struct SpyReader {
            reads: Arc<Mutex<u64>>,
        }
        impl Read for SpyReader {
            fn read(&mut self, _buffer: &mut [u8]) -> std::io::Result<usize> {
                *self.reads.lock().expect("spy lock") += 1;
                Ok(0)
            }
        }
        let reads = Arc::new(Mutex::new(0));
        let later = BoundedInput::new(
            InputRole::ProfileBlockLengths,
            AdmissionPhase::ProfileTableBlockLengths,
            PathBuf::from("/authoring/profiles/block_lengths.tsv"),
            1,
            SpyReader {
                reads: Arc::clone(&reads),
            },
        )
        .read_to_end_with_ledger(Some(&mut ledger))
        .expect_err("un ledger rechazado impide otra lectura");
        assert_eq!(later.code(), "GPS-RES-101");
        assert_eq!(*reads.lock().expect("spy lock"), 0, "cero I/O posterior");
    }

    #[test]
    fn physical_read_reservation_precedes_io_and_restores_unused_eof_capacity() {
        #[derive(Clone)]
        struct CountingReader {
            cursor: Cursor<Vec<u8>>,
            bytes_requested: Arc<Mutex<u64>>,
        }
        impl Read for CountingReader {
            fn read(&mut self, buffer: &mut [u8]) -> std::io::Result<usize> {
                *self.bytes_requested.lock().expect("counter lock") += buffer.len() as u64;
                self.cursor.read(buffer)
            }
        }

        let requested = Arc::new(Mutex::new(0));
        let mut exact_reader = CountingReader {
            cursor: Cursor::new(b"abcd".to_vec()),
            bytes_requested: Arc::clone(&requested),
        };
        let mut exact = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(5),
            WorkingMemoryPeak::new(1),
            OutputReservation::new(1),
        ));
        exact
            .charge(
                InputRole::GeneralToml,
                LedgerCharge::Operation(OperationCost::new(1)),
            )
            .expect("one already admitted future operation");
        exact
            .earmark_read_operation_escrow(InputRole::GeneralToml, 1)
            .expect("the future operation backs the EOF probe");
        let mut bytes = [0_u8; 4];
        let reservation = exact
            .reserve_read(InputRole::GeneralToml, bytes.len(), 1)
            .expect("the full physical block is admitted before read");
        let count = exact_reader.read(&mut bytes).expect("physical read");
        exact
            .settle_read(reservation, count)
            .expect("exact delivered block");
        let eof = exact
            .reserve_read(InputRole::GeneralToml, 1, 1)
            .expect("prospective EOF probe");
        let count = exact_reader.read(&mut bytes[..1]).expect("EOF probe");
        exact.settle_read(eof, count).expect("restore unused EOF");
        assert_eq!(count, 0);
        assert_eq!(exact.totals().operation_units().value(), 5);
        assert_eq!(*requested.lock().expect("counter lock"), 5);

        let delivered_reads = Arc::new(Mutex::new(0));
        let mut delivered_reader = CountingReader {
            cursor: Cursor::new(b"abcde".to_vec()),
            bytes_requested: Arc::clone(&delivered_reads),
        };
        let mut delivered = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(5),
            WorkingMemoryPeak::new(1),
            OutputReservation::new(1),
        ));
        delivered
            .charge(
                InputRole::GeneralToml,
                LedgerCharge::Operation(OperationCost::new(1)),
            )
            .expect("one already admitted future operation");
        delivered
            .earmark_read_operation_escrow(InputRole::GeneralToml, 1)
            .expect("one-byte read escrow");
        let first = delivered
            .reserve_read(InputRole::GeneralToml, 4, 1)
            .expect("four regular bytes are admitted");
        let count = delivered_reader.read(&mut bytes).expect("four bytes");
        delivered
            .settle_read(first, count)
            .expect("four bytes charged");
        let borrowed = delivered
            .reserve_read(InputRole::GeneralToml, 1, 1)
            .expect("the last physical byte is backed by charged escrow");
        let count = delivered_reader
            .read(&mut bytes[..1])
            .expect("fifth physical byte");
        delivered
            .settle_read(borrowed, count)
            .expect("the actual byte consumes charged escrow");
        assert_eq!(count, 1);
        assert_eq!(delivered.totals().operation_units().value(), 5);
        assert_eq!(*delivered_reads.lock().expect("counter lock"), 5);
        assert!(
            delivered
                .charge(
                    InputRole::GeneralToml,
                    LedgerCharge::Operation(OperationCost::new(1)),
                )
                .is_err(),
            "the displaced future operation cannot run without formula headroom"
        );

        let rejected_reads = Arc::new(Mutex::new(0));
        let _rejected_reader = CountingReader {
            cursor: Cursor::new(b"abcd".to_vec()),
            bytes_requested: Arc::clone(&rejected_reads),
        };
        let mut minus_one = F5ResourceLedger::new(LedgerLimits::new(
            OperationCost::new(3),
            WorkingMemoryPeak::new(1),
            OutputReservation::new(1),
        ));
        assert!(
            minus_one
                .reserve_read(InputRole::GeneralToml, bytes.len(), 1)
                .is_err()
        );
        assert_eq!(*rejected_reads.lock().expect("counter lock"), 0);
    }
}
