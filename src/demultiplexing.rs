//! Demultiplexado técnico transitorio de F4.5.1.
//!
//! Esta frontera recibe la instantánea posterior a contaminación biológica y
//! mueve cada pareja a lo sumo una vez. No publica FASTQ, manifiestos ni
//! procedencia por pareja; el llamador conserva solamente sus agregados y las
//! parejas finales mientras prepara la siguiente fase.

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostics::Diagnostic;
use crate::identity::{EventId, LibraryId, PreHopInstanceId};
use crate::randomness::{
    F4IndexHopInstance, F4TechnicalDecision, RandomSamplingError, RootSeed,
    derive_f4_technical_subseed, f4_bernoulli, f4_uniform_index,
};
use crate::technical_identity::{F4TechnicalIdentityDeriver, TechnicalPairIdentity};

/// UDI explícito en orientación de secuenciación. La concatenación conserva
/// literalmente `i7+i5`: F4 nunca complementa ni invierte `i5`.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct Udi {
    i7: String,
    i5: String,
}

impl Udi {
    pub fn new(i7: impl Into<String>, i5: impl Into<String>) -> Result<Self, DemultiplexingError> {
        let i7 = i7.into();
        let i5 = i5.into();
        if !is_udi_component(&i7) || !is_udi_component(&i5) {
            return Err(DemultiplexingError::InvalidUdi);
        }
        Ok(Self { i7, i5 })
    }

    #[must_use]
    pub fn i7(&self) -> &str {
        &self.i7
    }

    #[must_use]
    pub fn i5(&self) -> &str {
        &self.i5
    }

    #[must_use]
    pub fn literal(&self) -> String {
        format!("{}+{}", self.i7, self.i5)
    }
}

fn is_udi_component(value: &str) -> bool {
    value.len() == 8
        && value
            .bytes()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

/// Una biblioteca participante con el único UDI completo que puede recibir
/// una pareja tras el salto.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DemultiplexingLibrary {
    library_id: LibraryId,
    udi: Udi,
}

impl DemultiplexingLibrary {
    #[must_use]
    pub const fn new(library_id: LibraryId, udi: Udi) -> Self {
        Self { library_id, udi }
    }

    #[must_use]
    pub const fn library_id(&self) -> &LibraryId {
        &self.library_id
    }

    #[must_use]
    pub const fn udi(&self) -> &Udi {
        &self.udi
    }
}

/// Pareja de la instantánea posterior a contaminación y anterior al salto.
/// `payload` se mueve a la biblioteca final; no hay copia de contenido.
#[derive(Debug)]
pub struct PostContaminationPair<T> {
    pre_hop_library_id: LibraryId,
    pre_hop_instance_id: PreHopInstanceId,
    observed_udi: Udi,
    is_contaminant: bool,
    payload: T,
}

impl<T> PostContaminationPair<T> {
    #[must_use]
    pub const fn new(
        pre_hop_library_id: LibraryId,
        pre_hop_instance_id: PreHopInstanceId,
        observed_udi: Udi,
        is_contaminant: bool,
        payload: T,
    ) -> Self {
        Self {
            pre_hop_library_id,
            pre_hop_instance_id,
            observed_udi,
            is_contaminant,
            payload,
        }
    }
}

/// Pareja ya asignada a su biblioteca final. La identidad `TP_` solo aparece
/// después de que la transferencia simultánea haya fijado esa biblioteca.
#[derive(Debug)]
pub struct DemultiplexedPair<T> {
    final_library_id: LibraryId,
    final_udi: Udi,
    technical_identity: TechnicalPairIdentity,
    is_contaminant: bool,
    payload: T,
}

impl<T> DemultiplexedPair<T> {
    #[must_use]
    pub const fn final_library_id(&self) -> &LibraryId {
        &self.final_library_id
    }

    #[must_use]
    pub const fn final_udi(&self) -> &Udi {
        &self.final_udi
    }

    #[must_use]
    pub const fn technical_identity(&self) -> &TechnicalPairIdentity {
        &self.technical_identity
    }

    #[must_use]
    pub const fn is_contaminant(&self) -> bool {
        self.is_contaminant
    }

    #[must_use]
    pub fn payload(&self) -> &T {
        &self.payload
    }

    #[must_use]
    pub fn into_payload(self) -> T {
        self.payload
    }
}

/// Solo hechos agregados del demultiplexado. `unassigned_pairs` es el
/// subconjunto de descartes cuyo UDI no pertenece a la tanda declarada.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub enum DemultiplexingCause {
    IndexHopping,
    UnrecognizedUdi,
}

impl DemultiplexingCause {
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::IndexHopping => "index_hopping",
            Self::UnrecognizedUdi => "unrecognized_udi",
        }
    }
}

/// Un hecho agregado del demultiplexado. No contiene identificadores de
/// pareja ni UDI observados: solo biblioteca de origen, destino si existe,
/// causa y recuento.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct DemultiplexingAggregate {
    source_library_id: LibraryId,
    destination_library_id: Option<LibraryId>,
    cause: DemultiplexingCause,
    pair_count: u64,
}

impl DemultiplexingAggregate {
    #[must_use]
    pub const fn source_library_id(&self) -> &LibraryId {
        &self.source_library_id
    }

    #[must_use]
    pub const fn destination_library_id(&self) -> Option<&LibraryId> {
        self.destination_library_id.as_ref()
    }

    #[must_use]
    pub const fn cause(&self) -> DemultiplexingCause {
        self.cause
    }

    #[must_use]
    pub const fn pair_count(&self) -> u64 {
        self.pair_count
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct DemultiplexingReport {
    index_hopping_rate: f64,
    reassigned_pairs: u64,
    discarded_pairs: u64,
    unassigned_pairs: u64,
    reassigned_by_sample: Vec<DemultiplexingAggregate>,
    discarded_by_sample: Vec<DemultiplexingAggregate>,
    unassigned_by_sample: Vec<DemultiplexingAggregate>,
}

impl DemultiplexingReport {
    #[must_use]
    pub const fn index_hopping_rate(&self) -> f64 {
        self.index_hopping_rate
    }

    #[must_use]
    pub const fn reassigned_pairs(&self) -> u64 {
        self.reassigned_pairs
    }

    #[must_use]
    pub const fn discarded_pairs(&self) -> u64 {
        self.discarded_pairs
    }

    #[must_use]
    pub const fn unassigned_pairs(&self) -> u64 {
        self.unassigned_pairs
    }

    #[must_use]
    pub fn reassigned_by_sample(&self) -> &[DemultiplexingAggregate] {
        &self.reassigned_by_sample
    }

    #[must_use]
    pub fn discarded_by_sample(&self) -> &[DemultiplexingAggregate] {
        &self.discarded_by_sample
    }

    #[must_use]
    pub fn unassigned_by_sample(&self) -> &[DemultiplexingAggregate] {
        &self.unassigned_by_sample
    }
}

/// Resultado efímero agrupado por biblioteca final y el resumen que podrá
/// consumir F4.5.2. No escribe ni conserva un FASTQ limpio adicional.
#[derive(Debug)]
pub struct DemultiplexingOutcome<T> {
    final_libraries: BTreeMap<LibraryId, Vec<DemultiplexedPair<T>>>,
    report: DemultiplexingReport,
}

impl<T> DemultiplexingOutcome<T> {
    #[must_use]
    pub fn final_libraries(&self) -> &BTreeMap<LibraryId, Vec<DemultiplexedPair<T>>> {
        &self.final_libraries
    }

    #[must_use]
    pub const fn report(&self) -> &DemultiplexingReport {
        &self.report
    }

    /// Entrega las parejas ya clasificadas al único consumidor que publicará
    /// FASTQ. Esta transferencia evita clonar su contenido o conservar una
    /// instantánea de procedencia después de la salida final.
    #[must_use]
    pub fn into_final_libraries(self) -> BTreeMap<LibraryId, Vec<DemultiplexedPair<T>>> {
        self.final_libraries
    }
}

#[derive(Debug)]
pub enum DemultiplexingError {
    InvalidUdi,
    InvalidIndexHoppingRate,
    DuplicateLibrary,
    DuplicateI7,
    DuplicateI5,
    MissingAlternativeLibrary,
    UnknownPreHopLibrary,
    RandomSampling(RandomSamplingError),
    TechnicalIdentity(Diagnostic),
    CounterOverflow,
}

/// Consumidor de la instantánea post-contaminación. El registro de
/// identidades vive durante toda la tanda para detectar colisiones `TP_`.
#[derive(Debug)]
pub struct Demultiplexer {
    root_seed: RootSeed,
    event_id: EventId,
    index_hopping_rate: f64,
    libraries: BTreeMap<LibraryId, Udi>,
    udi_libraries: BTreeMap<Udi, LibraryId>,
    identity_deriver: F4TechnicalIdentityDeriver,
}

impl Demultiplexer {
    pub fn new(
        root_seed: RootSeed,
        event_id: EventId,
        index_hopping_rate: f64,
        libraries: Vec<DemultiplexingLibrary>,
    ) -> Result<Self, DemultiplexingError> {
        if !index_hopping_rate.is_finite() || !(0.0..=1.0).contains(&index_hopping_rate) {
            return Err(DemultiplexingError::InvalidIndexHoppingRate);
        }
        if index_hopping_rate > 0.0 && libraries.len() < 2 {
            return Err(DemultiplexingError::MissingAlternativeLibrary);
        }
        let mut by_library = BTreeMap::new();
        let mut by_udi = BTreeMap::new();
        let mut i7_values = BTreeSet::new();
        let mut i5_values = BTreeSet::new();
        for library in libraries {
            if !i7_values.insert(library.udi.i7.clone()) {
                return Err(DemultiplexingError::DuplicateI7);
            }
            if !i5_values.insert(library.udi.i5.clone()) {
                return Err(DemultiplexingError::DuplicateI5);
            }
            if by_library
                .insert(library.library_id.clone(), library.udi.clone())
                .is_some()
            {
                return Err(DemultiplexingError::DuplicateLibrary);
            }
            if by_udi.insert(library.udi, library.library_id).is_some() {
                return Err(DemultiplexingError::InvalidUdi);
            }
        }
        Ok(Self {
            root_seed,
            event_id: event_id.clone(),
            index_hopping_rate,
            libraries: by_library,
            udi_libraries: by_udi,
            identity_deriver: F4TechnicalIdentityDeriver::new(root_seed, event_id),
        })
    }

    /// Evalúa cada entrada una vez contra la misma instantánea post-biológica.
    /// Cada subsemilla depende de la biblioteca e instancia anteriores al
    /// salto, de manera que el orden de entrada no cambia ningún sorteo.
    pub fn demultiplex<T>(
        &mut self,
        pairs: Vec<PostContaminationPair<T>>,
    ) -> Result<DemultiplexingOutcome<T>, DemultiplexingError> {
        let mut final_libraries = self
            .libraries
            .keys()
            .cloned()
            .map(|library_id| (library_id, Vec::new()))
            .collect::<BTreeMap<_, _>>();
        let mut report = DemultiplexingReport {
            index_hopping_rate: self.index_hopping_rate,
            reassigned_pairs: 0,
            discarded_pairs: 0,
            unassigned_pairs: 0,
            reassigned_by_sample: Vec::new(),
            discarded_by_sample: Vec::new(),
            unassigned_by_sample: Vec::new(),
        };

        for pair in pairs {
            let pre_hop_udi = self
                .libraries
                .get(&pair.pre_hop_library_id)
                .ok_or(DemultiplexingError::UnknownPreHopLibrary)?;
            let instance = F4IndexHopInstance::parse(pair.pre_hop_instance_id.as_str())
                .map_err(DemultiplexingError::TechnicalIdentity)?;
            let Some(mut final_library_id) = self.udi_libraries.get(&pair.observed_udi).cloned()
            else {
                report.discarded_pairs = checked_increment(report.discarded_pairs)?;
                report.unassigned_pairs = checked_increment(report.unassigned_pairs)?;
                record_aggregate(
                    &mut report.discarded_by_sample,
                    &pair.pre_hop_library_id,
                    None,
                    DemultiplexingCause::UnrecognizedUdi,
                )?;
                record_aggregate(
                    &mut report.unassigned_by_sample,
                    &pair.pre_hop_library_id,
                    None,
                    DemultiplexingCause::UnrecognizedUdi,
                )?;
                continue;
            };
            if &pair.observed_udi != pre_hop_udi {
                return Err(DemultiplexingError::InvalidUdi);
            }
            if self.index_hopping_rate > 0.0 {
                let decision = F4TechnicalDecision::index_hop_bernoulli(
                    &self.event_id,
                    &pair.pre_hop_library_id,
                    instance,
                );
                let mut rng = derive_f4_technical_subseed(self.root_seed, decision)
                    .map_err(DemultiplexingError::TechnicalIdentity)?
                    .local_rng();
                if f4_bernoulli(&mut rng, self.index_hopping_rate)
                    .map_err(DemultiplexingError::RandomSampling)?
                {
                    let destination = F4TechnicalDecision::index_hop_destination(
                        &self.event_id,
                        &pair.pre_hop_library_id,
                        instance,
                    );
                    let mut destination_rng =
                        derive_f4_technical_subseed(self.root_seed, destination)
                            .map_err(DemultiplexingError::TechnicalIdentity)?
                            .local_rng();
                    let destination_index = usize::try_from(
                        f4_uniform_index(
                            &mut destination_rng,
                            u64::try_from(self.libraries.len() - 1)
                                .map_err(|_| DemultiplexingError::CounterOverflow)?,
                        )
                        .map_err(DemultiplexingError::RandomSampling)?,
                    )
                    .map_err(|_| DemultiplexingError::CounterOverflow)?;
                    final_library_id = self
                        .libraries
                        .keys()
                        .filter(|library_id| *library_id != &pair.pre_hop_library_id)
                        .nth(destination_index)
                        .cloned()
                        .ok_or(DemultiplexingError::MissingAlternativeLibrary)?;
                    report.reassigned_pairs = checked_increment(report.reassigned_pairs)?;
                    record_aggregate(
                        &mut report.reassigned_by_sample,
                        &pair.pre_hop_library_id,
                        Some(&final_library_id),
                        DemultiplexingCause::IndexHopping,
                    )?;
                }
            }
            let final_udi = self
                .libraries
                .get(&final_library_id)
                .cloned()
                .ok_or(DemultiplexingError::UnknownPreHopLibrary)?;
            let technical_pair_id = self
                .identity_deriver
                .derive_technical_pair(&final_library_id, &pair.pre_hop_instance_id)
                .map_err(DemultiplexingError::TechnicalIdentity)?;
            final_libraries
                .get_mut(&final_library_id)
                .ok_or(DemultiplexingError::UnknownPreHopLibrary)?
                .push(DemultiplexedPair {
                    final_library_id,
                    final_udi,
                    technical_identity: TechnicalPairIdentity::new(technical_pair_id),
                    is_contaminant: pair.is_contaminant,
                    payload: pair.payload,
                });
        }
        Ok(DemultiplexingOutcome {
            final_libraries,
            report,
        })
    }
}

fn record_aggregate(
    aggregates: &mut Vec<DemultiplexingAggregate>,
    source_library_id: &LibraryId,
    destination_library_id: Option<&LibraryId>,
    cause: DemultiplexingCause,
) -> Result<(), DemultiplexingError> {
    if let Some(existing) = aggregates.iter_mut().find(|entry| {
        entry.source_library_id == *source_library_id
            && entry.destination_library_id.as_ref() == destination_library_id
            && entry.cause == cause
    }) {
        existing.pair_count = checked_increment(existing.pair_count)?;
    } else {
        aggregates.push(DemultiplexingAggregate {
            source_library_id: source_library_id.clone(),
            destination_library_id: destination_library_id.cloned(),
            cause,
            pair_count: 1,
        });
        aggregates.sort_by(|left, right| {
            (
                left.source_library_id.as_str(),
                left.cause,
                left.destination_library_id.as_ref().map(LibraryId::as_str),
            )
                .cmp(&(
                    right.source_library_id.as_str(),
                    right.cause,
                    right.destination_library_id.as_ref().map(LibraryId::as_str),
                ))
        });
    }
    Ok(())
}

fn checked_increment(value: u64) -> Result<u64, DemultiplexingError> {
    value
        .checked_add(1)
        .ok_or(DemultiplexingError::CounterOverflow)
}
