//! Planificación y materialización efímera de la contaminación biológica F4.4.
//!
//! El plan recibe solo claves de parejas limpias. La materialización posterior
//! regenera cada fuente bajo demanda, la entrega al emisor y la descarta: no
//! conserva bases, calidades, FASTQ ni procedencia individual durable.

use std::collections::BTreeMap;
use std::fmt;

use crate::diagnostics::Diagnostic;
use crate::identity::{
    CompactPreHopInstanceId, EventId, IndividualId, LibraryId, PopulationId, PreHopInstanceId,
};
use crate::randomness::{
    F4TechnicalDecision, RootSeed, Subseed, derive_f4_technical_subseed, f4_partial_fisher_yates,
    f4_uniform_closed,
};
use crate::technical_identity::F4TechnicalIdentityDeriver;

/// Intervalo inclusivo de porcentajes y receptoras opcionales de F4.4.
#[derive(Clone, Debug, PartialEq)]
pub struct BiologicalContaminationSettings {
    minimum_percent: f64,
    maximum_percent: f64,
    configured_recipients: Option<Vec<IndividualId>>,
}

impl BiologicalContaminationSettings {
    pub fn new(
        minimum_percent: f64,
        maximum_percent: f64,
        configured_recipients: Option<Vec<IndividualId>>,
    ) -> Result<Self, BiologicalContaminationError> {
        if !minimum_percent.is_finite()
            || !maximum_percent.is_finite()
            || !(0.0..=100.0).contains(&minimum_percent)
            || !(0.0..=100.0).contains(&maximum_percent)
            || minimum_percent > maximum_percent
        {
            return Err(BiologicalContaminationError::InvalidPercentRange);
        }
        if let Some(recipients) = &configured_recipients
            && (recipients.is_empty() || recipients.windows(2).any(|pair| pair[0] >= pair[1]))
        {
            return Err(BiologicalContaminationError::InvalidRecipientSubset);
        }
        Ok(Self {
            minimum_percent,
            maximum_percent,
            configured_recipients,
        })
    }

    #[must_use]
    pub fn disabled() -> Self {
        Self {
            minimum_percent: 0.0,
            maximum_percent: 0.0,
            configured_recipients: None,
        }
    }

    #[must_use]
    pub const fn minimum_percent(&self) -> f64 {
        self.minimum_percent
    }

    #[must_use]
    pub const fn maximum_percent(&self) -> f64 {
        self.maximum_percent
    }

    /// `None` significa toda la población seleccionada; una lista explícita ya
    /// debe venir ordenada y sin repetidos desde la frontera de configuración.
    #[must_use]
    pub fn configured_recipients(&self) -> Option<&[IndividualId]> {
        self.configured_recipients.as_deref()
    }
}

/// Clave transitoria de una pareja propia anterior a toda contaminación.
#[derive(Clone, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct CleanPairKey {
    source_individual_id: IndividualId,
    /// La instancia se guarda por sus símbolos, sin asignación: el inventario
    /// retiene una clave por pareja de todo el evento y su texto era el mayor
    /// consumidor que crecía con la cobertura. El orden es el mismo, porque
    /// comparar los símbolos equivale a comparar el texto.
    pre_hop_instance_id: CompactPreHopInstanceId,
}

impl CleanPairKey {
    pub fn new(
        source_individual_id: IndividualId,
        pre_hop_instance_id: &PreHopInstanceId,
    ) -> Result<Self, Diagnostic> {
        Ok(Self {
            source_individual_id,
            pre_hop_instance_id: CompactPreHopInstanceId::new(pre_hop_instance_id)?,
        })
    }

    #[must_use]
    pub const fn source_individual_id(&self) -> &IndividualId {
        &self.source_individual_id
    }

    #[must_use]
    pub fn pre_hop_instance_id(&self) -> PreHopInstanceId {
        self.pre_hop_instance_id.to_instance_id()
    }

    /// La forma compacta, para quien solo necesita comparar o indexar.
    #[must_use]
    pub const fn compact_pre_hop_instance_id(&self) -> CompactPreHopInstanceId {
        self.pre_hop_instance_id
    }
}

/// Inventario de parejas propias de una única biblioteca antes de F4.4.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct CleanLibraryPairs {
    individual_id: IndividualId,
    population_id: PopulationId,
    library_id: LibraryId,
    own_pairs: Vec<CleanPairKey>,
}

impl CleanLibraryPairs {
    pub fn new(
        individual_id: IndividualId,
        population_id: PopulationId,
        library_id: LibraryId,
        mut own_pairs: Vec<CleanPairKey>,
    ) -> Result<Self, BiologicalContaminationError> {
        if own_pairs
            .iter()
            .any(|pair| pair.source_individual_id() != &individual_id)
        {
            return Err(BiologicalContaminationError::ForeignOwnPair);
        }
        own_pairs.sort();
        if own_pairs.windows(2).any(|pair| pair[0] == pair[1]) {
            return Err(BiologicalContaminationError::DuplicateCleanPair);
        }
        Ok(Self {
            individual_id,
            population_id,
            library_id,
            own_pairs,
        })
    }

    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn population_id(&self) -> &PopulationId {
        &self.population_id
    }

    #[must_use]
    pub const fn library_id(&self) -> &LibraryId {
        &self.library_id
    }

    #[must_use]
    pub fn own_pairs(&self) -> &[CleanPairKey] {
        &self.own_pairs
    }
}

/// Una sustitución futura: copia una pareja propia de otra persona y retira una
/// pareja propia de la receptora. Ninguna de las dos se materializa aquí.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BiologicalContaminationAssignment {
    source_pair: CleanPairKey,
    removed_recipient_pair: CleanPairKey,
}

impl BiologicalContaminationAssignment {
    #[must_use]
    pub const fn source_pair(&self) -> &CleanPairKey {
        &self.source_pair
    }

    #[must_use]
    pub const fn removed_recipient_pair(&self) -> &CleanPairKey {
        &self.removed_recipient_pair
    }
}

/// Hechos temporales por receptora que F4.4.2 consume para materializar cada
/// copia. No puede convertirse directamente en un registro de manifiesto:
/// este solo existe después de una materialización satisfactoria.
#[derive(Clone, Debug, PartialEq)]
pub struct RecipientBiologicalContaminationPlan {
    recipient_id: IndividualId,
    recipient_library_id: LibraryId,
    minimum_percent: f64,
    maximum_percent: f64,
    sampled_percent: f64,
    clean_pair_count: u64,
    target_pair_count: u64,
    percent_subseed: Subseed,
    donor_subseed: Subseed,
    removal_subseed: Subseed,
    donor_contributions: BTreeMap<IndividualId, u64>,
    assignments: Vec<BiologicalContaminationAssignment>,
}

impl RecipientBiologicalContaminationPlan {
    #[must_use]
    pub const fn recipient_id(&self) -> &IndividualId {
        &self.recipient_id
    }
    #[must_use]
    pub const fn recipient_library_id(&self) -> &LibraryId {
        &self.recipient_library_id
    }
    #[must_use]
    pub const fn minimum_percent(&self) -> f64 {
        self.minimum_percent
    }
    #[must_use]
    pub const fn maximum_percent(&self) -> f64 {
        self.maximum_percent
    }
    #[must_use]
    pub const fn sampled_percent(&self) -> f64 {
        self.sampled_percent
    }
    #[must_use]
    pub const fn clean_pair_count(&self) -> u64 {
        self.clean_pair_count
    }
    #[must_use]
    pub const fn target_pair_count(&self) -> u64 {
        self.target_pair_count
    }
    #[must_use]
    pub const fn percent_subseed(&self) -> Subseed {
        self.percent_subseed
    }
    #[must_use]
    pub const fn donor_subseed(&self) -> Subseed {
        self.donor_subseed
    }
    #[must_use]
    pub const fn removal_subseed(&self) -> Subseed {
        self.removal_subseed
    }
    #[must_use]
    pub const fn donor_contributions(&self) -> &BTreeMap<IndividualId, u64> {
        &self.donor_contributions
    }
    #[must_use]
    pub fn assignments(&self) -> &[BiologicalContaminationAssignment] {
        &self.assignments
    }

    /// Une los hechos del plan con los únicos contadores que prueba una
    /// emisión satisfactoria. El informe no se puede construir fuera del
    /// materializador, por lo que esta es la única frontera que crea el
    /// resumen publicable F4.4.2.
    #[must_use]
    pub fn realized_manifest_recipient(
        &self,
        report: &BiologicalContaminationMaterializationReport,
    ) -> Option<BiologicalContaminationManifestRecipient> {
        if report.recipient_id() != &self.recipient_id
            || report.recipient_library_id() != &self.recipient_library_id
        {
            return None;
        }
        let final_pair_count = report.final_pair_count();
        (report.removed_own_pairs() == self.target_pair_count
            && report.inserted_contaminant_pairs() == self.target_pair_count
            && report
                .retained_own_pairs()
                .checked_add(report.removed_own_pairs())
                == Some(self.clean_pair_count)
            && report
                .retained_own_pairs()
                .checked_add(report.inserted_contaminant_pairs())
                == Some(final_pair_count)
            && final_pair_count == self.clean_pair_count)
            .then(|| BiologicalContaminationManifestRecipient {
                recipient_id: self.recipient_id.clone(),
                recipient_library_id: self.recipient_library_id.clone(),
                minimum_percent: self.minimum_percent,
                maximum_percent: self.maximum_percent,
                sampled_percent: self.sampled_percent,
                clean_pair_count: self.clean_pair_count,
                target_pair_count: self.target_pair_count,
                removed_own_pair_count: report.removed_own_pairs(),
                inserted_contaminant_pair_count: report.inserted_contaminant_pairs(),
                retained_own_pair_count: report.retained_own_pairs(),
                final_pair_count,
                percent_subseed: self.percent_subseed,
                donor_subseed: self.donor_subseed,
                removal_subseed: self.removal_subseed,
                donor_contributions: self.donor_contributions.clone(),
            })
    }
}

#[derive(Clone, Debug, PartialEq)]
pub struct BiologicalContaminationPlan {
    recipients: Vec<RecipientBiologicalContaminationPlan>,
}

impl BiologicalContaminationPlan {
    #[must_use]
    pub fn recipients(&self) -> &[RecipientBiologicalContaminationPlan] {
        &self.recipients
    }
}

/// Pareja entregada inmediatamente al escritor final de una biblioteca.
///
/// `payload` se mueve desde la pareja propia o se regenera bajo demanda desde
/// una clave propia del donante. Esta vista no conserva el ID de origen:
/// `TechnicalPairId` se deriva después, en F4.5, desde la biblioteca final y
/// esta instancia `PH_` nueva.
#[derive(Debug)]
pub struct MaterializedBiologicalPair<T> {
    recipient_library_id: LibraryId,
    pre_hop_instance_id: PreHopInstanceId,
    is_contaminant: bool,
    payload: T,
}

impl<T> MaterializedBiologicalPair<T> {
    #[must_use]
    pub const fn recipient_library_id(&self) -> &LibraryId {
        &self.recipient_library_id
    }

    #[must_use]
    pub const fn pre_hop_instance_id(&self) -> &PreHopInstanceId {
        &self.pre_hop_instance_id
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

/// Recuentos de una materialización consumida, ligados a su receptora. No
/// contiene una clave, un ID de origen ni un enlace individual que pueda acabar
/// en el manifiesto.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BiologicalContaminationMaterializationReport {
    recipient_id: IndividualId,
    recipient_library_id: LibraryId,
    retained_own_pairs: u64,
    removed_own_pairs: u64,
    inserted_contaminant_pairs: u64,
}

impl BiologicalContaminationMaterializationReport {
    #[must_use]
    pub const fn recipient_id(&self) -> &IndividualId {
        &self.recipient_id
    }

    #[must_use]
    pub const fn recipient_library_id(&self) -> &LibraryId {
        &self.recipient_library_id
    }

    #[must_use]
    pub const fn retained_own_pairs(&self) -> u64 {
        self.retained_own_pairs
    }

    #[must_use]
    pub const fn removed_own_pairs(&self) -> u64 {
        self.removed_own_pairs
    }

    #[must_use]
    pub const fn inserted_contaminant_pairs(&self) -> u64 {
        self.inserted_contaminant_pairs
    }

    #[must_use]
    pub const fn final_pair_count(&self) -> u64 {
        self.retained_own_pairs + self.inserted_contaminant_pairs
    }
}

/// Error de la frontera de materialización. `Source` y `Emit` pertenecen al
/// productor regenerable y al escritor temporal del llamador, respectivamente.
#[derive(Debug)]
pub enum BiologicalContaminationMaterializationError<E> {
    RecipientInventoryMismatch,
    InvalidPlan,
    TechnicalIdentity(Diagnostic),
    Source(E),
    Emit(E),
    CounterOverflow,
}

/// Estado efímero que comparte el registro de colisiones durante F4.4.
///
/// Se inicializa con todas las parejas propias del evento, de modo que las
/// copias contaminantes no pueden colisionar con una instancia ya existente.
/// No recibe secuencias ni calidades y no retiene parejas entre emisiones.
#[derive(Debug)]
pub struct BiologicalContaminationMaterializer {
    identity_deriver: F4TechnicalIdentityDeriver,
}

impl BiologicalContaminationMaterializer {
    pub fn new(
        root_seed: RootSeed,
        event_id: EventId,
        libraries: &[CleanLibraryPairs],
    ) -> Result<Self, BiologicalContaminationMaterializationError<()>> {
        let mut identity_deriver = F4TechnicalIdentityDeriver::new(root_seed, event_id);
        for pair in libraries.iter().flat_map(|library| library.own_pairs()) {
            identity_deriver
                .register_existing_pre_hop_instance(&pair.pre_hop_instance_id())
                .map_err(BiologicalContaminationMaterializationError::TechnicalIdentity)?;
        }
        Ok(Self { identity_deriver })
    }

    /// Emite una receptora completamente sustituida sin acumular sus lecturas.
    ///
    /// Cada carga se regenera justo antes de emitirse. Las parejas propias no
    /// retiradas y cada contaminante se consumen de una en una; nunca se forma
    /// un mapa de cargas por pareja. Por tanto la memoria de cargas es O(1), no
    /// se clonan secuencias ni calidades, y no queda un vínculo por lectura al
    /// acabar.
    pub fn emit_materialized_recipient<T, E>(
        &mut self,
        recipient_plan: &RecipientBiologicalContaminationPlan,
        clean_recipient: &CleanLibraryPairs,
        mut regenerate_clean: impl FnMut(&CleanPairKey) -> Result<T, E>,
        mut emit: impl FnMut(MaterializedBiologicalPair<T>) -> Result<(), E>,
    ) -> Result<
        BiologicalContaminationMaterializationReport,
        BiologicalContaminationMaterializationError<E>,
    > {
        if recipient_plan.recipient_id() != clean_recipient.individual_id()
            || recipient_plan.recipient_library_id() != clean_recipient.library_id()
        {
            return Err(BiologicalContaminationMaterializationError::RecipientInventoryMismatch);
        }
        let removed = recipient_plan
            .assignments()
            .iter()
            .map(BiologicalContaminationAssignment::removed_recipient_pair)
            .collect::<std::collections::BTreeSet<_>>();
        if removed.len() != recipient_plan.assignments().len()
            || removed.len()
                != usize::try_from(recipient_plan.target_pair_count())
                    .map_err(|_| BiologicalContaminationMaterializationError::CounterOverflow)?
            || !removed.iter().all(|key| {
                key.source_individual_id() == clean_recipient.individual_id()
                    && clean_recipient.own_pairs().binary_search(key).is_ok()
            })
        {
            return Err(BiologicalContaminationMaterializationError::InvalidPlan);
        }

        let mut retained_own_pairs = 0_u64;
        for key in clean_recipient.own_pairs() {
            if removed.contains(key) {
                continue;
            }
            let payload = regenerate_clean(key)
                .map_err(BiologicalContaminationMaterializationError::Source)?;
            emit(MaterializedBiologicalPair {
                recipient_library_id: clean_recipient.library_id().clone(),
                pre_hop_instance_id: key.pre_hop_instance_id().clone(),
                is_contaminant: false,
                payload,
            })
            .map_err(BiologicalContaminationMaterializationError::Emit)?;
            retained_own_pairs = retained_own_pairs
                .checked_add(1)
                .ok_or(BiologicalContaminationMaterializationError::CounterOverflow)?;
        }

        let mut inserted_contaminant_pairs = 0_u64;
        for assignment in recipient_plan.assignments() {
            let payload = regenerate_clean(assignment.source_pair())
                .map_err(BiologicalContaminationMaterializationError::Source)?;
            let pre_hop_instance_id = self
                .identity_deriver
                .derive_contamination_pre_hop_instance(
                    clean_recipient.library_id(),
                    &assignment.source_pair().pre_hop_instance_id(),
                )
                .map_err(BiologicalContaminationMaterializationError::TechnicalIdentity)?;
            emit(MaterializedBiologicalPair {
                recipient_library_id: clean_recipient.library_id().clone(),
                pre_hop_instance_id,
                is_contaminant: true,
                payload,
            })
            .map_err(BiologicalContaminationMaterializationError::Emit)?;
            inserted_contaminant_pairs = inserted_contaminant_pairs
                .checked_add(1)
                .ok_or(BiologicalContaminationMaterializationError::CounterOverflow)?;
        }
        let removed_own_pairs = u64::try_from(removed.len())
            .map_err(|_| BiologicalContaminationMaterializationError::CounterOverflow)?;
        if retained_own_pairs.checked_add(removed_own_pairs)
            != Some(recipient_plan.clean_pair_count())
            || removed_own_pairs != inserted_contaminant_pairs
        {
            return Err(BiologicalContaminationMaterializationError::InvalidPlan);
        }
        Ok(BiologicalContaminationMaterializationReport {
            recipient_id: clean_recipient.individual_id().clone(),
            recipient_library_id: clean_recipient.library_id().clone(),
            retained_own_pairs,
            removed_own_pairs,
            inserted_contaminant_pairs,
        })
    }
}

/// Registro agregado por receptora que F4.4.2 entrega al manifiesto de evento.
/// No tiene claves de pareja, bases ni calidades.
#[derive(Clone, Debug, PartialEq)]
pub struct BiologicalContaminationManifestRecipient {
    recipient_id: IndividualId,
    recipient_library_id: LibraryId,
    minimum_percent: f64,
    maximum_percent: f64,
    sampled_percent: f64,
    clean_pair_count: u64,
    target_pair_count: u64,
    removed_own_pair_count: u64,
    inserted_contaminant_pair_count: u64,
    retained_own_pair_count: u64,
    final_pair_count: u64,
    percent_subseed: Subseed,
    donor_subseed: Subseed,
    removal_subseed: Subseed,
    donor_contributions: BTreeMap<IndividualId, u64>,
}

impl BiologicalContaminationManifestRecipient {
    #[must_use]
    pub const fn recipient_id(&self) -> &IndividualId {
        &self.recipient_id
    }
    #[must_use]
    pub const fn recipient_library_id(&self) -> &LibraryId {
        &self.recipient_library_id
    }
    #[must_use]
    pub const fn minimum_percent(&self) -> f64 {
        self.minimum_percent
    }
    #[must_use]
    pub const fn maximum_percent(&self) -> f64 {
        self.maximum_percent
    }
    #[must_use]
    pub const fn sampled_percent(&self) -> f64 {
        self.sampled_percent
    }
    #[must_use]
    pub const fn clean_pair_count(&self) -> u64 {
        self.clean_pair_count
    }
    #[must_use]
    pub const fn target_pair_count(&self) -> u64 {
        self.target_pair_count
    }
    #[must_use]
    pub const fn removed_own_pair_count(&self) -> u64 {
        self.removed_own_pair_count
    }
    #[must_use]
    pub const fn inserted_contaminant_pair_count(&self) -> u64 {
        self.inserted_contaminant_pair_count
    }
    #[must_use]
    pub const fn retained_own_pair_count(&self) -> u64 {
        self.retained_own_pair_count
    }
    #[must_use]
    pub const fn final_pair_count(&self) -> u64 {
        self.final_pair_count
    }
    #[must_use]
    pub const fn percent_subseed(&self) -> Subseed {
        self.percent_subseed
    }
    #[must_use]
    pub const fn donor_subseed(&self) -> Subseed {
        self.donor_subseed
    }
    #[must_use]
    pub const fn removal_subseed(&self) -> Subseed {
        self.removal_subseed
    }
    #[must_use]
    pub const fn donor_contributions(&self) -> &BTreeMap<IndividualId, u64> {
        &self.donor_contributions
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum BiologicalContaminationError {
    InvalidPercentRange,
    InvalidRecipientSubset,
    DuplicateLibrary,
    ForeignOwnPair,
    DuplicateCleanPair,
    UnknownRecipient(IndividualId),
    MissingDonor {
        recipient: IndividualId,
    },
    InsufficientDonorCapacity {
        recipient: IndividualId,
        requested: u64,
        available: u64,
    },
    CountDoesNotFitPlatform,
    Randomness,
}

impl fmt::Display for BiologicalContaminationError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidPercentRange => write!(
                formatter,
                "el rango de porcentaje debe ser finito, estar en 0..=100 y ser creciente"
            ),
            Self::InvalidRecipientSubset => write!(
                formatter,
                "recipient_ids debe ser una lista no vacía, ordenada y sin repetidos"
            ),
            Self::DuplicateLibrary => write!(
                formatter,
                "cada individuo solo puede aportar una biblioteca limpia"
            ),
            Self::ForeignOwnPair => {
                write!(formatter, "una pareja limpia propia tiene un origen ajeno")
            }
            Self::DuplicateCleanPair => write!(
                formatter,
                "una biblioteca repite una clave de pareja limpia"
            ),
            Self::UnknownRecipient(recipient) => {
                write!(formatter, "receptora desconocida: {recipient}")
            }
            Self::MissingDonor { recipient } => write!(
                formatter,
                "la receptora {recipient} no tiene donantes elegibles"
            ),
            Self::InsufficientDonorCapacity {
                recipient,
                requested,
                available,
            } => write!(
                formatter,
                "la receptora {recipient} solicita {requested} parejas y solo dispone de {available} donantes"
            ),
            Self::CountDoesNotFitPlatform => {
                write!(formatter, "el número de parejas no cabe en esta plataforma")
            }
            Self::Randomness => write!(formatter, "falló el sorteo reproducible de contaminación"),
        }
    }
}

impl std::error::Error for BiologicalContaminationError {}

/// Construye todo el plan sobre instantáneas propias limpias. El resultado es
/// estable frente al orden de `libraries` y no usa resultados de otra
/// receptora como candidatos.
pub fn plan_biological_contamination(
    root_seed: RootSeed,
    event_id: &EventId,
    settings: &BiologicalContaminationSettings,
    libraries: &[CleanLibraryPairs],
) -> Result<BiologicalContaminationPlan, BiologicalContaminationError> {
    let by_individual = libraries
        .iter()
        .try_fold(BTreeMap::new(), |mut values, library| {
            if values.insert(library.individual_id(), library).is_some() {
                Err(BiologicalContaminationError::DuplicateLibrary)
            } else {
                Ok(values)
            }
        })?;
    let recipient_ids = match settings.configured_recipients() {
        Some(recipients) => recipients.to_vec(),
        None => by_individual.keys().map(|id| (*id).clone()).collect(),
    };
    let mut recipients = Vec::with_capacity(recipient_ids.len());
    for recipient_id in recipient_ids {
        let recipient = by_individual
            .get(&recipient_id)
            .copied()
            .ok_or_else(|| BiologicalContaminationError::UnknownRecipient(recipient_id.clone()))?;
        let percent_subseed = derive_f4_technical_subseed(
            root_seed,
            F4TechnicalDecision::contamination_percent(
                event_id,
                recipient.library_id(),
                recipient.individual_id(),
            ),
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let donor_subseed = derive_f4_technical_subseed(
            root_seed,
            F4TechnicalDecision::contamination_donors(
                event_id,
                recipient.library_id(),
                recipient.individual_id(),
            ),
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let removal_subseed = derive_f4_technical_subseed(
            root_seed,
            F4TechnicalDecision::contamination_recipients(
                event_id,
                recipient.library_id(),
                recipient.individual_id(),
            ),
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let sampled_percent = f4_uniform_closed(
            &mut percent_subseed.local_rng(),
            settings.minimum_percent(),
            settings.maximum_percent(),
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let clean_pair_count = u64::try_from(recipient.own_pairs().len())
            .map_err(|_| BiologicalContaminationError::CountDoesNotFitPlatform)?;
        let target_pair_count = contamination_target_pair_count(clean_pair_count, sampled_percent)?;
        let selected_count = usize::try_from(target_pair_count)
            .map_err(|_| BiologicalContaminationError::CountDoesNotFitPlatform)?;
        let mut donor_candidates = by_individual
            .values()
            .filter(|library| {
                library.individual_id() != recipient.individual_id()
                    && library.population_id() == recipient.population_id()
            })
            .flat_map(|library| library.own_pairs().iter().cloned())
            .collect::<Vec<_>>();
        donor_candidates.sort();
        let available = u64::try_from(donor_candidates.len())
            .map_err(|_| BiologicalContaminationError::CountDoesNotFitPlatform)?;
        if target_pair_count > 0 && donor_candidates.is_empty() {
            return Err(BiologicalContaminationError::MissingDonor {
                recipient: recipient_id,
            });
        }
        if target_pair_count > available {
            return Err(BiologicalContaminationError::InsufficientDonorCapacity {
                recipient: recipient_id,
                requested: target_pair_count,
                available,
            });
        }
        f4_partial_fisher_yates(
            &mut donor_subseed.local_rng(),
            &mut donor_candidates,
            selected_count,
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let mut removed_candidates = recipient.own_pairs().to_vec();
        f4_partial_fisher_yates(
            &mut removal_subseed.local_rng(),
            &mut removed_candidates,
            selected_count,
        )
        .map_err(|_| BiologicalContaminationError::Randomness)?;
        let assignments = donor_candidates
            .into_iter()
            .take(selected_count)
            .zip(removed_candidates.into_iter().take(selected_count))
            .map(
                |(source_pair, removed_recipient_pair)| BiologicalContaminationAssignment {
                    source_pair,
                    removed_recipient_pair,
                },
            )
            .collect::<Vec<_>>();
        let donor_contributions =
            assignments
                .iter()
                .fold(BTreeMap::new(), |mut totals, assignment| {
                    *totals
                        .entry(assignment.source_pair().source_individual_id().clone())
                        .or_insert(0) += 1;
                    totals
                });
        recipients.push(RecipientBiologicalContaminationPlan {
            recipient_id,
            recipient_library_id: recipient.library_id().clone(),
            minimum_percent: settings.minimum_percent(),
            maximum_percent: settings.maximum_percent(),
            sampled_percent,
            clean_pair_count,
            target_pair_count,
            percent_subseed,
            donor_subseed,
            removal_subseed,
            donor_contributions,
            assignments,
        });
    }
    Ok(BiologicalContaminationPlan { recipients })
}

/// Convierte el porcentaje ya sorteado en sustituciones de parejas: entero más
/// próximo y empate hacia arriba. La función se comparte con los lectores de
/// manifiesto para que no acepten un resumen físicamente incoherente.
pub fn contamination_target_pair_count(
    clean_pair_count: u64,
    percent: f64,
) -> Result<u64, BiologicalContaminationError> {
    let raw = (clean_pair_count as f64) * percent / 100.0;
    if !raw.is_finite() || raw < 0.0 || raw > clean_pair_count as f64 {
        return Err(BiologicalContaminationError::CountDoesNotFitPlatform);
    }
    let floor = raw.floor();
    let rounded = if raw - floor >= 0.5 {
        floor + 1.0
    } else {
        floor
    };
    u64::try_from(rounded as u128)
        .map_err(|_| BiologicalContaminationError::CountDoesNotFitPlatform)
}
