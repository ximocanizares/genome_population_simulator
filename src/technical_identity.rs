//! Identidades técnicas mínimas y ámbitos RNG de F3.
//!
//! Esta frontera no produce ni selecciona lecturas. Solo fija las causas
//! estables de la única tanda, bibliotecas, grupos de lectura y plantillas del
//! único evento técnico de F3.

use std::collections::BTreeMap;

use crate::diagnostics::Diagnostic;
use crate::identity::{
    AutomaticIdRequest, BatchId, EventId, IdClass, IdentityRegistry, IndividualId, LibraryId,
    ReadGroupId, ReferenceId, TemplateId, compose_template_id,
};
use crate::randomness::{RootSeed, RuntimeSubseedScope, Subseed, derive_subseed};

/// Ordinal lógico de una plantilla dentro de una biblioteca. F3.1.3 deberá
/// asignarlo en su recorrido canónico de muestreo, no desde el orden físico de
/// individuos, hilos o tareas.
#[derive(Clone, Copy, Debug, Eq, Ord, PartialEq, PartialOrd)]
pub struct TemplateOrdinal(u64);

impl TemplateOrdinal {
    #[must_use]
    pub const fn new(value: u64) -> Self {
        Self(value)
    }

    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TechnicalLibraryIdentity {
    individual_id: IndividualId,
    library_id: LibraryId,
    read_group_id: ReadGroupId,
}

impl TechnicalLibraryIdentity {
    #[must_use]
    pub const fn individual_id(&self) -> &IndividualId {
        &self.individual_id
    }

    #[must_use]
    pub const fn library_id(&self) -> &LibraryId {
        &self.library_id
    }

    #[must_use]
    pub const fn read_group_id(&self) -> &ReadGroupId {
        &self.read_group_id
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TechnicalTemplateIdentity {
    template_id: TemplateId,
    template_ordinal: TemplateOrdinal,
    subseed: Subseed,
}

impl TechnicalTemplateIdentity {
    #[must_use]
    pub const fn template_id(&self) -> &TemplateId {
        &self.template_id
    }

    #[must_use]
    pub const fn template_ordinal(&self) -> TemplateOrdinal {
        self.template_ordinal
    }

    #[must_use]
    pub const fn subseed(&self) -> Subseed {
        self.subseed
    }
}

/// Plan inmutable de las entidades técnicas ya determinadas.
///
/// Las plantillas no se registran: `TemplateId` es único por construcción desde
/// una biblioteca única del evento y un ordinal único dentro de ella, por lo
/// que derivar millones de ordinales no incrementa el estado residente.
#[derive(Debug)]
pub struct TechnicalIdentityPlan {
    root_seed: RootSeed,
    reference_id: ReferenceId,
    event_id: EventId,
    batch_id: BatchId,
    event_subseed: Subseed,
    libraries: BTreeMap<IndividualId, TechnicalLibraryIdentity>,
}

impl TechnicalIdentityPlan {
    /// Referencia tipada que causó la identidad y la subsemilla del evento.
    ///
    /// Los consumidores deben compararla con la referencia del estado antes de
    /// derivar cualquier decisión aleatoria.
    #[must_use]
    pub const fn reference_id(&self) -> &ReferenceId {
        &self.reference_id
    }

    #[must_use]
    pub const fn event_id(&self) -> &EventId {
        &self.event_id
    }

    #[must_use]
    pub const fn batch_id(&self) -> &BatchId {
        &self.batch_id
    }

    #[must_use]
    pub const fn event_subseed(&self) -> Subseed {
        self.event_subseed
    }

    #[must_use]
    pub fn library(&self, individual_id: &IndividualId) -> Option<&TechnicalLibraryIdentity> {
        self.libraries.get(individual_id)
    }

    pub fn libraries(&self) -> impl Iterator<Item = &TechnicalLibraryIdentity> {
        self.libraries.values()
    }

    /// Deriva una plantilla sin consumir un flujo ni materializar una lectura.
    ///
    /// Repetir la misma causa devuelve la misma identidad; su unicidad procede
    /// de la composición canónica, no de retener identidades ya emitidas.
    pub fn derive_template(
        &self,
        individual_id: &IndividualId,
        template_ordinal: TemplateOrdinal,
    ) -> Result<TechnicalTemplateIdentity, Diagnostic> {
        let library = self.libraries.get(individual_id).ok_or_else(|| {
            Diagnostic::domain(
                "GPS-DOMAIN-001",
                "Un identificador aportado incumple la sintaxis o el ámbito de su categoría.",
            )
            .with("identifier", individual_id.as_str())
            .with("category", "individual")
            .with("scope", self.event_id.as_str())
        })?;
        let template_id = compose_template_id(library.library_id(), template_ordinal.get());
        let subseed = derive_subseed(
            self.root_seed,
            RuntimeSubseedScope::TechnicalTemplate {
                event_id: self.event_id.as_str(),
                library_id: library.library_id().as_str(),
                template_ordinal: template_ordinal.get(),
            },
        )?;
        Ok(TechnicalTemplateIdentity {
            template_id,
            template_ordinal,
            subseed,
        })
    }

    /// Deriva el flujo exclusivo que resuelve el número estocástico de
    /// plantillas de una biblioteca. No comparte estado ni subsemilla con las
    /// decisiones de origen y orientación de sus plantillas.
    pub fn derive_library_count_subseed(
        &self,
        individual_id: &IndividualId,
    ) -> Result<Subseed, Diagnostic> {
        let library = self.libraries.get(individual_id).ok_or_else(|| {
            Diagnostic::domain(
                "GPS-DOMAIN-001",
                "Un identificador aportado incumple la sintaxis o el ámbito de su categoría.",
            )
            .with("identifier", individual_id.as_str())
            .with("category", "individual")
            .with("scope", self.event_id.as_str())
        })?;
        derive_subseed(
            self.root_seed,
            RuntimeSubseedScope::TechnicalLibraryCount {
                event_id: self.event_id.as_str(),
                library_id: library.library_id().as_str(),
            },
        )
    }
}

/// Construye las identidades del único experimento F3 desde sus causas
/// estables. Los individuos se indexan en orden de identidad para que el orden
/// declarado no adquiera autoridad semántica.
///
/// La frontera requiere una referencia tipada; una ruta o una etiqueta no
/// pueden alcanzar la derivación de la subsemilla de evento.
///
/// ```compile_fail
/// use genome_population_simulator::identity::{EventId, IndividualId};
/// use genome_population_simulator::randomness::RootSeed;
/// use genome_population_simulator::technical_identity::derive_technical_identity_plan;
///
/// let root = RootSeed::parse(
///     "root",
///     "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
/// ).unwrap();
/// let event = EventId::parse("event-illumina", "test").unwrap();
/// let individuals = [IndividualId::parse("sample-A", "test").unwrap()];
/// let reference_path = "references/reference-main.fa";
/// let _ = derive_technical_identity_plan(root, reference_path, event, &individuals);
/// ```
pub fn derive_technical_identity_plan(
    root_seed: RootSeed,
    reference_id: &ReferenceId,
    event_id: EventId,
    individual_ids: &[IndividualId],
) -> Result<TechnicalIdentityPlan, Diagnostic> {
    let event_subseed = derive_subseed(
        root_seed,
        RuntimeSubseedScope::TechnicalEvent {
            event_id: event_id.as_str(),
            reference_id: reference_id.as_str(),
        },
    )?;
    let mut identities = IdentityRegistry::default();
    let batch_id = identities.derive_batch(
        root_seed,
        &AutomaticIdRequest::new(IdClass::Batch, event_id.as_str(), &[])?,
    )?;
    let mut libraries = BTreeMap::new();
    for individual_id in individual_ids {
        let library_id = identities.derive_library(
            root_seed,
            &AutomaticIdRequest::new(
                IdClass::Library,
                event_id.as_str(),
                &[("individual_id", individual_id.as_str())],
            )?,
        )?;
        let read_group_id = identities.derive_read_group(
            root_seed,
            &AutomaticIdRequest::new(
                IdClass::ReadGroup,
                event_id.as_str(),
                &[
                    ("batch_id", batch_id.as_str()),
                    ("library_id", library_id.as_str()),
                ],
            )?,
        )?;
        let identity = TechnicalLibraryIdentity {
            individual_id: individual_id.clone(),
            library_id,
            read_group_id,
        };
        if libraries.insert(individual_id.clone(), identity).is_some() {
            return Err(Diagnostic::domain(
                "GPS-DOMAIN-002",
                "La derivación de una identidad produjo una colisión.",
            )
            .with("identifier", individual_id.as_str())
            .with("category", "individual")
            .with("causes", "technical_event.individual_ids"));
        }
    }
    Ok(TechnicalIdentityPlan {
        root_seed,
        reference_id: reference_id.clone(),
        event_id,
        batch_id,
        event_subseed,
        libraries,
    })
}
