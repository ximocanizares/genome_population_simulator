//! Proyección canónica del estado F2.3 a `truth.jsonl` v1.
//!
//! Este módulo no reconstruye la verdad ni interpreta el fichero que escribe:
//! esas responsabilidades pertenecen a F2.4.3. La proyección conserva los
//! tipos del contrato externo y delega la visibilidad atómica del artefacto en
//! el workspace de publicación F1.

use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Write};
use std::path::Path;

use serde_json::{Value, json};
use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;
use crate::haplotype::HaplotypeCopy;
use crate::initial_individual_materialization::{
    InitialIndividualsMaterialization, MaterializedInitialIndividual,
};
use crate::materialization::{
    GenotypeQueryError, MaterializationError, materialize_copy, observed_snp_allele,
};
use crate::population_snvs::PreparedPopulationSnvs;
use crate::publication::ExecutionWorkspace;
use crate::reference::Reference;
use crate::results::ArtifactRecord;
use crate::segments::{
    CanonicalBase, CurrentIndelEvent, CurrentIndelPiece, CurrentSourcePiece, DeletionModification,
    InsertionModification, ModificationOrigin, SmallModification, SnpModification,
};
use crate::truth_schema::{
    AlleleCarrier, CopyIdentity, CopyProvenance, GenotypeCaseKind, GenotypeRecord,
    HaplotypeCopyRecord, MatrixProvenance, ModificationProvenance, ModificationProvenanceKind,
    NamedDigest, SegmentSource, TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND,
    TRUTH_SCHEMA_VERSION, TRUTH_SCHEMA_VERSION_WITH_REPEATS, TruthAllele, TruthContractError,
    TruthCurrentPiece, TruthCurrentSnp, TruthDocument, TruthHeader, TruthHistoryEvent, TruthLocus,
    TruthModification, TruthPieceLink, TruthRecord, TruthRepeatBlock, TruthSegment,
    TruthSourceKind, TruthSourceRef, TruthStructuralEvent, records_carry_repeat_alleles,
};
use crate::truth_validation::{
    InspectedValidationError, TruthValidationError, reference_content_sha256,
    validate_truth_with_copy_inspection,
};

pub const TRUTH_ARTIFACT_ID: &str = "truth-jsonl";
pub const TRUTH_ARTIFACT_KIND: &str = "truth_jsonl";
pub const TRUTH_ARTIFACT_PATH: &str = "outputs/truth.jsonl";

/// Entradas de procedencia que no pertenecen al estado haplotípico F2.3.
///
/// Los SHA-256 identifican los bytes de artefactos o grupos que alimentaron el
/// estado; no se inventan a partir de la proyección. La huella semántica de
/// `Reference` se calcula separadamente para `reference_content_sha256`.
/// El llamador debe incluir la referencia y cada fuente de modificación que
/// aparezca en las copias.
#[derive(Clone, Copy, Debug)]
pub struct TruthProjectionInput<'a> {
    pub execution_id: &'a str,
    pub reference: &'a Reference,
    pub individuals: &'a InitialIndividualsMaterialization,
    /// Catálogo canónico obligatorio de SNPs F2.3; también produce casos REF.
    pub matrix: &'a PreparedPopulationSnvs,
    pub source_digests: &'a [NamedDigest],
    pub group_digests: &'a [NamedDigest],
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum TruthProjectionError {
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
    MissingMatrix,
    MatrixReferenceMismatch,
    MatrixSourceMismatch,
    MissingDigest(TruthSourceRef),
    InvalidUtf8Insertion {
        event_id: String,
    },
    Genotype(GenotypeQueryError),
    Materialization(MaterializationError),
    /// El validador independiente no pudo reconstruir una copia publicada.
    IndependentReconstruction {
        identity: CopyIdentity,
        error: Box<TruthValidationError>,
    },
    /// La copia que reconstruye el validador independiente desde la verdad no
    /// es la que el productor materializa desde su estado interno. Es un
    /// defecto interno: la verdad no describiría el genoma simulado.
    ReconstructionMismatch {
        identity: CopyIdentity,
    },
    /// Hay una copia en el estado interno sin registro publicado, o al revés.
    CopySetMismatch,
    /// El juez semántico de la verdad rechaza el documento proyectado —un caso
    /// genotípico (locus, alelos, portadores o dosis) no concuerda con las
    /// copias reconstruidas— **o** su cabecera no concuerda con la referencia
    /// (identificador o huella; también si sus bases no se conservaron). Es un
    /// defecto interno del productor.
    GenotypeJudgement(Box<TruthValidationError>),
    DuplicateCaseId(String),
    /// El individuo no tiene copias de la secuencia que un locus nombra.
    /// Es imposible con una referencia coherente y se declara en vez de
    /// silenciarse.
    MissingRepeatSequence(String),
    /// Una copia lleva un alelo de repetición que el catálogo de su caso no
    /// contiene.
    ///
    /// **Es inalcanzable por construcción**: el catálogo se arma recorriendo
    /// las mismas copias que después se consultan, de modo que todo alelo
    /// observado está en él. Se declara en vez de silenciarse, y se dice aquí
    /// que ninguna prueba lo ejercita porque no hay forma de provocarlo desde
    /// fuera.
    RepeatAlleleOutsideItsCase(String),
    Contract(TruthContractError),
}

impl std::fmt::Display for TruthProjectionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self:?}")
    }
}

impl std::error::Error for TruthProjectionError {}

/// Proyecta sin E/S las copias F2.3 y el catálogo obligatorio de matriz.
/// Las colecciones resultantes se ordenan antes de validar el
/// contrato externo, por lo que el orden de entrada no afecta a los bytes.
pub fn project_truth(
    input: TruthProjectionInput<'_>,
) -> Result<TruthDocument, TruthProjectionError> {
    let reference_source = reference_source(input.reference);
    require_digest(input.source_digests, &reference_source)?;
    if input.matrix.reference_id() != input.reference.id() {
        return Err(TruthProjectionError::MatrixReferenceMismatch);
    }
    let matrix_artifact = matrix_artifact_source(input.matrix);
    require_digest(input.source_digests, &matrix_artifact)?;
    // El catálogo es obligatorio incluso para loci todo-REF. Verificar aquí
    // sus identidades impide que un caso sin ALT materializado oculte una
    // matriz suministrada/reutilizada sin procedencia publicable.
    for population_id in input.matrix.population_ids() {
        require_digest(
            input.group_digests,
            &TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: population_id.as_str().to_owned(),
            },
        )?;
    }
    for row in input.matrix.rows() {
        require_digest(
            input.source_digests,
            &TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: matrix_source_id(row.source_id()).to_owned(),
            },
        )?;
    }

    let mut records = Vec::new();
    for individual in input.individuals.individuals() {
        for sequence in individual.individual().sequences() {
            for copy in sequence.copies() {
                records.push(TruthRecord::HaplotypeCopy(project_copy(
                    individual,
                    copy,
                    input.source_digests,
                )?));
            }
        }
    }
    project_matrix_genotypes(
        &mut records,
        input.reference,
        input.individuals,
        input.matrix,
    )?;
    project_direct_genotypes(&mut records, input.reference, input.individuals)?;
    project_active_indel_genotypes(&mut records, input.reference, input.individuals)?;
    project_repeat_genotypes(&mut records, input.reference, input.individuals)?;
    records.sort_by(record_order);
    // **Una sola definición de «lleva alelos de repetición».** Estaba escrita
    // dos veces —aquí y en el contrato— y las dos versiones no miraban lo
    // mismo: esta se olvidaba del historial. Con dos respuestas distintas a la
    // misma pregunta, la versión que se publica y la que el contrato exige
    // podían discrepar, y el documento se habría rechazado a sí mismo.
    let carries_repeats = records_carry_repeat_alleles(&records);
    let copy_count = records
        .iter()
        .filter(|record| matches!(record, TruthRecord::HaplotypeCopy(_)))
        .count();
    let genotype_count = records
        .iter()
        .filter(|record| matches!(record, TruthRecord::Genotype(_)))
        .count();
    let document = TruthDocument {
        header: TruthHeader {
            document_kind: TRUTH_DOCUMENT_KIND.to_owned(),
            // La versión la decide **lo publicado**: sube a la 5 cuando algún
            // registro lleva un alelo de repetición y se queda en la 4 cuando
            // no lleva ninguno, aunque la configuración declare loci.
            schema_version: if carries_repeats {
                TRUTH_SCHEMA_VERSION_WITH_REPEATS.to_owned()
            } else {
                TRUTH_SCHEMA_VERSION.to_owned()
            },
            execution_id: input.execution_id.to_owned(),
            reference_id: input.reference.id().as_str().to_owned(),
            reference_content_sha256: reference_content_sha256(input.reference)
                .map_err(|_| TruthProjectionError::BasesNotMaterialized)?,
            coordinate_convention: TRUTH_COORDINATE_CONVENTION.to_owned(),
            copy_count: u64::try_from(copy_count).expect("record count fits u64"),
            genotype_count: u64::try_from(genotype_count).expect("record count fits u64"),
            source_digests: sorted_digests(input.source_digests),
            group_digests: sorted_digests(input.group_digests),
        },
        records,
    };
    #[cfg(test)]
    let document = {
        let mut document = document;
        test_tamper::apply(&mut document);
        document
    };
    document
        .validate_structure()
        .map_err(TruthProjectionError::Contract)?;
    verify_publication_against_simulation(&document, input.reference, input.individuals)?;
    Ok(document)
}

/// Gancho de prueba: altera el documento proyectado justo antes de que se
/// verifique, para ejercitar la verificación **a través de `project_truth`**
/// —y de la ejecución completa— con un productor defectuoso. Solo existe en
/// compilación de pruebas y es local a cada hilo.
#[cfg(test)]
pub(crate) mod test_tamper {
    use super::TruthDocument;
    use std::cell::RefCell;

    type Tamper = Box<dyn Fn(&mut TruthDocument)>;

    thread_local! {
        static TAMPER: RefCell<Option<Tamper>> = const { RefCell::new(None) };
    }

    /// Instala la alteración hasta que el guarda devuelto se suelte.
    pub(crate) fn install(tamper: impl Fn(&mut TruthDocument) + 'static) -> Guard {
        TAMPER.with(|slot| *slot.borrow_mut() = Some(Box::new(tamper)));
        Guard
    }

    pub(crate) struct Guard;

    impl Drop for Guard {
        fn drop(&mut self) {
            TAMPER.with(|slot| *slot.borrow_mut() = None);
        }
    }

    pub(super) fn apply(document: &mut TruthDocument) {
        TAMPER.with(|slot| {
            if let Some(tamper) = slot.borrow().as_ref() {
                tamper(document);
            }
        });
    }
}

/// **Lo publicado describe lo simulado, y sus casos genotípicos concuerdan con
/// sus copias.**
///
/// Hace dos contrastes con **una sola** reconstrucción por copia:
///
/// 1. Para **todas** las copias —lleven o no *indels*—, compara la huella
///    SHA-256 de la secuencia que el productor materializa desde su estado
///    interno con la de la que el validador independiente reconstruye desde el
///    registro publicado.
/// 2. Ejecuta el juez semántico de la verdad sobre el documento ya proyectado:
///    locus, alelos, portadores y dosis de cada caso genotípico contra las
///    copias reconstruidas.
///
/// Hasta el 30/09/2026 el primero lo hacía, de rebote y solo para las copias
/// con *indels*, la secuencia final que cada alelo estructural copiaba en
/// `final_haplotype_bases`; hasta el 3/10/2026 el segundo solo corría dentro de
/// la publicación del VCF, que está apagado por omisión, y así una ejecución
/// sin VCF publicó un locus de *indel* falso sin error (DT-019). No escribe
/// nada: es un contraste en memoria, una copia cada vez, que se libera antes de
/// pasar a la siguiente.
pub(crate) fn verify_publication_against_simulation(
    document: &TruthDocument,
    reference: &Reference,
    individuals: &InitialIndividualsMaterialization,
) -> Result<(), TruthProjectionError> {
    let simulated = individuals
        .individuals()
        .iter()
        .flat_map(|individual| individual.individual().sequences())
        .flat_map(|sequence| sequence.copies())
        .map(|copy| (copy_identity(copy), copy))
        .collect::<BTreeMap<_, _>>();
    let mut published = 0_usize;
    validate_truth_with_copy_inspection(document, reference, |record, reconstructed| {
        published += 1;
        let copy = simulated
            .get(&record.identity)
            .ok_or(TruthProjectionError::CopySetMismatch)?;
        let materialized: [u8; 32] = Sha256::digest(
            materialize_copy(copy, reference).map_err(TruthProjectionError::Materialization)?,
        )
        .into();
        if reconstructed != materialized {
            return Err(TruthProjectionError::ReconstructionMismatch {
                identity: record.identity.clone(),
            });
        }
        Ok(())
    })
    .map_err(|error| match error {
        InspectedValidationError::CopyReconstruction { identity, error } => {
            TruthProjectionError::IndependentReconstruction { identity, error }
        }
        InspectedValidationError::Semantics(error) => {
            TruthProjectionError::GenotypeJudgement(error)
        }
        InspectedValidationError::Inspection(error) => error,
    })?;
    if published != simulated.len() {
        return Err(TruthProjectionError::CopySetMismatch);
    }
    Ok(())
}

fn project_copy(
    individual: &MaterializedInitialIndividual,
    copy: &HaplotypeCopy,
    source_digests: &[NamedDigest],
) -> Result<HaplotypeCopyRecord, TruthProjectionError> {
    let identity = copy_identity(copy);
    let mut sources = BTreeSet::new();
    sources.insert(reference_source_from_id(copy.id().reference_id().as_str()));
    let mut segments = Vec::with_capacity(copy.segments().segments().len());
    for segment in copy.segments().segments() {
        let source = reference_source_from_id(segment.source().reference_id().as_str());
        require_digest(source_digests, &source)?;
        let modifications = segment
            .modifications()
            .iter()
            .map(|modification| project_modification(modification, source_digests))
            .collect::<Result<Vec<_>, _>>()?;
        for modification in &modifications {
            sources.insert(modification_source(modification).clone());
        }
        segments.push(TruthSegment {
            source: SegmentSource {
                source,
                reference_id: segment.source().reference_id().as_str().to_owned(),
                sequence_id: segment.source().sequence_id().as_str().to_owned(),
                start: segment.source().start(),
                end: segment.source().end(),
                sequence_length: segment.source().sequence_length(),
            },
            modifications,
        });
    }
    let current_pieces = copy
        .segments()
        .current_modifications()
        .iter()
        .map(|piece| project_current_piece(piece, source_digests))
        .collect::<Result<Vec<_>, _>>()?;
    for piece in &current_pieces {
        match piece {
            TruthCurrentPiece::Source { source, snps } => {
                sources.insert(source.source.clone());
                for snp in snps {
                    sources.insert(snp.provenance.source.clone());
                }
            }
            TruthCurrentPiece::Inserted { provenance, .. }
            | TruthCurrentPiece::RepeatAllele { provenance, .. } => {
                sources.insert(provenance.source.clone());
            }
        }
    }
    let history = copy
        .segments()
        .small_event_history()
        .iter()
        .map(|event| project_history_event(event, source_digests))
        .collect::<Result<Vec<_>, _>>()?;
    for event in &history {
        let provenance = match event {
            TruthHistoryEvent::Snp { provenance, .. }
            | TruthHistoryEvent::Insertion { provenance, .. }
            | TruthHistoryEvent::Deletion { provenance, .. }
            | TruthHistoryEvent::RepeatAllele { provenance, .. } => provenance,
        };
        sources.insert(provenance.source.clone());
    }
    Ok(HaplotypeCopyRecord {
        identity,
        phase: copy.phase().copy_index(),
        provenance: CopyProvenance {
            origin: format!("initial_population:{}", individual.population_id().as_str()),
            sources: sources.into_iter().collect(),
        },
        segments,
        current_pieces,
        history,
    })
}

fn project_current_piece(
    piece: &CurrentIndelPiece,
    source_digests: &[NamedDigest],
) -> Result<TruthCurrentPiece, TruthProjectionError> {
    match piece {
        CurrentIndelPiece::Source(source) => project_current_source(source, source_digests),
        CurrentIndelPiece::Inserted(inserted) => {
            let event_id = inserted.provenance().event_id().as_str().to_owned();
            let source = domain_source(inserted.provenance());
            require_digest(source_digests, &source)?;
            Ok(TruthCurrentPiece::Inserted {
                bases: bytes_to_truth(&event_id, inserted.bases())?,
                provenance: domain_provenance(inserted.provenance(), source),
            })
        }
        CurrentIndelPiece::RepeatAllele(allele) => {
            let event_id = allele.provenance().event_id().as_str().to_owned();
            let source = domain_source(allele.provenance());
            require_digest(source_digests, &source)?;
            let reference_source =
                reference_source_from_id(allele.interval().reference_id().as_str());
            require_digest(source_digests, &reference_source)?;
            Ok(TruthCurrentPiece::RepeatAllele {
                locus_id: allele.locus_id().as_str().to_owned(),
                source: segment_source(&reference_source, allele.interval()),
                bases: bytes_to_truth(&event_id, allele.tract().bases())?,
                blocks: repeat_blocks(allele.tract()),
                provenance: domain_provenance(allele.provenance(), source),
            })
        }
    }
}

/// Estructura en bloques publicada: motivo y **bases**, nunca unidades.
///
/// El número de unidades es `bases / motivo`, y se calcula al serializar. Si
/// se guardara junto a las bases podrían contradecirse.
fn repeat_blocks(tract: &crate::tandem_repeat::TandemRepeatTract) -> Vec<TruthRepeatBlock> {
    tract
        .blocks()
        .iter()
        .map(|block| TruthRepeatBlock {
            motif: block.motif().as_text(),
            bases: block.length(),
        })
        .collect()
}

fn segment_source(
    source: &TruthSourceRef,
    interval: &crate::coordinates::GenomicInterval,
) -> SegmentSource {
    SegmentSource {
        source: source.clone(),
        reference_id: interval.reference_id().as_str().to_owned(),
        sequence_id: interval.sequence_id().as_str().to_owned(),
        start: interval.start(),
        end: interval.end(),
        sequence_length: interval.sequence_length(),
    }
}

fn project_current_source(
    source: &CurrentSourcePiece,
    source_digests: &[NamedDigest],
) -> Result<TruthCurrentPiece, TruthProjectionError> {
    let source_ref = reference_source_from_id(source.source().reference_id().as_str());
    require_digest(source_digests, &source_ref)?;
    let snps = source
        .snps()
        .iter()
        .map(|snp| {
            Ok(TruthCurrentSnp {
                event_id: snp.provenance().event_id().as_str().to_owned(),
                position: snp.position().offset(),
                reference: base_text(snp.reference()),
                alternate: base_text(snp.alternate()),
                provenance: snp_provenance(snp, source_digests)?,
            })
        })
        .collect::<Result<Vec<_>, TruthProjectionError>>()?;
    Ok(TruthCurrentPiece::Source {
        source: SegmentSource {
            source: source_ref,
            reference_id: source.source().reference_id().as_str().to_owned(),
            sequence_id: source.source().sequence_id().as_str().to_owned(),
            start: source.source().start(),
            end: source.source().end(),
            sequence_length: source.source().sequence_length(),
        },
        snps,
    })
}

fn project_history_event(
    event: &CurrentIndelEvent,
    source_digests: &[NamedDigest],
) -> Result<TruthHistoryEvent, TruthProjectionError> {
    match event {
        CurrentIndelEvent::Snp { snp } => Ok(TruthHistoryEvent::Snp {
            event_id: snp.provenance().event_id().as_str().to_owned(),
            position: snp.position().offset(),
            provenance: snp_provenance(snp, source_digests)?,
        }),
        CurrentIndelEvent::Insertion {
            event_id,
            point,
            bases,
            provenance,
            ..
        } => {
            let source = domain_source(provenance);
            require_digest(source_digests, &source)?;
            Ok(TruthHistoryEvent::Insertion {
                event_id: event_id.as_str().to_owned(),
                point: *point,
                bases: bytes_to_truth(event_id.as_str(), bases)?,
                provenance: domain_provenance(provenance, source),
            })
        }
        CurrentIndelEvent::Deletion {
            event_id,
            start,
            end,
            provenance,
            ..
        } => {
            let source = domain_source(provenance);
            require_digest(source_digests, &source)?;
            Ok(TruthHistoryEvent::Deletion {
                event_id: event_id.as_str().to_owned(),
                start: *start,
                end: *end,
                provenance: domain_provenance(provenance, source),
            })
        }
        CurrentIndelEvent::RepeatAllele { allele } => {
            let event_id = allele.provenance().event_id().as_str().to_owned();
            let source = domain_source(allele.provenance());
            require_digest(source_digests, &source)?;
            Ok(TruthHistoryEvent::RepeatAllele {
                event_id: event_id.clone(),
                locus_id: allele.locus_id().as_str().to_owned(),
                start: allele.interval().start(),
                end: allele.interval().end(),
                bases: bytes_to_truth(&event_id, allele.tract().bases())?,
                blocks: repeat_blocks(allele.tract()),
                provenance: domain_provenance(allele.provenance(), source),
            })
        }
    }
}

fn project_modification(
    modification: &SmallModification,
    source_digests: &[NamedDigest],
) -> Result<TruthModification, TruthProjectionError> {
    match modification {
        SmallModification::Snp(snp) => project_snp(snp, source_digests),
        SmallModification::Insertion(insertion) => project_insertion(insertion, source_digests),
        SmallModification::Deletion(deletion) => project_deletion(deletion, source_digests),
    }
}

fn project_snp(
    snp: &SnpModification,
    source_digests: &[NamedDigest],
) -> Result<TruthModification, TruthProjectionError> {
    let provenance = snp_provenance(snp, source_digests)?;
    Ok(TruthModification::Snp {
        event_id: snp.provenance().event_id().as_str().to_owned(),
        position: snp.position().offset(),
        reference: base_text(snp.reference()),
        alternate: base_text(snp.alternate()),
        provenance,
    })
}

fn project_insertion(
    insertion: &InsertionModification,
    source_digests: &[NamedDigest],
) -> Result<TruthModification, TruthProjectionError> {
    let event_id = insertion.provenance().event_id().as_str().to_owned();
    let source = domain_source(insertion.provenance());
    require_digest(source_digests, &source)?;
    let bases = std::str::from_utf8(insertion.bases())
        .map_err(|_| TruthProjectionError::InvalidUtf8Insertion {
            event_id: event_id.clone(),
        })?
        .to_owned();
    Ok(TruthModification::Insertion {
        event_id,
        point: insertion.anchor().point().offset(),
        bases,
        provenance: domain_provenance(insertion.provenance(), source),
    })
}

fn project_deletion(
    deletion: &DeletionModification,
    source_digests: &[NamedDigest],
) -> Result<TruthModification, TruthProjectionError> {
    let event_id = deletion.provenance().event_id().as_str().to_owned();
    let source = domain_source(deletion.provenance());
    require_digest(source_digests, &source)?;
    Ok(TruthModification::Deletion {
        event_id,
        start: deletion.deletion().interval().start(),
        end: deletion.deletion().interval().end(),
        provenance: domain_provenance(deletion.provenance(), source),
    })
}

fn snp_provenance(
    snp: &SnpModification,
    source_digests: &[NamedDigest],
) -> Result<ModificationProvenance, TruthProjectionError> {
    match (
        snp.matrix_source_id(),
        snp.matrix_population_id(),
        snp.matrix_variant_id(),
    ) {
        (Some(source_id), Some(population_id), Some(variant_id)) => {
            let source_id = matrix_source_id(source_id);
            let source = TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: source_id.to_owned(),
            };
            require_digest(source_digests, &source)?;
            Ok(ModificationProvenance {
                kind: ModificationProvenanceKind::Matrix,
                version: "1".to_owned(),
                source: source.clone(),
                matrix: Some(MatrixProvenance {
                    matrix_source_id: source_id.to_owned(),
                    population_id: population_id.as_str().to_owned(),
                    variant_id: variant_id.as_str().to_owned(),
                }),
            })
        }
        (None, None, None) => {
            let source = domain_source(snp.provenance());
            require_digest(source_digests, &source)?;
            Ok(domain_provenance(snp.provenance(), source))
        }
        _ => Err(TruthProjectionError::Contract(
            TruthContractError::InvalidModificationProvenance {
                event_id: snp.provenance().event_id().as_str().to_owned(),
            },
        )),
    }
}

fn project_matrix_genotypes(
    records: &mut Vec<TruthRecord>,
    reference: &Reference,
    individuals: &InitialIndividualsMaterialization,
    matrix: &PreparedPopulationSnvs,
) -> Result<(), TruthProjectionError> {
    let mut loci = BTreeMap::<(String, String, u64, u8), Vec<_>>::new();
    for row in matrix.rows() {
        loci.entry((
            row.population_id().as_str().to_owned(),
            row.contig().as_str().to_owned(),
            row.position_1based(),
            row.reference(),
        ))
        .or_default()
        .push(row);
    }
    let mut cases = BTreeSet::new();
    for ((population_id, sequence_id, position_1based, reference_base), mut rows) in loci {
        rows.sort_by(|left, right| {
            left.alternate()
                .cmp(&right.alternate())
                .then(left.variant_id().cmp(right.variant_id()))
        });
        let mut alleles = vec![TruthAllele::Reference];
        for row in &rows {
            let allele = TruthAllele::Snp {
                base: char::from(row.alternate()).to_string(),
            };
            if !alleles.contains(&allele) {
                alleles.push(allele);
            }
        }
        for individual in individuals
            .individuals()
            .iter()
            .filter(|individual| individual.population_id().as_str() == population_id)
        {
            let case_id = case_id(
                &population_id,
                individual.individual().id().as_str(),
                &sequence_id,
                position_1based,
                &alleles,
            );
            if !cases.insert(case_id.clone()) {
                return Err(TruthProjectionError::DuplicateCaseId(case_id));
            }
            let sequence = reference
                .sequences()
                .iter()
                .find(|sequence| sequence.id().as_str() == sequence_id)
                .expect("matrix was prepared against this reference");
            let position = crate::coordinates::BasePosition::from_base_1(
                reference.id().clone(),
                sequence.id().clone(),
                position_1based,
                sequence.length(),
            )
            .expect("prepared matrix positions are valid");
            let reference_base = CanonicalBase::from_byte(reference_base)
                .expect("prepared matrix uses canonical bases");
            let mut observed = Vec::new();
            for copy in individual
                .individual()
                .copies_for(sequence.id())
                .expect("individual derives from this reference")
            {
                let allele = match observed_snp_allele(copy, &position, reference_base) {
                    Ok(allele) if allele == reference_base => TruthAllele::Reference,
                    Ok(allele) => TruthAllele::Snp {
                        base: char::from(allele.as_byte()).to_string(),
                    },
                    // Una deleción no es un alelo SNP. Su estado se publica en
                    // el caso estructural tipado que la ocasionó.
                    Err(GenotypeQueryError::LocusDeletedFromCopy { .. }) => continue,
                    Err(error) => return Err(TruthProjectionError::Genotype(error)),
                };
                if !alleles.contains(&allele) {
                    alleles.push(allele.clone());
                }
                observed.push((copy, allele));
            }
            let mut carriers = Vec::new();
            for (copy, allele) in observed {
                let allele_index =
                    alleles
                        .iter()
                        .position(|text| text == &allele)
                        .ok_or_else(|| {
                            TruthProjectionError::Genotype(
                                GenotypeQueryError::IncompatibleSnpInCopy {
                                    phase: copy.phase(),
                                },
                            )
                        })?;
                carriers.push(AlleleCarrier {
                    allele_index,
                    copy: copy_identity(copy),
                });
            }
            carriers.sort_by(|left, right| left.copy.cmp(&right.copy));
            let mut dosages = vec![0_u8; alleles.len()];
            for carrier in &carriers {
                dosages[carrier.allele_index] = dosages[carrier.allele_index]
                    .checked_add(1)
                    .expect("ploidy is at most ten");
            }
            records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id,
                reference_id: reference.id().as_str().to_owned(),
                individual_id: individual.individual().id().as_str().to_owned(),
                case_kind: GenotypeCaseKind::Snp,
                locus: TruthLocus {
                    sequence_id: sequence_id.clone(),
                    start: position.offset(),
                    end: position.offset() + 1,
                },
                alleles: alleles.clone(),
                allele_dosages: dosages,
                carriers,
            }));
        }
    }
    Ok(())
}

/// Añade los casos SNP vigentes que proceden de cambios directos y no estaban
/// ya cubiertos por una fila activa de matriz. El catálogo se toma de los
/// fragmentos vigentes, no del historial: un SNP eliminado sigue trazado pero
/// deja de ser un caso genotípico activo.
fn project_direct_genotypes(
    records: &mut Vec<TruthRecord>,
    reference: &Reference,
    individuals: &InitialIndividualsMaterialization,
) -> Result<(), TruthProjectionError> {
    let existing = records
        .iter()
        .filter_map(|record| match record {
            TruthRecord::Genotype(genotype) => Some((
                genotype.reference_id.clone(),
                genotype.individual_id.clone(),
                genotype.locus.clone(),
            )),
            TruthRecord::HaplotypeCopy(_) => None,
        })
        .collect::<BTreeSet<_>>();
    let mut cases = BTreeSet::new();
    for individual in individuals.individuals() {
        let mut loci = BTreeMap::<(String, u64, u8), BTreeSet<u8>>::new();
        for sequence in individual.individual().sequences() {
            for copy in sequence.copies() {
                for piece in copy.segments().current_modifications() {
                    if let CurrentIndelPiece::Source(source) = piece {
                        for snp in source.snps() {
                            loci.entry((
                                source.source().sequence_id().as_str().to_owned(),
                                snp.position().offset(),
                                snp.reference().as_byte(),
                            ))
                            .or_default()
                            .insert(snp.alternate().as_byte());
                        }
                    }
                }
            }
        }
        for ((sequence_id, position, reference_base), alternates) in loci {
            let locus = TruthLocus {
                sequence_id: sequence_id.clone(),
                start: position,
                end: position + 1,
            };
            if existing.contains(&(
                reference.id().as_str().to_owned(),
                individual.individual().id().as_str().to_owned(),
                locus.clone(),
            )) {
                continue;
            }
            let mut alleles = vec![TruthAllele::Reference];
            alleles.extend(alternates.into_iter().map(|base| TruthAllele::Snp {
                base: char::from(base).to_string(),
            }));
            let sequence = reference
                .sequences()
                .iter()
                .find(|sequence| sequence.id().as_str() == sequence_id)
                .expect("current source belongs to the individual reference");
            let position = crate::coordinates::BasePosition::from_base_1(
                reference.id().clone(),
                sequence.id().clone(),
                position + 1,
                sequence.length(),
            )
            .expect("current SNP position belongs to the source sequence");
            let reference_base = CanonicalBase::from_byte(reference_base)
                .expect("current SNP uses a canonical base");
            let mut observed = Vec::new();
            for copy in individual
                .individual()
                .copies_for(sequence.id())
                .expect("individual derives from this reference")
            {
                let allele = match observed_snp_allele(copy, &position, reference_base) {
                    Ok(allele) if allele == reference_base => TruthAllele::Reference,
                    Ok(allele) => TruthAllele::Snp {
                        base: char::from(allele.as_byte()).to_string(),
                    },
                    Err(GenotypeQueryError::LocusDeletedFromCopy { .. }) => continue,
                    Err(error) => return Err(TruthProjectionError::Genotype(error)),
                };
                if !alleles.contains(&allele) {
                    alleles.push(allele.clone());
                }
                observed.push((copy, allele));
            }
            let case_id = case_id(
                individual.population_id().as_str(),
                individual.individual().id().as_str(),
                &sequence_id,
                locus.start + 1,
                &alleles,
            );
            if !cases.insert(case_id.clone()) {
                return Err(TruthProjectionError::DuplicateCaseId(case_id));
            }
            let mut carriers = Vec::new();
            for (copy, allele) in observed {
                let allele_index = alleles
                    .iter()
                    .position(|value| value == &allele)
                    .expect("observed allele is inserted into the case catalog");
                carriers.push(AlleleCarrier {
                    allele_index,
                    copy: copy_identity(copy),
                });
            }
            carriers.sort_by(|left, right| left.copy.cmp(&right.copy));
            let mut dosages = vec![0_u8; alleles.len()];
            for carrier in &carriers {
                dosages[carrier.allele_index] = dosages[carrier.allele_index]
                    .checked_add(1)
                    .expect("ploidy is at most ten");
            }
            records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id,
                reference_id: reference.id().as_str().to_owned(),
                individual_id: individual.individual().id().as_str().to_owned(),
                case_kind: GenotypeCaseKind::Snp,
                locus,
                alleles,
                allele_dosages: dosages,
                carriers,
            }));
        }
    }
    Ok(())
}

/// Publica cada indel que sigue activo como un caso diploide/poliploide
/// estructural. Los portadores incluyen todas las copias: la copia a la que
/// se aplicó el evento lleva el alelo tipado y el resto lleva `reference`.
/// Así, una deleción que cubre un locus SNP nunca se convierte en el texto
/// artificial `"<deleted>"` dentro del caso SNP.
fn project_active_indel_genotypes(
    records: &mut Vec<TruthRecord>,
    reference: &Reference,
    individuals: &InitialIndividualsMaterialization,
) -> Result<(), TruthProjectionError> {
    for individual in individuals.individuals() {
        let mut active_insertions = BTreeSet::new();
        for sequence in individual.individual().sequences() {
            for copy in sequence.copies() {
                for piece in copy.segments().current_modifications() {
                    if let CurrentIndelPiece::Inserted(inserted) = piece {
                        active_insertions.insert((
                            copy.id().sequence_id().as_str().to_owned(),
                            inserted.provenance().event_id().as_str().to_owned(),
                        ));
                    }
                }
            }
        }
        let mut cases =
            BTreeMap::<(String, String), (TruthLocus, TruthAllele, BTreeSet<CopyIdentity>)>::new();
        for sequence in individual.individual().sequences() {
            for copy in sequence.copies() {
                for event in copy.segments().small_event_history() {
                    let sequence_id = copy.id().sequence_id().as_str().to_owned();
                    let (event_id, locus, allele, active) = match event {
                        // Los alelos de repetición tienen su propio caso, con
                        // su locus, su motivo y su número de unidades: no son
                        // una inserción ni una deleción anónima.
                        CurrentIndelEvent::Snp { .. } | CurrentIndelEvent::RepeatAllele { .. } => {
                            continue;
                        }
                        // El locus del caso va en coordenadas de la
                        // referencia; el evento fuente conserva las vigentes
                        // con las que se declaró y se replica el historial.
                        CurrentIndelEvent::Insertion {
                            event_id,
                            point,
                            bases,
                            reference_point,
                            ..
                        } => {
                            let event_id = event_id.as_str().to_owned();
                            let bases = bytes_to_truth(&event_id, bases)?;
                            (
                                event_id.clone(),
                                TruthLocus {
                                    sequence_id: sequence_id.clone(),
                                    start: *reference_point,
                                    end: *reference_point,
                                },
                                TruthAllele::Insertion {
                                    event_id: event_id.clone(),
                                    bases: bases.clone(),
                                    source_event: TruthStructuralEvent::Insertion {
                                        event_id: event_id.clone(),
                                        point: *point,
                                        bases,
                                    },
                                    piece_links: current_piece_links(copy, &event_id),
                                },
                                active_insertions.contains(&(sequence_id.clone(), event_id)),
                            )
                        }
                        CurrentIndelEvent::Deletion {
                            event_id,
                            start,
                            end,
                            reference_start,
                            reference_end,
                            ..
                        } => {
                            let event_id = event_id.as_str().to_owned();
                            (
                                event_id.clone(),
                                TruthLocus {
                                    sequence_id: sequence_id.clone(),
                                    start: *reference_start,
                                    end: *reference_end,
                                },
                                TruthAllele::Deletion {
                                    event_id: event_id.clone(),
                                    start: *start,
                                    end: *end,
                                    source_event: TruthStructuralEvent::Deletion {
                                        event_id: event_id.clone(),
                                        start: *start,
                                        end: *end,
                                    },
                                    piece_links: current_piece_links(copy, &event_id),
                                },
                                true,
                            )
                        }
                    };
                    if !active {
                        continue;
                    }
                    let entry = cases
                        .entry((sequence_id, event_id))
                        .or_insert_with(|| (locus, allele, BTreeSet::new()));
                    entry.2.insert(copy_identity(copy));
                }
            }
        }
        for ((_sequence_id, _event_id), (locus, alternate, alternate_copies)) in cases {
            let alleles = vec![TruthAllele::Reference, alternate];
            // A stable allele order is contractual even if a future domain
            // event happens to name an alternate that sorts before REF.
            debug_assert!(matches!(alleles[0], TruthAllele::Reference));
            let case_id = case_id(
                individual.population_id().as_str(),
                individual.individual().id().as_str(),
                &locus.sequence_id,
                locus.start + 1,
                &alleles,
            );
            let mut carriers = Vec::new();
            for copy in individual
                .individual()
                .copies_for(
                    &crate::identity::SequenceId::parse(&locus.sequence_id, "truth projection")
                        .expect("canonical sequence id"),
                )
                .expect("individual derives from this reference")
            {
                carriers.push(AlleleCarrier {
                    allele_index: usize::from(alternate_copies.contains(&copy_identity(copy))),
                    copy: copy_identity(copy),
                });
            }
            carriers.sort_by(|left, right| left.copy.cmp(&right.copy));
            let mut allele_dosages = vec![0_u8; alleles.len()];
            for carrier in &carriers {
                allele_dosages[carrier.allele_index] = allele_dosages[carrier.allele_index]
                    .checked_add(1)
                    .expect("ploidy is at most ten");
            }
            records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id,
                reference_id: reference.id().as_str().to_owned(),
                individual_id: individual.individual().id().as_str().to_owned(),
                case_kind: GenotypeCaseKind::Indel,
                locus,
                alleles,
                allele_dosages,
                carriers,
            }));
        }
    }
    Ok(())
}

/// Publica un caso por locus de repetición con al menos una copia que **no**
/// lleva el tracto de la referencia.
///
/// Un locus cuyas copias llevan todas el alelo de la referencia no produce
/// ningún caso: no hay variante que publicar, y por eso su verdad sigue siendo
/// byte a byte la de un escenario sin microsatélites.
fn project_repeat_genotypes(
    records: &mut Vec<TruthRecord>,
    reference: &Reference,
    individuals: &InitialIndividualsMaterialization,
) -> Result<(), TruthProjectionError> {
    for individual in individuals.individuals() {
        // Clave: secuencia e intervalo del locus, que es lo que identifica al
        // caso. El identificador del locus viaja dentro de cada alelo.
        let mut cases = BTreeMap::<(String, u64, u64), Vec<TruthAllele>>::new();
        for sequence in individual.individual().sequences() {
            for copy in sequence.copies() {
                for piece in copy.segments().current_modifications() {
                    let CurrentIndelPiece::RepeatAllele(allele) = piece else {
                        continue;
                    };
                    let key = (
                        copy.id().sequence_id().as_str().to_owned(),
                        allele.interval().start(),
                        allele.interval().end(),
                    );
                    let catalogued = TruthAllele::Repeat {
                        locus_id: allele.locus_id().as_str().to_owned(),
                        bases: bytes_to_truth(allele.locus_id().as_str(), allele.tract().bases())?,
                        blocks: repeat_blocks(allele.tract()),
                    };
                    let entry = cases
                        .entry(key)
                        .or_insert_with(|| vec![TruthAllele::Reference]);
                    if !entry.contains(&catalogued) {
                        entry.push(catalogued);
                    }
                }
            }
        }
        for ((sequence_id, start, end), alleles) in cases {
            let locus = TruthLocus {
                sequence_id: sequence_id.clone(),
                start,
                end,
            };
            let case_id = case_id(
                individual.population_id().as_str(),
                individual.individual().id().as_str(),
                &sequence_id,
                start + 1,
                &alleles,
            );
            let sequence_key = crate::identity::SequenceId::parse(&sequence_id, "truth projection")
                .map_err(|_| TruthProjectionError::MissingRepeatSequence(sequence_id.clone()))?;
            let copies = individual
                .individual()
                .copies_for(&sequence_key)
                .ok_or_else(|| TruthProjectionError::MissingRepeatSequence(sequence_id.clone()))?;
            let mut carriers = Vec::with_capacity(copies.len());
            for copy in copies {
                let observed = copy
                    .segments()
                    .current_modifications()
                    .iter()
                    .find_map(|piece| match piece {
                        CurrentIndelPiece::RepeatAllele(allele)
                            if allele.interval().start() == start
                                && allele.interval().end() == end =>
                        {
                            Some(allele)
                        }
                        CurrentIndelPiece::RepeatAllele(_)
                        | CurrentIndelPiece::Source(_)
                        | CurrentIndelPiece::Inserted(_) => None,
                    });
                let allele_index = match observed {
                    None => 0,
                    Some(allele) => {
                        let catalogued = TruthAllele::Repeat {
                            locus_id: allele.locus_id().as_str().to_owned(),
                            bases: bytes_to_truth(
                                allele.locus_id().as_str(),
                                allele.tract().bases(),
                            )?,
                            blocks: repeat_blocks(allele.tract()),
                        };
                        alleles
                            .iter()
                            .position(|value| value == &catalogued)
                            .ok_or_else(|| {
                                TruthProjectionError::RepeatAlleleOutsideItsCase(case_id.clone())
                            })?
                    }
                };
                carriers.push(AlleleCarrier {
                    allele_index,
                    copy: copy_identity(copy),
                });
            }
            carriers.sort_by(|left, right| left.copy.cmp(&right.copy));
            let mut allele_dosages = vec![0_u8; alleles.len()];
            for carrier in &carriers {
                allele_dosages[carrier.allele_index] = allele_dosages[carrier.allele_index]
                    .checked_add(1)
                    .expect("ploidy is at most ten");
            }
            records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id,
                reference_id: reference.id().as_str().to_owned(),
                individual_id: individual.individual().id().as_str().to_owned(),
                case_kind: GenotypeCaseKind::Repeat,
                locus,
                alleles,
                allele_dosages,
                carriers,
            }));
        }
    }
    Ok(())
}

/// Enlaces a los fragmentos vigentes que todavía proceden de una inserción.
/// Las coordenadas son offsets de la secuencia final de la copia y hacen
/// explícita una inserción dividida por cambios posteriores.
fn current_piece_links(copy: &HaplotypeCopy, event_id: &str) -> Vec<TruthPieceLink> {
    let mut cursor = 0_u64;
    copy.segments()
        .current_modifications()
        .iter()
        .enumerate()
        .filter_map(|(index, piece)| {
            let len = match piece {
                CurrentIndelPiece::Source(source) => {
                    source.source().end() - source.source().start()
                }
                CurrentIndelPiece::Inserted(inserted) => {
                    u64::try_from(inserted.bases().len()).expect("piece length fits u64")
                }
                CurrentIndelPiece::RepeatAllele(allele) => allele.tract().length(),
            };
            let start = cursor;
            cursor += len;
            match piece {
                CurrentIndelPiece::Inserted(inserted)
                    if inserted.provenance().event_id().as_str() == event_id =>
                {
                    Some(TruthPieceLink {
                        piece_index: u64::try_from(index).expect("piece index fits u64"),
                        start,
                        end: cursor,
                    })
                }
                _ => None,
            }
        })
        .collect()
}

/// Escribe JSON Lines UTF-8 canónico. Siempre añade LF tras cada registro,
/// incluida la última línea, y no cierra ni trunca el escritor prestado.
pub fn write_truth_jsonl(document: &TruthDocument, writer: &mut impl Write) -> io::Result<()> {
    write_jsonl_value(writer, &header_json(&document.header))?;
    for record in &document.records {
        write_jsonl_value(writer, &record_json(record))?;
    }
    Ok(())
}

fn write_jsonl_value(writer: &mut impl Write, value: &Value) -> io::Result<()> {
    serde_json::to_writer(&mut *writer, value).map_err(io::Error::other)?;
    writer.write_all(b"\n")
}

/// Escribe la verdad por flujo dentro del temporal de F1. Si el flujo falla,
/// `ExecutionWorkspace` elimina su temporal, no registra el artefacto y deja
/// intacto el marcador `INCOMPLETE` del provisional.
pub fn write_truth_artifact(
    workspace: &mut ExecutionWorkspace,
    document: &TruthDocument,
) -> Result<ArtifactRecord, Diagnostic> {
    document.validate_structure().map_err(contract_diagnostic)?;
    workspace.write_streamed_artifact(
        TRUTH_ARTIFACT_ID,
        TRUTH_ARTIFACT_KIND,
        Path::new(TRUTH_ARTIFACT_PATH),
        |file| {
            write_truth_jsonl(document, file).map_err(|error| {
                Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
                    .with("path", TRUTH_ARTIFACT_PATH)
                    .with("operation", "write_truth_jsonl")
                    .with("detail", error.to_string())
            })
        },
    )
}

fn header_json(header: &TruthHeader) -> Value {
    json!({
        "record_type": "header", "document_kind": header.document_kind,
        "schema_version": header.schema_version, "execution_id": header.execution_id,
        "reference_id": header.reference_id, "reference_content_sha256": header.reference_content_sha256,
        "coordinate_convention": header.coordinate_convention,
        "counts": {"copies": header.copy_count, "genotypes": header.genotype_count},
        "source_digests": header.source_digests.iter().map(digest_json).collect::<Vec<_>>(),
        "group_digests": header.group_digests.iter().map(digest_json).collect::<Vec<_>>(),
    })
}

fn record_json(record: &TruthRecord) -> Value {
    match record {
        TruthRecord::HaplotypeCopy(copy) => json!({
            "record_type":"haplotype_copy", "identity": copy_identity_json(&copy.identity), "phase":copy.phase,
            "provenance":{"origin":copy.provenance.origin,"sources":copy.provenance.sources.iter().map(source_json).collect::<Vec<_>>()},
            "segments": copy.segments.iter().map(|segment| json!({
                "source":{"source":source_json(&segment.source.source), "reference_id":segment.source.reference_id,
                    "sequence_id":segment.source.sequence_id,"start":segment.source.start,"end":segment.source.end,"sequence_length":segment.source.sequence_length},
                "modifications":segment.modifications.iter().map(modification_json).collect::<Vec<_>>()
            })).collect::<Vec<_>>(),
            "current_pieces":copy.current_pieces.iter().map(current_piece_json).collect::<Vec<_>>(),
            "history":copy.history.iter().map(history_event_json).collect::<Vec<_>>()
        }),
        TruthRecord::Genotype(genotype) => json!({
            "record_type":"genotype", "case_id":genotype.case_id, "reference_id":genotype.reference_id,
            "individual_id":genotype.individual_id, "case_kind":genotype_case_kind_json(genotype.case_kind),
            "locus":{"sequence_id":genotype.locus.sequence_id,"start":genotype.locus.start,"end":genotype.locus.end},
            "alleles":genotype.alleles.iter().map(allele_json).collect::<Vec<_>>(), "allele_dosages":genotype.allele_dosages,
            "carriers":genotype.carriers.iter().map(|carrier| json!({"allele_index":carrier.allele_index,"copy":copy_identity_json(&carrier.copy)})).collect::<Vec<_>>()
        }),
    }
}
fn genotype_case_kind_json(kind: GenotypeCaseKind) -> &'static str {
    match kind {
        GenotypeCaseKind::Snp => "snp",
        GenotypeCaseKind::Indel => "indel",
        GenotypeCaseKind::Repeat => "repeat",
    }
}
fn allele_json(allele: &TruthAllele) -> Value {
    match allele {
        TruthAllele::Reference => json!({"kind":"reference"}),
        TruthAllele::Snp { base } => json!({"kind":"snp","base":base}),
        TruthAllele::Insertion {
            event_id,
            bases,
            source_event,
            piece_links,
        } => json!({"kind":"insertion","event_id":event_id,"bases":bases,
            "source_event":structural_event_json(source_event),
            "piece_links":piece_links.iter().map(piece_link_json).collect::<Vec<_>>() }),
        TruthAllele::Deletion {
            event_id,
            start,
            end,
            source_event,
            piece_links,
        } => json!({"kind":"deletion","event_id":event_id,"start":start,"end":end,
            "source_event":structural_event_json(source_event),
            "piece_links":piece_links.iter().map(piece_link_json).collect::<Vec<_>>() }),
        TruthAllele::Repeat {
            locus_id,
            bases,
            blocks,
        } => json!({"kind":"repeat","locus_id":locus_id,"bases":bases,
            "blocks":blocks.iter().map(repeat_block_json).collect::<Vec<_>>() }),
    }
}

/// Un bloque publicado: su motivo, su longitud, sus bases y **su número de
/// unidades derivado de las dos**, que es el campo que el VCF 4.5 pide como
/// decimal porque un bloque puede terminar en una unidad parcial.
fn repeat_block_json(block: &TruthRepeatBlock) -> Value {
    let motif_length = block.motif.len() as u64;
    json!({
        "motif": block.motif,
        "motif_length": motif_length,
        "bases": block.bases,
        "units": block.bases as f64 / motif_length as f64,
    })
}
fn structural_event_json(event: &TruthStructuralEvent) -> Value {
    match event {
        TruthStructuralEvent::Insertion {
            event_id,
            point,
            bases,
        } => json!({"kind":"insertion","event_id":event_id,"point":point,"bases":bases}),
        TruthStructuralEvent::Deletion {
            event_id,
            start,
            end,
        } => json!({"kind":"deletion","event_id":event_id,"start":start,"end":end}),
    }
}
fn piece_link_json(link: &TruthPieceLink) -> Value {
    json!({"piece_index":link.piece_index,"start":link.start,"end":link.end})
}

fn current_piece_json(piece: &TruthCurrentPiece) -> Value {
    match piece {
        TruthCurrentPiece::Source { source, snps } => json!({
            "kind":"source",
            "source":{"source":source_json(&source.source),"reference_id":source.reference_id,"sequence_id":source.sequence_id,"start":source.start,"end":source.end,"sequence_length":source.sequence_length},
            "snps":snps.iter().map(|snp| json!({"event_id":snp.event_id,"position":snp.position,"reference":snp.reference,"alternate":snp.alternate,"provenance":provenance_json(&snp.provenance)})).collect::<Vec<_>>()
        }),
        TruthCurrentPiece::Inserted { bases, provenance } => {
            json!({"kind":"inserted","bases":bases,"provenance":provenance_json(provenance)})
        }
        TruthCurrentPiece::RepeatAllele {
            locus_id,
            source,
            bases,
            blocks,
            provenance,
        } => json!({
            "kind":"repeat_allele","locus_id":locus_id,
            "source":{"source":source_json(&source.source),"reference_id":source.reference_id,"sequence_id":source.sequence_id,"start":source.start,"end":source.end,"sequence_length":source.sequence_length},
            "bases":bases,
            "blocks":blocks.iter().map(repeat_block_json).collect::<Vec<_>>(),
            "provenance":provenance_json(provenance)
        }),
    }
}

fn history_event_json(event: &TruthHistoryEvent) -> Value {
    match event {
        TruthHistoryEvent::Snp {
            event_id,
            position,
            provenance,
        } => {
            json!({"kind":"snp","event_id":event_id,"position":position,"provenance":provenance_json(provenance)})
        }
        TruthHistoryEvent::Insertion {
            event_id,
            point,
            bases,
            provenance,
        } => {
            json!({"kind":"insertion","event_id":event_id,"point":point,"bases":bases,"provenance":provenance_json(provenance)})
        }
        TruthHistoryEvent::Deletion {
            event_id,
            start,
            end,
            provenance,
        } => {
            json!({"kind":"deletion","event_id":event_id,"start":start,"end":end,"provenance":provenance_json(provenance)})
        }
        TruthHistoryEvent::RepeatAllele {
            event_id,
            locus_id,
            start,
            end,
            bases,
            blocks,
            provenance,
        } => {
            json!({"kind":"repeat_allele","event_id":event_id,"locus_id":locus_id,"start":start,"end":end,
                "bases":bases,"blocks":blocks.iter().map(repeat_block_json).collect::<Vec<_>>(),
                "provenance":provenance_json(provenance)})
        }
    }
}

fn modification_json(modification: &TruthModification) -> Value {
    match modification {
        TruthModification::Snp {
            event_id,
            position,
            reference,
            alternate,
            provenance,
        } => {
            json!({"kind":"snp","event_id":event_id,"position":position,"reference":reference,"alternate":alternate,"provenance":provenance_json(provenance)})
        }
        TruthModification::Insertion {
            event_id,
            point,
            bases,
            provenance,
        } => {
            json!({"kind":"insertion","event_id":event_id,"point":point,"bases":bases,"provenance":provenance_json(provenance)})
        }
        TruthModification::Deletion {
            event_id,
            start,
            end,
            provenance,
        } => {
            json!({"kind":"deletion","event_id":event_id,"start":start,"end":end,"provenance":provenance_json(provenance)})
        }
    }
}

fn provenance_json(provenance: &ModificationProvenance) -> Value {
    let kind = match provenance.kind {
        ModificationProvenanceKind::Explicit => "explicit",
        ModificationProvenanceKind::Random => "random",
        ModificationProvenanceKind::Matrix => "matrix",
    };
    let mut value =
        json!({"kind":kind,"version":provenance.version,"source":source_json(&provenance.source)});
    if let Some(matrix) = &provenance.matrix {
        value["matrix"] = json!({"matrix_source_id":matrix.matrix_source_id,"population_id":matrix.population_id,"variant_id":matrix.variant_id});
    }
    value
}

fn digest_json(digest: &NamedDigest) -> Value {
    json!({"source":source_json(&digest.source),"sha256":digest.sha256})
}
fn source_json(source: &TruthSourceRef) -> Value {
    json!({"kind":source_kind_text(source.kind),"source_id":source.source_id})
}
fn copy_identity_json(identity: &CopyIdentity) -> Value {
    json!({"reference_id":identity.reference_id,"individual_id":identity.individual_id,"sequence_id":identity.sequence_id,"copy_index":identity.copy_index})
}
fn source_kind_text(kind: TruthSourceKind) -> &'static str {
    match kind {
        TruthSourceKind::Reference => "reference",
        TruthSourceKind::Founder => "founder",
        TruthSourceKind::Explicit => "explicit",
        TruthSourceKind::Random => "random",
        TruthSourceKind::Matrix => "matrix",
    }
}
fn base_text(base: CanonicalBase) -> String {
    char::from(base.as_byte()).to_string()
}
fn bytes_to_truth(event_id: &str, bases: &[u8]) -> Result<String, TruthProjectionError> {
    std::str::from_utf8(bases).map(str::to_owned).map_err(|_| {
        TruthProjectionError::InvalidUtf8Insertion {
            event_id: event_id.to_owned(),
        }
    })
}
fn reference_source(reference: &Reference) -> TruthSourceRef {
    reference_source_from_id(reference.id().as_str())
}
/// Identidad del recurso de matriz preparado que realmente consumió F2.3.
/// No reutiliza los `source_id` de filas, que pueden mezclar procedencias de
/// recursos suministrados y reutilizados.
fn matrix_artifact_source(matrix: &PreparedPopulationSnvs) -> TruthSourceRef {
    TruthSourceRef {
        kind: TruthSourceKind::Matrix,
        source_id: matrix.canonical_fingerprint(),
    }
}
fn matrix_source_id(source_id: &str) -> &str {
    if source_id.is_empty() {
        "MATRIX_DEFAULT"
    } else {
        source_id
    }
}
fn reference_source_from_id(id: &str) -> TruthSourceRef {
    TruthSourceRef {
        kind: TruthSourceKind::Reference,
        source_id: id.to_owned(),
    }
}
fn domain_source(provenance: &crate::segments::ModificationProvenance) -> TruthSourceRef {
    TruthSourceRef {
        kind: match provenance.origin() {
            ModificationOrigin::Explicit => TruthSourceKind::Explicit,
            ModificationOrigin::Random => TruthSourceKind::Random,
        },
        source_id: provenance.event_id().as_str().to_owned(),
    }
}
fn domain_provenance(
    provenance: &crate::segments::ModificationProvenance,
    source: TruthSourceRef,
) -> ModificationProvenance {
    ModificationProvenance {
        kind: match provenance.origin() {
            ModificationOrigin::Explicit => ModificationProvenanceKind::Explicit,
            ModificationOrigin::Random => ModificationProvenanceKind::Random,
        },
        version: "1".to_owned(),
        source,
        matrix: None,
    }
}
fn modification_source(modification: &TruthModification) -> &TruthSourceRef {
    match modification {
        TruthModification::Snp { provenance, .. }
        | TruthModification::Insertion { provenance, .. }
        | TruthModification::Deletion { provenance, .. } => &provenance.source,
    }
}
fn copy_identity(copy: &HaplotypeCopy) -> CopyIdentity {
    CopyIdentity {
        reference_id: copy.id().reference_id().as_str().to_owned(),
        individual_id: copy.id().individual_id().as_str().to_owned(),
        sequence_id: copy.id().sequence_id().as_str().to_owned(),
        copy_index: copy.id().copy_index(),
    }
}
fn require_digest(
    digests: &[NamedDigest],
    source: &TruthSourceRef,
) -> Result<(), TruthProjectionError> {
    if digests.iter().any(|digest| &digest.source == source) {
        Ok(())
    } else {
        Err(TruthProjectionError::MissingDigest(source.clone()))
    }
}
fn sorted_digests(digests: &[NamedDigest]) -> Vec<NamedDigest> {
    let mut sorted = digests.to_vec();
    sorted.sort_by(|left, right| left.source.cmp(&right.source));
    sorted
}
fn record_order(left: &TruthRecord, right: &TruthRecord) -> std::cmp::Ordering {
    match (left, right) {
        (TruthRecord::HaplotypeCopy(left), TruthRecord::HaplotypeCopy(right)) => {
            left.identity.cmp(&right.identity)
        }
        (TruthRecord::HaplotypeCopy(_), TruthRecord::Genotype(_)) => std::cmp::Ordering::Less,
        (TruthRecord::Genotype(_), TruthRecord::HaplotypeCopy(_)) => std::cmp::Ordering::Greater,
        (TruthRecord::Genotype(left), TruthRecord::Genotype(right)) => (
            &left.reference_id,
            &left.individual_id,
            &left.locus,
            &left.case_id,
        )
            .cmp(&(
                &right.reference_id,
                &right.individual_id,
                &right.locus,
                &right.case_id,
            )),
    }
}
pub(crate) fn case_id(
    population: &str,
    individual: &str,
    sequence: &str,
    position_1based: u64,
    alleles: &[TruthAllele],
) -> String {
    let mut hasher = Sha256::new();
    for value in [
        population,
        individual,
        sequence,
        &position_1based.to_string(),
        &alleles.iter().map(allele_key).collect::<Vec<_>>().join(","),
    ] {
        hasher.update(value.as_bytes());
        hasher.update([0]);
    }
    let digest = format!("{:x}", hasher.finalize());
    format!("CASE_{}", &digest[..59])
}
fn allele_key(allele: &TruthAllele) -> String {
    match allele {
        TruthAllele::Reference => "REF".to_owned(),
        TruthAllele::Snp { base } => format!("SNP:{base}"),
        TruthAllele::Insertion {
            event_id, bases, ..
        } => format!("INS:{event_id}:{bases}"),
        TruthAllele::Deletion {
            event_id,
            start,
            end,
            ..
        } => format!("DEL:{event_id}:{start}:{end}"),
        // La clave lleva **las bases**, no el número de unidades: dos alelos de
        // igual longitud con la interrupción en sitio distinto son alelos
        // distintos y no pueden compartir identificador de caso.
        TruthAllele::Repeat {
            locus_id, bases, ..
        } => format!("STR:{locus_id}:{bases}"),
    }
}
fn contract_diagnostic(error: TruthContractError) -> Diagnostic {
    Diagnostic::domain(
        "GPS-DOMAIN-003",
        "La proyección de verdad no cumple el contrato externo v1.",
    )
    .with("component", "truth_projection")
    .with("detail", format!("{error:?}"))
}

#[cfg(test)]
mod tests {
    use std::fs;
    use std::io::{self, Cursor, Write};
    use std::sync::atomic::{AtomicU64, Ordering};

    use sha2::{Digest, Sha256};

    use super::*;
    use crate::adapters::{admit_toml_header, parse_declared_config};
    use crate::application::build_effective_config;
    use crate::configuration::ExpectedConfigKind;
    use crate::coordinates::BasePosition;
    use crate::haplotype::HaplotypePhase;
    use crate::identity::{EventId, ReferenceId, SequenceId};
    use crate::initial_individual_integration::IncompatibleGenotypeAuthorization;
    use crate::initial_individual_materialization::materialize_initial_individuals;
    use crate::materialization::SnpLocus;
    use crate::population_snvs::{
        TSV_HEADER, parse_canonical_tsv, prepare_population_snvs_from_reference,
    };
    use crate::reference::{
        DnaSequence, FastaHeader, ReferenceOrigin, ReferenceSequence, SequenceClass,
    };
    use crate::segments::{CanonicalBase, ModificationProvenance as DomainModificationProvenance};
    use crate::small_event_composition::{
        DeclaredSmallEvent, TargetedCurrentDeletion, TargetedCurrentInsertion,
    };
    use crate::snp_application::{ApplicableSnp, TargetedSnp};
    use crate::truth_validation::{
        TruthValidationError, reference_content_sha256, validate_and_reconstruct_truth,
        validate_and_reconstruct_truth_jsonl, validate_truth_jsonl_streaming,
    };

    static NONCE: AtomicU64 = AtomicU64::new(0);

    fn document() -> TruthDocument {
        let source = TruthSourceRef {
            kind: TruthSourceKind::Reference,
            source_id: "REF_main".to_owned(),
        };
        TruthDocument {
            header: TruthHeader {
                document_kind: TRUTH_DOCUMENT_KIND.to_owned(),
                schema_version: TRUTH_SCHEMA_VERSION.to_owned(),
                execution_id: "RUN_01".to_owned(),
                reference_id: "REF_main".to_owned(),
                reference_content_sha256: "b".repeat(64),
                coordinate_convention: TRUTH_COORDINATE_CONVENTION.to_owned(),
                copy_count: 1,
                genotype_count: 0,
                source_digests: vec![NamedDigest {
                    source: source.clone(),
                    sha256: "a".repeat(64),
                }],
                group_digests: vec![],
            },
            records: vec![TruthRecord::HaplotypeCopy(HaplotypeCopyRecord {
                identity: CopyIdentity {
                    reference_id: "REF_main".to_owned(),
                    individual_id: "IND_01".to_owned(),
                    sequence_id: "chr1".to_owned(),
                    copy_index: 0,
                },
                phase: 0,
                provenance: CopyProvenance {
                    origin: "initial_population:P1".to_owned(),
                    sources: vec![source.clone()],
                },
                segments: vec![TruthSegment {
                    source: SegmentSource {
                        source,
                        reference_id: "REF_main".to_owned(),
                        sequence_id: "chr1".to_owned(),
                        start: 0,
                        end: 4,
                        sequence_length: 4,
                    },
                    modifications: vec![],
                }],
                current_pieces: vec![TruthCurrentPiece::Source {
                    source: SegmentSource {
                        source: TruthSourceRef {
                            kind: TruthSourceKind::Reference,
                            source_id: "REF_main".to_owned(),
                        },
                        reference_id: "REF_main".to_owned(),
                        sequence_id: "chr1".to_owned(),
                        start: 0,
                        end: 4,
                        sequence_length: 4,
                    },
                    snps: vec![],
                }],
                history: vec![],
            })],
        }
    }

    fn reference() -> Reference {
        Reference::new(
            ReferenceId::parse("REF_main", "test").expect("reference ID"),
            ReferenceOrigin::Generated,
            vec![
                ReferenceSequence::new(
                    SequenceId::parse("chr1", "test").expect("sequence ID"),
                    FastaHeader::parse(b"chr1").expect("header"),
                    SequenceClass::Chromosome,
                    &DnaSequence::try_from_bytes(b"ACGT").expect("DNA"),
                )
                .expect("reference sequence"),
            ],
        )
        .expect("reference")
    }

    fn f2_3_state(
        reference: &Reference,
        events: Vec<DeclaredSmallEvent>,
    ) -> (InitialIndividualsMaterialization, PreparedPopulationSnvs) {
        f2_3_state_with(
            reference,
            2,
            &["VAR_01\tchr1\t2\tC\tT\tPOP_01\t1\tMATRIX_01"],
            events,
        )
    }

    fn f2_3_state_with(
        reference: &Reference,
        ploidy: u8,
        rows: &[&str],
        events: Vec<DeclaredSmallEvent>,
    ) -> (InitialIndividualsMaterialization, PreparedPopulationSnvs) {
        let config = format!(
            r#"config_kind = "general"
schema_version = "2"
seed = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f"
[reference]
kind = "simulated"
id = "REF_main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 4
[[initial_populations]]
population_id = "POP_01"
individual_count = 1
individual_ids = ["IND_01"]
ploidy = {ploidy}
genotype_model = "hwe""#
        );
        let admitted = admit_toml_header(
            Path::new("truth-projection.toml"),
            config.as_bytes().to_vec(),
            ExpectedConfigKind::General,
        )
        .expect("header");
        let effective = build_effective_config(&parse_declared_config(admitted).expect("declared"))
            .expect("effective");
        let bytes = format!("{}\n{}", TSV_HEADER.join("\t"), rows.join("\n"));
        let matrix = prepare_population_snvs_from_reference(
            &parse_canonical_tsv(bytes.as_bytes()).expect("TSV"),
            reference,
            1e-9,
            true,
        )
        .expect("prepared matrix");
        let state = materialize_initial_individuals(
            &effective,
            reference,
            &matrix,
            IncompatibleGenotypeAuthorization::denied(),
            move |_| Ok(events.clone()),
        )
        .expect("F2.3 state");
        (state, matrix)
    }

    fn required_source_digests(
        reference: &Reference,
        matrix: &PreparedPopulationSnvs,
        additional: impl IntoIterator<Item = TruthSourceRef>,
    ) -> Vec<NamedDigest> {
        let mut sources = vec![reference_source(reference), matrix_artifact_source(matrix)];
        sources.extend(matrix.rows().iter().map(|row| TruthSourceRef {
            kind: TruthSourceKind::Matrix,
            source_id: matrix_source_id(row.source_id()).to_owned(),
        }));
        sources.extend(additional);
        sources
            .into_iter()
            .enumerate()
            .map(|(index, source)| NamedDigest {
                source,
                sha256: format!("{index:064x}"),
            })
            .collect()
    }

    fn required_group_digests(matrix: &PreparedPopulationSnvs) -> Vec<NamedDigest> {
        matrix
            .population_ids()
            .iter()
            .enumerate()
            .map(|(index, population_id)| NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Matrix,
                    source_id: population_id.as_str().to_owned(),
                },
                sha256: format!("{index:064x}"),
            })
            .collect()
    }

    #[test]
    fn projection_from_f2_3_state_is_ordered_and_counts_copies() {
        let reference = reference();
        let (state, matrix) = f2_3_state(&reference, Vec::new());
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_01",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &[
                NamedDigest {
                    source: reference_source(&reference),
                    sha256: "a".repeat(64),
                },
                NamedDigest {
                    source: matrix_artifact_source(&matrix),
                    sha256: "b".repeat(64),
                },
                NamedDigest {
                    source: TruthSourceRef {
                        kind: TruthSourceKind::Matrix,
                        source_id: "MATRIX_01".to_owned(),
                    },
                    sha256: "c".repeat(64),
                },
            ],
            group_digests: &[NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Matrix,
                    source_id: "POP_01".to_owned(),
                },
                sha256: "d".repeat(64),
            }],
        })
        .expect("projection");
        assert_eq!(document.header.copy_count, 2);
        assert_eq!(document.header.genotype_count, 1);
        assert!(
            document
                .records
                .windows(2)
                .all(|pair| record_order(&pair[0], &pair[1]).is_lt())
        );
    }

    #[test]
    fn projection_retains_nested_partial_indel_links_and_final_bytes() {
        let reference = reference();
        let phase = HaplotypePhase::from_copy_index(0).expect("phase");
        let sequence = SequenceId::parse("chr1", "test").expect("sequence");
        let insertion = |id: &str, point, bases: &[u8]| {
            DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
                sequence.clone(),
                phase,
                point,
                std::sync::Arc::<[u8]>::from(bases),
                DomainModificationProvenance::new(EventId::parse(id, "test").expect("event")),
            ))
        };
        let deletion = DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
            sequence.clone(),
            phase,
            1,
            5,
            DomainModificationProvenance::new(EventId::parse("DEL_01", "test").expect("event")),
        ));
        let (state, matrix) = f2_3_state(
            &reference,
            vec![
                insertion("INS_01", 1, b"ACGT"),
                insertion("INS_02", 2, b"T"),
                deletion,
            ],
        );
        let matrix_artifact = matrix_artifact_source(&matrix);
        let digests = vec![
            NamedDigest {
                source: reference_source(&reference),
                sha256: "a".repeat(64),
            },
            NamedDigest {
                source: matrix_artifact,
                sha256: "b".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Matrix,
                    source_id: "MATRIX_01".into(),
                },
                sha256: "c".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "INS_01".into(),
                },
                sha256: "d".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "INS_02".into(),
                },
                sha256: "e".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "DEL_01".into(),
                },
                sha256: "f".repeat(64),
            },
        ];
        let groups = [NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: "POP_01".into(),
            },
            sha256: "0".repeat(64),
        }];
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_02",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &digests,
            group_digests: &groups,
        })
        .expect("projection");
        let insertion = |event_id: &str| {
            document
                .records
                .iter()
                .find_map(|record| match record {
                    TruthRecord::Genotype(GenotypeRecord {
                        case_kind: GenotypeCaseKind::Indel,
                        alleles,
                        ..
                    }) => alleles.iter().find_map(|allele| match allele {
                        TruthAllele::Insertion {
                            event_id: observed,
                            source_event,
                            piece_links,
                            ..
                        } if observed == event_id => Some((source_event, piece_links)),
                        _ => None,
                    }),
                    _ => None,
                })
                .expect("insertion case")
        };
        let (source_event, links) = insertion("INS_01");
        assert!(
            matches!(source_event, TruthStructuralEvent::Insertion { event_id, point: 1, bases } if event_id == "INS_01" && bases == "ACGT")
        );
        assert_eq!(
            links,
            &[TruthPieceLink {
                piece_index: 1,
                start: 1,
                end: 2
            }]
        );
        assert!(document.records.iter().all(|record| !matches!(record,
            TruthRecord::Genotype(GenotypeRecord { alleles, .. })
            if alleles.iter().any(|allele| matches!(allele, TruthAllele::Insertion { event_id, .. } if event_id == "INS_02"))
        )));
        assert!(document.records.iter().any(|record| matches!(record,
            TruthRecord::HaplotypeCopy(copy)
            if copy.history.iter().any(|event| matches!(event, TruthHistoryEvent::Insertion { event_id, .. } if event_id == "INS_02"))
        )));
        // Los bytes finales ya no viajan en el alelo: se exigen a la copia que
        // el validador independiente reconstruye desde la verdad, y que sean
        // también los que el productor materializa desde su estado.
        let reconstructed =
            validate_and_reconstruct_truth(document, &reference).expect("independent replay");
        let (identity, bases) = reconstructed
            .copies
            .iter()
            .find(|(identity, _)| identity.copy_index == 0)
            .expect("copy 0");
        assert_eq!(bases.as_slice(), b"ATTGT");
        let simulated = state
            .individuals()
            .iter()
            .flat_map(|individual| individual.individual().sequences())
            .flat_map(|sequence| sequence.copies())
            .find(|copy| &copy_identity(copy) == identity)
            .expect("simulated copy 0");
        assert_eq!(
            materialize_copy(simulated, &reference).expect("materialized copy 0"),
            b"ATTGT"
        );
    }

    /// **La verdad publicada tiene que ser la copia simulada, copia a copia.**
    ///
    /// El contraste compara la reconstrucción independiente de cada registro
    /// con la materialización del estado interno. Se prueba con dos verdades
    /// que se reconstruyen sin error pero no describen el estado:
    ///
    /// - la de un estado cuya inserción lleva otras bases, contrastada con el
    ///   estado verdadero: falla la copia 0, la del *indel*;
    /// - la verdadera con el SNP de la **copia 1**, que no tiene ningún
    ///   *indel*, cambiado de forma coherente en sus segmentos y en sus
    ///   fragmentos vigentes: la copia 0 pasa y falla la 1. Así se ve que el
    ///   contraste no se limita a las copias con *indels*.
    ///
    /// Lo que esta prueba **no** puede ver es que `project_truth` deje de
    /// llamar a la verificación: sin una proyección defectuosa no hay
    /// discrepancia que cazar. Eso lo cazan, a través del gancho
    /// `test_tamper` y de la ejecución completa, las pruebas
    /// `shifted_case_locus_*` y `changed_case_dosage_*` de `execution.rs`.
    #[test]
    fn publication_rejects_a_truth_that_does_not_describe_the_simulated_copies() {
        let reference = reference();
        let phase = HaplotypePhase::from_copy_index(0).expect("phase");
        let sequence = SequenceId::parse("chr1", "test").expect("sequence");
        let insertion = |bases: &[u8]| {
            DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
                sequence.clone(),
                phase,
                1,
                std::sync::Arc::<[u8]>::from(bases),
                DomainModificationProvenance::new(EventId::parse("INS_01", "test").expect("event")),
            ))
        };
        let (state, matrix) = f2_3_state(&reference, vec![insertion(b"ACGT")]);
        let (other_state, _) = f2_3_state(&reference, vec![insertion(b"ACGG")]);
        let sources = required_source_digests(
            &reference,
            &matrix,
            [TruthSourceRef {
                kind: TruthSourceKind::Explicit,
                source_id: "INS_01".to_owned(),
            }],
        );
        let groups = required_group_digests(&matrix);
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_AGREEMENT",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &sources,
            group_digests: &groups,
        })
        .expect("the true projection agrees with its own state");

        assert!(matches!(
            verify_publication_against_simulation(&document, &reference, &other_state),
            Err(TruthProjectionError::ReconstructionMismatch { identity }) if identity.copy_index == 0
        ));

        let mut tampered = document.clone();
        let mut changed = 0;
        for record in &mut tampered.records {
            let TruthRecord::HaplotypeCopy(copy) = record else {
                continue;
            };
            if copy.identity.copy_index != 1 {
                continue;
            }
            assert!(
                copy.history
                    .iter()
                    .all(|event| matches!(event, TruthHistoryEvent::Snp { .. })),
                "copy 1 carries no indel"
            );
            for segment in &mut copy.segments {
                for modification in &mut segment.modifications {
                    if let TruthModification::Snp { alternate, .. } = modification {
                        *alternate = "G".to_owned();
                        changed += 1;
                    }
                }
            }
            for piece in &mut copy.current_pieces {
                if let TruthCurrentPiece::Source { snps, .. } = piece {
                    for snp in snps {
                        snp.alternate = "G".to_owned();
                        changed += 1;
                    }
                }
            }
        }
        assert_eq!(changed, 2, "the SNP of copy 1 must change in both views");
        assert!(matches!(
            verify_publication_against_simulation(&tampered, &reference, &state),
            Err(TruthProjectionError::ReconstructionMismatch { identity }) if identity.copy_index == 1
        ));
    }

    #[test]
    fn f2_3_projection_writer_and_independent_readers_replay_all_small_event_classes() {
        let reference = reference();
        let sequence = SequenceId::parse("chr1", "test").expect("sequence");
        let phase = HaplotypePhase::from_copy_index(0).expect("phase");
        let insertion = |id: &str, point, bases: &[u8]| {
            DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
                sequence.clone(),
                phase,
                point,
                std::sync::Arc::<[u8]>::from(bases),
                DomainModificationProvenance::new(EventId::parse(id, "test").expect("event")),
            ))
        };
        let snp = DeclaredSmallEvent::Snp(TargetedSnp::new(
            phase,
            ApplicableSnp::explicit(
                SnpLocus::new(
                    BasePosition::new(reference.id().clone(), sequence.clone(), 2, 4)
                        .expect("position"),
                    CanonicalBase::from_byte(b'G').expect("reference base"),
                    CanonicalBase::from_byte(b'A').expect("alternate base"),
                )
                .expect("SNP locus"),
                EventId::parse("SNP_DIRECT", "test").expect("event"),
            ),
        ));
        let deletion = DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
            sequence.clone(),
            phase,
            1,
            5,
            DomainModificationProvenance::new(EventId::parse("DEL_01", "test").expect("event")),
        ));
        let (state, matrix) = f2_3_state(
            &reference,
            vec![
                snp,
                insertion("INS_01", 1, b"ACGT"),
                insertion("INS_02", 2, b"T"),
                deletion,
            ],
        );
        let sources = required_source_digests(
            &reference,
            &matrix,
            ["SNP_DIRECT", "INS_01", "INS_02", "DEL_01"].map(|source_id| TruthSourceRef {
                kind: TruthSourceKind::Explicit,
                source_id: source_id.to_owned(),
            }),
        );
        let groups = required_group_digests(&matrix);
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_END_TO_END",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &sources,
            group_digests: &groups,
        })
        .expect("F2.3 projection");
        assert_eq!(
            document.header.reference_content_sha256,
            reference_content_sha256(&reference).expect("huella")
        );
        let mut bytes = Vec::new();
        write_truth_jsonl(&document, &mut bytes).expect("truth writer");
        let reconstructed =
            validate_and_reconstruct_truth_jsonl(&mut Cursor::new(&bytes), &reference)
                .expect("independent full validation");
        assert_eq!(reconstructed.document, document);
        let summary =
            validate_truth_jsonl_streaming(&mut Cursor::new(bytes), &reference, |_, _| {})
                .expect("independent streaming validation");
        assert_eq!(summary.copy_count, document.header.copy_count);
        assert_eq!(summary.genotype_count, document.header.genotype_count);

        let mut duplicated = document.clone();
        let TruthRecord::HaplotypeCopy(copy) = duplicated
            .records
            .iter_mut()
            .find(|record| matches!(record, TruthRecord::HaplotypeCopy(copy) if copy.identity.copy_index == 0))
            .expect("target copy")
        else {
            unreachable!("the target record is a copy");
        };
        let inherited_snp_id = copy
            .segments
            .iter()
            .flat_map(|segment| segment.modifications.iter())
            .find_map(|modification| match modification {
                TruthModification::Snp { event_id, .. } if event_id == "SNP_DIRECT" => {
                    Some(event_id.clone())
                }
                _ => None,
            })
            .expect("inherited SNP");
        copy.history.retain(|event| {
            !matches!(event, TruthHistoryEvent::Snp { event_id, .. } if event_id == &inherited_snp_id)
        });
        let TruthHistoryEvent::Insertion { event_id, .. } = copy
            .history
            .iter_mut()
            .find(|event| matches!(event, TruthHistoryEvent::Insertion { event_id, .. } if event_id == "INS_01"))
            .expect("active insertion")
        else {
            unreachable!("the selected event is an insertion");
        };
        *event_id = inherited_snp_id;
        duplicated
            .validate_structure()
            .expect("the cross-class collision is structurally representable");
        let mut duplicate_bytes = Vec::new();
        write_truth_jsonl(&duplicated, &mut duplicate_bytes).expect("mutated truth writer");
        assert!(matches!(
            validate_and_reconstruct_truth_jsonl(&mut Cursor::new(&duplicate_bytes), &reference),
            Err(TruthValidationError::HistoryDoesNotReplay { .. })
        ));
        assert!(matches!(
            validate_truth_jsonl_streaming(
                &mut Cursor::new(duplicate_bytes),
                &reference,
                |_, _| {}
            ),
            Err(TruthValidationError::HistoryDoesNotReplay { .. })
        ));
    }

    #[test]
    fn projection_requires_reference_matrix_row_and_group_digests() {
        let reference = reference();
        let (state, matrix) = f2_3_state(&reference, Vec::new());
        let reference_digest = NamedDigest {
            source: reference_source(&reference),
            sha256: "a".repeat(64),
        };
        let matrix_digest = NamedDigest {
            source: matrix_artifact_source(&matrix),
            sha256: "b".repeat(64),
        };
        let row_digest = NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: "MATRIX_01".into(),
            },
            sha256: "c".repeat(64),
        };
        let group_digest = NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: "POP_01".into(),
            },
            sha256: "d".repeat(64),
        };
        let project = |sources: &[NamedDigest], groups: &[NamedDigest]| {
            project_truth(TruthProjectionInput {
                execution_id: "RUN_03",
                reference: &reference,
                individuals: &state,
                matrix: &matrix,
                source_digests: sources,
                group_digests: groups,
            })
        };
        assert!(
            matches!(project(&[matrix_digest.clone(), row_digest.clone()], std::slice::from_ref(&group_digest)), Err(TruthProjectionError::MissingDigest(source)) if source == reference_digest.source)
        );
        assert!(
            matches!(project(&[reference_digest.clone(), row_digest.clone()], std::slice::from_ref(&group_digest)), Err(TruthProjectionError::MissingDigest(source)) if source == matrix_digest.source)
        );
        assert!(
            matches!(project(&[reference_digest.clone(), matrix_digest.clone()], std::slice::from_ref(&group_digest)), Err(TruthProjectionError::MissingDigest(source)) if source == row_digest.source)
        );
        assert!(
            matches!(project(&[reference_digest, matrix_digest, row_digest], &[]), Err(TruthProjectionError::MissingDigest(source)) if source == group_digest.source)
        );
    }

    #[test]
    fn projection_preserves_random_indel_provenance() {
        let reference = reference();
        let event = DeclaredSmallEvent::Insertion(TargetedCurrentInsertion::new(
            SequenceId::parse("chr1", "test").expect("sequence"),
            HaplotypePhase::from_copy_index(0).expect("phase"),
            1,
            std::sync::Arc::<[u8]>::from(&b"A"[..]),
            DomainModificationProvenance::random(EventId::parse("RAND_01", "test").expect("event")),
        ));
        let (state, matrix) = f2_3_state(&reference, vec![event]);
        let digests = vec![
            NamedDigest {
                source: reference_source(&reference),
                sha256: "a".repeat(64),
            },
            NamedDigest {
                source: matrix_artifact_source(&matrix),
                sha256: "b".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Matrix,
                    source_id: "MATRIX_01".into(),
                },
                sha256: "c".repeat(64),
            },
            NamedDigest {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Random,
                    source_id: "RAND_01".into(),
                },
                sha256: "d".repeat(64),
            },
        ];
        let groups = [NamedDigest {
            source: TruthSourceRef {
                kind: TruthSourceKind::Matrix,
                source_id: "POP_01".into(),
            },
            sha256: "e".repeat(64),
        }];
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_04",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &digests,
            group_digests: &groups,
        })
        .expect("projection");
        assert!(document.records.iter().any(|record| matches!(record,
            TruthRecord::HaplotypeCopy(copy)
            if copy.history.iter().any(|event| matches!(event, TruthHistoryEvent::Insertion { provenance, .. } if provenance.kind == ModificationProvenanceKind::Random))
        )));
    }

    #[test]
    fn projection_direct_snp_and_deletion_cases_cover_p1_p2_and_p10() {
        let reference = reference();
        let sequence = SequenceId::parse("chr1", "test").expect("sequence");
        for ploidy in [1_u8, 2, 10] {
            let phase = HaplotypePhase::from_copy_index(0).expect("phase");
            let direct = DeclaredSmallEvent::Snp(TargetedSnp::new(
                phase,
                ApplicableSnp::explicit(
                    SnpLocus::new(
                        BasePosition::new(reference.id().clone(), sequence.clone(), 2, 4)
                            .expect("position"),
                        CanonicalBase::from_byte(b'G').expect("reference base"),
                        CanonicalBase::from_byte(b'A').expect("alternate base"),
                    )
                    .expect("SNP locus"),
                    EventId::parse("SNP_DIRECT", "test").expect("event"),
                ),
            ));
            let (state, matrix) = f2_3_state_with(
                &reference,
                ploidy,
                &["VAR_01\tchr1\t2\tC\tT\tPOP_01\t1\tMATRIX_01"],
                vec![direct],
            );
            let sources = required_source_digests(
                &reference,
                &matrix,
                [TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "SNP_DIRECT".to_owned(),
                }],
            );
            let groups = required_group_digests(&matrix);
            let document = project_truth(TruthProjectionInput {
                execution_id: "RUN_DIRECT",
                reference: &reference,
                individuals: &state,
                matrix: &matrix,
                source_digests: &sources,
                group_digests: &groups,
            })
            .expect("direct SNP projection");
            let direct_case = document
                .records
                .iter()
                .find_map(|record| match record {
                    TruthRecord::Genotype(genotype)
                        if genotype.case_kind == GenotypeCaseKind::Snp
                            && genotype.locus.sequence_id == "chr1"
                            && genotype.locus.start == 2 =>
                    {
                        Some(genotype)
                    }
                    _ => None,
                })
                .expect("direct SNP case");
            assert_eq!(
                direct_case.alleles,
                vec![
                    TruthAllele::Reference,
                    TruthAllele::Snp { base: "A".into() }
                ]
            );
            assert_eq!(direct_case.carriers.len(), usize::from(ploidy));
            assert_eq!(direct_case.allele_dosages, vec![ploidy - 1, 1]);

            let deletion = DeclaredSmallEvent::Deletion(TargetedCurrentDeletion::new(
                sequence.clone(),
                phase,
                1,
                2,
                DomainModificationProvenance::new(
                    EventId::parse("DEL_PLOIDY", "test").expect("event"),
                ),
            ));
            let (state, matrix) = f2_3_state_with(
                &reference,
                ploidy,
                &["VAR_01\tchr1\t2\tC\tT\tPOP_01\t1\tMATRIX_01"],
                vec![deletion],
            );
            let sources = required_source_digests(
                &reference,
                &matrix,
                [TruthSourceRef {
                    kind: TruthSourceKind::Explicit,
                    source_id: "DEL_PLOIDY".to_owned(),
                }],
            );
            let groups = required_group_digests(&matrix);
            let document = project_truth(TruthProjectionInput {
                execution_id: "RUN_DELETION",
                reference: &reference,
                individuals: &state,
                matrix: &matrix,
                source_digests: &sources,
                group_digests: &groups,
            })
            .expect("deletion projection");
            let deletion_case = document
                .records
                .iter()
                .find_map(|record| {
                    match record {
                    TruthRecord::Genotype(genotype)
                        if genotype.case_kind == GenotypeCaseKind::Indel
                            && genotype.alleles.iter().any(|allele| matches!(
                                allele,
                                TruthAllele::Deletion { event_id, .. } if event_id == "DEL_PLOIDY"
                            )) =>
                    {
                        Some(genotype)
                    }
                    _ => None,
                }
                })
                .expect("deletion structural case");
            assert_eq!(deletion_case.carriers.len(), usize::from(ploidy));
            assert_eq!(deletion_case.allele_dosages, vec![ploidy - 1, 1]);
            assert!(matches!(
                deletion_case.alleles.as_slice(),
                [
                    TruthAllele::Reference,
                    TruthAllele::Deletion {
                        start: 1,
                        end: 2,
                        ..
                    }
                ]
            ));
        }
    }

    #[test]
    fn projection_accepts_mixed_row_sources_and_retains_catalog_alleles() {
        let reference = reference();
        let (state, matrix) = f2_3_state_with(
            &reference,
            2,
            &[
                "VAR_SUPPLIED\tchr1\t2\tC\tA\tPOP_01\t0.25\tSUPPLIED_01",
                "VAR_REUSED\tchr1\t2\tC\tT\tPOP_01\t0.25\tREUSED_01",
            ],
            Vec::new(),
        );
        let sources = required_source_digests(&reference, &matrix, []);
        let groups = required_group_digests(&matrix);
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_MIXED",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &sources,
            group_digests: &groups,
        })
        .expect("mixed supplied/reused matrix rows");
        assert!(document.header.source_digests.iter().any(|digest| {
            digest.source.kind == TruthSourceKind::Matrix
                && digest.source.source_id == "SUPPLIED_01"
        }));
        assert!(document.header.source_digests.iter().any(|digest| {
            digest.source.kind == TruthSourceKind::Matrix && digest.source.source_id == "REUSED_01"
        }));
        let case = document
            .records
            .iter()
            .find_map(|record| match record {
                TruthRecord::Genotype(genotype)
                    if genotype.case_kind == GenotypeCaseKind::Snp && genotype.locus.start == 1 =>
                {
                    Some(genotype)
                }
                _ => None,
            })
            .expect("matrix catalog case");
        assert_eq!(
            case.alleles,
            vec![
                TruthAllele::Reference,
                TruthAllele::Snp { base: "A".into() },
                TruthAllele::Snp { base: "T".into() },
            ]
        );
        assert_eq!(case.carriers.len(), 2);
        assert_eq!(
            case.allele_dosages
                .iter()
                .map(|dose| u16::from(*dose))
                .sum::<u16>(),
            2
        );
    }

    #[test]
    fn projection_keeps_authorized_empty_matrix_canonical_identity() {
        let reference = reference();
        let (state, matrix) = f2_3_state_with(&reference, 2, &[], Vec::new());
        assert!(matrix.empty_authorized());
        let matrix_artifact = matrix_artifact_source(&matrix);
        let sources = required_source_digests(&reference, &matrix, []);
        let groups = required_group_digests(&matrix);
        let document = project_truth(TruthProjectionInput {
            execution_id: "RUN_EMPTY",
            reference: &reference,
            individuals: &state,
            matrix: &matrix,
            source_digests: &sources,
            group_digests: &groups,
        })
        .expect("authorized empty matrix projection");
        assert!(groups.is_empty());
        assert!(
            document
                .header
                .source_digests
                .iter()
                .any(|digest| digest.source == matrix_artifact)
        );
        assert_eq!(document.header.genotype_count, 0);
        assert_eq!(document.header.copy_count, 2);
    }

    #[test]
    fn jsonl_bytes_are_deterministic_ordered_and_hashed() {
        let document = document();
        document.validate_structure().expect("valid fixture");
        let mut first = Vec::new();
        let mut second = Vec::new();
        write_truth_jsonl(&document, &mut first).expect("first write");
        write_truth_jsonl(&document, &mut second).expect("second write");
        assert_eq!(first, second);
        assert_eq!(first.iter().filter(|byte| **byte == b'\n').count(), 2);
        assert!(first.ends_with(b"\n"));
        assert_eq!(
            format!("{:x}", Sha256::digest(&first)),
            "4db6f76dc244b5bdefcd3fc140c8db784b841317776cbd8db32a91b56f578b07"
        );
        // Sin indels, lo único que cambió al pasar de la versión 1 a la 4 es
        // el número de versión: con él devuelto, los bytes son los de antes.
        let as_version_one = std::str::from_utf8(&first).expect("UTF-8").replacen(
            "\"schema_version\":\"4\"",
            "\"schema_version\":\"1\"",
            1,
        );
        assert_eq!(
            format!("{:x}", Sha256::digest(as_version_one.as_bytes())),
            "af91d7f472b3bb9db4369e0c383c687adae2a733f6d958a654a4d237c20c617e"
        );
        let lines = std::str::from_utf8(&first)
            .expect("UTF-8")
            .lines()
            .collect::<Vec<_>>();
        assert_eq!(
            serde_json::from_str::<Value>(lines[0]).expect("header")["record_type"],
            "header"
        );
        assert_eq!(
            serde_json::from_str::<Value>(lines[1]).expect("copy")["record_type"],
            "haplotype_copy"
        );
    }

    #[test]
    fn artifact_uses_f1_path_identity_and_complete_bytes() {
        let root = std::env::temp_dir().join(format!(
            "gps-truth-projection-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let document = document();
        let mut expected = Vec::new();
        write_truth_jsonl(&document, &mut expected).expect("direct writer");
        let record = write_truth_artifact(&mut workspace, &document).expect("artifact");
        assert_eq!(record.id, TRUTH_ARTIFACT_ID);
        assert_eq!(record.kind, TRUTH_ARTIFACT_KIND);
        assert_eq!(record.path, TRUTH_ARTIFACT_PATH);
        let bytes = fs::read(workspace.root().join(TRUTH_ARTIFACT_PATH)).expect("truth bytes");
        assert_eq!(bytes, expected);
        assert_eq!(record.size_bytes, u64::try_from(bytes.len()).expect("size"));
        assert_eq!(record.digest, format!("{:x}", Sha256::digest(&bytes)));
        drop(workspace);
        fs::remove_dir_all(&root).expect("cleanup");
    }

    struct FailingWriter {
        remaining: usize,
        written: Vec<u8>,
    }

    impl Write for FailingWriter {
        fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
            let accepted = self.remaining.min(bytes.len());
            self.written.extend_from_slice(&bytes[..accepted]);
            self.remaining -= accepted;
            if accepted == 0 {
                return Err(io::Error::new(io::ErrorKind::WriteZero, "injected failure"));
            }
            Ok(accepted)
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn writer_reports_truncation_or_failure_to_its_caller() {
        let mut writer = FailingWriter {
            remaining: 17,
            written: Vec::new(),
        };
        let error =
            write_truth_jsonl(&document(), &mut writer).expect_err("injected write failure");
        assert_ne!(error.kind(), io::ErrorKind::Interrupted);
        assert!(!writer.written.is_empty());
    }

    #[test]
    fn workspace_discards_partial_truth_when_its_actual_writer_fails() {
        let root = std::env::temp_dir().join(format!(
            "gps-truth-workspace-failure-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("root");
        let destination = root.join("result");
        let mut workspace = ExecutionWorkspace::create(&destination).expect("workspace");
        let error = workspace
            .write_streamed_artifact(
                TRUTH_ARTIFACT_ID,
                TRUTH_ARTIFACT_KIND,
                Path::new(TRUTH_ARTIFACT_PATH),
                |file| {
                    file.write_all(b"partial\n").map_err(|write_error| {
                        Diagnostic::io("GPS-IO-001", "fallo inyectado")
                            .with("detail", write_error.to_string())
                    })?;
                    Err(Diagnostic::io("GPS-IO-001", "fallo inyectado"))
                },
            )
            .expect_err("writer failure");
        assert_eq!(error.code(), "GPS-IO-001");
        assert!(!workspace.root().join(TRUTH_ARTIFACT_PATH).exists());
        assert!(workspace.artifacts().is_empty());
        drop(workspace);
        fs::remove_dir_all(&root).expect("cleanup");
    }

    #[test]
    fn write_truth_artifact_keeps_catalog_empty_when_finalization_fails() {
        let root = std::env::temp_dir().join(format!(
            "gps-truth-artifact-failure-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("result")).expect("workspace");
        // Fuerza que el rename final del escritor real de truth colisione con
        // un directorio: el temporal debe desaparecer y no hay catálogo.
        fs::create_dir_all(workspace.root().join(TRUTH_ARTIFACT_PATH)).expect("blocking path");
        write_truth_artifact(&mut workspace, &document()).expect_err("finalization failure");
        assert!(workspace.artifacts().is_empty());
        assert!(
            fs::read_dir(workspace.root().join("outputs"))
                .expect("payload entries")
                .all(|entry| !entry
                    .expect("entry")
                    .file_name()
                    .to_string_lossy()
                    .starts_with(".truth.jsonl.gps-tmp"))
        );
        drop(workspace);
        fs::remove_dir_all(&root).expect("cleanup");
    }
}
