//! Materialización física y publicación provisional del VCF de verdad F2.5.5.
//!
//! Este módulo recibe solo DTOs ya construidos por F2.5.2--F2.5.4. No
//! reconstruye ni valida semánticamente `truth.jsonl`: esa comparación
//! independiente es responsabilidad exclusiva de F2.5.6.

use std::collections::{BTreeMap, BTreeSet};
use std::{fs::File, io::BufRead, path::Path};

use crate::diagnostics::Diagnostic;
use crate::identity::ExecutionId;
use crate::publication::ExecutionWorkspace;
use crate::reference::Reference;
use crate::results::ArtifactRecord;
use crate::truth_schema::{TruthDocument, TruthRecord};
use crate::truth_validation::validate_truth_for_vcf;
use crate::vcf_contract::{
    VcfContig, VcfContractError, VcfHeaderContract, VcfIndexKind, VcfOutputMode, write_vcf_index,
    write_vcf_records,
};
use crate::vcf_grouping::group_projected_alleles_from_validated_truth;
use crate::vcf_projection::{VcfProjectionError, project_view_from_validated_truth};
use crate::vcf_representation::prepare_conditional_vcf_records;

pub const VCF_ARTIFACT_ID: &str = "truth-vcf";
pub const VCF_ARTIFACT_KIND: &str = "truth_vcf";
pub const VCF_TEXT_ARTIFACT_PATH: &str = "outputs/truth.vcf";
pub const VCF_BGZF_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz";
pub const VCF_INDEX_ARTIFACT_ID: &str = "truth-vcf-index";
pub const VCF_TBI_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz.tbi";
pub const VCF_CSI_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz.csi";

#[cfg(test)]
thread_local! {
    /// Hace fallar la escritura del índice dentro de la publicación, con el
    /// VCF ya escrito: es el fallo real del índice que el paquete fallido
    /// tiene que tratar como un VCF pedido y no publicado (DT-038).
    pub(crate) static FORCE_INDEX_WRITE_FAILURE: std::cell::Cell<bool> = const { std::cell::Cell::new(false) };
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfPublication {
    pub vcf: ArtifactRecord,
    pub index: Option<ArtifactRecord>,
    pub mode: VcfOutputMode,
    pub index_kind: Option<VcfIndexKind>,
    pub reference_id: String,
    pub samples: Vec<String>,
    /// Bases de secuencia de tracto que el VCF publicado lleva de verdad,
    /// contadas sobre los registros que se acaban de escribir.
    pub tract_bases: u64,
    /// Las mismas bases, repartidas por locus (DTC-097).
    pub tract_bases_by_locus: Vec<(String, u64)>,
}

/// El tope de bases de tracto que el usuario declaró, para aplicarlo sobre lo
/// que la ejecución **va a publicar**: la verdad ya construida más lo que el
/// VCF va a escribir.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct TractOutputCap {
    pub cap: u64,
}

/// Bases de tracto que la verdad lleva, por locus: la vista vigente de cada
/// copia, su historial y los alelos de cada caso.
fn truth_tract_bases_by_locus(truth: &TruthDocument) -> BTreeMap<String, u64> {
    use crate::truth_schema::{TruthAllele, TruthCurrentPiece, TruthHistoryEvent};
    let mut totals = BTreeMap::<String, u64>::new();
    let mut add = |locus: &str, bases: &str| {
        let entry = totals.entry(locus.to_owned()).or_insert(0);
        *entry = entry.saturating_add(bases.len() as u64);
    };
    for record in &truth.records {
        match record {
            TruthRecord::HaplotypeCopy(copy) => {
                for piece in &copy.current_pieces {
                    if let TruthCurrentPiece::RepeatAllele {
                        locus_id, bases, ..
                    } = piece
                    {
                        add(locus_id, bases);
                    }
                }
                for event in &copy.history {
                    if let TruthHistoryEvent::RepeatAllele {
                        locus_id, bases, ..
                    } = event
                    {
                        add(locus_id, bases);
                    }
                }
            }
            TruthRecord::Genotype(genotype) => {
                for allele in &genotype.alleles {
                    if let TruthAllele::Repeat {
                        locus_id, bases, ..
                    } = allele
                    {
                        add(locus_id, bases);
                    }
                }
            }
        }
    }
    totals
}

/// Reparte por locus las bases de tracto de los registros de repetición.
///
/// Cada registro se atribuye al locus cuyo intervalo de referencia contiene; la
/// verdad publica ese intervalo y su identificador en los casos de
/// repetición. Si un registro contiene varios, se atribuye entero al primero
/// por posición, de modo que la suma por loci es siempre la del VCF.
fn tract_bases_by_locus(
    truth: &TruthDocument,
    records: &[crate::vcf_grouping::GroupedVcfRecord],
) -> Result<Vec<(String, u64)>, Diagnostic> {
    let mut loci = BTreeMap::<(String, u64, u64), String>::new();
    for record in &truth.records {
        if let TruthRecord::Genotype(genotype) = record {
            for allele in &genotype.alleles {
                if let crate::truth_schema::TruthAllele::Repeat { locus_id, .. } = allele {
                    loci.insert(
                        (
                            genotype.locus.sequence_id.clone(),
                            genotype.locus.start,
                            genotype.locus.end,
                        ),
                        locus_id.clone(),
                    );
                }
            }
        }
    }
    let mut totals = BTreeMap::<String, u64>::new();
    for locus in loci.values() {
        totals.entry(locus.clone()).or_insert(0);
    }
    for record in records.iter().filter(|record| {
        record
            .alternate_repeats
            .iter()
            .any(|sequences| !sequences.is_empty())
    }) {
        let end = record.position + record.reference.len() as u64;
        // El primero en el orden del mapa, como al recorrerlo entero, pero
        // solo entre los loci de su secuencia que empiezan dentro del registro:
        // un locus contenido empieza en `[position, end]`. Recorrer todos los
        // loci por cada registro era cuadrático (F6.4.a.8a).
        let locus = loci
            .range((record.sequence_id.clone(), record.position, 0)..)
            .take_while(|((sequence, start, _), _)| {
                *sequence == record.sequence_id && *start <= end
            })
            .find(|((sequence, start, stop), _)| {
                *sequence == record.sequence_id && record.position <= *start && *stop <= end
            })
            .map(|(_, id)| id)
            .ok_or_else(|| {
                vcf_invariant(
                    "tract_bases_by_locus",
                    "un registro de repetición no contiene el intervalo de ningún locus",
                )
            })?;
        let bases = record
            .alternates
            .iter()
            .fold(record.reference.len() as u64, |total, alternate| {
                total.saturating_add(alternate.len() as u64)
            });
        let entry = totals.entry(locus.clone()).or_insert(0);
        *entry = entry.saturating_add(bases);
    }
    Ok(totals.into_iter().collect())
}

/// Bases de tracto de los registros de repetición: su `REF` y cada uno de sus
/// `ALT`, que desde la decisión S-09 llevan la secuencia completa.
///
/// Cuenta **lo que se publica**, no lo que se declaró. No comparte nada con la
/// cota, que se calcula en la frontera de configuración desde los alelos
/// declarados, los individuos y la ploidía.
#[must_use]
pub fn realized_vcf_tract_bases(records: &[crate::vcf_grouping::GroupedVcfRecord]) -> u64 {
    records
        .iter()
        .filter(|record| {
            record
                .alternate_repeats
                .iter()
                .any(|sequences| !sequences.is_empty())
        })
        .fold(0_u64, |total, record| {
            let alternates = record.alternates.iter().fold(0_u64, |sum, alternate| {
                sum.saturating_add(alternate.len() as u64)
            });
            total
                .saturating_add(record.reference.len() as u64)
                .saturating_add(alternates)
        })
}

/// Construye y cataloga exactamente un VCF del conjunto actual.
///
/// El índice se escribe solo después de cerrar el BGZF y ambos artefactos viven
/// en el `ExecutionWorkspace` provisional. Si la indexación falla, no se llega
/// al manifiesto ni al `renameat2` final de F1.
pub fn publish_truth_vcf(
    workspace: &mut ExecutionWorkspace,
    reference: &Reference,
    truth: &TruthDocument,
    mode: VcfOutputMode,
) -> Result<VcfPublication, Diagnostic> {
    publish_truth_vcf_within(workspace, reference, truth, mode, None)
}

/// Como [`publish_truth_vcf`], aplicando antes de escribir nada el tope de
/// bases de tracto que el usuario declaró (DTC-097).
///
/// Las bases de tracto del VCF se **miden** sobre los registros ya agrupados,
/// que son exactamente los que se van a escribir. Si las de la verdad ya
/// construida más lo medido superan el tope, se rechaza con un diagnóstico de
/// configuración que nombra el tope, la cifra y el locus que más aporta, y no
/// se escribe el VCF.
pub fn publish_truth_vcf_within(
    workspace: &mut ExecutionWorkspace,
    reference: &Reference,
    truth: &TruthDocument,
    mode: VcfOutputMode,
    cap: Option<&TractOutputCap>,
) -> Result<VcfPublication, Diagnostic> {
    let samples = samples_from_truth(truth);
    // Históricamente los fallos del contrato de cabecera tenían precedencia
    // sobre la comprobación semántica de verdad; se conserva ese diagnóstico.
    let header = header_contract(reference, truth, samples.clone())?;
    // Esta validación es anterior a cualquier escritor VCF. Reconstruye una
    // copia cada vez, contrasta todos los genotipos y solo retiene metadatos
    // compactos para que la ruta posterior no pueda observar verdad inválida.
    truth
        .validate_structure()
        .map_err(|error| vcf_invariant("validate_truth_structure", format!("{error:?}")))?;
    let validated = validate_truth_for_vcf(truth, reference)
        // Antes de la vista compacta, estos errores emergían de proyección.
        // Mantener el contexto público no retrasa la validación al escritor.
        .map_err(|error| vcf_invariant("project", error.to_string()))?;
    let view = project_view_from_validated_truth(validated, reference)
        .map_err(|error| match error {
            // No es una invariante rota sino un escenario que el VCF no sabe
            // anotar sin inventar una regla: se le dice a quien lo escribió.
            VcfProjectionError::UndeclaredCompoundTract {
                locus_id,
                sequence_id,
                start,
                end,
            } => Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", "reference.repeat_loci,direct_events")
            .with("locus_id", locus_id.clone())
            .with(
                "detail",
                format!(
                    "locus {locus_id} ({sequence_id}:{} a {end}): hay que anotar en el VCF sin bloques declarados un alelo de este microsatélite compuesto —una copia que conserva el tracto de la referencia dentro de un registro fundido con una inserción o deleción directa, o dos copias con las mismas bases y anotaciones distintas— y con varios motivos su estructura en bloques no se conoce sin ambigüedad",
                    start + 1
                ),
            ),
            VcfProjectionError::SeveralTractsInOneRecord { locus_ids } => {
                Diagnostic::configuration(
                    "GPS-CONF-008",
                    "La configuración contiene una combinación de campos incompatible.",
                )
                .with("fields", "reference.repeat_loci,direct_events")
                .with("locus_id", locus_ids.join(","))
                .with(
                    "detail",
                    format!(
                        "loci {}: una inserción o deleción directa funde en un solo registro VCF varios microsatélites, y la anotación de repetición de un alelo hay que derivarla de sus bases sin saber cuáles son de cada tracto",
                        locus_ids.join(", ")
                    ),
                )
            }
            other => vcf_invariant("project", other.to_string()),
        })?;
    let grouped = group_projected_alleles_from_validated_truth(validated, reference, &view)
        .map_err(|error| vcf_invariant("group", error.to_string()))?;
    let conditional = prepare_conditional_vcf_records(grouped)
        .map_err(|error| vcf_invariant("representation", error.to_string()))?;
    debug_assert!(
        conditional.complement().is_none(),
        "F2 cannot emit complement entries"
    );

    let tract_bases = realized_vcf_tract_bases(&conditional.records);
    let by_locus = tract_bases_by_locus(truth, &conditional.records)?;
    if let Some(cap) = cap {
        // Lo que la ejecución publicará de verdad: las bases de tracto que la
        // verdad ya construida lleva, más las que el VCF va a escribir.
        let truth_by_locus = truth_tract_bases_by_locus(truth);
        let truth_bases = truth_by_locus
            .values()
            .fold(0_u64, |total, bases| total.saturating_add(*bases));
        let total = truth_bases.saturating_add(tract_bases);
        if total > cap.cap {
            let mut per_locus = truth_by_locus.clone();
            for (locus, bases) in &by_locus {
                let entry = per_locus.entry(locus.clone()).or_insert(0);
                *entry = entry.saturating_add(*bases);
            }
            let heaviest = per_locus
                .iter()
                .map(|(id, bases)| (*bases, id.as_str()))
                .max()
                .map_or_else(
                    || "ninguno".to_owned(),
                    |(bases, id)| format!("{id} con {bases} bases"),
                );
            return Err(Diagnostic::configuration(
                "GPS-CONF-008",
                "La configuración contiene una combinación de campos incompatible.",
            )
            .with("fields", "reference.max_tract_output_bases")
            .with(
                "detail",
                format!(
                    "el tope declarado de bases de tracto es {} y la ejecución publicaría {total}: {truth_bases} en la verdad más {tract_bases} medidas en el VCF que se iba a escribir; el locus que más aporta es {heaviest}",
                    cap.cap
                ),
            ));
        }
    }

    let (relative_path, vcf) = match mode {
        VcfOutputMode::Text => {
            let path = Path::new(VCF_TEXT_ARTIFACT_PATH);
            let artifact = workspace.write_streamed_artifact(
                VCF_ARTIFACT_ID,
                VCF_ARTIFACT_KIND,
                path,
                |file| {
                    write_vcf_records(file, &header, &conditional.records)
                        .map(|_| ())
                        .map_err(|error| vcf_contract_error(path, "write_text", error))
                },
            )?;
            (path, artifact)
        }
        VcfOutputMode::Bgzf => {
            let path = Path::new(VCF_BGZF_ARTIFACT_PATH);
            let artifact = workspace.write_streamed_artifact(
                VCF_ARTIFACT_ID,
                VCF_ARTIFACT_KIND,
                path,
                |file| {
                    let bgzf = noodles_bgzf::io::Writer::new(file);
                    write_vcf_records(bgzf, &header, &conditional.records)
                        .map_err(|error| vcf_contract_error(path, "write_bgzf", error))?
                        .finish()
                        .map_err(|error| vcf_io_error(path, "finish_bgzf", error))?;
                    Ok(())
                },
            )?;
            (path, artifact)
        }
    };

    let (index, index_kind) = if mode == VcfOutputMode::Bgzf {
        let kind = header.selected_index();
        let index_path = match kind {
            VcfIndexKind::Tbi => Path::new(VCF_TBI_ARTIFACT_PATH),
            VcfIndexKind::Csi => Path::new(VCF_CSI_ARTIFACT_PATH),
        };
        let vcf_path = workspace.root().join(relative_path);
        let artifact = workspace.write_streamed_artifact(
            VCF_INDEX_ARTIFACT_ID,
            match kind {
                VcfIndexKind::Tbi => "vcf_tbi_index",
                VcfIndexKind::Csi => "vcf_csi_index",
            },
            index_path,
            |file| {
                #[cfg(test)]
                if FORCE_INDEX_WRITE_FAILURE.replace(false) {
                    return Err(vcf_io_detail(
                        index_path,
                        "write_index",
                        "injected index writer failure",
                    ));
                }
                let actual = write_vcf_index(&vcf_path, &header, file)
                    .map_err(|error| vcf_contract_error(index_path, "write_index", error))?;
                if actual != kind {
                    return Err(vcf_invariant("write_index", "tipo de índice incoherente"));
                }
                Ok(())
            },
        )?;
        (Some(artifact), Some(kind))
    } else {
        (None, None)
    };

    // La comprobación de esta frontera es estructural, no semántica: relee
    // completamente la vista ya cerrada y, si existe, abre y consulta su único
    // índice. La comparación independiente con truth.jsonl es F2.5.6.
    validate_materialized_vcf(
        &workspace.root().join(relative_path),
        mode,
        index
            .as_ref()
            .map(|artifact| workspace.root().join(&artifact.path)),
        index_kind,
        &header,
        &conditional.records,
    )?;

    Ok(VcfPublication {
        vcf,
        index,
        mode,
        index_kind,
        reference_id: reference.id().as_str().to_owned(),
        samples,
        tract_bases,
        tract_bases_by_locus: by_locus,
    })
}

fn validate_materialized_vcf(
    vcf_path: &Path,
    mode: VcfOutputMode,
    index_path: Option<std::path::PathBuf>,
    index_kind: Option<VcfIndexKind>,
    header: &VcfHeaderContract,
    records: &[crate::vcf_grouping::GroupedVcfRecord],
) -> Result<(), Diagnostic> {
    match mode {
        VcfOutputMode::Text => read_complete_vcf(
            std::io::BufReader::new(
                File::open(vcf_path).map_err(|error| vcf_io_error(vcf_path, "open_text", error))?,
            ),
            vcf_path,
        )?,
        VcfOutputMode::Bgzf => read_complete_vcf(
            noodles_bgzf::io::Reader::new(
                File::open(vcf_path).map_err(|error| vcf_io_error(vcf_path, "open_bgzf", error))?,
            ),
            vcf_path,
        )?,
    }

    if mode == VcfOutputMode::Bgzf {
        let (index_path, index_kind) = match (index_path, index_kind) {
            (Some(path), Some(kind)) => (path, kind),
            _ => {
                return Err(vcf_invariant(
                    "validate_index",
                    "un VCF BGZF materializado no tiene su índice único",
                ));
            }
        };
        // Incluso una vista BGZF sin variantes tiene que llevar un índice
        // legible de su tipo declarado. Para vistas no vacías, esta misma
        // frontera abre el lector indexado y recupera cada registro esperado.
        validate_all_indexed_records(vcf_path, &index_path, index_kind, header, records)?;
    }

    Ok(())
}

fn read_complete_vcf<R>(source: R, path: &Path) -> Result<(), Diagnostic>
where
    R: BufRead,
{
    let mut reader = noodles_vcf::io::Reader::new(source);
    let header = reader
        .read_header()
        .map_err(|error| vcf_io_error(path, "read_header", error))?;
    let mut record = noodles_vcf::variant::RecordBuf::default();
    while reader
        .read_record_buf(&header, &mut record)
        .map_err(|error| vcf_io_error(path, "read_record", error))?
        != 0
    {}
    Ok(())
}

fn validate_all_indexed_records(
    vcf_path: &Path,
    index_path: &Path,
    kind: VcfIndexKind,
    header: &VcfHeaderContract,
    records: &[crate::vcf_grouping::GroupedVcfRecord],
) -> Result<(), Diagnostic> {
    match kind {
        VcfIndexKind::Tbi => {
            let index = noodles_tabix::io::Reader::new(
                File::open(index_path)
                    .map_err(|error| vcf_io_error(index_path, "open_tbi", error))?,
            )
            .read_index()
            .map_err(|error| vcf_io_error(index_path, "read_tbi", error))?;
            query_expected_records(vcf_path, index_path, index, header, records)
        }
        VcfIndexKind::Csi => {
            let index = noodles_csi::io::Reader::new(
                File::open(index_path)
                    .map_err(|error| vcf_io_error(index_path, "open_csi", error))?,
            )
            .read_index()
            .map_err(|error| vcf_io_error(index_path, "read_csi", error))?;
            query_expected_records(vcf_path, index_path, index, header, records)
        }
    }
}

fn query_expected_records<I>(
    vcf_path: &Path,
    index_path: &Path,
    index: I,
    header_contract: &VcfHeaderContract,
    records: &[crate::vcf_grouping::GroupedVcfRecord],
) -> Result<(), Diagnostic>
where
    I: noodles_csi::BinningIndex + 'static,
{
    let source =
        File::open(vcf_path).map_err(|error| vcf_io_error(vcf_path, "open_bgzf", error))?;
    let mut reader = noodles_vcf::io::indexed_reader::IndexedReader::new(source, index);
    let header = reader
        .read_header()
        .map_err(|error| vcf_io_error(vcf_path, "read_header", error))?;
    for expected in records {
        let position = expected
            .position
            .checked_add(1)
            .and_then(|value| usize::try_from(value).ok())
            .and_then(|value| noodles_core::Position::try_from(value).ok())
            .ok_or_else(|| vcf_invariant("query_index", "posición VCF fuera de rango"))?;
        if !header_contract
            .contigs()
            .iter()
            .any(|contig| contig.name() == expected.sequence_id)
        {
            return Err(vcf_invariant("query_index", "contig del registro ausente"));
        }
        let region = noodles_core::Region::new(expected.sequence_id.as_str(), position..=position);
        let mut query = reader
            .query(&header, &region)
            .map_err(|error| vcf_io_error(index_path, "seek_query", error))?;
        let mut observed = noodles_vcf::Record::default();
        if query
            .read_record(&mut observed)
            .map_err(|error| vcf_io_error(vcf_path, "read_query_record", error))?
            == 0
        {
            return Err(vcf_io_detail(
                index_path,
                "query_coverage",
                "la consulta regional no recupera el registro VCF esperado",
            ));
        }
        let observed_position = observed
            .variant_start()
            .transpose()
            .map_err(|error| vcf_io_error(vcf_path, "read_query_position", error))?;
        if observed.reference_sequence_name() != expected.sequence_id
            || observed_position != Some(position)
        {
            return Err(vcf_io_detail(
                index_path,
                "query_coverage",
                "la consulta regional recupera un registro distinto al esperado",
            ));
        }
    }
    Ok(())
}

fn header_contract(
    reference: &Reference,
    truth: &TruthDocument,
    samples: Vec<String>,
) -> Result<VcfHeaderContract, Diagnostic> {
    if truth.header.reference_id != reference.id().as_str() {
        return Err(vcf_invariant(
            "header",
            "truth.jsonl y la referencia efectiva no comparten reference_id",
        ));
    }
    let execution_id = ExecutionId::parse(truth.header.execution_id.clone(), "VCF header")
        .map_err(|error| vcf_invariant("header", format!("{error:?}")))?;
    let contigs = reference
        .sequences()
        .iter()
        .map(|sequence| VcfContig::new(sequence.id().as_str(), sequence.length()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| vcf_invariant("header", error.to_string()))?;
    VcfHeaderContract::new(
        crate::vcf_contract::VcfHeaderIdentity::new(reference.id().clone(), execution_id),
        contigs,
        samples,
    )
    .map_err(|error| vcf_invariant("header", error.to_string()))
}

fn samples_from_truth(truth: &TruthDocument) -> Vec<String> {
    truth
        .records
        .iter()
        .filter_map(|record| match record {
            TruthRecord::HaplotypeCopy(copy) => Some(copy.identity.individual_id.clone()),
            TruthRecord::Genotype(_) => None,
        })
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn vcf_contract_error(path: &Path, operation: &str, error: VcfContractError) -> Diagnostic {
    match error {
        VcfContractError::Io(detail) => vcf_io_detail(path, operation, detail),
        error => vcf_invariant(operation, error.to_string()),
    }
}

fn vcf_io_error(path: &Path, operation: &str, error: std::io::Error) -> Diagnostic {
    vcf_io_detail(path, operation, error.to_string())
}

fn vcf_io_detail(path: &Path, operation: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::io("GPS-IO-001", "Una operación de lectura o escritura falló.")
        .with("path", path.display().to_string())
        .with("operation", operation)
        .with("detail", detail.into())
}

fn vcf_invariant(operation: &str, detail: impl Into<String>) -> Diagnostic {
    Diagnostic::internal(
        "GPS-INTERNAL-001",
        "Una invariante interna considerada imposible se ha roto.",
    )
    .with("component", "vcf_publication")
    .with("operation", operation)
    .with("detail", detail.into())
}

#[cfg(test)]
mod tests {
    use std::{
        fs::{self, File},
        io::{self, BufReader, Write},
        path::PathBuf,
        process::Command,
        sync::atomic::{AtomicU64, Ordering},
    };

    use super::*;
    use crate::identity::{ReferenceId, SequenceId};
    use crate::reference::{
        DnaSequence, FastaHeader, ReferenceOrigin, ReferenceSequence, SequenceClass,
    };
    use crate::truth_schema::{
        AlleleCarrier, CopyIdentity, CopyProvenance, GenotypeCaseKind, GenotypeRecord,
        HaplotypeCopyRecord, ModificationProvenance, ModificationProvenanceKind, NamedDigest,
        SegmentSource, TRUTH_COORDINATE_CONVENTION, TRUTH_DOCUMENT_KIND, TRUTH_SCHEMA_VERSION,
        TruthAllele, TruthCurrentPiece, TruthCurrentSnp, TruthHeader, TruthHistoryEvent,
        TruthLocus, TruthModification, TruthPieceLink, TruthRecord, TruthSegment, TruthSourceKind,
        TruthSourceRef, TruthStructuralEvent,
    };
    use crate::truth_validation::reference_content_sha256;
    use crate::vcf_contract::{TBI_MAXIMUM_POSITION_EXCLUSIVE, write_vcf_index, write_vcf_records};
    use crate::vcf_grouping::{GroupedVcfRecord, PhasedGenotype};
    use crate::vcf_semantic_validation::validate_vcf_against_truth;

    static NONCE: AtomicU64 = AtomicU64::new(0);

    struct FailingWriter;

    impl Write for FailingWriter {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("injected index failure"))
        }

        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }

    #[test]
    fn regional_validation_seeks_and_reads_all_records_for_tbi_and_csi() {
        for use_csi in [false, true] {
            let root = temporary(if use_csi {
                "regional-csi"
            } else {
                "regional-tbi"
            });
            fs::create_dir_all(&root).expect("root");
            let header = header(use_csi);
            let records = records();
            let vcf = root.join("truth.vcf.gz");
            write_vcf_records(
                noodles_bgzf::io::Writer::new(File::create(&vcf).expect("VCF")),
                &header,
                &records,
            )
            .expect("write VCF")
            .finish()
            .expect("finish BGZF");
            let index = root.join(if use_csi {
                "truth.vcf.gz.csi"
            } else {
                "truth.vcf.gz.tbi"
            });
            let kind = write_vcf_index(&vcf, &header, File::create(&index).expect("index"))
                .expect("write index");
            validate_all_indexed_records(&vcf, &index, kind, &header, &records)
                .expect("every regional query reads its BGZF record");
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    #[test]
    fn corrupted_index_or_missing_coverage_is_an_io_failure() {
        let root = temporary("corrupt-index");
        fs::create_dir_all(&root).expect("root");
        let header = header(false);
        let records = records();
        let vcf = root.join("truth.vcf.gz");
        write_vcf_records(
            noodles_bgzf::io::Writer::new(File::create(&vcf).expect("VCF")),
            &header,
            &records,
        )
        .expect("write VCF")
        .finish()
        .expect("finish BGZF");
        let index = root.join("truth.vcf.gz.tbi");
        let kind = write_vcf_index(&vcf, &header, File::create(&index).expect("index"))
            .expect("write index");

        let mut uncovered = records.clone();
        uncovered[1].position = 8;
        let coverage = validate_all_indexed_records(&vcf, &index, kind, &header, &uncovered)
            .expect_err("query must reject an index that does not cover an expected record");
        assert_eq!(coverage.code(), "GPS-IO-001");
        assert_eq!(coverage.exit_code(), 5);

        fs::write(&index, b"corrupt").expect("corrupt index");
        let corrupt = validate_all_indexed_records(&vcf, &index, kind, &header, &records)
            .expect_err("corrupt offset/index must fail before publication");
        assert_eq!(corrupt.code(), "GPS-IO-001");
        assert_eq!(corrupt.exit_code(), 5);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn empty_bgzf_still_rejects_a_corrupt_tbi_or_csi_before_publication() {
        for use_csi in [false, true] {
            let root = temporary(if use_csi {
                "empty-corrupt-csi"
            } else {
                "empty-corrupt-tbi"
            });
            fs::create_dir_all(&root).expect("root");
            let header = header(use_csi);
            let vcf = root.join("truth.vcf.gz");
            write_vcf_records(
                noodles_bgzf::io::Writer::new(File::create(&vcf).expect("VCF")),
                &header,
                &[],
            )
            .expect("write empty VCF")
            .finish()
            .expect("finish BGZF");
            let index = root.join(if use_csi {
                "truth.vcf.gz.csi"
            } else {
                "truth.vcf.gz.tbi"
            });
            let kind = write_vcf_index(&vcf, &header, File::create(&index).expect("index"))
                .expect("write empty index");

            validate_materialized_vcf(
                &vcf,
                VcfOutputMode::Bgzf,
                Some(index.clone()),
                Some(kind),
                &header,
                &[],
            )
            .expect("empty VCF has a valid readable index");

            fs::write(&index, b"corrupt").expect("corrupt index");
            let error = validate_materialized_vcf(
                &vcf,
                VcfOutputMode::Bgzf,
                Some(index),
                Some(kind),
                &header,
                &[],
            )
            .expect_err("corrupt empty index must fail before publication");
            assert_eq!(error.code(), "GPS-IO-001");
            assert_eq!(error.exit_code(), 5);
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    #[test]
    fn failing_index_writer_maps_to_the_public_io_contract() {
        let root = temporary("failing-index-writer");
        fs::create_dir_all(&root).expect("root");
        let header = header(false);
        let vcf = root.join("truth.vcf.gz");
        write_vcf_records(
            noodles_bgzf::io::Writer::new(File::create(&vcf).expect("VCF")),
            &header,
            &records(),
        )
        .expect("write VCF")
        .finish()
        .expect("finish BGZF");
        let error = write_vcf_index(&vcf, &header, FailingWriter).expect_err("injected error");
        let diagnostic = vcf_contract_error(Path::new(VCF_TBI_ARTIFACT_PATH), "write_index", error);
        assert_eq!(diagnostic.code(), "GPS-IO-001");
        assert_eq!(diagnostic.exit_code(), 5);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn corrupt_truth_is_rejected_before_creating_a_vcf_artifact() {
        let root = temporary("corrupt-truth");
        fs::create_dir_all(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("published")).expect("workspace");
        let truth = TruthDocument {
            header: TruthHeader {
                document_kind: "not_truth".to_owned(),
                schema_version: "1".to_owned(),
                execution_id: "RUN_01".to_owned(),
                reference_id: "REF_main".to_owned(),
                reference_content_sha256: "a".repeat(64),
                coordinate_convention: "zero_based_half_open_v1".to_owned(),
                copy_count: 0,
                genotype_count: 0,
                source_digests: Vec::new(),
                group_digests: Vec::new(),
            },
            records: Vec::new(),
        };
        let error = publish_truth_vcf(&mut workspace, &reference(), &truth, VcfOutputMode::Text)
            .expect_err("invalid truth must stop before VCF writing");

        assert_eq!(error.code(), "GPS-INTERNAL-001");
        assert_eq!(
            error.context().get("operation").map(String::as_str),
            Some("validate_truth_structure")
        );
        assert!(workspace.artifacts().is_empty());
        assert!(!workspace.root().join(VCF_TEXT_ARTIFACT_PATH).exists());
        drop(workspace);
        fs::remove_dir_all(root).expect("cleanup");
    }

    #[test]
    fn semantic_truth_failures_keep_project_diagnostic_and_write_nothing() {
        let mut digest = valid_truth();
        digest.header.reference_content_sha256 = "0".repeat(64);

        let mut history = valid_truth();
        let TruthRecord::HaplotypeCopy(copy) = &mut history.records[0] else {
            panic!("copy");
        };
        copy.segments[0].modifications = vec![TruthModification::Snp {
            event_id: "SNP_01".to_owned(),
            position: 1,
            reference: "C".to_owned(),
            alternate: "T".to_owned(),
            provenance: explicit_provenance(),
        }];
        copy.history = vec![TruthHistoryEvent::Snp {
            event_id: "SNP_01".to_owned(),
            position: 1,
            provenance: explicit_provenance(),
        }];

        let mut genotype = valid_truth();
        let TruthRecord::HaplotypeCopy(copy) = &mut genotype.records[0] else {
            panic!("copy");
        };
        copy.segments[0].modifications = vec![TruthModification::Snp {
            event_id: "SNP_01".to_owned(),
            position: 1,
            reference: "C".to_owned(),
            alternate: "T".to_owned(),
            provenance: explicit_provenance(),
        }];
        copy.current_pieces = vec![TruthCurrentPiece::Source {
            source: source(),
            snps: vec![TruthCurrentSnp {
                event_id: "SNP_01".to_owned(),
                position: 1,
                reference: "C".to_owned(),
                alternate: "T".to_owned(),
                provenance: explicit_provenance(),
            }],
        }];
        copy.history = vec![TruthHistoryEvent::Snp {
            event_id: "SNP_01".to_owned(),
            position: 1,
            provenance: explicit_provenance(),
        }];
        genotype.header.genotype_count = 1;
        genotype.records.push(TruthRecord::Genotype(GenotypeRecord {
            case_id: "CASE_01".to_owned(),
            reference_id: "REF_main".to_owned(),
            individual_id: "IND_01".to_owned(),
            case_kind: GenotypeCaseKind::Snp,
            locus: TruthLocus {
                sequence_id: "chr1".to_owned(),
                start: 1,
                end: 2,
            },
            alleles: vec![
                TruthAllele::Reference,
                TruthAllele::Snp {
                    base: "G".to_owned(),
                },
            ],
            allele_dosages: vec![0, 1],
            carriers: vec![AlleleCarrier {
                allele_index: 1,
                copy: copy_identity(),
            }],
        }));

        for (label, truth) in [
            ("digest", digest),
            ("history", history),
            ("genotype", genotype),
        ] {
            let root = temporary(label);
            fs::create_dir_all(&root).expect("root");
            let mut workspace =
                ExecutionWorkspace::create(&root.join("published")).expect("workspace");
            let error =
                publish_truth_vcf(&mut workspace, &reference(), &truth, VcfOutputMode::Text)
                    .expect_err("structurally valid semantic corruption must stop publication");

            assert_eq!(error.code(), "GPS-INTERNAL-001", "{label}");
            assert_eq!(
                error.context().get("operation").map(String::as_str),
                Some("project"),
                "{label}"
            );
            assert_eq!(
                error.context().get("component").map(String::as_str),
                Some("vcf_publication"),
                "{label}"
            );
            assert!(workspace.artifacts().is_empty(), "{label}");
            assert!(
                !workspace.root().join(VCF_TEXT_ARTIFACT_PATH).exists(),
                "{label}"
            );
            drop(workspace);
            fs::remove_dir_all(root).expect("cleanup");
        }
    }

    #[test]
    fn large_identical_direct_insertions_materialize_once_per_vcf_alt() {
        if let Ok(value) = std::env::var("GPS_VCF_HWM_COPY_COUNT") {
            let copies = value.parse::<u8>().expect("copy count");
            let (reference, truth) = direct_insertion_truth(copies, false);
            let peak_before = peak_resident_set_bytes().expect("Linux peak RSS");
            let root = temporary("large-direct-insertion");
            fs::create_dir_all(&root).expect("root");
            let mut workspace =
                ExecutionWorkspace::create(&root.join("published")).expect("workspace");
            publish_truth_vcf(&mut workspace, &reference, &truth, VcfOutputMode::Text)
                .expect("publish valid insertion truth");
            let delta = peak_resident_set_bytes()
                .expect("Linux peak RSS")
                .saturating_sub(peak_before);
            let vcf = workspace.root().join(VCF_TEXT_ARTIFACT_PATH);
            let text = fs::read_to_string(&vcf).expect("VCF");
            let record = text
                .lines()
                .find(|line| !line.starts_with('#'))
                .expect("one VCF record")
                .split('\t')
                .collect::<Vec<_>>();
            assert_eq!(record[0], "chr1");
            assert_eq!(record[1], "2");
            assert_eq!(record[3], "A");
            assert_eq!(record[4].len(), 1_000_001);
            assert!(record[4].as_bytes()[1..].iter().all(|base| *base == b'G'));
            assert_eq!(
                text.lines().filter(|line| !line.starts_with('#')).count(),
                1
            );
            validate_vcf_against_truth(
                &reference,
                &truth,
                &mut BufReader::new(File::open(&vcf).expect("VCF")),
            )
            .expect("independent semantic VCF oracle");
            drop(workspace);
            fs::remove_dir_all(root).expect("cleanup");
            println!("GPS_VCF_HWM_DELTA_BYTES={delta}");
            return;
        }

        let ten = measured_large_insertion_hwm(10);
        let twenty = measured_large_insertion_hwm(20);
        // La entrada ya contiene una cadena por copia. Se mide únicamente el
        // incremento de RSS máximo desde que está construida; la proyección y
        // agrupación no pueden añadir otra inserción de 1 MiB por portador.
        assert!(
            twenty <= ten.saturating_add(4 * 1024 * 1024),
            "20 copias añadió {twenty} bytes frente a {ten} para 10; el ALT se materializó por copia"
        );

        let (reference, truth) = direct_insertion_truth(2, true);
        let root = temporary("different-large-direct-insertions");
        fs::create_dir_all(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("published")).expect("workspace");
        publish_truth_vcf(&mut workspace, &reference, &truth, VcfOutputMode::Text)
            .expect("publish distinct valid insertion truth");
        let vcf = workspace.root().join(VCF_TEXT_ARTIFACT_PATH);
        let text = fs::read_to_string(&vcf).expect("VCF");
        let alternates = text
            .lines()
            .find(|line| !line.starts_with('#'))
            .expect("one VCF record")
            .split('\t')
            .nth(4)
            .expect("ALT")
            .split(',')
            .collect::<Vec<_>>();
        assert_eq!(alternates.len(), 2, "distinct insertions stay distinct ALT");
        validate_vcf_against_truth(
            &reference,
            &truth,
            &mut BufReader::new(File::open(&vcf).expect("VCF")),
        )
        .expect("independent semantic VCF oracle for distinct insertions");
        drop(workspace);
        fs::remove_dir_all(root).expect("cleanup");
    }

    fn measured_large_insertion_hwm(copies: u8) -> u64 {
        let output = Command::new(std::env::current_exe().expect("test executable"))
            .args([
                "--exact",
                "vcf_publication::tests::large_identical_direct_insertions_materialize_once_per_vcf_alt",
                "--nocapture",
                "--test-threads=1",
            ])
            .env("GPS_VCF_HWM_COPY_COUNT", copies.to_string())
            .output()
            .expect("isolated RSS measurement process");
        assert!(
            output.status.success(),
            "isolated measurement failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        String::from_utf8(output.stdout)
            .expect("measurement output")
            .lines()
            .find_map(|line| {
                line.split_once("GPS_VCF_HWM_DELTA_BYTES=")
                    .map(|(_, value)| value)
            })
            .expect("RSS measurement")
            .trim()
            .parse::<u64>()
            .expect("RSS bytes")
    }

    fn peak_resident_set_bytes() -> Option<u64> {
        fs::read_to_string("/proc/self/status")
            .ok()?
            .lines()
            .find_map(|line| {
                line.strip_prefix("VmHWM:")?
                    .split_whitespace()
                    .next()?
                    .parse::<u64>()
                    .ok()
            })
            .and_then(|kibibytes| kibibytes.checked_mul(1024))
    }

    fn direct_insertion_truth(copy_count: u8, distinct: bool) -> (Reference, TruthDocument) {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "large insertion VCF test").expect("sequence ID"),
            FastaHeader::parse(b"chr1").expect("FASTA header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(b"CAAAA").expect("DNA"),
        )
        .expect("reference sequence");
        let reference = Reference::new(
            ReferenceId::parse("REF_main", "large insertion VCF test").expect("reference ID"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference");
        let shared_inserted = "G".repeat(1_000_000);
        let mut copy_records = Vec::with_capacity(usize::from(copy_count));
        let mut genotype_records = Vec::with_capacity(usize::from(copy_count));
        for index in 0..copy_count {
            let individual = format!("IND_{index:02}");
            let identity = CopyIdentity {
                reference_id: "REF_main".to_owned(),
                individual_id: individual.clone(),
                sequence_id: "chr1".to_owned(),
                copy_index: 0,
            };
            let source = |start, end| SegmentSource {
                source: TruthSourceRef {
                    kind: TruthSourceKind::Reference,
                    source_id: "REF_main".to_owned(),
                },
                reference_id: "REF_main".to_owned(),
                sequence_id: "chr1".to_owned(),
                start,
                end,
                sequence_length: 5,
            };
            let event_id = format!("INS_{index:02}");
            let provenance = explicit_provenance();
            let inserted = if distinct && index % 2 == 1 {
                "T".repeat(1_000_000)
            } else {
                shared_inserted.clone()
            };
            copy_records.push(TruthRecord::HaplotypeCopy(HaplotypeCopyRecord {
                identity: identity.clone(),
                phase: 0,
                provenance: CopyProvenance {
                    origin: "test".to_owned(),
                    sources: vec![TruthSourceRef {
                        kind: TruthSourceKind::Reference,
                        source_id: "REF_main".to_owned(),
                    }],
                },
                segments: vec![TruthSegment {
                    source: source(0, 5),
                    modifications: Vec::new(),
                }],
                current_pieces: vec![
                    TruthCurrentPiece::Source {
                        source: source(0, 2),
                        snps: Vec::new(),
                    },
                    TruthCurrentPiece::Inserted {
                        bases: inserted.clone(),
                        provenance: provenance.clone(),
                    },
                    TruthCurrentPiece::Source {
                        source: source(2, 5),
                        snps: Vec::new(),
                    },
                ],
                history: vec![TruthHistoryEvent::Insertion {
                    event_id: event_id.clone(),
                    point: 2,
                    bases: inserted.clone(),
                    provenance: provenance.clone(),
                }],
            }));
            genotype_records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id: format!("CASE_{index:02}"),
                reference_id: "REF_main".to_owned(),
                individual_id: individual,
                case_kind: GenotypeCaseKind::Indel,
                locus: TruthLocus {
                    sequence_id: "chr1".to_owned(),
                    start: 2,
                    end: 2,
                },
                alleles: vec![
                    TruthAllele::Reference,
                    TruthAllele::Insertion {
                        event_id: event_id.clone(),
                        bases: inserted.clone(),
                        source_event: TruthStructuralEvent::Insertion {
                            event_id,
                            point: 2,
                            bases: inserted.clone(),
                        },
                        piece_links: vec![TruthPieceLink {
                            piece_index: 1,
                            start: 2,
                            end: 1_000_002,
                        }],
                    },
                ],
                allele_dosages: vec![0, 1],
                carriers: vec![AlleleCarrier {
                    allele_index: 1,
                    copy: identity,
                }],
            }));
        }
        copy_records.extend(genotype_records);
        let reference_digest = reference_content_sha256(&reference).expect("huella");
        (
            reference,
            TruthDocument {
                header: TruthHeader {
                    document_kind: TRUTH_DOCUMENT_KIND.to_owned(),
                    schema_version: TRUTH_SCHEMA_VERSION.to_owned(),
                    execution_id: "RUN_01".to_owned(),
                    reference_id: "REF_main".to_owned(),
                    reference_content_sha256: reference_digest,
                    coordinate_convention: TRUTH_COORDINATE_CONVENTION.to_owned(),
                    copy_count: u64::from(copy_count),
                    genotype_count: u64::from(copy_count),
                    source_digests: vec![
                        NamedDigest {
                            source: TruthSourceRef {
                                kind: TruthSourceKind::Reference,
                                source_id: "REF_main".to_owned(),
                            },
                            sha256: "a".repeat(64),
                        },
                        NamedDigest {
                            source: TruthSourceRef {
                                kind: TruthSourceKind::Explicit,
                                source_id: "EXPLICIT_01".to_owned(),
                            },
                            sha256: "b".repeat(64),
                        },
                    ],
                    group_digests: Vec::new(),
                },
                records: copy_records,
            },
        )
    }

    #[test]
    fn header_contract_keeps_precedence_over_a_reference_id_mismatch() {
        let mut truth = valid_truth();
        truth.header.reference_id = "REF_other".to_owned();
        let root = temporary("reference-id-precedence");
        fs::create_dir_all(&root).expect("root");
        let mut workspace = ExecutionWorkspace::create(&root.join("published")).expect("workspace");

        let error = publish_truth_vcf(&mut workspace, &reference(), &truth, VcfOutputMode::Text)
            .expect_err("header reference mismatch");

        assert_eq!(error.code(), "GPS-INTERNAL-001");
        assert_eq!(
            error.context().get("operation").map(String::as_str),
            Some("header")
        );
        assert!(workspace.artifacts().is_empty());
        drop(workspace);
        fs::remove_dir_all(root).expect("cleanup");
    }

    fn header(use_csi: bool) -> VcfHeaderContract {
        let first_length = if use_csi {
            TBI_MAXIMUM_POSITION_EXCLUSIVE
        } else {
            16
        };
        VcfHeaderContract::new(
            crate::vcf_contract::VcfHeaderIdentity::new(
                crate::identity::ReferenceId::parse("ref-main", "VCF publication test")
                    .expect("reference ID"),
                crate::identity::ExecutionId::parse("run-main", "VCF publication test")
                    .expect("execution ID"),
            ),
            vec![
                VcfContig::new("chr1", first_length).expect("chr1"),
                VcfContig::new("chr2", 16).expect("chr2"),
            ],
            vec!["sample-1".to_owned()],
        )
        .expect("header")
    }

    fn records() -> Vec<GroupedVcfRecord> {
        [("chr1", 1, "CASE_1"), ("chr2", 4, "CASE_2")]
            .into_iter()
            .map(|(sequence_id, position, case)| GroupedVcfRecord {
                sequence_id: sequence_id.to_owned(),
                position,
                reference: "A".to_owned(),
                alternates: vec!["C".to_owned()],
                samples: vec![(
                    "sample-1".to_owned(),
                    PhasedGenotype {
                        alleles: vec![1, 0],
                        allele_dosages: vec![1, 1],
                    },
                )],
                gps_cases: vec![case.to_owned()],
                alternate_cases: vec![vec![case.to_owned()]],
                alternate_repeats: vec![Vec::new()],
            })
            .collect()
    }

    fn reference() -> Reference {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "VCF publication test").expect("sequence ID"),
            FastaHeader::parse(b"chr1").expect("FASTA header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(b"ACGT").expect("DNA"),
        )
        .expect("reference sequence");
        Reference::new(
            ReferenceId::parse("REF_main", "VCF publication test").expect("reference ID"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("reference")
    }

    fn valid_truth() -> TruthDocument {
        let reference = reference();
        TruthDocument {
            header: TruthHeader {
                document_kind: TRUTH_DOCUMENT_KIND.to_owned(),
                schema_version: TRUTH_SCHEMA_VERSION.to_owned(),
                execution_id: "RUN_01".to_owned(),
                reference_id: "REF_main".to_owned(),
                reference_content_sha256: reference_content_sha256(&reference).expect("huella"),
                coordinate_convention: TRUTH_COORDINATE_CONVENTION.to_owned(),
                copy_count: 1,
                genotype_count: 0,
                source_digests: vec![
                    NamedDigest {
                        source: TruthSourceRef {
                            kind: TruthSourceKind::Reference,
                            source_id: "REF_main".to_owned(),
                        },
                        sha256: "a".repeat(64),
                    },
                    NamedDigest {
                        source: TruthSourceRef {
                            kind: TruthSourceKind::Explicit,
                            source_id: "EXPLICIT_01".to_owned(),
                        },
                        sha256: "b".repeat(64),
                    },
                ],
                group_digests: Vec::new(),
            },
            records: vec![TruthRecord::HaplotypeCopy(HaplotypeCopyRecord {
                identity: copy_identity(),
                phase: 0,
                provenance: CopyProvenance {
                    origin: "test".to_owned(),
                    sources: vec![TruthSourceRef {
                        kind: TruthSourceKind::Reference,
                        source_id: "REF_main".to_owned(),
                    }],
                },
                segments: vec![TruthSegment {
                    source: source(),
                    modifications: Vec::new(),
                }],
                current_pieces: vec![TruthCurrentPiece::Source {
                    source: source(),
                    snps: Vec::new(),
                }],
                history: Vec::new(),
            })],
        }
    }

    fn copy_identity() -> CopyIdentity {
        CopyIdentity {
            reference_id: "REF_main".to_owned(),
            individual_id: "IND_01".to_owned(),
            sequence_id: "chr1".to_owned(),
            copy_index: 0,
        }
    }

    fn source() -> SegmentSource {
        SegmentSource {
            source: TruthSourceRef {
                kind: TruthSourceKind::Reference,
                source_id: "REF_main".to_owned(),
            },
            reference_id: "REF_main".to_owned(),
            sequence_id: "chr1".to_owned(),
            start: 0,
            end: 4,
            sequence_length: 4,
        }
    }

    fn explicit_provenance() -> ModificationProvenance {
        ModificationProvenance {
            kind: ModificationProvenanceKind::Explicit,
            version: "1".to_owned(),
            source: TruthSourceRef {
                kind: TruthSourceKind::Explicit,
                source_id: "EXPLICIT_01".to_owned(),
            },
            matrix: None,
        }
    }

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-f2-5-5-vcf-publication-{}-{}-{label}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }
}
