//! Materialización física y publicación provisional del VCF de verdad F2.5.5.
//!
//! Este módulo recibe solo DTOs ya construidos por F2.5.2--F2.5.4. No
//! reconstruye ni valida semánticamente `truth.jsonl`: esa comparación
//! independiente es responsabilidad exclusiva de F2.5.6.

use std::collections::BTreeSet;
use std::{fs::File, io::BufRead, path::Path};

use crate::diagnostics::Diagnostic;
use crate::publication::ExecutionWorkspace;
use crate::reference::Reference;
use crate::results::ArtifactRecord;
use crate::truth_schema::{TruthDocument, TruthRecord};
use crate::vcf_contract::{
    VcfContig, VcfContractError, VcfHeaderContract, VcfIndexKind, VcfOutputMode, write_vcf_index,
    write_vcf_records,
};
use crate::vcf_grouping::group_normalized_alleles;
use crate::vcf_projection::project_normalized_alleles;
use crate::vcf_representation::prepare_conditional_vcf_records;

pub const VCF_ARTIFACT_ID: &str = "truth-vcf";
pub const VCF_ARTIFACT_KIND: &str = "truth_vcf";
pub const VCF_TEXT_ARTIFACT_PATH: &str = "outputs/truth.vcf";
pub const VCF_BGZF_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz";
pub const VCF_INDEX_ARTIFACT_ID: &str = "truth-vcf-index";
pub const VCF_TBI_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz.tbi";
pub const VCF_CSI_ARTIFACT_PATH: &str = "outputs/truth.vcf.gz.csi";

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct VcfPublication {
    pub vcf: ArtifactRecord,
    pub index: Option<ArtifactRecord>,
    pub mode: VcfOutputMode,
    pub index_kind: Option<VcfIndexKind>,
    pub reference_id: String,
    pub samples: Vec<String>,
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
    let samples = samples_from_truth(truth);
    let header = header_contract(reference, samples.clone())?;
    let normalized = project_normalized_alleles(truth, reference)
        .map_err(|error| vcf_invariant("project", error.to_string()))?;
    let grouped = group_normalized_alleles(truth, reference, &normalized)
        .map_err(|error| vcf_invariant("group", error.to_string()))?;
    let conditional = prepare_conditional_vcf_records(grouped)
        .map_err(|error| vcf_invariant("representation", error.to_string()))?;
    debug_assert!(
        conditional.complement().is_none(),
        "F2 cannot emit complement entries"
    );

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
    samples: Vec<String>,
) -> Result<VcfHeaderContract, Diagnostic> {
    let contigs = reference
        .sequences()
        .iter()
        .map(|sequence| VcfContig::new(sequence.id().as_str(), sequence.length()))
        .collect::<Result<Vec<_>, _>>()
        .map_err(|error| vcf_invariant("header", error.to_string()))?;
    VcfHeaderContract::new(contigs, samples)
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
        io::{self, Write},
        path::PathBuf,
        sync::atomic::{AtomicU64, Ordering},
    };

    use super::*;
    use crate::vcf_contract::{TBI_MAXIMUM_POSITION_EXCLUSIVE, write_vcf_index, write_vcf_records};
    use crate::vcf_grouping::{GroupedVcfRecord, PhasedGenotype};

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

    fn header(use_csi: bool) -> VcfHeaderContract {
        let first_length = if use_csi {
            TBI_MAXIMUM_POSITION_EXCLUSIVE
        } else {
            16
        };
        VcfHeaderContract::new(
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
            })
            .collect()
    }

    fn temporary(label: &str) -> PathBuf {
        std::env::temp_dir().join(format!(
            "gps-f2-5-5-vcf-publication-{}-{}-{label}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ))
    }
}
