//! Oráculo semántico independiente para la vista VCF de verdad F2.
//!
//! Esta frontera lee los tres artefactos publicados. No llama a
//! `vcf_projection`, `vcf_grouping`, `vcf_representation` ni al escritor VCF:
//! la verdad se valida y se interpreta antes de abrir el VCF.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::identity::{ReferenceId, SequenceId};
use crate::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use crate::truth_schema::{CopyIdentity, TruthCurrentPiece, TruthDocument, TruthRecord};
use crate::truth_validation::{
    TruthValidationError, read_truth_jsonl, validate_and_reconstruct_truth,
};

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum VcfSemanticValidationError {
    /// Se pidieron bases de una referencia validada sin conservarlas.
    BasesNotMaterialized,
    Io {
        path: String,
        detail: String,
    },
    Fasta {
        detail: String,
    },
    Truth(TruthValidationError),
    Vcf {
        line: u64,
        detail: String,
    },
}

impl std::fmt::Display for VcfSemanticValidationError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(formatter, "{self:?}")
    }
}

impl std::error::Error for VcfSemanticValidationError {}

/// Valida un VCF de texto o BGZF publicado frente a sus tres autoridades.
///
/// El JSONL se lee y valida completamente antes de abrir el VCF. De este modo
/// ningún dato VCF puede completar, corregir ni alimentar la verdad.
pub fn validate_published_vcf_semantics(
    fasta_path: &Path,
    truth_path: &Path,
    vcf_path: &Path,
) -> Result<(), VcfSemanticValidationError> {
    let truth_file = File::open(truth_path).map_err(|error| io_error(truth_path, error))?;
    let truth = read_truth_jsonl(&mut BufReader::new(truth_file))
        .map_err(VcfSemanticValidationError::Truth)?;
    let reference = read_oracle_fasta(fasta_path, &truth.header.reference_id)?;
    // La validación completa de la verdad y de la referencia precede incluso
    // a `File::open(vcf_path)`: un VCF ausente o corrupto nunca puede ocultar
    // una relación truth--FASTA rota.
    validate_declared_copy_phases(&truth)?;
    validate_and_reconstruct_truth(truth.clone(), &reference)
        .map_err(VcfSemanticValidationError::Truth)?;

    if vcf_path
        .extension()
        .is_some_and(|extension| extension == "gz")
    {
        let file = File::open(vcf_path).map_err(|error| io_error(vcf_path, error))?;
        validate_vcf_against_truth(&reference, &truth, &mut noodles_bgzf::io::Reader::new(file))
    } else {
        let file = File::open(vcf_path).map_err(|error| io_error(vcf_path, error))?;
        validate_vcf_against_truth(&reference, &truth, &mut BufReader::new(file))
    }
}

/// Comprueba una vista VCF ya abierta contra la referencia y la verdad.
///
/// Está expuesta también para corpus sintéticos, sin necesitar publicación ni
/// rutas. El lector de VCF es deliberadamente pequeño: aplica el subconjunto
/// F2 (GT faseado y `GPS_CASE`) sin compartir el parser productor.
pub fn validate_vcf_against_truth(
    reference: &Reference,
    truth: &TruthDocument,
    reader: &mut impl BufRead,
) -> Result<(), VcfSemanticValidationError> {
    // Este reconstructor valida que el JSONL sea autosuficiente antes de que
    // el oráculo examine un solo byte VCF.
    validate_declared_copy_phases(truth)?;
    validate_and_reconstruct_truth(truth.clone(), reference)
        .map_err(VcfSemanticValidationError::Truth)?;
    let expected = expected_records(reference, truth)?;
    let expected_samples = samples_from_truth(truth);
    let sequence_order = reference
        .sequences()
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.id().as_str(), index))
        .collect::<BTreeMap<_, _>>();

    let mut line = String::new();
    let mut line_number = 0_u64;
    let mut samples = None;
    let mut seen = BTreeSet::new();
    let mut previous = None;
    let mut vcf_45 = false;
    let mut reference_declarations = 0_u8;
    let mut execution_declarations = 0_u8;
    let mut qual_declarations = 0_u8;
    let mut pass_filter_declarations = 0_u8;
    let mut gt_declarations = 0_u8;
    let mut gps_case_declarations = 0_u8;
    let mut records_started = false;
    loop {
        line.clear();
        let count =
            reader
                .read_line(&mut line)
                .map_err(|error| VcfSemanticValidationError::Io {
                    path: "VCF".to_owned(),
                    detail: error.to_string(),
                })?;
        if count == 0 {
            break;
        }
        line_number = line_number.saturating_add(1);
        let content = line.trim_end_matches(['\n', '\r']);
        if content.starts_with("##") {
            if records_started {
                return vcf_error(
                    line_number,
                    "una declaración ## no puede aparecer tras #CHROM",
                );
            }
            if line_number == 1 && content == "##fileformat=VCFv4.5" {
                vcf_45 = true;
            }
            if content.starts_with("##reference=") {
                reference_declarations = reference_declarations.saturating_add(1);
                if content != format!("##reference={}", truth.header.reference_id) {
                    return vcf_error(
                        line_number,
                        "declaración reference distinta de la identidad truth F2",
                    );
                }
            }
            if content.starts_with("##gps_execution_id=") {
                execution_declarations = execution_declarations.saturating_add(1);
                if content != format!("##gps_execution_id={}", truth.header.execution_id) {
                    return vcf_error(
                        line_number,
                        "declaración gps_execution_id distinta de la identidad truth F2",
                    );
                }
            }
            if content.starts_with("##gps_qual=") {
                qual_declarations = qual_declarations.saturating_add(1);
                if content != "##gps_qual=." {
                    return vcf_error(line_number, "declaración gps_qual distinta de QUAL=. en F2");
                }
            }
            if content.starts_with("##FILTER=<") {
                pass_filter_declarations = pass_filter_declarations.saturating_add(1);
                if content != "##FILTER=<ID=PASS,Description=\"All filters passed\">" {
                    return vcf_error(
                        line_number,
                        "declaración FILTER distinta de PASS fijo en F2",
                    );
                }
            }
            if content.starts_with("##FORMAT=<") {
                gt_declarations = gt_declarations.saturating_add(1);
                if content != "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" {
                    return vcf_error(line_number, "declaración FORMAT distinta de GT fijo en F2");
                }
            }
            if is_gps_case_declaration(content) {
                gps_case_declarations = gps_case_declarations.saturating_add(1);
                if !is_exact_gps_case_header(content) {
                    return vcf_error(
                        line_number,
                        "declaración INFO/GPS_CASE distinta del contrato canónico F2",
                    );
                }
            }
            if content.starts_with("##INFO=<ID=GPS_LOSS,")
                || content.starts_with("##INFO=<ID=GPS_COMPLEMENT,")
            {
                return vcf_error(line_number, "F2 prohíbe declarar GPS_LOSS y GPS_COMPLEMENT");
            }
            continue;
        }
        if content.starts_with('#') {
            if records_started {
                return vcf_error(line_number, "una línea # no puede aparecer tras #CHROM");
            }
            if content.starts_with("#CHROM\t") {
                if samples.is_some() {
                    return vcf_error(line_number, "cabecera #CHROM duplicada");
                }
                let fields = content.split('\t').collect::<Vec<_>>();
                if fields.len() < 9
                    || fields[..9]
                        != [
                            "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
                        ]
                {
                    return vcf_error(
                        line_number,
                        "cabecera de columnas VCF distinta del contrato F2",
                    );
                }
                let parsed = fields[9..]
                    .iter()
                    .map(|value| (*value).to_owned())
                    .collect::<Vec<_>>();
                if parsed != expected_samples || !strictly_sorted(&parsed) {
                    return vcf_error(
                        line_number,
                        "muestras VCF no coinciden con la verdad en orden canónico",
                    );
                }
                samples = Some(parsed);
                records_started = true;
            } else {
                return vcf_error(line_number, "línea de cabecera VCF no reconocida");
            }
            continue;
        }
        let samples = samples
            .as_ref()
            .ok_or_else(|| VcfSemanticValidationError::Vcf {
                line: line_number,
                detail: "registro VCF antes de #CHROM".to_owned(),
            })?;
        let fields = content.split('\t').collect::<Vec<_>>();
        if fields.len() != 9 + samples.len() {
            return vcf_error(
                line_number,
                "número de columnas VCF incompatible con las muestras",
            );
        }
        if fields[8] != "GT" {
            return vcf_error(line_number, "FORMAT debe ser exactamente GT en F2");
        }
        if fields[2] != "." || fields[5] != "." || fields[6] != "PASS" {
            return vcf_error(
                line_number,
                "ID, QUAL y FILTER deben ser respectivamente '.', '.' y PASS en F2",
            );
        }
        if fields[7].contains("GPS_LOSS=") || fields[7].contains("GPS_COMPLEMENT=") {
            return vcf_error(line_number, "F2 prohíbe GPS_LOSS y GPS_COMPLEMENT");
        }
        let position = fields[1]
            .parse::<u64>()
            .ok()
            .filter(|position| *position > 0)
            .ok_or_else(|| VcfSemanticValidationError::Vcf {
                line: line_number,
                detail: "POS debe ser un entero positivo".to_owned(),
            })?;
        let key = (fields[0].to_owned(), position);
        let Some(expected_record) = expected.get(&key) else {
            return vcf_error(line_number, "registro CHROM+POS ausente de la verdad");
        };
        if !seen.insert(key.clone()) {
            return vcf_error(line_number, "CHROM+POS está duplicado");
        }
        let order = (
            *sequence_order
                .get(fields[0])
                .ok_or_else(|| VcfSemanticValidationError::Vcf {
                    line: line_number,
                    detail: "CHROM no existe en el FASTA autoritativo".to_owned(),
                })?,
            position,
        );
        if previous.is_some_and(|previous| order <= previous) {
            return vcf_error(
                line_number,
                "los registros no siguen el orden FASTA/POS canónico",
            );
        }
        previous = Some(order);
        if fields[3] != expected_record.reference
            || fields[4] != expected_record.alternates.join(",")
        {
            return vcf_error(
                line_number,
                format!(
                    "REF o ALT no coincide con la derivación independiente (esperado {}>{})",
                    expected_record.reference,
                    expected_record.alternates.join(","),
                ),
            );
        }
        validate_reference_allele(reference, fields[0], position, fields[3], line_number)?;
        let gps_cases = parse_gps_case(fields[7], line_number)?;
        if gps_cases != expected_record.gps_cases || !strictly_sorted(&gps_cases) {
            return vcf_error(
                line_number,
                "GPS_CASE no es la unión canónica de los casos fuente ALT",
            );
        }
        for (sample, observed) in samples.iter().zip(&fields[9..]) {
            let expected_gt = expected_gt(expected_record, sample)?;
            if *observed != expected_gt {
                return vcf_error(
                    line_number,
                    "GT, fase o dosis no coincide con las copias de verdad",
                );
            }
        }
    }
    if samples.is_none() {
        return vcf_error(line_number, "falta la cabecera #CHROM");
    }
    if !vcf_45 {
        return vcf_error(
            1,
            "la primera línea debe ser exactamente ##fileformat=VCFv4.5",
        );
    }
    if reference_declarations != 1 || execution_declarations != 1 {
        return vcf_error(
            0,
            "debe existir exactamente una declaración canónica de referencia y ejecución",
        );
    }
    if qual_declarations != 1 || pass_filter_declarations != 1 || gt_declarations != 1 {
        return vcf_error(
            0,
            "deben existir exactamente las declaraciones QUAL=., FILTER=PASS y FORMAT/GT de F2",
        );
    }
    if gps_case_declarations != 1 {
        return vcf_error(
            0,
            "debe existir exactamente una declaración INFO/GPS_CASE canónica",
        );
    }
    if seen.len() != expected.len() {
        return vcf_error(line_number, "faltan registros VCF que exige la verdad");
    }
    Ok(())
}

#[derive(Clone, Debug)]
struct ExpectedRecord {
    reference: String,
    alternates: Vec<String>,
    gps_cases: Vec<String>,
    copies: BTreeMap<String, Vec<(u8, CopyIdentity, String)>>,
}

fn expected_records(
    reference: &Reference,
    truth: &TruthDocument,
) -> Result<BTreeMap<(String, u64), ExpectedRecord>, VcfSemanticValidationError> {
    let mut paths = Vec::new();
    for record in &truth.records {
        let TruthRecord::HaplotypeCopy(copy) = record else {
            continue;
        };
        let sequence = reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id().as_str() == copy.identity.sequence_id)
            .ok_or_else(|| oracle_error("copia truth nombra una secuencia FASTA ausente"))?;
        paths.push((
            copy.phase,
            oracle_path(
                copy.identity.clone(),
                &copy.current_pieces,
                sequence
                    .sequence()
                    .as_bytes()
                    .map_err(|_| VcfSemanticValidationError::BasesNotMaterialized)?,
            )?,
        ));
    }
    validate_oracle_phases(&paths)?;
    let mut spans = BTreeMap::<String, Vec<OracleSpan>>::new();
    for (_, path) in &paths {
        spans
            .entry(path.sequence_id.clone())
            .or_default()
            .extend(path.changes());
    }
    let mut alleles = Vec::new();
    for (sequence_id, changes) in spans {
        let sequence = reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id().as_str() == sequence_id)
            .ok_or_else(|| oracle_error("secuencia FASTA ausente durante composición"))?;
        let components = merge_oracle_spans(changes);
        for (_, path) in paths
            .iter()
            .filter(|(_, path)| path.sequence_id == sequence_id)
        {
            alleles.extend(normalize_oracle_components(
                path,
                sequence
                    .sequence()
                    .as_bytes()
                    .map_err(|_| VcfSemanticValidationError::BasesNotMaterialized)?,
                &components,
            )?);
        }
    }
    let mut by_record = BTreeMap::<(String, u64), Vec<OracleAllele>>::new();
    for allele in alleles {
        by_record
            .entry((allele.sequence_id.clone(), allele.position + 1))
            .or_default()
            .push(allele);
    }
    by_record
        .into_iter()
        .map(|(key, alleles)| {
            let common_reference = alleles
                .iter()
                .map(|allele| allele.reference.as_str())
                .max_by_key(|value| value.len())
                .ok_or_else(|| oracle_error("registro esperado sin alelos"))?
                .to_owned();
            let mut per_copy = BTreeMap::<CopyIdentity, (String, OracleSpan)>::new();
            for allele in &alleles {
                if !common_reference.starts_with(&allele.reference) {
                    return Err(oracle_error(
                        "REF incompatibles tras composición independiente",
                    ));
                }
                let alternate = format!(
                    "{}{}",
                    allele.alternate,
                    &common_reference[allele.reference.len()..]
                );
                if let Some((previous, _)) =
                    per_copy.insert(allele.copy.clone(), (alternate.clone(), allele.source))
                    && previous != alternate
                {
                    return Err(oracle_error(
                        "una copia produce ALT incompatibles al mismo CHROM+POS",
                    ));
                }
            }
            let alternates = per_copy
                .values()
                .map(|(alternate, _)| alternate.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let mut alt_cases = BTreeMap::<String, BTreeSet<String>>::new();
            for genotype in truth.records.iter().filter_map(|record| {
                if let TruthRecord::Genotype(value) = record {
                    Some(value)
                } else {
                    None
                }
            }) {
                for carrier in genotype
                    .carriers
                    .iter()
                    .filter(|carrier| carrier.allele_index != 0)
                {
                    if let Some((alternate, span)) = per_copy.get(&carrier.copy)
                        && carrier.copy.sequence_id == key.0
                        && oracle_span_contains(*span, genotype.locus.start, genotype.locus.end)
                    {
                        alt_cases
                            .entry(alternate.clone())
                            .or_default()
                            .insert(genotype.case_id.clone());
                    }
                }
            }
            if alternates
                .iter()
                .any(|alternate| alt_cases.get(alternate).is_none_or(BTreeSet::is_empty))
            {
                return Err(oracle_error("ALT compuesto sin caso truth fuente"));
            }
            let gps_cases = alt_cases
                .values()
                .flatten()
                .cloned()
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            let mut copies = BTreeMap::<String, Vec<(u8, CopyIdentity, String)>>::new();
            for (phase, path) in &paths {
                if path.sequence_id == key.0 {
                    let alternate = per_copy
                        .get(&path.identity)
                        .map_or_else(String::new, |(alternate, _)| alternate.clone());
                    copies
                        .entry(path.identity.individual_id.clone())
                        .or_default()
                        .push((*phase, path.identity.clone(), alternate));
                }
            }
            for values in copies.values_mut() {
                values.sort_by_key(|(phase, copy, _)| (*phase, copy.clone()));
            }
            Ok((
                key,
                ExpectedRecord {
                    reference: common_reference,
                    alternates,
                    gps_cases,
                    copies,
                },
            ))
        })
        .collect()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
struct OracleSpan {
    start: u64,
    end: u64,
}
#[derive(Clone, Debug)]
struct OracleAllele {
    sequence_id: String,
    position: u64,
    reference: String,
    alternate: String,
    source: OracleSpan,
    copy: CopyIdentity,
}
#[derive(Clone, Debug)]
struct OraclePath {
    identity: CopyIdentity,
    sequence_id: String,
    pieces: Vec<OraclePiece>,
}
#[derive(Clone, Debug)]
enum OraclePiece {
    Source {
        start: u64,
        end: u64,
        bases: Vec<u8>,
        snps: Vec<u64>,
    },
    Inserted {
        anchor: u64,
        bases: Vec<u8>,
    },
}

fn oracle_path(
    identity: CopyIdentity,
    current: &[TruthCurrentPiece],
    reference: &[u8],
) -> Result<OraclePath, VcfSemanticValidationError> {
    let mut pieces = Vec::new();
    let mut prior = None;
    for piece in current {
        match piece {
            TruthCurrentPiece::Source { source, snps } => {
                if source.start >= source.end
                    || source.end > u64::try_from(reference.len()).unwrap_or(u64::MAX)
                    || prior.is_some_and(|end| source.start < end)
                {
                    return Err(oracle_error("piezas source truth no colineales"));
                }
                let start = usize::try_from(source.start)
                    .map_err(|_| oracle_error("coordenada source demasiado grande"))?;
                let end = usize::try_from(source.end)
                    .map_err(|_| oracle_error("coordenada source demasiado grande"))?;
                let mut bases = reference[start..end].to_vec();
                let mut positions = Vec::new();
                for snp in snps {
                    let offset = usize::try_from(
                        snp.position
                            .checked_sub(source.start)
                            .ok_or_else(|| oracle_error("SNP fuera de source"))?,
                    )
                    .map_err(|_| oracle_error("coordenada SNP demasiado grande"))?;
                    if offset >= bases.len() {
                        return Err(oracle_error("SNP fuera de source"));
                    }
                    bases[offset] = snp.alternate.as_bytes()[0];
                    positions.push(snp.position);
                }
                pieces.push(OraclePiece::Source {
                    start: source.start,
                    end: source.end,
                    bases,
                    snps: positions,
                });
                prior = Some(source.end);
            }
            TruthCurrentPiece::Inserted { bases, .. } => {
                pieces.push(OraclePiece::Inserted {
                    anchor: prior
                        .ok_or_else(|| oracle_error("inserción truth sin source precedente"))?,
                    bases: bases.as_bytes().to_vec(),
                });
            }
        }
    }
    Ok(OraclePath {
        sequence_id: identity.sequence_id.clone(),
        identity,
        pieces,
    })
}

impl OraclePath {
    fn changes(&self) -> Vec<OracleSpan> {
        let mut out = Vec::new();
        let mut prior = None;
        for piece in &self.pieces {
            match piece {
                OraclePiece::Source {
                    start, end, snps, ..
                } => {
                    if let Some(previous) = prior
                        && previous < *start
                    {
                        out.push(OracleSpan {
                            start: previous,
                            end: *start,
                        });
                    }
                    out.extend(snps.iter().map(|position| OracleSpan {
                        start: *position,
                        end: position + 1,
                    }));
                    prior = Some(*end);
                }
                OraclePiece::Inserted { anchor, .. } => out.push(OracleSpan {
                    start: *anchor,
                    end: *anchor,
                }),
            }
        }
        out
    }
    fn window(&self, start: u64, end: u64) -> Result<Vec<u8>, VcfSemanticValidationError> {
        let mut out = Vec::new();
        for piece in &self.pieces {
            match piece {
                OraclePiece::Source {
                    start: left,
                    end: right,
                    bases,
                    ..
                } if *left < end && start < *right => {
                    let from = usize::try_from((*left).max(start) - left)
                        .map_err(|_| oracle_error("coordenada window demasiado grande"))?;
                    let to = usize::try_from((*right).min(end) - left)
                        .map_err(|_| oracle_error("coordenada window demasiado grande"))?;
                    out.extend_from_slice(&bases[from..to]);
                }
                OraclePiece::Inserted { anchor, bases } if start <= *anchor && *anchor <= end => {
                    out.extend_from_slice(bases)
                }
                _ => {}
            }
        }
        Ok(out)
    }
}

fn normalize_oracle_components(
    path: &OraclePath,
    reference: &[u8],
    initial: &[OracleSpan],
) -> Result<Vec<OracleAllele>, VcfSemanticValidationError> {
    let mut components = initial.to_vec();
    loop {
        let normalized = components
            .iter()
            .copied()
            .map(|span| oracle_normalize_window(path, reference, span).map(|value| (span, value)))
            .collect::<Result<Vec<_>, _>>()?;
        let mut merge = Vec::new();
        for (index, (span, allele)) in normalized.iter().enumerate() {
            let Some(allele) = allele else { continue };
            let end = allele.position
                + u64::try_from(allele.reference.len())
                    .map_err(|_| oracle_error("REF demasiado largo"))?;
            let traversed = OracleSpan {
                start: span.start.min(allele.position),
                end: span.end.max(end),
            };
            for (other_index, other) in components.iter().enumerate() {
                if other_index == index {
                    continue;
                }
                if (allele.position < other.end && other.start < end)
                    || (traversed.start <= other.start && other.end <= traversed.end)
                {
                    merge.push((*span, *other));
                }
            }
        }
        if merge.is_empty() {
            return Ok(normalized
                .into_iter()
                .filter_map(|(span, allele)| {
                    allele.map(|mut allele| {
                        allele.source = span;
                        allele
                    })
                })
                .collect());
        }
        for (left, right) in merge {
            let envelope = OracleSpan {
                start: left.start.min(right.start),
                end: left.end.max(right.end),
            };
            components.retain(|item| *item != left && *item != right);
            components.push(envelope);
        }
        components = merge_oracle_spans(components);
    }
}

fn oracle_normalize_window(
    path: &OraclePath,
    reference_bases: &[u8],
    span: OracleSpan,
) -> Result<Option<OracleAllele>, VcfSemanticValidationError> {
    let mut start = span.start;
    let end = span.end;
    let mut reference = reference_bases[usize::try_from(start)
        .map_err(|_| oracle_error("inicio fuera de rango"))?
        ..usize::try_from(end).map_err(|_| oracle_error("fin fuera de rango"))?]
        .to_vec();
    let mut alternate = path.window(start, end)?;
    if reference == alternate {
        return Ok(None);
    }
    if reference.is_empty() || alternate.is_empty() {
        start = start
            .checked_sub(1)
            .ok_or_else(|| oracle_error("indel sin ancla izquierda"))?;
        reference = reference_bases[usize::try_from(start)
            .map_err(|_| oracle_error("inicio fuera de rango"))?
            ..usize::try_from(end).map_err(|_| oracle_error("fin fuera de rango"))?]
            .to_vec();
        alternate = path.window(start, end)?;
    }
    while reference.len() > 1 && alternate.len() > 1 && reference.last() == alternate.last() {
        reference.pop();
        alternate.pop();
    }
    while reference.len() > 1 && alternate.len() > 1 && reference.first() == alternate.first() {
        reference.remove(0);
        alternate.remove(0);
        start += 1;
    }
    while reference.len() != alternate.len() && start > 0 && reference.last() == alternate.last() {
        start -= 1;
        let base = reference_bases
            [usize::try_from(start).map_err(|_| oracle_error("posición fuera de rango"))?];
        reference.insert(0, base);
        alternate.insert(0, base);
        reference.pop();
        alternate.pop();
    }
    Ok(Some(OracleAllele {
        sequence_id: path.sequence_id.clone(),
        position: start,
        reference: String::from_utf8(reference).map_err(|_| oracle_error("REF no UTF-8"))?,
        alternate: String::from_utf8(alternate).map_err(|_| oracle_error("ALT no UTF-8"))?,
        source: span,
        copy: path.identity.clone(),
    }))
}

fn merge_oracle_spans(mut spans: Vec<OracleSpan>) -> Vec<OracleSpan> {
    spans.sort();
    let mut output = Vec::<OracleSpan>::new();
    let mut right_boundary_bridge = None;
    for span in spans {
        if let Some(last) = output.last_mut() {
            let previous_end = last.end;
            let overlaps = if last.start == last.end || span.start == span.end {
                last.start <= span.end && span.start <= last.end
            } else {
                last.start < span.end && span.start < last.end
            };
            let joins_through_bridge =
                right_boundary_bridge == Some(previous_end) && span.start == previous_end;
            if overlaps || joins_through_bridge {
                last.end = last.end.max(span.end);
                right_boundary_bridge = if span.start == span.end && span.start == previous_end {
                    Some(previous_end)
                } else if last.end > previous_end {
                    None
                } else {
                    right_boundary_bridge
                };
                continue;
            }
        }
        output.push(span);
        right_boundary_bridge = None;
    }
    output
}

fn is_exact_gps_case_header(content: &str) -> bool {
    let Some(body) = content
        .strip_prefix("##INFO=<")
        .and_then(|value| value.strip_suffix('>'))
    else {
        return false;
    };
    let mut id = None;
    let mut number = None;
    let mut value_type = None;
    let mut description = None;
    for field in body.split(',') {
        let Some((name, value)) = field.split_once('=') else {
            return false;
        };
        match name {
            "ID" if id.replace(value).is_none() => {}
            "Number" if number.replace(value).is_none() => {}
            "Type" if value_type.replace(value).is_none() => {}
            "Description" if description.replace(value).is_none() => {}
            "ID" | "Number" | "Type" | "Description" => return false,
            _ => {}
        }
    }
    id == Some("GPS_CASE")
        && number == Some(".")
        && value_type == Some("String")
        && description == Some("\"Canonical truth case identifiers contributing an ALT allele\"")
}

fn is_gps_case_declaration(content: &str) -> bool {
    content
        .strip_prefix("##INFO=<")
        .and_then(|value| value.strip_suffix('>'))
        .is_some_and(|body| body.split(',').any(|field| field == "ID=GPS_CASE"))
}
fn oracle_span_contains(span: OracleSpan, start: u64, end: u64) -> bool {
    if start == end {
        span.start <= start && start <= span.end
    } else {
        span.start < end && start < span.end
    }
}
fn validate_oracle_phases(paths: &[(u8, OraclePath)]) -> Result<(), VcfSemanticValidationError> {
    let mut values = BTreeMap::<(String, String), Vec<u8>>::new();
    for (phase, path) in paths {
        values
            .entry((
                path.identity.individual_id.clone(),
                path.sequence_id.clone(),
            ))
            .or_default()
            .push(*phase);
    }
    for ((sample, sequence), mut phases) in values {
        phases.sort_unstable();
        if !(1..=10).contains(&phases.len())
            || phases
                .iter()
                .enumerate()
                .any(|(index, phase)| usize::from(*phase) != index)
        {
            return Err(oracle_error(format!(
                "fases de copia no contiguas para {sample}/{sequence}"
            )));
        }
    }
    Ok(())
}

fn validate_declared_copy_phases(truth: &TruthDocument) -> Result<(), VcfSemanticValidationError> {
    let mut phases = BTreeMap::<(String, String), Vec<u8>>::new();
    for copy in truth.records.iter().filter_map(|record| match record {
        TruthRecord::HaplotypeCopy(copy) => Some(copy),
        TruthRecord::Genotype(_) => None,
    }) {
        phases
            .entry((
                copy.identity.individual_id.clone(),
                copy.identity.sequence_id.clone(),
            ))
            .or_default()
            .push(copy.phase);
    }
    for ((sample, sequence), mut observed) in phases {
        observed.sort_unstable();
        if !(1..=10).contains(&observed.len())
            || observed
                .iter()
                .enumerate()
                .any(|(index, phase)| usize::from(*phase) != index)
        {
            return Err(oracle_error(format!(
                "fases declaradas no contiguas para {sample}/{sequence}"
            )));
        }
    }
    Ok(())
}
fn oracle_error(detail: impl Into<String>) -> VcfSemanticValidationError {
    VcfSemanticValidationError::Vcf {
        line: 0,
        detail: detail.into(),
    }
}

fn expected_gt(
    record: &ExpectedRecord,
    sample: &str,
) -> Result<String, VcfSemanticValidationError> {
    let copies = record
        .copies
        .get(sample)
        .ok_or_else(|| VcfSemanticValidationError::Vcf {
            line: 0,
            detail: "muestra VCF no tiene copias truth para el contig".to_owned(),
        })?;
    copies
        .iter()
        .map(|(_, _, alternate)| {
            if alternate.is_empty() {
                Ok("0".to_owned())
            } else {
                record
                    .alternates
                    .iter()
                    .position(|candidate| candidate == alternate)
                    .and_then(|index| index.checked_add(1))
                    .map(|index| index.to_string())
                    .ok_or_else(|| VcfSemanticValidationError::Vcf {
                        line: 0,
                        detail: "ALT truth no aparece en el registro esperado".to_owned(),
                    })
            }
        })
        .collect::<Result<Vec<_>, _>>()
        .map(|alleles| alleles.join("|"))
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

fn parse_gps_case(info: &str, line: u64) -> Result<Vec<String>, VcfSemanticValidationError> {
    let Some(value) = info.strip_prefix("GPS_CASE=") else {
        return vcf_error(line, "INFO debe declarar GPS_CASE");
    };
    if value.is_empty() || value.contains(';') || value.split(',').any(str::is_empty) {
        return vcf_error(
            line,
            "INFO debe contener exclusivamente GPS_CASE canónico sin valores vacíos",
        );
    }
    Ok(value.split(',').map(str::to_owned).collect())
}

fn validate_reference_allele(
    reference: &Reference,
    chrom: &str,
    position: u64,
    allele: &str,
    line: u64,
) -> Result<(), VcfSemanticValidationError> {
    let sequence = reference
        .sequences()
        .iter()
        .find(|sequence| sequence.id().as_str() == chrom)
        .ok_or_else(|| VcfSemanticValidationError::Vcf {
            line,
            detail: "CHROM no existe en el FASTA autoritativo".to_owned(),
        })?;
    let start = position
        .checked_sub(1)
        .and_then(|value| usize::try_from(value).ok());
    let end = start.and_then(|start| start.checked_add(allele.len()));
    if start.zip(end).and_then(|(start, end)| {
        sequence
            .sequence()
            .as_bytes()
            .ok()
            .and_then(|bases| bases.get(start..end))
    }) != Some(allele.as_bytes())
    {
        return vcf_error(line, "REF no coincide con el FASTA autoritativo");
    }
    Ok(())
}

fn read_oracle_fasta(
    path: &Path,
    reference_id: &str,
) -> Result<Reference, VcfSemanticValidationError> {
    let file = File::open(path).map_err(|error| io_error(path, error))?;
    let mut records = Vec::new();
    let mut line = String::new();
    let mut reader = BufReader::new(file);
    let mut current: Option<(String, Vec<u8>)> = None;
    loop {
        line.clear();
        if reader
            .read_line(&mut line)
            .map_err(|error| io_error(path, error))?
            == 0
        {
            break;
        }
        let content = line.trim_end_matches(['\n', '\r']);
        if let Some(header) = content.strip_prefix('>') {
            if let Some((header, bases)) = current.take() {
                records.push(oracular_sequence(&header, bases, reference_id)?);
            }
            current = Some((header.to_owned(), Vec::new()));
        } else if let Some((_, bases)) = current.as_mut() {
            for base in content.bytes() {
                bases.push(base.to_ascii_uppercase());
            }
        } else {
            return Err(VcfSemanticValidationError::Fasta {
                detail: "bases antes de la primera cabecera FASTA".to_owned(),
            });
        }
    }
    if let Some((header, bases)) = current {
        records.push(oracular_sequence(&header, bases, reference_id)?);
    }
    Reference::new(
        ReferenceId::parse(reference_id, "oráculo VCF").map_err(|error| {
            VcfSemanticValidationError::Fasta {
                detail: format!("{error:?}"),
            }
        })?,
        ReferenceOrigin::Generated,
        records,
    )
    .map_err(|error| VcfSemanticValidationError::Fasta {
        detail: format!("{error:?}"),
    })
}

fn oracular_sequence(
    header: &str,
    bases: Vec<u8>,
    reference_id: &str,
) -> Result<ReferenceSequence, VcfSemanticValidationError> {
    let fasta = FastaHeader::parse(header.as_bytes()).map_err(|error| {
        VcfSemanticValidationError::Fasta {
            detail: format!("{error:?}"),
        }
    })?;
    let id = SequenceId::parse(fasta.id_text(), reference_id).map_err(|error| {
        VcfSemanticValidationError::Fasta {
            detail: format!("{error:?}"),
        }
    })?;
    let dna =
        DnaSequence::try_from_bytes(&bases).map_err(|error| VcfSemanticValidationError::Fasta {
            detail: format!("{error:?}"),
        })?;
    ReferenceSequence::new(id, fasta, SequenceClass::Chromosome, &dna).map_err(|error| {
        VcfSemanticValidationError::Fasta {
            detail: format!("{error:?}"),
        }
    })
}

fn strictly_sorted(values: &[String]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

fn io_error(path: &Path, error: std::io::Error) -> VcfSemanticValidationError {
    VcfSemanticValidationError::Io {
        path: path.display().to_string(),
        detail: error.to_string(),
    }
}

fn vcf_error<T>(line: u64, detail: impl Into<String>) -> Result<T, VcfSemanticValidationError> {
    Err(VcfSemanticValidationError::Vcf {
        line,
        detail: detail.into(),
    })
}

#[cfg(test)]
mod tests {
    use super::merge_oracle_spans;

    #[test]
    fn independently_merges_overlapping_or_inserted_components() {
        assert_eq!(
            merge_oracle_spans(vec![
                super::OracleSpan { start: 2, end: 3 },
                super::OracleSpan { start: 3, end: 3 },
            ]),
            [super::OracleSpan { start: 2, end: 3 }]
        );
    }
}
