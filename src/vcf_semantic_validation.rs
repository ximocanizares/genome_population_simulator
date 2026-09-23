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
use crate::truth_schema::{
    CopyIdentity, TruthCurrentPiece, TruthDocument, TruthRecord, TruthRepeatBlock,
};
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
    let mut tandem_repeat_declarations = BTreeSet::<&'static str>::new();
    let mut records_started = false;
    loop {
        line.clear();
        let count = reader
            .read_line(&mut line)
            .map_err(|error| io_reject("VCF", error.to_string()))?;
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
            if let Some(key) = tandem_repeat_declaration_key(content) {
                if content != canonical_tandem_repeat_declaration(key) {
                    return vcf_error(
                        line_number,
                        "declaración INFO de repetición distinta de la del estándar VCFv4.5",
                    );
                }
                if !tandem_repeat_declarations.insert(key) {
                    return vcf_error(
                        line_number,
                        "declaración INFO de repetición duplicada en la cabecera",
                    );
                }
            }
            continue;
        }
        if content.starts_with('#') {
            if records_started {
                return vcf_error(line_number, "una línea # no puede aparecer tras #CHROM");
            }
            if content.starts_with("#CHROM\t") {
                // Ninguna prueba ejercita esta rama y se cree inalcanzable
                // porque `samples` y `records_started` se asignan en la misma
                // línea de más abajo y en ningún otro sitio: un segundo
                // `#CHROM` ya ha caído arriba por «una línea # no puede
                // aparecer tras #CHROM». Se conserva porque es la condición
                // que el contrato enuncia, y desaparecería sin aviso si
                // alguien separase las dos asignaciones.
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
            .ok_or_else(|| vcf_reject(line_number, "registro VCF antes de #CHROM"))?;
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
            .ok_or_else(|| vcf_reject(line_number, "POS debe ser un entero positivo"))?;
        let key = (fields[0].to_owned(), position);
        let Some(expected_record) = expected.get(&key) else {
            return vcf_error(line_number, "registro CHROM+POS ausente de la verdad");
        };
        if !seen.insert(key.clone()) {
            return vcf_error(line_number, "CHROM+POS está duplicado");
        }
        let order = (
            // Inalcanzable: un CHROM que el FASTA no tiene no llega hasta aquí,
            // porque la búsqueda del registro esperado de unas líneas más
            // arriba ya ha rechazado por «registro CHROM+POS ausente de la
            // verdad». Ninguna prueba la ejercita.
            *sequence_order.get(fields[0]).ok_or_else(|| {
                vcf_reject(line_number, "CHROM no existe en el FASTA autoritativo")
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
        validate_tandem_repeat_info(expected_record, fields[7], line_number)?;
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
    // La cabecera declara **cada campo usado y ninguno más**: los cinco si
    // algún registro publica un tracto, ninguno si no lo publica ninguno.
    let expects_tandem_repeats = expected.values().any(|record| {
        record
            .alternate_repeats
            .iter()
            .any(|blocks| !blocks.is_empty())
    });
    let expected_declarations = if expects_tandem_repeats {
        TANDEM_REPEAT_INFO_KEYS.len()
    } else {
        0
    };
    if tandem_repeat_declarations.len() != expected_declarations {
        return vcf_error(
            0,
            "la cabecera no declara exactamente los campos de repetición que los registros usan",
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
    /// Paralela a `alternates`: los bloques que cada ALT debe publicar.
    alternate_repeats: Vec<Vec<TruthRepeatBlock>>,
}

fn expected_records(
    reference: &Reference,
    truth: &TruthDocument,
) -> Result<BTreeMap<(String, u64), ExpectedRecord>, VcfSemanticValidationError> {
    // Las dos lecturas de bases de esta función rechazan si la referencia no
    // las conserva. **Ninguna prueba las ejercita y se creen inalcanzables por
    // la frontera con rutas**, que construye la referencia leyendo el FASTA y
    // por tanto siempre con bases; solo un llamador de
    // `validate_vcf_against_truth` podría entregar una referencia validada sin
    // ellas, y la batería no tiene ninguna.
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
                    .map_err(|_| bases_not_materialized())?,
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
            // Ninguna prueba ejercita este rechazo y se cree inalcanzable
            // porque las claves de `spans` salen de los caminos ya
            // construidos, y construir cada uno exigió encontrar su secuencia
            // en el FASTA unas líneas más arriba.
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
                    .map_err(|_| bases_not_materialized())?,
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
                // Ninguna prueba ejercita este rechazo y se cree inalcanzable
                // porque cada entrada de `by_record` nace al empujar su primer
                // alelo: el vector nunca está vacío y `max_by_key` nunca da
                // `None`.
                .ok_or_else(|| oracle_error("registro esperado sin alelos"))?
                .to_owned();
            let mut per_copy = BTreeMap::<CopyIdentity, OracleCarrier>::new();
            let mut repeats_by_alternate = BTreeMap::<String, Vec<TruthRepeatBlock>>::new();
            for allele in &alleles {
                // Ninguna prueba ejercita este rechazo y se cree inalcanzable
                // porque todos los alelos de un registro comparten `position`
                // y cada `reference` es la rebanada del FASTA que empieza
                // ahí: el más corto es siempre prefijo del más largo, que es
                // el que `common_reference` elige. Se conserva porque es la
                // premisa de la línea siguiente, que rebana
                // `common_reference` por la longitud del alelo.
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
                if let Some(previous) =
                    repeats_by_alternate.insert(alternate.clone(), allele.repeat_blocks.clone())
                    && previous != allele.repeat_blocks
                {
                    return Err(oracle_error(
                        "dos copias dan el mismo ALT con tractos distintos",
                    ));
                }
                if let Some(previous) = per_copy.insert(
                    allele.copy.clone(),
                    OracleCarrier {
                        alternate: alternate.clone(),
                        source: allele.source,
                    },
                ) && previous.alternate != alternate
                {
                    return Err(oracle_error(
                        "una copia produce ALT incompatibles al mismo CHROM+POS",
                    ));
                }
            }
            let alternates = per_copy
                .values()
                .map(|carrier| carrier.alternate.clone())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect::<Vec<_>>();
            let alternate_repeats = alternates
                .iter()
                .map(|alternate| {
                    repeats_by_alternate
                        .get(alternate)
                        .cloned()
                        .unwrap_or_default()
                })
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
                    if let Some(carrier_allele) = per_copy.get(&carrier.copy)
                        && carrier.copy.sequence_id == key.0
                        && oracle_span_contains(
                            carrier_allele.source,
                            genotype.locus.start,
                            genotype.locus.end,
                        )
                    {
                        alt_cases
                            .entry(carrier_allele.alternate.clone())
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
                        .map_or_else(String::new, |carrier| carrier.alternate.clone());
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
                    alternate_repeats,
                },
            ))
        })
        .collect()
}

/// Lo que una copia aporta a un registro: su ALT y el tramo que lo compuso.
#[derive(Clone, Debug)]
struct OracleCarrier {
    alternate: String,
    source: OracleSpan,
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
    /// Bloques declarados de los tractos contenidos en la ventana.
    repeat_blocks: Vec<TruthRepeatBlock>,
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
    /// Tracto de repetición: sustituye `[start, end)` por sus propias bases.
    Replaced {
        start: u64,
        end: u64,
        bases: Vec<u8>,
        /// Estructura en bloques que la verdad declara para este tracto.
        blocks: Vec<TruthRepeatBlock>,
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
                // Las conversiones de coordenada de esta función y de las tres
                // siguientes —`window`, `normalize_oracle_components` y
                // `oracle_normalize_window`— **no las ejercita ninguna prueba y
                // se creen inalcanzables** en los objetivos de 64 bits que el
                // proyecto compila: la guarda de colinealidad de arriba acota
                // ya cada coordenada por `reference.len()`, que es un `usize`.
                // Se conservan porque el tipo del documento es `u64` y en un
                // objetivo de 32 bits dejarían de serlo.
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
            TruthCurrentPiece::RepeatAllele {
                source,
                bases,
                blocks,
                ..
            } => {
                if source.start >= source.end
                    || source.end > u64::try_from(reference.len()).unwrap_or(u64::MAX)
                    || prior.is_some_and(|end| source.start < end)
                {
                    return Err(oracle_error("piezas source truth no colineales"));
                }
                pieces.push(OraclePiece::Replaced {
                    start: source.start,
                    end: source.end,
                    bases: bases.as_bytes().to_vec(),
                    blocks: blocks.clone(),
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
                OraclePiece::Replaced { start, end, .. } => {
                    if let Some(previous) = prior
                        && previous < *start
                    {
                        out.push(OracleSpan {
                            start: previous,
                            end: *start,
                        });
                    }
                    out.push(OracleSpan {
                        start: *start,
                        end: *end,
                    });
                    prior = Some(*end);
                }
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
                OraclePiece::Replaced {
                    start: left,
                    end: right,
                    bases,
                    ..
                } if *left < end && start < *right => {
                    if *left < start || end < *right {
                        return Err(oracle_error("la ventana parte un tracto de repetición"));
                    }
                    out.extend_from_slice(bases);
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
    // **El registro de repetición no se recorta**, igual que en el productor y
    // por la misma razón: la decisión S-09 manda publicar la secuencia
    // completa en `REF` y `ALT`, y DTC-064 sigue mandando en todos los demás
    // registros. El oráculo aplica la regla por su cuenta —esta es su copia,
    // no la del productor—, que es lo que hace que la comparación signifique
    // algo.
    let repeat_blocks = path.repeats_in(span.start, span.end);
    if repeat_blocks.is_empty() {
        while reference.len() > 1 && alternate.len() > 1 && reference.last() == alternate.last() {
            reference.pop();
            alternate.pop();
        }
        while reference.len() > 1 && alternate.len() > 1 && reference.first() == alternate.first() {
            reference.remove(0);
            alternate.remove(0);
            start += 1;
        }
        while reference.len() != alternate.len()
            && start > 0
            && reference.last() == alternate.last()
        {
            start -= 1;
            let base = reference_bases
                [usize::try_from(start).map_err(|_| oracle_error("posición fuera de rango"))?];
            reference.insert(0, base);
            alternate.insert(0, base);
            reference.pop();
            alternate.pop();
        }
    }
    Ok(Some(OracleAllele {
        sequence_id: path.sequence_id.clone(),
        position: start,
        // Ninguna prueba ejercita estos dos rechazos y se creen inalcanzables,
        // pero **no** por el alfabeto: este lector no comprueba en ninguna
        // parte qué bases trae una pieza insertada de la verdad. Lo que los
        // hace inalcanzables es el recorte de arriba, que solo quita un byte
        // cuando coincide con uno del FASTA —que sí es ASCII—, de modo que
        // nunca puede partir un carácter de varios bytes: si la entrada era
        // texto válido, la salida también lo es.
        reference: String::from_utf8(reference).map_err(|_| oracle_error("REF no UTF-8"))?,
        alternate: String::from_utf8(alternate).map_err(|_| oracle_error("ALT no UTF-8"))?,
        source: span,
        copy: path.identity.clone(),
        repeat_blocks,
    }))
}

impl OraclePath {
    /// Bloques declarados de los tractos contenidos en la ventana, en el orden
    /// en que aparecen en la secuencia.
    fn repeats_in(&self, start: u64, end: u64) -> Vec<TruthRepeatBlock> {
        self.pieces
            .iter()
            .filter_map(|piece| match piece {
                OraclePiece::Replaced {
                    start: left,
                    end: right,
                    blocks,
                    ..
                } if *left < end && start < *right => Some(blocks.iter().cloned()),
                _ => None,
            })
            .flatten()
            .collect()
    }
}

/// Los cinco campos de repetición en tándem, en el orden en que la cabecera
/// los declara.
const TANDEM_REPEAT_INFO_KEYS: [&str; 5] = ["RN", "RUS", "RUL", "RUC", "RB"];

/// Las declaraciones canónicas, **transcritas del estándar** y no importadas
/// del productor.
///
/// Compartir la constante con quien escribe el VCF convertiría esta
/// comprobación en un espejo: el oráculo aceptaría cualquier cosa que el
/// productor decidiera emitir. Por eso están copiadas aquí, igual que ya lo
/// estaban las de `FILTER`, `FORMAT` y `GPS_CASE`.
fn canonical_tandem_repeat_declaration(key: &str) -> &'static str {
    match key {
        "RN" => {
            "##INFO=<ID=RN,Number=A,Type=Integer,Description=\"Total number of repeat sequences in this allele\">"
        }
        "RUS" => {
            "##INFO=<ID=RUS,Number=.,Type=String,Description=\"Repeat unit sequence of the corresponding repeat sequence\">"
        }
        "RUL" => {
            "##INFO=<ID=RUL,Number=.,Type=Integer,Description=\"Repeat unit length of the corresponding repeat sequence\">"
        }
        "RUC" => {
            "##INFO=<ID=RUC,Number=.,Type=Float,Description=\"Repeat unit count of corresponding repeat sequence\">"
        }
        "RB" => {
            "##INFO=<ID=RB,Number=.,Type=Integer,Description=\"Total number of bases in the corresponding repeat sequence\">"
        }
        _ => unreachable!("solo se consultan los cinco campos del estándar"),
    }
}

fn tandem_repeat_declaration_key(content: &str) -> Option<&'static str> {
    TANDEM_REPEAT_INFO_KEYS
        .into_iter()
        .find(|key| content.starts_with(&format!("##INFO=<ID={key},")))
}

/// Valor de una clave `INFO` de la línea, sin interpretar.
fn info_value<'line>(info: &'line str, key: &str) -> Option<&'line str> {
    info.split(';').find_map(|field| {
        field
            .strip_prefix(key)
            .and_then(|rest| rest.strip_prefix('='))
    })
}

/// Comprueba los campos de repetición de un registro contra la verdad.
///
/// # Qué discrimina y qué no
///
/// Las cuatro listas por secuencia repetida se comparan **posición a
/// posición** contra los bloques que este oráculo ha reconstruido para cada
/// alelo alternativo, en el orden de los alelos. Eso es lo que caza el error
/// silencioso más probable de estos campos —una lista desplazada respecto de
/// la de alelos—, y lo caza lista por lista: desplazar solo el número de
/// copias, o solo las bases, falla igual.
///
/// **No se comprueba aparte que las bases publicadas sumen la longitud del
/// tracto.** Esa suma ya está implicada por la comparación bloque a bloque y
/// por la validación de la verdad, que exige que los bloques sumen las bases
/// del alelo: una aserción más no podría fallar, y aquí no se finge que
/// discrimine.
///
/// # Por qué aquí no se trocea, y qué queda sin ejercitar
///
/// Esta función **no reparte las listas entre los alelos**: compara el
/// recuento de cada alelo contra los bloques que la verdad le da, y después
/// las cuatro listas aplanadas posición a posición. Las dos comprobaciones
/// juntas fijan lo mismo que un troceado y no pueden desalinearse entre sí.
///
/// Lo que la producción no ha llegado a construir es un registro con
/// **recuentos distintos entre alelos**: los alelos de un locus comparten su
/// lista de motivos, de modo que todos traen el mismo número de bloques. Para
/// que ese camino no quedara sin recorrer hay una prueba unitaria en este
/// mismo módulo que lo alimenta a mano con dos bloques, ninguno y uno.
fn validate_tandem_repeat_info(
    expected: &ExpectedRecord,
    info: &str,
    line: u64,
) -> Result<(), VcfSemanticValidationError> {
    let expected_blocks = expected
        .alternate_repeats
        .iter()
        .flatten()
        .collect::<Vec<_>>();
    let present = TANDEM_REPEAT_INFO_KEYS
        .into_iter()
        .filter(|key| info_value(info, key).is_some())
        .count();
    if expected_blocks.is_empty() {
        if present != 0 {
            return vcf_error(line, "un registro sin tracto publica campos de repetición");
        }
        return Ok(());
    }
    if present != TANDEM_REPEAT_INFO_KEYS.len() {
        return vcf_error(
            line,
            "un registro con tracto no publica los cinco campos de repetición",
        );
    }
    let counts = parse_u64_list(info_value(info, "RN"), line, "RN")?;
    if counts.len() != expected.alternates.len() {
        return vcf_error(
            line,
            "RN no trae exactamente un valor por alelo alternativo",
        );
    }
    for (observed, blocks) in counts.iter().zip(&expected.alternate_repeats) {
        if usize::try_from(*observed).ok() != Some(blocks.len()) {
            return vcf_error(
                line,
                "RN no cuenta las secuencias repetidas que la verdad da a ese alelo",
            );
        }
    }
    let motifs = info_value(info, "RUS")
        .unwrap_or_default()
        .split(',')
        .collect::<Vec<_>>();
    let motif_lengths = parse_u64_list(info_value(info, "RUL"), line, "RUL")?;
    let unit_counts = parse_f64_list(info_value(info, "RUC"), line, "RUC")?;
    let bases = parse_u64_list(info_value(info, "RB"), line, "RB")?;
    let total = expected_blocks.len();
    if motifs.len() != total
        || motif_lengths.len() != total
        || unit_counts.len() != total
        || bases.len() != total
    {
        return vcf_error(
            line,
            "las listas de repetición no suman los recuentos declarados en RN",
        );
    }
    for (index, block) in expected_blocks.iter().enumerate() {
        let motif_length = block.motif.len() as u64;
        if motifs[index] != block.motif {
            return vcf_error(line, "RUS no es el motivo que la verdad declara");
        }
        if motif_lengths[index] != motif_length {
            return vcf_error(line, "RUL no es la longitud del motivo que la acompaña");
        }
        if bases[index] != block.bases {
            return vcf_error(line, "RB no son las bases que la verdad declara");
        }
        let derived = block.bases as f64 / motif_length as f64;
        if (unit_counts[index] - derived).abs() > 1e-9 * derived.max(1.0) {
            return vcf_error(line, "RUC no son las bases entre la longitud del motivo");
        }
    }
    Ok(())
}

fn parse_u64_list(
    value: Option<&str>,
    line: u64,
    key: &str,
) -> Result<Vec<u64>, VcfSemanticValidationError> {
    value
        .unwrap_or_default()
        .split(',')
        .map(|item| item.parse::<u64>().ok())
        .collect::<Option<Vec<_>>>()
        .map_or_else(
            || vcf_error(line, format!("{key} trae un valor que no es un entero")),
            Ok,
        )
}

fn parse_f64_list(
    value: Option<&str>,
    line: u64,
    key: &str,
) -> Result<Vec<f64>, VcfSemanticValidationError> {
    value
        .unwrap_or_default()
        .split(',')
        .map(|item| item.parse::<f64>().ok().filter(|value| value.is_finite()))
        .collect::<Option<Vec<_>>>()
        .map_or_else(
            || vcf_error(line, format!("{key} trae un valor que no es un decimal")),
            Ok,
        )
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
// ---------------------------------------------------------------------------
// Las cinco constructoras de rechazo, y por qué son las únicas
//
// Todo rechazo que **este** verificador decide nace en una de ellas. No es una
// manía de estilo: son el punto por el que pasa el barrido de guardas de
// [`crate::guard_trace`], de modo que un rechazo construido a mano queda
// **fuera de la puerta** y podría romperse sin que nadie lo notase. Ocurrió:
// trece rechazos se construían a mano y el denominador que la puerta publicaba
// se quedaba corto, uno de ellos dentro de la propia función que ata la
// columna de referencia al genoma real.
//
// La prueba léxica `the_vcf_oracle_builds_every_rejection_through_its_constructors`
// de `tests/vcf_semantic_validation.rs` impide que el agujero se reabra.
//
// La única variante que no tiene constructora es `Truth`, y a propósito: no es
// un rechazo que decida este verificador, sino el veredicto del lector
// independiente de la verdad, que se reenvía tal cual con sus propias guardas.
// ---------------------------------------------------------------------------

/// Rechazo del contrato VCF, con la línea examinada.
///
/// La anotación `#[track_caller]` y la llamada a `record` son el barrido de
/// guardas: anotan **desde qué línea** se rechazó, no esta.
#[track_caller]
fn vcf_reject(line: u64, detail: impl Into<String>) -> VcfSemanticValidationError {
    crate::guard_trace::record();
    VcfSemanticValidationError::Vcf {
        line,
        detail: detail.into(),
    }
}

/// Rechazo del contrato VCF sin línea, para lo que no se juzga sobre una.
#[track_caller]
fn oracle_error(detail: impl Into<String>) -> VcfSemanticValidationError {
    vcf_reject(0, detail)
}

/// Rechazo del FASTA autoritativo.
#[track_caller]
fn fasta_reject(detail: impl Into<String>) -> VcfSemanticValidationError {
    crate::guard_trace::record();
    VcfSemanticValidationError::Fasta {
        detail: detail.into(),
    }
}

/// Rechazo de entrada y salida, con el fichero al que se refiere.
#[track_caller]
fn io_reject(path: impl Into<String>, detail: impl Into<String>) -> VcfSemanticValidationError {
    crate::guard_trace::record();
    VcfSemanticValidationError::Io {
        path: path.into(),
        detail: detail.into(),
    }
}

/// Rechazo por pedir bases de una referencia validada que no las conserva.
#[track_caller]
fn bases_not_materialized() -> VcfSemanticValidationError {
    crate::guard_trace::record();
    VcfSemanticValidationError::BasesNotMaterialized
}

/// El genotipo que la verdad exige para una muestra en un registro.
///
/// Sus dos rechazos —la muestra sin copias en el contig y el alelo alternativo
/// que no aparece entre los del registro— **no los ejercita ninguna prueba y no
/// se ha demostrado que sean inalcanzables**. Los dos exigirían una verdad en
/// la que un individuo no tuviera copias en una secuencia que otro sí publica,
/// o en la que un portador llevara un alelo que la composición no ha recogido,
/// y no he construido ninguna de las dos. Queda pendiente.
fn expected_gt(
    record: &ExpectedRecord,
    sample: &str,
) -> Result<String, VcfSemanticValidationError> {
    let copies = record
        .copies
        .get(sample)
        .ok_or_else(|| oracle_error("muestra VCF no tiene copias truth para el contig"))?;
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
                    .ok_or_else(|| oracle_error("ALT truth no aparece en el registro esperado"))
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

/// `GPS_CASE` es el **primer** campo de `INFO` y, tras él, solo pueden venir
/// los campos de repetición en tándem del estándar.
///
/// La comprobación era antes «`INFO` es exactamente `GPS_CASE`», que rechazaba
/// una clave repetida por el mero hecho de que sobraba. Al ampliarla a una
/// lista cerrada esa propiedad se perdió sin querer: una línea con dos veces el
/// mismo campo pasaba y el lector se quedaba con la primera aparición. Aquí se
/// exige otra vez, porque el estándar no define qué significa una clave
/// duplicada y quedarse con una de las dos es elegir en silencio. Cualquier
/// otra clave sigue siendo un rechazo, que es lo que impide que un productor
/// cuele un campo propio sin declararlo.
fn parse_gps_case(info: &str, line: u64) -> Result<Vec<String>, VcfSemanticValidationError> {
    let mut seen = BTreeSet::new();
    for field in info.split(';') {
        let key = field.split_once('=').map_or(field, |(key, _)| key);
        if !seen.insert(key) {
            return vcf_error(line, "INFO repite una clave");
        }
    }
    let mut fields = info.split(';');
    let Some(value) = fields
        .next()
        .and_then(|first| first.strip_prefix("GPS_CASE="))
    else {
        return vcf_error(line, "INFO debe declarar GPS_CASE en primer lugar");
    };
    if value.is_empty() || value.split(',').any(str::is_empty) {
        return vcf_error(
            line,
            "INFO debe contener GPS_CASE canónico sin valores vacíos",
        );
    }
    for field in fields {
        let key = field.split_once('=').map_or(field, |(key, _)| key);
        if !TANDEM_REPEAT_INFO_KEYS.contains(&key) {
            return vcf_error(
                line,
                "INFO solo admite GPS_CASE y los campos de repetición del estándar",
            );
        }
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
        .ok_or_else(|| vcf_reject(line, "CHROM no existe en el FASTA autoritativo"))?;
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
            return Err(fasta_reject("bases antes de la primera cabecera FASTA"));
        }
    }
    if let Some((header, bases)) = current {
        records.push(oracular_sequence(&header, bases, reference_id)?);
    }
    Reference::new(
        // Inalcanzable por la frontera con rutas, y demostrado ejecutándolo:
        // la identidad de referencia sale de la cabecera de la verdad, que se
        // lee y valida **antes** de abrir el FASTA, y el lector de la verdad
        // rechaza una identidad que no lo sea. Se conserva porque esta función
        // no se fía de quién le pasa el nombre.
        ReferenceId::parse(reference_id, "oráculo VCF")
            .map_err(|error| fasta_reject(format!("{error:?}")))?,
        ReferenceOrigin::Generated,
        records,
    )
    // Este otro sí se ejercita: un FASTA con dos secuencias del mismo nombre
    // produce piezas individualmente válidas y un conjunto que no lo es.
    .map_err(|error| fasta_reject(format!("{error:?}")))
}

/// Traduce un registro FASTA leído a una secuencia de referencia.
///
/// De sus cuatro rechazos, la batería ejercita dos —la cabecera ilegible y las
/// bases que no son ADN—. Los otros dos, el de la identidad de secuencia y el
/// de la construcción final, **no los ejercita ninguna prueba y no se ha
/// demostrado que sean inalcanzables**: no he conseguido una cabecera que el
/// analizador acepte y cuya identidad rechace después, ni una secuencia cuyas
/// tres piezas sean válidas y cuya construcción falle. Queda pendiente, y se
/// dice aquí en vez de fingir que están cubiertos.
fn oracular_sequence(
    header: &str,
    bases: Vec<u8>,
    reference_id: &str,
) -> Result<ReferenceSequence, VcfSemanticValidationError> {
    let fasta = FastaHeader::parse(header.as_bytes())
        .map_err(|error| fasta_reject(format!("{error:?}")))?;
    let id = SequenceId::parse(fasta.id_text(), reference_id)
        .map_err(|error| fasta_reject(format!("{error:?}")))?;
    let dna =
        DnaSequence::try_from_bytes(&bases).map_err(|error| fasta_reject(format!("{error:?}")))?;
    ReferenceSequence::new(id, fasta, SequenceClass::Chromosome, &dna)
        .map_err(|error| fasta_reject(format!("{error:?}")))
}

fn strictly_sorted(values: &[String]) -> bool {
    values.windows(2).all(|pair| pair[0] < pair[1])
}

#[track_caller]
fn io_error(path: &Path, error: std::io::Error) -> VcfSemanticValidationError {
    io_reject(path.display().to_string(), error.to_string())
}

/// Igual que [`oracle_error`], pero devolviendo ya el `Err` y conservando el
/// número de línea del VCF examinado. También participa en el barrido.
#[track_caller]
fn vcf_error<T>(line: u64, detail: impl Into<String>) -> Result<T, VcfSemanticValidationError> {
    Err(vcf_reject(line, detail))
}

#[cfg(test)]
mod tests {
    use super::{
        ExpectedRecord, OraclePath, OraclePiece, OracleSpan, VcfSemanticValidationError,
        expected_records, merge_oracle_spans, oracle_normalize_window, oracle_path,
        validate_oracle_phases, validate_reference_allele, validate_tandem_repeat_info,
    };
    use crate::identity::{ReferenceId, SequenceId};
    use crate::reference::{
        DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
    };
    use crate::truth_schema::{
        CopyIdentity, ModificationProvenance, ModificationProvenanceKind, SegmentSource,
        TruthCurrentPiece, TruthCurrentSnp, TruthDocument, TruthRecord, TruthRepeatBlock,
        TruthSourceKind, TruthSourceRef,
    };
    use crate::truth_validation::read_truth_jsonl;
    use std::collections::BTreeMap;
    use std::io::Cursor;

    // -----------------------------------------------------------------------
    // Guardas de **segunda opinión**
    //
    // Las comprobaciones que siguen no las alcanza ninguna corrupción del VCF
    // publicado, y no porque falten pruebas: **otra comprobación anterior caza
    // antes la discrepancia**. Son el segundo par de ojos del oráculo sobre su
    // propia reconstrucción, y el único modo honrado de ejercitarlas es
    // llamarlas directamente con una entrada que solo ellas puedan rechazar.
    // Cada prueba dice, en su documentación, quién caza antes.
    // -----------------------------------------------------------------------

    /// Referencia mínima de una sola secuencia `chr1`.
    fn tiny_reference(bases: &[u8]) -> Reference {
        let sequence = ReferenceSequence::new(
            SequenceId::parse("chr1", "prueba del oráculo VCF").expect("identidad de secuencia"),
            FastaHeader::parse(b"chr1").expect("cabecera FASTA"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(bases).expect("ADN"),
        )
        .expect("secuencia de referencia");
        Reference::new(
            ReferenceId::parse("REF_main", "prueba del oráculo VCF")
                .expect("identidad de referencia"),
            ReferenceOrigin::Generated,
            vec![sequence],
        )
        .expect("referencia")
    }

    fn source_ref() -> TruthSourceRef {
        TruthSourceRef {
            kind: TruthSourceKind::Reference,
            source_id: "REF_main".to_owned(),
        }
    }

    fn provenance() -> ModificationProvenance {
        ModificationProvenance {
            kind: ModificationProvenanceKind::Explicit,
            version: "1".to_owned(),
            source: source_ref(),
            matrix: None,
        }
    }

    fn span(start: u64, end: u64) -> SegmentSource {
        SegmentSource {
            source: source_ref(),
            reference_id: "REF_main".to_owned(),
            sequence_id: "chr1".to_owned(),
            start,
            end,
            sequence_length: 8,
        }
    }

    fn identity(copy_index: u8) -> CopyIdentity {
        CopyIdentity {
            reference_id: "REF_main".to_owned(),
            individual_id: "IND_01".to_owned(),
            sequence_id: "chr1".to_owned(),
            copy_index,
        }
    }

    /// El motivo exacto del rechazo, o el pánico de la prueba.
    fn detail(error: &VcfSemanticValidationError) -> String {
        match error {
            VcfSemanticValidationError::Vcf { detail, .. } => detail.clone(),
            other => panic!("se esperaba un rechazo del contrato VCF y llegó {other:?}"),
        }
    }

    /// **La columna REF está atada al genoma real, no a la derivación.**
    ///
    /// *Quién caza antes, y por eso esta prueba es unitaria.* De extremo a
    /// extremo, `validate_vcf_against_truth` compara la columna REF contra su
    /// propia composición —«REF o ALT no coincide con la derivación
    /// independiente»— **antes** de llamar aquí, y esa composición sale de
    /// rebanar el FASTA. Las dos comparaciones no pueden discrepar sin un
    /// defecto **dentro** del oráculo, que es justo lo que esta guarda vigila.
    /// Ninguna corrupción del VCF la alcanza; llamarla directamente sí.
    ///
    /// *Qué hace que pueda fallar:* el primer caso fija que la columna correcta
    /// se acepta —sin él, una guarda que rechazara todo pasaría el resto—, y
    /// los tres siguientes cambian una base, desplazan la posición y se salen
    /// del final de la secuencia.
    #[test]
    fn the_reference_column_is_tied_to_the_real_genome() {
        let reference = tiny_reference(b"ACGTACGT");
        validate_reference_allele(&reference, "chr1", 3, "GT", 7)
            .expect("la columna que el FASTA dice se acepta");

        for (reason, position, allele) in [
            ("una base cambiada", 3_u64, "GA"),
            ("la misma base una posición más allá", 4, "GT"),
            ("una columna que se sale del final", 8, "TT"),
        ] {
            let error = validate_reference_allele(&reference, "chr1", position, allele, 7)
                .expect_err(reason);
            assert_eq!(
                detail(&error),
                "REF no coincide con el FASTA autoritativo",
                "«{reason}» se rechazó por otra razón"
            );
        }

        let error = validate_reference_allele(&reference, "chr2", 3, "GT", 7)
            .expect_err("una secuencia que el FASTA no tiene");
        assert_eq!(detail(&error), "CHROM no existe en el FASTA autoritativo");
    }

    /// **Las fases que el oráculo recompone también tienen que ser contiguas.**
    ///
    /// *Quién caza antes.* `validate_declared_copy_phases` hace la misma
    /// comprobación sobre las fases **declaradas** y corre antes, de modo que
    /// ningún documento llega aquí con fases rotas. Esta es la comprobación
    /// sobre las fases que el oráculo ha vuelto a leer de la verdad: existe
    /// para que una reconstrucción que perdiera o duplicara una copia no
    /// pasara inadvertida.
    #[test]
    fn the_recomputed_copy_phases_must_be_contiguous() {
        let path = |copy_index| OraclePath {
            identity: identity(copy_index),
            sequence_id: "chr1".to_owned(),
            pieces: Vec::new(),
        };
        validate_oracle_phases(&[(0, path(0)), (1, path(1))]).expect("dos fases contiguas");

        for (reason, phases) in [
            ("un hueco entre las fases", vec![0_u8, 2]),
            ("una fase repetida", vec![0, 0]),
            ("una fase que no empieza en cero", vec![1, 2]),
        ] {
            let paths = phases
                .iter()
                .enumerate()
                .map(|(index, phase)| (*phase, path(u8::try_from(index).expect("índice pequeño"))))
                .collect::<Vec<_>>();
            let error = validate_oracle_phases(&paths).expect_err(reason);
            assert_eq!(
                detail(&error),
                "fases de copia no contiguas para IND_01/chr1",
                "«{reason}» se rechazó por otra razón"
            );
        }
    }

    /// **La reconstrucción independiente rechaza una verdad mal enhebrada.**
    ///
    /// *Quién caza antes.* `validate_and_reconstruct_truth` valida el documento
    /// entero contra el FASTA antes de que `expected_records` lo mire, así que
    /// ninguna de estas formas llega hasta aquí en una ejecución completa. El
    /// oráculo las recomprueba de todos modos porque **no se fía de quien
    /// escribió la verdad**, y esa desconfianza es su razón de ser.
    ///
    /// *Qué hace que pueda fallar:* el primer caso fija que la forma correcta
    /// se acepta y reproduce las bases esperadas.
    #[test]
    fn the_independent_path_rejects_a_badly_threaded_truth() {
        let reference = b"ACGTACGT";
        let good = oracle_path(
            identity(0),
            &[
                TruthCurrentPiece::Source {
                    source: span(0, 4),
                    snps: vec![TruthCurrentSnp {
                        event_id: "EVT_01".to_owned(),
                        position: 2,
                        reference: "G".to_owned(),
                        alternate: "T".to_owned(),
                        provenance: provenance(),
                    }],
                },
                TruthCurrentPiece::Inserted {
                    bases: "CC".to_owned(),
                    provenance: provenance(),
                },
                TruthCurrentPiece::Source {
                    source: span(4, 8),
                    snps: Vec::new(),
                },
            ],
            reference,
        )
        .expect("la forma correcta se acepta");
        assert_eq!(
            good.window(0, 8).expect("ventana entera"),
            b"ACTTCCACGT".to_vec(),
            "el SNP y la inserción tienen que aparecer en la secuencia reconstruida"
        );

        for (reason, expected_detail, pieces) in [
            (
                "un tramo source vacío",
                "piezas source truth no colineales",
                vec![TruthCurrentPiece::Source {
                    source: span(4, 4),
                    snps: Vec::new(),
                }],
            ),
            (
                "un tramo source que se sale del FASTA",
                "piezas source truth no colineales",
                vec![TruthCurrentPiece::Source {
                    source: span(0, 9),
                    snps: Vec::new(),
                }],
            ),
            (
                "dos tramos source que retroceden",
                "piezas source truth no colineales",
                vec![
                    TruthCurrentPiece::Source {
                        source: span(2, 6),
                        snps: Vec::new(),
                    },
                    TruthCurrentPiece::Source {
                        source: span(1, 4),
                        snps: Vec::new(),
                    },
                ],
            ),
            (
                "un SNP anterior al tramo que lo contiene",
                "SNP fuera de source",
                vec![TruthCurrentPiece::Source {
                    source: span(4, 8),
                    snps: vec![TruthCurrentSnp {
                        event_id: "EVT_01".to_owned(),
                        position: 1,
                        reference: "C".to_owned(),
                        alternate: "T".to_owned(),
                        provenance: provenance(),
                    }],
                }],
            ),
            (
                "un SNP posterior al tramo que lo contiene",
                "SNP fuera de source",
                vec![TruthCurrentPiece::Source {
                    source: span(0, 4),
                    snps: vec![TruthCurrentSnp {
                        event_id: "EVT_01".to_owned(),
                        position: 6,
                        reference: "C".to_owned(),
                        alternate: "T".to_owned(),
                        provenance: provenance(),
                    }],
                }],
            ),
            (
                "un tracto de repetición que retrocede",
                "piezas source truth no colineales",
                vec![
                    TruthCurrentPiece::Source {
                        source: span(0, 4),
                        snps: Vec::new(),
                    },
                    TruthCurrentPiece::RepeatAllele {
                        locus_id: "LOC_01".to_owned(),
                        source: span(2, 6),
                        bases: "CACACA".to_owned(),
                        blocks: vec![TruthRepeatBlock {
                            motif: "CA".to_owned(),
                            bases: 6,
                        }],
                        provenance: provenance(),
                    },
                ],
            ),
            (
                "una inserción sin ningún tramo delante",
                "inserción truth sin source precedente",
                vec![TruthCurrentPiece::Inserted {
                    bases: "CC".to_owned(),
                    provenance: provenance(),
                }],
            ),
        ] {
            let error = oracle_path(identity(0), &pieces, reference).expect_err(reason);
            assert_eq!(
                detail(&error),
                expected_detail,
                "«{reason}» se rechazó por otra razón"
            );
        }
    }

    /// **Una ventana no puede partir un tracto de repetición realizado.**
    ///
    /// *Quién caza antes, y por qué esta prueba es unitaria.* Las ventanas que
    /// el oráculo se pide a sí mismo salen de `changes`, que empuja el tracto
    /// **entero** como un solo tramo, y de fusiones que solo pueden agrandarlo.
    /// Ninguna ejecución pide, pues, media ventana de un tracto. La guarda
    /// existe porque leer la mitad de un tracto realizado daría bases que no
    /// son de nadie, y se ejercita llamando a `window` con la ventana que la
    /// producción no construye.
    #[test]
    fn a_window_may_not_split_a_realized_repeat_tract() {
        let path = OraclePath {
            identity: identity(0),
            sequence_id: "chr1".to_owned(),
            pieces: vec![OraclePiece::Replaced {
                start: 2,
                end: 6,
                bases: b"CACACA".to_vec(),
                blocks: vec![TruthRepeatBlock {
                    motif: "CA".to_owned(),
                    bases: 6,
                }],
            }],
        };
        assert_eq!(
            path.window(2, 6).expect("el tracto entero se lee"),
            b"CACACA".to_vec()
        );
        assert_eq!(
            path.window(0, 8)
                .expect("una ventana que lo envuelve también"),
            b"CACACA".to_vec()
        );
        for (reason, start, end) in [
            ("una ventana que empieza dentro", 3_u64, 6_u64),
            ("una ventana que termina dentro", 2, 5),
        ] {
            let error = path.window(start, end).expect_err(reason);
            assert_eq!(
                detail(&error),
                "la ventana parte un tracto de repetición",
                "«{reason}» se rechazó por otra razón"
            );
        }
    }

    /// **Una inserción pegada al primer par de bases no tiene ancla izquierda.**
    ///
    /// *Quién caza antes.* La producción no publica una inserción antes de la
    /// primera base: el formato VCF exige un ancla, y la admisión la rechaza
    /// antes de simular. La guarda existe porque el oráculo lee documentos que
    /// no ha escrito él, y sin ella el cálculo del ancla restaría por debajo de
    /// cero. Se ejercita normalizando la ventana vacía en la posición cero.
    #[test]
    fn an_insertion_at_the_very_start_has_no_left_anchor() {
        let reference = b"ACGTACGT";
        let inserted_at = |anchor| OraclePath {
            identity: identity(0),
            sequence_id: "chr1".to_owned(),
            pieces: vec![
                OraclePiece::Source {
                    start: 0,
                    end: 4,
                    bases: b"ACGT".to_vec(),
                    snps: Vec::new(),
                },
                OraclePiece::Inserted {
                    anchor,
                    bases: b"TT".to_vec(),
                },
            ],
        };
        let path = inserted_at(0);
        let error = oracle_normalize_window(&path, reference, OracleSpan { start: 0, end: 0 })
            .expect_err("una inserción en la posición cero no tiene ancla");
        assert_eq!(detail(&error), "indel sin ancla izquierda");

        // Control positivo: la misma inserción una base más allá sí normaliza,
        // y lo hace anclando en la base anterior.
        let anchored = inserted_at(1);
        let allele = oracle_normalize_window(&anchored, reference, OracleSpan { start: 1, end: 1 })
            .expect("con ancla sí normaliza")
            .expect("hay alelo");
        assert_eq!(allele.position, 0);
        assert_eq!(allele.reference, "A");
        assert_eq!(allele.alternate, "ATT");
    }

    fn block(motif: &str, bases: u64) -> TruthRepeatBlock {
        TruthRepeatBlock {
            motif: motif.to_owned(),
            bases,
        }
    }

    /// Registro con **recuentos distintos por alelo**, incluido uno sin
    /// tracto: dos bloques, ninguno y uno.
    fn heterogeneous_record() -> ExpectedRecord {
        ExpectedRecord {
            reference: "A".to_owned(),
            alternates: vec!["AA".to_owned(), "AC".to_owned(), "AG".to_owned()],
            gps_cases: Vec::new(),
            copies: BTreeMap::new(),
            alternate_repeats: vec![
                vec![block("CA", 10), block("GGG", 9)],
                Vec::new(),
                vec![block("TTTT", 12)],
            ],
        }
    }

    /// **Recuentos distintos entre alelos, que la producción no construye.**
    ///
    /// Los alelos de un locus comparten su lista de motivos, de modo que todos
    /// traen el mismo número de bloques: el camino en el que un alelo tiene
    /// dos, otro ninguno y otro uno no lo recorre ninguna prueba de extremo a
    /// extremo. Aquí se alimenta a mano. La primera aserción fija que la forma
    /// correcta se acepta —sin ella, un oráculo que rechazara todo pasaría las
    /// demás—.
    #[test]
    fn slicing_the_lists_per_allele_catches_a_shift_with_heterogeneous_counts() {
        let expected = heterogeneous_record();
        let correct = "GPS_CASE=C;RN=2,0,1;RUS=CA,GGG,TTTT;RUL=2,3,4;RUC=5,3,3;RB=10,9,12";
        validate_tandem_repeat_info(&expected, correct, 1).expect("la forma correcta se acepta");

        for (reason, info) in [
            (
                "recuentos desplazados entre alelos",
                "GPS_CASE=C;RN=1,0,2;RUS=CA,GGG,TTTT;RUL=2,3,4;RUC=5,3,3;RB=10,9,12",
            ),
            (
                "bases desplazadas dentro de la lista aplanada",
                "GPS_CASE=C;RN=2,0,1;RUS=CA,GGG,TTTT;RUL=2,3,4;RUC=5,3,3;RB=9,12,10",
            ),
            (
                "un alelo sin tracto al que se le da uno",
                "GPS_CASE=C;RN=2,1,1;RUS=CA,GGG,TTTT,TTTT;RUL=2,3,4,4;RUC=5,3,3,3;RB=10,9,12,12",
            ),
        ] {
            assert!(
                validate_tandem_repeat_info(&expected, info, 1).is_err(),
                "el oráculo tiene que rechazar «{reason}»"
            );
        }
    }

    /// La verdad mínima de F2 que el corpus independiente ya usa.
    fn fixture_truth() -> TruthDocument {
        read_truth_jsonl(&mut Cursor::new(include_str!(
            "../tests/fixtures/truth-v1/valid.jsonl"
        )))
        .expect("verdad de referencia")
    }

    /// **La composición de los registros esperados no se fía de la verdad.**
    ///
    /// *Quién caza antes.* `validate_and_reconstruct_truth` corre antes y
    /// rechazaría una verdad que nombrase una secuencia ausente, de modo que
    /// ninguna ejecución completa llega aquí con ella. La segunda mitad —un
    /// alelo alternativo que ningún caso de la verdad respalda— tampoco se
    /// alcanza desde fuera: la verdad publicada siempre trae sus genotipos.
    ///
    /// *Qué hace que pueda fallar:* el control positivo compone los registros
    /// de la verdad íntegra, y cada caso exige el motivo exacto.
    #[test]
    fn the_expected_records_do_not_trust_the_truth_they_read() {
        let reference = tiny_reference(b"AAAAAAAA");
        let truth = fixture_truth();
        assert!(
            !expected_records(&reference, &truth)
                .expect("la verdad íntegra compone sus registros")
                .is_empty(),
            "el control positivo tiene que producir algún registro"
        );

        let mut renamed = truth.clone();
        for record in &mut renamed.records {
            if let TruthRecord::HaplotypeCopy(copy) = record {
                copy.identity.sequence_id = "chr9".to_owned();
            }
        }
        let error = expected_records(&reference, &renamed)
            .expect_err("una copia que nombra una secuencia ausente");
        assert_eq!(
            detail(&error),
            "copia truth nombra una secuencia FASTA ausente"
        );

        let mut without_genotypes = truth.clone();
        without_genotypes
            .records
            .retain(|record| matches!(record, TruthRecord::HaplotypeCopy(_)));
        let error = expected_records(&reference, &without_genotypes)
            .expect_err("un ALT que ningún caso respalda");
        assert_eq!(detail(&error), "ALT compuesto sin caso truth fuente");

        // Dos copias que realizan el **mismo** tracto y lo descomponen en
        // bloques distintos: el ALT coincide y la estructura declarada no.
        let mut contradictory = truth.clone();
        let mut blocks = [
            vec![TruthRepeatBlock {
                motif: "AC".to_owned(),
                bases: 4,
            }],
            vec![TruthRepeatBlock {
                motif: "ACAC".to_owned(),
                bases: 4,
            }],
        ]
        .into_iter();
        for record in &mut contradictory.records {
            if let TruthRecord::HaplotypeCopy(copy) = record {
                copy.current_pieces = vec![
                    TruthCurrentPiece::Source {
                        source: span(0, 2),
                        snps: Vec::new(),
                    },
                    TruthCurrentPiece::RepeatAllele {
                        locus_id: "LOC_01".to_owned(),
                        source: span(2, 6),
                        bases: "ACAC".to_owned(),
                        blocks: blocks.next().expect("dos copias, dos descomposiciones"),
                        provenance: provenance(),
                    },
                    TruthCurrentPiece::Source {
                        source: span(6, 8),
                        snps: Vec::new(),
                    },
                ];
            }
        }
        let error = expected_records(&reference, &contradictory)
            .expect_err("el mismo ALT con dos descomposiciones");
        assert_eq!(
            detail(&error),
            "dos copias dan el mismo ALT con tractos distintos"
        );

        // La **misma copia** declarada dos veces con contenidos que no
        // coinciden: dos alelos distintos en el mismo CHROM+POS para una sola
        // identidad. Las fases siguen siendo contiguas, de modo que la
        // comprobación de fases no lo caza antes.
        let mut duplicated = truth.clone();
        let mut twin = duplicated
            .records
            .iter()
            .find_map(|record| match record {
                TruthRecord::HaplotypeCopy(copy) if copy.phase == 0 => Some(copy.clone()),
                _ => None,
            })
            .expect("la verdad trae la copia de fase cero");
        twin.phase = 1;
        for piece in &mut twin.current_pieces {
            if let TruthCurrentPiece::Source { snps, .. } = piece {
                for snp in snps {
                    snp.alternate = "T".to_owned();
                }
            }
        }
        duplicated.records.retain(
            |record| !matches!(record, TruthRecord::HaplotypeCopy(copy) if copy.phase == 1),
        );
        duplicated.records.push(TruthRecord::HaplotypeCopy(twin));
        let error = expected_records(&reference, &duplicated)
            .expect_err("una copia declarada dos veces con contenidos distintos");
        assert_eq!(
            detail(&error),
            "una copia produce ALT incompatibles al mismo CHROM+POS"
        );
    }

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
