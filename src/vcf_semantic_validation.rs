//! Oráculo semántico independiente para la vista VCF de verdad F2.
//!
//! Esta frontera lee los tres artefactos publicados. No llama a
//! `vcf_projection`, `vcf_grouping`, `vcf_representation` ni al escritor VCF:
//! la verdad se valida y se interpreta antes de abrir el VCF.

use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;

use sha2::{Digest, Sha256};
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
        let mut file = File::open(vcf_path).map_err(|error| io_error(vcf_path, error))?;
        // Un bloque que declara más de 64 KiB haría entrar en pánico al lector
        // (DT-030): se rechaza antes de dárselo.
        crate::bgzf_blocks::reject_oversized_bgzf_blocks(&mut file)
            .map_err(|fault| io_reject(vcf_path.display().to_string(), fault.to_string()))?;
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
    let mut covered = None::<(String, u64, u64)>;
    let mut rebuilt = rebuild_states(reference, truth, &expected_samples);
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
        // **Dos registros que se pisan son un error del productor** (DTC-064),
        // igual que dos con el mismo `CHROM+POS`. Se juzga sobre el fichero,
        // antes de buscar el registro esperado, y sin pasar por la fusión de
        // este oráculo: así también lo caza si productor y oráculo dejaran de
        // fundir a la vez. Un registro con el mismo `POS` o anterior lo juzgan
        // las guardas de duplicado y de orden de más abajo.
        let reference_end = (position - 1).saturating_add(fields[3].len() as u64);
        if let Some((chrom, previous_position, previous_end)) = &covered
            && chrom == fields[0]
            && position > *previous_position
            && position - 1 < *previous_end
        {
            return vcf_error(
                line_number,
                "el intervalo REF del registro pisa el de un registro anterior",
            );
        }
        covered = match covered.take() {
            Some((chrom, previous_position, previous_end)) if chrom == fields[0] => Some((
                chrom,
                previous_position.max(position),
                previous_end.max(reference_end),
            )),
            _ => Some((fields[0].to_owned(), position, reference_end)),
        };
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
        // La reconstrucción de cada copia avanza con lo que el registro
        // **dice**, antes de compararlo con lo esperado: es la comprobación
        // que no depende de que este oráculo derive bien los registros.
        if let Some(sequence) = rebuilt.get_mut(fields[0]) {
            // El REF ya se ató al FASTA unas líneas más arriba, de modo que la
            // posición cabe en la secuencia; y los registros llegan ordenados
            // y sin pisarse —lo exigen las guardas de arriba—, de modo que el
            // hueco desde el anterior nunca es negativo.
            let start = (position - 1) as usize;
            let gap = &sequence.bases[sequence.cursor..start];
            let alternates = fields[4].split(',').collect::<Vec<_>>();
            for copy in &mut sequence.copies {
                let Some(text) = fields[9 + copy.sample].split('|').nth(copy.phase) else {
                    return vcf_error(
                        line_number,
                        "GT no trae un alelo para cada fase de la copia",
                    );
                };
                let Ok(index) = text.parse::<usize>() else {
                    return vcf_error(line_number, "GT trae un alelo que no es un índice");
                };
                let allele = if index == 0 {
                    fields[3]
                } else if let Some(alternate) = alternates.get(index - 1) {
                    alternate
                } else {
                    return vcf_error(line_number, "GT nombra un ALT que el registro no tiene");
                };
                copy.hasher.update(gap);
                copy.hasher.update(allele.as_bytes());
                copy.length += (gap.len() + allele.len()) as u64;
            }
            sequence.cursor = start + fields[3].len();
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
    compare_rebuilt_copies(
        &truth_copy_digests(reference, truth),
        finish_rebuild(rebuilt),
    )
}

// ---------------------------------------------------------------------------
// La reconstrucción de cada copia desde el VCF
//
// Todo lo anterior compara el VCF con los registros que **este oráculo**
// deriva de la verdad. Si el oráculo y el productor compartieran un error de
// derivación —ocurrió dos veces: una inserción contada dos veces en el borde
// de una ventana, y una copia anotada con la referencia en una base que había
// borrado—, el VCF falso coincidiría con lo esperado y pasaría. Esto no
// deriva nada: aplica el `GT` de cada registro sobre la referencia, rehace
// cada copia y la compara con la que describen las piezas de `truth.jsonl`.
// Se compara por SHA-256 y longitud para no tener dos copias de cada genoma
// en memoria.
// ---------------------------------------------------------------------------

struct RebuiltCopy {
    sample: usize,
    phase: usize,
    individual: String,
    hasher: Sha256,
    length: u64,
}

struct RebuiltSequence<'reference> {
    bases: &'reference [u8],
    cursor: usize,
    copies: Vec<RebuiltCopy>,
}

/// Las bases de cada secuencia del FASTA. Una referencia sin bases ya la ha
/// rechazado la derivación de lo esperado antes de llegar aquí.
fn sequence_bases(reference: &Reference) -> BTreeMap<&str, &[u8]> {
    reference
        .sequences()
        .iter()
        .filter_map(|sequence| {
            sequence
                .sequence()
                .as_bytes()
                .ok()
                .map(|bases| (sequence.id().as_str(), bases))
        })
        .collect()
}

/// El estado inicial de la reconstrucción: una copia por individuo y fase en
/// cada secuencia, en el orden de las muestras de la cabecera.
fn rebuild_states<'reference>(
    reference: &'reference Reference,
    truth: &TruthDocument,
    samples: &[String],
) -> BTreeMap<String, RebuiltSequence<'reference>> {
    let mut phases = BTreeMap::<(String, String), usize>::new();
    for record in &truth.records {
        if let TruthRecord::HaplotypeCopy(copy) = record {
            *phases
                .entry((
                    copy.identity.sequence_id.clone(),
                    copy.identity.individual_id.clone(),
                ))
                .or_default() += 1;
        }
    }
    sequence_bases(reference)
        .into_iter()
        .map(|(id, bases)| {
            let id = id.to_owned();
            let copies = samples
                .iter()
                .enumerate()
                .flat_map(|(sample, individual)| {
                    let count = phases
                        .get(&(id.clone(), individual.clone()))
                        .copied()
                        .unwrap_or(0);
                    (0..count).map(move |phase| RebuiltCopy {
                        sample,
                        phase,
                        individual: individual.clone(),
                        hasher: Sha256::new(),
                        length: 0,
                    })
                })
                .collect();
            (
                id,
                RebuiltSequence {
                    bases,
                    cursor: 0,
                    copies,
                },
            )
        })
        .collect()
}

type CopyDigests = BTreeMap<(String, String, usize), ([u8; 32], u64)>;

/// Cierra la reconstrucción con la cola de cada secuencia tras su último
/// registro.
fn finish_rebuild(rebuilt: BTreeMap<String, RebuiltSequence<'_>>) -> CopyDigests {
    let mut digests = CopyDigests::new();
    for (sequence_id, state) in rebuilt {
        let tail = &state.bases[state.cursor..];
        for mut copy in state.copies {
            copy.hasher.update(tail);
            copy.length += tail.len() as u64;
            digests.insert(
                (sequence_id.clone(), copy.individual, copy.phase),
                (copy.hasher.finalize().into(), copy.length),
            );
        }
    }
    digests
}

/// Cada copia tal como la describen las piezas de `truth.jsonl`.
fn truth_copy_digests(reference: &Reference, truth: &TruthDocument) -> CopyDigests {
    let all_bases = sequence_bases(reference);
    let mut digests = CopyDigests::new();
    for record in &truth.records {
        let TruthRecord::HaplotypeCopy(copy) = record else {
            continue;
        };
        let Some(bases) = all_bases.get(copy.identity.sequence_id.as_str()) else {
            continue;
        };
        let mut hasher = Sha256::new();
        let mut length = 0_u64;
        for piece in &copy.current_pieces {
            match piece {
                TruthCurrentPiece::Source { source, snps } => {
                    // La verdad ya se validó y reconstruyó entera antes de
                    // abrir el VCF: los intervalos y los cambios de base caben
                    // en la referencia.
                    let mut cursor = source.start as usize;
                    let mut ordered = snps.iter().collect::<Vec<_>>();
                    ordered.sort_by_key(|snp| snp.position);
                    for snp in ordered {
                        let at = snp.position as usize;
                        hasher.update(&bases[cursor..at]);
                        hasher.update(snp.alternate.as_bytes());
                        cursor = at + 1;
                    }
                    hasher.update(&bases[cursor..source.end as usize]);
                    length += source.end - source.start;
                }
                TruthCurrentPiece::Inserted {
                    bases: inserted, ..
                } => {
                    hasher.update(inserted.as_bytes());
                    length += inserted.len() as u64;
                }
                TruthCurrentPiece::RepeatAllele { bases: tract, .. } => {
                    hasher.update(tract.as_bytes());
                    length += tract.len() as u64;
                }
            }
        }
        digests.insert(
            (
                copy.identity.sequence_id.clone(),
                copy.identity.individual_id.clone(),
                usize::from(copy.phase),
            ),
            (hasher.finalize().into(), length),
        );
    }
    digests
}

/// **El VCF tiene que describir exactamente cada copia de la verdad.**
fn compare_rebuilt_copies(
    expected: &CopyDigests,
    rebuilt: CopyDigests,
) -> Result<(), VcfSemanticValidationError> {
    for (copy, digest) in expected {
        if rebuilt.get(copy) != Some(digest) {
            return vcf_error(
                0,
                format!(
                    "la copia {} fase {} de {} que describe el VCF no es la de la verdad",
                    copy.1, copy.2, copy.0
                ),
            );
        }
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
    // Los cambios de cada copia, aunque no tenga ninguno: si ninguno toca el
    // locus de un caso, la copia es la referencia ahí.
    let mut copy_changes = BTreeMap::<CopyIdentity, Vec<OracleSpan>>::new();
    for (_, path) in &paths {
        let path_changes = path.changes();
        copy_changes.insert(path.identity.clone(), path_changes.clone());
        spans
            .entry(path.sequence_id.clone())
            .or_default()
            .extend(path_changes);
    }
    let mut alleles = Vec::new();
    // Ventanas en que cada copia queda igual a la referencia: las únicas en
    // que un caso puede no aportar ALT (excepción de REQ-VER-25).
    let mut unchanged = BTreeMap::<CopyIdentity, Vec<OracleSpan>>::new();
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
        let sequence_bases = sequence
            .sequence()
            .as_bytes()
            .map_err(|_| bases_not_materialized())?;
        let sequence_paths = paths
            .iter()
            .filter(|(_, path)| path.sequence_id == sequence_id)
            .map(|(_, path)| path)
            .collect::<Vec<_>>();
        // Los tractos que alguna copia realiza, con sus bloques, en el orden
        // de la secuencia. Son los únicos que tienen registro de repetición.
        let realized = OracleTracts::new(&sequence_paths);
        let mut sequence_alleles = Vec::new();
        for path in &sequence_paths {
            let (path_alleles, path_unchanged) =
                normalize_oracle_components(path, sequence_bases, &components, &realized)?;
            sequence_alleles.extend(path_alleles);
            unchanged
                .entry(path.identity.clone())
                .or_default()
                .extend(path_unchanged);
        }
        alleles.extend(oracle_fuse(sequence_alleles, sequence_bases, &realized)?);
    }
    every_case_has_an_alternate_or_an_unchanged_copy(truth, &alleles, &unchanged, &copy_changes)?;
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
                // Dos copias con el mismo ALT traen siempre la misma
                // anotación: la reconstrucción del registro ya da a las dos la
                // derivada de sus bases cuando la verdad les declara otras
                // (DTC-064), de modo que aquí no queda nada que decidir.
                repeats_by_alternate
                    .entry(alternate.clone())
                    .or_insert_with(|| allele.repeat_blocks.clone());
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

/// **Cada portador de un caso aporta un ALT, o su copia queda igual a la
/// referencia alrededor del locus.**
///
/// Un portador cuyo caso no cae en el tramo de ningún alelo de su copia solo
/// es legítimo si su efecto neto es nulo (excepción de REQ-VER-25), es decir,
/// si su copia es la referencia alrededor del locus. Hay dos maneras de
/// saberlo:
///
/// - en esa copia hay una ventana, de las que la normalización compuso, que
///   contiene el locus entero y cuya secuencia es la de la referencia;
/// - ningún cambio de la copia —hueco, punto de inserción, cambio de base o
///   tracto sustituido, sacados de sus piezas finales— toca el locus, extremos
///   incluidos. Es el caso de una inserción que una deleción posterior borra
///   entera: la copia queda intacta, no produce ningún cambio y, si ninguna
///   otra copia cambia ahí, tampoco se compone ninguna ventana.
///
/// Cualquier otro es un ALT perdido, y sin esta guarda el oráculo lo
/// admitiría en silencio: los registros esperados se derivan de las copias,
/// no de los casos.
///
/// *Quién caza antes.* Con una verdad válida no se alcanza de extremo a
/// extremo: el juez de la verdad exige que el locus de cada caso sea el que su
/// evento ocupa en la copia, y ese locus queda siempre dentro de un componente
/// de cambios de la copia, que o produce un alelo cuyo tramo lo contiene o es
/// una ventana igual a la referencia. La ejercita la prueba unitaria
/// `a_case_without_an_alternate_needs_an_unchanged_copy`, que llama a
/// [`expected_records`] con una verdad sin validar.
fn every_case_has_an_alternate_or_an_unchanged_copy(
    truth: &TruthDocument,
    alleles: &[OracleAllele],
    unchanged: &BTreeMap<CopyIdentity, Vec<OracleSpan>>,
    copy_changes: &BTreeMap<CopyIdentity, Vec<OracleSpan>>,
) -> Result<(), VcfSemanticValidationError> {
    for genotype in truth.records.iter().filter_map(|record| match record {
        TruthRecord::Genotype(genotype) => Some(genotype),
        TruthRecord::HaplotypeCopy(_) => None,
    }) {
        let locus = &genotype.locus;
        for carrier in genotype
            .carriers
            .iter()
            .filter(|carrier| carrier.allele_index != 0)
        {
            let contributes = alleles.iter().any(|allele| {
                allele.copy == carrier.copy
                    && allele.sequence_id == locus.sequence_id
                    && oracle_span_contains(allele.source, locus.start, locus.end)
            });
            let in_window = unchanged.get(&carrier.copy).is_some_and(|windows| {
                windows
                    .iter()
                    .any(|window| window.start <= locus.start && locus.end <= window.end)
            });
            let untouched = copy_changes.get(&carrier.copy).is_some_and(|changes| {
                !changes
                    .iter()
                    .any(|change| change.start <= locus.end && locus.start <= change.end)
            });
            let unchanged_copy =
                carrier.copy.sequence_id == locus.sequence_id && (in_window || untouched);
            if !contributes && !unchanged_copy {
                return Err(oracle_error(
                    "un caso sin ALT cuya copia no queda igual a la referencia",
                ));
            }
        }
    }
    Ok(())
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
        self.window_from(start, end, start)
    }

    /// La ventana, copiando solo las inserciones con punto en
    /// `[insertions_from, end]`: al anclar, la base que la ventana gana por la
    /// izquierda no trae consigo la inserción de ese punto, que es de la
    /// ventana de su propio componente (DTC-064).
    fn window_from(
        &self,
        start: u64,
        end: u64,
        insertions_from: u64,
    ) -> Result<Vec<u8>, VcfSemanticValidationError> {
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
                OraclePiece::Inserted { anchor, bases }
                    if insertions_from <= *anchor && *anchor <= end =>
                {
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
    tracts: &OracleTracts<'_>,
) -> Result<(Vec<OracleAllele>, Vec<OracleSpan>), VcfSemanticValidationError> {
    let mut components = initial.to_vec();
    loop {
        let normalized = components
            .iter()
            .copied()
            .map(|span| {
                // DTC-064: el alelo se alinea libre; solo si acabaría **dentro**
                // de un tracto realizado se repite con el límite en su última
                // base. El que cruza el tracto entero lo cruza. «Dentro» se
                // mide en la primera base que cambia, la siguiente al ancla.
                let free = oracle_normalize_window(path, reference, span, 0)?;
                let mut floor = 0_u64;
                if let Some(allele) = &free {
                    for (tract, _) in tracts.ending_after(allele.position, allele.position + 1) {
                        if tract.end <= span.start
                            && tract.start <= allele.position
                            && allele.position + 1 < tract.end
                        {
                            floor = floor.max(tract.end - 1);
                        }
                    }
                }
                if floor == 0 {
                    Ok((span, free))
                } else {
                    oracle_normalize_window(path, reference, span, floor).map(|value| (span, value))
                }
            })
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
                // Solo choca con un componente en el que **esta** copia tiene
                // alelo: un tramo que la copia conserva como la referencia no
                // la obliga a nada. Importa desde DTC-064, porque un alelo
                // anclado en la última base de un tracto que otra copia
                // realiza pisa ese componente sin que esta copia lo cambie, y
                // fundirlo aquí lo dejaría alinearse hacia dentro del tracto.
                // El componente atravesado al alinear sí se absorbe siempre.
                // Y choca también con el alelo de otro componente de la copia
                // cuyo `REF` pisa el suyo, aunque los componentes no se toquen.
                let other_reference_overlaps =
                    normalized[other_index]
                        .1
                        .as_ref()
                        .is_some_and(|other_allele| {
                            other_allele.position < end
                                && allele.position
                                    < other_allele.position + other_allele.reference.len() as u64
                        });
                if (normalized[other_index].1.is_some()
                    && allele.position < other.end
                    && other.start < end)
                    || other_reference_overlaps
                    || (traversed.start <= other.start && other.end <= traversed.end)
                {
                    merge.push((*span, *other));
                }
            }
        }
        if merge.is_empty() {
            // Un componente sin alelo es una ventana en que la copia queda
            // igual a la referencia; se devuelve aparte.
            let mut alleles = Vec::new();
            let mut unchanged = Vec::new();
            for (span, allele) in normalized {
                match allele {
                    Some(mut allele) => {
                        allele.source = span;
                        alleles.push(allele);
                    }
                    None => unchanged.push(span),
                }
            }
            return Ok((alleles, unchanged));
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

/// Los tractos que realizan las copias de una secuencia, **ordenados por
/// intervalo**, para que cada componente y cada borrador miren solo los que
/// pueden tocarlos (F6.4.a.8b). Antes, cada uno recorría todos los tractos de
/// todas las copias: borradores × copias × loci. Escrito aquí por su cuenta,
/// sin el índice del productor.
struct OracleTracts<'a> {
    /// `(intervalo, bloques)` de cada tracto realizado, en orden estable de
    /// intervalo: los del mismo intervalo conservan el orden de las copias.
    spans: Vec<(OracleSpan, &'a Vec<TruthRepeatBlock>)>,
    /// `reach[i]` es el mayor fin de `spans[..=i]`.
    reach: Vec<u64>,
    /// La primera copia con cada identidad, por su posición en la lista.
    paths: BTreeMap<&'a CopyIdentity, usize>,
    /// Los bloques del primer tracto de cada copia con cada intervalo.
    own: BTreeMap<(usize, OracleSpan), &'a Vec<TruthRepeatBlock>>,
}

impl<'a> OracleTracts<'a> {
    fn new(paths: &[&'a OraclePath]) -> Self {
        let mut spans = Vec::new();
        let mut own = BTreeMap::new();
        let mut by_identity = BTreeMap::new();
        for (index, path) in paths.iter().enumerate() {
            by_identity.entry(&path.identity).or_insert(index);
            for piece in &path.pieces {
                if let OraclePiece::Replaced {
                    start, end, blocks, ..
                } = piece
                {
                    let span = OracleSpan {
                        start: *start,
                        end: *end,
                    };
                    own.entry((index, span)).or_insert(blocks);
                    spans.push((span, blocks));
                }
            }
        }
        Self::indexed(spans, by_identity, own)
    }

    /// Ordena los tractos y calcula el fin acumulado.
    fn indexed(
        mut spans: Vec<(OracleSpan, &'a Vec<TruthRepeatBlock>)>,
        paths: BTreeMap<&'a CopyIdentity, usize>,
        own: BTreeMap<(usize, OracleSpan), &'a Vec<TruthRepeatBlock>>,
    ) -> Self {
        spans.sort_by_key(|(span, _)| *span);
        let mut furthest = 0;
        let reach = spans
            .iter()
            .map(|(span, _)| {
                furthest = span.end.max(furthest);
                furthest
            })
            .collect();
        Self {
            spans,
            reach,
            paths,
            own,
        }
    }

    /// Los que empiezan en `[from, to]`.
    fn starting_in(&self, from: u64, to: u64) -> &[(OracleSpan, &'a Vec<TruthRepeatBlock>)] {
        let first = self.spans.partition_point(|(span, _)| span.start < from);
        let last = self.spans.partition_point(|(span, _)| span.start <= to);
        self.spans.get(first..last).unwrap_or_default()
    }

    /// Los que empiezan antes de `to` y no están entre los que, todos, acaban
    /// en `from` o antes: un superconjunto de los que pisan `[from, to)`.
    fn reaching(&self, from: u64, to: u64) -> &[(OracleSpan, &'a Vec<TruthRepeatBlock>)] {
        let first = self.reach.partition_point(|reach| *reach <= from);
        let last = self.spans.partition_point(|(span, _)| span.start < to);
        self.spans.get(first..last).unwrap_or_default()
    }

    /// Los que empiezan en `position` o antes y no están entre los que,
    /// todos, acaban en `after` o antes: un superconjunto de los que contienen
    /// `position` y acaban después de `after`.
    fn ending_after(
        &self,
        position: u64,
        after: u64,
    ) -> &[(OracleSpan, &'a Vec<TruthRepeatBlock>)] {
        let first = self.reach.partition_point(|reach| *reach <= after);
        let last = self
            .spans
            .partition_point(|(span, _)| span.start <= position);
        self.spans.get(first..last).unwrap_or_default()
    }
}

/// Intervalo `[POS, POS+len(REF))` de un alelo del oráculo, ancla incluida.
fn oracle_vcf_span(allele: &OracleAllele) -> OracleSpan {
    OracleSpan {
        start: allele.position,
        end: allele.position + allele.reference.len() as u64,
    }
}

/// Funde los alelos de copias distintas que se pisan, **por su cuenta**.
///
/// Es la regla de DTC-064 escrita otra vez, sin mirar la del productor: dos
/// alelos se pisan si sus intervalos VCF, ancla incluida, comparten una base,
/// y se funden por transitividad. Un borrador cuyos alelos comparten `POS` se
/// deja como está, porque el agrupador por `CHROM+POS` ya los integra. Un
/// borrador con varias posiciones se reconstruye sobre la unión de sus
/// intervalos, aplicando a la referencia los alelos de cada copia, y después
/// se normaliza como un registro entero salvo que lleve un tracto. Se repite
/// hasta que ningún borrador pisa a otro, porque la alineación puede llevar un
/// registro fundido encima de su vecino.
fn oracle_fuse(
    alleles: Vec<OracleAllele>,
    reference: &[u8],
    realized: &OracleTracts<'_>,
) -> Result<Vec<OracleAllele>, VcfSemanticValidationError> {
    // Borradores iniciales: los alelos cuyos intervalos se encadenan.
    let mut sorted = alleles;
    sorted.sort_by_key(|allele| (allele.position, oracle_vcf_span(allele).end));
    let mut drafts = Vec::<Vec<OracleAllele>>::new();
    let mut reach = 0_u64;
    for allele in sorted {
        let span = oracle_vcf_span(&allele);
        match drafts.last_mut() {
            Some(draft) if span.start < reach => {
                reach = reach.max(span.end);
                draft.push(allele);
            }
            _ => {
                reach = span.end;
                drafts.push(vec![allele]);
            }
        }
    }
    loop {
        let mut built = Vec::with_capacity(drafts.len());
        for draft in &drafts {
            built.push(oracle_build_draft(draft, reference, realized)?);
        }
        // Los registros construidos que se pisan forman componentes; si todos
        // son de un solo borrador, es el punto fijo.
        let mut published = built
            .iter()
            .enumerate()
            .filter_map(|(index, (span, _))| span.map(|span| (span, index)))
            .collect::<Vec<_>>();
        published.sort();
        let mut components = Vec::<Vec<usize>>::new();
        let mut reach = 0_u64;
        for (span, index) in published {
            match components.last_mut() {
                Some(component) if span.start < reach => {
                    reach = reach.max(span.end);
                    component.push(index);
                }
                _ => {
                    reach = span.end;
                    components.push(vec![index]);
                }
            }
        }
        if components.iter().all(|component| component.len() == 1) {
            return Ok(built.into_iter().flat_map(|(_, alleles)| alleles).collect());
        }
        let mut taken = drafts.into_iter().map(Some).collect::<Vec<_>>();
        let mut next = Vec::new();
        for component in components {
            next.push(
                component
                    .iter()
                    .filter_map(|index| taken[*index].take())
                    .flatten()
                    .collect::<Vec<_>>(),
            );
        }
        // Un borrador que no publica nada no pisa a nadie y sigue igual.
        next.extend(taken.into_iter().flatten());
        drafts = next;
    }
}

/// Construye un borrador: su intervalo publicado, si publica algo, y los
/// alelos que aporta.
fn oracle_build_draft(
    draft: &[OracleAllele],
    reference: &[u8],
    realized: &OracleTracts<'_>,
) -> Result<(Option<OracleSpan>, Vec<OracleAllele>), VcfSemanticValidationError> {
    let position = draft[0].position;
    let has_tract = draft.iter().any(|allele| !allele.repeat_blocks.is_empty());
    // Un borrador con tracto se reconstruye siempre, aunque comparta `POS`,
    // para resolver dos copias con las mismas bases y anotaciones distintas.
    if !has_tract && draft.iter().all(|allele| allele.position == position) {
        let end = draft
            .iter()
            .map(|allele| oracle_vcf_span(allele).end)
            .max()
            .unwrap_or(position);
        return Ok((
            Some(OracleSpan {
                start: position,
                end,
            }),
            draft.to_vec(),
        ));
    }
    let start = draft
        .iter()
        .map(|allele| allele.position)
        .min()
        .unwrap_or(0);
    let mut end = draft
        .iter()
        .map(|allele| oracle_vcf_span(allele).end)
        .max()
        .unwrap_or(0);
    // DTC-064: el alelo que cruzó entero al alinearse un tracto realizado
    // del registro lleva el registro hasta el final de su tramo original, para
    // que el `ALT` de su copia lleve el tracto entero.
    if has_tract {
        let inside_end = end;
        for allele in draft {
            let crossed = realized
                .starting_in(start, inside_end)
                .iter()
                .any(|(span, _)| {
                    start <= span.start
                        && span.end <= inside_end
                        && allele.position < span.start
                        && span.end < allele.source.end
                });
            if crossed {
                end = end.max(allele.source.end);
            }
        }
    }
    // Los tramos que se rebanan salen de intervalos de alelos ya reconstruidos
    // contra esta misma referencia, de modo que siempre caben en ella.
    let window = |from: u64, to: u64| &reference[from as usize..to as usize];
    let mut copies = BTreeMap::<CopyIdentity, Vec<&OracleAllele>>::new();
    for allele in draft {
        copies.entry(allele.copy.clone()).or_default().push(allele);
    }
    let reference_window = window(start, end).to_vec();
    let mut carriers = Vec::<(CopyIdentity, Vec<u8>, OracleSpan)>::new();
    for (copy, mut own) in copies {
        own.sort_by_key(|allele| allele.position);
        let mut bases = Vec::new();
        let mut at = start;
        let mut source = own[0].source;
        for allele in own {
            if allele.position < at {
                return Err(oracle_error("dos alelos de una copia se pisan"));
            }
            bases.extend_from_slice(window(at, allele.position));
            bases.extend_from_slice(allele.alternate.as_bytes());
            at = oracle_vcf_span(allele).end;
            source = OracleSpan {
                start: source.start.min(allele.source.start),
                end: source.end.max(allele.source.end),
            };
        }
        bases.extend_from_slice(window(at, end));
        if bases != reference_window {
            carriers.push((copy, bases, source));
        }
    }
    if carriers.is_empty() {
        return Ok((None, Vec::new()));
    }

    let mut published_start = start;
    let mut alleles = vec![reference_window];
    alleles.extend(carriers.iter().map(|(_, bases, _)| bases.clone()));
    if !has_tract {
        // Normalización de un registro multialélico: se quita la última base
        // común; si un alelo se vacía, se antepone la base anterior de la
        // referencia a todos; y al final se quita la primera base común
        // mientras todos conserven dos.
        loop {
            if alleles.iter().any(Vec::is_empty) {
                if published_start == 0 {
                    return Err(oracle_error("registro fundido sin ancla izquierda"));
                }
                published_start -= 1;
                let base = reference[published_start as usize];
                for allele in &mut alleles {
                    allele.insert(0, base);
                }
            } else if alleles
                .iter()
                .all(|allele| allele.last() == alleles[0].last())
            {
                for allele in &mut alleles {
                    allele.pop();
                }
            } else {
                break;
            }
        }
        while alleles.iter().all(|allele| allele.len() > 1)
            && alleles.iter().all(|allele| allele[0] == alleles[0][0])
        {
            for allele in &mut alleles {
                allele.remove(0);
            }
            published_start += 1;
        }
    }
    let published_reference = alleles[0].clone();
    let published = OracleSpan {
        start: published_start,
        end: published_start + published_reference.len() as u64,
    };
    let tracts = if has_tract {
        realized
            .reaching(start, end)
            .iter()
            .filter(|(span, _)| span.start < end && start < span.end)
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    if tracts
        .iter()
        .any(|(span, _)| span.start < start || end < span.end)
    {
        return Err(oracle_error("el registro fundido parte un tracto"));
    }
    let mut distinct_tracts = Vec::<OracleSpan>::new();
    for (span, _) in &tracts {
        if !distinct_tracts.contains(span) {
            distinct_tracts.push(*span);
        }
    }
    distinct_tracts.sort();
    // La anotación de cada ALT, por posición y tracto a tracto (DTC-064): la
    // declarada si la copia realiza el tracto; si conserva el de la
    // referencia, un bloque de su motivo con todas las bases del tracto. Y la
    // derivada de las bases si otra copia con las mismas bases trae otra.
    let alternates = alleles[1..].to_vec();
    let mut annotations = Vec::with_capacity(carriers.len());
    for (copy, _, _) in &carriers {
        let path = realized.paths.get(copy).copied();
        let mut blocks = Vec::new();
        for tract in &distinct_tracts {
            let own = path.and_then(|path| realized.own.get(&(path, *tract)).copied());
            if let Some(own) = own {
                blocks.extend(own.iter().cloned());
                continue;
            }
            let motifs = tracts
                .iter()
                .filter(|(span, _)| span == tract)
                .flat_map(|(_, declared)| declared.iter().map(|block| block.motif.as_str()))
                .collect::<BTreeSet<_>>();
            let [motif] = motifs.into_iter().collect::<Vec<_>>()[..] else {
                return Err(oracle_error(
                    "una copia conserva un tracto compuesto dentro de un registro fundido",
                ));
            };
            blocks.push(TruthRepeatBlock {
                motif: motif.to_owned(),
                bases: tract.end - tract.start,
            });
        }
        annotations.push(blocks);
    }
    let mut derive = Vec::new();
    for (index, alternate) in alternates.iter().enumerate() {
        if alternates
            .iter()
            .zip(&annotations)
            .any(|(other, blocks)| other == alternate && *blocks != annotations[index])
        {
            derive.push(index);
        }
    }
    for index in derive {
        annotations[index] = oracle_derived_blocks(
            &alternates[index],
            &published_reference,
            published_start,
            &distinct_tracts,
            &tracts,
        )?;
    }
    let mut output = Vec::with_capacity(carriers.len());
    for (((copy, _, source), alternate), repeat_blocks) in
        carriers.into_iter().zip(alternates).zip(annotations)
    {
        output.push(OracleAllele {
            sequence_id: draft[0].sequence_id.clone(),
            position: published_start,
            // Las bases salen del FASTA y de las piezas de la verdad, que son
            // ADN en ASCII: la conversión no puede perder nada.
            reference: String::from_utf8_lossy(&published_reference).into_owned(),
            alternate: String::from_utf8_lossy(&alternate).into_owned(),
            source,
            copy,
            repeat_blocks,
        });
    }
    Ok((Some(published), output))
}

/// La anotación de un ALT **derivada de sus bases** (DTC-064), escrita aquí
/// por su cuenta, solo para dos copias con las mismas bases y anotaciones
/// distintas.
///
/// Se descuentan las bases del ALT que coinciden con el REF del registro por
/// delante —como mucho las que el registro tiene antes del tracto— y por
/// detrás —como mucho las que tiene después—. En lo que queda se buscan las
/// copias completas del motivo de izquierda a derecha, sin solapar, y el
/// bloque va de la primera a la última: solo lo que de verdad repite el
/// motivo. Solo hay regla con un tracto y un motivo.
fn oracle_derived_blocks(
    alternate: &[u8],
    record_reference: &[u8],
    record_start: u64,
    tracts: &[OracleSpan],
    realized: &[&(OracleSpan, &Vec<TruthRepeatBlock>)],
) -> Result<Vec<TruthRepeatBlock>, VcfSemanticValidationError> {
    let [tract] = tracts else {
        return Err(oracle_error(
            "hay que derivar la anotación de un ALT en un registro con varios tractos",
        ));
    };
    let motifs = realized
        .iter()
        .filter(|(span, _)| span == tract)
        .flat_map(|(_, blocks)| blocks.iter().map(|block| block.motif.as_str()))
        .collect::<BTreeSet<_>>();
    let [motif] = motifs.into_iter().collect::<Vec<_>>()[..] else {
        return Err(oracle_error(
            "hay que derivar la anotación de un ALT de un tracto compuesto",
        ));
    };
    let room_before = (tract.start - record_start) as usize;
    let room_after = (record_start + record_reference.len() as u64 - tract.end) as usize;
    let mut lead = 0;
    while lead < room_before && lead < alternate.len() && alternate[lead] == record_reference[lead]
    {
        lead += 1;
    }
    let mut trail = 0;
    while trail < room_after
        && lead + trail < alternate.len()
        && lead + trail < record_reference.len()
        && alternate[alternate.len() - 1 - trail]
            == record_reference[record_reference.len() - 1 - trail]
    {
        trail += 1;
    }
    let middle = &alternate[lead..alternate.len() - trail];
    let unit = motif.as_bytes();
    let mut first = None;
    let mut last_end = 0;
    let mut at = 0;
    while at + unit.len() <= middle.len() {
        if &middle[at..at + unit.len()] == unit {
            first.get_or_insert(at);
            at += unit.len();
            last_end = at;
        } else {
            at += 1;
        }
    }
    Ok(match first {
        None => Vec::new(),
        Some(first) => vec![TruthRepeatBlock {
            motif: motif.to_owned(),
            bases: (last_end - first) as u64,
        }],
    })
}

fn oracle_normalize_window(
    path: &OraclePath,
    reference_bases: &[u8],
    span: OracleSpan,
    floor: u64,
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
        // DTC-064: si el ancla es una base de un tracto que **esta** copia
        // realiza, la ventana arranca donde arranca el tracto. Anclar dentro
        // lo partiría, y la base de ancla no es la de la referencia.
        for piece in &path.pieces {
            if let OraclePiece::Replaced {
                start: left,
                end: right,
                ..
            } = piece
                && *left <= start
                && start < *right
            {
                start = *left;
            }
        }
        reference = reference_bases[usize::try_from(start)
            .map_err(|_| oracle_error("inicio fuera de rango"))?
            ..usize::try_from(end).map_err(|_| oracle_error("fin fuera de rango"))?]
            .to_vec();
        alternate = path.window_from(start, end, span.start)?;
    }
    // **El registro de repetición no se recorta**, igual que en el productor y
    // por la misma razón: la decisión S-09 manda publicar la secuencia
    // completa en `REF` y `ALT`, y DTC-064 sigue mandando en todos los demás
    // registros. El oráculo aplica la regla por su cuenta —esta es su copia,
    // no la del productor—, que es lo que hace que la comparación signifique
    // algo.
    let repeat_blocks = path.repeats_in(start, end);
    // `floor` es la posición más a la izquierda a la que puede llegar la
    // alineación: cero, o la última base de un tracto en el que acabaría.
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
            && start > floor
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
/// Sus dos rechazos los ejercita una prueba unitaria que llama a esta función
/// con un registro armado a mano (`the_expected_genotype_refuses_a_sample_without_copies_and_an_unknown_alternate`),
/// y de extremo a extremo se sabe esto de cada uno:
///
/// - **El alelo alternativo que no aparece entre los del registro es
///   inalcanzable desde `expected_records`**, que es quien construye el
///   registro: `alternates` es el conjunto de los ALT de cada copia, y cada
///   tupla de `copies` lleva el ALT de su copia o vacío, de modo que todo ALT no
///   vacío está en `alternates`.
/// - **La muestra sin copias en el contig sigue sin demostrarse inalcanzable**:
///   exigiría una verdad en la que un individuo no publicara copias en una
///   secuencia donde otro sí las publica, y ninguna verdad de la batería lo
///   hace. No se ha construido una de extremo a extremo.
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
/// De sus cuatro rechazos, tres se ejercitan: la cabecera ilegible, las bases
/// que no son ADN y la identidad de secuencia que el analizador de cabeceras
/// acepta y `SequenceId::parse` rechaza (`chr#1`, o más de 64 caracteres; la
/// prueba unitaria `a_fasta_header_the_parser_accepts_but_the_identity_type_refuses_is_rejected`).
/// El cuarto, la construcción final, es **inalcanzable**: el único error de
/// `ReferenceSequence::new` es que la identidad difiera del texto de la
/// cabecera, y la identidad acaba de salir de `SequenceId::parse` sobre ese
/// mismo texto, que la conserva tal cual.
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
    // Inalcanzable: `ReferenceSequence::new` solo rechaza que `id` no sea el texto
    // de la cabecera, y `id` acaba de salir de `SequenceId::parse` sobre ese
    // texto, que conserva tal cual. Se conserva porque esta función no se fía de
    // lo que haga `new` mañana.
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
        ExpectedRecord, OracleAllele, OraclePath, OraclePiece, OracleSpan, OracleTracts,
        VcfSemanticValidationError, expected_gt, expected_records, merge_oracle_spans,
        oracle_build_draft, oracle_normalize_window, oracle_path, oracular_sequence,
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

    fn draft_allele(copy: u8, position: u64, reference: &str, alternate: &str) -> OracleAllele {
        OracleAllele {
            sequence_id: "chr1".to_owned(),
            position,
            reference: reference.to_owned(),
            alternate: alternate.to_owned(),
            source: OracleSpan {
                start: position,
                end: position + reference.len() as u64,
            },
            copy: identity(copy),
            repeat_blocks: Vec::new(),
        }
    }

    /// **Las guardas internas del registro fundido (DTC-064).**
    ///
    /// *Quién caza antes, y por eso esta prueba es unitaria.* Ninguna verdad
    /// las alcanza de extremo a extremo, porque cada una vigila una premisa que
    /// la reconstrucción del propio oráculo ya garantiza:
    ///
    /// - **dos alelos de una copia que se pisan**: la normalización por copia
    ///   funde antes los alelos solapados de una misma copia;
    /// - **un registro que habría que anclar antes de la primera base**: la
    ///   normalización por copia ya rechaza el *indel* sin ancla izquierda;
    /// - **un tracto que el registro fundido parte**: un tracto realizado es el
    ///   alelo de su copia, y ese alelo pisa el registro y entra en él entero;
    /// - **una anotación que habría que derivar de las bases con varios tractos
    ///   en el registro**: hacen falta dos microsatélites en un mismo registro
    ///   fundido y, además, dos copias con las mismas bases y anotaciones
    ///   distintas; ningún escenario de la batería lo hace, y se llama aquí
    ///   directamente a la derivación. El productor tiene su propio rechazo
    ///   con diagnóstico para ese caso.
    ///
    /// Se llama directamente al constructor con la única entrada que cada
    /// guarda puede rechazar, y se exige su motivo. El control fija que dos
    /// alelos de copias distintas sin nada de eso se funden sin error.
    fn no_tracts() -> OracleTracts<'static> {
        OracleTracts::indexed(Vec::new(), BTreeMap::new(), BTreeMap::new())
    }

    #[test]
    fn the_fused_record_builder_refuses_what_its_own_reconstruction_never_produces() {
        let control = oracle_build_draft(
            &[draft_allele(0, 1, "CG", "C"), draft_allele(1, 2, "G", "T")],
            b"ACGTACGT",
            &no_tracts(),
        )
        .expect("control: dos copias distintas se funden");
        assert_eq!(control.0, Some(OracleSpan { start: 1, end: 3 }));
        assert_eq!(control.1.len(), 2);

        let overlapping = oracle_build_draft(
            &[draft_allele(0, 1, "CG", "C"), draft_allele(0, 2, "G", "T")],
            b"ACGTACGT",
            &no_tracts(),
        )
        .expect_err("una copia no puede llevar dos alelos que se pisan");
        assert_eq!(detail(&overlapping), "dos alelos de una copia se pisan");

        // Las dos copias borran un `CA` del `ACAC` inicial: el registro se
        // vacía al recortarlo y habría que anclarlo antes de la base 0.
        let unanchored = oracle_build_draft(
            &[
                draft_allele(0, 0, "ACA", "A"),
                draft_allele(1, 1, "CAC", "C"),
            ],
            b"ACACGTGT",
            &no_tracts(),
        )
        .expect_err("no hay base anterior a la primera");
        assert_eq!(detail(&unanchored), "registro fundido sin ancla izquierda");

        let block = vec![TruthRepeatBlock {
            motif: "CG".to_owned(),
            bases: 2,
        }];
        let mut with_tract = draft_allele(0, 1, "CG", "C");
        with_tract.repeat_blocks = block.clone();
        let cut = oracle_build_draft(
            &[with_tract.clone(), draft_allele(1, 2, "G", "T")],
            b"ACGTACGT",
            &OracleTracts::indexed(
                vec![(OracleSpan { start: 2, end: 6 }, &block)],
                BTreeMap::new(),
                BTreeMap::new(),
            ),
        )
        .expect_err("el registro no puede partir un tracto");
        assert_eq!(detail(&cut), "el registro fundido parte un tracto");

        let tracts = [
            (OracleSpan { start: 1, end: 3 }, &block),
            (OracleSpan { start: 4, end: 6 }, &block),
        ];
        let several = super::oracle_derived_blocks(
            b"ATTTACGT",
            b"ACGTACGT",
            0,
            &[tracts[0].0, tracts[1].0],
            &tracts.iter().collect::<Vec<_>>(),
        )
        .expect_err("con dos tractos no se sabe qué bases son de cuál");
        assert_eq!(
            detail(&several),
            "hay que derivar la anotación de un ALT en un registro con varios tractos"
        );
    }

    /// **La comparación final de cada copia rechaza la que no coincide.**
    ///
    /// *Por qué es unitaria, dicho sin adornos.* De extremo a extremo no se
    /// alcanza con un oráculo correcto: un VCF que supera todas las demás
    /// comprobaciones coincide con los registros que este oráculo deriva, y si
    /// la derivación es correcta esos registros describen las copias. Esta
    /// guarda existe para el día en que productor y oráculo compartan un error
    /// de derivación, que ya ocurrió dos veces; lo que demuestra que caza ese
    /// caso es el mutante que reintroduce uno de ellos en los dos. Aquí se fija
    /// que una huella distinta, o una copia ausente, se rechaza, y que las
    /// mismas huellas se aceptan.
    #[test]
    fn the_rebuilt_copies_must_be_the_truth_copies() {
        let key = ("chr1".to_owned(), "IND_01".to_owned(), 0_usize);
        let expected = super::CopyDigests::from([(key.clone(), ([7_u8; 32], 8))]);
        super::compare_rebuilt_copies(&expected, expected.clone())
            .expect("control: las mismas huellas se aceptan");
        let wrong = super::CopyDigests::from([(key, ([8_u8; 32], 8))]);
        let error = super::compare_rebuilt_copies(&expected, wrong)
            .expect_err("una copia con otra huella se rechaza");
        assert_eq!(
            detail(&error),
            "la copia IND_01 fase 0 de chr1 que describe el VCF no es la de la verdad"
        );
        let error = super::compare_rebuilt_copies(&expected, super::CopyDigests::new())
            .expect_err("una copia que el VCF no describe se rechaza");
        assert!(detail(&error).contains("no es la de la verdad"));
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
        let error = oracle_normalize_window(&path, reference, OracleSpan { start: 0, end: 0 }, 0)
            .expect_err("una inserción en la posición cero no tiene ancla");
        assert_eq!(detail(&error), "indel sin ancla izquierda");

        // Control positivo: la misma inserción una base más allá sí normaliza,
        // y lo hace anclando en la base anterior.
        let anchored = inserted_at(1);
        let allele =
            oracle_normalize_window(&anchored, reference, OracleSpan { start: 1, end: 1 }, 0)
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
            "../tests/fixtures/truth-v4/valid.jsonl"
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
        // Desde el 2 de octubre de 2026 la anotación de ese ALT se deriva de
        // sus bases (DTC-064), y con dos motivos no hay regla: se rechaza.
        let error = expected_records(&reference, &contradictory)
            .expect_err("el mismo ALT con dos descomposiciones de motivos distintos");
        assert_eq!(
            detail(&error),
            "hay que derivar la anotación de un ALT de un tracto compuesto"
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

    /// **Un caso sin ALT exige que su copia sea la referencia alrededor del
    /// locus**, y el oráculo no puede rechazar uno que lo es.
    ///
    /// *Quién caza antes.* Con una verdad validada nadie llega aquí con un
    /// caso cuya copia cambie en el locus sin ALT: el juez de la verdad fija
    /// el locus de cada caso donde su evento cambió la copia (véase
    /// [`super::every_case_has_an_alternate_or_an_unchanged_copy`]). Por eso
    /// se llama a `expected_records` con una verdad sin validar.
    ///
    /// *Qué hace que pueda fallar:*
    ///
    /// - se admite un caso cuya copia es la referencia en la ventana que la
    ///   copia vecina obliga a componer;
    /// - se admite un caso en un sitio que **ninguna** copia cambia, de modo
    ///   que no hay ventana: es la forma de una inserción que una deleción
    ///   posterior borra entera. Hasta el 3 de octubre de 2026 el oráculo la
    ///   rechazaba en falso con el motivo de esta guarda;
    /// - se rechaza, con su motivo exacto, un caso cuyo locus no está en
    ///   ningún alelo de su copia y que toca un cambio de ella.
    #[test]
    fn a_case_without_an_alternate_needs_an_unchanged_copy() {
        use crate::truth_schema::{
            AlleleCarrier, GenotypeCaseKind, GenotypeRecord, TruthAllele, TruthLocus,
        };
        let reference = tiny_reference(b"AAAAAAAA");
        // La copia 0 lleva el cambio de base en 2; la copia 1 es la
        // referencia, pero la ventana [2, 3) se compone para las dos.
        let with_case = |case_id: &str, copy_index: u8, start: u64| {
            let mut truth = fixture_truth();
            truth.records.push(TruthRecord::Genotype(GenotypeRecord {
                case_id: case_id.to_owned(),
                reference_id: "REF_main".to_owned(),
                individual_id: "IND_01".to_owned(),
                case_kind: GenotypeCaseKind::Snp,
                locus: TruthLocus {
                    sequence_id: "chr1".to_owned(),
                    start,
                    end: start + 1,
                },
                alleles: vec![
                    TruthAllele::Reference,
                    TruthAllele::Snp {
                        base: "G".to_owned(),
                    },
                ],
                allele_dosages: vec![1, 1],
                carriers: vec![
                    AlleleCarrier {
                        allele_index: usize::from(copy_index == 0),
                        copy: identity(0),
                    },
                    AlleleCarrier {
                        allele_index: usize::from(copy_index == 1),
                        copy: identity(1),
                    },
                ],
            }));
            truth
        };

        for (copy_index, start, why) in [
            (
                1,
                2,
                "la copia 1 queda igual a la referencia en la ventana [2, 3)",
            ),
            (
                0,
                5,
                "ninguna copia cambia en 5: sin ventana, la copia 0 es la referencia",
            ),
            (
                1,
                5,
                "ninguna copia cambia en 5: sin ventana, la copia 1 es la referencia",
            ),
        ] {
            let records = expected_records(&reference, &with_case("CASE_02", copy_index, start))
                .unwrap_or_else(|error| panic!("{why}: {error:?}"));
            assert!(
                records
                    .values()
                    .all(|record| !record.gps_cases.contains(&"CASE_02".to_owned())),
                "{why}: un caso sin ALT no entra en ningún GPS_CASE"
            );
        }

        let error = expected_records(&reference, &with_case("CASE_03", 0, 3))
            .expect_err("la copia 0 cambia en 2, que toca el locus [3, 4) por su extremo");
        assert_eq!(
            detail(&error),
            "un caso sin ALT cuya copia no queda igual a la referencia"
        );
    }

    /// **Las dos pruebas de que una copia es la referencia, una por una.**
    ///
    /// *Por qué se llama a la guarda directamente.* Ninguna verdad combina a
    /// la vez un caso sin ALT y una copia que cambia junto a su locus (véase
    /// la prueba anterior), y aquí hace falta fijar por separado cada
    /// condición: la ventana tiene que ser **de la misma copia** y
    /// **contener** el locus entero, no solo tocarlo; y un cambio de la copia
    /// que toca el locus por cualquiera de sus extremos ya impide la segunda
    /// prueba.
    #[test]
    fn each_proof_of_an_unchanged_copy_is_checked_on_its_own() {
        use super::every_case_has_an_alternate_or_an_unchanged_copy;
        use crate::truth_schema::{
            AlleleCarrier, GenotypeCaseKind, GenotypeRecord, TruthAllele, TruthLocus,
        };
        let mut truth = fixture_truth();
        truth
            .records
            .retain(|record| matches!(record, TruthRecord::HaplotypeCopy(_)));
        // Un caso de la copia 0 con locus [3, 5).
        truth.records.push(TruthRecord::Genotype(GenotypeRecord {
            case_id: "CASE_X".to_owned(),
            reference_id: "REF_main".to_owned(),
            individual_id: "IND_01".to_owned(),
            case_kind: GenotypeCaseKind::Snp,
            locus: TruthLocus {
                sequence_id: "chr1".to_owned(),
                start: 3,
                end: 5,
            },
            alleles: vec![
                TruthAllele::Reference,
                TruthAllele::Snp {
                    base: "G".to_owned(),
                },
            ],
            allele_dosages: vec![1, 1],
            carriers: vec![
                AlleleCarrier {
                    allele_index: 1,
                    copy: identity(0),
                },
                AlleleCarrier {
                    allele_index: 0,
                    copy: identity(1),
                },
            ],
        }));
        let span = |start, end| OracleSpan { start, end };
        let check = |windows: Vec<(u8, OracleSpan)>, changes: Vec<(u8, Vec<OracleSpan>)>| {
            let mut unchanged = BTreeMap::new();
            for (copy, window) in windows {
                unchanged
                    .entry(identity(copy))
                    .or_insert_with(Vec::new)
                    .push(window);
            }
            let copy_changes = changes
                .into_iter()
                .map(|(copy, spans)| (identity(copy), spans))
                .collect::<BTreeMap<_, _>>();
            every_case_has_an_alternate_or_an_unchanged_copy(&truth, &[], &unchanged, &copy_changes)
        };
        // La copia 0 cambia en 4, dentro del locus, para que solo una ventana
        // pueda excusar el caso.
        let touching = vec![(0, vec![span(4, 5)])];
        check(vec![(0, span(3, 5))], touching.clone())
            .expect("control: una ventana de la copia 0 que contiene el locus");
        check(vec![], vec![(0, vec![span(6, 7), span(1, 2)])])
            .expect("control: ningún cambio de la copia 0 toca [3, 5]");
        for (windows, changes, why) in [
            (
                vec![],
                touching.clone(),
                "sin ventana y con un cambio dentro",
            ),
            (
                vec![(1, span(3, 5))],
                touching.clone(),
                "con la ventana de otra copia",
            ),
            (
                vec![(0, span(3, 4))],
                touching.clone(),
                "con una ventana que no cubre el final",
            ),
            (
                vec![(0, span(4, 6))],
                touching.clone(),
                "con una ventana que no cubre el inicio",
            ),
            (
                vec![(0, span(1, 3))],
                touching.clone(),
                "con una ventana que solo toca el inicio",
            ),
            (
                vec![],
                vec![(0, vec![span(5, 6)])],
                "con un cambio que toca el extremo final",
            ),
            (
                vec![],
                vec![(0, vec![span(2, 3)])],
                "con un cambio que toca el extremo inicial",
            ),
            (
                vec![],
                vec![(0, vec![span(3, 3)])],
                "con una inserción en el extremo inicial",
            ),
            (vec![], vec![], "con una copia cuyos cambios no se conocen"),
        ] {
            let error = check(windows, changes).expect_err(why);
            assert_eq!(
                detail(&error),
                "un caso sin ALT cuya copia no queda igual a la referencia",
                "{why}"
            );
        }
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

    /// **El genotipo esperado de una muestra: qué exige y qué rechaza.**
    ///
    /// *Quién caza antes, y por qué esta prueba es unitaria.* `expected_gt` solo
    /// la llama la comparación de cada columna de muestra, con un
    /// `ExpectedRecord` que construye `expected_records`, y esa construcción
    /// hace que sus dos rechazos no se alcancen desde un VCF:
    ///
    /// - **El ALT ausente del registro.** `alternates` es el conjunto de los ALT
    ///   de `per_copy`, y cada tupla de `copies` lleva exactamente el ALT que
    ///   `per_copy` da a su copia, o vacío. Todo ALT no vacío de `copies` está
    ///   por tanto en `alternates`, y `position` nunca da `None`. Es
    ///   **inalcanzable desde `expected_records`**: solo un registro armado a
    ///   mano, como el de aquí, lo dispara.
    /// - **La muestra sin copias en el contig.** Un individuo que no publicara
    ///   ninguna copia en una secuencia donde otro sí las publica no lo
    ///   construye ninguna verdad de la batería. Queda sin demostrarse que sea
    ///   inalcanzable de extremo a extremo; aquí se ejercita el rechazo mismo.
    ///
    /// *Qué hace que pueda fallar:* el caso íntegro fija el genotipo que se
    /// acepta —una guarda que rechazara todo pasaría el resto—, y cada rechazo
    /// exige su motivo exacto.
    #[test]
    fn the_expected_genotype_refuses_a_sample_without_copies_and_an_unknown_alternate() {
        let record = |alternate_of_second_copy: &str| ExpectedRecord {
            reference: "A".to_owned(),
            alternates: vec!["T".to_owned()],
            gps_cases: Vec::new(),
            copies: BTreeMap::from([(
                "IND_01".to_owned(),
                vec![
                    (0, identity(0), String::new()),
                    (1, identity(1), alternate_of_second_copy.to_owned()),
                ],
            )]),
            alternate_repeats: vec![Vec::new()],
        };
        assert_eq!(
            expected_gt(&record("T"), "IND_01").expect("un genotipo de dos copias"),
            "0|1"
        );
        let error =
            expected_gt(&record("T"), "IND_02").expect_err("una muestra sin copias en el contig");
        assert_eq!(
            detail(&error),
            "muestra VCF no tiene copias truth para el contig"
        );
        let error = expected_gt(&record("G"), "IND_01")
            .expect_err("un ALT que el registro esperado no tiene");
        assert_eq!(
            detail(&error),
            "ALT truth no aparece en el registro esperado"
        );
    }

    /// **La identidad de una secuencia FASTA que el analizador de cabeceras
    /// acepta y el tipo de identidad rechaza.**
    ///
    /// `FastaHeader::parse` toma como identidad cualquier texto sin espacios;
    /// `SequenceId::parse` exige de 1 a 64 caracteres ASCII, alfanumérico el
    /// primero y `.`, `_` o `-` después. Entre los dos cabe, por ejemplo,
    /// `chr#1`, y esa diferencia es lo que esta guarda vigila. La medición
    /// anterior no había encontrado una cabecera así: basta una con un carácter
    /// fuera del alfabeto de identidades.
    ///
    /// El rechazo lleva el diagnóstico de identidad (`GPS-DOMAIN-001`), que es
    /// lo que la distingue de los otros tres rechazos de esta función: cabecera
    /// ilegible y bases que no son ADN.
    #[test]
    fn a_fasta_header_the_parser_accepts_but_the_identity_type_refuses_is_rejected() {
        oracular_sequence("chr1 descripción", b"ACGT".to_vec(), "REF_main")
            .expect("una cabecera con identidad válida se acepta");
        for header in ["chr#1", "-chr", &"c".repeat(65)] {
            let error = oracular_sequence(header, b"ACGT".to_vec(), "REF_main")
                .expect_err("una identidad fuera del alfabeto o del largo admitido");
            let VcfSemanticValidationError::Fasta { detail } = error else {
                panic!("se esperaba un rechazo del FASTA");
            };
            assert!(
                detail.contains("GPS-DOMAIN-001"),
                "«{header}» se rechazó por otra razón: {detail}"
            );
        }
    }

    /// **F6.4.a.8b. El índice de tractos del juez no pierde el que empieza
    /// antes de la ventana y acaba dentro.**
    ///
    /// Un tracto largo `[100, 400)` y dos cortos `[450, 460)` de dos copias.
    /// La ventana `[300, 455)` la pisa el largo, que empieza **antes**; un
    /// índice que buscara por inicio lo perdería. Igual para el límite de la
    /// alineación: la base 300 está dentro del largo, que acaba después de la
    /// 301. Se exige el conjunto exacto tras aplicar la condición de siempre,
    /// y que los del mismo intervalo conserven el orden de las copias.
    #[test]
    fn the_oracle_tract_index_keeps_a_tract_that_starts_before_the_window() {
        let long = vec![TruthRepeatBlock {
            motif: "CA".to_owned(),
            bases: 300,
        }];
        let short_first = vec![TruthRepeatBlock {
            motif: "GT".to_owned(),
            bases: 10,
        }];
        let short_second = vec![TruthRepeatBlock {
            motif: "GT".to_owned(),
            bases: 12,
        }];
        let span = |start: u64, end: u64| OracleSpan { start, end };
        let tracts = OracleTracts::indexed(
            vec![
                (span(450, 460), &short_first),
                (span(100, 400), &long),
                (span(450, 460), &short_second),
            ],
            BTreeMap::new(),
            BTreeMap::new(),
        );
        let found = |slice: &[(OracleSpan, &Vec<TruthRepeatBlock>)],
                     keep: &dyn Fn(OracleSpan) -> bool| {
            slice
                .iter()
                .filter(|(span, _)| keep(*span))
                .map(|(span, blocks)| (span.start, span.end, blocks[0].bases))
                .collect::<Vec<_>>()
        };
        assert_eq!(
            found(tracts.reaching(300, 455), &|span| span.start < 455
                && 300 < span.end),
            vec![(100, 400, 300), (450, 460, 10), (450, 460, 12)]
        );
        assert_eq!(
            found(tracts.ending_after(300, 301), &|span| span.start <= 300
                && 301 < span.end),
            vec![(100, 400, 300)]
        );
        assert_eq!(
            found(tracts.starting_in(440, 460), &|_| true),
            vec![(450, 460, 10), (450, 460, 12)]
        );
    }

    // -----------------------------------------------------------------------
    // El índice de tractos del juez contra el recorrido completo
    //
    // El juez ya no recorre todos los tractos: repite el algoritmo del
    // productor (ordenar, fin máximo acumulado, búsqueda binaria), y con ello
    // pierde independencia. Estas pruebas se la devuelven: el oráculo es el
    // **recorrido completo**, un filtro lineal sin índice.
    // -----------------------------------------------------------------------

    /// Generador xorshift de semilla fija: las pruebas son deterministas.
    struct XorShift(u64);

    impl XorShift {
        fn below(&mut self, bound: u64) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0 % bound
        }
    }

    #[test]
    fn oracle_tracts_queries_match_a_full_scan_on_random_inputs() {
        let block = vec![TruthRepeatBlock {
            motif: "CA".into(),
            bases: 2,
        }];
        let mut rng = XorShift(0x1234_5678_9ABC_DEF1);
        for _round in 0..2500 {
            let spans = (0..rng.below(12))
                .map(|_| {
                    let start = rng.below(80);
                    // Alguno largo, para que el fin acumulado difiera del propio.
                    let long = rng.below(4) == 0;
                    let length = 1 + rng.below(if long { 60 } else { 8 });
                    (
                        OracleSpan {
                            start,
                            end: start + length,
                        },
                        &block,
                    )
                })
                .collect::<Vec<_>>();
            let all = spans.iter().map(|(span, _)| *span).collect::<Vec<_>>();
            let tracts = OracleTracts::indexed(spans, BTreeMap::new(), BTreeMap::new());
            for _query in 0..20 {
                // Un tercio de las consultas arranca justo en el inicio de un
                // tracto: el borde que un `<` en vez de `<=` perdería.
                let from = if !all.is_empty() && rng.below(3) == 0 {
                    all[rng.below(all.len() as u64) as usize].start
                } else {
                    rng.below(100)
                };
                let to = from + rng.below(30);

                let mut got = tracts
                    .reaching(from, to)
                    .iter()
                    .map(|(span, _)| *span)
                    .filter(|span| span.start < to && from < span.end)
                    .collect::<Vec<_>>();
                let mut want = all
                    .iter()
                    .copied()
                    .filter(|span| span.start < to && from < span.end)
                    .collect::<Vec<_>>();
                got.sort_unstable();
                want.sort_unstable();
                assert_eq!(got, want, "reaching({from}, {to})");

                let mut got = tracts
                    .ending_after(from, from + 1)
                    .iter()
                    .map(|(span, _)| *span)
                    .filter(|span| span.start <= from && from + 1 < span.end)
                    .collect::<Vec<_>>();
                let mut want = all
                    .iter()
                    .copied()
                    .filter(|span| span.start <= from && from + 1 < span.end)
                    .collect::<Vec<_>>();
                got.sort_unstable();
                want.sort_unstable();
                assert_eq!(got, want, "ending_after({from}, {})", from + 1);

                let mut got = tracts
                    .starting_in(from, to)
                    .iter()
                    .map(|(span, _)| *span)
                    .collect::<Vec<_>>();
                let mut want = all
                    .iter()
                    .copied()
                    .filter(|span| from <= span.start && span.start <= to)
                    .collect::<Vec<_>>();
                got.sort_unstable();
                want.sort_unstable();
                assert_eq!(got, want, "starting_in({from}, {to})");
            }
        }
    }

    #[test]
    fn ending_after_keeps_a_tract_that_starts_exactly_at_the_anchor_base() {
        // Un tracto [10, 20) que empieza **en** `position`, la base de
        // anclaje: `ending_after` debe devolverlo, y con `<` en lugar de `<=`
        // en su búsqueda binaria desaparecería. Los tractos [0, 5) y [2, 9)
        // van antes y su fin acumulado no pasa de `after`.
        let block = vec![TruthRepeatBlock {
            motif: "CA".into(),
            bases: 2,
        }];
        let spans = [(0, 5), (2, 9), (10, 20)]
            .into_iter()
            .map(|(start, end)| (OracleSpan { start, end }, &block))
            .collect::<Vec<_>>();
        let tracts = OracleTracts::indexed(spans, BTreeMap::new(), BTreeMap::new());
        let found = tracts
            .ending_after(10, 11)
            .iter()
            .map(|(span, _)| *span)
            .filter(|span| span.start <= 10 && 11 < span.end)
            .collect::<Vec<_>>();
        assert_eq!(found, vec![OracleSpan { start: 10, end: 20 }]);
        // Un paso a la izquierda, el tracto ya no cuenta como contenedor.
        let found = tracts
            .ending_after(9, 10)
            .iter()
            .filter(|(span, _)| span.start <= 9 && 10 < span.end)
            .count();
        assert_eq!(found, 0);
    }
}
