//! **Juez independiente del catálogo de microsatélites publicado**
//! (REQ-VAL-64, DTC-104; paso F6.4.a.15 del plan de detección).
//!
//! # Qué juzga
//!
//! El catálogo aparte que publica una ejecución con la sección de
//! microsatélites (REQ-MAN-98): la tabla `repeat_catalog.tsv`, su metadato
//! `repeat_catalog.metadata.json`, la tabla de alelos de los detectados
//! `repeat_alleles.tsv` y los enlaces con huella de la configuración efectiva
//! (`repeat_detection.catalog` y `repeat_detection.alleles_table`). Lo relee **contra el FASTA final
//! publicado** y exige:
//!
//! - **identidad**: que las huellas de la tabla, del metadato y del FASTA final
//!   sean las de los ficheros y las que dicen el manifiesto, la efectiva y el
//!   metadato; que los criterios del metadato sean los de la efectiva;
//! - **de cada tramo detectado** —los loci detectados de la tabla, los tramos
//!   que un declarado sustituye y los retirados por una variante—: que está
//!   dentro de su secuencia con una base A/C/G/T a cada lado; que su motivo
//!   tiene de 1 a 6 bases, es primitivo y es la lectura de las bases del
//!   inicio del tramo; que sus copias completas alcanzan el suelo de su
//!   longitud de motivo; que su pureza y su tirada más larga alcanzan los
//!   mínimos; que no rebasa el tope de unidades; que cabe entero en una región
//!   si las hay; y que su puntuación es la de sus propias comparaciones;
//! - que el identificador de cada detectado y de cada retirado sale de sus
//!   coordenadas por la fórmula `SSR-<secuencia en base uno>-<inicio>-<fin>-<motivo>`;
//! - que los loci de la tabla van por secuencia y por inicio y **no se solapan**,
//!   y que los detectados, los sustituidos y los retirados tampoco entre sí;
//! - que **cada declarado coincide** con un tramo sustituido de su secuencia
//!   cuyo motivo es una rotación de uno de los suyos (REQ-REF-40), salvo los
//!   de motivos de 7 bases o más, y que cada sustituido toca a un declarado;
//! - que **los recuentos cuadran**: por longitud de motivo, intervalos del
//!   detector = descartes de las ocho etapas + aceptados, y aceptados =
//!   sustituidos + retirados + sin copias haplotípicas + publicados, con cada
//!   sumando igual a lo que listan la tabla y el metadato;
//! - que **cada retirado contiene de verdad su variante**: la posición cae
//!   dentro del tramo y la variante existe ahí, en la matriz de SNV o en los
//!   cambios directos publicados;
//! - que las filas declaradas de la tabla son los loci de `repeat_loci` de la
//!   efectiva, y las detectadas, los loci de la tabla de alelos
//!   (`repeat_alleles.tsv`), en el mismo orden (entrega 4 bis);
//! - y de la tabla de alelos, enlazada con huella desde la efectiva
//!   (`repeat_detection.alleles_table`) y listada en el manifiesto: su forma
//!   (cabecera, trece columnas, enteros canónicos, proporciones en `[0, 1]`),
//!   índices de alelo de 0 en 0, origen `generated`, secuencias de A, C, G, T
//!   y N distintas dentro de cada locus, pasos aplicados más no aplicados
//!   iguales a los sorteados en cada fila, copias realizadas que suman las
//!   copias haplotípicas publicadas y cada frecuencia realizada igual, bit a
//!   bit, a sus copias entre las del locus.
//!
//! El abanico de cada detectado no lo rehace este juez sino el de REQ-VAL-60
//! (`validate_published_detected_fans`, en `src/results.rs`), que lee la tabla
//! de alelos con el mismo lector estricto de aquí.
//!
//! # Por qué no es un espejo
//!
//! No llama al detector (`src/repeat_detection.rs`), ni a la medida de pureza
//! del productor (`observe_tract`), ni a la conversión de tramo en locus, ni a
//! nada que los use: la descomposición en copias, la pureza, la tirada, la
//! puntuación, la primitividad y la rotación se rehacen aquí con aritmética
//! propia. Si llamara al productor compararía el catálogo consigo mismo y no
//! podría fallar (REQ-VAL-64). Lo vigila
//! `tests/repeat_catalog_validation.rs`, que también comprueba con una copia
//! del árbol que un productor con la pureza estropeada muere aquí.
//!
//! # Lo que no puede comprobar, dicho sin adornos
//!
//! **La exhaustividad**: que no falte ningún tramo que cumpla los criterios.
//! Comprobarlo exige un detector, y este juez no lo tiene a propósito. Un
//! tracto puro de muchas copias absorbido en un segmento impuro que falle la
//! pureza desaparece sin aviso, y este juez no lo ve (riesgo R-7 del plan).
//! Tampoco puede rehacer los descartes por etapa: solo exige que cuadren.
//! Las dos cosas las ata la reconstrucción del juez de los recursos
//! publicados (`src/results.rs`), que sí llama al productor, y las cubren las
//! pruebas del detector con verdad plantada y con un buscador por fuerza bruta.
//!
//! Cada rechazo pasa por [`catalog_reject`], anotada para la puerta de
//! cobertura de los jueces.

use std::collections::BTreeMap;
use std::path::Path;

use serde_json::{Map, Value, json};
use sha2::{Digest, Sha256};

use crate::diagnostics::Diagnostic;

const METADATA_ID: &str = "repeat-catalog-metadata";
const METADATA_KIND: &str = "repeat_catalog_metadata";
const METADATA_PATH: &str = "repeat_catalog.metadata.json";
const TABLE_ID: &str = "repeat-catalog-tsv";
const TABLE_KIND: &str = "repeat_catalog_tsv";
const TABLE_PATH: &str = "repeat_catalog.tsv";
const TABLE_HEADER: &str = "locus_id\tsequence_id\tstart\tend\tmotif\tscore\torigin";
const ALLELES_ID: &str = "repeat-alleles-tsv";
const ALLELES_KIND: &str = "repeat_alleles_tsv";
/// Ruta publicada de la tabla de alelos de los detectados (DTC-104).
pub const ALLELES_PATH: &str = "repeat_alleles.tsv";
const ALLELES_HEADER: &str = "locus_id\tallele_index\tsequence\torigin\tdeclared_frequency\t\
effective_frequency\trealized_copies\trealized_frequency\tdrawn_events\tapplied_events\t\
unapplied_events\ttruth_tract_bases_bound\tvcf_tract_bases_bound";
/// Las ocho etapas de descarte del barrido, en su orden (DTC-103).
const STAGES: [&str; 8] = [
    "non_canonical_motif",
    "min_copies",
    "min_purity",
    "min_longest_run",
    "max_units",
    "edge",
    "outside_regions",
    "overlap",
];
/// Longitud máxima de un motivo detectado (REQ-REF-38).
const MAX_MOTIF_LENGTH: usize = 6;
/// La puntuación fija del detector, que DTC-103 declara constante del
/// algoritmo: +2 por base igual a la de un periodo antes, −7 si no.
const MATCH_REWARD: i64 = 2;
const MISMATCH_PENALTY: i64 = 7;

/// Un artefacto de la lista del manifiesto general.
#[derive(Clone, Copy, Debug)]
pub struct PublishedArtifact<'a> {
    pub id: &'a str,
    pub kind: &'a str,
    pub path: &'a str,
    pub size_bytes: u64,
    pub digest: &'a str,
}

/// Una secuencia del FASTA final publicado, con sus bases.
#[derive(Clone, Copy, Debug)]
pub struct PublishedSequence<'a> {
    pub id: &'a str,
    pub bases: &'a [u8],
}

/// Lo que el juez necesita de un paquete publicado.
#[derive(Clone, Copy, Debug)]
pub struct PublishedCatalogInputs<'a> {
    /// Raíz del paquete.
    pub root: &'a Path,
    /// La configuración efectiva publicada.
    pub effective: &'a Map<String, Value>,
    /// La lista de artefactos del manifiesto general.
    pub artifacts: &'a [PublishedArtifact<'a>],
    /// Identificador de la referencia final.
    pub reference_id: &'a str,
    /// Las secuencias del FASTA final publicado, en su orden.
    pub sequences: &'a [PublishedSequence<'a>],
    /// Las filas de la matriz de SNV publicada: `(secuencia, posición en base
    /// uno, identificador)`.
    pub matrix_variants: &'a [(&'a str, u64, &'a str)],
}

/// El único constructor de rechazo del juez: anota la línea que rechaza.
///
/// Usa el código de un documento producido que no cumple su contrato
/// (`GPS-FMT-006`), cuyo texto fijo del catálogo de diagnósticos habla de un
/// documento JSON también cuando el fallo es de la tabla; `path` y `schema`
/// nombran el fichero de verdad: el metadato y su esquema, o la tabla y la
/// decisión que fija sus columnas (DTC-104), que no tiene esquema JSON.
#[track_caller]
fn catalog_reject(file: &str, detail: impl Into<String>) -> Diagnostic {
    crate::guard_trace::record();
    Diagnostic::format(
        "GPS-FMT-006",
        "El catálogo de microsatélites publicado no satisface su contrato.",
    )
    .with("path", file)
    .with(
        "schema",
        match file {
            TABLE_PATH => "repeat_catalog.tsv (columnas de DTC-104)",
            ALLELES_PATH => "repeat_alleles.tsv (columnas de DTC-104)",
            _ => "repeat-catalog-metadata-v1.schema.json",
        },
    )
    .with("detail", detail.into())
}

/// Juzga el catálogo de microsatélites publicado. Sin la sección de
/// microsatélites en la efectiva solo exige que no se haya publicado ningún
/// catálogo.
pub fn validate_published_repeat_catalog(
    inputs: &PublishedCatalogInputs<'_>,
) -> Result<(), Diagnostic> {
    let artifact = |id: &str| inputs.artifacts.iter().find(|artifact| artifact.id == id);
    let Some(detection) = inputs.effective.get("repeat_detection") else {
        if artifact(METADATA_ID).is_some()
            || artifact(TABLE_ID).is_some()
            || artifact(ALLELES_ID).is_some()
        {
            return Err(catalog_reject(
                METADATA_PATH,
                "hay un catálogo publicado y la configuración efectiva no declara la sección de microsatélites",
            ));
        }
        return Ok(());
    };

    // --- Identidad: enlace, artefactos y huellas -------------------------
    let link = object(detection.get("catalog"), "repeat_detection.catalog")?;
    if text(link.get("metadata_path"), "catalog.metadata_path")? != METADATA_PATH {
        return Err(catalog_reject(
            METADATA_PATH,
            "repeat_detection.catalog.metadata_path no es la ruta del metadato",
        ));
    }
    let linked_sha = text(link.get("metadata_sha256"), "catalog.metadata_sha256")?;
    let listed = |id: &str, kind: &str, path: &str| {
        artifact(id)
            .filter(|artifact| artifact.kind == kind && artifact.path == path)
            .ok_or_else(|| {
                catalog_reject(
                    path,
                    format!("el manifiesto no lista {id} con clase {kind} y ruta {path}"),
                )
            })
    };
    let metadata_artifact = listed(METADATA_ID, METADATA_KIND, METADATA_PATH)?;
    let table_artifact = listed(TABLE_ID, TABLE_KIND, TABLE_PATH)?;
    let read = |path: &str| {
        std::fs::read(inputs.root.join(path))
            .map_err(|error| catalog_reject(path, format!("no se puede leer: {error}")))
    };
    let metadata_bytes = read(METADATA_PATH)?;
    let table_bytes = read(TABLE_PATH)?;
    let metadata_sha = sha256_hex(&metadata_bytes);
    let table_sha = sha256_hex(&table_bytes);
    if metadata_sha != linked_sha || metadata_sha != metadata_artifact.digest {
        return Err(catalog_reject(
            METADATA_PATH,
            "la huella del metadato no es la del enlace de la efectiva o la del manifiesto",
        ));
    }
    if table_sha != table_artifact.digest {
        return Err(catalog_reject(
            TABLE_PATH,
            "la huella de la tabla no es la del manifiesto",
        ));
    }
    if metadata_bytes.len() as u64 != metadata_artifact.size_bytes {
        return Err(catalog_reject(
            METADATA_PATH,
            "el tamaño del metadato no es el del manifiesto",
        ));
    }
    if table_bytes.len() as u64 != table_artifact.size_bytes {
        return Err(catalog_reject(
            TABLE_PATH,
            "el tamaño de la tabla no es el del manifiesto",
        ));
    }
    // La tabla de alelos de los detectados (entrega 4 bis): enlace con huella
    // desde la efectiva y artefacto del manifiesto, como el par del catálogo.
    let alleles_link = object(
        detection.get("alleles_table"),
        "repeat_detection.alleles_table",
    )?;
    if text(alleles_link.get("path"), "alleles_table.path")? != ALLELES_PATH {
        return Err(catalog_reject(
            ALLELES_PATH,
            "repeat_detection.alleles_table.path no es la ruta de la tabla de alelos",
        ));
    }
    let alleles_linked_sha = text(alleles_link.get("sha256"), "alleles_table.sha256")?;
    let alleles_artifact = listed(ALLELES_ID, ALLELES_KIND, ALLELES_PATH)?;
    let alleles_bytes = read(ALLELES_PATH)?;
    let alleles_sha = sha256_hex(&alleles_bytes);
    if alleles_sha != alleles_linked_sha || alleles_sha != alleles_artifact.digest {
        return Err(catalog_reject(
            ALLELES_PATH,
            "la huella de la tabla de alelos no es la del enlace de la efectiva o la del manifiesto",
        ));
    }
    if alleles_bytes.len() as u64 != alleles_artifact.size_bytes {
        return Err(catalog_reject(
            ALLELES_PATH,
            "el tamaño de la tabla de alelos no es el del manifiesto",
        ));
    }
    let metadata: Value = serde_json::from_slice(&metadata_bytes)
        .map_err(|error| catalog_reject(METADATA_PATH, format!("JSON no legible: {error}")))?;
    let metadata = object(Some(&metadata), "metadato")?;
    let header = json!({
        "document_kind": metadata.get("document_kind"),
        "schema_version": metadata.get("schema_version"),
        "coordinate_system": metadata.get("coordinate_system"),
        "algorithm_version": metadata.get("algorithm_version"),
        "table": metadata.get("table"),
    });
    let expected_header = json!({
        "document_kind": "repeat_catalog_metadata",
        "schema_version": "1",
        "coordinate_system": "base1_closed",
        "algorithm_version": detection.get("algorithm_version"),
        "table": {"path": TABLE_PATH, "sha256": table_sha},
    });
    if header != expected_header {
        return Err(catalog_reject(
            METADATA_PATH,
            "la clase, la versión, el sistema de coordenadas, la versión del algoritmo o la tabla \
             (ruta y huella) del metadato no son los esperados",
        ));
    }
    let final_fasta_id = inputs
        .effective
        .get("inserted_reference")
        .or_else(|| inputs.effective.get("reference"))
        .and_then(|reference| reference.get("fasta_artifact_id"))
        .and_then(Value::as_str);
    let final_fasta_digest = final_fasta_id
        .and_then(artifact)
        .map(|artifact| artifact.digest);
    if metadata.get("reference")
        != Some(&json!({
            "reference_id": inputs.reference_id,
            "fasta_sha256": final_fasta_digest,
        }))
    {
        return Err(catalog_reject(
            METADATA_PATH,
            "la referencia del metadato no es la referencia final publicada (identificador y huella del FASTA)",
        ));
    }
    let admission = inputs.effective.get("repeat_admission");
    let expected_criteria = json!({
        "min_copies": detection.get("min_copies"),
        "min_purity": admission.and_then(|value| value.get("min_purity")),
        "min_longest_run_units": admission.and_then(|value| value.get("min_longest_run_units")),
        "max_units": detection.get("max_units"),
        "regions": detection.get("regions"),
    });
    let criteria_value = metadata.get("criteria");
    if criteria_value != Some(&expected_criteria) {
        return Err(catalog_reject(
            METADATA_PATH,
            "los criterios del metadato no son los de la configuración efectiva",
        ));
    }
    let criteria = Criteria::read(object(criteria_value, "criteria")?)?;

    // --- La tabla ----------------------------------------------------------
    let sequences: BTreeMap<&str, (usize, &[u8])> = inputs
        .sequences
        .iter()
        .enumerate()
        .map(|(index, sequence)| (sequence.id, (index, sequence.bases)))
        .collect();
    let rows = parse_table(&table_bytes)?;
    for row in &rows {
        let (index, bases) = locate(&sequences, row.sequence_id, row.start, row.end, TABLE_PATH)?;
        match row.origin {
            Origin::Detected => {
                let [motif] = row.motifs.as_slice() else {
                    return Err(catalog_reject(
                        TABLE_PATH,
                        format!("{}: un locus detectado tiene un solo motivo", row.id),
                    ));
                };
                let tract = Tract {
                    sequence_id: row.sequence_id,
                    start: row.start,
                    end: row.end,
                    motif,
                    score: row.score.unwrap_or_default(),
                };
                check_tract(&tract, bases, &criteria, row.id, TABLE_PATH)?;
                check_identifier(row.id, index, &tract, TABLE_PATH)?;
            }
            Origin::Declared => {
                if row.motifs.iter().any(|motif| !is_acgt(motif)) {
                    return Err(catalog_reject(
                        TABLE_PATH,
                        format!(
                            "{}: un motivo declarado no es una cadena de A, C, G y T",
                            row.id
                        ),
                    ));
                }
            }
        }
    }
    let order = |sequence_id: &str| {
        sequences
            .get(sequence_id)
            .map_or(usize::MAX, |entry| entry.0)
    };
    for pair in rows.windows(2) {
        let (left, right) = (&pair[0], &pair[1]);
        let (left_key, right_key) = (
            (order(left.sequence_id), left.start),
            (order(right.sequence_id), right.start),
        );
        if left_key >= right_key {
            return Err(catalog_reject(
                TABLE_PATH,
                format!(
                    "{} y {} no van por secuencia y por inicio",
                    left.id, right.id
                ),
            ));
        }
        if left.sequence_id == right.sequence_id && right.start <= left.end {
            return Err(catalog_reject(
                TABLE_PATH,
                format!("{} y {} se solapan", left.id, right.id),
            ));
        }
    }

    // --- Tramos sustituidos y retirados --------------------------------------
    // Los declarados aparte: son pocos (los escribe el usuario) y cada
    // sustituido toca a alguno, de modo que los dos recorridos de abajo son
    // declarados × sustituidos y no crecen con los detectados.
    let declared: Vec<&Row<'_>> = rows
        .iter()
        .filter(|row| row.origin == Origin::Declared)
        .collect();
    let replaced = array(metadata.get("replaced_tracts"), "replaced_tracts")?
        .iter()
        .map(|value| read_tract(value, "replaced_tracts"))
        .collect::<Result<Vec<_>, _>>()?;
    for tract in &replaced {
        let (_, bases) = locate(
            &sequences,
            tract.sequence_id,
            tract.start,
            tract.end,
            METADATA_PATH,
        )?;
        check_tract(
            tract,
            bases,
            &criteria,
            "un tramo sustituido",
            METADATA_PATH,
        )?;
        let touches_declared = declared.iter().any(|row| {
            row.sequence_id == tract.sequence_id
                && overlaps(row.start, row.end, tract.start, tract.end)
        });
        if !touches_declared {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "el tramo sustituido {}:{}..{} no toca a ningún locus declarado",
                    tract.sequence_id, tract.start, tract.end
                ),
            ));
        }
    }
    for row in &declared {
        let exempt = row
            .motifs
            .iter()
            .all(|motif| motif.len() > MAX_MOTIF_LENGTH);
        let matched = replaced.iter().any(|tract| {
            tract.sequence_id == row.sequence_id
                && overlaps(row.start, row.end, tract.start, tract.end)
                && row
                    .motifs
                    .iter()
                    .any(|motif| is_rotation(tract.motif, motif))
        });
        if !exempt && !matched {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "el locus declarado {} no coincide con ningún tramo sustituido de un motivo \
                     rotación del suyo",
                    row.id
                ),
            ));
        }
    }
    // Una búsqueda por retirado: con una matriz suministrada de un VCF real
    // puede haber muchas filas y muchos retirados.
    let matrix: std::collections::BTreeSet<(&str, u64, &str)> =
        inputs.matrix_variants.iter().copied().collect();
    let direct_events = inputs
        .effective
        .get("direct_events")
        .and_then(Value::as_array)
        .map_or(&[][..], Vec::as_slice);
    let mut retired = Vec::new();
    for value in array(metadata.get("retired_tracts"), "retired_tracts")? {
        let entry = object(Some(value), "retired_tracts")?;
        let tract = read_tract(value, "retired_tracts")?;
        let id = text(entry.get("locus_id"), "retired_tracts.locus_id")?;
        let (index, bases) = locate(
            &sequences,
            tract.sequence_id,
            tract.start,
            tract.end,
            METADATA_PATH,
        )?;
        check_tract(&tract, bases, &criteria, id, METADATA_PATH)?;
        check_identifier(id, index, &tract, METADATA_PATH)?;
        let variant = object(entry.get("variant"), "retired_tracts.variant")?;
        let source = text(variant.get("source"), "variant.source")?;
        let variant_id = text(variant.get("id"), "variant.id")?;
        let position = unsigned(variant.get("position"), "variant.position")?;
        if position < tract.start || position > tract.end {
            return Err(catalog_reject(
                METADATA_PATH,
                format!("{id}: la variante {variant_id} en {position} no cae dentro del tramo"),
            ));
        }
        let exists = match source {
            "population_snvs" => matrix.contains(&(tract.sequence_id, position, variant_id)),
            "direct_event" => direct_events.iter().any(|event| {
                event.get("kind").and_then(Value::as_str) == Some("snp")
                    && event.get("event_id").and_then(Value::as_str) == Some(variant_id)
                    && event.get("contig").and_then(Value::as_str) == Some(tract.sequence_id)
                    && event.get("position_1based").and_then(Value::as_u64) == Some(position)
            }),
            _ => false,
        };
        if !exists {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "{id}: no hay ninguna variante {variant_id} de origen {source} en {}:{position}",
                    tract.sequence_id
                ),
            ));
        }
        retired.push(tract);
    }
    // Los aceptados del barrido son disjuntos: los detectados que quedan, los
    // que un declarado sustituye y los retirados no pueden solaparse.
    let mut accepted: Vec<(usize, u64, u64)> = rows
        .iter()
        .filter(|row| row.origin == Origin::Detected)
        .map(|row| (order(row.sequence_id), row.start, row.end))
        .chain(
            replaced
                .iter()
                .chain(&retired)
                .map(|tract| (order(tract.sequence_id), tract.start, tract.end)),
        )
        .collect();
    accepted.sort_unstable();
    if accepted
        .windows(2)
        .any(|pair| pair[0].0 == pair[1].0 && pair[1].1 <= pair[0].2)
    {
        return Err(catalog_reject(
            METADATA_PATH,
            "dos tramos aceptados (detectados, sustituidos o retirados) se solapan",
        ));
    }

    // --- Recuentos -------------------------------------------------------
    let counts = array(
        metadata.get("counts_by_motif_length"),
        "counts_by_motif_length",
    )?;
    if counts.len() != MAX_MOTIF_LENGTH {
        return Err(catalog_reject(
            METADATA_PATH,
            "counts_by_motif_length no tiene una entrada por longitud de motivo de 1 a 6",
        ));
    }
    for (index, entry) in counts.iter().enumerate() {
        let length = index + 1;
        let entry = object(Some(entry), "counts_by_motif_length")?;
        let read = |key: &str| unsigned(entry.get(key), key);
        if read("motif_length")? != length as u64 {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "la entrada {length} de counts_by_motif_length no es la de motivos de {length} pb"
                ),
            ));
        }
        let discarded = object(entry.get("discarded"), "discarded")?;
        if discarded.len() != STAGES.len()
            || !STAGES.iter().all(|stage| discarded.contains_key(*stage))
        {
            return Err(catalog_reject(
                METADATA_PATH,
                format!("motivos de {length} pb: los descartes no son los de las ocho etapas"),
            ));
        }
        let mut discarded_total: u64 = 0;
        for stage in STAGES {
            discarded_total =
                discarded_total.saturating_add(unsigned(discarded.get(stage), stage)?);
        }
        let accepted_count = read("accepted")?;
        if read("detector_intervals")? != discarded_total.saturating_add(accepted_count) {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: los intervalos del detector no son los descartes más \
                     los aceptados"
                ),
            ));
        }
        let replaced_count = read("replaced_by_declared")?;
        let retired_count = read("retired_by_variant")?;
        let without = read("without_haplotype_copies")?;
        let published = read("published")?;
        if accepted_count
            != replaced_count
                .saturating_add(retired_count)
                .saturating_add(without)
                .saturating_add(published)
        {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: los aceptados no son los sustituidos, los retirados, \
                     los que quedan sin copias y los publicados"
                ),
            ));
        }
        let listed_by_length = |tracts: &[Tract<'_>]| {
            tracts
                .iter()
                .filter(|tract| tract.motif.len() == length)
                .count() as u64
        };
        if replaced_count != listed_by_length(&replaced) {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: replaced_by_declared no es lo que lista replaced_tracts"
                ),
            ));
        }
        if retired_count != listed_by_length(&retired) {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: retired_by_variant no es lo que lista retired_tracts"
                ),
            ));
        }
        let detected_rows = rows
            .iter()
            .filter(|row| row.origin == Origin::Detected && row.motifs[0].len() == length)
            .count() as u64;
        if published != detected_rows {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: published dice {published} y la tabla tiene \
                     {detected_rows} loci detectados"
                ),
            ));
        }
        if without > 0 && (published > 0 || retired_count > 0) {
            return Err(catalog_reject(
                METADATA_PATH,
                format!(
                    "motivos de {length} pb: sin copias haplotípicas no entra ningún detectado, \
                     y hay publicados o retirados"
                ),
            ));
        }
    }

    // --- Las filas declaradas son `repeat_loci` ---------------------------
    let loci = array(inputs.effective.get("repeat_loci"), "repeat_loci")?;
    if loci.len() != declared.len() {
        return Err(catalog_reject(
            TABLE_PATH,
            format!(
                "la tabla tiene {} loci declarados y repeat_loci de la efectiva {}",
                declared.len(),
                loci.len()
            ),
        ));
    }
    for (row, locus) in declared.iter().zip(loci) {
        let expected = json!({
            "id": row.id,
            "origin": row.origin.as_str(),
            "sequence_id": row.sequence_id,
            "start": row.start,
            "end": row.end,
            "motifs": row
                .motifs
                .iter()
                .map(|motif| String::from_utf8_lossy(motif))
                .collect::<Vec<_>>(),
        });
        let published = json!({
            "id": locus.get("id"),
            "origin": locus.get("origin"),
            "sequence_id": locus.get("sequence_id"),
            "start": locus.get("start"),
            "end": locus.get("end"),
            "motifs": locus.get("motifs"),
        });
        if expected != published {
            return Err(catalog_reject(
                TABLE_PATH,
                format!(
                    "{}: la fila no es el locus de repeat_loci de la efectiva",
                    row.id
                ),
            ));
        }
    }

    // --- La tabla de alelos: los detectados, su forma y su aritmética -------
    let alleles = parse_alleles_table(&alleles_bytes)?;
    let detected: Vec<&str> = rows
        .iter()
        .filter(|row| row.origin == Origin::Detected)
        .map(|row| row.id)
        .collect();
    let mut grouped: Vec<(&str, Vec<&AlleleRow<'_>>)> = Vec::new();
    for row in &alleles {
        match grouped.last_mut() {
            Some((locus, rows)) if *locus == row.locus_id => rows.push(row),
            _ => grouped.push((row.locus_id, vec![row])),
        }
    }
    if grouped
        .iter()
        .map(|(locus, _)| *locus)
        .ne(detected.iter().copied())
    {
        return Err(catalog_reject(
            ALLELES_PATH,
            "los loci de la tabla de alelos no son los detectados del catálogo, en su orden",
        ));
    }
    let haplotype_copies = published_haplotype_copies(inputs.effective).ok_or_else(|| {
        catalog_reject(
            ALLELES_PATH,
            "la efectiva no publica las copias haplotípicas de sus poblaciones iniciales",
        )
    })?;
    for (locus, rows) in &grouped {
        let mut sequences_seen = std::collections::BTreeSet::new();
        let mut copies: u64 = 0;
        for (index, row) in rows.iter().enumerate() {
            if row.allele_index != index as u64 {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!("{locus}: los índices de alelo no van de 0 en 0 sin huecos"),
                ));
            }
            if row.origin != "generated" {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!("{locus}: un alelo de un detectado no es generated"),
                ));
            }
            // Una `N` puede estar dentro de un tramo detectado (cuenta como
            // una base que no encaja en la pureza), y los alelos la heredan.
            if row.sequence.is_empty()
                || !row
                    .sequence
                    .bytes()
                    .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T' | b'N'))
            {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!("{locus}: la secuencia de un alelo no es una cadena de A, C, G, T y N"),
                ));
            }
            if !sequences_seen.insert(row.sequence) {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!("{locus}: dos alelos tienen la misma secuencia"),
                ));
            }
            if row.applied_events.checked_add(row.unapplied_events) != Some(row.drawn_events) {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!("{locus}: los pasos aplicados y no aplicados no suman los sorteados"),
                ));
            }
            copies = copies.saturating_add(row.realized_copies);
        }
        if copies != haplotype_copies {
            return Err(catalog_reject(
                ALLELES_PATH,
                format!(
                    "{locus}: las copias realizadas suman {copies} y las copias haplotípicas \
                     publicadas son {haplotype_copies}"
                ),
            ));
        }
        for row in rows {
            let expected = (copies > 0).then(|| row.realized_copies as f64 / copies as f64);
            if row.realized_frequency.map(f64::to_bits) != expected.map(f64::to_bits) {
                return Err(catalog_reject(
                    ALLELES_PATH,
                    format!(
                        "{locus}: la frecuencia realizada de un alelo no es sus copias entre las \
                         del locus, o «.» sin copias"
                    ),
                ));
            }
        }
    }
    Ok(())
}

/// Las copias haplotípicas de las poblaciones iniciales **publicadas**:
/// individuos por ploidía, con aritmética comprobada; `None` si algo falta o
/// desborda.
fn published_haplotype_copies(effective: &Map<String, Value>) -> Option<u64> {
    effective
        .get("initial_populations")?
        .as_array()?
        .iter()
        .try_fold(0_u64, |total, population| {
            let individuals = population.get("individual_count")?.as_u64()?;
            let ploidy = population.get("ploidy")?.as_u64()?;
            total.checked_add(individuals.checked_mul(ploidy)?)
        })
}

// ---------------------------------------------------------------------------
// Lectura
// ---------------------------------------------------------------------------

fn sha256_hex(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn object<'a>(value: Option<&'a Value>, what: &str) -> Result<&'a Map<String, Value>, Diagnostic> {
    value
        .and_then(Value::as_object)
        .ok_or_else(|| catalog_reject(METADATA_PATH, format!("{what} no es un objeto")))
}

fn array<'a>(value: Option<&'a Value>, what: &str) -> Result<&'a Vec<Value>, Diagnostic> {
    value
        .and_then(Value::as_array)
        .ok_or_else(|| catalog_reject(METADATA_PATH, format!("{what} no es una lista")))
}

fn text<'a>(value: Option<&'a Value>, what: &str) -> Result<&'a str, Diagnostic> {
    value
        .and_then(Value::as_str)
        .ok_or_else(|| catalog_reject(METADATA_PATH, format!("{what} no es un texto")))
}

fn unsigned(value: Option<&Value>, what: &str) -> Result<u64, Diagnostic> {
    value
        .and_then(Value::as_u64)
        .ok_or_else(|| catalog_reject(METADATA_PATH, format!("{what} no es un entero no negativo")))
}

/// Los criterios del metadato, ya comparados con los de la efectiva.
struct Criteria {
    min_copies: Vec<u64>,
    min_purity: f64,
    min_longest_run_units: u64,
    max_units: u64,
    /// `(secuencia, inicio, fin)` en base uno; `None` es la referencia entera.
    regions: Option<Vec<(String, u64, u64)>>,
}

impl Criteria {
    fn read(criteria: &Map<String, Value>) -> Result<Self, Diagnostic> {
        let min_copies = array(criteria.get("min_copies"), "criteria.min_copies")?
            .iter()
            .map(|value| unsigned(Some(value), "criteria.min_copies"))
            .collect::<Result<Vec<_>, _>>()?;
        let min_purity = criteria
            .get("min_purity")
            .and_then(Value::as_f64)
            .ok_or_else(|| catalog_reject(METADATA_PATH, "criteria.min_purity no es un número"))?;
        if min_copies.len() != MAX_MOTIF_LENGTH {
            return Err(catalog_reject(
                METADATA_PATH,
                "criteria.min_copies no tiene seis suelos",
            ));
        }
        let regions = match criteria.get("regions") {
            Some(Value::Null) => None,
            value => Some(
                array(value, "criteria.regions")?
                    .iter()
                    .map(|region| {
                        let region = object(Some(region), "criteria.regions")?;
                        Ok((
                            text(region.get("sequence_id"), "region.sequence_id")?.to_owned(),
                            unsigned(region.get("start"), "region.start")?,
                            unsigned(region.get("end"), "region.end")?,
                        ))
                    })
                    .collect::<Result<Vec<_>, Diagnostic>>()?,
            ),
        };
        Ok(Self {
            min_copies,
            min_purity,
            min_longest_run_units: unsigned(
                criteria.get("min_longest_run_units"),
                "criteria.min_longest_run_units",
            )?,
            max_units: unsigned(criteria.get("max_units"), "criteria.max_units")?,
            regions,
        })
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Origin {
    Detected,
    Declared,
}

impl Origin {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Detected => "detected",
            Self::Declared => "declared",
        }
    }
}

struct Row<'a> {
    id: &'a str,
    sequence_id: &'a str,
    start: u64,
    end: u64,
    motifs: Vec<&'a [u8]>,
    score: Option<i64>,
    origin: Origin,
}

/// Un entero escrito en su forma canónica: sin signo, sin ceros a la izquierda.
fn canonical_integer(field: &str) -> Option<u64> {
    let value = field.parse::<u64>().ok()?;
    (value.to_string() == field).then_some(value)
}

fn parse_table(bytes: &[u8]) -> Result<Vec<Row<'_>>, Diagnostic> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| catalog_reject(TABLE_PATH, "la tabla no es UTF-8"))?;
    let body = text
        .strip_suffix('\n')
        .ok_or_else(|| catalog_reject(TABLE_PATH, "la tabla no termina en un salto de línea"))?;
    let mut lines = body.split('\n');
    if lines.next() != Some(TABLE_HEADER) {
        return Err(catalog_reject(
            TABLE_PATH,
            "la cabecera de la tabla no es la del contrato",
        ));
    }
    let mut rows = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let [id, sequence_id, start, end, motifs, score, origin] = fields.as_slice() else {
            return Err(catalog_reject(
                TABLE_PATH,
                format!("una fila no tiene siete columnas: {line:?}"),
            ));
        };
        let (Some(start), Some(end)) = (canonical_integer(start), canonical_integer(end)) else {
            return Err(catalog_reject(
                TABLE_PATH,
                format!("{id}: el inicio o el fin no es un entero canónico"),
            ));
        };
        let origin = match *origin {
            "detected" => Origin::Detected,
            "declared" => Origin::Declared,
            other => {
                return Err(catalog_reject(
                    TABLE_PATH,
                    format!("{id}: origen desconocido {other:?}"),
                ));
            }
        };
        let score = match (origin, *score) {
            (Origin::Declared, ".") => None,
            (Origin::Detected, value)
                if value
                    .parse::<i64>()
                    .is_ok_and(|parsed| parsed.to_string() == value) =>
            {
                value.parse::<i64>().ok()
            }
            _ => {
                return Err(catalog_reject(
                    TABLE_PATH,
                    format!(
                        "{id}: la puntuación de un detectado es un entero y la de un declarado es «.»"
                    ),
                ));
            }
        };
        rows.push(Row {
            id,
            sequence_id,
            start,
            end,
            motifs: motifs.split(',').map(str::as_bytes).collect(),
            score,
            origin,
        });
    }
    Ok(rows)
}

/// Una fila de la tabla de alelos de los detectados (DTC-104).
pub(crate) struct AlleleRow<'a> {
    pub(crate) locus_id: &'a str,
    pub(crate) allele_index: u64,
    pub(crate) sequence: &'a str,
    pub(crate) origin: &'a str,
    pub(crate) declared_frequency: f64,
    pub(crate) realized_copies: u64,
    /// `None` es el `.` de un locus sin copias sorteadas.
    pub(crate) realized_frequency: Option<f64>,
    pub(crate) drawn_events: u64,
    pub(crate) applied_events: u64,
    pub(crate) unapplied_events: u64,
    pub(crate) vcf_tract_bases_bound: u64,
}

/// Una proporción escrita como decimal: finita y en `[0, 1]`.
fn proportion(field: &str) -> Option<f64> {
    field
        .parse::<f64>()
        .ok()
        .filter(|value| value.is_finite() && (0.0..=1.0).contains(value))
}

/// Lee la tabla de alelos con su forma estricta: cabecera exacta, trece
/// columnas, enteros canónicos, proporciones decimales en `[0, 1]` y `.` solo
/// en la frecuencia realizada. La usan también los jueces de `src/results.rs`
/// que la leen, para que la tabla se lea de una sola manera.
pub(crate) fn parse_alleles_table(bytes: &[u8]) -> Result<Vec<AlleleRow<'_>>, Diagnostic> {
    let text = std::str::from_utf8(bytes)
        .map_err(|_| catalog_reject(ALLELES_PATH, "la tabla de alelos no es UTF-8"))?;
    let body = text.strip_suffix('\n').ok_or_else(|| {
        catalog_reject(
            ALLELES_PATH,
            "la tabla de alelos no termina en un salto de línea",
        )
    })?;
    let mut lines = body.split('\n');
    if lines.next() != Some(ALLELES_HEADER) {
        return Err(catalog_reject(
            ALLELES_PATH,
            "la cabecera de la tabla de alelos no es la del contrato",
        ));
    }
    let mut rows = Vec::new();
    for line in lines {
        let fields: Vec<&str> = line.split('\t').collect();
        let [
            locus_id,
            allele_index,
            sequence,
            origin,
            declared_frequency,
            effective_frequency,
            realized_copies,
            realized_frequency,
            drawn_events,
            applied_events,
            unapplied_events,
            truth_tract_bases_bound,
            vcf_tract_bases_bound,
        ] = fields.as_slice()
        else {
            return Err(catalog_reject(
                ALLELES_PATH,
                format!(
                    "{}: una fila de la tabla de alelos no tiene trece columnas",
                    fields[0]
                ),
            ));
        };
        let integers = [
            allele_index,
            realized_copies,
            drawn_events,
            applied_events,
            unapplied_events,
            truth_tract_bases_bound,
            vcf_tract_bases_bound,
        ]
        .map(|field| canonical_integer(field));
        let [
            Some(allele_index),
            Some(realized_copies),
            Some(drawn_events),
            Some(applied_events),
            Some(unapplied_events),
            Some(_truth_tract_bases_bound),
            Some(vcf_tract_bases_bound),
        ] = integers
        else {
            return Err(catalog_reject(
                ALLELES_PATH,
                format!("{locus_id}: una columna entera no es un entero canónico"),
            ));
        };
        let (Some(declared_frequency), Some(_effective_frequency)) = (
            proportion(declared_frequency),
            proportion(effective_frequency),
        ) else {
            return Err(catalog_reject(
                ALLELES_PATH,
                format!("{locus_id}: una frecuencia no es un decimal en [0, 1]"),
            ));
        };
        let realized_frequency = match *realized_frequency {
            "." => None,
            value => Some(proportion(value).ok_or_else(|| {
                catalog_reject(
                    ALLELES_PATH,
                    format!(
                        "{locus_id}: la frecuencia realizada no es un decimal en [0, 1] ni «.»"
                    ),
                )
            })?),
        };
        rows.push(AlleleRow {
            locus_id,
            allele_index,
            sequence,
            origin,
            declared_frequency,
            realized_copies,
            realized_frequency,
            drawn_events,
            applied_events,
            unapplied_events,
            vcf_tract_bases_bound,
        });
    }
    Ok(rows)
}

/// Un locus de la tabla del catálogo publicada, como lo leen los jueces de
/// `src/results.rs` que atribuyen registros a loci o rehacen el abanico de un
/// detectado.
pub(crate) struct CatalogLocus<'a> {
    pub(crate) id: &'a str,
    pub(crate) sequence_id: &'a str,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) motifs: Vec<&'a [u8]>,
    pub(crate) detected: bool,
}

/// Los loci de la tabla del catálogo publicada, leídos con su forma estricta.
pub(crate) fn published_catalog_loci(bytes: &[u8]) -> Result<Vec<CatalogLocus<'_>>, Diagnostic> {
    Ok(parse_table(bytes)?
        .into_iter()
        .map(|row| CatalogLocus {
            id: row.id,
            sequence_id: row.sequence_id,
            start: row.start,
            end: row.end,
            motifs: row.motifs,
            detected: row.origin == Origin::Detected,
        })
        .collect())
}

/// Un tramo detectado: un locus detectado de la tabla, un sustituido o un
/// retirado. Base uno, cerrado.
struct Tract<'a> {
    sequence_id: &'a str,
    start: u64,
    end: u64,
    motif: &'a [u8],
    score: i64,
}

fn read_tract<'a>(value: &'a Value, what: &str) -> Result<Tract<'a>, Diagnostic> {
    let entry = object(Some(value), what)?;
    let score = entry.get("score").and_then(Value::as_i64).ok_or_else(|| {
        catalog_reject(
            METADATA_PATH,
            format!("{what}: la puntuación no es un entero"),
        )
    })?;
    Ok(Tract {
        sequence_id: text(entry.get("sequence_id"), "sequence_id")?,
        start: unsigned(entry.get("start"), "start")?,
        end: unsigned(entry.get("end"), "end")?,
        motif: text(entry.get("motif"), "motif")?.as_bytes(),
        score,
    })
}

/// La secuencia del tramo y su índice, con el tramo dentro de ella.
fn locate<'a>(
    sequences: &BTreeMap<&str, (usize, &'a [u8])>,
    sequence_id: &str,
    start: u64,
    end: u64,
    file: &str,
) -> Result<(usize, &'a [u8]), Diagnostic> {
    let &(index, bases) = sequences.get(sequence_id).ok_or_else(|| {
        catalog_reject(
            file,
            format!("la secuencia {sequence_id} no está en el FASTA final publicado"),
        )
    })?;
    if start == 0 || start > end || end > bases.len() as u64 {
        return Err(catalog_reject(
            file,
            format!("{sequence_id}:{start}..{end} no cabe en su secuencia"),
        ));
    }
    Ok((index, bases))
}

// ---------------------------------------------------------------------------
// Aritmética propia
// ---------------------------------------------------------------------------

fn is_acgt(motif: &[u8]) -> bool {
    !motif.is_empty()
        && motif
            .iter()
            .all(|base| matches!(base, b'A' | b'C' | b'G' | b'T'))
}

/// Primitivo: no es la repetición de un trozo más corto (`AA` no lo es).
fn is_primitive(motif: &[u8]) -> bool {
    let length = motif.len();
    !(1..length).any(|divisor| {
        length.is_multiple_of(divisor)
            && motif
                .chunks(divisor)
                .all(|chunk| chunk == &motif[..divisor])
    })
}

/// Rotación: misma longitud y aparece dentro de `motif` repetido dos veces.
fn is_rotation(candidate: &[u8], motif: &[u8]) -> bool {
    candidate.len() == motif.len()
        && [motif, motif]
            .concat()
            .windows(candidate.len())
            .any(|window| window == candidate)
}

/// Dos intervalos cerrados que comparten al menos una base.
const fn overlaps(start: u64, end: u64, other_start: u64, other_end: u64) -> bool {
    start <= other_end && other_start <= end
}

/// Descomposición voraz de izquierda a derecha en copias completas del
/// motivo: `(bases cubiertas, tirada más larga en copias)`. Una copia es el
/// motivo exacto; lo que no encaja avanza una base y corta la tirada.
fn decompose(tract: &[u8], motif: &[u8]) -> (u64, u64) {
    let period = motif.len();
    let (mut offset, mut covered, mut run, mut longest) = (0, 0_u64, 0_u64, 0_u64);
    while offset < tract.len() {
        if tract.get(offset..offset + period) == Some(motif) {
            covered += period as u64;
            run += 1;
            longest = longest.max(run);
            offset += period;
        } else {
            run = 0;
            offset += 1;
        }
    }
    (covered, longest)
}

/// La puntuación de un tramo con sus propias comparaciones: cada base desde
/// la segunda copia, contra la de un periodo antes.
fn own_score(tract: &[u8], period: usize) -> i64 {
    (period..tract.len())
        .map(|offset| {
            let (here, back) = (tract[offset], tract[offset - period]);
            if here == back && is_acgt(&[here]) {
                MATCH_REWARD
            } else {
                -MISMATCH_PENALTY
            }
        })
        .sum()
}

/// Cada criterio de la detección sobre un tramo, contra las bases publicadas.
fn check_tract(
    tract: &Tract<'_>,
    bases: &[u8],
    criteria: &Criteria,
    name: &str,
    file: &str,
) -> Result<(), Diagnostic> {
    let at = format!(
        "{name} ({}:{}..{})",
        tract.sequence_id, tract.start, tract.end
    );
    let motif = tract.motif;
    let period = motif.len();
    if period > MAX_MOTIF_LENGTH || !is_acgt(motif) {
        return Err(catalog_reject(
            file,
            format!("{at}: el motivo no tiene de 1 a 6 bases A, C, G o T"),
        ));
    }
    if !is_primitive(motif) {
        return Err(catalog_reject(
            file,
            format!("{at}: el motivo no es primitivo"),
        ));
    }
    let (begin, end) = ((tract.start - 1) as usize, tract.end as usize);
    let span = &bases[begin..end];
    if !span.starts_with(motif) {
        return Err(catalog_reject(
            file,
            format!("{at}: el motivo no es la lectura de las bases del inicio del tramo"),
        ));
    }
    let flanked = begin > 0
        && end < bases.len()
        && is_acgt(&bases[begin - 1..begin])
        && is_acgt(&bases[end..=end]);
    if !flanked {
        return Err(catalog_reject(
            file,
            format!("{at}: no tiene una base A, C, G o T a cada lado"),
        ));
    }
    if span.len() as u64 > criteria.max_units.saturating_mul(period as u64) {
        return Err(catalog_reject(
            file,
            format!("{at}: rebasa el tope de {} unidades", criteria.max_units),
        ));
    }
    let (covered, longest) = decompose(span, motif);
    let floor = criteria.min_copies[period - 1];
    let copies = covered / period as u64;
    if copies < floor {
        return Err(catalog_reject(
            file,
            format!(
                "{at}: tiene {copies} copias completas y el suelo de motivos de {period} pb es {floor}"
            ),
        ));
    }
    if (covered as f64) / (span.len() as f64) < criteria.min_purity {
        return Err(catalog_reject(
            file,
            format!(
                "{at}: su pureza {covered}/{} no alcanza {}",
                span.len(),
                criteria.min_purity
            ),
        ));
    }
    if longest < criteria.min_longest_run_units {
        return Err(catalog_reject(
            file,
            format!(
                "{at}: su tirada más larga es de {longest} copias y el mínimo es {}",
                criteria.min_longest_run_units
            ),
        ));
    }
    if let Some(regions) = &criteria.regions
        && !regions.iter().any(|(sequence_id, start, end)| {
            sequence_id == tract.sequence_id && *start <= tract.start && tract.end <= *end
        })
    {
        return Err(catalog_reject(
            file,
            format!("{at}: no cabe entero en ninguna región de detección"),
        ));
    }
    let score = own_score(span, period);
    if score != tract.score {
        return Err(catalog_reject(
            file,
            format!(
                "{at}: su puntuación publicada es {} y la de sus comparaciones {score}",
                tract.score
            ),
        ));
    }
    Ok(())
}

/// El identificador de un detectado sale de sus coordenadas.
fn check_identifier(
    id: &str,
    index: usize,
    tract: &Tract<'_>,
    file: &str,
) -> Result<(), Diagnostic> {
    let expected = format!(
        "SSR-{}-{}-{}-{}",
        index + 1,
        tract.start,
        tract.end,
        String::from_utf8_lossy(tract.motif)
    );
    if id != expected {
        return Err(catalog_reject(
            file,
            format!("{id}: el identificador de un detectado es {expected}"),
        ));
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn primitive_rotation_and_decomposition_are_computed_here() {
        assert!(is_primitive(b"CA"));
        assert!(!is_primitive(b"CACA"));
        assert!(!is_primitive(b"AA"));
        assert!(is_primitive(b"AAC"));
        assert!(is_rotation(b"CA", b"AC"));
        assert!(!is_rotation(b"CC", b"AC"));
        assert_eq!(decompose(b"CACAGCACA", b"CA"), (8, 2));
        assert_eq!(own_score(b"AAAA", 1), 6);
        assert_eq!(own_score(b"AANA", 1), 2 - 7 - 7);
    }
}
