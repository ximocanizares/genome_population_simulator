//! **La sección de microsatélites para las pruebas escritas antes de la
//! detección** (entrega 3 del plan de detección, paso F6.4.a.13b).
//!
//! Desde el 29/10/2026 un locus declarado solo existe con la sección
//! `[reference.repeat_detection]` (REQ-CONF-142), y entonces tiene que
//! coincidir con un tramo detectado (REQ-REF-40). Las pruebas de antes
//! declaraban sus loci sin ella. Esta función les añade la sección de modo que
//! **el catálogo sea exactamente el de antes**: los loci declarados, con sus
//! fronteras, y ninguno más. Dos cosas lo procuran:
//!
//! 1. **Una región por locus declarado** (en coordenadas de la referencia
//!    final; con inserción, ya traducidas) que cubre su tramo **ensanchado
//!    `MARGIN` bases a cada lado**. Hace falta margen porque el detector da la
//!    repetición que **hay** en la referencia, y en varias pruebas es algo más
//!    larga que el locus declarado: una copia parcial en el flanco (`(CA)×10`
//!    seguido de una `C` sale de 21 bases), o una o dos copias más que el fondo
//!    al azar o el escenario de la prueba ponen al lado (los escenarios de la
//!    revisión de la fusión del VCF declaran a propósito loci más cortos que
//!    su tracto). Seis bases cubren todos los casos de la batería. El locus
//!    declarado coincide con esa repetición y **mandan sus fronteras**
//!    (REQ-REF-40). Las regiones que se tocan o se solapan se funden en una.
//! 2. **Suelos bajos en las longitudes de motivo que declara algún locus**:
//!    la mitad de las copias del locus de un motivo más corto de esa longitud,
//!    y una copia si la declara un compuesto. Varias pruebas declaran loci
//!    cortos a propósito (por debajo del suelo por omisión) y lo que prueban
//!    no es la detección. Las demás longitudes llevan un suelo
//!    **inalcanzable** (`UNUSED_FLOOR`), para que ningún tramo de otro periodo
//!    entre en el catálogo por el margen.
//!
//! **Lo que no garantiza, dicho sin adornos:** un tramo ajeno de una longitud
//! de motivo declarada que quepa entero en el margen, o entre dos loci cuyas
//! regiones se funden, y alcance ese suelo entraría en el catálogo como
//! detectado (con suelo de una copia, la `G×4` que separa el `A×15` y el
//! `(CA)×17` de la evidencia de DT-009 lo hacía). Se comprobó el 30/10/2026
//! ejecutando la batería entera con un registro de cada locus detectado que
//! queda en el catálogo; donde pasaba, la prueba declara su propia sección.
//!
//! El deslizamiento común no lo usa ningún locus de esas pruebas salvo que lo
//! pida: todas declaran sus alelos o su generador.
//!
//! No toca un documento sin loci, uno que ya declare la sección, uno de una
//! versión que no sea la 4 (que no la conoce) ni un texto que no sea TOML.
#![allow(dead_code)]

/// Suelo de copias de las longitudes de motivo que no declara ningún locus:
/// mayor que el tope de unidades de un detectado (200), de modo que ningún
/// tramo de esa longitud lo alcanza.
pub const UNUSED_FLOOR: u64 = 1_000;

/// Bases de margen a cada lado de cada locus declarado (ver la cabecera).
pub const MARGIN: u64 = 6;

/// Añade la sección (ver la cabecera).
pub fn with_detection_section(toml: &str) -> String {
    use toml::de::{DeTable, DeValue};
    // Lo que no es TOML (un FASTA, una matriz) pasa sin tocar: así la función
    // puede envolver cualquier escritura de una prueba.
    let Ok(document) = DeTable::parse(toml) else {
        return toml.to_owned();
    };
    let document = document.get_ref();
    // Solo el general v4 conoce la sección; en un documento de otra versión,
    // que la prueba escribe para ver su rechazo, no se toca nada.
    if document
        .get("schema_version")
        .and_then(|value| value.get_ref().as_str())
        .is_some_and(|version| version != "4")
    {
        return toml.to_owned();
    }
    let Some(reference) = document
        .get("reference")
        .and_then(|value| value.get_ref().as_table())
    else {
        return toml.to_owned();
    };
    if reference.contains_key("repeat_detection") {
        return toml.to_owned();
    }
    let Some(loci) = reference
        .get("repeat_loci")
        .and_then(|value| value.get_ref().as_array())
        .filter(|loci| !loci.is_empty())
    else {
        return toml.to_owned();
    };
    let as_u64 = |value: &DeValue<'_>| {
        value
            .as_integer()
            .and_then(|integer| u64::from_str_radix(integer.as_str(), integer.radix()).ok())
    };
    let integer =
        |table: &DeTable<'_>, key: &str| table.get(key).and_then(|value| as_u64(value.get_ref()));
    let text = |table: &DeTable<'_>, key: &str| {
        table
            .get(key)
            .and_then(|value| value.get_ref().as_str())
            .map(str::to_owned)
    };
    let motifs = |table: &DeTable<'_>| -> Vec<String> {
        text(table, "motif").map_or_else(
            || {
                // Sin motivos legibles (un locus mal formado a propósito), ninguno.
                table
                    .get("motifs")
                    .and_then(|value| value.get_ref().as_array())
                    .map(|motifs| {
                        motifs
                            .iter()
                            .filter_map(|motif| motif.get_ref().as_str().map(str::to_owned))
                            .collect()
                    })
                    .unwrap_or_default()
            },
            |motif| vec![motif],
        )
    };
    let tables: Vec<&DeTable<'_>> = loci
        .iter()
        .filter_map(|locus| locus.get_ref().as_table())
        .collect();
    // Los puntos de inserción: (secuencia, base de partida, bases insertadas).
    let mut points: Vec<(String, u64, u64)> = Vec::new();
    for locus in &tables {
        if let Some(anchor) = integer(locus, "insert_after") {
            // Sin `reference_units` legibles (una inserción mal formada a
            // propósito) no se sabe cuánto mide: no hay región que declarar, y
            // la prueba espera el rechazo del lector.
            let Some(units) = locus
                .get("reference_units")
                .and_then(|value| value.get_ref().as_array())
                .and_then(|units| {
                    units
                        .iter()
                        .map(|value| as_u64(value.get_ref()))
                        .collect::<Option<Vec<u64>>>()
                })
            else {
                continue;
            };
            let length = motifs(locus)
                .iter()
                .zip(units)
                .map(|(motif, units)| motif.len() as u64 * units)
                .sum();
            let Some(sequence) = text(locus, "sequence") else {
                continue;
            };
            points.push((sequence, anchor, length));
        }
    }
    let shift = |sequence: &str, before: u64| -> u64 {
        points
            .iter()
            .filter(|(id, anchor, _)| id == sequence && *anchor < before)
            .map(|(_, _, length)| length)
            .sum()
    };
    // El suelo de cada longitud de motivo: la **mitad** de las copias del
    // locus de un motivo más corto de esa longitud (al menos una), y una sola
    // copia si la declara un compuesto. Basta para que cada locus encuentre
    // su tramo, incluso impuro, y deja fuera las tiras cortas de los flancos.
    let mut min_copies = [UNUSED_FLOOR; 6];
    for locus in &tables {
        let motifs = motifs(locus);
        let span = match integer(locus, "insert_after") {
            Some(_) => locus
                .get("reference_units")
                .and_then(|value| value.get_ref().as_array())
                .and_then(|units| {
                    units
                        .iter()
                        .zip(&motifs)
                        .map(|(units, motif)| Some(as_u64(units.get_ref())? * motif.len() as u64))
                        .sum::<Option<u64>>()
                }),
            None => integer(locus, "start")
                .zip(integer(locus, "end"))
                .and_then(|(start, end)| (end + 1).checked_sub(start)),
        };
        for motif in &motifs {
            let length = motif.len();
            if (1..=6).contains(&length) {
                let floor = match (motifs.len(), span) {
                    (1, Some(span)) => (span / length as u64 / 2).max(1),
                    _ => 1,
                };
                min_copies[length - 1] = min_copies[length - 1].min(floor);
            }
        }
    }
    let mut section = format!(
        "\n[reference.repeat_detection]\nmin_copies = {min_copies:?}\n\n\
         [reference.repeat_detection.slippage]\nmodel = \"stepwise\"\nmean_events = 1.0\n\n"
    );
    // Longitud de cada secuencia en la referencia final: la declarada más lo
    // que se le inserta.
    let lengths: Vec<(String, u64)> = reference
        .get("sequences")
        .and_then(|value| value.get_ref().as_array())
        .map(|sequences| {
            sequences
                .iter()
                .filter_map(|sequence| sequence.get_ref().as_table())
                .filter_map(|sequence| Some((text(sequence, "id")?, integer(sequence, "length")?)))
                .map(|(id, length)| {
                    let inserted: u64 = points
                        .iter()
                        .filter(|(point, _, _)| *point == id)
                        .map(|(_, _, length)| length)
                        .sum();
                    (id, length + inserted)
                })
                .collect()
        })
        .unwrap_or_default();
    let mut regions: Vec<(String, u64, u64)> = Vec::new();
    for locus in &tables {
        let Some(sequence) = text(locus, "sequence") else {
            continue;
        };
        let (start, end) = match integer(locus, "insert_after") {
            Some(anchor) => {
                let Some(length) = points
                    .iter()
                    .find(|(id, point, _)| *id == sequence && *point == anchor)
                    .map(|(_, _, length)| *length)
                else {
                    continue;
                };
                let start = anchor + shift(&sequence, anchor) + 1;
                (start, start + length - 1)
            }
            None => {
                let (Some(start), Some(end)) = (integer(locus, "start"), integer(locus, "end"))
                else {
                    // Un locus mal formado a propósito: la prueba espera el
                    // rechazo del lector, y no hay región que declarar.
                    continue;
                };
                (start + shift(&sequence, start), end + shift(&sequence, end))
            }
        };
        let length = lengths
            .iter()
            .find(|(id, _)| *id == sequence)
            .map_or(u64::MAX, |(_, length)| *length);
        regions.push((
            sequence,
            start.saturating_sub(MARGIN).max(1),
            (end + MARGIN).min(length),
        ));
    }
    regions.sort();
    let mut merged: Vec<(String, u64, u64)> = Vec::new();
    for (sequence, start, end) in regions {
        match merged.last_mut() {
            Some(last) if last.0 == sequence && start <= last.2 + 1 => last.2 = last.2.max(end),
            _ => merged.push((sequence, start, end)),
        }
    }
    // Sin ningún locus con intervalo (uno mal formado a propósito) no hay
    // región que declarar y la sección va sin lista, que es la referencia
    // entera: los suelos inalcanzables de las longitudes que nadie declara
    // la dejan sin nada que detectar fuera de lo declarado, y la prueba ve el
    // rechazo del lector que esperaba.
    for (sequence, start, end) in merged {
        section.push_str(&format!(
            "[[reference.repeat_detection.regions]]\nsequence = \"{sequence}\"\nstart = {start}\n\
             end = {end}\n\n"
        ));
    }
    let mut out = toml.to_owned();
    if !out.ends_with('\n') {
        out.push('\n');
    }
    out.push_str(&section);
    out
}

/// **Un fondo de `G` y `T` sin tiras de más de dos bases iguales**,
/// reproducible (xorshift64 con la semilla dada), para las pruebas que
/// rodeaban su tracto de una tira larga de una sola base.
///
/// # Por qué hace falta
///
/// El detector no ve un microsatélite **corto entre dos tiras largas de una
/// misma base**: con `G×116 (CA)×13 G×122` el tramo de periodo 1 cruza el
/// `(CA)×13` como una interrupción, la suma de las dos tiras paga de sobra sus
/// discrepancias, y el `(CA)×13` no sale como tramo propio (comprobado el
/// 30/10/2026 con el núcleo de `fa215c23`; con `G×60` a la izquierda sí
/// sale). Un locus declarado ahí no puede coincidir con nada. Este fondo
/// conserva lo que esas pruebas usaban de la tira de `G`: no lleva ni `A` ni
/// `C`, de modo que un motivo con esas bases solo aparece en el tracto, y
/// contar sus copias en una lectura sigue contando solo las del tracto; y no
/// lleva tres bases iguales seguidas, de modo que un marcador como `TTTT`
/// colocado junto al tracto sigue siendo el único de la secuencia.
pub fn gt_background(seed: u64, length: usize) -> String {
    let mut state = seed;
    let mut bases: Vec<u8> = Vec::with_capacity(length);
    for _ in 0..length {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let mut base = if state >> 63 == 0 { b'G' } else { b'T' };
        if let [.., before, last] = bases.as_slice()
            && before == last
            && *last == base
        {
            base = if base == b'G' { b'T' } else { b'G' };
        }
        bases.push(base);
    }
    String::from_utf8(bases).expect("ASCII")
}
