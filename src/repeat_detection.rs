//! Detección de microsatélites (repeticiones en tándem de periodo corto) sobre
//! una secuencia de bases, y conversión de lo detectado en candidatos a locus.
//!
//! # Procedencia y licencia
//!
//! El **núcleo** de este módulo —[`PeriodRange`], [`ScanParams`],
//! [`RepeatInterval`], la tabla `CANONICAL`, `is_primitive_motif`, el pase de
//! Ruzzo y Tompa y [`find_tandem_repeats`], con sus pruebas— está **copiado**
//! del fichero `src/tandem_repeat.rs` del proyecto `pop_var_caller`
//! (repositorio `join_per_sample_vcfs`), commit de origen
//! `950909f9962f463ca157264fa36193703ed2b467`.
//!
//! - **Licencia del original:** MIT, «Copyright (c) 2026 Jose Blanca & Joaquin
//!   Cañizares», los mismos autores que este proyecto. Se copia sin depender del
//!   *crate* (decisión 8 del responsable; DTC-103).
//! - **Nota de sala limpia del original, literal:**
//!
//!   > **Clean-room note:** derived only from Benson (1999)'s published idea
//!   > (compare to a period-shifted copy; grow runs of matches) and our own code
//!   > — **not** from the AGPL-v3 `TRF-mod` source, which was not read.
//!
//!   Es decir: derivado solo de la idea publicada por Benson (1999) —comparar
//!   con una copia desplazada un periodo y hacer crecer las tiradas de
//!   coincidencias— y de código propio; **no** del código AGPL-v3 de
//!   `TRF-mod`, que no se leyó.
//! - **Qué se quitó del original:** la dependencia `thiserror` (el error de
//!   [`PeriodRange`] escribe su `Display` a mano), el lector de FASTA del otro
//!   proyecto (`crate::fasta`) y su error de lectura (`ScanError`), toda la
//!   parte de regiones y ventanas (`SegmentOptions`, `Region`,
//!   `RegionScanner`, `scan_windowed` y compañía) y la función muerta
//!   `bases_match`. Los comentarios se han traducido; el algoritmo y sus
//!   pruebas no se han tocado.
//!
//! # La conversión en candidatos (propia de este proyecto)
//!
//! [`detect_repeat_candidates`] aplica a los tramos del detector los criterios
//! de REQ-REF-39 y cuenta, **por etapa y por longitud de motivo**, dónde cae
//! cada tramo descartado. Las reglas son las recomendaciones de DP-SSR-02 a
//! DP-SSR-08 del plan de detección y son **provisionales**: se fijan en la
//! pausa de decisión de la entrega 1, con las cifras medidas. Por eso los
//! umbrales son parámetros ([`DetectionCriteria`]) y no constantes.
//!
//! Este módulo **no está conectado** a ninguna ejecución: lo usa solo la
//! medición de la entrega 1 (`tests/repeat_detection_measurement.rs`).

use std::collections::BTreeMap;
use std::fmt;

use crate::application::{DEFAULT_REPEAT_MAX_UNITS, EffectiveRepeatAdmission};
use crate::tandem_repeat::{RepeatMotif, observe_tract};

// ---------------------------------------------------------------------
// Entradas: el alcance y la puntuación (copiado)
// ---------------------------------------------------------------------

/// Intervalo inclusivo de **periodos** (longitudes de motivo, en pb) que se
/// buscan.
///
/// Lo elige quien llama; el algoritmo no impone techo propio. Es el único tipo
/// **restringido** del núcleo —un periodo 0 no significa nada y un intervalo
/// vacío (`min > max`) es un error de quien llama—, así que sus campos son
/// privados y todo valor pasa por [`PeriodRange::new`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PeriodRange {
    min: u8,
    max: u8,
}

impl PeriodRange {
    /// Construye el intervalo y rechaza `min == 0` y `min > max`. Son errores de
    /// quien llama (el periodo 0 no tiene sentido; un intervalo vacío no busca
    /// nada), así que fallan en voz alta en vez de corregirse en silencio.
    pub fn new(min: u8, max: u8) -> Result<Self, PeriodRangeError> {
        if min == 0 {
            return Err(PeriodRangeError::ZeroMin);
        }
        if min > max {
            return Err(PeriodRangeError::MinExceedsMax { min, max });
        }
        Ok(Self { min, max })
    }

    /// El periodo más pequeño que se busca (inclusivo).
    #[inline]
    #[must_use]
    pub fn min(self) -> u8 {
        self.min
    }

    /// El periodo más grande que se busca (inclusivo).
    #[inline]
    #[must_use]
    pub fn max(self) -> u8 {
        self.max
    }
}

/// No se pudo construir un [`PeriodRange`]: siempre un error de quien llama,
/// nunca una condición de los datos.
#[non_exhaustive]
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PeriodRangeError {
    /// `min` valía 0; el periodo 0 no significa nada.
    ZeroMin,
    /// `min` superaba a `max`; el intervalo no buscaría nada.
    MinExceedsMax { min: u8, max: u8 },
}

impl fmt::Display for PeriodRangeError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ZeroMin => {
                formatter.write_str("period range min must be >= 1 (period 0 is meaningless)")
            }
            Self::MinExceedsMax { min, max } => {
                write!(formatter, "period range min ({min}) exceeds max ({max})")
            }
        }
    }
}

impl std::error::Error for PeriodRangeError {}

/// Intervalo de periodos general del detector, de 1 a 6: una constante con
/// nombre y no un literal suelto.
pub const DEFAULT_PERIODS: (u8, u8) = (1, 6);

/// Recompensa por coincidencia: una coincidencia a distancia `p` suma esto a la
/// puntuación.
pub const DEFAULT_MATCH_REWARD: i32 = 2;
/// Penalización por discrepancia: una discrepancia a distancia `p` (o una base
/// que no sea A, C, G o T) resta esto. El cociente `penalización / recompensa`
/// es la tolerancia a discrepancias: un cociente `r` sostiene un tracto hasta
/// una pureza de `r / (1 + r)`; 7/2 = 3,5 lo sostiene hasta ≈ 0,78.
pub const DEFAULT_MISMATCH_PENALTY: i32 = 7;
/// Mínimo de copias por omisión para emitir un segmento: una repetición
/// necesita al menos dos copias para existir.
pub const DEFAULT_MIN_COPIES: u32 = 2;

/// La puntuación de la autocomparación a distancia `p` y el suelo mínimo de
/// copias para emitir.
///
/// Una repetición en tándem de periodo `p` aparece como una tirada de alta
/// puntuación de la señal que suma [`match_reward`](Self::match_reward) donde
/// `seq[j] == seq[j-p]` (las dos A/C/G/T, sin distinguir mayúsculas) y resta
/// [`mismatch_penalty`](Self::mismatch_penalty) en otro caso. `Default` son los
/// pesos de partida del detector (`2 / 7 / 2`). En este proyecto son
/// **constantes del algoritmo** (DTC-103): la conversión en candidatos usa
/// siempre los de omisión y no los expone como parámetro.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ScanParams {
    /// Puntuación sumada por coincidencia. Se espera `> 0`.
    pub match_reward: i32,
    /// Puntuación restada por discrepancia (o base no ACGT). Se espera `> 0`.
    pub mismatch_penalty: i32,
    /// Un segmento más corto que `min_copies · periodo` no se emite. Es un suelo
    /// permisivo: el consumidor aplica después sus suelos por periodo.
    pub min_copies: u32,
}

impl Default for ScanParams {
    fn default() -> Self {
        Self {
            match_reward: DEFAULT_MATCH_REWARD,
            mismatch_penalty: DEFAULT_MISMATCH_PENALTY,
            min_copies: DEFAULT_MIN_COPIES,
        }
    }
}

// ---------------------------------------------------------------------
// Salida del detector (copiado)
// ---------------------------------------------------------------------

/// Un intervalo de repetición en tándem detectado.
///
/// Coordenadas **en base cero y semiabiertas** (`[start, end)`), que indexan
/// directamente la secuencia que recibió [`find_tandem_repeats`]; `period` es
/// la distancia a la que se encontró; `score` es el total del segmento de
/// Ruzzo y Tompa, una medida conjunta de longitud y pureza. **No lleva
/// motivo**: el motivo es `seq[start..start + period]`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct RepeatInterval {
    /// Inicio del tramo, inclusivo (base cero).
    pub start: u64,
    /// Fin del tramo, exclusivo (base cero).
    pub end: u64,
    /// Periodo de la repetición (longitud del motivo, en pb).
    pub period: u8,
    /// Puntuación del segmento de Ruzzo y Tompa.
    pub score: i32,
}

// ---------------------------------------------------------------------
// El buscador de intervalos: puntuación a distancia p + Ruzzo–Tompa (copiado)
// ---------------------------------------------------------------------

/// Si `motif` es **primitivo**: no se forma repitiendo un periodo más corto que
/// divida su longitud. `AAA` no lo es (se forma con `A`); `ATAT` tampoco (con
/// `AT`); `AT` y `CAG` sí. Un motivo de periodo 1 siempre lo es. No distingue
/// mayúsculas, igual que la comparación del barrido.
///
/// Es lo que hace que [`find_tandem_repeats`] respete su intervalo de periodos
/// sobre el periodo **verdadero** de cada tracto: un motivo no primitivo es un
/// alias de una repetición de periodo menor, y la menor es la que se pidió.
fn is_primitive_motif(motif: &[u8]) -> bool {
    let p = motif.len();
    for q in 1..p {
        if !p.is_multiple_of(q) {
            continue; // solo los *divisores* propios de p pueden formarlo
        }
        if (q..p).all(|i| CANONICAL[motif[i] as usize] == CANONICAL[motif[i - q] as usize]) {
            return false; // se forma con sus q primeras bases → no es primitivo
        }
    }
    true
}

/// Base canónica de cada byte: `A`/`C`/`G`/`T` en cualquier caja pasan a la
/// mayúscula, y **todo lo demás pasa a `0`**, un centinela al que ninguna base
/// canónica puede igualarse.
///
/// El centinela sirve a los dos usos, cuya semántica difiere:
///
/// - `is_primitive_motif` compara dos bytes no ACGT como **iguales** (`0 == 0`).
/// - la puntuación de [`find_tandem_repeats`] exige `h != 0` además de la
///   igualdad, de modo que un byte no ACGT **no coincide con nada**, ni
///   siquiera con otro igual: una `N` junto a una `N` es ausencia de pruebas,
///   no prueba de una repetición.
static CANONICAL: [u8; 256] = {
    let mut table = [0u8; 256];
    table[b'A' as usize] = b'A';
    table[b'a' as usize] = b'A';
    table[b'C' as usize] = b'C';
    table[b'c' as usize] = b'C';
    table[b'G' as usize] = b'G';
    table[b'g' as usize] = b'G';
    table[b'T' as usize] = b'T';
    table[b't' as usize] = b'T';
    table
};

/// Un segmento abierto en la pila de trabajo de Ruzzo y Tompa: `l`/`r` son la
/// puntuación acumulada estrictamente antes del segmento y hasta su final (la
/// puntuación propia es `r - l`); `start`/`end` son índices inclusivos en la
/// secuencia de puntuaciones.
struct RtSeg {
    l: i64,
    r: i64,
    start: usize,
    end: usize,
    /// **Enlace al anterior menor**: `1 +` el índice del segmento más cercano
    /// *por debajo* de este cuyo `l` es menor, o `0` si no hay ninguno. El
    /// sesgo de `+ 1` permite representar «ninguno» sin `Option` y hace que la
    /// comprobación de límite y el salto del bucle de búsqueda sean la misma
    /// comparación.
    ///
    /// Es lo que permite a la búsqueda saltar en vez de recorrer. Mantenerlo no
    /// cuesta nada: el enlace que un segmento necesita es exactamente el `j`
    /// que ya encontró la búsqueda que lo colocó.
    prev_smaller: usize,
}

/// Encuentra todas las **subsecuencias de puntuación máxima** de una secuencia
/// de puntuaciones (Ruzzo y Tompa, 1999) y llama a `emit(start, end, score)`
/// una vez por segmento —`start`/`end` inclusivos, `score = r - l > 0`—, en
/// orden creciente de `start`.
///
/// Una subsecuencia de puntuación máxima es una tirada contigua de total
/// positivo que no puede alargarse ni recortarse sin bajar su puntuación; las
/// sustituciones y las interrupciones cortas quedan dentro de un segmento
/// mientras la puntuación de alrededor las compense.
///
/// Es el pase fuera de línea de manual: una pila de segmentos abiertos que se
/// vacía al final. Su corrección la fija una prueba de propiedad contra fuerza
/// bruta (`ruzzo_tompa_matches_brute_force_on_pseudo_random`). La pila puede
/// crecer hasta O(n) sobre una entrada patológica de una secuencia entera.
fn maximal_scoring_subsequences(
    scores: impl Iterator<Item = i64>,
    mut emit: impl FnMut(usize, usize, i64),
) {
    let mut stack: Vec<RtSeg> = Vec::new();
    let mut total: i64 = 0;
    for (i, s) in scores.enumerate() {
        let l = total;
        total += s;
        let r = total;
        if s <= 0 {
            continue; // una puntuación no positiva no abre segmento; solo avanza el total
        }
        let mut cur = RtSeg {
            l,
            r,
            start: i,
            end: i,
            prev_smaller: 0,
        };
        loop {
            // El segmento apilado más a la derecha cuyo `l` es menor que
            // `cur.l`, hallado **saltando por la cadena de anteriores menores**
            // en vez de recorrer.
            //
            // Corrección: la cadena se salta exactamente los segmentos que no
            // pueden ser la respuesta. Todo lo que hay entre `t` y
            // `stack[t].prev_smaller` tiene `l >= l[t]`, por definición del
            // enlace; así que una vez `l[t] >= cur.l`, ninguno está por debajo
            // de `cur.l`, y el primer `t` con `l[t] < cur.l` es el de más a la
            // derecha. La misma respuesta que el recorrido al que sustituye.
            //
            // **No es una búsqueda binaria** (se probó y se revirtió en el
            // original, 2026-07-20): los `l` de la pila **no** están ordenados,
            // y la regla 2 es lo que los desordena —cuando no hay ninguno
            // menor, `cur` se apila igualmente y deja un `l` menor encima de
            // otros mayores—. `partition_point` ahí compila, pasa un
            // `debug_assert` que en optimizado desaparece y encuentra en
            // silencio **cero** loci sobre secuencia real.
            let mut idx = stack.len();
            let mut found = None;
            while idx > 0 {
                let t = idx - 1;
                if stack[t].l < cur.l {
                    found = Some(t);
                    break;
                }
                idx = stack[t].prev_smaller;
            }
            match found {
                None => {
                    cur.prev_smaller = 0;
                    break; // Regla 2: sin segmento menor a la izquierda → `cur` va aparte.
                }
                Some(j) => {
                    if stack[j].r >= cur.r {
                        // Regla 3: `cur` sigue aparte a la derecha, y `j`, la
                        // respuesta que acaba de dar la búsqueda, es justo el
                        // enlace que `cur` debe llevar.
                        cur.prev_smaller = j + 1;
                        break;
                    }
                    // Regla 4: `cur` absorbe `stack[j..]`: amplía su borde
                    // izquierdo y vuelve a buscar. Truncar no invalida ningún
                    // enlace superviviente: todos apuntan estrictamente por
                    // debajo de su propio segmento.
                    cur.l = stack[j].l;
                    cur.start = stack[j].start;
                    stack.truncate(j);
                }
            }
        }
        stack.push(cur);
    }
    for seg in stack.drain(..) {
        emit(seg.start, seg.end, seg.r - seg.l);
    }
}

/// Encuentra todo intervalo de repetición en tándem de `seq` cuyo periodo esté
/// en `periods`.
///
/// Para cada periodo `p` puntúa la posición `j` (`j >= p`) con `+match_reward`
/// cuando `seq[j]` y `seq[j-p]` son la misma base ACGT y con
/// `-mismatch_penalty` en otro caso (una base no ACGT nunca coincide), toma
/// cada segmento de puntuación máxima de Ruzzo y Tompa, lo traduce al tramo
/// semiabierto en base cero que certifica y lo emite si sus copias implícitas
/// alcanzan `params.min_copies`. Pura y total sobre bytes arbitrarios;
/// determinista.
///
/// **Cada repetición se informa una vez, en su periodo verdadero**, con dos
/// mecanismos de pura geometría:
/// - un **motivo no primitivo** se descarta al emitir (`AA` se reduce a `A`), de
///   modo que un homopolímero es de periodo 1 y nunca un alias de periodo 2;
/// - una **re-detección del mismo tracto con un periodo múltiplo** se descarta
///   en un pase posterior: un `(GAA)n` de periodo 3 también encaja, mal, como
///   un periodo 6 `AAGATG`, que es un 6-mero primitivo, pero el de periodo 3
///   **cubre la mayor parte** de él: es el mismo tracto visto peor.
///
/// La prueba de «cubre la mayor parte» es lo que la mantiene fuera de la
/// política: una repetición corta que solo *está dentro* de otra larga (la
/// `AAA` dentro de un `(GAAA)n`) cubre una astilla y **no** se toma por
/// re-detección; se conservan las dos.
#[must_use]
pub fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> Vec<RepeatInterval> {
    let mut out = Vec::new();
    let n = seq.len();
    // Se canonicaliza UNA vez, no dos por periodo.
    let mut canonical: Vec<u8> = Vec::with_capacity(n);
    canonical.extend(seq.iter().map(|&b| CANONICAL[b as usize]));
    for period in periods.min()..=periods.max() {
        let p = period as usize;
        if p >= n {
            continue; // ninguna posición tiene pareja p bases atrás
        }
        // El índice de puntuación `k` (base cero sobre `p..n`) corresponde a la
        // posición `j = k + p`.
        let reward = i64::from(params.match_reward);
        let penalty = -i64::from(params.mismatch_penalty);
        let back = &canonical[..n - p];
        let here = &canonical[p..];
        let scores = back
            .iter()
            .zip(here.iter())
            .map(|(&b, &h)| if h != 0 && h == b { reward } else { penalty });
        maximal_scoring_subsequences(scores, |k0, k1, score| {
            // Segmento [k0, k1] → tramo [k0, k1 + p + 1): la primera base
            // implicada es `j0 - p = k0` y la última `j1 = k1 + p`, así que el
            // fin exclusivo es `k1 + p + 1`.
            let tract_start = k0 as u64;
            let tract_end = (k1 + p + 1) as u64;
            let copies = (tract_end - tract_start) / u64::from(period);
            // Solo repeticiones de periodo PRIMITIVO: `AAAAAA` se detecta en
            // periodo 2, pero su motivo `AA` se reduce a `A`; es un
            // homopolímero, no un dinucleótido, e informarlo como periodo 2
            // sería incumplir el contrato del detector.
            let motif = &seq[tract_start as usize..tract_start as usize + p];
            if copies >= u64::from(params.min_copies) && is_primitive_motif(motif) {
                out.push(RepeatInterval {
                    start: tract_start,
                    end: tract_end,
                    period,
                    score: i32::try_from(score).unwrap_or(i32::MAX),
                });
            }
        });
    }

    // Descarta la re-detección con periodo múltiplo del mismo tracto. `iv` (el
    // periodo mayor) sobra si un periodo menor `a` lo divide y **cubre al menos
    // la mitad** de él: «cubre la mayor parte» es lo que dice *mismo tracto* y
    // no solo *solapado*. Ninguna política de copias aquí: solo cuánto de `iv`
    // explica el menor.
    out.sort_by_key(|iv| (iv.period, iv.start));

    // Supervivientes agrupados por periodo, que es lo que convierte la prueba
    // de redundancia en una búsqueda binaria. Tres hechos la hacen exacta:
    //
    // 1. **Solo un periodo menor que lo divida puede eliminar a `iv`**: a lo
    //    sumo `{1, 2, 3}`, para el periodo 6.
    // 2. **Esas listas ya están completas** al llegar a `iv`: `out` está
    //    ordenado primero por periodo.
    // 3. **Dentro de un periodo los tramos están ordenados por inicio Y por
    //    fin**: los segmentos de Ruzzo y Tompa son disjuntos, y cada tramo
    //    alarga su segmento en el mismo `p + 1`. Es la precondición de
    //    `partition_point`.
    let max_period = out.iter().map(|iv| iv.period).max().unwrap_or(0) as usize;
    let mut by_period: Vec<Vec<RepeatInterval>> = vec![Vec::new(); max_period + 1];
    let mut kept_len = 0usize;
    for iv in out {
        let iv_len = iv.end - iv.start;
        let mut redundant = false;
        for shorter in 1..iv.period {
            if !iv.period.is_multiple_of(shorter) {
                continue;
            }
            let candidates = &by_period[shorter as usize];
            // Los tramos que solapan `iv` son una tirada contigua: se saltan los
            // que acaban antes de que empiece y se para en el primero que
            // empieza donde acaba o después.
            let first = candidates.partition_point(|a| a.end <= iv.start);
            for a in &candidates[first..] {
                if a.start >= iv.end {
                    break;
                }
                // `a` debe cubrir la mitad de `iv`; uno más corto que eso no
                // puede, se coloque como se coloque.
                if (a.end - a.start) * 2 < iv_len {
                    continue;
                }
                let overlap = a.end.min(iv.end).saturating_sub(a.start.max(iv.start));
                if overlap * 2 >= iv_len {
                    redundant = true;
                    break;
                }
            }
            if redundant {
                break;
            }
        }
        if !redundant {
            kept_len += 1;
            by_period[iv.period as usize].push(iv);
        }
    }
    let mut kept: Vec<RepeatInterval> = Vec::with_capacity(kept_len);
    for group in by_period {
        kept.extend(group);
    }
    kept
}

// ---------------------------------------------------------------------
// Conversión de tramos en candidatos (propia; provisional)
// ---------------------------------------------------------------------

/// Longitudes de motivo que ve la detección: de 1 a 6 pb (REQ-REF-38).
pub const DETECTION_MOTIF_LENGTHS: usize = 6;

/// Mínimo de copias completas por longitud de motivo, de 1 a 6 pb
/// (REQ-REF-39, decisión 5 del responsable).
pub const DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH: [u64; DETECTION_MOTIF_LENGTHS] = [8, 6, 6, 6, 5, 4];

/// Pureza mínima por omisión: la general de la admisión de loci
/// (`[reference.repeat_admission]`), sin una segunda perilla (DP-SSR-02). Se
/// toma de la propia admisión para que no puedan separarse.
pub const DEFAULT_DETECTION_MIN_PURITY: f64 = EffectiveRepeatAdmission::DEFAULT.min_purity();

/// Tirada mínima por omisión, en copias contiguas: la general de la admisión.
pub const DEFAULT_DETECTION_MIN_LONGEST_RUN_UNITS: u64 =
    EffectiveRepeatAdmission::DEFAULT.min_longest_run_units();

/// Tope de unidades por omisión para un detectado: el de todo locus
/// (`DEFAULT_REPEAT_MAX_UNITS`); por encima es un satélite (DP-SSR-07).
pub const DEFAULT_DETECTION_MAX_UNITS: u64 = DEFAULT_REPEAT_MAX_UNITS;

/// Criterios **provisionales** que convierten un tramo detectado en candidato
/// (DP-SSR-02 a DP-SSR-08). Son parámetros para poder medirlos; la pausa de la
/// entrega 1 decide cuáles se exponen.
#[derive(Clone, Debug, PartialEq)]
pub struct DetectionCriteria {
    /// Mínimo de **copias completas** de la descomposición voraz del simulador,
    /// por longitud de motivo (índice 0 = 1 pb). DP-SSR-02, opción (b).
    pub min_copies: [u64; DETECTION_MOTIF_LENGTHS],
    /// Pureza mínima sobre bases, comparada igual que la admisión de loci.
    pub min_purity: f64,
    /// Tirada más larga mínima, en copias contiguas.
    pub min_longest_run_units: u64,
    /// Tope de unidades: un tramo con más es un satélite y se descarta.
    pub max_units: u64,
}

impl Default for DetectionCriteria {
    fn default() -> Self {
        Self {
            min_copies: DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH,
            min_purity: DEFAULT_DETECTION_MIN_PURITY,
            min_longest_run_units: DEFAULT_DETECTION_MIN_LONGEST_RUN_UNITS,
            max_units: DEFAULT_DETECTION_MAX_UNITS,
        }
    }
}

/// Etapa en la que se descarta un tramo, **en el orden en que se aplican**.
/// Un tramo se cuenta solo en la primera etapa que no supera.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub enum DetectionStage {
    /// El motivo leído al inicio del tramo tiene una base que no es A, C, G o
    /// T en mayúscula, y el simulador no puede representarlo. El detector lo
    /// admite en un periodo de 5 o 6 si la discrepancia cae al final del
    /// primer periodo.
    NonCanonicalMotif,
    /// Menos copias completas que el mínimo de su longitud de motivo.
    MinCopies,
    /// Pureza por debajo del mínimo.
    MinPurity,
    /// Tirada más larga por debajo del mínimo.
    MinLongestRun,
    /// Más unidades que el tope: satélite.
    MaxUnits,
    /// Sin una base A/C/G/T a cada lado: borde de secuencia o `N` contigua
    /// (DP-SSR-08, opción (a)).
    Edge,
    /// Solapa con otro candidato que gana por la regla de DP-SSR-05.
    Overlap,
}

impl DetectionStage {
    /// Todas las etapas, en orden de aplicación.
    pub const ALL: [Self; 7] = [
        Self::NonCanonicalMotif,
        Self::MinCopies,
        Self::MinPurity,
        Self::MinLongestRun,
        Self::MaxUnits,
        Self::Edge,
        Self::Overlap,
    ];

    /// Nombre estable en inglés, para tablas de medición.
    #[must_use]
    pub const fn name(self) -> &'static str {
        match self {
            Self::NonCanonicalMotif => "non_canonical_motif",
            Self::MinCopies => "min_copies",
            Self::MinPurity => "min_purity",
            Self::MinLongestRun => "min_longest_run",
            Self::MaxUnits => "max_units",
            Self::Edge => "edge",
            Self::Overlap => "overlap",
        }
    }

    const fn index(self) -> usize {
        self as usize
    }
}

/// Recuento por etapa y por longitud de motivo.
///
/// Cumple siempre `detected = descartados en todas las etapas + accepted`
/// para cada longitud de motivo; [`DetectionCounts::is_balanced`] lo comprueba.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct DetectionCounts {
    detected: [u64; DETECTION_MOTIF_LENGTHS],
    discarded: [[u64; DETECTION_MOTIF_LENGTHS]; DetectionStage::ALL.len()],
    accepted: [u64; DETECTION_MOTIF_LENGTHS],
}

impl DetectionCounts {
    /// Intervalos que entregó el detector con esa longitud de motivo (1 a 6).
    #[must_use]
    pub fn detected(&self, motif_length: usize) -> u64 {
        self.detected[motif_length - 1]
    }

    /// Descartados en esa etapa con esa longitud de motivo.
    #[must_use]
    pub fn discarded(&self, stage: DetectionStage, motif_length: usize) -> u64 {
        self.discarded[stage.index()][motif_length - 1]
    }

    /// Candidatos aceptados con esa longitud de motivo.
    #[must_use]
    pub fn accepted(&self, motif_length: usize) -> u64 {
        self.accepted[motif_length - 1]
    }

    /// Total de intervalos del detector.
    #[must_use]
    pub fn detected_total(&self) -> u64 {
        self.detected.iter().sum()
    }

    /// Total de descartados en una etapa.
    #[must_use]
    pub fn discarded_total(&self, stage: DetectionStage) -> u64 {
        self.discarded[stage.index()].iter().sum()
    }

    /// Total de candidatos aceptados.
    #[must_use]
    pub fn accepted_total(&self) -> u64 {
        self.accepted.iter().sum()
    }

    /// Si cada longitud de motivo cuadra: detectados = descartados + aceptados.
    #[must_use]
    pub fn is_balanced(&self) -> bool {
        (0..DETECTION_MOTIF_LENGTHS).all(|index| {
            let discarded: u64 = self.discarded.iter().map(|stage| stage[index]).sum();
            self.detected[index] == discarded + self.accepted[index]
        })
    }

    /// Suma otro recuento a este (para acumular secuencias).
    pub fn add(&mut self, other: &Self) {
        for index in 0..DETECTION_MOTIF_LENGTHS {
            self.detected[index] += other.detected[index];
            self.accepted[index] += other.accepted[index];
            for stage in 0..DetectionStage::ALL.len() {
                self.discarded[stage][index] += other.discarded[stage][index];
            }
        }
    }

    fn discard(&mut self, stage: DetectionStage, motif_length: usize) {
        self.discarded[stage.index()][motif_length - 1] += 1;
    }
}

/// Un tramo que supera los criterios por sí solo. Coordenadas en base cero y
/// semiabiertas sobre la secuencia barrida; el motivo es el de la fase del
/// inicio del tramo (DP-SSR-03, opción (a)) y las fronteras son las del
/// detector, sin recortar (DP-SSR-04, opción (a)).
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RepeatCandidate {
    start: u64,
    end: u64,
    motif: RepeatMotif,
    score: i32,
    covered_bases: u64,
    longest_run_units: u64,
}

impl RepeatCandidate {
    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn end(&self) -> u64 {
        self.end
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.end - self.start
    }

    #[must_use]
    pub fn motif(&self) -> &RepeatMotif {
        &self.motif
    }

    /// Puntuación de Ruzzo y Tompa del detector.
    #[must_use]
    pub const fn score(&self) -> i32 {
        self.score
    }

    /// Copias completas de la descomposición voraz (bases cubiertas entre la
    /// longitud del motivo).
    #[must_use]
    pub fn complete_copies(&self) -> u64 {
        self.covered_bases / self.motif.length()
    }

    /// Bases cubiertas por copias completas.
    #[must_use]
    pub const fn covered_bases(&self) -> u64 {
        self.covered_bases
    }

    #[must_use]
    pub const fn longest_run_units(&self) -> u64 {
        self.longest_run_units
    }

    /// Clave de prioridad de DP-SSR-05: mayor puntuación, más bases, periodo
    /// menor, más a la izquierda. Menor clave = gana.
    fn priority(&self) -> (std::cmp::Reverse<i32>, std::cmp::Reverse<u64>, u64, u64) {
        (
            std::cmp::Reverse(self.score),
            std::cmp::Reverse(self.length()),
            self.motif.length(),
            self.start,
        )
    }
}

/// Resultado de barrer una secuencia: los candidatos, sin solapes y ordenados
/// por inicio, y el recuento por etapa.
#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct SequenceDetection {
    candidates: Vec<RepeatCandidate>,
    counts: DetectionCounts,
}

impl SequenceDetection {
    #[must_use]
    pub fn candidates(&self) -> &[RepeatCandidate] {
        &self.candidates
    }

    #[must_use]
    pub const fn counts(&self) -> &DetectionCounts {
        &self.counts
    }
}

/// Barre una secuencia y devuelve sus candidatos a microsatélite con el
/// recuento de descartes por etapa y longitud de motivo.
///
/// Usa el detector con su puntuación fija (2 y 7) y periodos de 1 a 6, aplica
/// a cada intervalo [`classify_interval`] y resuelve después los solapes entre
/// los supervivientes con la regla de DP-SSR-05. Pura y determinista.
#[must_use]
pub fn detect_repeat_candidates(bases: &[u8], criteria: &DetectionCriteria) -> SequenceDetection {
    let periods = PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1)
        .unwrap_or_else(|_| unreachable!("1..=6 es un intervalo válido"));
    let intervals = find_tandem_repeats(bases, periods, &ScanParams::default());
    let mut counts = DetectionCounts::default();
    let mut survivors = Vec::new();
    for interval in &intervals {
        let motif_length = usize::from(interval.period);
        counts.detected[motif_length - 1] += 1;
        match classify_interval(bases, interval, criteria) {
            Ok(candidate) => survivors.push(candidate),
            Err(stage) => counts.discard(stage, motif_length),
        }
    }
    let (candidates, overlapped) = resolve_overlaps(survivors);
    for loser in &overlapped {
        counts.discard(DetectionStage::Overlap, loser.motif.length() as usize);
    }
    for candidate in &candidates {
        counts.accepted[candidate.motif.length() as usize - 1] += 1;
    }
    SequenceDetection { candidates, counts }
}

/// Juzga un intervalo por sí solo, sin mirar a los demás, y devuelve el
/// candidato o la **primera** etapa que no supera.
///
/// Las etapas, en orden: motivo representable, copias completas, pureza,
/// tirada, tope de unidades y bordes. La pureza y la tirada se miden con
/// [`observe_tract`], la misma descomposición voraz con la que el simulador
/// admite un locus declarado, con el motivo en la fase del inicio del tramo.
///
/// # Panics
///
/// Si el intervalo no cabe en `bases` o tiene periodo 0 o mayor que 6: el
/// detector nunca los produce y una prueba que los fabrique es un error suyo.
pub fn classify_interval(
    bases: &[u8],
    interval: &RepeatInterval,
    criteria: &DetectionCriteria,
) -> Result<RepeatCandidate, DetectionStage> {
    let period = usize::from(interval.period);
    assert!(
        (1..=DETECTION_MOTIF_LENGTHS).contains(&period),
        "periodo {period} fuera de 1..=6"
    );
    let start = interval.start as usize;
    let end = interval.end as usize;
    let tract = &bases[start..end];
    let motif =
        RepeatMotif::new(&tract[..period]).map_err(|_| DetectionStage::NonCanonicalMotif)?;
    let observed = observe_tract(tract, std::slice::from_ref(&motif));
    let complete_copies = observed.covered_bases() / period as u64;
    if complete_copies < criteria.min_copies[period - 1] {
        return Err(DetectionStage::MinCopies);
    }
    if observed.purity_as_f64() < criteria.min_purity {
        return Err(DetectionStage::MinPurity);
    }
    if observed.longest_run_units() < criteria.min_longest_run_units {
        return Err(DetectionStage::MinLongestRun);
    }
    // Unidades = bases / longitud del motivo, la medida con la que el
    // generador de deslizamiento compara su tope (`units_as_f64` de un solo
    // bloque). `bases / p > tope` equivale exactamente a `bases > tope · p`.
    if tract.len() as u64 > criteria.max_units.saturating_mul(period as u64) {
        return Err(DetectionStage::MaxUnits);
    }
    let flanked = start > 0
        && end < bases.len()
        && CANONICAL[bases[start - 1] as usize] != 0
        && CANONICAL[bases[end] as usize] != 0;
    if !flanked {
        return Err(DetectionStage::Edge);
    }
    Ok(RepeatCandidate {
        start: interval.start,
        end: interval.end,
        motif,
        score: interval.score,
        covered_bases: observed.covered_bases(),
        longest_run_units: observed.longest_run_units(),
    })
}

/// Resuelve los solapes con la regla de DP-SSR-05, opción (a): ordena los
/// candidatos por prioridad —mayor puntuación del detector; a igual
/// puntuación, más bases; luego periodo menor; luego el de la izquierda— y
/// los recorre en ese orden aceptando cada uno que **no solape a uno ya
/// aceptado**.
///
/// No es «de cada grupo solapado gana uno»: en una cadena A–B–C donde A gana
/// a B y C solo solapa a B, B se descarta y **quedan A y C**, porque C ya no
/// solapa a nadie aceptado. El resultado **no depende del orden de
/// entrada**: la clave ordena del todo (dos candidatos con la misma clave
/// tendrían el mismo intervalo y periodo, y el detector no los emite).
/// Dos tramos que solo se tocan (uno acaba donde empieza el otro) no se
/// solapan. Devuelve (aceptados ordenados por inicio, descartados).
fn resolve_overlaps(
    mut survivors: Vec<RepeatCandidate>,
) -> (Vec<RepeatCandidate>, Vec<RepeatCandidate>) {
    survivors.sort_by_key(RepeatCandidate::priority);
    // Aceptados por inicio. Son disjuntos, así que también están ordenados por
    // fin y basta mirar el último que empieza antes del fin del aspirante.
    let mut accepted: BTreeMap<u64, RepeatCandidate> = BTreeMap::new();
    let mut discarded = Vec::new();
    for candidate in survivors {
        let overlaps = accepted
            .range(..candidate.end)
            .next_back()
            .is_some_and(|(_, previous)| previous.end > candidate.start);
        if overlaps {
            discarded.push(candidate);
        } else {
            accepted.insert(candidate.start, candidate);
        }
    }
    (accepted.into_values().collect(), discarded)
}

#[cfg(test)]
mod tests {
    use super::*;

    // ---- PeriodRange (copiadas) --------------------------------------------

    #[test]
    fn period_range_accepts_valid_ranges() {
        assert_eq!(
            PeriodRange::new(1, 6).map(|r| (r.min(), r.max())),
            Ok((1, 6))
        );
        // Un intervalo de un solo periodo es válido.
        assert_eq!(
            PeriodRange::new(6, 6).map(|r| (r.min(), r.max())),
            Ok((6, 6))
        );
        assert_eq!(
            PeriodRange::new(2, 6).map(|r| (r.min(), r.max())),
            Ok((2, 6))
        );
    }

    #[test]
    fn period_range_rejects_zero_min() {
        assert_eq!(PeriodRange::new(0, 6), Err(PeriodRangeError::ZeroMin));
        // El mínimo cero se comprueba antes que min > max: (0, 0) es ZeroMin.
        assert_eq!(PeriodRange::new(0, 0), Err(PeriodRangeError::ZeroMin));
    }

    #[test]
    fn period_range_rejects_min_above_max() {
        assert_eq!(
            PeriodRange::new(3, 2),
            Err(PeriodRangeError::MinExceedsMax { min: 3, max: 2 })
        );
    }

    /// Propia: el `Display` escrito a mano conserva los textos de `thiserror`.
    #[test]
    fn period_range_error_display_keeps_the_original_texts() {
        assert_eq!(
            PeriodRangeError::ZeroMin.to_string(),
            "period range min must be >= 1 (period 0 is meaningless)"
        );
        assert_eq!(
            PeriodRangeError::MinExceedsMax { min: 3, max: 2 }.to_string(),
            "period range min (3) exceeds max (2)"
        );
    }

    // ---- Valores por omisión (copiadas) ------------------------------------

    #[test]
    fn default_periods_is_one_to_six() {
        assert_eq!(DEFAULT_PERIODS, (1, 6));
        assert!(PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1).is_ok());
    }

    #[test]
    fn scan_params_default_is_two_seven_two() {
        let p = ScanParams::default();
        assert_eq!(p.match_reward, 2);
        assert_eq!(p.mismatch_penalty, 7);
        assert_eq!(p.min_copies, 2);
    }

    // ---- Ruzzo–Tompa contra un oráculo de fuerza bruta (copiadas) ----------

    /// Recoge `maximal_scoring_subsequences` en un vector `(start, end, score)`
    /// ordenado por inicio.
    fn rt_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let mut out = Vec::new();
        maximal_scoring_subsequences(scores.iter().copied(), |a, b, s| out.push((a, b, s)));
        out.sort_unstable();
        out
    }

    /// Todas las subsecuencias de puntuación máxima por fuerza bruta, desde la
    /// definición de Ruzzo y Tompa (O(n⁴), solo para entradas pequeñas: es el
    /// oráculo independiente). Un intervalo inclusivo `[i, j]` vale si (1) su
    /// puntuación es positiva; (2) todo subintervalo **propio** puntúa
    /// estrictamente menos; y (3) no está contenido propiamente en otro que
    /// cumpla (1) y (2).
    fn brute_segments(scores: &[i64]) -> Vec<(usize, usize, i64)> {
        let n = scores.len();
        let mut pre = vec![0i64; n + 1];
        for i in 0..n {
            pre[i + 1] = pre[i] + scores[i];
        }
        let score = |i: usize, j: usize| pre[j + 1] - pre[i]; // [i, j] inclusivo

        let mut cands: Vec<(usize, usize, i64)> = Vec::new();
        for i in 0..n {
            for j in i..n {
                let sc = score(i, j);
                if sc <= 0 {
                    continue;
                }
                let internally_optimal =
                    (i..=j).all(|i2| (i2..=j).all(|j2| (i2, j2) == (i, j) || score(i2, j2) < sc));
                if internally_optimal {
                    cands.push((i, j, sc));
                }
            }
        }

        let mut res: Vec<(usize, usize, i64)> = cands
            .iter()
            .copied()
            .filter(|&(i, j, _)| {
                !cands
                    .iter()
                    .any(|&(a, b, _)| a <= i && j <= b && (a, b) != (i, j))
            })
            .collect();
        res.sort_unstable();
        res
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_crafted_cases() {
        let cases: &[&[i64]] = &[
            &[],
            &[-7, -7, -7],
            &[2, 2, 2, 2],
            &[2, -7, 2],  // dos sueltos: la interrupción cuesta demasiado
            &[2, -3, 2],  // un solo segmento: el bache se salva
            &[5, -10, 8], // el bache baja del L izquierdo → dos segmentos
            &[5, -3, 8],  // el bache se queda por encima → uno
            &[2, 2, -7, 2, 2],
            &[-7, 2, 2, -7, 2, -7, 2, 2, -7],
            &[1, 1, 1, -2, 1, 1, 1, -10, 5, 5],
        ];
        for case in cases {
            assert_eq!(rt_segments(case), brute_segments(case), "case {case:?}");
        }
    }

    #[test]
    fn ruzzo_tompa_matches_brute_force_on_pseudo_random() {
        // LCG determinista sobre puntuaciones de tipo ±; sin dependencia de azar.
        let mut state: u64 = 0x9E3779B97F4A7C15;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        for _ in 0..2000 {
            let len = (next() % 160) as usize;
            let scores: Vec<i64> = (0..len)
                .map(|_| if next() % 4 == 0 { 2 } else { -7 }) // 25 % de coincidencias
                .collect();
            assert_eq!(
                rt_segments(&scores),
                brute_segments(&scores),
                "scores {scores:?}"
            );
        }
    }

    // ---- find_tandem_repeats (copiadas) ------------------------------------

    fn p(min: u8, max: u8) -> PeriodRange {
        PeriodRange::new(min, max).unwrap()
    }

    /// Un (CAG)×8 limpio sale en periodo 3 como un único intervalo exacto.
    #[test]
    fn finds_a_clean_perfect_tract() {
        let seq = b"CAGCAGCAGCAGCAGCAGCAGCAG"; // 8 copias, 24 pb
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(
            got,
            vec![RepeatInterval {
                start: 0,
                end: 24,
                period: 3,
                score: 2 * (24 - 3), // recompensa × coincidencias
            }]
        );
    }

    /// Una sustitución dentro de un tracto lo deja en **un** intervalo.
    #[test]
    fn a_substitution_keeps_one_interval() {
        let mut seq = Vec::new();
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        seq.extend_from_slice(b"CAT"); // la copia interrumpida (G→T)
        for _ in 0..5 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "one span across the substitution: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u64);
    }

    /// Un indel aislado (una ráfaga acotada de discrepancias) también deja el
    /// tracto en **un** intervalo: la comparación es local y las copias de
    /// después recuperan la fase.
    #[test]
    fn a_single_indel_keeps_one_interval() {
        let mut seq = Vec::new();
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        seq.push(b'T'); // inserción de 1 pb
        for _ in 0..6 {
            seq.extend_from_slice(b"CAG");
        }
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        let period3: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(
            period3.len(),
            1,
            "the single indel is a bounded burst, not a split: {got:?}"
        );
        assert_eq!(period3[0].start, 0);
        assert_eq!(period3[0].end, seq.len() as u64);
    }

    /// Una tirada de `N` no produce nada: `N == N` no es coincidencia.
    #[test]
    fn an_n_run_emits_nothing() {
        let seq = b"NNNNNNNNNNNNNNNNNNNN";
        assert!(find_tandem_repeats(seq, p(1, 6), &ScanParams::default()).is_empty());
    }

    /// Un tracto en minúsculas se encuentra igual: la comparación no distingue
    /// mayúsculas.
    #[test]
    fn a_soft_masked_tract_is_found() {
        let seq = b"cagcagcagcagcagcagcagcag";
        let got = find_tandem_repeats(seq, p(3, 3), &ScanParams::default());
        assert_eq!(got.len(), 1);
        assert_eq!((got[0].start, got[0].end, got[0].period), (0, 24, 3));
    }

    /// Un (AT)×10 coincide en periodo 2 y en sus múltiplos 4 y 6, pero su
    /// periodo primitivo es 2 y solo se emite ahí.
    #[test]
    fn only_the_primitive_period_is_emitted() {
        let seq = b"ATATATATATATATATATAT"; // (AT)×10, 20 pb
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        let periods: std::collections::BTreeSet<u8> = got.iter().map(|r| r.period).collect();
        assert_eq!(
            periods,
            std::collections::BTreeSet::from([2]),
            "only the primitive period 2, no 4/6 aliases (nor 3/5 non-divisors): {got:?}"
        );
    }

    /// Un periodo 4 primitivo de verdad (`ACGT`) sale en 4 y no se confunde con
    /// un alias.
    #[test]
    fn a_genuine_higher_primitive_period_is_kept() {
        let seq = b"ACGTACGTACGTACGT"; // (ACGT)×4
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 4),
            "ACGT is a primitive period-4 motif and must survive: {got:?}"
        );
        assert!(
            got.iter().all(|r| r.period != 2 && r.period != 3),
            "it does not tile under any shorter period, so no 2/3 detection: {got:?}"
        );
    }

    /// Un homopolímero es un tracto de periodo 1 y nada más.
    #[test]
    fn a_homopolymer_is_only_period_one() {
        let seq = b"AAAAAAAAAAAAAAAAAAAA"; // poli-A, 20 pb
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        let periods: std::collections::BTreeSet<u8> = got.iter().map(|r| r.period).collect();
        assert_eq!(
            periods,
            std::collections::BTreeSet::from([1]),
            "a poly-A is period 1 only — the whole aliasing problem, gone at the source: {got:?}"
        );
        // Y pedido solo de 2 en adelante, no devuelve nada: es el contrato.
        let got2 = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got2.is_empty(),
            "asked for period 2..6, a homopolymer is not any of them: {got2:?}"
        );
    }

    /// **Una re-detección con periodo múltiplo y fase desplazada se descarta;
    /// una repetición que solo está dentro de otra se conserva.**
    ///
    /// Nota propia (riesgo R-2 del plan): aquí el intervalo de periodo 6 lo
    /// elimina el **filtro de motivo primitivo**, no la geometría: su motivo
    /// sale no primitivo y se descarta antes del pase geométrico. Anular el
    /// descarte geométrico deja esta prueba en verde. La geometría la vigila
    /// `a_primitive_period_six_re_detection_is_dropped_by_geometry`.
    #[test]
    fn a_period_multiple_re_detection_is_dropped_by_geometry() {
        let seq = b"ATGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAG";
        let got = find_tandem_repeats(seq, p(2, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 3),
            "the primitive period-3 GAA survives: {got:?}"
        );
        assert!(
            !got.iter().any(|r| r.period == 6),
            "the period-6 re-detection of the same tract is dropped: {got:?}"
        );

        // Un homopolímero DENTRO de un motivo de periodo mayor es otra
        // repetición, más pequeña, que cubre una astilla: no debe eliminar a su
        // padre.
        let seq = b"GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA"; // (GAAA)×8
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        assert!(
            got.iter().any(|r| r.period == 4),
            "the period-4 (GAAA) survives its own internal AAA homopolymer runs: {got:?}"
        );
    }

    /// **La geometría, sola.** `GAA GTA (GAA)×10`: la sustitución de la segunda
    /// copia hace que el segmento de periodo 6 empiece en la base 0 con el
    /// motivo `GAAGTA`, que **es primitivo** (lo comprueba la prueba), así que
    /// el filtro de motivos no lo toca. El de periodo 3 empieza después de la
    /// sustitución y cubre más de la mitad del de periodo 6: solo el descarte
    /// geométrico puede quitarlo.
    #[test]
    fn a_primitive_period_six_re_detection_is_dropped_by_geometry() {
        let mut seq = b"GAAGTA".to_vec();
        for _ in 0..10 {
            seq.extend_from_slice(b"GAA");
        }
        assert!(
            is_primitive_motif(&seq[..6]),
            "el motivo de periodo 6 es primitivo"
        );
        let got = find_tandem_repeats(&seq, p(3, 6), &ScanParams::default());
        let three: Vec<_> = got.iter().filter(|r| r.period == 3).collect();
        assert_eq!(three.len(), 1, "{got:?}");
        assert!(
            three[0].start > 0 && 2 * (three[0].end - three[0].start) >= seq.len() as u64,
            "el de periodo 3 empieza tras la sustitución y cubre la mitad: {got:?}"
        );
        assert!(
            !got.iter().any(|r| r.period == 6),
            "la re-detección primitiva de periodo 6 se descarta por geometría: {got:?}"
        );
    }

    /// El suelo `min_copies` descarta segmentos cortos (p. ej. una coincidencia
    /// casual suelta).
    #[test]
    fn the_min_copies_floor_drops_short_segments() {
        let seq = b"ACGTACAGTCTGAC";
        let got = find_tandem_repeats(seq, p(2, 2), &ScanParams::default());
        assert!(
            got.iter().all(|r| (r.end - r.start) / 2 >= 2),
            "every emitted interval has >= min_copies copies: {got:?}"
        );
    }

    // ---- Conversión en candidatos (propias) --------------------------------

    /// Generador determinista de bases **sin repeticiones**: cada base difiere
    /// de las dos anteriores, así que no hay coincidencias a distancia 1 ni 2,
    /// y las de distancia 3 a 6 son sueltas. La ausencia de candidatos en el
    /// fondo se comprueba con `exact_runs`, un buscador independiente.
    fn background(seed: u64, length: usize) -> Vec<u8> {
        let mut state = seed;
        let mut out: Vec<u8> = Vec::with_capacity(length);
        while out.len() < length {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            let base = b"ACGT"[(state >> 62) as usize];
            let n = out.len();
            if (n >= 1 && out[n - 1] == base) || (n >= 2 && out[n - 2] == base) {
                continue;
            }
            out.push(base);
        }
        out
    }

    /// Buscador de tiradas **exactas** por fuerza bruta, ajeno al detector: para
    /// cada periodo y cada inicio, cuántas copias exactas consecutivas del
    /// motivo leído ahí hay. Devuelve las tiradas maximales que alcanzan
    /// `floor[p-1]` copias, como `(inicio, copias, periodo)`.
    fn exact_runs(seq: &[u8], floor: &[u64; 6]) -> Vec<(usize, usize, usize)> {
        let mut out = Vec::new();
        for period in 1..=6usize {
            for start in 0..seq.len().saturating_sub(period) {
                let motif = &seq[start..start + period];
                let mut copies = 1;
                while start + (copies + 1) * period <= seq.len()
                    && &seq[start + copies * period..start + (copies + 1) * period] == motif
                {
                    copies += 1;
                }
                let maximal_left = start < period || &seq[start - period..start] != motif;
                if maximal_left && copies as u64 >= floor[period - 1] {
                    out.push((start, copies, period));
                }
            }
        }
        out
    }

    /// Añade un tracto puro con flancos **saneados**: cada base de flanco
    /// difiere de la que está un periodo más allá en dirección al tracto y de
    /// su vecina, de modo que a distancia `p` todo el flanco es discrepancia y
    /// el detector no puede alargar el tracto ni salvando un bache. Solo dos
    /// restricciones, a propósito: con tres (también la vecina segunda) la
    /// elección quedaba forzada y el flanco salía periódico, es decir, otra
    /// repetición.
    fn plant(
        sequence: &mut Vec<u8>,
        motif: &[u8],
        copies: usize,
        seed: u64,
        flank: usize,
    ) -> (u64, u64) {
        let start = sequence.len();
        let tract: Vec<u8> = motif
            .iter()
            .copied()
            .cycle()
            .take(motif.len() * copies)
            .collect();
        // Flanco izquierdo ya escrito: se sanea de dentro afuera.
        let p = motif.len();
        let mut joined = sequence.clone();
        joined.extend_from_slice(&tract);
        for offset in 1..=flank.min(start) {
            let position = start - offset;
            let partner = joined[position + p];
            let mut choice = joined[position];
            let mut state = seed ^ (offset as u64);
            while choice == partner || joined.get(position + 1) == Some(&choice) {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                choice = b"ACGT"[(state >> 62) as usize];
            }
            joined[position] = choice;
        }
        // Flanco derecho, saneado de dentro afuera.
        let mut state = seed;
        for _ in 0..flank {
            let position = joined.len();
            loop {
                state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let choice = b"ACGT"[(state >> 62) as usize];
                if choice != joined[position - p] && choice != joined[position - 1] {
                    joined.push(choice);
                    break;
                }
            }
        }
        *sequence = joined;
        (start as u64, (start + tract.len()) as u64)
    }

    fn candidate_triples(detection: &SequenceDetection) -> Vec<(u64, u64, Vec<u8>)> {
        detection
            .candidates()
            .iter()
            .map(|c| (c.start(), c.end(), c.motif().bases().to_vec()))
            .collect()
    }

    /// **Verdad plantada (R-1, oráculo 1).** Un fondo sin repeticiones con
    /// tractos puros en posiciones escritas aquí: los candidatos son
    /// exactamente esos tractos, con su motivo en la fase plantada. El fondo
    /// solo no da ningún candidato, y el buscador de tiradas exactas confirma
    /// que solo hay las plantadas.
    #[test]
    fn planted_pure_tracts_are_exactly_the_candidates() {
        let criteria = DetectionCriteria::default();
        let base = background(0x5EED_0001, 400);
        assert!(
            detect_repeat_candidates(&base, &criteria)
                .candidates()
                .is_empty(),
            "el fondo no debe tener candidatos"
        );
        // Cada tracto alcanza justo su suelo, o uno más.
        let plan: [(&[u8], usize); 6] = [
            (b"A", 8),
            (b"CA", 6),
            (b"GAT", 7),
            (b"ACTG", 6),
            (b"AACTG", 5),
            (b"AGCTTC", 4),
        ];
        let mut sequence = background(0x5EED_0002, 60);
        let mut expected = Vec::new();
        for (index, (motif, copies)) in plan.iter().enumerate() {
            let (start, end) = plant(&mut sequence, motif, *copies, 0xF1A7 + index as u64, 40);
            expected.push((start, end, motif.to_vec()));
            sequence.extend(background(0x5EED_0100 + index as u64, 60));
        }
        // El buscador exacto ve tiradas que alcanzan el suelo solo dentro de
        // las plantadas (también sus rotaciones), y una al inicio de cada una.
        let runs = exact_runs(&sequence, &DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH);
        for &(start, copies, period) in &runs {
            let (start, end) = (start as u64, (start + copies * period) as u64);
            assert!(
                expected.iter().any(|(s, e, _)| *s <= start && end <= *e),
                "tirada exacta fuera de lo plantado: {start}..{end}"
            );
        }
        for (start, _, motif) in &expected {
            assert!(
                runs.iter()
                    .any(|&(s, _, period)| s as u64 == *start && period == motif.len()),
                "el buscador exacto no ve lo plantado en {start}"
            );
        }
        let detection = detect_repeat_candidates(&sequence, &criteria);
        assert_eq!(candidate_triples(&detection), expected);
        assert!(detection.counts().is_balanced());
        for (motif, _) in plan {
            assert_eq!(detection.counts().accepted(motif.len()), 1);
        }
    }

    /// Un tracto plantado **por debajo** de su suelo no sale; el mismo con una
    /// copia más sí. Ejercita el suelo de cada longitud de motivo de punta a
    /// punta, por la función pública.
    #[test]
    fn a_planted_tract_one_copy_below_its_floor_is_not_a_candidate() {
        let criteria = DetectionCriteria::default();
        let motifs: [&[u8]; 6] = [b"T", b"GA", b"CAT", b"AGTC", b"ACGTT", b"AACGTC"];
        for motif in motifs {
            let floor = DEFAULT_MIN_COPIES_BY_MOTIF_LENGTH[motif.len() - 1] as usize;
            for (copies, expected) in [(floor - 1, 0), (floor, 1)] {
                let mut sequence = background(0xAB, 50);
                plant(&mut sequence, motif, copies, 0xCD, 30);
                sequence.extend(background(0xEF, 50));
                let detection = detect_repeat_candidates(&sequence, &criteria);
                assert_eq!(
                    detection.candidates().len(),
                    expected,
                    "{} × {copies}",
                    String::from_utf8_lossy(motif)
                );
            }
        }
    }

    fn interval(start: usize, end: usize, period: u8, score: i32) -> RepeatInterval {
        RepeatInterval {
            start: start as u64,
            end: end as u64,
            period,
            score,
        }
    }

    /// Envuelve un tramo entre dos bases de flanco que no encajan con nada,
    /// para que el borde no decida.
    fn flanked(tract: &[u8]) -> (Vec<u8>, RepeatInterval) {
        let mut bases = vec![b'G'];
        bases.extend_from_slice(tract);
        bases.push(b'G');
        (bases, interval(1, 1 + tract.len(), 0, 0))
    }

    fn classify(tract: &[u8], period: u8) -> Result<RepeatCandidate, DetectionStage> {
        let (bases, mut iv) = flanked(tract);
        iv.period = period;
        classify_interval(&bases, &iv, &DetectionCriteria::default())
    }

    /// **Saturado para el suelo:** dos tramos de **la misma longitud** (12 pb,
    /// 6 copias por longitud), ambos por encima de la pureza y la tirada
    /// mínimas; uno tiene 6 copias completas y el otro 5. Solo el recuento por
    /// copias completas los distingue: contando por longitud, los dos pasarían.
    #[test]
    fn min_copies_counts_complete_copies_not_length() {
        let six = b"CACACACACACA";
        let five = b"CACACACACATT";
        assert_eq!(six.len(), five.len());
        let accepted = classify(six, 2).expect("seis copias completas");
        assert_eq!(accepted.complete_copies(), 6);
        assert_eq!(classify(five, 2), Err(DetectionStage::MinCopies));
        // El de cinco supera por sí solo pureza (10/12) y tirada (5): si no
        // fuera por el suelo, sería candidato.
        let mut lax = DetectionCriteria::default();
        lax.min_copies[1] = 5;
        let (bases, mut iv) = flanked(five);
        iv.period = 2;
        assert!(classify_interval(&bases, &iv, &lax).is_ok());
    }

    /// **Saturado para la pureza:** mismas copias completas (30) y misma tirada
    /// (10), pureza 30/38 ≈ 0,789 frente a 30/37 ≈ 0,811. Y el límite exacto,
    /// 0,8, pasa: la comparación es la de la admisión (`<` rechaza).
    #[test]
    fn min_purity_is_the_only_rule_that_separates_these_tracts() {
        let ten = &b"AAAAAAAAAA"[..];
        let low = [ten, b"CC", ten, b"CCC", ten, b"CCC"].concat();
        let high = [ten, b"CC", ten, b"CC", ten, b"CCC"].concat();
        let exact = [ten, b"CC", ten, b"CCC", ten, b"CCC", ten, b"CC"].concat();
        assert_eq!((low.len(), high.len()), (38, 37));
        assert_eq!(classify(&low, 1), Err(DetectionStage::MinPurity));
        let passed = classify(&high, 1).expect("0,81");
        assert_eq!(
            (passed.complete_copies(), passed.longest_run_units()),
            (30, 10)
        );
        // 40 de 50 = 0,8 justo.
        assert_eq!(exact.len(), 50);
        assert!(classify(&exact, 1).is_ok());
    }

    /// **Saturado para la tirada:** 12 copias completas y la misma pureza
    /// (24/26) en los dos; la tirada más larga es 2 en uno y 3 en el otro.
    #[test]
    fn min_longest_run_is_the_only_rule_that_separates_these_tracts() {
        // Seis tiradas de 2 frente a 3+2+2+2+2+1: en los dos, 12 copias
        // completas, 5 bases sueltas y 29 bases en total.
        let short_runs = b"CACATCACATCACATCACATCACATCACA";
        let long_runs = b"CACACATCACATCACATCACATCACATCA";
        assert_eq!(short_runs.len(), long_runs.len());
        assert_eq!(classify(short_runs, 2), Err(DetectionStage::MinLongestRun));
        let accepted = classify(long_runs, 2).expect("tirada 3");
        assert_eq!(accepted.longest_run_units(), 3);
        let (bases, mut iv) = flanked(short_runs);
        iv.period = 2;
        let lax = DetectionCriteria {
            min_longest_run_units: 2,
            ..DetectionCriteria::default()
        };
        let relaxed = classify_interval(&bases, &iv, &lax).expect("con tirada mínima 2");
        assert_eq!(relaxed.longest_run_units(), 2);
        assert_eq!(
            (relaxed.complete_copies(), relaxed.covered_bases()),
            (accepted.complete_copies(), accepted.covered_bases()),
            "mismas copias y misma pureza: solo la tirada los separa"
        );
    }

    /// **Saturado para el tope:** dos homopolímeros puros, de 200 y 201
    /// unidades. 200 es el tope y se admite; 201 es satélite.
    #[test]
    fn max_units_admits_the_cap_and_rejects_one_more() {
        assert!(classify(&[b'A'; 200], 1).is_ok());
        assert_eq!(classify(&[b'A'; 201], 1), Err(DetectionStage::MaxUnits));
        // En un dinucleótido el tope cuenta unidades, no bases.
        assert!(classify(&b"CA".repeat(200), 2).is_ok());
        assert_eq!(
            classify(&b"CA".repeat(201), 2),
            Err(DetectionStage::MaxUnits)
        );
    }

    /// **Saturado para el borde:** el mismo tramo puro, con flanco a los dos
    /// lados, pegado al inicio, pegado al final y pegado a una `N`.
    #[test]
    fn edge_requires_one_acgt_base_on_each_side() {
        let tract = b"CACACACACACACACA";
        let criteria = DetectionCriteria::default();
        let iv = |start: usize| interval(start, start + tract.len(), 2, 0);
        let with = |left: &[u8], right: &[u8]| [left, tract, right].concat();
        assert!(classify_interval(&with(b"G", b"G"), &iv(1), &criteria).is_ok());
        assert_eq!(
            classify_interval(&with(b"", b"G"), &iv(0), &criteria),
            Err(DetectionStage::Edge)
        );
        assert_eq!(
            classify_interval(&with(b"G", b""), &iv(1), &criteria),
            Err(DetectionStage::Edge)
        );
        assert_eq!(
            classify_interval(&with(b"N", b"G"), &iv(1), &criteria),
            Err(DetectionStage::Edge)
        );
        assert_eq!(
            classify_interval(&with(b"G", b"N"), &iv(1), &criteria),
            Err(DetectionStage::Edge)
        );
    }

    /// Un motivo con una base que el simulador no puede representar se cuenta
    /// en su propia etapa, antes que todas las demás.
    #[test]
    fn a_motif_with_a_non_acgt_base_is_counted_apart() {
        let tract = b"ACGTNACGTNACGTNACGTNACGTN";
        assert_eq!(classify(tract, 5), Err(DetectionStage::NonCanonicalMotif));
        assert_eq!(
            classify(b"cacacacacacacaca", 2),
            Err(DetectionStage::NonCanonicalMotif)
        );
    }

    fn candidate(start: u64, end: u64, motif: &[u8], score: i32) -> RepeatCandidate {
        RepeatCandidate {
            start,
            end,
            motif: RepeatMotif::new(motif).unwrap(),
            score,
            covered_bases: end - start,
            longest_run_units: (end - start) / motif.len() as u64,
        }
    }

    /// Resuelve en los dos órdenes de entrada y comprueba que el ganador es el
    /// mismo: es lo que distingue una regla fija de «el primero que llegue».
    fn winner_in_both_orders(a: RepeatCandidate, b: RepeatCandidate) -> RepeatCandidate {
        let (forward, lost_forward) = resolve_overlaps(vec![a.clone(), b.clone()]);
        let (backward, lost_backward) = resolve_overlaps(vec![b, a]);
        assert_eq!(forward, backward, "el ganador no depende del orden");
        assert_eq!(
            (forward.len(), lost_forward.len(), lost_backward.len()),
            (1, 1, 1)
        );
        forward.into_iter().next().unwrap()
    }

    /// **Desempate de solapes, cada nivel saturado:** en cada pareja todos los
    /// criterios anteriores empatan y solo el del nivel decide, en los dos
    /// órdenes de entrada.
    #[test]
    fn overlap_ties_follow_the_fixed_rule_in_any_input_order() {
        // 1. Mayor puntuación, aunque tenga menos bases.
        let won =
            winner_in_both_orders(candidate(10, 30, b"CA", 50), candidate(20, 50, b"GAT", 40));
        assert_eq!(won.start(), 10);
        // 2. A igual puntuación, más bases.
        let won =
            winner_in_both_orders(candidate(10, 30, b"CA", 40), candidate(20, 50, b"GAT", 40));
        assert_eq!(won.start(), 20);
        // 3. A igual puntuación y bases, periodo menor.
        let won =
            winner_in_both_orders(candidate(10, 34, b"GATC", 40), candidate(20, 44, b"CA", 40));
        assert_eq!(won.motif().bases(), b"CA");
        // 4. A igual puntuación, bases y periodo, el de la izquierda.
        let won = winner_in_both_orders(candidate(20, 44, b"CA", 40), candidate(10, 34, b"GT", 40));
        assert_eq!(won.start(), 10);
    }

    /// Dos tramos que solo se tocan no se solapan; y un perdedor no bloquea a
    /// un tercero que solo lo solapaba a él.
    #[test]
    fn touching_tracts_are_kept_and_losers_do_not_block() {
        let (kept, lost) = resolve_overlaps(vec![
            candidate(10, 20, b"CA", 10),
            candidate(20, 30, b"GT", 10),
        ]);
        assert_eq!((kept.len(), lost.len()), (2, 0));
        // A gana a B; C solo solapa a B y queda.
        let (kept, lost) = resolve_overlaps(vec![
            candidate(18, 40, b"GAT", 20),
            candidate(0, 20, b"CA", 30),
            candidate(38, 60, b"AC", 10),
        ]);
        assert_eq!(
            kept.iter().map(RepeatCandidate::start).collect::<Vec<_>>(),
            vec![0, 38]
        );
        assert_eq!(lost.len(), 1);
    }

    /// El recuento cuadra y atribuye cada descarte a su etapa, por la función
    /// pública: suelo, tope y borde con los criterios por omisión.
    #[test]
    fn counts_attribute_each_discard_to_its_stage() {
        let mut sequence = background(0x77, 40);
        // Suelo: (CA)×5.
        plant(&mut sequence, b"CA", 5, 1, 20);
        sequence.extend(background(0x78, 40));
        // Tope: poli-A de 201.
        plant(&mut sequence, b"A", 201, 2, 20);
        sequence.extend(background(0x79, 40));
        // Aceptado: (GAT)×7.
        plant(&mut sequence, b"GAT", 7, 3, 20);
        sequence.extend(background(0x7A, 40));
        // Borde: (TC)×10 pegado al final.
        sequence.extend(b"TC".repeat(10));
        let detection = detect_repeat_candidates(&sequence, &DetectionCriteria::default());
        let counts = detection.counts();
        assert!(counts.is_balanced());
        assert_eq!(counts.accepted_total(), 1);
        assert_eq!(counts.accepted(3), 1);
        assert_eq!(counts.discarded(DetectionStage::MaxUnits, 1), 1);
        assert_eq!(counts.discarded(DetectionStage::Edge, 2), 1);
        assert!(counts.discarded(DetectionStage::MinCopies, 2) >= 1);
    }

    /// Recuentos exactos de **pureza, tirada y solape** por la función
    /// pública, cada uno en una longitud de motivo propia:
    ///
    /// - pureza, periodo 3: `((GAT)×3 GCT)×4`, 12 copias completas de 16 y
    ///   pureza 0,75; su tirada (3) tampoco alcanza el mínimo de 5 de esta
    ///   prueba, de modo que **fija también que la pureza va antes que la
    ///   tirada**;
    /// - tirada, periodo 4: `((ACTG)×4 ATTG)×3 (ACTG)×4`, pureza 0,84 y
    ///   tirada 4;
    /// - solape, periodo 2 contra 3: `(CA)×10 G (CAG)×7`; el tramo de periodo
    ///   3 empieza en el último `CA` del dinucleótido (`CAG` con la `G`) y lo
    ///   pisa dos bases; el `(CA)×10` (puntuación 36) pierde ante el `(CAG)×8`
    ///   (puntuación 42). Un poli-A delante de un `(CA)n` **no** sirve de
    ///   solape: el tramo de periodo 2 se alarga hacia el poli-A, su motivo
    ///   sale `AA` y el filtro de motivos lo descarta entero (pendiente de
    ///   decisión en la pausa).
    #[test]
    fn counts_purity_run_and_overlap_through_the_detector() {
        let criteria = DetectionCriteria {
            min_longest_run_units: 5,
            ..DetectionCriteria::default()
        };
        let mut sequence = background(0x101, 40);
        let purity_tract = [&b"GATGATGATGCT"[..]; 4].concat();
        sequence.extend_from_slice(&purity_tract);
        sequence.extend(background(0x102, 40));
        let run_block = [&b"ACTGACTGACTGACTG"[..], b"ATTG"].concat();
        let run_tract = [&run_block[..], &run_block, &run_block, b"ACTGACTGACTGACTG"].concat();
        sequence.extend_from_slice(&run_tract);
        sequence.extend(background(0x103, 40));
        sequence.extend(b"CA".repeat(10));
        sequence.push(b'G');
        sequence.extend(b"CAG".repeat(7));
        sequence.extend(background(0x104, 40));
        let detection = detect_repeat_candidates(&sequence, &criteria);
        let counts = detection.counts();
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!(
            counts.discarded(DetectionStage::MinPurity, 3),
            1,
            "{counts:?}"
        );
        assert_eq!(
            counts.discarded_total(DetectionStage::MinPurity),
            1,
            "{counts:?}"
        );
        assert_eq!(
            counts.discarded(DetectionStage::MinLongestRun, 4),
            1,
            "{counts:?}"
        );
        assert_eq!(
            counts.discarded_total(DetectionStage::MinLongestRun),
            1,
            "{counts:?}"
        );
        assert_eq!(
            counts.discarded(DetectionStage::Overlap, 2),
            1,
            "{counts:?}"
        );
        assert_eq!(
            counts.discarded_total(DetectionStage::Overlap),
            1,
            "{counts:?}"
        );
        assert_eq!(
            candidate_triples(&detection)
                .into_iter()
                .map(|(_, _, motif)| motif)
                .collect::<Vec<_>>(),
            vec![b"CAG".to_vec()],
            "solo queda el (CAG)n, que gana el solape"
        );
    }

    /// **El orden de las etapas.** Cada tramo falla **dos** reglas a la vez y
    /// se cuenta en la primera: suelo antes que borde, suelo antes que tope,
    /// pureza antes que tirada y tope antes que borde. Por `classify_interval`,
    /// con intervalos escritos a mano; y el primer caso, también por la
    /// función pública.
    #[test]
    fn a_tract_failing_two_rules_is_counted_at_the_earlier_stage() {
        let criteria = DetectionCriteria::default();
        // Suelo y borde: (CA)×5 al inicio de la secuencia.
        let bases = [&b"CACACACACA"[..], b"GG"].concat();
        assert_eq!(
            classify_interval(&bases, &interval(0, 10, 2, 0), &criteria),
            Err(DetectionStage::MinCopies)
        );
        // Suelo y tope: (CA)×5 con un tope de 4 unidades.
        let capped = DetectionCriteria {
            max_units: 4,
            ..DetectionCriteria::default()
        };
        assert_eq!(
            classify(b"CACACACACA", 2).map(|_| ()),
            Err(DetectionStage::MinCopies)
        );
        let (flanked_bases, mut iv) = flanked(b"CACACACACA");
        iv.period = 2;
        assert_eq!(
            classify_interval(&flanked_bases, &iv, &capped).map(|_| ()),
            Err(DetectionStage::MinCopies)
        );
        // Pureza y tirada: (CACATT)×6, pureza 4/6 y tirada 2.
        assert_eq!(
            classify(&b"CACATT".repeat(6), 2).map(|_| ()),
            Err(DetectionStage::MinPurity)
        );
        // Tope y borde: (CA)×201 al inicio.
        let long = [&b"CA".repeat(201)[..], b"GG"].concat();
        assert_eq!(
            classify_interval(&long, &interval(0, 402, 2, 0), &criteria).map(|_| ()),
            Err(DetectionStage::MaxUnits)
        );
        // Por la función pública: (CA)×5 pegado al final cuenta como suelo.
        let mut sequence = background(0x201, 40);
        let tail = sequence.len();
        sequence.extend(b"CA".repeat(5));
        let detection = detect_repeat_candidates(&sequence, &criteria);
        assert!(detection.counts().discarded(DetectionStage::MinCopies, 2) >= 1);
        assert_eq!(detection.counts().discarded(DetectionStage::Edge, 2), 0);
        assert!(sequence[tail..].starts_with(b"CA"));
    }
}
