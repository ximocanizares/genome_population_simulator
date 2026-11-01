//! Detección de microsatélites (repeticiones en tándem de periodo corto) sobre
//! una secuencia de bases, y conversión de lo detectado en candidatos a locus.
//!
//! # Procedencia y licencia
//!
//! El **núcleo** de este módulo —[`PeriodRange`], [`ScanParams`],
//! [`RepeatInterval`], la tabla `CANONICAL`, `is_primitive_motif`, el pase de
//! Ruzzo y Tompa, el corte por cambio de motivo (`split_at_motif_changes`),
//! [`find_tandem_repeats`] y [`find_tandem_repeats_for_copy_floor`], con sus
//! pruebas— está **copiado** del fichero `src/tandem_repeat.rs` del proyecto
//! `pop_var_caller` (repositorio `join_per_sample_vcfs`), commit de origen
//! `fa215c2365ca05ed15fac056285b63b9771c6c49` (30/10/2026).
//!
//! La primera copia (28/10/2026) se hizo de
//! `950909f9962f463ca157264fa36193703ed2b467`. Aquel núcleo fundía en un solo
//! intervalo dos repeticiones pegadas de motivos distintos (`A×12 T×12`, o
//! `(CA)×10 (GT)×10`), con el motivo de la primera, y el tramo fundido no
//! pasaba la pureza: se perdían las dos. El 30/10/2026 se sustituyó por el de
//! `fa215c23`, que trae del original cuatro cambios: el corte por cambio de
//! motivo (`ccd57ea2` y sus correcciones de revisión, `306e3bd3`), la emisión
//! de un tramo unido junto con los trozos de los que se unió (`c237378b`) y
//! la revisión de rendimiento (`cdce162c`: la pila de Ruzzo y Tompa emite los
//! segmentos en cuanto son definitivos, los segmentos cortos no se parten y
//! el barrido con suelo de copias).
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
//!   `RegionScanner`, `scan_windowed` y compañía, con el tipo `RegionSpan`:
//!   el pase de periodos múltiplos guarda aquí el tramo como una tupla
//!   `(inicio, fin)`) y la función muerta `bases_match`. Tampoco se trae
//!   quien elige, en el original, entre un tramo unido y sus trozos
//!   (`segment_criteria::choose_readings`, que no es del buscador): aquí las
//!   dos lecturas entran en [`detect_repeat_candidates`] como intervalos
//!   corrientes, cada una se juzga sola y la regla de solapes se queda con
//!   una. Los comentarios se han traducido y las pruebas copiadas no se han
//!   tocado.
//! - **Única desviación del algoritmo (decisión del responsable, 29/10/2026):
//!   un tramo recortado empieza en la primera fase primitiva fuera de la tira
//!   que lo precede.** Nació con el núcleo de `950909f`: allí, si el motivo
//!   leído al inicio del segmento no era primitivo, el intervalo **entero** se
//!   descartaba, y eso perdía un microsatélite pegado a un homopolímero de su
//!   misma base (en `A×8 (GA)×10` el segmento de periodo 2 empieza dentro del
//!   poli-A, lee `AA` y el `(GA)n` desaparecía). La regla: el inicio avanza
//!   base a base hasta la primera fase cuyo motivo de periodo `p` sea
//!   primitivo **y cuya base sea distinta de la anterior**; solo se descarta
//!   si no queda ninguna antes de que el tramo baje del mínimo de copias. La
//!   segunda condición saca el inicio de la tira de una sola base que precede
//!   al microsatélite: `A×8 (GA)×10` da `[8, 28)` con motivo `GA`, que solo
//!   toca al poli-A `[0, 8)`, y se conservan los dos.
//!
//!   **Con el núcleo de `fa215c23` sigue haciendo falta, y cambia dónde se
//!   aplica.** El corte por cambio de motivo ya separa el poli-A del `(GA)n`,
//!   pero el trozo nuevo empieza en la primera coincidencia tras la unión,
//!   que es la **última base de la tira**: sin la desviación, `A×8 (GA)×10`
//!   da `[7, 28)` con motivo `AG`, `A×6 (GA)×10 C` da `[5, 26)` y `T×12
//!   (CTT)×7` da `[10, 33)` con `TTC`, y la regla de solapes deja uno de los
//!   dos. Así que la regla se aplica a todo tramo **recortado**, que ahora
//!   son dos: el que empieza con un motivo no primitivo, como antes, y el
//!   trozo nacido de un corte, que no empieza donde su segmento. Un tramo que
//!   empieza donde su segmento y con motivo primitivo sale exactamente igual
//!   que en el original, y uno puro de un motivo no primitivo (un alias de un
//!   homopolímero o de un dinucleótido) no tiene ninguna fase primitiva y se
//!   sigue descartando. La puntuación del tramo recortado **se recalcula**
//!   restando la del prefijo quitado, de modo que siempre es la suma de las
//!   comparaciones del propio tramo (la usa el desempate de solapes). Ver el
//!   cierre `emit` de [`find_tandem_repeats_for_copy_floor`].
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
//! **Los solapes (DP-SSR-05, reescrita el 01/11/2026).** Un microsatélite suele
//! acabar con una copia parcial del motivo que coincide con el principio de
//! lo que viene detrás, y el buscador entrega los dos tramos pisándose esa
//! copia. La regla ordena por prioridad (puntuación, longitud, periodo menor,
//! posición) y, de cada aspirante: lo que no pisa a un aceptado se acepta; lo
//! que pisa por cada lado **una copia de su propio motivo o menos** (solape ≤
//! periodo) **se recorta** a la parte que no solapa y vuelve a pasar por todos
//! los criterios; lo que pisa más de una copia por algún lado, o contiene o
//! queda cubierto por un aceptado, se descarta entero. Un recortado que cae se
//! cuenta en «solape». Ver `resolve_overlaps`.
//!
//! Desde la entrega 3 (29/10/2026) la usa la ejecución:
//! `EffectiveConfig::with_detected_repeat_loci` barre cada secuencia de la
//! referencia final con [`detect_repeat_candidates_in_regions`] cuando el
//! documento declara la sección `[reference.repeat_detection]`. Los criterios
//! son los del apartado 11.7 del plan, decididos por el responsable; siguen
//! siendo parámetros porque la sección los declara.

use std::cell::Cell;
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
/// bruta (`ruzzo_tompa_matches_brute_force_on_pseudo_random`). La pila guarda
/// solo los segmentos abiertos desde el último apilado por la regla 2, que es
/// cuando todos los de debajo pasan a ser definitivos y se emiten (ver esa
/// rama).
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
                    // Regla 2: sin segmento menor a la izquierda → `cur` va
                    // aparte. **Y todo lo que hay debajo es definitivo**: todo
                    // `l` apilado es al menos `cur.l`, así que una búsqueda
                    // posterior que fuera a dar con uno de ellos da antes con
                    // `cur`, que está encima, y que `cur` se encuentre o se
                    // absorba no los toca. Emitirlos ya deja la pila en los
                    // segmentos desde el último apilado así, en vez de los de
                    // una secuencia entera (en el original, 1,34 GB en el
                    // cromosoma 1 del tomate a periodo 1 cuando esperaban al
                    // final). Salen en el orden en que los habría emitido el
                    // vaciado final.
                    for seg in &stack {
                        emit(seg.start, seg.end, seg.r - seg.l);
                    }
                    stack.clear();
                    cur.prev_smaller = 0;
                    break;
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

/// Parte un segmento de puntuación máxima allí donde **cambia el motivo** y
/// emite cada trozo como `(k0, k1, score)`, en la misma forma de índices de
/// puntuación que emite [`maximal_scoring_subsequences`].
///
/// La comparación a distancia `p` solo pregunta «¿es esta base la de `p`
/// atrás?», lo que es igual de cierto dentro de una tira de `A` que dentro de
/// una de `T`. Así, `A×12 T×12` puntúa 11 coincidencias, una discrepancia y 11
/// coincidencias en periodo 1, y el segmento cruza la unión: un solo tramo, de
/// motivo `A`, con la mitad de `T`. La pureza lo rechaza después y se pierden
/// **las dos** repeticiones. Lo mismo pasa con `(CA)×10 (GT)×10` en periodo 2
/// o con `(GA)×12` pegado a `C×12`.
///
/// La prueba que separa un cambio de motivo de las interrupciones que deben
/// quedar dentro es esta. Un **bloque sólido** es una tirada de coincidencias
/// consecutivas a distancia `p` lo bastante larga para ser una repetición según
/// el propio suelo del detector —al menos `min_copies` copias, y nunca menos
/// de dos, para que su motivo quede determinado del todo—. Su motivo se
/// compara salvo rotación:
///
/// - una sustitución deja menos de `p` coincidencias entre sus dos
///   discrepancias, así que nunca forma un bloque sólido, y los bloques de los
///   dos lados comparten motivo;
/// - un indel desplaza la fase, así que el bloque de después lleva una
///   **rotación** del motivo de antes: el mismo motivo;
/// - una repetición distinta pegada a esta lleva un motivo que no es rotación
///   suya: un corte.
///
/// Los bloques sólidos consecutivos del mismo motivo forman un trozo; las
/// bases entre dos trozos (las discrepancias de la unión) no son de ninguno.
///
/// **Una tira corta de otro motivo dentro de una repetición es una
/// interrupción, no una segunda repetición.** `ATATATATACACATATATAT` es un
/// `AT` con dos copias cambiadas a `AC`; la pareja `AC` es un bloque sólido,
/// pero cortar ahí deja trozos demasiado cortos para un suelo de copias y la
/// repetición se pierde (en el original, 876 segmentos de esta forma en el
/// tomate). Por eso, cuando el motivo de un trozo vuelve, los trozos de en
/// medio se **absorben**: el trozo anterior más cercano del mismo motivo, todo
/// lo que le sigue y el nuevo trozo pasan a ser uno, siempre que el tramo
/// unido, comparado base a base con ese motivo repetido en fase, sea al menos
/// tan puro como la puntuación exige a un tramo (un cociente `r` sostiene un
/// tramo hasta una pureza de `r / (1 + r)`: 7/9 ≈ 0,78 con los pesos por
/// omisión). Ese ejemplo coincide con `AT` en 18 de 20 bases, 0,9, y queda
/// entero; también `(TAC)×2 (TGC)×4 (TAC)×3`, 0,85, cuyo motivo interior
/// difiere del exterior en una base, y `(AAT)×4 (GAT)×3 (GGT)×2 (AAT)×4`,
/// 0,82, con dos interrupciones distintas seguidas. Base a base y no trozo a
/// trozo, porque así mide la pureza quien lo lee después. Ninguna perilla
/// nueva: la tolerancia es la que ya fija la puntuación.
///
/// La fase es la del **último** bloque sólido del trozo exterior, de modo que
/// un indel anterior en ese trozo hace que las bases de antes cuenten como
/// impuras, y la interrupción siguiente puede no absorberse. Es deliberado: la
/// pureza se mide en una sola fase desde el inicio del tramo, así que un tramo
/// unido con un indel dentro se rechazaría entero, como antes de que existiera
/// este corte; cortar al menos conserva los trozos puros.
///
/// **Un trozo de cualquiera de los dos extremos nunca se absorbe**, por corto
/// que sea: no le sigue nada del primer motivo, así que es donde una
/// repetición acaba y empieza otra. `G×60 A×15` es `G` en más de 0,8 de sus
/// bases, y siguen siendo dos repeticiones.
///
/// Un segmento de un solo motivo —todo segmento que no sea una unión— y uno
/// cuyos trozos se absorbieron todos se emiten sin cambios, puntuación
/// incluida. Quien llama solo envía segmentos con alguna discrepancia (uno sin
/// ninguna no puede cambiar de motivo) y lo bastante largos para quedarse.
///
/// **Un tramo unido se emite con los trozos de los que se unió**, anidados
/// dentro de él en el mismo periodo: las dos lecturas, para que quien lo lea
/// elija por su propia pureza mínima. En este proyecto elige la cadena de
/// [`detect_repeat_candidates`]: cada lectura se juzga sola y la regla de
/// solapes se queda con una.
///
/// **Fuera de línea a propósito**: en línea dentro del cierre por segmento del
/// barrido, ese cierre crecía demasiado para meterse en el bucle del barrido y
/// toda la construcción pagaba una llamada por segmento (medido en el
/// original sobre el cromosoma 1 del tomate).
#[inline(never)]
fn split_at_motif_changes(
    canonical: &[u8],
    p: usize,
    k0: usize,
    k1: usize,
    score: i64,
    params: &ScanParams,
    emit: &mut impl FnMut(usize, usize, i64),
) {
    let is_match = |k: usize| {
        let h = canonical[k + p];
        h != 0 && h == canonical[k]
    };
    // Un bloque de `len` coincidencias consecutivas cubre `len + p` bases, así
    // que `min_copies` copias necesitan `(min_copies - 1) · p` coincidencias;
    // al menos `p`, para que el bloque tenga dos copias enteras.
    let solid_len = p * (params.min_copies.max(2) as usize - 1);
    let same_motif = |a: usize, b: usize| {
        let (x, y) = (&canonical[a..a + p], &canonical[b..b + p]);
        (0..p).any(|r| (0..p).all(|i| x[i] == y[(i + r) % p]))
    };

    let reward = i64::from(params.match_reward);
    let penalty = i64::from(params.mismatch_penalty);
    // Cuántas bases de `from..to` coinciden con el motivo leído en `motif_at`,
    // repetido en la fase de ese bloque —antes y después de él—, de modo que un
    // trozo que empieza con un fragmento más corto que un bloque sólido se mide
    // contra el motivo que de verdad repite.
    let in_phase = |motif_at: usize, from: usize, to: usize| {
        (from..to)
            .filter(|&x| canonical[x] == canonical[motif_at + (x % p + p - motif_at % p) % p])
            .count() as i64
    };
    // Todos los trozos tal como se cortaron, `(primer índice de puntuación,
    // último)`: la segunda lectura de un tramo unido. Vacío —y nunca reservado—
    // para un segmento de un solo motivo.
    let mut split: Vec<(usize, usize)> = Vec::new();
    // Los trozos tras la absorción, de izquierda a derecha, como `(primer
    // índice de puntuación, último, dónde se lee su motivo, bases de su tramo
    // que coinciden con ese motivo en fase, el primero de `split` que
    // contiene)`. El tramo de un trozo son sus índices ensanchados en `p + 1`,
    // igual que el del segmento.
    let mut pieces: Vec<(usize, usize, usize, i64, usize)> = Vec::new();
    // Apila un trozo y absorbe toda interrupción que cierre: el trozo anterior
    // más cercano del mismo motivo, con todo lo de en medio, cuando el tramo
    // unido coincide con ese motivo en fase con la tolerancia de la
    // puntuación —coincidencias / bases >= penalización / (recompensa +
    // penalización), en enteros—. Solo se comparan las bases nuevas, así que
    // una tirada larga de interrupciones absorbidas cuesta un pase y no uno por
    // unión.
    let push = |pieces: &mut Vec<(usize, usize, usize, i64, usize)>,
                split: &mut Vec<(usize, usize)>,
                (a, b, m): (usize, usize, usize)| {
        pieces.push((a, b, m, in_phase(m, a, b + p + 1), split.len()));
        split.push((a, b));
        loop {
            let right = pieces[pieces.len() - 1];
            let Some(j) = (0..pieces.len() - 1)
                .rev()
                .find(|&j| same_motif(pieces[j].2, right.2))
            else {
                return;
            };
            let left = pieces[j];
            let (left_end, right_end) = (left.1 + p + 1, right.1 + p + 1);
            let matches = left.3 + in_phase(left.2, left_end, right_end.max(left_end));
            if matches * (reward + penalty) < penalty * (right_end - left.0) as i64 {
                return;
            }
            pieces.truncate(j);
            pieces.push((left.0, right.1, left.2, matches, left.4));
        }
    };
    // El trozo abierto: dónde empieza, el final de su último bloque sólido y el
    // inicio de ese bloque (cuyas `p` primeras bases son el motivo actual del
    // trozo).
    let mut piece_start = k0;
    let mut motif_at: Option<usize> = None;
    let mut last_solid_end = k0;
    let mut k = k0;
    while k <= k1 {
        if !is_match(k) {
            k += 1;
            continue;
        }
        let b0 = k;
        while k <= k1 && is_match(k) {
            k += 1;
        }
        let b1 = k - 1;
        if b1 + 1 - b0 < solid_len {
            continue;
        }
        match motif_at {
            Some(m) if !same_motif(m, b0) => {
                push(&mut pieces, &mut split, (piece_start, last_solid_end, m));
                piece_start = b0;
            }
            _ => {}
        }
        motif_at = Some(b0);
        last_solid_end = b1;
    }
    let (Some(last_motif), false) = (motif_at, pieces.is_empty()) else {
        emit(k0, k1, score);
        return;
    };
    push(&mut pieces, &mut split, (piece_start, k1, last_motif));
    // La puntuación propia de cada trozo, sobre sus propios índices, como la
    // del segmento; la del segmento para un trozo que es el segmento entero.
    let score_of = |a: usize, b: usize| {
        if (a, b) == (k0, k1) {
            score
        } else {
            (a..=b)
                .map(|k| if is_match(k) { reward } else { -penalty })
                .sum()
        }
    };
    // Un trozo unido se emite **con los trozos de los que se unió**: las dos
    // lecturas, y la pureza de quien lee elige entre ellas.
    for (i, &(a, b, _, _, first_split)) in pieces.iter().enumerate() {
        emit(a, b, score_of(a, b));
        let end_split = pieces.get(i + 1).map_or(split.len(), |next| next.4);
        if end_split - first_split > 1 {
            for &(sa, sb) in &split[first_split..end_split] {
                emit(sa, sb, score_of(sa, sb));
            }
        }
    }
}

/// Encuentra todo intervalo de repetición en tándem de `seq` cuyo periodo esté
/// en `periods`.
///
/// Para cada periodo `p` puntúa la posición `j` (`j >= p`) con `+match_reward`
/// cuando `seq[j]` y `seq[j-p]` son la misma base ACGT y con
/// `-mismatch_penalty` en otro caso (una base no ACGT nunca coincide), toma
/// cada segmento de puntuación máxima de Ruzzo y Tompa, lo parte donde cambia
/// el motivo (`split_at_motif_changes`), traduce cada trozo al tramo
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
///
/// La única excepción a «una vez» son las dos lecturas de un tramo unido, que
/// salen anidadas en el mismo periodo (ver `split_at_motif_changes`).
#[must_use]
pub fn find_tandem_repeats(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
) -> Vec<RepeatInterval> {
    find_tandem_repeats_for_copy_floor(seq, periods, params, |_| 0)
}

/// [`find_tandem_repeats`] para quien va a descartar todo tramo con menos de
/// `copy_floor(period)` copias: todo tramo que devuelve en ese suelo o por
/// encima es exactamente uno de los que devuelve `find_tandem_repeats`, y de
/// los que quedan por debajo omite los demasiado cortos para importar.
///
/// **No es simplemente el suelo, porque un tramo por debajo puede quitar a uno
/// por encima.** El pase de periodos múltiplos quita un tramo del que un tramo
/// de periodo menor cubre la mitad, así que un tramo de periodo 2 de 9 bases,
/// por debajo de un suelo de 5 copias, puede quitar a uno de periodo 4 de 16
/// que supera un suelo de 4. Solo se omite lo que no puede hacer eso: un tramo
/// de periodo `p` hace falta cuando supera su propio suelo, o cuando mide al
/// menos la mitad del tramo más corto que hace falta en algún múltiplo de `p`
/// del intervalo; esas longitudes se calculan del periodo más largo hacia
/// abajo. Un tramo más corto no quita nada que haga falta, y los tramos a los
/// que habría protegido son más cortos todavía, así que omitirlo no cambia
/// nada de lo que supera el suelo.
///
/// En el original existe por la memoria: en el cromosoma 1 del tomate su
/// catálogo conserva 774 585 de los 23,6 millones de tramos que devuelve
/// `find_tandem_repeats`. Aquí la conversión en candidatos usa el barrido
/// entero (suelo 0), porque cuenta en qué etapa cae **cada** tramo detectado.
pub fn find_tandem_repeats_for_copy_floor(
    seq: &[u8],
    periods: PeriodRange,
    params: &ScanParams,
    copy_floor: impl Fn(u8) -> u32,
) -> Vec<RepeatInterval> {
    // `needed[p]`: el tramo más corto de periodo `p` que puede importar.
    let mut needed = vec![0u64; usize::from(periods.max()) + 1];
    for p in (periods.min()..=periods.max()).rev() {
        let own = u64::from(copy_floor(p)) * u64::from(p);
        needed[usize::from(p)] = (2 * usize::from(p)..=usize::from(periods.max()))
            .step_by(usize::from(p))
            .map(|m| needed[m].div_ceil(2))
            .fold(own, u64::min);
    }

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
        // El tramo más corto que `emit` conserva en este periodo. Todo trozo de
        // un corte cae dentro de su segmento, así que un segmento más corto no
        // da nada, se parta como se parta.
        let shortest = (u64::from(params.min_copies) * p as u64).max(needed[p]) as usize;
        // Dónde empieza el segmento que se está emitiendo: un tramo que empieza
        // ahí es el segmento o su primer trozo; uno que empieza después es un
        // trozo nacido de un corte. Lo usa la desviación propia de `emit`.
        let segment_start = Cell::new(0usize);
        let mut emit = |k0: usize, k1: usize, score: i64| {
            // Segmento [k0, k1] → tramo [k0, k1 + p + 1): la primera base
            // implicada es `j0 - p = k0` y la última `j1 = k1 + p`, así que el
            // fin exclusivo es `k1 + p + 1`.
            let tract_end = k1 + p + 1;
            // Solo repeticiones de periodo PRIMITIVO: `AAAAAA` se detecta en
            // periodo 2, pero su motivo `AA` se reduce a `A`; es un
            // homopolímero, no un dinucleótido, e informarlo como periodo 2
            // sería incumplir el contrato del detector.
            //
            // **Desviación respecto al original** (ver la cabecera): el inicio
            // avanza hasta la primera fase cuyo motivo sea primitivo **y cuya
            // base sea distinta de la anterior**, mientras el tramo conserve el
            // mínimo de copias, cuando el motivo del inicio no es primitivo o
            // cuando el tramo es un trozo nacido de un corte. La segunda
            // condición saca el inicio de la tira de una sola base que lo
            // precede, para que tira y microsatélite solo se toquen. Sin
            // ninguna posición así, se descarta como en el original. Un tramo
            // que empieza donde su segmento y con motivo primitivo no cambia.
            let at_segment_start = k0 == segment_start.get();
            let Some(start) = (k0..=k1)
                .take_while(|&start| (tract_end - start) >= shortest)
                .find(|&start| {
                    is_primitive_motif(&seq[start..start + p])
                        && ((at_segment_start && start == k0)
                            || canonical[start] != canonical[start - 1])
                })
            else {
                return;
            };
            // La puntuación del tramo recortado es la del trozo menos la del
            // prefijo quitado: la suma de las comparaciones del propio tramo,
            // igual que un tramo sin recortar.
            let removed: i64 = (k0..start)
                .map(|k| {
                    if here[k] != 0 && here[k] == back[k] {
                        reward
                    } else {
                        penalty
                    }
                })
                .sum();
            out.push(RepeatInterval {
                start: start as u64,
                end: tract_end as u64,
                period,
                score: i32::try_from(score - removed).unwrap_or(i32::MAX),
            });
        };
        maximal_scoring_subsequences(scores, |k0, k1, score| {
            if k1 + p + 1 - k0 < shortest {
                return;
            }
            segment_start.set(k0);
            // Sin discrepancia no hay cambio de motivo: la puntuación de un
            // segmento es toda recompensas justo cuando todas sus comparaciones
            // coincidieron, y se emite entero sin recorrerlo.
            if score == reward * (k1 - k0 + 1) as i64 {
                emit(k0, k1, score);
            } else {
                split_at_motif_changes(&canonical, p, k0, k1, score, params, &mut emit);
            }
        });
    }

    // Descarta la re-detección con periodo múltiplo del mismo tracto. `iv` (el
    // periodo mayor) sobra si un periodo menor `a` lo divide y **cubre al menos
    // la mitad** de él: «cubre la mayor parte» es lo que dice *mismo tracto* y
    // no solo *solapado*. Ninguna política de copias aquí: solo cuánto de `iv`
    // explica el menor.
    out.sort_by_key(|iv| (iv.period, iv.start, std::cmp::Reverse(iv.end)));

    // Supervivientes agrupados por periodo, que es lo que convierte la prueba
    // de redundancia en una búsqueda binaria. Tres hechos la hacen exacta:
    //
    // 1. **Solo un periodo menor que lo divida puede eliminar a `iv`**: a lo
    //    sumo `{1, 2, 3}`, para el periodo 6.
    // 2. **Esas listas ya están completas** al llegar a `iv`: `out` está
    //    ordenado primero por periodo.
    // 3. **Dentro de un periodo los tramos buscados están ordenados por inicio
    //    Y por fin**: los segmentos de Ruzzo y Tompa son disjuntos, y cada
    //    tramo alarga su segmento en el mismo `p + 1`. Es la precondición de
    //    `partition_point`. La única excepción son los trozos de un tramo
    //    unido, emitidos junto a él y anidados dentro: se dejan fuera de las
    //    listas buscadas, sin perder nada, porque el tramo que encierra a un
    //    trozo solapa a cualquier otro al menos tanto como el trozo. La
    //    ordenación pone primero al que encierra, así que un trozo se reconoce
    //    porque no acaba más allá del alcance acumulado.
    let max_period = out.iter().map(|iv| iv.period).max().unwrap_or(0) as usize;
    // Solo un periodo con algún múltiplo en el intervalo se busca (hecho 1),
    // así que solo esos llevan lista, y la búsqueda solo lee el tramo `(inicio,
    // fin)`. Cada lista se reserva una vez con el recuento de su periodo en
    // `out`.
    let mut by_period: Vec<Vec<(u64, u64)>> = (0..=max_period / 2)
        .map(|period| {
            let from = out.partition_point(|iv| usize::from(iv.period) < period);
            let to = out.partition_point(|iv| usize::from(iv.period) <= period);
            Vec::with_capacity(to - from)
        })
        .collect();
    // El fin más lejano de los tramos buscados de este periodo hasta ahora, y
    // el periodo al que pertenece. Solo lo mueve un tramo conservado: un tramo
    // unido descartado aquí deja buscables a sus trozos.
    let mut reach = (0u8, 0u64);
    // Los supervivientes se compactan en su sitio: `retain` visita en orden y
    // una vez cada uno, que es todo lo que necesita el pase.
    out.retain(|&iv| {
        let nested = reach.0 == iv.period && iv.end <= reach.1;
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
            let first = candidates.partition_point(|a| a.1 <= iv.start);
            for a in &candidates[first..] {
                if a.0 >= iv.end {
                    break;
                }
                // `a` debe cubrir la mitad de `iv`; uno más corto que eso no
                // puede, se coloque como se coloque.
                if (a.1 - a.0) * 2 < iv_len {
                    continue;
                }
                let overlap = a.1.min(iv.end).saturating_sub(a.0.max(iv.start));
                if overlap * 2 >= iv_len {
                    redundant = true;
                    break;
                }
            }
            if redundant {
                break;
            }
        }
        if !redundant && !nested {
            reach = (iv.period, iv.end);
            if let Some(list) = by_period.get_mut(usize::from(iv.period)) {
                list.push((iv.start, iv.end));
            }
        }
        !redundant
    });
    out
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
    /// No cabe **entero** en ninguna de las regiones declaradas (DP-SSR-18).
    /// Se aplica antes de los solapes: lo que queda fuera no puede quitarle
    /// el sitio a lo de dentro.
    OutsideRegions,
    /// Solapa con otro candidato que gana por la regla de DP-SSR-05 y no se
    /// salva recortándose, o es un recortado que ya no supera los criterios.
    Overlap,
}

impl DetectionStage {
    /// Todas las etapas, en orden de aplicación.
    pub const ALL: [Self; 8] = [
        Self::NonCanonicalMotif,
        Self::MinCopies,
        Self::MinPurity,
        Self::MinLongestRun,
        Self::MaxUnits,
        Self::Edge,
        Self::OutsideRegions,
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
            Self::OutsideRegions => "outside_regions",
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
/// detector, sin recortar (DP-SSR-04, opción (a)), salvo el recorte por solape
/// de una copia o menos (DP-SSR-05).
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
    detect_repeat_candidates_in_regions(bases, criteria, None)
}

/// Como [`detect_repeat_candidates`], pero conservando solo los tramos que
/// caben **enteros** en alguna de las `regions` (DP-SSR-18).
///
/// Las regiones van en base cero y semiabiertas sobre la misma secuencia,
/// ordenadas por inicio y sin solaparse entre sí (quien las resuelve lo
/// garantiza). `None` es la secuencia entera. El barrido es siempre el de la
/// secuencia entera, de modo que las fronteras de un tramo no dependen de
/// dónde corte una región: la región solo decide qué tramos se quedan. Un
/// tramo que se queda fuera se cuenta en [`DetectionStage::OutsideRegions`],
/// que va **antes** de los solapes.
#[must_use]
pub fn detect_repeat_candidates_in_regions(
    bases: &[u8],
    criteria: &DetectionCriteria,
    regions: Option<&[(u64, u64)]>,
) -> SequenceDetection {
    let intervals = scan_all_periods(bases);
    let mut counts = DetectionCounts::default();
    let mut survivors = Vec::new();
    for interval in &intervals {
        let motif_length = usize::from(interval.period);
        counts.detected[motif_length - 1] += 1;
        match classify_interval(bases, interval, criteria) {
            Ok(candidate) if !inside_regions(regions, candidate.start, candidate.end) => {
                counts.discard(DetectionStage::OutsideRegions, motif_length);
            }
            Ok(candidate) => survivors.push(candidate),
            Err(stage) => counts.discard(stage, motif_length),
        }
    }
    let (candidates, overlapped) = resolve_overlaps(bases, criteria, survivors);
    for motif_length in overlapped {
        counts.discard(DetectionStage::Overlap, motif_length);
    }
    for candidate in &candidates {
        counts.accepted[candidate.motif.length() as usize - 1] += 1;
    }
    SequenceDetection { candidates, counts }
}

/// Los intervalos del detector con su puntuación fija (2 y 7) y periodos de
/// 1 a 6.
fn scan_all_periods(bases: &[u8]) -> Vec<RepeatInterval> {
    let periods = PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1)
        .unwrap_or_else(|_| unreachable!("1..=6 es un intervalo válido"));
    find_tandem_repeats(bases, periods, &ScanParams::default())
}

/// Si `[start, end)` cabe entero en alguna región; sin regiones, siempre.
fn inside_regions(regions: Option<&[(u64, u64)]>, start: u64, end: u64) -> bool {
    let Some(regions) = regions else {
        return true;
    };
    // Las regiones no se solapan y van por inicio: la única que puede
    // contener el tramo es la última que empieza en su inicio o antes.
    let index = regions.partition_point(|(region_start, _)| *region_start <= start);
    index > 0 && regions[index - 1].1 >= end
}

/// Si `candidate` es una rotación de `motif`: misma longitud y aparece dentro
/// de `motif` repetido dos veces. `CA` es rotación de `AC`, y de sí mismo.
#[must_use]
pub fn is_rotation(candidate: &[u8], motif: &[u8]) -> bool {
    candidate.len() == motif.len()
        && !motif.is_empty()
        && [motif, motif]
            .concat()
            .windows(candidate.len())
            .any(|window| window == candidate)
}

/// Por qué un tramo de la secuencia que se solapa con `[start, end)` y cuyo
/// motivo es una rotación de alguno de `motifs` **no** quedó entre los
/// candidatos: la etapa en la que cayó y sus fronteras, del tramo más largo
/// de los que lo solapan. `None` si el detector no ve ninguno así.
///
/// Solo la usa el diagnóstico de un locus declarado que no coincide con
/// ningún detectado (REQ-REF-40): vuelve a barrer la secuencia entera, que es
/// caro, pero ocurre una vez y justo antes de rechazar.
#[must_use]
pub fn explain_missing_tract(
    bases: &[u8],
    criteria: &DetectionCriteria,
    regions: Option<&[(u64, u64)]>,
    start: u64,
    end: u64,
    motifs: &[RepeatMotif],
) -> Option<MissingTract> {
    scan_all_periods(bases)
        .into_iter()
        .filter(|interval| interval.start < end && start < interval.end)
        .filter(|interval| {
            let period = usize::from(interval.period);
            let first = &bases[interval.start as usize..interval.start as usize + period];
            motifs.iter().any(|motif| is_rotation(first, motif.bases()))
        })
        .max_by_key(|interval| {
            (
                interval.end - interval.start,
                std::cmp::Reverse(interval.start),
            )
        })
        .map(|interval| {
            let stage = match classify_interval(bases, &interval, criteria) {
                Err(stage) => stage,
                Ok(candidate) if !inside_regions(regions, candidate.start, candidate.end) => {
                    DetectionStage::OutsideRegions
                }
                // Pasó todo lo que se juzga solo: si no está entre los
                // candidatos, lo quitó un solape.
                Ok(_) => DetectionStage::Overlap,
            };
            MissingTract {
                start: interval.start,
                end: interval.end,
                period: interval.period,
                stage,
            }
        })
}

/// Un tramo del detector que no llegó a candidato, con la etapa que no superó.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct MissingTract {
    /// Base cero, semiabierto.
    pub start: u64,
    pub end: u64,
    pub period: u8,
    pub stage: DetectionStage,
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

/// Resuelve los solapes con la regla de DP-SSR-05: ordena los candidatos por
/// prioridad —mayor puntuación del detector; a igual puntuación, más bases;
/// luego periodo menor; luego el de la izquierda— y los recorre en ese orden.
/// De cada aspirante:
///
/// - **sin solape** con los ya aceptados, se acepta;
/// - con solape de **una copia de su motivo o menos por cada lado** (solape ≤
///   periodo del aspirante) y sin ningún aceptado dentro de él ni que lo cubra
///   entero, **se recorta** a la parte que no solapa, y ese tramo recortado
///   **vuelve a pasar por todos los criterios** ([`classify_interval`]): entra
///   solo si los supera;
/// - en cualquier otro caso (algún lado con más de una copia de solape, un
///   aceptado dentro de él o que lo cubre) **se descarta entero**.
///
/// Un microsatélite suele acabar con una copia parcial del motivo que coincide
/// con el principio de lo que viene detrás: `(CA)×10 A×12` da `CA` y `A`
/// pisándose una base, y el recorte conserva los dos.
///
/// **El recortado cae, si cae, en «solape»** (devuelve su longitud de motivo
/// entre los descartados), y no en la etapa del criterio que ya no supera:
/// superó todos por sí solo, y lo pierde el solape. Así las columnas de suelo,
/// pureza y demás siguen contando solo lo que cae por sí solo.
///
/// **La puntuación del recortado se recalcula** con las comparaciones de su
/// propio tramo, como el recorte del inicio del buscador: es lo que todo
/// candidato lleva, recortado o no. No decide nada aquí, porque el recortado
/// se acepta en el acto y no vuelve a competir.
///
/// **El recortado no necesita comprobar las regiones**: es un subconjunto del
/// aspirante, que ya cabía en la suya (las regiones se aplican antes).
///
/// No es «de cada grupo solapado gana uno»: en una cadena A–B–C donde A gana
/// a B y C solo solapa a B, B se recorta o se descarta y C compite contra lo
/// aceptado. El resultado **no depende del orden de entrada**: la clave
/// ordena del todo (dos candidatos con la misma clave tendrían el mismo
/// intervalo y periodo, y el detector no los emite). Dos tramos que solo se
/// tocan no se solapan. Devuelve (aceptados ordenados por inicio, longitudes de
/// motivo de los descartados).
fn resolve_overlaps(
    bases: &[u8],
    criteria: &DetectionCriteria,
    mut survivors: Vec<RepeatCandidate>,
) -> (Vec<RepeatCandidate>, Vec<usize>) {
    survivors.sort_by_key(RepeatCandidate::priority);
    // Aceptados por inicio. Son disjuntos, así que también están ordenados por
    // fin y los que solapan al aspirante son los últimos que empiezan antes de
    // su fin.
    let mut accepted: BTreeMap<u64, RepeatCandidate> = BTreeMap::new();
    let mut discarded = Vec::new();
    for candidate in survivors {
        let (start, end) = (candidate.start, candidate.end);
        let period = candidate.motif.length();
        // Solape por la izquierda y por la derecha; `whole` si algún aceptado
        // está dentro del aspirante o lo cubre entero (no hay nada que recortar).
        let (mut left, mut right, mut whole) = (0, 0, false);
        for (_, previous) in accepted.range(..end).rev() {
            if previous.end <= start {
                break;
            }
            if previous.start <= start && previous.end < end {
                left = previous.end - start;
            } else if previous.start > start && previous.end >= end {
                right = end - previous.start;
            } else {
                whole = true;
                break;
            }
        }
        if whole || left > period || right > period {
            discarded.push(period as usize);
            continue;
        }
        if left == 0 && right == 0 {
            accepted.insert(start, candidate);
            continue;
        }
        let (trimmed_start, trimmed_end) = (start + left, end - right);
        // `classify_interval` lee un motivo entero: un resto más corto que él
        // no es un tramo.
        if trimmed_end - trimmed_start < period {
            discarded.push(period as usize);
            continue;
        }
        let trimmed = RepeatInterval {
            start: trimmed_start,
            end: trimmed_end,
            period: period as u8,
            score: own_score(bases, trimmed_start, trimmed_end, period),
        };
        match classify_interval(bases, &trimmed, criteria) {
            Ok(rejudged) => {
                accepted.insert(rejudged.start, rejudged);
            }
            Err(_) => discarded.push(period as usize),
        }
    }
    (accepted.into_values().collect(), discarded)
}

/// Puntuación de un tramo con las comparaciones **de su propio tramo**: la
/// recompensa por cada base igual a la que está `period` antes dentro del
/// tramo, y la penalización por cada una distinta o no ACGT. Es la que el
/// buscador da a un tramo sin recortar.
fn own_score(bases: &[u8], start: u64, end: u64, period: u64) -> i32 {
    let params = ScanParams::default();
    let (start, end, period) = (start as usize, end as usize, period as usize);
    (start + period..end)
        .map(|j| {
            let (here, back) = (
                CANONICAL[bases[j] as usize],
                CANONICAL[bases[j - period] as usize],
            );
            if here != 0 && here == back {
                params.match_reward
            } else {
                -params.mismatch_penalty
            }
        })
        .sum()
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

    /// Lo que conserva quien aplica un suelo de copias es lo mismo con el
    /// barrido que conoce el suelo que con el barrido entero: sobre secuencia
    /// densa en repeticiones, donde las eliminaciones por periodo múltiplo son
    /// frecuentes, y con suelos de permisivos a estrictos.
    #[test]
    fn the_floor_aware_scan_keeps_what_the_floor_keeps() {
        let mut state: u64 = 0x2545F4914F6CDD1D;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (state >> 33) as u32
        };
        let floors: [[u32; 6]; 5] = [
            [5, 5, 4, 4, 4, 3],
            [2, 2, 2, 2, 2, 2],
            [3, 7, 2, 5, 2, 2],
            [12, 2, 6, 2, 3, 2],
            [9, 9, 9, 9, 9, 9],
        ];
        let mut sub_floor_kept = 0usize;
        for _ in 0..300 {
            let mut seq = Vec::new();
            while seq.len() < 3_000 {
                let period = 1 + (next() % 6) as usize;
                let motif: Vec<u8> = (0..period)
                    .map(|_| b"ACGT"[(next() % 4) as usize])
                    .collect();
                for _ in 0..1 + next() % 10 {
                    for &b in &motif {
                        seq.push(if next() % 12 == 0 {
                            b"ACGT"[(next() % 4) as usize]
                        } else {
                            b
                        });
                    }
                }
            }
            let whole = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
            for floor in &floors {
                let at = |period: u8| floor[usize::from(period) - 1];
                let clears = |iv: &&RepeatInterval| {
                    (iv.end - iv.start) / u64::from(iv.period) >= u64::from(at(iv.period))
                };
                let want: Vec<_> = whole.iter().filter(clears).collect();
                let got_all =
                    find_tandem_repeats_for_copy_floor(&seq, p(1, 6), &ScanParams::default(), at);
                let got: Vec<_> = got_all.iter().filter(clears).collect();
                assert_eq!(got, want, "floor {floor:?}");
                sub_floor_kept += got_all.len().saturating_sub(got.len());
            }
        }
        assert!(
            sub_floor_kept > 0,
            "the sequences exercise sub-floor tracts"
        );
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

    // ---- Corte por cambio de motivo (copiadas de fa215c23) -----------------

    /// Dos repeticiones pegadas de motivos distintos salen como **dos**
    /// intervalos, cada uno con su motivo, y no como uno con el primer motivo
    /// y la mitad del segundo, que la pureza rechazaría entero. Las cuatro
    /// formas de dos repeticiones que un simulador encontró fundidas.
    #[test]
    fn abutting_repeats_of_different_motifs_are_split() {
        /// Una secuencia y los intervalos `(start, end, period)` que debe dar.
        type Case = (Vec<u8>, Vec<(u64, u64, u8)>);
        let rep = |unit: &[u8], n: usize| unit.repeat(n);
        let cases: Vec<Case> = vec![
            (
                [rep(b"A", 12), rep(b"T", 12)].concat(),
                vec![(0, 12, 1), (12, 24, 1)],
            ),
            (
                [rep(b"G", 60), rep(b"A", 15)].concat(),
                vec![(0, 60, 1), (60, 75, 1)],
            ),
            (
                [rep(b"GA", 12), rep(b"C", 12)].concat(),
                vec![(0, 24, 2), (24, 36, 1)],
            ),
            (
                [rep(b"CA", 10), rep(b"GT", 10)].concat(),
                vec![(0, 20, 2), (20, 40, 2)],
            ),
        ];
        for (seq, want) in cases {
            let got = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
            let mut spans: Vec<(u64, u64, u8)> =
                got.iter().map(|r| (r.start, r.end, r.period)).collect();
            spans.sort_unstable();
            for w in &want {
                assert!(
                    spans.contains(w),
                    "{} lacks {w:?}: {spans:?}",
                    String::from_utf8_lossy(&seq)
                );
            }
            // Nada cruza la unión en el periodo del primer motivo.
            let junction = want[0].1;
            assert!(
                !got.iter().any(|r| r.start < junction - 1
                    && r.end > junction + 1
                    && r.period == want[0].2),
                "{} still has a merged interval: {spans:?}",
                String::from_utf8_lossy(&seq)
            );
        }
    }

    /// Los intervalos de `got` como `(start, end)`, tras comprobar que el
    /// primero es `joined` y que todos los demás caen dentro de él: la forma
    /// de un tramo unido emitido con sus trozos.
    fn joined_with_pieces(got: &[RepeatInterval], joined: (u64, u64)) -> Vec<(u64, u64)> {
        let spans: Vec<(u64, u64)> = got.iter().map(|r| (r.start, r.end)).collect();
        assert_eq!(spans.first(), Some(&joined), "{got:?}");
        assert!(spans.len() > 1, "the pieces come with it: {got:?}");
        assert!(
            spans[1..]
                .iter()
                .all(|&(s, e)| s >= joined.0 && e <= joined.1),
            "every piece lies inside the joined tract: {got:?}"
        );
        spans
    }

    /// Dos copias de otro motivo **dentro** de una repetición son una
    /// interrupción, no una segunda repetición: `(AT)×4 A (CA)×2 (TA)×4 T` sale
    /// como un intervalo `AT` sobre todo él, como antes de que existiera el
    /// corte, y con los trozos de los que se unió, para que quien lo lea pueda
    /// quedarse con ellos si el tramo unido no es lo bastante puro. Coincide
    /// con `AT` en 18 de 20 bases, por encima del 7/9 que exige la puntuación
    /// por omisión.
    #[test]
    fn a_short_run_of_another_motif_inside_a_repeat_is_not_split() {
        let seq = b"ATATATATACACATATATAT";
        let got = find_tandem_repeats(seq, p(2, 2), &ScanParams::default());
        assert_eq!(
            joined_with_pieces(&got, (0, 20)),
            vec![(0, 20), (0, 9), (8, 13), (12, 20)]
        );
    }

    /// Dos interrupciones **distintas** seguidas también se absorben cuando
    /// vuelve el motivo exterior: `(AAT)×4 (GAT)×3 (GGT)×2 (AAT)×4` coincide
    /// con `AAT` en 32 de 39 bases, 0,82. Una prueba que solo mirara a los
    /// vecinos no lo veía: ningún trozo tiene ahí el mismo motivo a los dos
    /// lados.
    #[test]
    fn two_different_interruptions_in_a_row_are_absorbed() {
        let seq = [
            b"AAT".repeat(4),
            b"GAT".repeat(3),
            b"GGT".repeat(2),
            b"AAT".repeat(4),
        ]
        .concat();
        let got = find_tandem_repeats(&seq, p(3, 3), &ScanParams::default());
        joined_with_pieces(&got, (0, 39));
    }

    /// La misma interrupción, demasiado larga para la tolerancia:
    /// `(CA)×10 (GT)×10 (CA)×10` tiene `CA` en dos tercios de sus bases, por
    /// debajo de 7/9, así que salen las tres repeticiones.
    #[test]
    fn a_long_run_of_another_motif_inside_a_repeat_is_split() {
        let seq = [b"CA".repeat(10), b"GT".repeat(10), b"CA".repeat(10)].concat();
        let got = find_tandem_repeats(&seq, p(2, 2), &ScanParams::default());
        assert_eq!(
            got.iter().map(|r| (r.start, r.end)).collect::<Vec<_>>(),
            vec![(0, 20), (20, 40), (40, 60)],
            "{got:?}"
        );
    }

    /// Un homopolímero solo no lo toca el corte.
    #[test]
    fn a_lone_homopolymer_is_one_interval() {
        let seq = b"AAAAAAAAAAAA";
        let got = find_tandem_repeats(seq, p(1, 6), &ScanParams::default());
        assert_eq!(
            got.iter()
                .map(|r| (r.start, r.end, r.period))
                .collect::<Vec<_>>(),
            vec![(0, 12, 1)]
        );
    }

    // ---- Desviación: inicio avanzado hasta la primera fase primitiva ------

    /// El intervalo de periodo `p` de `seq` que contiene `[from, to)`.
    fn covering(seq: &[u8], period: u8, from: u64, to: u64) -> Option<RepeatInterval> {
        find_tandem_repeats(seq, p(1, 6), &ScanParams::default())
            .into_iter()
            .find(|r| r.period == period && r.start <= from && to <= r.end)
    }

    /// Puntuación calculada aparte, desde la definición: +2 por cada base del
    /// tramo igual a la que está `p` antes **dentro del tramo**, −7 por cada
    /// una distinta.
    fn own_score(seq: &[u8], interval: &RepeatInterval) -> i32 {
        let p = usize::from(interval.period);
        (interval.start as usize + p..interval.end as usize)
            .map(|j| {
                let (a, b) = (seq[j].to_ascii_uppercase(), seq[j - p].to_ascii_uppercase());
                if a == b && b"ACGT".contains(&a) {
                    2
                } else {
                    -7
                }
            })
            .sum()
    }

    /// Inicio y fin relativos al núcleo, motivo y puntuación.
    type CoreCandidate = (u64, u64, Vec<u8>, i32);

    /// Pasa `core` por `detect_repeat_candidates` entre dos fondos sin
    /// repeticiones y unas bases de flanco elegidas para no alargar nada.
    /// Devuelve los candidatos con coordenadas relativas a `core` y el
    /// recuento.
    fn detect_core(
        left: &[u8],
        core: &[u8],
        right: &[u8],
    ) -> (Vec<CoreCandidate>, DetectionCounts) {
        let mut sequence = background(0xD1, 30);
        sequence.extend_from_slice(left);
        let offset = sequence.len() as u64;
        sequence.extend_from_slice(core);
        sequence.extend_from_slice(right);
        sequence.extend(background(0xD2, 30));
        let detection = detect_repeat_candidates(&sequence, &DetectionCriteria::default());
        let candidates = detection
            .candidates()
            .iter()
            .map(|c| {
                (
                    c.start() - offset,
                    c.end() - offset,
                    c.motif().bases().to_vec(),
                    c.score(),
                )
            })
            .collect();
        (candidates, detection.counts().clone())
    }

    /// **Regresión: `A×8 (GA)×10`.** Con el núcleo original solo salía el
    /// poli-A. El segmento de periodo 2 empieza dentro del poli-A (motivo
    /// `AA`, no primitivo); el inicio avanza a la primera fase primitiva **fuera
    /// de la tira de `A`**: `[8, 28)`, motivo `GA`, puntuación 36 (18
    /// coincidencias propias, no las del segmento con el poli-A dentro). El
    /// poli-A `[0, 8)` y el `(GA)n` solo se tocan, y `detect_repeat_candidates`
    /// acepta los dos sin ningún descarte por solape.
    #[test]
    fn a_microsatellite_after_a_long_homopolymer_of_its_base_is_found() {
        let seq = [&b"AAAAAAAA"[..], &b"GA".repeat(10)].concat();
        let got = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
        assert!(got.contains(&RepeatInterval {
            start: 0,
            end: 8,
            period: 1,
            score: 14
        }));
        assert!(got.contains(&RepeatInterval {
            start: 8,
            end: 28,
            period: 2,
            score: 36
        }));
        let (candidates, counts) = detect_core(b"C", &seq, b"C");
        assert_eq!(
            candidates,
            vec![(0, 8, b"A".to_vec(), 14), (8, 28, b"GA".to_vec(), 36)]
        );
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
    }

    /// **Regresión: `T×12 (CTT)×7`.** Con el núcleo original solo salía el
    /// poli-T. Ahora el de periodo 3 empieza en 12, motivo `CTT` (la primera
    /// fase primitiva, `TTC` en 10, está dentro de la tira de `T`), y llega al
    /// final: `[12, 33)`, puntuación 36. Los dos se aceptan, sin solape. Y la
    /// re-detección de periodo 6, que también empieza dentro del poli-T y se
    /// recorta, la sigue quitando el descarte geométrico.
    #[test]
    fn a_trinucleotide_after_a_homopolymer_of_its_base_is_found() {
        let seq = [&b"TTTTTTTTTTTT"[..], &b"CTT".repeat(7)].concat();
        let got = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
        assert!(got.contains(&RepeatInterval {
            start: 0,
            end: 12,
            period: 1,
            score: 22
        }));
        assert!(got.contains(&RepeatInterval {
            start: 12,
            end: 33,
            period: 3,
            score: 36
        }));
        assert!(
            got.iter().all(|r| r.period != 2 && r.period != 6),
            "{got:?}"
        );
        let (candidates, counts) = detect_core(b"G", &seq, b"G");
        assert_eq!(
            candidates,
            vec![(0, 12, b"T".to_vec(), 22), (12, 33, b"CTT".to_vec(), 36)]
        );
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
    }

    /// **Regresión: `A×6 (GA)×10 C`.** El `(GA)n` sale desde 6: `[6, 26)`,
    /// motivo `GA`, puntuación 36. El poli-A de 6 también sale del detector,
    /// pero no alcanza su suelo de 8 copias y se descarta en esa etapa: el
    /// único candidato es el `(GA)n`.
    #[test]
    fn a_microsatellite_after_a_short_homopolymer_of_its_base_is_found() {
        let seq = [&b"AAAAAA"[..], &b"GA".repeat(10), b"C"].concat();
        assert_eq!(
            covering(&seq, 2, 6, 26),
            Some(RepeatInterval {
                start: 6,
                end: 26,
                period: 2,
                score: 36
            })
        );
        let (candidates, counts) = detect_core(b"C", &seq, b"T");
        assert_eq!(candidates, vec![(6, 26, b"GA".to_vec(), 36)]);
        assert!(counts.discarded(DetectionStage::MinCopies, 1) >= 1);
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
    }

    /// **El caso que ya funcionaba, `A×4 (GA)×10`, no cambia:** el poli-A es
    /// demasiado corto para que el segmento de periodo 2 empiece dentro de él
    /// con un motivo no primitivo; el segmento empieza ya en fase `AG` en la
    /// última `A`, no hay recorte y no se aplica ninguna de las dos reglas.
    #[test]
    fn a_tract_whose_start_is_already_primitive_is_unchanged() {
        let seq = [&b"AAAA"[..], &b"GA".repeat(10)].concat();
        assert_eq!(
            covering(&seq, 2, 4, 24),
            Some(RepeatInterval {
                start: 3,
                end: 24,
                period: 2,
                score: 38
            })
        );
    }

    /// **El recorte respeta el mínimo de copias del detector**, contado sobre
    /// el tramo final y no sobre el segmento: con un mínimo de 10 copias,
    /// `A×8 (GA)×10` deja `[8, 28)`, 20 bases (10 copias), y sale; `A×8
    /// (GA)×9` deja 18 (9 copias) y no sale, aunque el segmento entero tenga 26
    /// bases (13 copias).
    #[test]
    fn trimming_keeps_the_detector_minimum_copies() {
        let params = ScanParams {
            min_copies: 10,
            ..ScanParams::default()
        };
        let ten = [&b"AAAAAAAA"[..], &b"GA".repeat(10)].concat();
        let nine = [&b"AAAAAAAA"[..], &b"GA".repeat(9)].concat();
        let period_two =
            |seq: &[u8]| -> Vec<RepeatInterval> { find_tandem_repeats(seq, p(2, 2), &params) };
        assert_eq!(
            period_two(&ten),
            vec![RepeatInterval {
                start: 8,
                end: 28,
                period: 2,
                score: 36
            }]
        );
        assert_eq!(period_two(&nine), vec![]);
    }

    /// **La puntuación de todo intervalo es la de sus propias comparaciones**,
    /// recortado o no, contra un cálculo hecho aparte. Sobre una secuencia con
    /// microsatélites de periodos 1 a 6 pegados a homopolímeros de su base, a
    /// homopolímeros de otra base y con interrupciones.
    #[test]
    fn every_interval_score_is_the_score_of_its_own_comparisons() {
        let mut seq = background(0x51, 30);
        let blocks: [(&[u8], usize, &[u8]); 8] = [
            (b"AAAAAAAA", 10, b"GA"),
            (b"TTTTTTTTTTTT", 7, b"CTT"),
            (b"CCCCCCCCC", 6, b"ACCG"),
            (b"GGGGGGGGGG", 5, b"GGATC"),
            (b"AAAAAAAAA", 4, b"AATCGA"),
            (b"CCCC", 9, b"CA"),
            (b"TTTTTTTT", 8, b"TA"),
            (b"", 6, b"CAGCAT"),
        ];
        let mut trimmed = 0;
        for (index, (homopolymer, copies, motif)) in blocks.iter().enumerate() {
            seq.extend_from_slice(homopolymer);
            seq.extend(motif.repeat(*copies));
            seq.extend(background(0x60 + index as u64, 30));
        }
        let got = find_tandem_repeats(&seq, p(1, 6), &ScanParams::default());
        for interval in &got {
            assert_eq!(interval.score, own_score(&seq, interval), "{interval:?}");
            let period = usize::from(interval.period);
            let start = interval.start as usize;
            assert!(is_primitive_motif(&seq[start..start + period]));
            // Probablemente recortado: de periodo 2 o más y justo detrás de
            // una tira de al menos cuatro bases iguales. Solo sirve para
            // comprobar que la secuencia ejercita el caso; no distingue un
            // tramo recortado de uno que ya empezaba ahí (`CCCC (CA)×9` empieza
            // sin recorte en la primera `C` del dinucleótido).
            if interval.period > 1
                && start >= 4
                && seq[start - 4..start]
                    .iter()
                    .all(|&base| base == seq[start - 1])
            {
                trimmed += 1;
            }
        }
        assert!(
            trimmed >= 3,
            "la secuencia debe ejercitar el recorte: {got:?}"
        );
    }

    // ---- Regresiones: tramos pegados de motivos distintos (propias) --------
    //
    // Los seis casos que el núcleo de `950909f` fundía en un solo tramo. Lo
    // esperado está escrito a mano: inicio y fin relativos al núcleo, motivo y
    // puntuación (+2 por coincidencia propia del tramo, −7 por discrepancia).
    // Todas, salvo la de control, fallan con aquel núcleo.

    /// Atajo: los candidatos de `core` entre dos flancos, sin el recuento.
    fn candidates_of(left: &[u8], core: &[u8], right: &[u8]) -> Vec<CoreCandidate> {
        detect_core(left, core, right).0
    }

    /// **Control: `A×12` sola** sale entera, `[0, 12)`, con 11 coincidencias.
    /// Ya salía bien con el núcleo viejo; está para que las demás comparen
    /// contra una tira que se sabe detectar.
    #[test]
    fn a_lone_homopolymer_is_a_candidate() {
        let core = b"A".repeat(12);
        assert_eq!(
            candidates_of(b"C", &core, b"C"),
            vec![(0, 12, b"A".to_vec(), 22)]
        );
    }

    /// **Regresión: `A×12 T×12`.** El núcleo viejo daba un solo tramo de
    /// periodo 1 `[0, 24)` con motivo `A` y la mitad de `T`, que no pasaba la
    /// pureza: ningún candidato. Ahora salen las dos tiras, que solo se tocan.
    #[test]
    fn two_abutting_homopolymers_are_two_candidates() {
        let core = [b"A".repeat(12), b"T".repeat(12)].concat();
        assert_eq!(
            candidates_of(b"C", &core, b"G"),
            vec![(0, 12, b"A".to_vec(), 22), (12, 24, b"T".to_vec(), 22)]
        );
    }

    /// **Regresión: `G×60 A×15`.** El núcleo viejo daba un solo candidato
    /// `[0, 75)` con motivo `G` (pureza 60/75 = 0,8, justo en el mínimo) y la
    /// tira de `A` desaparecía dentro de él. Ahora, dos candidatos que solo
    /// se tocan.
    #[test]
    fn a_long_homopolymer_and_a_short_one_are_two_candidates() {
        let core = [b"G".repeat(60), b"A".repeat(15)].concat();
        assert_eq!(
            candidates_of(b"C", &core, b"C"),
            vec![(0, 60, b"G".to_vec(), 118), (60, 75, b"A".to_vec(), 28)]
        );
    }

    /// **Regresión: `(GA)×12 C×12`.** El núcleo viejo daba un tramo de periodo
    /// 2 `[0, 36)` de pureza 24/36 y solo quedaba la tira de `C`. Ahora el
    /// `(GA)n` sale con sus 22 coincidencias y toca a la tira.
    #[test]
    fn a_dinucleotide_before_a_homopolymer_is_kept() {
        let core = [b"GA".repeat(12), b"C".repeat(12)].concat();
        assert_eq!(
            candidates_of(b"T", &core, b"A"),
            vec![(0, 24, b"GA".to_vec(), 44), (24, 36, b"C".to_vec(), 22)]
        );
    }

    /// **Regresión: `(CA)×10 (GT)×10`.** El núcleo viejo daba un tramo de
    /// periodo 2 `[0, 40)` con motivo `CA` de pureza 0,5: ningún candidato.
    /// Ahora, los dos dinucleótidos.
    #[test]
    fn two_abutting_dinucleotides_are_two_candidates() {
        let core = [b"CA".repeat(10), b"GT".repeat(10)].concat();
        assert_eq!(
            candidates_of(b"T", &core, b"C"),
            vec![(0, 20, b"CA".to_vec(), 36), (20, 40, b"GT".to_vec(), 36)]
        );
    }

    /// **Regresión: el compuesto `(TCTA)×8 (TCTG)×8 (TCTA)×8`.** El núcleo
    /// viejo daba un solo tramo de periodo 4 `[0, 96)`, de pureza 64/96 con
    /// `TCTA`, y no quedaba nada. El nuevo emite ese tramo unido (el corte lo
    /// une porque coincide con `TCTA` en fase en 88 de 96 bases, por encima de
    /// 7/9) **y** sus tres trozos: `[0, 35)` `TCTA`, `[32, 67)` `TCTG` y
    /// `[64, 96)` `TCTA`. El unido cae en la pureza; los trozos pasan solos.
    ///
    /// Cada trozo acaba con la copia parcial `TCT` que comparte con el motivo
    /// siguiente, así que el del centro solapa 3 bases con cada vecino, menos
    /// de una copia (4). Antes la regla de solapes descartaba el `(TCTG)n`
    /// entero; ahora **se recorta** y se queda. El primero (empate a 62 y 35
    /// bases: gana el de la izquierda) se acepta tal cual; el del centro pierde
    /// por la izquierda 3 bases y queda `[35, 67)`, de **motivo `GTCT`** (el
    /// motivo se lee al inicio del tramo, una rotación de `TCTG`) y puntuación
    /// 56 (28 comparaciones propias); el tercero ya solapa 3 bases con el
    /// centro recortado y queda `[67, 96)`, de motivo `ATCT` y puntuación 50.
    /// Los tres se vuelven a juzgar y superan los criterios; no hay descartes
    /// por solape.
    ///
    /// **Debe fallar contra** la regla de antes (descartar entero): daría solo
    /// los dos tramos de los extremos, `[0, 35)` y `[64, 96)`.
    #[test]
    fn a_compound_tract_yields_all_three_pieces() {
        let core = [b"TCTA".repeat(8), b"TCTG".repeat(8), b"TCTA".repeat(8)].concat();
        let (candidates, counts) = detect_core(b"G", &core, b"G");
        assert_eq!(
            candidates,
            vec![
                (0, 35, b"TCTA".to_vec(), 62),
                (35, 67, b"GTCT".to_vec(), 56),
                (67, 96, b"ATCT".to_vec(), 50)
            ]
        );
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!(counts.discarded(DetectionStage::MinPurity, 4), 1);
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
        assert_eq!(counts.accepted(4), 3);
    }

    /// **Regresión: `(CA)×10 A×12`.** El buscador da `CA [0, 20)` y la tira de
    /// `A [19, 32)` (13 bases), que se pisan **una base** (la última `A` del
    /// dinucleótido). Una base es una copia del motivo de la tira, justo el
    /// límite que el responsable fijó el 01/11/2026 (solape ≤ periodo): la
    /// tira se recorta a `[20, 32)`, 12 bases, se vuelve a juzgar (suelo de 8
    /// copias) y se queda, con la puntuación de su propio tramo (11
    /// comparaciones, 22). Cada una cuenta en su longitud de motivo y no hay
    /// descarte por solape.
    ///
    /// **Debe fallar contra** la regla de antes: solo saldría el `(CA)×10`.
    #[test]
    fn overlap_trim_keeps_a_homopolymer_sharing_one_base_with_a_dinucleotide() {
        let core = [b"CA".repeat(10), b"A".repeat(12)].concat();
        let (candidates, counts) = detect_core(b"G", &core, b"G");
        assert_eq!(
            candidates,
            vec![(0, 20, b"CA".to_vec(), 36), (20, 32, b"A".to_vec(), 22)]
        );
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!((counts.accepted(1), counts.accepted(2)), (1, 1));
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
    }

    /// **`(CA)×10 G (CAG)×7` da los dos.** Es el caso que la regla de antes
    /// resolvía descartando el `(CA)×10`: el `(CAG)×8` empieza en el último
    /// `CA` (`CAG` con la `G`), `[18, 42)` con puntuación 42, y gana al
    /// `(CA)×10` `[0, 20)` (36). El `(CA)×10` lo pisa **dos bases**, una copia
    /// de su motivo: se recorta a `[0, 18)` (`(CA)×9`, 16 comparaciones, 32) y
    /// supera el suelo de 6 copias. Son dos microsatélites reales pegados.
    ///
    /// **Debe fallar contra** la regla de antes: solo saldría el `(CAG)×8`.
    #[test]
    fn overlap_trim_keeps_a_dinucleotide_overlapped_by_one_copy_by_a_trinucleotide() {
        let core = [b"CA".repeat(10), b"G".to_vec(), b"CAG".repeat(7)].concat();
        let (candidates, counts) = detect_core(b"T", &core, b"T");
        assert_eq!(
            candidates,
            vec![(0, 18, b"CA".to_vec(), 32), (18, 42, b"CAG".to_vec(), 42)]
        );
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
    }

    /// **Un solape de más de una copia del que pierde se sigue descartando
    /// entero.** `(GC)×8 (ACGC)×6`: el `(ACGC)` entra tres bases en el
    /// dinucleótido (`[13, 16)`), más que sus dos bases de periodo, y gana por
    /// puntuación; el `(GC)×8` se descarta y se cuenta en «solape», sin
    /// recortarlo a `[0, 13)` (seis copias completas de `GC`), que pasaría el
    /// suelo.
    ///
    /// **Debe fallar contra** recortar con cualquier solape: saldría también
    /// ese `GC` de seis copias.
    #[test]
    fn an_overlap_of_more_than_one_copy_still_discards_the_loser_whole() {
        let core = [b"GC".repeat(8), b"ACGC".repeat(6)].concat();
        let (candidates, counts) = detect_core(b"T", &core, b"T");
        assert_eq!(candidates, vec![(13, 40, b"CGCA".to_vec(), 46)]);
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!(counts.discarded(DetectionStage::Overlap, 2), 1);
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 1);
    }

    /// **El límite exacto, por `resolve_overlaps` y con las bases a la vista.**
    /// Un `(GAT)×8` con un ganador escrito a mano que lo pisa por la
    /// izquierda: con 3 bases (una copia) se recorta y se queda en `[33, 54)`;
    /// con 4 (una copia y una base) se descarta entero. Es lo que distingue
    /// `≤` de `<` en el umbral, por cada lado: la prueba repite el caso por la
    /// derecha, con el ganador a la derecha.
    #[test]
    fn the_trim_threshold_is_exactly_one_copy_on_each_side() {
        let mut bases = background(0x41, 30);
        bases.extend(b"GAT".repeat(8));
        bases.extend(background(0x42, 30));
        let criteria = DetectionCriteria::default();
        let aspirant = || {
            classify_interval(&bases, &interval(30, 54, 3, 46), &criteria)
                .expect("el (GAT)×8 pasa solo")
        };
        // Por la izquierda.
        for (overlap, kept) in [(3u64, Some((33u64, 54u64))), (4, None)] {
            let winner = candidate(20, 30 + overlap, b"CA", 1_000);
            let (accepted, lost) =
                resolve_overlaps(&bases, &criteria, vec![winner.clone(), aspirant()]);
            let got: Vec<(u64, u64)> = accepted
                .iter()
                .filter(|c| c.motif().length() == 3)
                .map(|c| (c.start(), c.end()))
                .collect();
            assert_eq!(
                got,
                kept.into_iter().collect::<Vec<_>>(),
                "izquierda {overlap}"
            );
            assert_eq!(
                lost.len(),
                usize::from(kept.is_none()),
                "izquierda {overlap}"
            );
        }
        // Por la derecha.
        for (overlap, kept) in [(3u64, Some((30u64, 51u64))), (4, None)] {
            let winner = candidate(54 - overlap, 64, b"CA", 1_000);
            let (accepted, lost) =
                resolve_overlaps(&bases, &criteria, vec![winner.clone(), aspirant()]);
            let got: Vec<(u64, u64)> = accepted
                .iter()
                .filter(|c| c.motif().length() == 3)
                .map(|c| (c.start(), c.end()))
                .collect();
            assert_eq!(
                got,
                kept.into_iter().collect::<Vec<_>>(),
                "derecha {overlap}"
            );
            assert_eq!(lost.len(), usize::from(kept.is_none()), "derecha {overlap}");
        }
    }

    /// **Un solape por los dos lados, uno dentro del umbral y otro no, se
    /// descarta entero**; y un aceptado **dentro** del aspirante o que lo
    /// **cubre** tampoco se recorta (no hay un resto con sentido). Escrito con
    /// ganadores a mano sobre un `(GAT)×8` en `[30, 54)`.
    #[test]
    fn overlaps_on_both_sides_or_around_the_aspirant_are_not_trimmed() {
        let mut bases = background(0x43, 30);
        bases.extend(b"GAT".repeat(8));
        bases.extend(background(0x44, 30));
        let criteria = DetectionCriteria::default();
        let aspirant = classify_interval(&bases, &interval(30, 54, 3, 46), &criteria)
            .expect("el (GAT)×8 pasa solo");
        let cases: [(&str, Vec<RepeatCandidate>, usize); 5] = [
            // Tres bases por la izquierda (dentro) y cuatro por la derecha.
            (
                "izquierda ok, derecha excede",
                vec![
                    candidate(20, 33, b"CA", 1_000),
                    candidate(50, 64, b"CA", 999),
                ],
                0,
            ),
            // Tres por cada lado: se recorta por los dos y quedan 18 bases.
            (
                "los dos lados dentro",
                vec![
                    candidate(20, 33, b"CA", 1_000),
                    candidate(51, 64, b"CA", 999),
                ],
                1,
            ),
            // Un aceptado dentro del aspirante.
            ("dentro", vec![candidate(40, 44, b"CA", 1_000)], 0),
            // Un aceptado que lo cubre entero.
            ("cubre", vec![candidate(20, 64, b"CA", 1_000)], 0),
            // Un aceptado que lo cubre entero y empieza justo en su inicio.
            (
                "cubre desde el inicio",
                vec![candidate(30, 60, b"CA", 1_000)],
                0,
            ),
        ];
        for (label, winners, expected_trimmed) in cases {
            let mut input = winners.clone();
            input.push(aspirant.clone());
            let (accepted, lost) = resolve_overlaps(&bases, &criteria, input);
            let trimmed = accepted.iter().filter(|c| c.motif().length() == 3).count();
            assert_eq!(trimmed, expected_trimmed, "{label}");
            assert_eq!(lost.len(), 1 - expected_trimmed, "{label}");
            assert_eq!(accepted.len(), winners.len() + expected_trimmed, "{label}");
        }
    }

    /// **Un recortado que ya no llega al suelo cae**, y cae en «solape», no en
    /// la etapa del suelo. `(CA)×8 A×7`: la tira de `A` mide 8 bases `[17, 25)`
    /// y pasa el suelo de 8 por sí sola; pierde la base que comparte con el
    /// dinucleótido, queda en 7 y **al volver a juzgarla** no lo alcanza. El
    /// `(CA)×8`, que gana, se queda.
    ///
    /// Esta prueba **no distingue** la regla de hoy de la nueva por el
    /// resultado (con las dos cae); fija que el recortado se **vuelve a
    /// juzgar** y no se acepta tal cual (aceptarlo daría `[18, 25)`), y dónde
    /// se cuenta.
    #[test]
    fn a_trimmed_tract_below_the_floor_falls_and_counts_as_overlap() {
        let core = [b"CA".repeat(8), b"A".repeat(7)].concat();
        let (candidates, counts) = detect_core(b"G", &core, b"G");
        assert_eq!(candidates, vec![(0, 16, b"CA".to_vec(), 28)]);
        assert!(counts.is_balanced(), "{counts:?}");
        assert_eq!(counts.discarded(DetectionStage::Overlap, 1), 1);
        // Cae en «solape», no en el suelo que ya no supera.
        assert_eq!(counts.accepted(1), 0);
    }

    /// **Un tramo que no cabe en su región se descarta antes de los solapes y
    /// no recorta a nadie.** `(CA)×10 A×12` con la región justo sobre la tira
    /// de `A` (`[19, 32)`, sus 13 bases): el `(CA)×10` no cabe y cae en
    /// «fuera de las regiones»; la tira se queda **entera**, con sus 13 bases,
    /// y no recortada a 12 como si el dinucleótido la hubiera pisado. Con la
    /// región sobre el `(CA)×10`, la tira cae fuera y el dinucleótido no
    /// cambia.
    ///
    /// **Debe fallar contra** aplicar las regiones **después** de los solapes:
    /// la tira saldría recortada, `[20, 32)`.
    #[test]
    fn a_tract_outside_the_regions_is_dropped_before_overlaps_and_trims_nobody() {
        let core = [b"CA".repeat(10), b"A".repeat(12)].concat();
        let mut sequence = background(0xD1, 30);
        let offset = sequence.len() as u64;
        sequence.push(b'G');
        sequence.extend_from_slice(&core);
        sequence.push(b'G');
        sequence.extend(background(0xD2, 30));
        let start = offset + 1;
        let criteria = DetectionCriteria::default();
        let run = |regions: &[(u64, u64)]| {
            let detection =
                detect_repeat_candidates_in_regions(&sequence, &criteria, Some(regions));
            let kept: Vec<(u64, u64)> = detection
                .candidates()
                .iter()
                .map(|c| (c.start() - start, c.end() - start))
                .collect();
            (kept, detection.counts().clone())
        };
        let (kept, counts) = run(&[(start + 19, start + 32)]);
        assert_eq!(kept, vec![(19, 32)]);
        assert_eq!(counts.discarded(DetectionStage::OutsideRegions, 2), 1);
        assert_eq!(counts.discarded_total(DetectionStage::Overlap), 0);
        assert!(counts.is_balanced(), "{counts:?}");
        let (kept, counts) = run(&[(start, start + 20)]);
        assert_eq!(kept, vec![(0, 20)]);
        assert_eq!(counts.discarded(DetectionStage::OutsideRegions, 1), 1);
        assert!(counts.is_balanced(), "{counts:?}");
    }

    /// **Un recorte demasiado corto para tener un motivo no cuelga ni
    /// entra.** Un aspirante de 6 bases de periodo 3 con dos bases que le
    /// quitan por un lado y una por el otro (todas dentro de una copia) deja
    /// 3 bases: justo un motivo, y `classify_interval` lo juzga. Con 2 por
    /// cada lado deja 2, menos que el periodo, y cae sin leer el motivo (que
    /// no cabe).
    #[test]
    fn a_trim_shorter_than_one_motif_falls_without_panicking() {
        let bases = [&b"CAGT"[..], &b"GAT".repeat(2), &b"CAGT"[..]].concat();
        let criteria = DetectionCriteria {
            min_copies: [1; 6],
            min_purity: 0.0,
            min_longest_run_units: 1,
            ..DetectionCriteria::default()
        };
        let aspirant = candidate(4, 10, b"GAT", 2);
        for (left, right, kept) in [(2u64, 1u64, true), (2, 2, false), (3, 3, false)] {
            let winners = vec![
                candidate(0, 4 + left, b"CAGT", 1_000),
                candidate(10 - right, 14, b"CAGT", 999),
            ];
            let mut input = winners.clone();
            input.push(aspirant.clone());
            let (accepted, lost) = resolve_overlaps(&bases, &criteria, input);
            assert_eq!(accepted.len(), 2 + usize::from(kept), "{left} {right}");
            assert_eq!(lost.len(), usize::from(!kept), "{left} {right}");
        }
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
    /// Los solapes de estas parejas son de más de una copia del que pierde, así
    /// que se descarta entero y las bases no se miran.
    fn winner_in_both_orders(a: RepeatCandidate, b: RepeatCandidate) -> RepeatCandidate {
        let bases = background(0x31, 100);
        let criteria = DetectionCriteria::default();
        let (forward, lost_forward) =
            resolve_overlaps(&bases, &criteria, vec![a.clone(), b.clone()]);
        let (backward, lost_backward) = resolve_overlaps(&bases, &criteria, vec![b, a]);
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
    /// un tercero que solo lo solapaba a él. Los solapes son de más de una
    /// copia del que pierde, de modo que se descarta entero (el recorte lo
    /// prueban las pruebas de `overlap_trim_*`).
    #[test]
    fn touching_tracts_are_kept_and_losers_do_not_block() {
        let bases = background(0x31, 100);
        let criteria = DetectionCriteria::default();
        let (kept, lost) = resolve_overlaps(
            &bases,
            &criteria,
            vec![candidate(10, 20, b"CA", 10), candidate(20, 30, b"GT", 10)],
        );
        assert_eq!((kept.len(), lost.len()), (2, 0));
        // A gana a B (solape de 10 bases, más de una copia de `GAT`); C solo
        // solapa a B y queda.
        let (kept, lost) = resolve_overlaps(
            &bases,
            &criteria,
            vec![
                candidate(10, 40, b"GAT", 20),
                candidate(0, 20, b"CA", 30),
                candidate(38, 60, b"AC", 10),
            ],
        );
        assert_eq!(
            kept.iter().map(RepeatCandidate::start).collect::<Vec<_>>(),
            vec![0, 38]
        );
        assert_eq!(lost, vec![3]);
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
    /// - solape, periodo 2 contra 4: `(GC)×8 (ACGC)×6`; el tramo de periodo 4
    ///   (`CGCA`) empieza tres bases antes de acabar el dinucleótido, **más de
    ///   una copia de este (dos bases)**, y gana por puntuación; el `(GC)×8`
    ///   se descarta entero. Con un solape de una copia o menos se recortaría
    ///   en vez de descartarse (`(CA)×10 G (CAG)×7` da los dos; lo fija
    ///   `overlap_trim_keeps_a_dinucleotide_overlapped_by_one_copy_by_a_trinucleotide`),
    ///   así que este caso es el que ejercita la etapa «solape» con el
    ///   descarte entero;
    /// - sin solape tras recortar, periodo 1 contra 2: `(CA)×40 A×8`. El
    ///   poli-A `[403, 412)` pisa **una base** al `(CA)×40` `[324, 404)`, una
    ///   copia de su motivo: se recorta a `[404, 412)`, 8 bases, justo el
    ///   suelo, y se queda; **no** cuenta en «solape»;
    /// - un poli-A **delante** de un `(CA)n` tampoco pierde: desde que el
    ///   inicio recortado sale de la tira de una sola base, los dos solo se
    ///   tocan (lo prueban las regresiones de `A×8 (GA)×10` y `T×12
    ///   (CTT)×7`).
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
        sequence.extend(b"GC".repeat(8));
        sequence.extend(b"ACGC".repeat(6));
        sequence.extend(background(0x105, 40));
        sequence.extend(b"CA".repeat(40));
        sequence.extend_from_slice(&[b'A'; 8]);
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
            counts.discarded(DetectionStage::Overlap, 1),
            0,
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
            vec![b"CGCA".to_vec(), b"CA".to_vec(), b"A".to_vec()],
            "quedan el (ACGC)n, el (CA)n y el poli-A recortado a una copia de solape"
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
