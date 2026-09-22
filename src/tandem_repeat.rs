//! Modelo puro del locus de repetición en tándem y de su tracto realizado.
//!
//! # La regla que gobierna el módulo entero
//!
//! Nada de lo anotado se almacena junto a la secuencia: **todo se deriva de
//! ella** cada vez que se pide, de modo que no puede contradecirla. La pureza
//! es siempre una **salida** calculada, nunca una entrada declarada. Esa es la
//! decisión S-02 del plan de repeticiones y tartamudeo, y la razón por la que
//! aquí no existe ningún campo `purity` ni ningún contador de unidades.
//!
//! # El intervalo no tiene por qué ser múltiplo del motivo
//!
//! Las interrupciones reales suelen ser de **una sola base** —una `T` dentro de
//! un `(CA)n`—, con lo que la longitud del tracto deja de ser múltiplo de la
//! del motivo y el número de unidades tiene parte fraccionaria. El propio
//! estándar VCF 4.5 trae como ejemplo `(CAG)5(CA)1(CAG)4`, que un modelo de
//! múltiplos exactos no podría representar.
//!
//! # Definiciones exactas
//!
//! - **Número de unidades:** la longitud del tracto dividida por la del motivo.
//!   Puede tener parte fraccionaria. Con un solo bloque se expone como fracción
//!   exacta (`exact_units`) además de como `f64`; con varios, solo como `f64`,
//!   porque cada bloque tiene su propio denominador.
//! - **Pureza:** la proporción de **bases** del tracto que pertenecen a copias
//!   completas y exactas del motivo. No es la proporción de unidades: contada
//!   sobre unidades, `CACACATCACA` daría 5/6, y sobre bases da 10/11, que es lo
//!   correcto.
//! - **Tirada ininterrumpida más larga:** el mayor número de copias completas
//!   **consecutivas y contiguas** del motivo. En caso de empate manda la
//!   primera por la izquierda, para que la operación de paso sea determinista.
//!
//! # Un locus puede tener **varios motivos**, en bloques
//!
//! Un microsatélite compuesto como `(TCTA)n(TCTG)n(TCTA)n` —la forma de muchos
//! marcadores forenses— no es un tracto de un solo motivo. Declararlo con uno
//! solo lo hace parecer muy impuro y lo rechaza la guarda de admisión, cuando en
//! realidad es perfectamente regular: lo que pasa es que su regularidad tiene
//! tres tramos.
//!
//! Por eso un tracto es una **lista ordenada de bloques**, cada uno con su
//! motivo y sus bases. Un tracto de un solo motivo es el caso de **un bloque**,
//! y se comporta exactamente igual que antes.
//!
//! **Esto acerca el modelo al estándar en vez de alejarlo.** El VCF 4.5 trae
//! como ejemplo literal `(CAG)10(TG)7(CAGG)3` y tiene campos pensados justo para
//! esto: el número de secuencias repetidas del alelo y **listas paralelas** de
//! secuencia del motivo, longitud del motivo, número de copias y bases totales.
//! La anotación de este módulo expone exactamente esas magnitudes por bloque.
//!
//! ## Qué se deriva y qué se declara
//!
//! La **estructura en bloques se declara**; no puede derivarse. Dadas unas bases
//! sueltas no hay forma de saber dónde acaba un motivo y empieza el siguiente:
//! esa es precisamente la información que el bloque aporta. Todo lo **numérico**
//! —copias completas, pureza, interrupciones, tiradas— se sigue derivando de las
//! bases, ahora dentro del tramo de cada bloque.
//!
//! ## Las tres reglas de un tracto con varios bloques
//!
//! 1. **La tirada ininterrumpida más larga se mide sobre todos los bloques y
//!    gana la más larga, sea del motivo que sea.** La polimerasa resbala donde
//!    hay más repetición seguida, no donde el usuario escribió el primer motivo.
//! 2. **El paso de unidades gana o pierde una unidad del bloque donde ocurre el
//!    resbalón**, es decir, del motivo de esa tirada.
//! 3. **La pureza se sigue midiendo sobre bases** y sobre el tracto completo,
//!    sumando las copias completas de todos los bloques.
//!
//! # La descomposición en copias, que hay que fijar porque no es única
//!
//! Las tres magnitudes anteriores se leen de una única descomposición del
//! tracto en copias completas del motivo: **voraz, de izquierda a derecha y sin
//! solapamiento**. Se avanza por el tracto; si en la posición actual empiezan
//! exactamente las bases del motivo se anota una copia y se salta su longitud
//! entera; si no, se avanza una sola base. Las bases que ninguna copia cubre
//! son las **interrupciones**.
//!
//! **La dirección es una decisión, no una consecuencia.** Una redacción anterior
//! de esta nota afirmaba que el barrido de izquierda a derecha es «la única
//! lectura que reproduce los ejemplos del plan». **Es falso**, y la revisión
//! adversarial del 20 de septiembre de 2026 lo demostró: el barrido de derecha a
//! izquierda reproduce los cuatro ejemplos con **idénticos números e idénticas
//! posiciones de interrupción**, porque todos sus motivos son no solapantes
//! consigo mismos. Lo que los ejemplos fijan es que la descomposición sea
//! **máxima y sin solapamiento**, no por dónde se empieza.
//!
//! Se elige la izquierda por coherencia con el resto del proyecto —el desempate
//! de la tirada más larga también manda a la izquierda, y la proyección al VCF
//! alinea los *indels* a la izquierda— y la diferencia **sí es observable**: con
//! un motivo de período interno, como el homopolímero `AA` sobre el tracto
//! `AAA`, la interrupción cae en el desplazamiento 2 por la izquierda y en el 0
//! por la derecha. Hay un fixture que lo fija, precisamente porque los
//! homopolímeros son la clase más común.

use std::sync::Arc;

use crate::coordinates::GenomicInterval;
use crate::identity::RepeatLocusId;

/// Motivo de repetición: no vacío y solo A, C, G o T.
///
/// `N` queda excluido a propósito: un motivo ambiguo no define ninguna copia
/// exacta, de modo que la pureza y la tirada más larga no serían computables.
#[derive(Clone, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct RepeatMotif {
    bases: Arc<[u8]>,
}

impl RepeatMotif {
    /// Valida el motivo en mayúsculas canónicas.
    pub fn new(bases: impl AsRef<[u8]>) -> Result<Self, TandemRepeatError> {
        let bases = bases.as_ref();
        if bases.is_empty() {
            return Err(TandemRepeatError::EmptyMotif);
        }
        for &base in bases {
            if !matches!(base, b'A' | b'C' | b'G' | b'T') {
                return Err(TandemRepeatError::NonCanonicalMotifBase { base });
            }
        }
        Ok(Self {
            bases: Arc::from(bases),
        })
    }

    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    /// Longitud del motivo en pares de base. Nunca es cero.
    #[must_use]
    pub fn length(&self) -> u64 {
        // El constructor rechaza el motivo vacío y ningún camino lo modifica.
        self.bases.len() as u64
    }

    /// Texto ASCII del motivo, que es lo que publican el VCF y la verdad.
    #[must_use]
    pub fn as_text(&self) -> String {
        String::from_utf8(self.bases.to_vec())
            .unwrap_or_else(|_| unreachable!("el constructor solo admite A, C, G y T"))
    }
}

/// Una interrupción: bases contiguas del tracto que ninguna copia completa
/// cubre. Los desplazamientos son relativos al inicio del tracto.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct Interruption {
    start: u64,
    length: u64,
}

impl Interruption {
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    #[must_use]
    pub const fn length(self) -> u64 {
        self.length
    }
}

/// Una tirada: copias completas del motivo consecutivas y contiguas.
///
/// El desplazamiento es **absoluto dentro del tracto**, no relativo al bloque,
/// para que comparar tiradas de bloques distintos no exija recolocarlas.
#[derive(Clone, Copy, Debug, Eq, Hash, Ord, PartialEq, PartialOrd)]
pub struct RepeatRun {
    start: u64,
    units: u64,
    block_index: usize,
}

impl RepeatRun {
    /// Desplazamiento en bases del comienzo de la tirada dentro del tracto.
    #[must_use]
    pub const fn start(self) -> u64 {
        self.start
    }

    /// Copias completas y contiguas que la componen.
    #[must_use]
    pub const fn units(self) -> u64 {
        self.units
    }

    /// Bloque al que pertenece. Es el bloque cuyo motivo gana o pierde una
    /// unidad cuando el resbalón cae en esta tirada.
    #[must_use]
    pub const fn block_index(self) -> usize {
        self.block_index
    }
}

/// Anotación de **un bloque**: su motivo y lo que sus bases dicen de él.
///
/// Es la unidad que el VCF 4.5 llama «secuencia repetida» del alelo, y de aquí
/// salen sus listas paralelas: la secuencia del motivo, su longitud, el número
/// de copias y las bases totales del bloque.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct BlockAnnotation {
    motif: RepeatMotif,
    start: u64,
    length: u64,
    complete_copies: u64,
    runs: Vec<RepeatRun>,
    interruptions: Vec<Interruption>,
}

impl BlockAnnotation {
    #[must_use]
    pub const fn motif(&self) -> &RepeatMotif {
        &self.motif
    }

    /// Desplazamiento del bloque dentro del tracto, en base cero.
    #[must_use]
    pub const fn start(&self) -> u64 {
        self.start
    }

    /// Bases del bloque. Es el campo de bases totales del VCF.
    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }

    #[must_use]
    pub const fn complete_copies(&self) -> u64 {
        self.complete_copies
    }

    /// Número de copias del bloque, que puede no ser entero. Es el campo de
    /// número de copias del VCF, que por eso es decimal y no entero.
    #[must_use]
    pub fn units_as_f64(&self) -> f64 {
        self.length as f64 / self.motif.length() as f64
    }

    #[must_use]
    pub fn runs(&self) -> &[RepeatRun] {
        &self.runs
    }

    /// Interrupciones **de este bloque**, con desplazamiento absoluto en el
    /// tracto y cortadas en la frontera del bloque.
    ///
    /// Para la cuenta biológica —tramos maximales, sin partir en la frontera—
    /// hay que usar [`TractAnnotation::interruptions`], que es la que se
    /// publica.
    #[must_use]
    pub fn interruptions(&self) -> &[Interruption] {
        &self.interruptions
    }
}

/// Anotación completa **derivada** de las bases del tracto, bloque a bloque.
///
/// No se construye desde fuera: la única forma de obtenerla es pedírsela a un
/// [`TandemRepeatTract`], que la recalcula desde sus bases. Así no puede existir
/// una anotación que contradiga la secuencia que dice describir.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TractAnnotation {
    total_bases: u64,
    blocks: Vec<BlockAnnotation>,
}

impl TractAnnotation {
    /// Bases del tracto entero.
    #[must_use]
    pub const fn total_bases(&self) -> u64 {
        self.total_bases
    }

    /// Los bloques, en orden. Su número es el campo del VCF que cuenta cuántas
    /// secuencias repetidas tiene el alelo.
    #[must_use]
    pub fn blocks(&self) -> &[BlockAnnotation] {
        &self.blocks
    }

    /// Fracción exacta del número de unidades, **y solo cuando hay un bloque**.
    ///
    /// Devuelve `None` con varios bloques, porque entonces la fracción única no
    /// existe: cada bloque tiene su propio denominador. Antes había dos métodos
    /// que con varios bloques devolvían las bases totales y la longitud del
    /// motivo del **primer** bloque: un número sin significado y sin aviso, y
    /// una trampa puesta para quien fuera a publicar la verdad y el VCF.
    #[must_use]
    pub fn exact_units(&self) -> Option<(u64, u64)> {
        match self.blocks.as_slice() {
            [only] => Some((self.total_bases, only.motif.length())),
            _ => None,
        }
    }

    /// Número **total** de unidades: la suma de las de cada bloque.
    ///
    /// Con un solo bloque coincide exactamente con las bases entre la longitud
    /// del motivo, que es lo que valía antes de que hubiera bloques.
    #[must_use]
    pub fn units_as_f64(&self) -> f64 {
        self.blocks.iter().map(BlockAnnotation::units_as_f64).sum()
    }

    /// Copias completas y exactas halladas en todos los bloques.
    #[must_use]
    pub fn complete_copies(&self) -> u64 {
        self.blocks.iter().map(|block| block.complete_copies).sum()
    }

    /// Numerador exacto de la pureza: **bases** cubiertas por copias completas,
    /// sumadas sobre todos los bloques.
    #[must_use]
    pub fn purity_numerator(&self) -> u64 {
        self.blocks
            .iter()
            .map(|block| block.complete_copies * block.motif.length())
            .sum()
    }

    /// Denominador exacto de la pureza: las bases del tracto entero.
    #[must_use]
    pub const fn purity_denominator(&self) -> u64 {
        self.total_bases
    }

    /// Pureza como decimal en `[0, 1]`.
    #[must_use]
    pub fn purity_as_f64(&self) -> f64 {
        self.purity_numerator() as f64 / self.total_bases as f64
    }

    /// Las interrupciones del tracto, **maximales**, en orden de izquierda a
    /// derecha.
    ///
    /// Una interrupción que cruza la frontera entre dos bloques es **una**, no
    /// dos. La contabilidad por bloque la parte en el punto donde el usuario
    /// decidió cortar la estructura, y eso no es una propiedad de la secuencia:
    /// el número de interrupciones de un microsatélite tiene significado
    /// biológico y no puede depender de dónde se haya puesto la frontera. Por
    /// eso aquí los tramos adyacentes se funden, y es esta lista —no la de cada
    /// bloque— la que se publica.
    ///
    /// La lista de cada bloque, en [`BlockAnnotation::interruptions`], sigue
    /// siendo la contabilidad **de ese bloque**: sirve para mirar un bloque por
    /// separado, y ahí la frontera sí es el límite.
    #[must_use]
    pub fn interruptions(&self) -> Vec<Interruption> {
        let mut merged: Vec<Interruption> = Vec::new();
        for interruption in self
            .blocks
            .iter()
            .flat_map(|block| block.interruptions.iter())
        {
            match merged.last_mut() {
                Some(previous) if previous.start + previous.length == interruption.start => {
                    previous.length += interruption.length;
                }
                _ => merged.push(*interruption),
            }
        }
        merged
    }

    /// Todas las tiradas del tracto, en orden de izquierda a derecha.
    #[must_use]
    pub fn runs(&self) -> Vec<RepeatRun> {
        self.blocks
            .iter()
            .flat_map(|block| block.runs.iter().copied())
            .collect()
    }

    /// Tirada ininterrumpida más larga **de todo el tracto**, mire al motivo que
    /// mire, con desempate por la izquierda.
    ///
    /// Es la regla 1 de los tractos con varios bloques: la polimerasa resbala
    /// donde hay más repetición seguida, no donde el usuario escribió el primer
    /// motivo. Devuelve `None` solo cuando ningún bloque contiene una copia
    /// completa de su motivo.
    #[must_use]
    pub fn longest_run(&self) -> Option<RepeatRun> {
        let mut best: Option<RepeatRun> = None;
        for block in &self.blocks {
            for &run in &block.runs {
                // Estrictamente mayor: el empate conserva la primera por la
                // izquierda, y como los bloques se recorren en orden y sus
                // tiradas también, «la primera» es la de menor desplazamiento.
                if best.is_none_or(|current| run.units > current.units) {
                    best = Some(run);
                }
            }
        }
        best
    }

    /// Unidades de la tirada más larga; cero si no hay ninguna copia completa.
    #[must_use]
    pub fn longest_run_units(&self) -> u64 {
        self.longest_run().map_or(0, RepeatRun::units)
    }
}

/// Un bloque del tracto: su motivo y cuántas bases ocupa.
///
/// No guarda las bases: viven una sola vez, en el búfer del tracto. El bloque
/// dice qué motivo rige ese tramo y dónde acaba.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TractBlock {
    motif: RepeatMotif,
    length: u64,
}

impl TractBlock {
    pub fn new(motif: RepeatMotif, length: u64) -> Result<Self, TandemRepeatError> {
        if length == 0 {
            return Err(TandemRepeatError::EmptyBlock);
        }
        Ok(Self { motif, length })
    }

    #[must_use]
    pub const fn motif(&self) -> &RepeatMotif {
        &self.motif
    }

    #[must_use]
    pub const fn length(&self) -> u64 {
        self.length
    }
}

/// Tracto realizado: su estructura en bloques y la **secuencia concreta** que
/// la realiza.
///
/// Las bases son la autoridad de todo lo numérico. La estructura en bloques es
/// declarada, porque no puede derivarse de las bases: dadas unas bases sueltas
/// no hay forma de saber dónde acaba un motivo y empieza el siguiente.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TandemRepeatTract {
    blocks: Arc<[TractBlock]>,
    bases: Arc<[u8]>,
}

impl TandemRepeatTract {
    /// Construye el tracto de **un solo bloque**, que es el caso de siempre.
    ///
    /// Se admite `N` en el tracto, aunque no en el motivo: una referencia real
    /// puede traerlas, y una `N` sencillamente nunca pertenece a una copia
    /// exacta, de modo que cuenta como interrupción sin ninguna regla extra.
    pub fn new(motif: RepeatMotif, bases: impl AsRef<[u8]>) -> Result<Self, TandemRepeatError> {
        Self::from_owned_bases(motif, bases.as_ref().to_vec())
    }

    /// Igual que [`Self::new`], pero **consume** el búfer en vez de leerlo por
    /// referencia, de modo que no se recorre dos veces un tracto que el
    /// llamador acaba de construir.
    ///
    /// No reduce el pico a la mitad: `Arc<[u8]>` no puede adoptar la
    /// asignación de un `Vec`, así que la conversión sigue asignando y
    /// copiando una vez. Lo que ahorra es la pasada extra y la vida útil
    /// solapada más larga. **El pico de verdad lo evita rechazar antes de
    /// materializar**, que es lo que hace la frontera de configuración con el
    /// tope en unidades.
    pub fn from_owned_bases(motif: RepeatMotif, bases: Vec<u8>) -> Result<Self, TandemRepeatError> {
        validate_tract_alphabet(&bases)?;
        let block = TractBlock::new(motif, bases.len() as u64)?;
        Ok(Self {
            blocks: Arc::from(vec![block]),
            bases: Arc::from(bases),
        })
    }

    /// Construye un tracto de **varios bloques** desde su estructura y sus
    /// bases. La suma de las longitudes de los bloques tiene que ser
    /// exactamente la del búfer.
    pub fn from_blocks(blocks: Vec<TractBlock>, bases: Vec<u8>) -> Result<Self, TandemRepeatError> {
        validate_tract_alphabet(&bases)?;
        if blocks.is_empty() {
            return Err(TandemRepeatError::NoBlocks);
        }
        let declared: u64 = blocks.iter().map(TractBlock::length).sum();
        let observed = bases.len() as u64;
        if declared != observed {
            return Err(TandemRepeatError::BlockLengthsDoNotCoverTract { declared, observed });
        }
        Ok(Self {
            blocks: Arc::from(blocks),
            bases: Arc::from(bases),
        })
    }

    /// Materializa un tracto puro de `units` copias exactas de un solo motivo.
    pub fn pure(motif: &RepeatMotif, units: u64) -> Result<Self, TandemRepeatError> {
        Self::pure_blocks(&[(motif.clone(), units)])
    }

    /// Materializa un tracto puro de varios bloques, cada uno con su motivo y
    /// su número de copias.
    pub fn pure_blocks(specification: &[(RepeatMotif, u64)]) -> Result<Self, TandemRepeatError> {
        if specification.is_empty() {
            return Err(TandemRepeatError::NoBlocks);
        }
        let mut total = 0_usize;
        for (motif, units) in specification {
            if *units == 0 {
                return Err(TandemRepeatError::EmptyBlock);
            }
            let motif_length =
                usize::try_from(motif.length()).map_err(|_| TandemRepeatError::TractTooLong {
                    requested: u64::MAX,
                })?;
            let block = usize::try_from(*units)
                .ok()
                .and_then(|units| units.checked_mul(motif_length))
                .ok_or(TandemRepeatError::TractTooLong {
                    requested: units.saturating_mul(motif.length()),
                })?;
            total = total
                .checked_add(block)
                .ok_or(TandemRepeatError::TractTooLong {
                    requested: u64::MAX,
                })?;
        }
        let mut bases = Vec::with_capacity(total);
        let mut blocks = Vec::with_capacity(specification.len());
        for (motif, units) in specification {
            for _ in 0..*units {
                bases.extend_from_slice(motif.bases());
            }
            blocks.push(TractBlock::new(motif.clone(), units * motif.length())?);
        }
        Self::from_blocks(blocks, bases)
    }

    /// Los bloques declarados, en orden.
    #[must_use]
    pub fn blocks(&self) -> &[TractBlock] {
        &self.blocks
    }

    /// Motivo del tracto **cuando solo tiene un bloque**.
    ///
    /// Con varios bloques no hay «el motivo» del tracto: cada bloque tiene el
    /// suyo, y hay que preguntarle a [`Self::blocks`].
    #[must_use]
    pub fn motif(&self) -> Option<&RepeatMotif> {
        match self.blocks.as_ref() {
            [only] => Some(only.motif()),
            _ => None,
        }
    }

    #[must_use]
    pub fn bases(&self) -> &[u8] {
        &self.bases
    }

    #[must_use]
    pub fn length(&self) -> u64 {
        self.bases.len() as u64
    }

    /// Texto ASCII del tracto, que es lo que van a llevar `REF` y `ALT`.
    #[must_use]
    pub fn as_text(&self) -> String {
        String::from_utf8(self.bases.to_vec())
            .unwrap_or_else(|_| unreachable!("el constructor solo admite A, C, G, T y N"))
    }

    /// Recalcula la anotación **desde las bases**, bloque a bloque, cada vez
    /// que se pide.
    #[must_use]
    pub fn annotation(&self) -> TractAnnotation {
        let mut blocks = Vec::with_capacity(self.blocks.len());
        let mut offset = 0_u64;
        for (block_index, block) in self.blocks.iter().enumerate() {
            let start = offset;
            let end = offset + block.length;
            blocks.push(annotate_block(
                &self.bases[start as usize..end as usize],
                block.motif(),
                start,
                block_index,
            ));
            offset = end;
        }
        TractAnnotation {
            total_bases: self.bases.len() as u64,
            blocks,
        }
    }

    /// Aplica un paso de `step` unidades **en la tirada ininterrumpida más
    /// larga del tracto entero**, con desempate por la primera por la izquierda.
    ///
    /// Con varios bloques, la unidad que se gana o se pierde es **la del motivo
    /// de esa tirada**, y solo su bloque cambia de longitud; los demás quedan
    /// byte a byte iguales. El resto del tracto también, de modo que la
    /// diferencia respecto del original es una inserción o una deleción
    /// contigua y limpia. Las unidades se añaden o se quitan en el extremo
    /// derecho de la tirada; como todas sus copias son idénticas, da lo mismo
    /// dónde dentro de ella, y fijarlo mantiene la operación determinista.
    ///
    /// # Decisión de modelización: una contracción no puede vaciar su tirada
    ///
    /// La contracción debe dejar **al menos una** copia completa en la tirada
    /// sobre la que actúa. Es una **decisión propia**, no algo que el plan pida,
    /// y se deja escrita aquí con sus tres consecuencias porque un comentario de
    /// paso no basta para algo que cambia la biología del modelo:
    ///
    /// 1. **Es por tirada, no por tracto.** Un tracto con tiradas de 5 y 3
    ///    unidades puede perder hasta 4 de la primera; lo que no puede es
    ///    quedarse sin ella.
    /// 2. **El motivo es que vaciar una tirada mueve la pureza por una decisión
    ///    nuestra y no por biología**, que es justo lo que la decisión S-07
    ///    retira al abandonar la convención del extremo derecho.
    /// 3. **Tiene una frontera más grande de lo que parece, y hay que verla.**
    ///    Una redacción anterior decía solo que «el número de tiradas nunca
    ///    disminuye». Se queda corta. La consecuencia real es que **un tracto
    ///    cuya tirada más larga sea de una sola copia no puede contraerse nunca,
    ///    mientras que expandirse sí puede**: solo puede crecer. Y eso puede
    ///    pasarle a un tracto con muchas copias completas, si están repartidas
    ///    en tiradas de una. Es decir, **una clase entera de tractos queda
    ///    congelada frente a la contracción con la expansión libre**, lo que
    ///    sesga el modelo hacia el alargamiento.
    ///
    ///    Hasta la decisión S-31 esa clase no llegaba a entrar, porque la tirada
    ///    mínima de admisión general es 3. Desde que **un locus puede declarar la
    ///    suya**, uno con tirada mínima 1 la admite. El sesgo dejó de ser teórico
    ///    el día que los mínimos pasaron a ser por locus.
    ///
    ///    **Salida acordada por el responsable el 21 de septiembre de 2026:**
    ///    el modelo **no se toca**, y en su lugar la admisión exigirá que la
    ///    tirada mínima efectiva de un locus sea **de al menos 2 unidades
    ///    cuando ese locus vaya a sufrir tartamudeo**. Con eso desaparece la
    ///    clase de tractos que no pueden contraerse nunca, y la asimetría que
    ///    queda —con una tirada de dos copias siguen prohibidos los pasos de
    ///    dos y tres— es **físicamente correcta**: en la realidad tampoco se
    ///    pueden perder tres unidades de una tirada que solo tiene dos.
    ///
    ///    **Esa guarda no vive aquí y todavía no existe.** El microsatélite se
    ///    declara junto a la referencia y el tartamudeo se configura en el
    ///    evento técnico, de modo que la condición solo puede comprobarse
    ///    cuando las dos configuraciones están resueltas, no al leer el locus.
    ///    Está planificada como el paso **F4.8.1.1** de la entrega 2 del plan
    ///    `design/plan-implementacion-repeticiones-pcr-0.1.md`, con su
    ///    razonamiento completo y sus variantes defectuosas.
    pub fn apply_unit_step(&self, step: i64) -> Result<Self, TandemRepeatError> {
        if step == 0 {
            return Err(TandemRepeatError::ZeroStep);
        }
        let annotation = self.annotation();
        let run = annotation
            .longest_run()
            .ok_or(TandemRepeatError::NoCompleteUnit)?;
        let block = &self.blocks[run.block_index()];
        let motif_length = block.motif().length();
        let run_end = run.start() + run.units() * motif_length;
        let run_end_index = usize::try_from(run_end)
            .map_err(|_| TandemRepeatError::TractTooLong { requested: run_end })?;

        let mut bases = Vec::with_capacity(self.bases.len() + motif_length as usize);
        let delta = if step > 0 {
            let added = u64::try_from(step).unwrap_or(0);
            bases.extend_from_slice(&self.bases[..run_end_index]);
            for _ in 0..added {
                bases.extend_from_slice(block.motif().bases());
            }
            bases.extend_from_slice(&self.bases[run_end_index..]);
            i64::try_from(added * motif_length).unwrap_or(i64::MAX)
        } else {
            let removed = step.unsigned_abs();
            if removed >= run.units() {
                return Err(TandemRepeatError::StepEmptiesLongestRun {
                    run_units: run.units(),
                    requested: removed,
                });
            }
            let cut = usize::try_from(run_end - removed * motif_length)
                .map_err(|_| TandemRepeatError::TractTooLong { requested: run_end })?;
            bases.extend_from_slice(&self.bases[..cut]);
            bases.extend_from_slice(&self.bases[run_end_index..]);
            -i64::try_from(removed * motif_length).unwrap_or(i64::MAX)
        };
        // Solo cambia de longitud el bloque donde cayó el resbalón.
        let mut blocks = self.blocks.to_vec();
        let affected = &mut blocks[run.block_index()];
        let updated =
            affected
                .length
                .checked_add_signed(delta)
                .ok_or(TandemRepeatError::TractTooLong {
                    requested: affected.length,
                })?;
        *affected = TractBlock::new(affected.motif.clone(), updated)?;
        Self::from_blocks(blocks, bases)
    }
}

/// Descompone un bloque en copias completas de su motivo: **voraz, de izquierda
/// a derecha y sin solapamiento**, dentro del tramo del bloque.
fn annotate_block(
    bases: &[u8],
    motif: &RepeatMotif,
    block_start: u64,
    block_index: usize,
) -> BlockAnnotation {
    let motif_bases = motif.bases();
    let motif_length = motif_bases.len();
    let total = bases.len();
    let mut copy_starts: Vec<u64> = Vec::new();
    let mut offset = 0_usize;
    while offset + motif_length <= total {
        if &bases[offset..offset + motif_length] == motif_bases {
            copy_starts.push(block_start + offset as u64);
            offset += motif_length;
        } else {
            offset += 1;
        }
    }

    let motif_length_u64 = motif_length as u64;
    let mut runs: Vec<RepeatRun> = Vec::new();
    for &start in &copy_starts {
        match runs.last_mut() {
            Some(run) if run.start + run.units * motif_length_u64 == start => {
                run.units += 1;
            }
            _ => runs.push(RepeatRun {
                start,
                units: 1,
                block_index,
            }),
        }
    }

    let mut interruptions: Vec<Interruption> = Vec::new();
    let mut cursor = block_start;
    for run in &runs {
        if run.start > cursor {
            interruptions.push(Interruption {
                start: cursor,
                length: run.start - cursor,
            });
        }
        cursor = run.start + run.units * motif_length_u64;
    }
    let block_end = block_start + total as u64;
    if cursor < block_end {
        interruptions.push(Interruption {
            start: cursor,
            length: block_end - cursor,
        });
    }

    BlockAnnotation {
        motif: motif.clone(),
        start: block_start,
        length: total as u64,
        complete_copies: copy_starts.len() as u64,
        runs,
        interruptions,
    }
}

fn validate_tract_alphabet(bases: &[u8]) -> Result<(), TandemRepeatError> {
    if bases.is_empty() {
        return Err(TandemRepeatError::EmptyTract);
    }
    for &base in bases {
        if !matches!(base, b'A' | b'C' | b'G' | b'T' | b'N') {
            return Err(TandemRepeatError::NonCanonicalTractBase { base });
        }
    }
    Ok(())
}

/// Lo que se **observa** en un tramo de una referencia real, sin estructura en
/// bloques declarada.
///
/// # Por qué hace falta un tipo aparte
///
/// Un alelo declarado trae su estructura en bloques: el usuario dice qué motivo
/// rige cada tramo. Un tracto **leído de una referencia** no la trae. La guarda
/// de admisión tiene que decidir si en la posición declarada hay realmente un
/// tracto de los motivos declarados, y para eso hay que mirar las bases sin
/// saber dónde estarían las fronteras.
///
/// # La convención de atribución, que hay que fijar porque no es única
///
/// Se recorre el tramo **de izquierda a derecha**, y en cada posición se prueban
/// los motivos declarados **en el orden en que se declararon**; el primero que
/// encaje se anota como copia y se salta su longitud; si no encaja ninguno, se
/// avanza una base. Es la misma descomposición voraz de siempre, extendida a
/// varios motivos, y **con un solo motivo es exactamente la de siempre**.
///
/// Una **tirada** es una cadena de copias completas contiguas **del mismo
/// motivo**: es donde la polimerasa resbala, y dos copias de motivos distintos
/// pegadas no son una tirada de dos.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub struct ObservedTract {
    total_bases: u64,
    covered_bases: u64,
    longest_run_units: u64,
}

impl ObservedTract {
    #[must_use]
    pub const fn total_bases(self) -> u64 {
        self.total_bases
    }

    /// Bases que pertenecen a una copia completa de alguno de los motivos.
    #[must_use]
    pub const fn covered_bases(self) -> u64 {
        self.covered_bases
    }

    /// Pureza **observada**, sobre bases, como en un alelo declarado.
    #[must_use]
    pub fn purity_as_f64(self) -> f64 {
        if self.total_bases == 0 {
            return 0.0;
        }
        self.covered_bases as f64 / self.total_bases as f64
    }

    /// Copias contiguas del mismo motivo en la tirada más larga observada.
    #[must_use]
    pub const fn longest_run_units(self) -> u64 {
        self.longest_run_units
    }
}

/// Observa un tramo de referencia contra los motivos de un locus.
///
/// No valida el alfabeto: las bases vienen de una referencia ya validada, y una
/// `N` sencillamente no encaja con ningún motivo y cuenta como no cubierta.
#[must_use]
pub fn observe_tract(bases: &[u8], motifs: &[RepeatMotif]) -> ObservedTract {
    let total = bases.len();
    let mut covered = 0_u64;
    let mut longest = 0_u64;
    let mut current_run = 0_u64;
    let mut current_motif: Option<&[u8]> = None;
    let mut offset = 0_usize;
    while offset < total {
        let matched = motifs.iter().map(RepeatMotif::bases).find(|motif| {
            offset + motif.len() <= total && &bases[offset..offset + motif.len()] == *motif
        });
        match matched {
            Some(motif) => {
                if current_motif == Some(motif) {
                    current_run += 1;
                } else {
                    current_motif = Some(motif);
                    current_run = 1;
                }
                longest = longest.max(current_run);
                covered += motif.len() as u64;
                offset += motif.len();
            }
            None => {
                current_motif = None;
                current_run = 0;
                offset += 1;
            }
        }
    }
    ObservedTract {
        total_bases: total as u64,
        covered_bases: covered,
        longest_run_units: longest,
    }
}

/// Locus de repetición **declarado**: dónde está el tracto y con qué motivos.
///
/// El locus no lleva secuencia. La secuencia es del alelo, y por eso vive en
/// [`TandemRepeatTract`]: dos individuos comparten el locus y no el tracto. Lo
/// que el locus sí fija es la **estructura**: la lista ordenada de motivos que
/// todos sus alelos comparten.
#[derive(Clone, Debug, Eq, PartialEq)]
pub struct TandemRepeatLocus {
    id: RepeatLocusId,
    interval: GenomicInterval,
    motifs: Arc<[RepeatMotif]>,
}

impl TandemRepeatLocus {
    /// El intervalo declarado **no** tiene que ser múltiplo de los motivos, pero
    /// sí tiene que caber al menos una copia de cada uno: sin ella no hay tracto
    /// que anotar ni tirada donde resbalar.
    pub fn new(
        id: RepeatLocusId,
        interval: GenomicInterval,
        motifs: Vec<RepeatMotif>,
    ) -> Result<Self, TandemRepeatError> {
        if motifs.is_empty() {
            return Err(TandemRepeatError::NoBlocks);
        }
        let minimum: u64 = motifs.iter().map(RepeatMotif::length).sum();
        if interval.len() < minimum {
            return Err(TandemRepeatError::IntervalShorterThanMotif {
                interval_length: interval.len(),
                motif_length: minimum,
            });
        }
        Ok(Self {
            id,
            interval,
            motifs: Arc::from(motifs),
        })
    }

    #[must_use]
    pub const fn id(&self) -> &RepeatLocusId {
        &self.id
    }

    #[must_use]
    pub const fn interval(&self) -> &GenomicInterval {
        &self.interval
    }

    /// Los motivos del locus, en orden. Un locus de un solo motivo tiene uno.
    #[must_use]
    pub fn motifs(&self) -> &[RepeatMotif] {
        &self.motifs
    }

    /// El motivo **cuando solo hay uno**. Con varios no existe «el motivo» del
    /// locus y hay que preguntar por la lista.
    #[must_use]
    pub fn single_motif(&self) -> Option<&RepeatMotif> {
        match self.motifs.as_ref() {
            [only] => Some(only),
            _ => None,
        }
    }
}

/// Errores del modelo puro. La traducción a diagnóstico público, con el
/// identificador del locus, la hace la frontera de configuración.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum TandemRepeatError {
    EmptyMotif,
    NoBlocks,
    EmptyBlock,
    BlockLengthsDoNotCoverTract {
        declared: u64,
        observed: u64,
    },
    NonCanonicalMotifBase {
        base: u8,
    },
    EmptyTract,
    NonCanonicalTractBase {
        base: u8,
    },
    TractTooLong {
        requested: u64,
    },
    IntervalShorterThanMotif {
        interval_length: u64,
        motif_length: u64,
    },
    ZeroStep,
    NoCompleteUnit,
    StepEmptiesLongestRun {
        run_units: u64,
        requested: u64,
    },
}

impl TandemRepeatError {
    /// Detalle estable en inglés para el contexto del diagnóstico público.
    #[must_use]
    pub fn detail(self) -> String {
        match self {
            Self::EmptyMotif => "the repeat motif is empty".to_owned(),
            Self::NoBlocks => "the repeat tract declares no block".to_owned(),
            Self::EmptyBlock => "a repeat block declares no base".to_owned(),
            Self::BlockLengthsDoNotCoverTract { declared, observed } => {
                format!("the declared blocks cover {declared} bases but the tract has {observed}")
            }
            Self::NonCanonicalMotifBase { base } => {
                format!("the repeat motif contains 0x{base:02x}, which is not A, C, G or T")
            }
            Self::EmptyTract => "the repeat tract is empty".to_owned(),
            Self::NonCanonicalTractBase { base } => {
                format!("the repeat tract contains 0x{base:02x}, which is not A, C, G, T or N")
            }
            Self::TractTooLong { requested } => {
                format!("the repeat tract of {requested} bases exceeds the address space")
            }
            Self::IntervalShorterThanMotif {
                interval_length,
                motif_length,
            } => format!(
                "the declared interval of {interval_length} bases is shorter than the {motif_length}-base motif"
            ),
            Self::ZeroStep => "a slippage step of zero units changes nothing".to_owned(),
            Self::NoCompleteUnit => "the tract contains no complete copy of the motif".to_owned(),
            Self::StepEmptiesLongestRun {
                run_units,
                requested,
            } => format!(
                "removing {requested} units would empty the longest run of {run_units} units"
            ),
        }
    }
}
