# Diseño de bibliotecas preparadas con PCR 0.1

- **Fecha:** 19 de septiembre de 2026.
- **Estado:** propuesta de diseño. No es normativa. Ninguna cifra, requisito
  ni decisión de este documento está aplicada; todo cambio a otros documentos
  se enumera en la sección 9 y queda pendiente de aprobación.
- **Autoridad afectada:** `especificaciones-0.1.md` (F4, §3 y §4.4),
  `registro-ampliaciones-futuras-0.1.md` (DTF-056), `hoja-de-ruta-0.1.md`
  (§F4), `plan-implementacion-f4-0.1.md`, `decisiones-tecnicas-0.1.md`,
  `trazabilidad-requisitos-0.1.md`, `calibracion-tartamudeo-0.1.md` y
  `expediente-decisiones-repeticiones-tandem-0.1.md` (ficha D-12).
- **Decisión de partida, ya tomada por el usuario y no reabierta aquí:** el
  simulador debe poder producir lecturas de una biblioteca WGS Illumina
  **amplificada por PCR**, con sus duplicados, su tartamudeo clonal y su sesgo
  de amplificación, además de la biblioteca sin PCR que ya existe.

---

## 1. Lo esencial en una página

**Qué se propone.** Un modo de preparación de biblioteca, `pcr`, alternativo
al único que hoy existe (`pcr_free`). Con él, el recorrido técnico de F4 deja
de producir una pareja por fragmento y pasa a producir **familias de
amplificación**: grupos de parejas que descienden de la misma molécula
original. Todo lo que la PCR hace observable —duplicados, sobrerrepresentación
de unas moléculas frente a otras, sesgo por composición de bases y tartamudeo
que se hereda dentro de la familia— sale de esa estructura de familias.

**Cómo se resuelve el problema duro.** No se simula la población de moléculas
que produce la amplificación, que crece como `2^ciclos` y no se puede acotar.
Se muestrea **hacia atrás desde las parejas que se van a publicar**: como el
número de parejas `N` de una biblioteca ya es determinista desde la cobertura,
se generan familias una tras otra hasta cubrir esas `N` parejas, y cada
familia cuesta un número de operaciones acotado por constantes de contrato.
Dentro de cada familia, lo que en la realidad es un árbol exponencial se
sustituye por dos piezas de forma cerrada —la distribución del número de
copias tras los primeros ciclos y la distribución del tamaño de familia— y por
una genealogía pequeña, de tamaño igual al de la familia, que es lo único que
hace falta para que el tartamudeo sea clonal. La cota resultante es lineal en
`N` y está en la sección 6.

**Qué invariante no se toca.** La verdad genotípica publicada —`truth.jsonl`,
el VCF de verdad y sus manifiestos— no cambia con ningún parámetro de PCR. La
PCR solo afecta a los FASTQ y a los recuentos agregados del manifiesto de
evento. Con `pcr_free`, que sigue siendo el valor predeterminado, todos los
artefactos publicados son **byte a byte** los de hoy.

> **Alcance exacto de esa congelación (corregido el 19 de septiembre de 2026).**
> Vale para los **escenarios que no declaran loci de repetición**, que son
> todos los que existen hoy. Un escenario que sí los declare **debe** producir
> tartamudeo incluso con `pcr_free`, porque el fenómeno existe sin
> amplificación y está medido: es la columna sin PCR de
> `design/calibracion-tartamudeo-0.1.md`. Leer la congelación como absoluta
> dejaría el caso sin PCR en cero tartamudeo, que contradice la calibración del
> propio proyecto, y **rompería además la aritmética del caso con PCR**, como
> se explica en la sección 5.7.

**Qué cambia de contrato.** Cinco sitios de la especificación declaran que F4
es «exclusivamente sin PCR» (REQ-TEC-07, REQ-TEC-08, REQ-CONF-88, REQ-DAT-46 y
ACE-49, más ACE-68 y el párrafo de alcance de §3). Se reescriben, y la
ampliación aplazada DTF-056 se divide: la preparación amplificada, los
duplicados de PCR y el sesgo de composición **configurable** salen de ella; el
perfil de sesgo GC **calibrado con datos reales**, los duplicados ópticos y el
marcado de duplicados en BAM/CRAM siguen aplazados. Se documenta como cambio
de contrato posterior al cierre de F4, con la misma forma que el plan de F5
usó el 13 de septiembre de 2026 para los techos de recursos.

**Qué depende de las repeticiones en tándem y qué no.** La amplificación, los
duplicados, el sesgo de composición y la sobrerrepresentación de los primeros
ciclos no dependen de nada de F6 y pueden implementarse antes. El tartamudeo
clonal depende de que existan loci de repetición declarados en la verdad
(hito F6.4, fichas D-01 a D-06 del expediente) y del tartamudeo sin PCR de la
parte C del expediente, del que reutiliza el modelo de paso y la aplicación
sobre la molécula.

**Qué queda por decidir.** Once decisiones, en la sección 11, con opciones y
recomendación. Las dos que más condicionan: cómo se acepta el sesgo de
composición sin romper la cota (rechazo con tope de intentos, recomendado,
frente a muestreo ponderado por ventanas) y si la cota de recursos se
**publica** en el manifiesto o además se **impone** como rechazo previo, lo
que F4 nunca ha hecho.

---

## 2. Hechos comprobados en el repositorio

Cada afirmación de esta sección se ha leído en el fichero citado el 19 de
septiembre de 2026. El diseño se apoya en ellas; si alguna dejara de ser
cierta, hay que revisar la parte que la usa.

### 2.1. El contrato cerrado de F4

La especificación declara la preparación sin PCR en estos sitios:

| Sitio | Texto literal o resumen |
|---|---|
| §3, párrafo de alcance de F4 | «La preparación de biblioteca será exclusivamente sin PCR; la preparación WGS amplificada y sus duplicados o sesgos quedan aplazados en DTF-056.» El mismo §3 dice también que «los duplicados de PCR de una biblioteca WGS amplificada sí pertenecen al alcance: son un artefacto de su preparación, no una secuenciación dirigida por PCR». |
| REQ-CONF-88 | «F4 admite exclusivamente preparación WGS sin PCR. Rechaza un perfil amplificado y cualquier tasa de duplicados PCR u ópticos; el modelo correspondiente queda aplazado en DTF-056.» |
| REQ-CONF-89 | «El perfil inicial declara explícitamente duplicados PCR no aplicables y ópticos a cero. F4 no permite sustituir esos valores.» |
| REQ-CONF-94/95 | GC neutro; rechaza curvas o intensidades no nulas; el perfil con sesgo GC calibrado queda en DTF-056. |
| REQ-TEC-07 | «F4 usa exclusivamente preparación WGS sin PCR y la registra como tal. Debe rechazar una preparación con amplificación e indicar DTF-056.» |
| REQ-TEC-08 | «El perfil WGS amplificado, sus duplicados PCR y sus sesgos propios quedan aplazados en DTF-056; F4 no acepta sus parámetros ni los representa como secuenciación dirigida por PCR.» |
| REQ-DAT-46 | «F4 no genera duplicados: PCR no aplica al perfil inicial y la tasa óptica es cero. Los modelos de duplicación quedan aplazados en DTF-056.» |
| REQ-DAT-47/48 | Duplicación de parejas completas y procedencia individual de duplicados, aplazadas en DTF-056. |
| REQ-DAT-52/53 | Sin sesgo GC; relación no neutra aplazada en DTF-056. |
| REQ-MAN-45/48, REQ-MAN-60 | El manifiesto declara PCR no aplicable, ópticos a cero, GC neutro y «la preparación sin PCR». |
| ACE-49, ACE-51, ACE-68 | Criterios de aceptación que repiten lo anterior. |
| DTF-056 | «Preparación amplificada, duplicados y sesgo GC no neutro […] No se inventan tasas, curvas ni modelos sin calibración aprobada.» |

En el código, el contrato está fijado en tres capas: el perfil incorporado
`gps-hiseqx-pcrfree-pe150-v1` es el único que existe y su tipo
`LibraryPreparation` tiene una sola variante, `PcrFree`
(`src/technical_profile.rs`, líneas 28–30); el evento efectivo devuelve
`pcr_duplicates_applicable() == false` y ópticos `0` como constantes
(`src/application.rs`, líneas 776–786); y el manifiesto de evento y la
configuración efectiva exigen literalmente `pcr.applicable == false`,
`gc.mode == "neutral"` y `preparation == "pcr_free"`
(`src/event_manifest.rs`, líneas 1516–1531; `schemas/manifest-event-v2.schema.json`
y `schemas/effective-config-v4.schema.json`, con `const`). El TOML general v2
(`schemas/config-general-v2.schema.json`, `$defs/technical_event`) es un
objeto cerrado con once claves; ninguna nombra la PCR, de modo que hoy una
clave de PCR se rechaza como campo desconocido, no como capacidad aplazada.

### 2.2. El recorrido de una pareja y dónde caben los pasos nuevos

El recorrido real de una pareja, en `WgsPairedEndSampler::sample_next_pair`
(`src/wgs_sampling.rs`, líneas 756–826), es: derivar la plantilla por su
ordinal; sortear inserto, copia haplotípica e inicio con claves
`fragment`; **materializar** la secuencia con `materialize_copy_interval`;
construir `PairedEndFragment` —la molécula—; emitir R1/R2 con
`emit_paired_end`; y decidir la quimera con `chimera::transform_paired_end`.
Después, en `src/execution.rs` (líneas 1270–1276), cada pareja pasa por
`error_quality::transform_paired_end`, y más tarde por contaminación y
demultiplexado.

**Se confirma lo que el encargo daba por cierto:** la formación de la
quimera vuelve a materializar el prefijo del fragmento desde la copia
haplotípica y lo compara byte a byte (`src/chimera.rs`, líneas 232–238);
si no coincide devuelve `ChimeraError::PairSourceMismatch`. Por tanto ningún
paso que modifique la secuencia de la molécula puede ir antes de la quimera.
La transformación de calidad y errores solo usa las secuencias de R1 y R2 y
las identidades (`src/error_quality.rs`, líneas 341–373), no las coordenadas
del fragmento. El sitio correcto para el tartamudeo, con o sin PCR, es
**después de la quimera y antes de calidad y errores**, exactamente como
concluyó la ficha D-13 del expediente.

### 2.3. El número de parejas es determinista

`resolve_wgs_paired_end_plan` (`src/wgs_sampling.rs`, líneas 218–279) deriva
el número de parejas de una biblioteca como
`round(cobertura × longitud_haploide / (2 × longitud_lectura))`, con empate
hacia arriba, sin ningún sorteo. Es el hecho que hace posible una cota cerrada
para todo lo que sigue: `N` se conoce antes de abrir el primer ámbito
aleatorio.

### 2.4. Memoria y pasadas del recorrido actual

`publish_technical_event` (`src/execution.rs`, líneas 647–1120) recorre cada
biblioteca **dos veces**: una primera pasada conserva solo una clave compacta
por pareja del evento entero, y una segunda, `RegeneratedPairs::build`
(líneas 1183–1290), regenera la biblioteca completa y retiene en un `BTreeMap`
**todas sus parejas observadas** antes de escribirlas. El pico de memoria de
F4 ya es, hoy, proporcional al número de parejas de la biblioteca más grande.
La escritura final sí es incremental (`FinalPairedFastqStream::write_pair`,
`src/fastq_output.rs`, línea 288). Este diseño no empeora ese orden de
memoria y hereda el factor dos de las pasadas.

### 2.5. F4 no tiene presupuesto de recursos

Las constantes `max_operation_units` y `max_working_memory_bytes` aparecen en
`src/application.rs` únicamente en la vía F5 (líneas 531, 1110–1114). La
hoja de ruta reserva «F9.4. Presupuesto de recursos» para estimar y rechazar
antes de simular, y el plan F4 §14 dice literalmente que «los presupuestos,
advertencias de recursos y la evaluación a escala permanecen en F9.4». El
nivel de rigor que el encargo pide —cota superior demostrada sin ejecutar
nada— es el de DTC-089, que hoy solo rige F5. La sección 6 da esa cota; la
sección 11 deja al usuario decidir si además se impone.

### 2.6. Cómo acota el proyecto un recuento aleatorio

DTC-089 y el commit `a1f8104` fijan la técnica: un recuento Poisson se
**pinza** en un intervalo declarado en vez de rechazarse con reintentos,
«porque el rechazo no acota el número de sorteos y el presupuesto de
operaciones sí». El mismo commit deja dicho, con honestidad, que el
muestreador Poisson por rechazo «no tiene cota superior de palabras» RNG y
que el presupuesto de operaciones no contabiliza esos sorteos. Este diseño
sigue ambas reglas: pinza todo recuento y declara qué no acota.

### 2.7. Catálogo aleatorio y primitivas disponibles

El catálogo técnico vigente es `f4-technical-v3` (DTC-086), definido en
`src/randomness.rs` (líneas 1333–1470) con etapas `fragment`, `chimera`,
`error`, `contamination_*` e `index_hop`, cada decisión con un ChaCha12
propio derivado de su clave y sin consumo acumulado (DTC-081). Las primitivas
normativas de F4 son `uniform-index-v1`, `bernoulli-v1` —que no consume para
probabilidad 0 ni 1—, `uniform-open-53-v1`, `uniform-closed-53-bit-v1`,
`partial-fisher-yates-v1` y `normal-truncated-box-muller-v1` (DTC-088). La
familia general aporta además `categorical-f64-v1` (DTC-072) y el único
recuento normativo `poisson-count-v1` (DTC-067). Este diseño no introduce
ninguna primitiva nueva; construye sus distribuciones con `categorical-f64-v1`
sobre pesos calculados de forma cerrada y con `bernoulli-v1`.

### 2.8. La calibración del tartamudeo con PCR ya está documentada

`design/calibracion-tartamudeo-0.1.md` recoge la tabla 1 de Gymrek (2016),
CC BY 4.0, con la fracción de lecturas afectadas por tartamudeo **con** y
**sin** PCR por longitud de motivo, y establece que entre protocolos cambia
sobre todo la frecuencia, no la forma de la distribución del tamaño del
error. La columna «con PCR», hasta ahora descartada, es la que ahora aplica.
Sus cifras son de genoma humano en Illumina y el documento advierte que no son
universales.

---

## 3. El modelo biológico, en llano

Esta sección explica qué fenómenos se simulan y con qué fidelidad. Distingue
en cada punto lo que está medido en la literatura de lo que es una elección de
diseño. Las fuentes de fuera del repositorio se citan de memoria y **deben
verificarse antes de adoptar de ellas cualquier cifra**; la sección 13 las
enumera con ese aviso.

### 3.1. Qué hace la amplificación

Una biblioteca sin PCR es una colección de moléculas originales de las que
cada una, si se secuencia, produce como mucho una pareja de lecturas. Una
biblioteca amplificada pasa esas moléculas por `c` ciclos de PCR: en cada
ciclo, cada molécula presente se copia con una probabilidad `e`, la
eficiencia. Tras `c` ciclos la molécula original está representada por un
número de copias que, en promedio, vale `(1+e)^c`, pero que varía mucho de
una molécula a otra. La secuenciación toma entonces `N` parejas de esa mezcla.

Tres consecuencias observables, y son las que este diseño produce:

1. **Duplicados de PCR.** Dos parejas que proceden de copias de la misma
   molécula original tienen el mismo fragmento —mismo inicio, misma longitud,
   misma copia haplotípica— y difieren solo en sus errores de secuenciación.
   La fracción de parejas duplicadas depende de cuántas parejas se secuencian
   (`N`) frente a cuántas moléculas distintas hay en la biblioteca (`C`, la
   complejidad). Depende del material de partida porque menos DNA de entrada
   significa menos moléculas distintas. **El número de ciclos no cambia por
   sí mismo esa fracción**; la cambia de forma indirecta, porque cuantos más
   ciclos, más desigual es la representación de las moléculas y más
   concentradas en pocas familias quedan las parejas.
2. **Sobrerrepresentación estocástica de los primeros ciclos.** Si en el
   primer ciclo una molécula se copia y su vecina no, la primera parte con el
   doble de representación para todo lo que sigue, y esa ventaja no se
   corrige después: los ciclos tardíos afectan a muchas copias a la vez y sus
   fluctuaciones se promedian. Es el proceso de ramificación de Galton–Watson
   aplicado a la PCR (Sun 1995; Weiss y von Haeseler 1995). La varianza
   relativa del número final de copias tiene forma cerrada y está dominada por
   los primeros ciclos; Best y col. (2015) midieron además que la eficiencia
   difiere entre moléculas, lo que ensancha esa desigualdad.
3. **Sesgo por composición de bases.** La eficiencia de copia depende del
   contenido en G+C del fragmento: los fragmentos muy pobres o muy ricos en
   GC se amplifican peor en cada ciclo y, compuesto sobre `c` ciclos, quedan
   muy infrarrepresentados. Aird y col. (2011) identificaron la PCR como la
   causa principal del sesgo GC de las bibliotecas Illumina, y Benjamini y
   Speed (2012) mostraron que el efecto se mide mejor sobre el GC del
   **fragmento** completo que sobre el de la lectura; la curva es unimodal,
   con un óptimo intermedio y caída en ambos extremos. Kozarewa y col. (2009)
   mostraron que eliminar la PCR reduce ese sesgo y los duplicados.

### 3.2. Qué es el tartamudeo clonal y por qué es distinto

El tartamudeo es el cambio del número de unidades observado en una repetición
en tándem por deslizamiento de la polimerasa al copiar el tracto. Sin PCR
ocurre en la amplificación en agrupación de la celda y en la síntesis de la
secuenciación, una vez por molécula secuenciada, y ese es el modelo de la
parte C del expediente.

Con PCR se añade una fuente cualitativamente distinta: el deslizamiento puede
ocurrir en **cualquier ciclo**, y una copia que resbala transmite su número
de unidades a **toda su descendencia**. Un resbalón en el ciclo 1 alcanza a
la mitad de la familia; uno en el último ciclo, a una sola copia. El
resultado es que una familia de duplicados puede mostrar dos o tres alelos
distintos en proporciones que imitan a un individuo heterocigoto o a una
mezcla, cuando en realidad procede de una sola molécula. Por eso el
tartamudeo con PCR **no es independiente por lectura**: está correlacionado
con la estructura de duplicados, y un modelo que sortee el tartamudeo por
lectura no lo reproduce.

Gymrek (2016) midió el efecto neto de todo lo anterior: la fracción de
lecturas afectadas con PCR es entre 3,5 y 6,3 veces mayor que sin PCR para
motivos de 1 a 5 pb (la cifra de 6 pb es ruidosa, como advierte el autor), y
la distribución del tamaño del error es parecida en ambos protocolos.

### 3.3. Lo que se deja fuera a propósito

- **Errores de polimerasa heredables.** Las polimerasas de PCR introducen
  sustituciones con tasas del orden de `10^-6` a `10^-5` por base y ciclo,
  también heredables dentro de la familia; en datos reales aparecen como
  variantes de baja frecuencia soportadas por duplicados. Es un modelo
  independiente del tartamudeo, con su propia calibración, y no lo pide el
  encargo. Se registra como ampliación futura (sección 12).
- **Duplicados ópticos y de agrupación.** No son PCR; siguen a cero.
- **Quimeras formadas durante la PCR y heredadas por la familia.** La quimera
  actual se decide por pareja e independientemente. Con PCR sería más fiel
  formarla una vez por familia en un ciclo dado; este diseño conserva la
  decisión por pareja y lo declara como simplificación.
- **Preferencia de la PCR por fragmentos cortos** (sesgo de longitud). No se
  modela; el inserto sigue la distribución del perfil.
- **Marcado de duplicados en BAM/CRAM.** Sigue en DTF-055/DTF-056; no hay
  mapeo en F4.

---

## 4. El problema duro: tres formas de acotar una amplificación

El proyecto exige que las operaciones y la memoria se acoten antes de
ejecutar. Una amplificación fiel es un proceso de ramificación: `C`
moléculas que tras `c` ciclos son del orden de `C × (1+e)^c`. Tres
alternativas:

| Alternativa | Qué hace | Coste | Qué reproduce | Veredicto |
|---|---|---|---|---|
| **A. Simular la población de moléculas** | Materializa cada molécula y cada copia ciclo a ciclo; secuencia `N` de entre todas | Operaciones y memoria `O(C × (1+e)^c)`: con `C ≈ N` y `c = 10`, mil veces `N`; exponencial en `c` | Todo, con fidelidad máxima | **Rechazada.** No admite cota razonable y casi todo lo que calcula no se observa nunca: la mayoría de las copias no se secuencian |
| **B. Muestrear hacia atrás desde las parejas** | Genera familias de amplificación una tras otra hasta cubrir las `N` parejas; para cada familia sortea su molécula original, su peso de amplificación, su tamaño y, solo si hace falta, la genealogía de sus `K` parejas | `O(N × constantes de contrato)`: cada familia cuesta a lo sumo un número fijo de operaciones y `Σ K = N` | Duplicados con su distribución de tamaños, sobrerrepresentación de los primeros ciclos, sesgo de composición y tartamudeo clonal correlacionado con la familia | **Elegida.** Solo calcula lo que se observa; la cota es lineal en `N` y cerrada |
| **C. Aproximación analítica de los efectos** | Aplica una fracción fija de duplicados copiando parejas ya observadas, un adelgazamiento por GC y una tasa de tartamudeo mayor por lectura | `O(N)` con constantes mínimas | Tasas medias correctas | **Rechazada como diseño único.** No produce la cola pesada de tamaños de familia ni el tartamudeo clonal, y copiar parejas observadas duplica también los errores de secuenciación, que en la realidad son independientes entre duplicados. Sus fórmulas cerradas sí se usan **dentro** de B para las distribuciones |

La elección B tiene una consecuencia que conviene decir con claridad: las
parejas duplicadas se producen **antes** de calidad y errores, porque lo que
comparten es la molécula —fragmento, resbalones de PCR—, no la secuenciación.
La alternativa C, que es lo que muchos simuladores hacen, se equivoca en esto.

Dentro de B, los tres objetos que en A serían un árbol exponencial se
reducen a formas cerradas:

- el **número de copias tras los primeros ciclos**, que es donde nace la
  desigualdad, se calcula exactamente por programación dinámica sobre a lo
  sumo `2^c_early` estados una vez por evento, y se sortea con un categórico;
- el **tamaño de familia** se sortea de una Poisson truncada y pinzada con
  un solo uniforme;
- la **genealogía** de las `K` parejas de una familia se sustituye por un
  árbol recursivo aleatorio de `K` nodos con un ciclo de divergencia por
  nodo, que basta para que los resbalones se hereden.

---

## 5. Diseño

### 5.1. Vocabulario

- **Molécula original:** fragmento de DNA de la biblioteca antes de
  amplificar: copia haplotípica, inicio, longitud de inserto y secuencia.
  Es lo que hoy construye `PairedEndFragment`.
- **Familia de amplificación:** conjunto de parejas publicadas que descienden
  de la misma molécula original. Tamaño `K ≥ 1`. Una familia de tamaño 1 no
  tiene duplicados.
- **Peso de amplificación `W`:** representación relativa de una molécula en
  la mezcla amplificada. Producto de tres factores: el número de copias tras
  los primeros ciclos (`Z`), el factor de composición (`G`) y, si se declara,
  la clase de eficiencia de la molécula.
- **Linaje:** camino de copias que lleva de la molécula original a una pareja
  concreta, con un paso por ciclo. Dos parejas de la misma familia comparten
  el linaje hasta su **ciclo de divergencia** `d` y son independientes
  después.
- **Resbalón:** deslizamiento de la polimerasa en un tracto en un paso del
  linaje; cambia el número de unidades en un paso entero del modelo de
  tartamudeo (parte C del expediente, `stutter-v1`).

### 5.2. Parámetros y su origen

| Parámetro | Símbolo | Dominio | Predeterminado | Origen |
|---|---|---|---|---|
| Preparación | — | `pcr_free` \| `pcr` | `pcr_free` | Contrato vigente |
| Ciclos de PCR | `c` | entero en `1..=C_MAX` (`C_MAX = 30`) | **sin valor**: obligatorio con `pcr` | Los protocolos Illumina con PCR usan del orden de 6 a 12 ciclos según el DNA de entrada; cifra de memoria, sin fuente en el repositorio, por eso no se fija |
| Eficiencia por ciclo | `e` | `f64` en `(0, 1]` | **sin valor**: obligatorio | Diseño; la literatura da valores altos, típicamente por encima de 0,8, para PCR de biblioteca |
| Ciclos explícitos | `c_early` | entero en `0..=E_MAX` (`E_MAX = 6`), `c_early ≤ c` | `min(3, c)` | Diseño; la varianza relativa del proceso de ramificación satura en pocos ciclos (fórmula en 5.5) |
| Parejas por molécula | `m` | `f64` en `[0, M_MAX]` (`M_MAX = 64`) | **sin valor**: obligatorio | `m = N / C`; es el cociente entre parejas secuenciadas y moléculas distintas, que fija la fracción de duplicados |
| Tamaño máximo de familia | `k_max` | entero en `1..=K_MAX` (`K_MAX = 256`) | `64` | Diseño; tope de pinzamiento |
| Curva de composición | `g(GC)` | tabla por tramos lineales de `(fracción GC, eficiencia relativa)` con eficiencia en `(0, 1]` y nudos crecientes en `[0, 1]` | neutra: `[[0, 1], [1, 1]]` | Diseño de la forma; ninguna curva calibrada se incorpora (queda en DTF-056) |
| Intentos de aceptación | `A` | entero en `1..=A_MAX` (`A_MAX = 16`) | `8` | Diseño; tope de rechazo de 5.4 |
| Clases de eficiencia | `{(e_i, w_i)}` | lista de a lo sumo `Q_MAX = 4` pares, pesos normalizados DTC-049 | una sola clase `(e, 1)` | Best y col. (2015) para el fenómeno; sin valores |
| Tasa de resbalón por ciclo, por longitud de motivo | `p_ℓ` | `[0, 1]` | derivada de la columna con PCR de la calibración (5.7) | Gymrek (2016), CC BY 4.0, humano/Illumina |
| Distribución del ciclo de divergencia | `q` | derivada de `e` | `q = 1/(1+e)` | Diseño, aproximación del coalescente del proceso de ramificación |

Las constantes en mayúsculas son **constantes de contrato**: viven en una
sola constante de Rust cada una, atada a los `maximum` de los esquemas por
una prueba juez, como hizo el plan F5 §11.4 con los techos de recursos.
Sirven para que la cota de la sección 6 sea cerrada.

### 5.3. El muestreo por familias, paso a paso

El muestreador PE de una biblioteca, `WgsPairedEndSampler`, conserva su
contrato externo —un iterador que produce exactamente `N` parejas con
ordinales de plantilla `0..N-1`— y añade un estado de familia entre llamadas.
Para cada familia `f` con ordinal de familia igual al **ordinal de plantilla
de su primera pareja**, `o_f`:

1. **Clase de eficiencia.** Si hay más de una clase, un
   `categorical-f64-v1` sobre los pesos declarados elige `e_f`. Con una sola
   clase no se consume nada.
2. **Copias tras los primeros ciclos, `Z_f`.** Un `categorical-f64-v1` sobre
   la distribución **sesgada por tamaño** del número de copias tras `c_early`
   ciclos con eficiencia `e_f` (5.5). Con `c_early = 0` o `m = 0`, `Z_f = 1`
   sin consumo.
3. **Molécula original, con aceptación por composición.** Para el intento
   `a = 0`, se sortean inserto, copia haplotípica e inicio con las claves
   `fragment` **de hoy**, indexadas por `o_f`, y se materializa la secuencia.
   Se calcula la fracción GC del fragmento sobre A/C/G/T (las `N` no cuentan
   ni en numerador ni en denominador; un fragmento sin bases A/C/G/T recibe
   `g = 1`) y el factor `G_f = ((1 + e_f·g(GC)) / (1 + e_f))^c`. Un
   `bernoulli-v1` con probabilidad `G_f` acepta el fragmento; con curva neutra
   `G_f = 1` y no se consume nada. Si rechaza y `a + 1 < A`, se repite el
   sorteo del fragmento con las mismas claves más un componente `attempt`;
   si se agota `A`, se **acepta el último intento** y se incrementa el
   recuento agregado `gc_attempts_exhausted_families`. Este es el
   pinzamiento de esta etapa: la cota es `A` materializaciones por familia,
   y la distorsión que introduce se publica, no se oculta.
4. **Tamaño de familia, `K_f`.** Con `λ_f = m × (Z_f / E[Z]) × G_f`, donde
   `E[Z] = (1+e_f)^c_early`, un `categorical-f64-v1` sobre la Poisson
   truncada a cero y pinzada en `k_max` (5.6). Si `λ_f = 0`, `K_f = 1` sin
   consumo. Si la familia excede las parejas que quedan por producir
   (`o_f + K_f > N`), se pinza a `N − o_f` y se cuenta en
   `family_truncated_at_end` (a lo sumo una vez por biblioteca).
5. **Parejas de la familia.** Para `i` en `0..K_f`, con ordinal de plantilla
   `o = o_f + i`:
   1. la plantilla y la identidad previa al salto se derivan **como hoy**
      desde `o`; cada duplicado es una pareja con identidad propia, lo que
      exige REQ-DAT-18;
   2. se emite R1/R2 desde la molécula original sin modificar;
   3. la **quimera** se decide como hoy, por pareja, con las claves `chimera`
      indexadas por `o`, sobre la molécula sin modificar: la comparación byte a
      byte se conserva;
   4. el **tartamudeo con PCR** (5.7) calcula el paso neto de la pareja a
      partir del linaje de la familia y lo aplica sobre la molécula —quimérica
      si lo es— antes de reemitir R1/R2; después se aplica el tartamudeo sin
      PCR de la parte C del expediente, con las claves `stutter` indexadas por
      `o`;
   5. calidad y errores, contaminación y demultiplexado siguen sin cambios
      aguas abajo, por pareja.

El iterador devuelve las parejas en orden de ordinal, de modo que las
parejas de una familia son consecutivas en la generación. Como el orden de
escritura del FASTQ final es por identificador de instancia (`CleanLibraryPairs`
ordena sus claves, `src/execution.rs`, comentario en la línea 1183 y
siguientes), los duplicados **no** aparecen adyacentes en el fichero
publicado, como en un FASTQ real.

### 5.4. Por qué la composición tiene que actuar al elegir la molécula

Una tentación es meter el sesgo de composición solo en `λ_f`, es decir, en
el tamaño de familia, y evitar así el rechazo. No funciona: con `m` pequeño
—el caso normal, con pocos duplicados— casi todas las familias tienen
tamaño 1 tanto en regiones ricas como pobres en GC, y el número de parejas
que caen en cada región queda casi igual. La infrarrepresentación real de una
región pobre en GC viene sobre todo de las moléculas que **no llegan a
observarse** (`K = 0`), y en un modelo que produce familias observadas eso
solo puede expresarse al **elegir** qué molécula da lugar a la familia, con
probabilidad proporcional a su peso. Por eso la aceptación `bernoulli-v1(G_f)`
está en el paso 3 y no en el 4.

Alternativa evaluada: muestreo ponderado por ventanas. Se recorrería una vez
cada copia haplotípica para calcular el GC por ventanas de tamaño fijo,
se guardarían sumas prefijo de `g` por ventana y cada familia elegiría
ventana por búsqueda dicotómica y luego inicio uniforme. Ventajas: sin
rechazo ni tope; coste `O(log ventanas)` por familia. Inconvenientes: el
sesgo se aplica al GC de la ventana, no al del fragmento, que es lo que la
literatura dice que importa; exige un consumidor aleatorio nuevo con sus
vectores; y su memoria, `8 bytes × Σ_copias ⌈longitud/ventana⌉`, crece con
el genoma (para un tomate diploide de 800 Mb y ventanas de 100 pb, unos
128 MB). Se recomienda el rechazo con tope para la primera versión y se deja
la alternativa registrada (decisión 2 de la sección 11).

### 5.5. Los primeros ciclos: forma cerrada de la desigualdad

Con eficiencia `e` y `z` copias presentes, el siguiente ciclo produce
`z + Binomial(z, e)` copias. Partiendo de una copia, la distribución de `Z`
tras `c_early` ciclos se calcula **exactamente** por programación dinámica
sobre los estados `1..=2^c_early`: coste `O(c_early × 4^c_early)` una vez por
evento y clase, que para `E_MAX = 6` son unas 25 000 operaciones. Las
probabilidades binomiales se acumulan por multiplicaciones sucesivas de `e` y
`1−e` en orden canónico, sin `libm`. La distribución que se sortea es la
**sesgada por tamaño**, `p*(z) ∝ z × p(z)`, porque una familia observada es
más probable cuanto más copias tenga su molécula; es la razón por la que
se sortea `Z` antes que `K`.

Por qué bastan pocos ciclos explícitos: para el proceso de ramificación con
descendencia 1 o 2, la varianza relativa del número final de copias vale
`(1−e)/(1+e) × (1 − (1+e)^-c)`, que con `e = 0,9` converge a 0,053 y con
`c = 3` ya alcanza el 85 % de ese límite; los ciclos posteriores apenas
añaden desigualdad y se toman como el factor determinista `(1+e)^(c−c_early)`.
Que `c_early` sea configurable hasta `E_MAX` permite al usuario comprobarlo,
y la cota de la sección 6 se calcula con `2^c_early`, no con `2^c`.

### 5.6. Tamaño de familia y fracción de duplicados

`K_f` sigue una Poisson de media `λ_f` **condicionada a `K ≥ 1`** y pinzada
en `k_max`: para `1 ≤ k < k_max`, `p_k = (λ^k / k!) × e^-λ / (1 − e^-λ)`, y
`p_(k_max)` recoge toda la masa restante. Los pesos se calculan por
multiplicación sucesiva (`p_1 = λ`, `p_(k+1) = p_k × λ/(k+1)`), con `e^-λ`
de `libm` como en DTC-067, y se normalizan con la suma compensada de DTC-049.
Un solo `categorical-f64-v1`, es decir, un solo uniforme abierto por
familia. No se usa `poisson-count-v1` porque su rechazo no acota los sorteos y
su cola no está pinzada; la decisión técnica nueva debe decirlo para que nadie
lea DTC-067 como una prohibición.

Fracción de duplicados. Para una biblioteca homogénea (`Z = 1`, `G = 1`), la
fracción esperada de parejas duplicadas es `1 − (1 − e^-m)/m`: con `m = 0,02`
vale el 1 %; con `m = 0,1`, el 4,8 %; con `m = 0,5`, el 21 %. La
configuración efectiva publica ese valor como `expected_duplicate_fraction`
para que el usuario elija `m` con criterio, y el manifiesto publica la
fracción **realizada**, `1 − familias/N`. Con desigualdad (`Z`, `G`, clases)
la realizada es mayor que la homogénea, y esa diferencia es precisamente el
efecto de los primeros ciclos y de la composición.

### 5.7. Tartamudeo clonal: la genealogía mínima

Depende de la parte C del expediente: reutiliza el modelo `stutter-v1` (paso
entero en unidades, pesos de contracción y expansión, mínimo de unidades
contenidas, guardas de longitud) y la operación «aplicar un paso `k` al tracto
contenido en una molécula». Lo que añade es **de dónde sale `k`**.

**Linaje de una pareja.** Cada pareja `i` de una familia tiene un linaje de
`c` pasos, uno por ciclo. En cada paso y para cada tracto declarado contenido
en la molécula, hay un resbalón con probabilidad `p_ℓ` (por longitud de
motivo `ℓ`), decidido con `bernoulli-v1` y, si ocurre, un paso sorteado con el
categórico de `stutter-v1`. El paso neto de la pareja es la suma de sus
resbalones. Con `p_ℓ = 0` no se consume nada.

**Herencia.** La pareja `0` de la familia tiene un linaje propio completo.
La pareja `i ≥ 1` elige un **progenitor** `j` uniforme en `0..i`
(`uniform-index-v1(i)`) y un **ciclo de divergencia** `d` en `0..=c` con
`P(d = j) = (1−q) × q^j` para `j < c` y el resto en `d = c`, con
`q = 1/(1+e_f)`; hereda los resbalones del progenitor en los ciclos
`1..=d` y sortea los suyos en `d+1..=c`. Con `d = 0` no comparte nada: es
el caso, muy frecuente con eficiencia alta (`P(d=0) = e/(1+e) ≈ 0,47`), de
dos parejas que descienden de productos distintos del primer ciclo. Así, un
resbalón en el ciclo 1 del linaje de una pareja lo heredan todas las que
divergen de ella —o de sus descendientes— en `d ≥ 1`, y uno en el último
ciclo no lo hereda nadie. Es la propiedad que el encargo pide y la que la
prueba 4 de la sección 10 discrimina.

Esta genealogía es una **aproximación declarada** del coalescente del
proceso de ramificación: la forma geométrica de `d` refleja que dos
descendientes de una amplificación eficiente divergen casi siempre pronto, y
el progenitor uniforme entre las parejas anteriores produce un árbol
recursivo aleatorio. Kebschull y Zador (2015) estudiaron estas distorsiones
por simulación del árbol completo; no conozco una forma cerrada del
coalescente exacto y no se afirma tenerla.

**Calibración.** Gymrek mide fracciones de lecturas afectadas de extremo a
extremo, no tasas por ciclo. Se derivan así: sea `r_con` la fracción con
PCR y `r_sin` la fracción sin PCR de la tabla de la calibración, para un
motivo dado. La parte atribuible a la PCR es `r_pcr = 1 − (1−r_con)/(1−r_sin)`,
y la tasa por paso de linaje que la reproduce en `c` ciclos es
`p_ℓ = 1 − (1 − r_pcr)^(1/c)`, con la potencia de `libm`. Con `c = 10`, a
título ilustrativo y redondeado a tres cifras: motivo 1 pb, `p = 0,0158`;
2 pb, `0,00276`; 3 pb, `0,000815`; 4 pb, `0,000532`; 5 pb, `0,000593`; 6 pb,
`0,00131` (el menos fiable, como advierte la fuente). El simulador calcula
`p_ℓ` desde `r_con`, `r_sin` y `c` en la configuración efectiva, publica la
fórmula y los tres insumos, y admite sustituir `p_ℓ` directamente. Si
`r_con ≤ r_sin`, `p_ℓ = 0`. La dependencia con la longitud del tracto, sin
fuente abierta, queda como en la calibración: configurable, sin valor
derivado.

**Las dos columnas de Gymrek son totales, no componentes.** Conviene dejarlo
escrito porque invita al error contrario. Los dos conjuntos que comparó pasaron
por las **mismas** células, la **misma** amplificación en puente y la **misma**
síntesis; lo único que los distingue es el paso de PCR de preparación. Por
tanto la columna sin PCR ya contiene el mosaicismo somático de la muestra, el
deslizamiento del puente y el de la síntesis, y la columna con PCR contiene
todo eso **más** la PCR de biblioteca. **Sumar el puente y la síntesis al
modelo de PCR sería contarlos dos veces.** Por decisión del proyecto, el
mosaicismo somático se queda dentro de la tasa sin PCR y **no** se modela
aparte: hacerlo obligaría a descontarlo de esta cifra, que no permite
separarlo.

**El caso con PCR solo cierra si el caso sin PCR está implementado.** La
descomposición de arriba despeja **únicamente el incremento** de la PCR de
biblioteca. Si se implanta la preparación amplificada sin la parte C, las
bibliotecas con PCR se quedan en ese incremento y **subestiman** el
tartamudeo: un 14 % menos en homopolímeros, 28 % en dinucleótidos, 26 % en
trinucleótidos, 23 % en tetranucleótidos y 20 % en pentanucleótidos. Con las
dos partes, la composición devuelve exactamente la fracción medida
(`1 − (1−r_sin)(1−r_pcr) = r_con`). Es una dependencia de orden, no una
preferencia: figura en la sección 7.

**Orden con la quimera y con el tartamudeo sin PCR.** Sobre la molécula ya
quimérica si lo es, primero el paso neto de PCR y después la decisión sin PCR
(cluster y síntesis), que es por pareja y sigue el diseño de la parte C. Los
tractos se procesan en orden canónico de inicio sobre la molécula, cada uno
con su `locus_id`. Las guardas de viabilidad (al menos una unidad, molécula
no más corta que la lectura) se aplican al paso **neto acumulado**: un
resbalón cuyo paso dejaría el tracto inviable se excluye del categórico
antes de sortear, sin reintento, como en D-14.

### 5.8. Identidades y ámbitos aleatorios

Ningún identificador nuevo: cada pareja de una familia es una plantilla con
su ordinal, su `TemplateId`, su `PreHopInstanceId` propio y su
`TechnicalPairId` final, como hoy (DTC-087). El manifiesto **no** publica
qué parejas forman familia (REQ-DAT-48 sigue aplazado: sin procedencia por
pareja); un validador puede rederivarlo desde la semilla.

Ámbitos, en una versión nueva del catálogo técnico (la que corresponda al
implementar: `f4-technical-v4` si se hace junto con `stutter`, la siguiente
si después), con el patrón literal de DTC-086
`[event_id, library_id, stage, family_ordinal, decision, …]` y etapa
`amplification`:

| Decisión | Componentes tras `decision` | Consumidor | Cuándo no consume |
|---|---|---|---|
| `efficiency_class` | — | `categorical-f64-v1` | una sola clase |
| `early_copies` | — | `categorical-f64-v1` sobre la DP sesgada por tamaño | `c_early = 0` o `m = 0` |
| `gc_accept` | `attempt` | `bernoulli-v1(G_f)` | `G_f = 1` |
| `family_size` | — | `categorical-f64-v1` sobre la Poisson truncada y pinzada | `λ_f = 0` |
| `lineage_parent` | `copy_index` | `uniform-index-v1(i)` | `i = 1` (un solo candidato) |
| `lineage_divergence` | `copy_index` | `categorical-f64-v1` sobre `P(d)` | `c = 0` no es admisible; siempre consume si `K ≥ 2` y hay tractos |
| `lineage_slip` | `copy_index`, `cycle`, `locus_id` | `bernoulli-v1(p_ℓ)` | `p_ℓ = 0` |
| `lineage_step` | `copy_index`, `cycle`, `locus_id` | categórico de `stutter-v1` | solo si hubo resbalón |

Los sorteos de `fragment` de los intentos `a ≥ 1` añaden un componente
`attempt` tras `decision`, igual que `second_source` añade `breakpoint`; el
intento `0` usa las claves de hoy sin cambio. Las decisiones de linaje se
indexan por `locus_id`, no por posición, por la misma razón que D-11:
añadir un locus no desplaza las decisiones de los demás. Cada clave abre su
ChaCha12 y no hay consumo acumulado (DTC-081), de modo que activar la PCR
no puede mover ninguna decisión de las etapas existentes.

### 5.9. Configuración

Propuesta de forma, en el evento técnico del TOML general:

```toml
[technical_event]
library_preparation = "pcr"          # predeterminado "pcr_free"

[technical_event.pcr]
cycles = 10
efficiency = 0.9
early_cycles = 3
mean_pairs_per_molecule = 0.1
max_family_size = 64
gc_max_attempts = 8
gc_relative_efficiency = [[0.0, 0.3], [0.3, 0.9], [0.45, 1.0], [0.6, 0.9], [1.0, 0.3]]
# efficiency_classes = [{ efficiency = 0.95, weight = 0.7 }, { efficiency = 0.6, weight = 0.3 }]

[technical_event.pcr.stutter]           # solo tiene efecto con loci de repetición declarados
cycles_source = "profile"               # las tasas por ciclo se derivan de la calibración con PCR
# per_cycle_slip = { "1" = 0.0158, "2" = 0.00276 }   # sustitución directa por longitud de motivo
```

Reglas: la tabla `[technical_event.pcr]` solo se admite con
`library_preparation = "pcr"` y es obligatoria en ese caso; `cycles`,
`efficiency` y `mean_pairs_per_molecule` no tienen predeterminado; toda
clave desconocida, todo valor fuera de dominio y toda tabla presente con
`pcr_free` se rechazan antes de simular, con parámetro y valor. La interfaz
Marimo añade los controles correspondientes y muestra la fracción esperada de
duplicados calculada desde `m` (regla de cierre Marimo de la hoja de ruta).

### 5.10. Perfil, configuración efectiva y manifiesto

**Perfil.** Un recurso hermano `resources/f4/gps-hiseqx-pcr-pe150-v1.json`,
idéntico al vigente en lectura, inserto, calidad y error, con
`"preparation": "pcr"` y un bloque `amplification` con los dominios y las
constantes de contrato. Mantiene la regla de REQ-CONF-82 —un perfil
identificado y versionado por tecnología— y deja intacto el recurso PCR-free
y su SHA-256. `LibraryPreparation` gana la variante `Pcr`. Alternativa: un
solo recurso con ambas preparaciones; se descarta porque el nombre
`pcrfree` del vigente es literal en esquemas y manifiestos publicados y su
hash está congelado.

**Configuración efectiva.** `technical_profile.preparation` pasa a admitir
`pcr_free | pcr`; `technical_artifacts.pcr` pasa de `{applicable: false}` a
`{applicable: bool, cycles, efficiency, early_cycles, mean_pairs_per_molecule,
max_family_size, gc_max_attempts, gc_relative_efficiency, efficiency_classes,
expected_duplicate_fraction, stutter: {per_cycle_slip por motivo, fuente}}`,
con `applicable: false` y ningún otro campo cuando es `pcr_free`. `gc.mode`
admite `neutral | curve`. Los campos nuevos no entran en `required`, y con
`pcr_free` el documento es byte a byte el de hoy.

**Manifiesto de evento.** `artifacts.pcr` añade, con `applicable: true`, los
recuentos agregados: `families`, `pairs`, `duplicate_pairs`,
`realized_duplicate_fraction`, `family_size_histogram` (índices `1..=k_max`),
`family_truncated_at_end`, `gc_attempts_total`,
`gc_attempts_exhausted_families`, `lineage_slips_by_locus` y
`lineage_inherited_slips_by_locus`. La cobertura sigue el desglose de
REQ-WGS-11: la básica antes de artefactos no cambia (`N` parejas); se añade
`unique_molecule_x`, la cobertura de moléculas distintas, que es la que un
analista vería tras eliminar duplicados. Nada por pareja (REQ-DAT-56).

**Versiones de esquema.** El plan F6 §2 promete no reabrir esquemas
publicados y crear versiones nuevas (general v4, efectivo v6, manifiestos).
La PCR entra en la versión que esté abierta cuando se implemente; si se
implementa antes que F6.1.1, abre ella esas versiones y F6 las completa
dentro de la misma versión mientras no estén publicadas externamente, que es
la regla que el plan F5 §1 usó para `resource_preflight`. Decisión 8 de la
sección 11.

### 5.11. Contaminación, demultiplexado y validación

Sin cambios de contrato: una pareja duplicada es una pareja más para la
contaminación biológica y el *index hopping*; una donante puede aportar un
duplicado a una receptora. La validación estructural del FASTQ final
(`src/fastq_validation.rs`, `validate_final_paired_fastq`) no cambia. La
validación de cierre añade el oráculo de la sección 10.

---

## 6. Cota de recursos

La cota se deriva de la configuración y de la geometría ya conocida antes de
abrir el primer ámbito aleatorio, con aritmética entera comprobada y nombres
dimensionalmente distintos, al estilo de `f5-resource-accounting-v9`. Se
expresa **por encima del coste que F4 ya tiene** para producir `N` parejas
(`X_pair`: emisión, quimera, calidad y errores por pareja), que este diseño no
modifica.

### 6.1. Magnitudes

- `N`: parejas de la biblioteca (2.3), determinista.
- `P`: número de copias haplotípicas del individuo; `L_max = 450`, inserto
  máximo del perfil.
- Constantes de contrato: `C_MAX = 30`, `E_MAX = 6`, `K_MAX = 256`,
  `A_MAX = 16`, `Q_MAX = 4`.
- Valores configurados, todos `≤` su constante: `c`, `c_early`, `k_max`, `A`,
  `Q` (clases).
- `T`: número de loci de repetición declarados que caben, total o
  parcialmente, en un inserto; cota grosera y válida: `T ≤ min(L_max, loci
  declarados)`.

### 6.2. Operaciones por familia

Cada familia cuesta a lo sumo:

- selección de clase y de `Z`: `Q + 2^c_early` unidades (recorrido lineal
  del categórico);
- aceptación por composición: `A × (P + 3·L_max)` unidades —enumerar
  candidatas, materializar el inserto, contar GC y emitir— más `A` sorteos
  Bernoulli;
- tamaño de familia: `k_max` unidades para construir y recorrer los pesos;
- por cada una de sus `K` parejas: `X_pair`, más el linaje,
  `c × T × 2` unidades (resbalón y, como mucho, un paso por ciclo y tracto),
  más `K` unidades para el progenitor y `c + 1` para la divergencia.

Como `Σ_f K_f = N` y el número de familias `F ≤ N`, la cota total es

```
ops_pcr(N) ≤ N × [ Q + 2^c_early + A × (P + 3·L_max) + k_max + 2·c·T + c + 1 + k_max ]
             + N × X_pair
             + Q × c_early × 4^c_early                     (DP de los primeros ciclos)
```

y, sustituyendo por las constantes de contrato, una cota que solo depende de
`N`, `P` y `T`:

```
ops_pcr(N) ≤ N × [ 4 + 64 + 16 × (P + 1350) + 256 + 60·T + 31 + 256 ] + N × X_pair + 4 × 6 × 4096
          =  N × [ 22 211 + 16·P + 60·T ] + N × X_pair + 98 304
```

Como el recorrido F4 regenera cada biblioteca dos veces (2.4), el coste
efectivo es el doble; la cota lo declara en vez de suponer una pasada.

### 6.3. Memoria

Por encima del pico actual de F4 (2.4), el estado añadido es constante:

```
mem_pcr ≤ L_max                                (molécula original)
        + 8 × (Q × 2^c_early + k_max + c + 1)  (tablas de pesos)
        + K_MAX × T × 8                        (pasos netos por pareja y tracto de la familia en curso)
        + 8 × (k_max + 1)                      (histograma)
        ≤ 450 + 8 × (256 + 256 + 31) + 256 × T × 8 + 2056
```

Es independiente de `N` y del genoma. Las tablas de pesos de la Poisson
truncada se reconstruyen por familia sobre el mismo búfer.

### 6.4. Lo que no se acota, dicho expresamente

- El número de **palabras RNG** de los sorteos por rechazo que F4 ya usa —la
  normal truncada del inserto y `uniform-index-v1`— no tiene cota superior
  en el peor caso; DTC-088 y el commit `a1f8104` ya lo reconocen y no lo
  contabilizan. Este diseño no añade ningún sorteo por rechazo: todas sus
  distribuciones son categóricos de una palabra o Bernoulli. Sí multiplica
  por `A` como máximo los sorteos de inserto de una familia.
- La cota usa `T ≤ min(L_max, loci)`; un usuario con miles de loci de
  repetición declarados obtiene una cota holgada pero válida.

### 6.5. Cómo se demuestra y cómo se impone

**Demostración.** Un contador de unidades instrumentado en el muestreador PE
(solo en pruebas) cuenta las unidades reales; una prueba ejecuta un evento
con `A`, `k_max`, `c_early` y `T` en sus máximos de contrato y comprueba
`unidades ≤ ops_pcr(N)`, con control positivo de que el contador contó, y un
**mutante** que omita el término `A × (P + 3·L_max)` de la fórmula debe
hacer fallar la prueba en un fixture donde la curva rechace casi siempre
(por ejemplo `g ≡ 0,01` fuera de un óptimo estrecho y `A = 16`). Sin ese
fixture saturado, la prueba no discrimina: con curva neutra `A` no se
ejercita.

**Imposición.** Hoy F4 no rechaza por presupuesto (2.5). Dos opciones,
decisión 9 de la sección 11: (a) publicar `ops_pcr(N)` y `mem_pcr` en el
manifiesto de evento como `resource_bound`, sin rechazo, y dejar la
imposición para F9.4; (b) añadir a F4 una admisión con `max_operation_units`
reutilizando las constantes de DTC-089, lo que reabre contratos de F5
(«subir de techo reabre DTC-089», plan F6 §2). Se recomienda (a).

---

## 7. Dónde encaja cada pieza y de qué depende

| Pieza | Depende de | Puede ir antes de F6.4 |
|---|---|---|
| Perfil `gps-hiseqx-pcr-pe150-v1`, variante `Pcr`, configuración y esquemas | Nada de F6 | **Sí** |
| Muestreo por familias, `Z`, `K`, aceptación por composición | Nada de F6 | **Sí** |
| Recuentos del manifiesto, fracción esperada y realizada de duplicados, `unique_molecule_x` | Lo anterior | **Sí** |
| Catálogo `amplification` (sin las decisiones `lineage_*`) | Lo anterior | **Sí** |
| Oráculo de duplicados por bytes (10, prueba 3) | Lo anterior | **Sí** |
| Genealogía y `lineage_*` | Loci de repetición en la verdad (F6.4.1–F6.4.3, fichas D-01 a D-06) y el tartamudeo sin PCR (parte C, D-11 a D-16) | **No** |
| Derivación de `p_ℓ` desde la calibración | Lo anterior | No |

Consecuencia práctica: el trabajo se puede cortar en dos entregas. La primera
—amplificación, duplicados, composición, primeros ciclos— cierra el cambio
de contrato de F4 y es verificable por sí sola. La segunda añade el
tartamudeo clonal cuando F6.4 y el tartamudeo sin PCR existan, y solo entonces
sube el catálogo técnico con las decisiones `lineage_*`. Si se hace la
primera entrega sin la segunda, el manifiesto declara
`stutter: not_applicable_no_repeat_loci` y el requisito lo dice.

---

## 8. Orden normativo de F4, reescrito

El párrafo de §4.4 de la especificación pasaría a:

> F4 sigue este orden normativo: calcula las parejas de cobertura; sortea la
> molécula original de cada familia de amplificación —una familia por
> pareja cuando la preparación es sin PCR— y, con preparación amplificada,
> su tamaño; para cada pareja de la familia, forma una quimera cuando
> corresponda, aplica el tartamudeo cuando corresponda —primero el heredado
> por el linaje de PCR, después el de agrupación y síntesis—, realiza calidad
> y errores de R1/R2, aplica contaminación biológica y finalmente
> demultiplexa.

REQ-DAT-17 enumera el mismo orden y cambia igual.

---

## 9. Cambio de contrato: qué se modifica y cómo queda

Ninguna modificación está aplicada. Las redacciones son propuestas. La
numeración nueva se elige **después** de la que el expediente de repeticiones
propone para el tartamudeo (REQ-DAT-74 a 80, REQ-CONF-131 a 133, REQ-MAN-93,
REQ-VAL-57, ACE-118, y una decisión técnica DTC-090); si aquella cambia,
esta se desplaza.

### 9.1. Requisitos que cambian

| Requisito | Redacción propuesta |
|---|---|
| **§3, alcance de F4** | «La preparación de biblioteca predeterminada es sin PCR. F4 admite además, por elección expresa, una preparación WGS amplificada por PCR con sus duplicados, su sobrerrepresentación estocástica, su sesgo de composición configurable y su tartamudeo heredable; el perfil de sesgo GC calibrado con datos reales, los duplicados ópticos y el marcado de duplicados quedan aplazados en DTF-056 y DTF-055.» |
| **REQ-CONF-88** | «F4 admite dos preparaciones WGS, `pcr_free`, predeterminada, y `pcr`. Con `pcr_free` rechaza cualquier parámetro de amplificación, de duplicados o de composición no neutra. Con `pcr` exige ciclos, eficiencia y parejas por molécula, y admite los parámetros de REQ-CONF-134 a 138. Las tasas de duplicados ópticos siguen rechazadas.» |
| **REQ-CONF-89** | «Cada perfil declara su preparación. El perfil `pcr_free` declara duplicados PCR no aplicables; el perfil `pcr` declara los dominios y constantes de contrato de su amplificación. Los duplicados ópticos valen cero en ambos y no pueden sustituirse.» |
| **REQ-CONF-94** | «Con `pcr_free`, F4 usa una relación GC neutra y rechaza curvas o intensidades no nulas. Con `pcr`, admite una curva por tramos lineales de eficiencia relativa por fracción GC del fragmento, neutra por omisión. Un perfil con curva calibrada queda aplazado en DTF-056.» |
| **REQ-CONF-95** | «La composición actúa exclusivamente al elegir la molécula original de cada familia y con un número máximo de intentos declarado; no existen excepciones regionales ni combinación con la cobertura objetivo, que conserva su valor.» |
| **REQ-TEC-07** | «F4 registra la preparación de cada biblioteca, `pcr_free` o `pcr`. Con `pcr_free` rechaza todo parámetro de amplificación e indica DTF-056 solo para lo que sigue aplazado.» |
| **REQ-TEC-08** | «El perfil WGS amplificado `pcr` es un perfil de preparación de biblioteca WGS; F4 no lo representa como secuenciación dirigida por PCR, que sigue en DTF-002.» |
| **REQ-DAT-46** | «Con `pcr_free`, F4 no genera duplicados. Con `pcr`, genera duplicados de PCR como parejas de una misma familia de amplificación, antes de calidad y errores; la tasa óptica es cero en ambos casos.» |
| **REQ-DAT-47** | «Una pareja duplicada es una pareja completa con identidad técnica propia, mismo fragmento que su familia y errores, calidades, quimera y contaminación propios.» |
| **REQ-DAT-52** | «Con `pcr`, la probabilidad de que un fragmento origine una familia se multiplica por `((1 + e·g(GC)) / (1 + e))^c`; con `pcr_free` toda pareja elegible conserva la probabilidad del muestreo WGS básico.» |
| **REQ-DAT-53** | «La cobertura objetivo no se corrige por composición; el manifiesto registra la cobertura básica, la de moléculas distintas y la final.» |
| **REQ-MAN-45** | «El manifiesto declara la preparación y, con `pcr`, los parámetros efectivos y los recuentos agregados de familias, duplicados, histograma de tamaños, intentos de composición y resbalones heredados; ópticos a cero.» |
| **REQ-MAN-48** | «El manifiesto declara `neutral` o `curve` con la curva efectiva.» |
| **REQ-MAN-60** | «…la preparación de la biblioteca, el objetivo original de cobertura, el número de parejas derivado y realizado, la cobertura de moléculas distintas…» |
| **ACE-49** | «Con `pcr_free` F4 declara PCR no aplicable y no genera duplicados. Con `pcr`, cada duplicado comparte fragmento con su familia y difiere en errores; la fracción realizada se publica y coincide con la que un oráculo cuenta sobre los FASTQ con errores desactivados; ópticos a cero.» |
| **ACE-51** | «Con `pcr`, una curva no neutra reduce las familias en las regiones de composición desfavorable en la proporción configurada, dentro del tope de intentos declarado.» |
| **ACE-68** | «F4 acepta exclusivamente WGS de DNA genómico con Illumina PE, en preparación `pcr_free` o `pcr`. Rechaza …» (resto igual). |

### 9.2. Requisitos nuevos

- **REQ-CONF-134.** La configuración del evento admite `library_preparation`
  con valores `pcr_free`, predeterminado, y `pcr`.
- **REQ-CONF-135.** Con `pcr`, la tabla de amplificación debe declarar
  ciclos, eficiencia por ciclo y parejas por molécula, sin predeterminados, y
  puede declarar ciclos explícitos, tamaño máximo de familia, intentos de
  composición, curva de composición y clases de eficiencia dentro de los
  dominios y constantes de contrato del perfil.
- **REQ-CONF-136.** Una configuración con `pcr_free` y cualquier parámetro
  de amplificación, o con `pcr` y un parámetro fuera de dominio, se rechaza
  antes de simular identificando parámetro y valor.
- **REQ-CONF-137.** La configuración efectiva publica la fracción esperada de
  duplicados para una biblioteca homogénea, calculada desde las parejas por
  molécula, como orientación y no como objetivo.
- **REQ-CONF-138.** Con `pcr` y loci de repetición declarados, las tasas de
  resbalón por ciclo se derivan por longitud de motivo de las fracciones con y
  sin PCR de la calibración y del número de ciclos, y pueden sustituirse
  directamente. Sin loci declarados, el tartamudeo heredable no aplica y se
  declara así.
- **REQ-DAT-81.** Con `pcr`, las `N` parejas de una biblioteca se organizan
  en familias de amplificación consecutivas por ordinal de plantilla. Cada
  familia sortea una molécula original, un número de copias tras los ciclos
  explícitos sesgado por tamaño, y un tamaño de familia de una Poisson
  condicionada a uno y pinzada en el máximo declarado. La molécula original
  se acepta con probabilidad igual a su factor de composición, con un número
  máximo de intentos; agotado, se acepta el último y se cuenta.
- **REQ-DAT-82.** El tartamudeo heredable se decide por familia y linaje:
  cada pareja tiene un linaje de un paso por ciclo, hereda de un progenitor
  de su familia los resbalones anteriores a su ciclo de divergencia y sortea
  los posteriores. Su paso neto se aplica a la molécula antes de reemitir
  R1/R2 y antes del tartamudeo de agrupación y síntesis.
- **REQ-DAT-83.** Ningún parámetro de amplificación altera el estado
  genotípico, la verdad publicada ni el VCF de verdad. Dos ejecuciones
  idénticas salvo en la preparación o en sus parámetros publican verdad, VCF y
  manifiestos de verdad byte a byte iguales.
- **REQ-DAT-84.** Con `pcr_free`, todos los artefactos publicados son byte a
  byte iguales a los de una ejecución equivalente anterior a la existencia de
  la preparación amplificada. Con `pcr` y todos sus parámetros neutros
  —parejas por molécula cero, curva neutra, tasas de resbalón cero—, los FASTQ
  finales son byte a byte iguales a los de `pcr_free` y solo difieren los
  campos declarativos del manifiesto.
- **REQ-DAT-85.** Toda decisión de amplificación y linaje se deriva de un
  ámbito propio del catálogo técnico, indexado por evento, biblioteca,
  ordinal de familia y, cuando corresponda, índice de copia, ciclo, intento y
  locus; nunca por posición ni por orden físico.
- **REQ-MAN-94.** El manifiesto de evento registra los recuentos agregados de
  REQ-MAN-45 y la cota de operaciones y memoria derivada de la configuración
  (sección 6), como valores calculados sin ejecutar y no como medidas.
- **REQ-VAL-58.** La validación de cierre debe contar sobre los FASTQ
  finales, con errores desactivados, quimera y tartamudeo a cero, el número
  de parejas con secuencias R1/R2 idénticas y comprobar que coincide con los
  duplicados realizados del manifiesto; y debe rederivar desde la semilla, sin
  reutilizar el productor, las familias y los pasos netos de linaje.
- **ACE-119.** Criterio de aceptación que resume 9.1 y 9.2: familias,
  duplicados antes de errores, composición al elegir la molécula, herencia
  del tartamudeo por linaje, verdad intacta, `pcr_free` idéntico byte a byte,
  cota publicada.

### 9.3. Decisión técnica nueva

Una fila **DTC-091** (provisional) que fije: el algoritmo
`pcr-amplification-v1` con los pasos de 5.3; la DP de los primeros ciclos y
su sesgo por tamaño; la Poisson truncada y pinzada como categórico de una
palabra, con la aclaración expresa de que DTC-067 no prohíbe construirla así;
el rechazo con tope y aceptación del último intento; la genealogía de 5.7 y
sus fórmulas; las constantes de contrato y su juez contra los esquemas; la
versión del catálogo técnico y el orden literal de componentes de 5.8; y la
fórmula de derivación de `p_ℓ`.

### 9.4. Registro de ampliaciones

DTF-056 se **divide**:

- Se retira de DTF-056 lo que este diseño cubre: preparación amplificada,
  duplicados de PCR, sesgo de composición configurable, sobrerrepresentación
  estocástica y tartamudeo heredable. La fila lo anota con fecha y remite a
  este documento y a los REQ nuevos.
- Permanece en DTF-056: perfil de sesgo GC **calibrado** con datos reales,
  duplicados ópticos y de agrupación, y combinación del sesgo con excepciones
  regionales.
- Se añade **DTF-059**: errores de polimerasa heredables por familia;
  quimeras formadas durante la PCR y heredadas; sesgo de longitud de la
  amplificación; genealogía exacta del proceso de ramificación. Dependencias:
  calibración de tasas, modelo de herencia común con el tartamudeo,
  procedencia agregada.
- DTF-055 no cambia: el marcado de duplicados sigue esperando a BAM/CRAM.

### 9.5. Hoja de ruta, plan de F4 y trazabilidad

- **Hoja de ruta §F4.3.** «Simular quimeras y, por elección expresa, una
  preparación amplificada con duplicados, sobrerrepresentación, composición y
  tartamudeo heredable; rechazar adaptadores, sesgo GC calibrado y anomalías
  *paired-end*, aplazados en DTF-054–056.» Y una nota fechada: el contrato de
  F4 se reabrió el día que se apruebe, con constancia en el plan de F4.
- **Plan de F4, sección 16 nueva**, con el título y la forma de la sección
  11 del plan de F5: «Cambio de contrato posterior al cierre — preparación
  amplificada (fecha)». Contenido: qué no se reabre (los hitos F4.0–F4.5 y sus
  evidencias), qué cambia y por decisión de quién, dónde vive cada literal
  (`LibraryPreparation`, los `const` de los tres esquemas, `ui/gps_app.py`),
  las pruebas nuevas con sus mutantes, los fixtures que se reabren y la
  repetición de la evidencia de cierre de F4 con `pcr_free`, que debe dar los
  mismos bytes.
- **Plan de F6.** F6.4 gana un paso, después de F6.4.9: «Tartamudeo heredable
  con PCR», que depende de F6.4.3 y de la parte C del expediente. La tabla de
  decisiones DP-F6-04 anota que el catálogo técnico nuevo incluye las etapas
  `stutter` y `amplification`.
- **Trazabilidad.** Filas de REQ-CONF-88/89/94/95, REQ-TEC-07/08,
  REQ-DAT-46/47/52/53, REQ-MAN-45/48/60, ACE-49/51/68 pasan de «F4.3/DTF-056»
  a «F4 §16»; filas nuevas para 9.2; REQ-DAT-48 y REQ-CONF-90 no cambian.
- **Calibración del tartamudeo §2.** «La columna con PCR aplica a la
  preparación `pcr` y se convierte en tasa por ciclo con la fórmula de este
  diseño».
- **Expediente, ficha D-12.** La opción (b) «reabrir DTF-056» es la que el
  usuario ha decidido para la preparación amplificada; la opción (a) sigue
  siendo la forma correcta de nombrar el tartamudeo sin PCR, que existe en
  ambos perfiles. La ficha se anota, no se reescribe.

---

## 10. Pruebas, y contra qué variante defectuosa debe fallar cada una

Cada prueba que exija mutar código sigue la disciplina de `CLAUDE.md` §4:
copia de seguridad fuera de git, SHA antes y después de mutar, restauración
comprobada byte a byte.

| # | Qué comprueba | Debe fallar contra | Cuidado |
|---|---|---|---|
| 1 | **`pcr_free` es hoy, en escenarios sin loci de repetición.** `sha256sum` de FASTQ, verdad, VCF y manifiestos con el binario nuevo y `pcr_free` frente a los del binario anterior. **Un escenario con loci declarados queda fuera de esta prueba a propósito**: ahí `pcr_free` sí debe cambiar, y lo cubre la parte C | Un lector que materialice un predeterminado nuevo en el efectivo; una versión de esquema publicada con `pcr_free` | Los hashes de referencia se miden **antes** de tocar el código y se guardan fuera de git; medirlos después es circular (precedente de `a1f8104`) |
| 2 | **Reducción al caso sin PCR.** `pcr` con `m = 0`, curva neutra, tasas cero: FASTQ byte a byte iguales a `pcr_free`; manifiesto distinto solo en los campos declarativos | Una implementación que consuma RNG con parámetros neutros o que desplace ordinales de plantilla | Comprobar además con un RNG instrumentado que el consumo es cero en `amplification` |
| 3 | **Duplicados por bytes.** `sequencing_errors = disabled`, quimera 0, sin loci: el número de parejas cuyas secuencias R1 y R2 coinciden con otra, contado por un oráculo que solo lee el FASTQ, es igual a `duplicate_pairs` del manifiesto; y con `m = 0,5` es mayor que cero | Un manifiesto que publique el recuento esperado en vez del realizado; una implementación que sortee `K` pero emita una sola pareja | Exigir `> 0` evita que pase por el motivo equivocado. Con errores activados la prueba no discrimina, porque los duplicados dejan de ser idénticos |
| 4 | **Tartamudeo clonal.** Fixture saturado: un locus dentro de la molécula, `m = 8` (familias grandes), divergencia forzada a `d = c` para todos (`q` degenerado en prueba), `p_ℓ = 1` en el ciclo 1 y `0` en los demás (sustitución directa por ciclo, solo en prueba). Todas las parejas de cada familia muestran el **mismo** número de unidades, distinto del verdadero | Un tartamudeo sorteado **por pareja** o por lectura; una herencia que copie los resbalones **posteriores** a `d` en vez de los anteriores | Sin forzar `d = c` y sin `p = 1`, algunas discordancias son legítimas y la aserción no discrimina. Docstring obligatoria diciéndolo. Una segunda variante con `d = 0` forzado exige que las parejas sean **independientes** |
| 5 | **La verdad no cambia.** Misma semilla, `pcr_free` y `pcr` con `m = 0,5` y curva dura: `truth.jsonl`, VCF y manifiestos de verdad con el mismo SHA; FASTQ distintos | Una implementación que ajuste el tracto en la copia haplotípica en vez de en la molécula | Exigir que los FASTQ **difieran** |
| 6 | **Composición al elegir la molécula.** Referencia con una mitad al 30 % GC y otra al 70 %, curva que da `g = 1` a la primera y `g = 0,2` a la segunda, `c = 10`, `m = 0,01`: la fracción de familias en la segunda mitad cae hasta el orden de `G = ((1+0,18)/(1+0,9))^10 ≈ 0,009`, muy por debajo de 0,5, con `A = 16` | Una implementación que aplique la composición solo a `λ` (argumento de 5.4): con `m = 0,01` daría ≈ 0,5 | Fijar `m` pequeño; con `m` grande la variante defectuosa también produce sesgo y no se distingue |
| 7 | **Tope de intentos.** Curva `g ≡ 0,01` salvo un óptimo estrecho y `A = 2`: `gc_attempts_exhausted_families` es alto y las familias aceptadas al agotar aparecen contadas; con `A = 16`, menor | Un rechazo sin tope (la prueba no termina o excede la cota de la prueba 10) | Ejecutar con límite de tiempo externo |
| 8 | **Primeros ciclos.** `c_early = 3`, `e = 0,5`: la DP produce la distribución exacta de `Z` en `1..8`, comparada con una enumeración a mano de los `2^3` caminos; la sesgada por tamaño tiene media `E[Z²]/E[Z]` | Una DP que olvide el término «no se copia» o que no sesgue por tamaño | Enumerar a mano, no reusar la DP |
| 9 | **Poisson truncada.** Para `λ = 0,3` y `k_max = 4`, los pesos coinciden con los calculados a mano y `p_4 = P(K ≥ 4 \| K ≥ 1)`; `λ = 0` no consume | Pesos normalizados sobre `1..k_max` sin pinzar la cola; un `λ = 0` que consuma | Contar consumo con RNG instrumentado |
| 10 | **Cota.** Contador de unidades en pruebas; fixture con `A`, `k_max`, `c_early`, `T` en máximos y curva casi siempre rechazada: `unidades ≤ ops_pcr(N)` con control positivo de conteo | Un mutante que suprima el término `A × (P + 3·L_max)` de la fórmula | Sin curva de rechazo, `A` no se ejercita y el mutante sobrevive |
| 11 | **Juez de constantes.** `C_MAX`, `E_MAX`, `K_MAX`, `A_MAX`, `Q_MAX` iguales a los `maximum` de los tres esquemas y a los de la UI | Un `maximum` cambiado en un esquema | Nueve comparaciones con control positivo, como en el plan F5 §11.4 |
| 12 | **Reproducibilidad y ámbitos.** Misma semilla, mismo SHA; añadir un locus en otro cromosoma no cambia los resbalones de los loci previos; añadir una biblioteca no cambia las familias de las demás | Un ámbito indexado por posición o por orden de biblioteca | — |
| 13 | **Validador independiente.** El oráculo rederiva familias, `Z`, `K` y pasos netos desde la semilla y la configuración efectiva sin llamar al productor, y coincide con los recuentos del manifiesto; `./scripts/check-independent-validators.sh` sobre un directorio nuevo | Un oráculo que lea el histograma publicado en vez de rederivarlo | La independencia algebraica de la Poisson truncada es la de sus pesos a mano (prueba 9), no la del oráculo |
| 14 | **Rechazos de configuración.** `pcr_free` con tabla `pcr`; `pcr` sin `cycles`; `cycles = 31`; `efficiency = 0`; `mean_pairs_per_molecule = -1`; curva con eficiencia 0 o nudos decrecientes; clases con pesos que no suman 1 fuera de tolerancia | Un esquema con `additionalProperties: true` o un dominio abierto | Cada caso cita parámetro y valor en el diagnóstico |

---

## 11. Decisiones abiertas para el usuario

| # | Decisión | Opciones | Recomendación |
|---|---|---|---|
| 1 | **Cómo se expresa la cantidad de material de partida** | (a) `mean_pairs_per_molecule` `m`, con la fracción esperada de duplicados publicada como orientación; (b) una fracción objetivo de duplicados que el sistema invierta numéricamente; (c) número absoluto de moléculas `C` | **(a).** Es la magnitud física, la fórmula es cerrada y no exige una inversión numérica cuya reproducibilidad habría que fijar. La UI muestra la fracción esperada |
| 2 | **Aceptación por composición** | (a) rechazo con tope `A` y aceptación del último intento, contada; (b) muestreo ponderado por ventanas con un recorrido previo del genoma | **(a)** para la primera versión: GC del fragmento, sin consumidor nuevo ni memoria por genoma; la distorsión del tope se publica. (b) si el coste `A×` resulta prohibitivo en genomas grandes; entonces F9 |
| 3 | **Genealogía dentro de la familia** | (a) la de 5.7: progenitor uniforme y divergencia geométrica; (b) estrella: todas las parejas divergen del linaje de la pareja 0; (c) coalescente exacto del proceso de ramificación por simulación | **(a).** (b) no representa resbalones compartidos por una rama sin la pareja 0; (c) es la alternativa A en pequeño y no tiene forma cerrada conocida |
| 4 | **Heterogeneidad de eficiencia entre moléculas** | (a) incluir las clases de eficiencia opcionales desde el principio, sin valores; (b) dejarlas a DTF-059 | **(a).** Cuesta un categórico y `Q` tablas; sin valores predeterminados no afirma nada |
| 5 | **Quimeras con PCR** | (a) por pareja como hoy; (b) una vez por familia, heredada | **(a)** ahora, (b) en DTF-059 con los errores de polimerasa, que comparten el mecanismo de herencia |
| 6 | **Tartamudeo heredable: cómo se dan las tasas** | (a) derivarlas de las columnas con y sin PCR de la calibración y de `c`, con sustitución directa; (b) solo tasas por ciclo declaradas por el usuario | **(a).** Es lo único que hace usable la columna con PCR, y la derivación es una fórmula publicada; la etiqueta «humano, Illumina» viaja con los valores |
| 7 | **Perfil** | (a) recurso hermano `gps-hiseqx-pcr-pe150-v1`; (b) un solo recurso con dos preparaciones | **(a)**, por el hash congelado y los `const` del vigente |
| 8 | **Versiones de esquema** | (a) la PCR abre las versiones nuevas que F6 planea y F6 las completa; (b) la PCR espera a que F6.1.1 las abra | **(a)** si la primera entrega de la sección 7 va antes que F6; en ambos casos, misma versión mientras no esté publicada externamente |
| 9 | **La cota de recursos** | (a) se publica en el manifiesto y se prueba, sin rechazo; (b) se impone con una admisión F4 que reutilice las constantes de DTC-089 | **(a).** F4 nunca ha rechazado por presupuesto y (b) reabre contratos de F5; F9.4 es el sitio de la imposición |
| 10 | **Cobertura de moléculas distintas en el manifiesto** | (a) añadir `unique_molecule_x`; (b) solo la fracción de duplicados | **(a).** Es lo que un analista compara tras eliminar duplicados y cuesta una división |
| 11 | **Entregas** | (a) dos entregas según la sección 7; (b) una sola tras F6.4 | **(a).** La primera es verificable sola y no depende de las dieciséis decisiones abiertas de F6 |

---

## 12. Lo que se deja fuera a propósito

Ya enumerado en 3.3 y en 9.4: errores de polimerasa heredables, quimeras
heredadas, sesgo de longitud, genealogía exacta, duplicados ópticos y de
agrupación, marcado de duplicados, perfil GC calibrado. Y además:

- **Calibrar nada.** Ningún valor de este documento, salvo las fracciones de
  tartamudeo de Gymrek ya documentadas, se propone como predeterminado; los
  parámetros que fijan el fenómeno (`c`, `e`, `m`) son obligatorios para que
  el usuario los declare. Es la regla de la ficha D-14 y la de DTF-056.
- **Reabrir el orden del recorrido más allá de lo necesario.** Contaminación,
  demultiplexado y escritura no cambian.
- **Una vía de rendimiento.** El coste `A ×` de la aceptación y el factor dos
  de las dos pasadas se declaran y no se optimizan; F9.

---

## 13. Referencias

Todas **verificadas el 19 de septiembre de 2026** contra la fuente publicada:
existencia, revista, año y **que dicen lo que aquí se les atribuye**. Antes de
esa fecha siete de las ocho estaban citadas de memoria.

- **Gymrek M. (2016).** *PCR-free library preparation greatly reduces stutter
  noise at short tandem repeats.* bioRxiv 043448. **CC BY 4.0.** En el
  repositorio, `design/calibracion-tartamudeo-0.1.md`. Fuente de las fracciones
  de tartamudeo con y sin PCR. **Verificada.**

- **Aird D., Ross M. G., Chen W.-S. y col. (2011).** *Analyzing and minimizing
  PCR amplification bias in Illumina sequencing libraries.* Genome Biology
  12:R18. <https://doi.org/10.1186/gb-2011-12-2-r18>
  **Verificada** para la PCR como causa principal del sesgo de composición: el
  artículo la identifica como «a principal source of bias».
  **No verificada** la forma unimodal de la curva que este documento le
  atribuye; si esa forma se usa para fijar algún valor, hay que leer el
  artículo completo antes.

- **Benjamini Y., Speed T. P. (2012).** *Summarizing and correcting the GC
  content bias in high-throughput sequencing.* Nucleic Acids Research
  40(10):e72. <https://doi.org/10.1093/nar/gks001>
  **Verificada, y sostiene una decisión del diseño:** es el contenido de bases
  **del fragmento completo**, no el de la lectura, lo que más influye en el
  recuento. Respalda aplicar la composición al **elegir la molécula** y no al
  repartir copias.

- **Kozarewa I. y col. (2009).** *Amplification-free Illumina
  sequencing-library preparation facilitates improved mapping and assembly of
  (G+C)-biased genomes.* Nature Methods 6:291–295. PMID 19287394.
  **Verificada:** sin amplificación hay menos duplicados y menos sesgo. Además
  describe que el protocolo sin PCR **usa la amplificación en agrupación sobre
  la celda** en lugar de la PCR, que es el mecanismo por el que el tartamudeo
  sigue existiendo en bibliotecas sin amplificar.

- **Best K., Oakes T., Heather J. M., Shawe-Taylor J., Chain B. (2015).**
  *Computational analysis of stochastic heterogeneity in PCR amplification
  efficiency revealed by single molecule barcoding.* Scientific Reports
  5:14629. PMID 26459131.
  **Verificada, con un matiz que refuerza el diseño:** la heterogeneidad de
  eficiencia es **independiente de la secuencia del cebador y de las
  condiciones globales**, y se atribuye tanto a diferencias heredadas entre
  moléculas plantilla como a la estocasticidad propia de la reacción. Es decir,
  modelarla como una clase aleatoria por molécula no es arbitrario: la
  secuencia no la explica.

- **Kebschull J. M., Zador A. M. (2015).** *Sources of PCR-induced distortions
  in high-throughput sequencing data sets.* Nucleic Acids Research 43(21):e143.
  <https://doi.org/10.1093/nar/gkv717>
  **Verificada, con dos precisiones.** Primera: es un trabajo
  **experimental** —marcado de moléculas individuales— con modelado, no una
  simulación del árbol como decía la versión anterior de esta ficha. Segunda,
  y **sostiene lo que aquí se aplaza**: encuentra que la estocasticidad es la
  fuerza principal que distorsiona la representación, mientras que los errores
  de polimerasa, aun siendo frecuentes en los ciclos tardíos, **apenas afectan
  a la distribución**, y los cambios de plantilla —las quimeras de PCR— son
  **raros y limitados a números de copia bajos**. Aplazar ambos a DTF-059 queda
  así respaldado por medición y no solo por conveniencia.

- **Sun F. (1995).** *The polymerase chain reaction and branching processes.*
  Journal of Computational Biology 2(1):63–86. PMID 7497121.
  **Weiss G., von Haeseler A. (1995).** *Modeling the polymerase chain
  reaction.* Journal of Computational Biology 2(1):49–61. PMID 7497120.
  **Ambas verificadas**, con volumen, número y páginas. Son el origen del
  proceso de ramificación como modelo de la PCR y de los estimadores de tasa de
  mutación y eficiencia.

- **Número de ciclos y tasas típicas de duplicados de los protocolos Illumina
  con PCR: sigue sin fuente.** Por eso este diseño **no fija ningún
  predeterminado** para esos parámetros y los declara obligatorios en la
  configuración. No es un descuido: es lo que impide publicar como valor de
  referencia un número que nadie ha medido aquí.
