# Plan de implementación — repeticiones en tándem y tartamudeo, con y sin PCR 0.1

## 1. Estado, autoridad y alcance

- **Fecha:** 19 de septiembre de 2026.
- **Estado:** **ejecutable bajo supuestos declarados.** Este plan es una
  excepción a la regla que prohíbe a un plan decidir diseño. Las veintiocho
  fichas de decisión abiertas que reparten los dos documentos de origen se han
  resuelto **adoptando la opción que cada documento recomienda**, y cada
  adopción está marcada como **supuesto** en la tabla de la sección 2 y en el
  paso que la usa. El usuario puede leer esa tabla, confirmar o cambiar
  cualquier fila, y el plan sigue siendo válido: cambiar un supuesto cambia los
  pasos que lo citan y ningún otro.
- **Por qué la excepción.** El plan de la fase de variación compleja
  (`plan-implementacion-f6-0.1.md`) se paró ante cada decisión abierta y dejó
  21 de sus 46 pasos bloqueados. Un plan hueco no sirve para trabajar. Aquí se
  adopta lo recomendado, se dice en voz alta y se deja revocable.
- **Fuentes normativas y de diseño:**
  - `design/expediente-decisiones-repeticiones-tandem-0.1.md` — corte mínimo
    para tener microsatélites que varíen entre individuos, dieciocho fichas de
    decisión y el borrador del requisito del tartamudeo.
  - `design/diseno-bibliotecas-pcr-0.1.md` — muestreo por familias de
    amplificación, duplicados, sesgo por composición, tartamudeo clonal, cota
    de recursos, cambio de contrato y once decisiones.
  - `design/calibracion-tartamudeo-0.1.md` — valores medidos con y sin PCR,
    **declarados humanos**, con sus fuentes y su situación de licencia.
  - `design/especificaciones-0.1.md`, `design/decisiones-tecnicas-0.1.md`,
    `design/hoja-de-ruta-0.1.md` (§7, reglas de avance),
    `design/plan-implementacion-f6-0.1.md` §8 (hito de repeticiones),
    `design/plan-implementacion-f4-0.1.md` (fase de realismo técnico, cerrada).
  - Requisitos existentes que el trabajo toca: REQ-CONF-55/57/59/60/63,
    REQ-CONF-88/89/94/95, REQ-IND-31/32/33/34, REQ-REF-19/23, REQ-TEC-07/08,
    REQ-DAT-17/46/47/52/53/56, REQ-MAN-29/45/48/60, REQ-VER-20/22,
    ACE-31/49/51/68/103/110, DTC-049, DTC-056, DTC-061, DTC-064, DTC-067,
    DTC-072, DTC-081, DTC-086, DTC-087, DTC-088, DTC-089, DTF-055, DTF-056.
- **Modelo de forma:** el plan de la fase de población fundadora
  (`plan-implementacion-f5-0.1.md`), y en particular su sección 11, que es el
  precedente de cómo este proyecto documenta un **cambio de contrato posterior
  al cierre de una fase**. La sección 9 de este plan copia esa forma.
- **Qué entrega.** Dos cosas que el usuario pidió, en este orden:
  1. **Microsatélites que varían entre individuos**: loci de repetición
     declarados en coordenadas de la referencia y alelos —números de unidades—
     distintos entre individuos y entre las dos copias de un mismo individuo,
     publicados en la verdad y en el VCF.
  2. **Tartamudeo visible en las lecturas**: el cambio del número de unidades
     observado por deslizamiento de la polimerasa, primero en la biblioteca sin
     amplificar y después en la amplificada por PCR, con sus duplicados y su
     herencia por linaje.
- **Qué no entrega, y dónde queda:** sección 12.
- **Prohibiciones que este plan hereda y no relaja.** No autoriza *commit* ni
  publicación externa. No revierte trabajo ajeno. No se etiqueta una puerta con
  un comando que no se haya ejecutado. Las copias de seguridad de cualquier
  mutación van **fuera de git**; nunca `git checkout --`.

---

## 2. Tabla de supuestos

**Cómo se lee.** Cada fila es una decisión que los documentos de origen dejaron
abierta y que este plan ha resuelto adoptando lo que ese documento recomienda.
La columna «Opción adoptada» es lo que el plan da por decidido. La columna
«Pasos que la usan» dice exactamente qué hay que rehacer si el usuario cambia
la fila.

Las diecisiete primeras proceden del expediente de repeticiones (fichas D-01 a
D-18, menos D-10, que es una confirmación y no una decisión). Las once últimas,
de la sección 11 del diseño de bibliotecas con PCR.

| # | Qué se decide, en llano | Opción adoptada | Ficha | Pasos que la usan |
|---|---|---|---|---|
| **S-01** | Dónde se declaran los loci y los alelos, y qué versión de fichero de configuración nace | Esquema general **v4** que extiende la vía directa: tabla `[[tandem_repeat_loci]]` y clase `kind = "tandem_repeat"` dentro de `direct_events`. Sin escenario enlazado nuevo | D-01 (a) | F6.1.1.1, F6.4.8 |
| **S-02** | Cómo se guarda el alelo dentro del estado genotípico | **Modificación pequeña tipada**: variante nueva `SmallModification::TandemRepeat` junto a `Snp`, `Insertion` y `Deletion`, que materializa como la inserción o la deleción equivalente. **No** se crea ninguna clase de segmento estructural | D-02 (a) | F6.4.3, F6.4.3.2 |
| **S-03** | De dónde sale el tracto que hay en la referencia | Las dos vías, **empezando por la suministrada**: se admite un FASTA que ya contiene el tracto y se comprueba el intervalo declarado; después se siembra el tracto en la referencia simulada | D-03 (c), con (a) primero | F6.4.2 |
| **S-04** | Cómo se asignan los alelos a los individuos | **Explícito** en el fichero general ahora; frecuencias alélicas por locus y población, después. La forma del TOML de ahora no debe impedir añadir las frecuencias luego | D-04 (a) ahora, (c) después | F6.1.1.1, F6.4.3 |
| **S-05** | Qué pasa con un alelo que se separa 50 pb o más de la referencia | **Rechazo antes de simular**, identificando locus, individuo y copia. El corte queda íntegramente en territorio de *indels* pequeños | D-05 (a) | F6.4.3.1, F6.4.8 |
| **S-06** | Pureza e interrupciones: qué se declara y qué se calcula | **Interrupciones declaradas** como pares (desplazamiento dentro del tracto, base); **pureza calculada y publicada** como valor realizado. Una configuración que declare una pureza objetivo se rechaza | D-06 (a) | F6.4.1, F6.4.8 |
| **S-07** | Por qué extremo se ganan o se pierden unidades | **Extremo derecho del tracto**, unidades completas; el resto del tracto queda byte a byte igual al de la referencia, y la diferencia es siempre un *indel* limpio | D-07 (a) | F6.4.1.1, F6.4.3 |
| **S-08** | Forma de la verdad publicada | **Caso *indel* con metadatos de repetición**: mismos tipos de registro, más `locus_id`, `motif` y `units` en el alelo, en una versión nueva del esquema de verdad (v3). `truth-v1` y `truth-v2` no cambian ni un byte | D-08 (a) | F6.4.3.2 |
| **S-09** | Campos del VCF y qué es «tamaño manejable» | **Campos de repetición en tándem del estándar 4.5**, declarados en la cabecera, con alelos explícitos en `REF`/`ALT`. El umbral de «manejable» es **el mismo** que separa *indel* de variante estructural (50 pb por omisión). **Sin fichero complementario** | D-09 (a) | F6.5.1.1 |
| **S-11** | Ámbitos de aleatoriedad e identificadores | Identificador del locus **declarado por el usuario** y validado por unicidad, sin clase de identificador automático nueva. Para el tartamudeo, catálogo técnico nuevo **`f4-technical-v4`** con etapa `stutter` y componentes `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`, indexada por `locus_id` y no por posición | D-11 | F4.8.1 |
| **S-12** | Si cabe un artefacto de PCR en una fase que declara no usar PCR | **Nombrarlo por su mecanismo** para el caso sin amplificar: deslizamiento de la polimerasa en la amplificación en puente y en la síntesis, que existe en el perfil vigente. Para la biblioteca amplificada, el usuario ya decidió reabrir DTF-056 (ver contradicción C-3) | D-12 (a), ampliada | F4.8.1, sección 9 |
| **S-13** | En qué punto exacto del recorrido se aplica el tartamudeo | **Después de la quimera y antes de calidad y errores**, con un `stutter::transform_paired_end` de la misma forma que los dos que ya existen | D-13 | F4.8.4 |
| **S-14** | El modelo de probabilidad | **`stutter-v1`**: ensayo de Bernoulli con `p = min(p_max, max(0, a + b·(u − u₀)))` por longitud de motivo, y sorteo categórico sobre un conjunto declarado de pasos enteros con pesos. Sin reintentos: los pasos inviables se excluyen **antes** de sortear | D-14 | F4.8.2, F4.8.3 |
| **S-15** | Qué pasa con la longitud del fragmento | **La molécula cambia de longitud** en `k·(longitud del motivo)`, y con ella el inserto observado. No se compensa ni se descarta la molécula. Las dos lecturas conservan la longitud del perfil | D-15 (a) | F4.8.5 |
| **S-16** | Compatibilidad con lo ya publicado | Congelación byte a byte **acotada a los escenarios que no declaran loci de repetición** (ver contradicción C-1, que corrige la ficha) | D-16, corregida | F4.8.6, F4.9.10 |
| **S-17** | Si la interfaz Marimo es condición de cierre | **Sí para el tartamudeo y sí, en versión mínima, para los loci** (ver contradicción C-6, que corrige la ficha: sin control de loci, el control de tartamudeo no sirve para nada) | D-17, corregida | F6.4.9.1, F4.8.9, F4.9.12 |
| **S-18** | Si esto es un hito propio o adelanta el hito de repeticiones | **Adelanta pasos existentes** de la fase de variación compleja, con la corrección de dependencia escrita, y añade hitos nuevos a la fase de realismo técnico con su cambio de contrato documentado | D-18 (b) | Sección 4, sección 9 |
| **S-19** | Cómo se expresa la cantidad de material de partida de la PCR | **Parejas por molécula** (`mean_pairs_per_molecule`), con la fracción esperada de duplicados publicada como orientación | PCR 1 (a) | F4.9.3, F4.9.5 |
| **S-20** | Cómo actúa el sesgo por composición | **Rechazo con tope de intentos** y aceptación del último, contada. Se aplica al **elegir la molécula**, no al repartir copias | PCR 2 (a) | F4.9.6 |
| **S-21** | Genealogía dentro de la familia | **Progenitor uniforme entre las parejas anteriores y ciclo de divergencia geométrico** con `q = 1/(1+e)` | PCR 3 (a) | F4.10.2 |
| **S-22** | Heterogeneidad de eficiencia entre moléculas | **Incluir las clases de eficiencia desde el principio**, opcionales y **sin valores predeterminados** | PCR 4 (a) | F4.9.3, F4.9.7 |
| **S-23** | Quimeras con PCR | **Por pareja, como hoy.** La quimera heredada por la familia queda para más adelante | PCR 5 (a) | F4.9.7 |
| **S-24** | Cómo se dan las tasas de resbalón por ciclo | **Derivadas** de las columnas con y sin PCR de la calibración y del número de ciclos, con sustitución directa admitida | PCR 6 (a) | F4.10.3 |
| **S-25** | Perfil tecnológico | **Recurso hermano** `resources/f4/gps-hiseqx-pcr-pe150-v1.json`, dejando intacto el recurso sin PCR y su SHA-256 | PCR 7 (a) | F4.9.1 |
| **S-26** | Versiones de esquema | Las versiones nuevas (general v4, efectivo v6, manifiesto de evento v3, verdad v3) las **abre este trabajo** y la fase de variación compleja las completa dentro de la misma versión mientras no estén publicadas fuera | PCR 8 (a) | F6.1.1.1, F4.9.3 |
| **S-27** | La cota de recursos | **Se publica** en el manifiesto y se prueba, **sin rechazo**. La imposición pertenece a la fase de presupuesto de recursos | PCR 9 (a) | F4.9.9 |
| **S-28** | Cobertura de moléculas distintas en el manifiesto | **Se añade** `unique_molecule_x`, que es lo que un analista compara tras eliminar duplicados | PCR 10 (a) | F4.9.8 |
| **S-29** | Reparto en entregas | **Dos entregas** en el diseño; este plan las convierte en **cuatro**, por la razón aritmética de la sección 3 | PCR 11 (a), refinada | Sección 3 |

**Por qué falta S-10.** La numeración de los supuestos sigue a la de las fichas
del expediente: S-01 corresponde a D-01, S-02 a D-02, y así hasta S-18. **No
hay S-10** porque la ficha D-10 no es una decisión, sino la confirmación que
sigue. Los supuestos S-19 a S-29 corresponden, en orden, a las once decisiones
numeradas 1 a 11 de la sección 11 del diseño de bibliotecas con PCR.

**Confirmación registrada, que no es un supuesto.** La ficha D-10 no decide
nada: DTC-064, la decisión que fija cómo se normaliza la proyección al VCF, ya
obliga a alinear los *indels* a la izquierda, de modo que un alelo de
repetición quedará anclado en el extremo izquierdo del tracto, y ya obliga a un
único registro por cromosoma y posición con varios `ALT`. Un diploide con 5 y 7
unidades produce **un** registro con dos alternativas y genotipo faseado. No
hay que tocar nada. Se anota para que no sorprenda en F6.5.1.1.

**Aviso de recuento.** El encargo cuenta veintiséis decisiones abiertas: quince
del expediente más once del diseño con PCR. El expediente anuncia en su parte B
que «tres» de sus dieciocho fichas son confirmaciones, pero **solo marca una**,
D-10. Por prudencia, este plan trata como supuesto todo lo que no está marcado
como confirmación, y de ahí salen veintiocho filas en vez de veintiséis. Si el
usuario identifica las otras dos confirmaciones, esas filas se bajan de la
tabla sin más consecuencia.

### 2.1. Decisiones **sin recomendación**, que este plan deja pendientes

No se inventan. El plan trabaja alrededor de ellas y dice cómo.

| # | Qué falta | Por qué no hay recomendación | Cómo se planifica alrededor |
|---|---|---|---|
| **P-1** | **Cómo crece la tasa de tartamudeo con la longitud del tracto.** Son los coeficientes `b` y `u₀` del modelo | La calibración lo dice sin rodeos: Fungtammasan y colaboradores (2015) lo describen, pero el artículo está tras muro de pago y **no se ha podido extraer la curva ajustada**; «no debe inventarse uno» | Los parámetros existen y son configurables, **con `b = 0` como valor de partida**, con lo que la probabilidad se reduce a la fracción medida por longitud de motivo, que sí está calibrada. El manifiesto publica `tract_length_dependence: not_calibrated`. F4.8.2 prueba que la función **sí** depende de `u` cuando `b > 0`; F4.8.3 prueba que con `b = 0` reproduce la tabla |
| **P-2** | **Qué versión de `bcftools` y `tabix` vale para la puerta externa del VCF.** El entorno tiene 1.21 y la decisión vigente exige 1.24 | Es la decisión pendiente DP-F6-14 del plan de variación compleja, que enumera opciones y **no recomienda ninguna** | F6.5.1.1 ejecuta la comprobación con la versión instalada y **registra la evidencia como no contractual, diciendo la versión**. La puerta externa **no se etiqueta como superada**. El paso no se bloquea, porque el oráculo semántico propio del proyecto sí discrimina |
| **P-3** | **Los pesos por omisión del conjunto de pasos** —cuánto pesa perder una unidad frente a perder dos o ganar una | La calibración establece el hecho cualitativo («sesgo fuerte hacia perder repeticiones», «la distribución del tamaño es parecida con y sin PCR») pero **no publica los pesos**, y el expediente escribe `{−2, −1, +1}` como ejemplo | El conjunto de pasos y sus pesos son **obligatorios** en la configuración cuando el tartamudeo está activo: sin valor predeterminado, igual que el diseño con PCR hace con los ciclos y la eficiencia. F4.8.3 rechaza una configuración de tartamudeo sin pesos declarados |

### 2.2. Contradicciones encontradas al adoptar las recomendaciones

Se declaran aquí porque son hallazgos, no obstáculos. Las tres primeras cambian
lo que hay que construir.

**C-1. «Tartamudeo desactivado por omisión» contra «un escenario con loci debe
tartamudear sin PCR».** La ficha D-16 y los requisitos propuestos REQ-DAT-80 y
REQ-CONF-132 piden el tartamudeo desactivado por omisión y la identidad byte a
byte con la ejecución anterior. El diseño con PCR corrige esa congelación el
mismo día en que se escribe: vale para los escenarios que **no declaran loci de
repetición**, que son todos los que existen hoy; uno que sí los declare debe
tartamudear también sin amplificar, porque el fenómeno existe sin PCR y está
medido en la columna sin PCR de la calibración. Y hay una razón aritmética
además de una empírica: si el caso sin PCR se deja en cero, el caso con PCR
queda solo con el incremento despejado y **subestima** el tartamudeo entre un
14 % y un 28 % según el motivo.

*Resolución adoptada:* manda la corrección. La congelación se acota a los
escenarios sin loci declarados. Con loci declarados, el tartamudeo de base está
**activo** con las tasas de la calibración, y lo que sigue desactivado por
omisión es el **tartamudeo heredable de PCR**, que solo existe con preparación
amplificada. REQ-DAT-80 y REQ-CONF-132 deben reescribirse con esa acotación
antes de aprobarse; la redacción está en la sección 10.

**C-2. Dos formas distintas del mismo modelo.** La calibración §5 propone la
forma de tres parámetros de HipSTR —probabilidad de añadir, probabilidad de
quitar y un parámetro geométrico del tamaño—. La ficha D-14 propone
`stutter-v1`: un ensayo de Bernoulli y un sorteo categórico sobre pasos
declarados. Además, D-14 afirma que «no existe en el repositorio ninguna
calibración de tartamudeo» y por eso recomienda valores predeterminados
**cero**; esa premisa **ya no es cierta**: `design/calibracion-tartamudeo-0.1.md`
existe, con fuentes verificadas y licencia CC BY 4.0.

*Resolución adoptada:* se conserva el mecanismo de D-14, porque un categórico
sobre pasos declarados **puede** expresar la forma geométrica de HipSTR sin
atarse a ella, y porque encaja con las primitivas normativas que el proyecto ya
tiene. Se abandona la parte de D-14 que recomendaba valores cero: las tasas por
longitud de motivo salen de la calibración, con su etiqueta «humano, Illumina»
viajando con ellas. Los pesos del tamaño del paso siguen sin valor (P-3).

**C-3. D-12 recomienda no reabrir la ampliación aplazada; el usuario la ha
reabierto.** La ficha D-12 recomienda nombrar el tartamudeo por su mecanismo
para evitar tocar DTF-056, la ampliación que aplaza la preparación amplificada.
El diseño con PCR parte de que el usuario ya decidió simular bibliotecas
amplificadas, y por tanto **divide** DTF-056. No es un choque real, sino un
solapamiento de fechas: la recomendación (a) de D-12 sigue siendo la forma
correcta de nombrar el tartamudeo **sin** PCR, que existe en los dos perfiles, y
la reapertura cubre la preparación amplificada. El plan usa las dos cosas y lo
dice en la sección 9.

**C-4. La ficha D-18 nombra menos coordenadas de las que el corte toca.** Dice
«adelantar F6.4.1 y F6.4.3». El corte necesita además el paso de la referencia
con el tracto (F6.4.2), el del registro en configuración efectiva y manifiesto
(F6.4.8), el de la interfaz (F6.4.9), un subconjunto del paso de contratos
(F6.1.1) y un subconjunto del paso de proyección VCF (F6.5.1). Se resuelve
nombrándolos todos en la sección 4.

**C-5. Un límite conocido que la opción adoptada en S-04 alcanza pronto.**
DTC-089 fija un techo duro de 65536 bytes al fichero general. Declarar cada
alelo a mano cuesta de 70 a 90 bytes; con 100 individuos diploides y 5 loci son
entre 70 y 90 kB, es decir, **se rompe el techo antes de llegar a un estudio de
tamaño normal**. No es una contradicción con otra decisión adoptada, pero sí un
límite que el plan declara y comprueba (F6.1.1.1), y la razón por la que la
forma del TOML debe admitir después las frecuencias alélicas sin cambiar de
versión de esquema otra vez.

**C-6. La asimetría de D-17 deja la interfaz inservible para lo que pide el
usuario.** D-17 recomienda control Marimo para el tartamudeo y lo difiere para
los loci, con el argumento de que el catálogo de loci se escribe una vez y los
parámetros del tartamudeo se tocan muchas. El argumento es bueno pero el
resultado no funciona: **el tartamudeo solo actúa sobre loci declarados**, de
modo que una interfaz con controles de tartamudeo y sin forma de declarar un
locus no permite producir ni una sola lectura tartamudeada. El plan corrige la
ficha: la interfaz incluye una tabla mínima de loci y de alelos (F6.4.9.1), y
la tabla completa —parámetros aleatorios y regionales— sigue difiriéndose al
cierre del hito de repeticiones de la fase de variación compleja.

---

## 3. El orden es forzoso, y por qué. El reparto en entregas

### 3.1. La aritmética manda el orden

Las dos columnas de la calibración —«con PCR» y «sin PCR»— son **totales, no
componentes**. Los dos conjuntos que la fuente comparó pasaron por las mismas
células, la misma amplificación en puente sobre la celda y la misma síntesis de
secuenciación; lo único que los distingue es el paso de PCR de preparación. Por
tanto la columna sin PCR **ya contiene** el mosaicismo somático de la muestra,
el deslizamiento del puente y el de la síntesis, y la columna con PCR contiene
todo eso **más** la PCR de biblioteca.

El diseño con PCR no suma nada a la columna con PCR: **despeja** el incremento
propio de la amplificación,

```
r_pcr = 1 − (1 − r_con) / (1 − r_sin)
```

y lo compone después con la parte de base, de modo que
`1 − (1−r_sin)(1−r_pcr) = r_con`, exactamente la fracción medida.

**Consecuencia que fija el orden de este plan:** una preparación amplificada
entregada **sin** el tartamudeo de base se queda solo en el incremento
despejado y subestima el tartamudeo en un 14 % en homopolímeros, 28 % en
dinucleótidos, 26 % en trinucleótidos, 23 % en tetranucleótidos y 20 % en
pentanucleótidos. No es una preferencia de orden: es que la aritmética no
cierra. La preparación amplificada va **después** del tartamudeo de base.

### 3.2. Qué no depende de que existan los microsatélites, y qué sí

| Pieza | ¿Necesita loci de repetición? |
|---|---|
| Perfil amplificado, variante `Pcr`, configuración, esquemas y rechazos | **No** |
| Muestreo por familias, copias de los primeros ciclos, tamaño de familia | **No** |
| Duplicados y su oráculo por bytes | **No** |
| Aceptación por composición de bases | **No** |
| Recuentos del manifiesto, fracción esperada y realizada, cobertura de moléculas distintas, cota de recursos | **No** |
| Tartamudeo de base (agrupación y síntesis) | **Sí** |
| Genealogía de la familia y tartamudeo heredable | **Sí**, y además necesita el tartamudeo de base |

### 3.3. El reparto: cuatro entregas, no dos

El diseño con PCR recomienda **dos** entregas: primero la amplificación con sus
duplicados y su composición, que no depende de nada de la fase de variación
compleja; después el tartamudeo clonal. Este plan **respeta esa separación pero
reordena y subdivide**, y quedan cuatro:

| Entrega | Contenido | Coordenadas | Depende de |
|---|---|---|---|
| **1** | Microsatélites que varían entre individuos: loci, alelos por copia, verdad, VCF, manifiesto, interfaz | F6.4.0–F6.4.9.1 (adelanto) | Nada |
| **2** | Tartamudeo de base, sin amplificar | F4.8.1–F4.8.9 | Entrega 1 |
| **3** | Preparación amplificada: familias, duplicados, composición, cota | F4.9.1–F4.9.12 | Nada técnicamente; **ordenada aquí por la aritmética** |
| **4** | Tartamudeo clonal heredable por linaje | F4.10.1–F4.10.6 | Entregas 1, 2 y 3 |

**Por qué cuatro y no dos.**

1. **Porque lo que el usuario pidió cabe en las dos primeras.** Microsatélites
   que varían entre individuos y tartamudeo visible en las lecturas se entregan
   al terminar la entrega 2, sin haber tocado el contrato de duplicados ni de
   sesgo por composición. Una entrega única que empezara por la amplificación
   retrasaría el resultado pedido detrás de un cambio de contrato grande.
2. **Porque la aritmética prohíbe el orden contrario.** Poner la amplificación
   antes del tartamudeo de base dejaría existir, aunque fuera unos días, una
   biblioteca amplificada que subestima el tartamudeo entre un 14 % y un 28 %.
3. **Porque la primera versión del esquema la abre quien la necesita.** El
   supuesto S-26 dice que las versiones nuevas las abre el trabajo que llegue
   primero. Con este orden, las abre la entrega 1 —que es la que tiene el
   contrato de entrada más grande— y las entregas 3 y 4 las completan dentro de
   la misma versión mientras no estén publicadas fuera.
4. **Porque la entrega 3 sigue siendo verificable sola**, que es la virtud que
   el diseño reclamaba para ella: sus doce pasos no citan ni un locus.

**Qué se pierde y cómo se compensa.** Si el usuario prefiriera el orden del
diseño —amplificación primero—, el plan sigue siendo válido intercambiando las
entregas 3 y 2, **con una condición inexcusable**: mientras falte el tartamudeo
de base, la configuración debe **rechazar** la combinación de preparación
amplificada con loci de repetición declarados, con un diagnóstico que nombre el
paso pendiente. Ese rechazo está planificado de todos modos en F4.9.10, porque
también hace falta en el orden recomendado mientras la entrega 4 no exista.

---

## 4. Cómo se numera, y puertas comunes

### 4.1. Numeración

Este trabajo cruza dos fases, una cerrada y otra sin empezar. Siguiendo la
recomendación de la ficha D-18 (supuesto S-18), **no** se crea una numeración
nueva:

- Los pasos de la entrega 1 **adelantan pasos que ya existen** en el plan de la
  fase de variación compleja, con su coordenada intacta: F6.4.1, F6.4.2,
  F6.4.3, F6.4.8 y F6.4.9. Cuando el corte necesita solo una parte de un paso
  mayor de aquel plan, se usa una subdivisión excepcional: **F6.1.1.1** es el
  subconjunto de repeticiones del paso de contratos externos, y **F6.5.1.1** el
  subconjunto de repeticiones del paso de proyección VCF. **F6.4.0** es nuevo y
  es la línea base de evidencias.
- Los pasos de las entregas 2, 3 y 4 son **hitos nuevos de la fase de realismo
  técnico**, numerados a continuación de los que existen, que llegan hasta
  F4.7: **F4.8** (tartamudeo de base), **F4.9** (preparación amplificada) y
  **F4.10** (tartamudeo clonal). Su cambio de contrato se documenta en la
  sección 9, con la forma de la sección 11 del plan de la fase de población
  fundadora.

**Correcciones de dependencia que este plan declara** y que hay que aplicar al
plan de la fase de variación compleja (sección 10):

- **F6.4.3 ya no depende del primer paso del hito de variantes estructurales.**
  La cadena de segmentos (`src/segments.rs`) ya sostiene modificaciones tipadas
  por copia haplotípica —`Snp`, `Insertion`, `Deletion`— y
  `src/materialization.rs` ya las materializa. Un alelo de repetición es,
  mecánicamente, una inserción o una deleción de unidades completas con
  etiqueta.
- **F6.4.2 ya no depende de derivar una referencia estructural.** El generador
  de referencias simuladas (`src/reference_generation.rs`) y el escritor FASTA
  canónico (`src/reference_fasta.rs`) existen desde la primera fase.
- **F6.5.1.1 no arrastra el fichero complementario del VCF.** El proyecto emite
  `VCFv4.5`, cuyo estándar ya define campos propios de repetición en tándem; si
  el formato lo expresa, no hay pérdida que complementar.

### 4.2. Puertas

Cada paso ejecuta sus pruebas específicas. Antes de cada pausa de hito se
ejecuta la batería completa, con los comandos literales:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
./scripts/check-independent-validators.sh <directorio-nuevo>
MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
uv run --frozen python -m unittest discover -s tests_ui
```

- **Base de partida:** `cargo test --workspace --all-targets --locked` da hoy
  **1126 pasan, 0 fallan, 12 ignoradas**. Toda entrega declara su recuento
  nuevo y la diferencia.
- `./scripts/check-independent-validators.sh` **falla si el directorio ya
  existe** y debe dar **5 PASS**. Cada hito usa un directorio propio:
  `evidence/repeats-f6-4`, `evidence/stutter-f4-8`, `evidence/pcr-f4-9`,
  `evidence/clonal-f4-10`.
- `uv run --frozen python -m unittest discover -s tests_ui` da hoy **96
  ejecutadas con 2 fallos preexistentes y ajenos** (interfaz del VCF de la fase
  2 y del evento técnico de la fase 4). No se persiguen. Tampoco las doce
  pruebas ignoradas históricas.
- **`pytest` no es puerta.** No está instalado ni declarado en
  `pyproject.toml`. Ningún paso de este plan lo cita.
- **Regla de cierre de la interfaz** (hoja de ruta §7): toda opción
  configurable que este plan implemente exige control en Marimo y equivalencia
  interfaz → TOML → configuración efectiva antes de cerrar su hito.

### 4.3. Disciplina de pruebas

Cada paso dice **contra qué variante defectuosa concreta** debe fallar la
prueba que pide. Es la respuesta al defecto recurrente número uno del proyecto:
pruebas que no pueden fallar por el defecto que dicen cubrir, que apareció diez
veces durante la fase de población fundadora y que hace pocos días dejó al
descubierto dos guardas del verificador frente a una batería de 1124 pruebas.

Reglas, sin excepción:

1. Toda regresión se verifica **fallando** contra la variante defectuosa.
2. Se comprueba que **la mutación mutó**: SHA del fichero antes y después.
3. La copia de seguridad va **fuera de git**. Nunca `git checkout --`.
4. Si una aserción **no puede ejercitarse**, se dice en su docstring, sin
   adornos, en vez de fingir que discrimina.
5. Cuando un fixture ordinario no ejercita la guarda —porque otra comprobación
   caza antes la discrepancia—, hace falta un fixture **saturado**, que haga
   indistinguible todo lo demás. Los pasos que lo necesitan lo dicen.

---

## 5. Entrega 1 — Microsatélites que varían entre individuos (adelanto de F6.4)

Resultado: loci de repetición en tándem declarados en coordenadas de una
referencia, alelos por copia haplotípica que difieren entre individuos, y su
publicación en la verdad, en el VCF y en el manifiesto, con control en la
interfaz. Sin tartamudeo todavía.

- ☐ **F6.4.0.** Línea base de evidencias, medida **antes** de tocar una sola
  línea de código y guardada **fuera de git**.
  - Depende de: ninguno.
  - Procede de: CLAUDE.md §4 y la prueba 1 de la sección 10 del diseño con PCR,
    que recuerda que medir los hashes después de tocar el código es circular.
  - Cambia: nada del repositorio. Crea
    `/tmp/gps-baseline-<fecha>/` con: `sha256sum` de los FASTQ finales, la
    verdad, el VCF de verdad, la configuración efectiva y los dos manifiestos
    de al menos dos ejecuciones de referencia; la salida completa de
    `cargo test --workspace --all-targets --locked`; y el SHA del `HEAD`
    (`11d26e5`).
  - Se comprueba con: el fichero de hashes existe, tiene una línea por
    artefacto y **ninguna** ruta bajo el árbol de git. **Debe fallar contra:**
    una línea base capturada después del primer cambio —se detecta porque el
    `HEAD` registrado no coincide con el de la captura.

- ☐ **F6.4.1.** Modelo de dominio del locus de repetición, puro: motivo
  A/C/G/T no vacío y sin `N`, longitud de unidad coincidente con la del motivo,
  número de copias de referencia, intervalo en la referencia, patrón de
  interrupciones **declarado** como pares (desplazamiento dentro del tracto,
  base sustituta) y **pureza calculada** como valor realizado.
  - Depende de: ninguno. Es código puro sobre secuencias en memoria; no hay
    FASTA, ni TOML, ni aleatoriedad.
  - Procede de: REQ-CONF-57/59/60, REQ-IND-32/33/34, REQ-REF-23, ACE-31/110;
    supuesto **S-06**; paso F6.4.1 del plan de la fase de variación compleja,
    que este paso **desbloquea** al cerrar su decisión pendiente sobre la
    realización de interrupciones.
  - Cambia: módulo nuevo `src/tandem_repeat.rs`.
  - Se comprueba con: `(ACG)×5` produce quince bases; longitud de unidad
    declarada 4 para el motivo `ACG` rechazada; motivo vacío, motivo `ACN` y
    motivo en minúsculas mal normalizado rechazados, cada uno identificando el
    locus; una configuración que declare una **pureza objetivo** rechazada con
    el parámetro nombrado; sin interrupciones la pureza realizada vale
    exactamente 1; con dos bases interrumpidas de veinte la pureza realizada
    vale exactamente 0,9; la secuencia realizada solo contiene A/C/G/T. **Debe
    fallar contra:** un validador que compare la longitud declarada con la del
    motivo en mayúsculas pero acepte `acn`; una pureza realizada que cuente las
    interrupciones sobre las unidades en vez de sobre las bases del tracto; un
    lector que acepte silenciosamente la pureza como entrada y la ignore.

- ☐ **F6.4.1.1.** Realización del tracto con un número de unidades distinto del
  de la referencia: las unidades se añaden o se quitan **por el extremo
  derecho** del tracto, en unidades completas, y el resto del tracto queda byte
  a byte igual al de la referencia, interrupciones incluidas.
  - Depende de: F6.4.1.
  - Procede de: supuesto **S-07**; ficha D-07, que advierte de que este es
    «exactamente el tipo de detalle que dos implementaciones resuelven distinto
    sin que ninguna prueba lo note».
  - Cambia: `src/tandem_repeat.rs`.
  - Se comprueba con: un tracto con una interrupción en la tercera unidad y
    alelos de 5, 7 y 3 unidades; las bases anteriores al extremo derecho son
    idénticas en los tres, comparadas byte a byte, y la diferencia frente a la
    referencia es exactamente una inserción o una deleción contigua de
    `k·(longitud del motivo)` bases. **Debe fallar contra:** una realización
    que **regenere** el tracto entero desde el motivo y las interrupciones —que
    es la opción (c) de la ficha—, porque entonces la interrupción se recoloca
    y la diferencia deja de ser un *indel* limpio. *Cuidado:* con un tracto
    **sin** interrupciones la regeneración da el mismo resultado y la prueba no
    discrimina; el fixture **debe** llevar interrupción, y su docstring debe
    decirlo.

- ☐ **F6.1.1.1.** Contrato de entrada del corte: esquema general **v4** que
  extiende la vía directa con la tabla `[[tandem_repeat_loci]]` —`locus_id`,
  `contig`, `start_1based`, `end_1based`, `motif`, interrupciones— y la clase
  `kind = "tandem_repeat"` dentro de `direct_events`, con `locus_id`, `phase` y
  `units`; y configuración efectiva **v6**.
  - Depende de: F6.4.1.
  - Procede de: REQ-CONF-60, REQ-MAN-29, DTC-089; supuestos **S-01**, **S-04**,
    **S-26**; corrección C-5.
  - Cambia: `schemas/config-general-v4.schema.json` y
    `schemas/effective-config-v6.schema.json` nuevos;
    `src/configuration.rs`, `src/adapters.rs`, `src/results.rs`. **Ningún
    esquema v1 a v5 cambia, ni `truth-v1`, ni `truth-v2`.**
  - Se comprueba con: pareja mínima válida admitida; clave desconocida dentro
    de la tabla de loci rechazada (`additionalProperties: false`);
    identificador de locus repetido rechazado por unicidad; un alelo que cita un
    `locus_id` inexistente rechazado nombrando el evento; validación externa con
    `python3 -m jsonschema` en Draft 2020-12, que es ajena al productor;
    comparación **byte a byte** de los esquemas heredados contra los de
    `11d26e5`; un documento legítimo de 65536 bytes admitido y el de 65537
    rechazado; y —esta es la que importa— **un fichero general v2 sin loci
    sigue publicando `effective-config-v4`, no v6**. **Debe fallar contra:** un
    lector que promueva todo documento a v6; un esquema heredado con una
    propiedad añadida; un esquema nuevo con `additionalProperties: true`.
  - *Límite declarado (C-5):* con la declaración explícita de alelos, 100
    individuos diploides y 5 loci ocupan entre 70 y 90 kB y **no caben** en el
    techo de 65536 bytes de DTC-089. El paso lo comprueba con un documento que
    lo rebasa y lo registra como límite conocido, no como defecto. La vía por
    frecuencias alélicas (S-04, segunda mitad) es la salida, y la forma del
    TOML de este paso no debe impedirla.

- ☐ **F6.4.2.** El tracto en la referencia, por las dos vías: (a) **referencia
  suministrada** que ya lo contiene, con una comprobación de admisión que
  verifique que el intervalo declarado es el motivo repetido un número entero
  de veces con sus interrupciones declaradas, y (b) **referencia simulada** con
  el tracto sembrado.
  - Depende de: F6.4.1.1, F6.1.1.1.
  - Procede de: REQ-REF-19/23, REQ-CONF-56; supuesto **S-03**. Corrección de
    dependencia: **ya no depende** de derivar una referencia estructural.
  - Cambia: `src/reference.rs` (admisión) y `src/reference_generation.rs`
    (siembra); el escritor FASTA canónico `src/reference_fasta.rs` no cambia.
  - Se comprueba con: un FASTA de juguete escrito a mano con `(AC)×10` en las
    coordenadas declaradas es admitido; el mismo con una base cambiada dentro
    del tracto es rechazado nombrando el locus y el desplazamiento; un
    intervalo declarado que no es múltiplo entero del motivo es rechazado; en
    la vía simulada, el FASTA contiene la secuencia realizada en la coordenada
    declarada, localizada con `grep -ob` externo, y misma semilla da mismo
    FASTA. **Debe fallar contra:** una escritura que coloque el locus una base
    más allá —el `grep -ob` externo la caza—; una admisión que compare solo la
    longitud del intervalo y no su contenido.

- ☐ **F6.4.3.** El alelo por copia haplotípica: variante nueva
  `SmallModification::TandemRepeat` con `locus_id`, número de unidades y ancla,
  que se materializa como la inserción o la deleción equivalente.
  - Depende de: F6.4.2.
  - Procede de: REQ-IND-32/33, REQ-CONF-58/59, ACE-103, DTC-061; supuesto
    **S-02**. Corrección de dependencia: **ya no depende** del primer paso del
    hito de variantes estructurales.
  - Cambia: `src/segments.rs` (la variante y su validación de extensión y
    orden), `src/materialization.rs`, `src/small_event_composition.rs`.
  - Se comprueba con: las dos copias de un diploide con 5 y 7 unidades
    materializan longitudes distintas y secuencias distintas; tres individuos
    con 4, 5 y 9 unidades producen tres secuencias distintas en el mismo locus;
    dos modificaciones sobre el mismo tracto se rechazan por solapamiento, como
    ya hace `SegmentError::OverlappingModifications`; la reconstrucción por un
    lector independiente devuelve el número de unidades declarado. **Debe
    fallar contra:** una materialización que use el número de copias **de la
    referencia** para ambas copias haplotípicas —el caso 5/7 la caza—; una
    variante que guarde el *indel* pero pierda el número de unidades, que es lo
    que REQ-IND-33 prohíbe.

- ☐ **F6.4.3.1.** Umbral de 50 pb: rechazo antes de simular de todo alelo cuya
  diferencia con la referencia **alcance** el umbral configurable que separa
  *indels* pequeños de variantes estructurales, identificando locus, individuo
  y copia.
  - Depende de: F6.4.3.
  - Procede de: REQ-CONF-55, REQ-IND-31; supuesto **S-05**.
  - Cambia: `src/configuration.rs`, `src/tandem_repeat.rs`.
  - Se comprueba con: motivo de 4 pb y alelo con 12 unidades de más (48 pb)
    admitido; con 13 (52 pb) rechazado nombrando locus, individuo y copia; con
    el umbral configurado a 20 pb, el de 48 pb pasa a rechazarse, lo que
    demuestra que se lee el valor configurado y no una constante. **Debe fallar
    contra:** un umbral aplicado solo al generar y no al validar la entrada —el
    caso del umbral configurado a 20 lo caza—; una comparación estricta
    `> umbral` donde REQ-IND-31 exige `>=`.

- ☐ **F6.4.3.2.** Verdad publicada, versión **v3**: el caso sigue siendo de
  tipo *indel* y el alelo gana `locus_id`, `motif` y `units`. `truth-v1` y
  `truth-v2` no cambian ni un byte, y una ejecución sin loci de repetición
  publica exactamente lo que publica hoy.
  - Depende de: F6.4.3.1.
  - Procede de: REQ-VER-22 —que prohíbe expresamente dejar el detalle restante
    únicamente en el manifiesto—, REQ-IND-33; supuesto **S-08**. Cierra la
    decisión pendiente DP-F6-03 en su parte de repeticiones, y solo en esa.
  - Cambia: `schemas/truth-v3.schema.json` nuevo, `src/truth_schema.rs`,
    `src/truth_projection.rs`, `src/truth_validation.rs`.
  - Se comprueba con: el registro de un diploide 5/7 lleva los dos alelos con
    su fase, su `locus_id`, su motivo y su número de unidades; el validador de
    verdad existente sigue aceptando un documento v2 sin cambios; `sha256sum`
    de la verdad de una ejecución **sin** loci, comparado con el de F6.4.0,
    idéntico; los dos esquemas heredados comparados byte a byte. **Debe fallar
    contra:** una proyección que publique el *indel* y omita el número de
    unidades; un cambio que promueva a v3 las ejecuciones sin repeticiones —el
    `sha256sum` contra la línea base lo caza.

- ☐ **F6.5.1.1.** Proyección al VCF de los alelos de repetición: campos de
  repetición en tándem **del estándar `VCFv4.5`** declarados en la cabecera,
  con alelos explícitos en `REF` y `ALT`, y el umbral de «tamaño manejable»
  fijado en el **mismo** umbral que separa *indel* de variante estructural.
  - Depende de: F6.4.3.2.
  - Procede de: REQ-VER-20/22, ACE-31, DTC-064 y DTC-066; supuesto **S-09**;
    confirmación D-10. Cierra la decisión pendiente DP-F6-06 y, al hacerlo,
    **mantiene fuera del corte el fichero complementario del VCF, su catálogo
    de códigos de pérdida y su esquema**.
  - Cambia: `src/vcf_contract.rs` (cabecera), `src/vcf_projection.rs`,
    `src/vcf_representation.rs`, `src/vcf_semantic_validation.rs`.
  - **Primero, verificar el estándar.** Antes de escribir una sola línea de
    cabecera, comprobar contra la especificación publicada de `VCFv4.5` los
    nombres, números y tipos exactos de los campos de repetición en tándem que
    la versión 4.4 introdujo y la 4.5 conserva, y **dejar la cita en el
    código**. El proyecto ya pagó una vez el precio de dar por buena una cifra
    no comprobada; no se repite.
  - Se comprueba con: un diploide 5/7 produce **un** registro con dos
    alternativas y genotipo faseado, alineado a la izquierda del tracto, tal
    como DTC-064 obliga; la cabecera declara cada campo de repetición usado y
    ninguno más; el oráculo semántico propio reconstruye el número de unidades
    desde `REF`/`ALT` y lo compara con la verdad; y la comprobación externa con
    las herramientas instaladas. **Debe fallar contra:** una cabecera que use
    un campo no declarado; una proyección que emita dos registros para el mismo
    sitio en vez de uno multialélico; un campo de repetición cuyo valor
    contradiga el número de unidades de la verdad.
  - **Pendiente P-2, dicho sin adornos:** la puerta externa del VCF exige
    `bcftools` y `tabix` en la versión 1.24, y el entorno tiene 1.21. Este paso
    ejecuta la comprobación con 1.21, **registra la versión en la evidencia** y
    **no etiqueta la puerta externa como superada**. El oráculo semántico
    propio sí discrimina y es el que sostiene el paso.

- ☐ **F6.4.8.** Registro en configuración efectiva y manifiesto: loci y
  parámetros configurados, procedencia concreta o aleatoria, motivo, unidad,
  número de copias o longitud, y los **valores realizados** de pureza e
  interrupciones.
  - Depende de: F6.5.1.1.
  - Procede de: REQ-MAN-29, REQ-VAL-48; supuestos **S-01**, **S-05**, **S-06**.
  - Cambia: `src/event_manifest.rs`, `src/results.rs`, `schemas/manifest-general-v6.schema.json`.
  - Se comprueba con: validación externa de los dos documentos con
    `python3 -m jsonschema`; la pureza publicada coincide con la calculada a
    mano sobre el tracto del fixture; el umbral efectivo publicado coincide con
    el configurado, no con el predeterminado, cuando la configuración lo
    sustituye. **Debe fallar contra:** un manifiesto que registre la pureza
    declarada como realizada —imposible aquí, porque no se declara, y por eso
    la aserción que lo vigila es la del **umbral sustituido**, que sí puede
    fallar; se dice así en su docstring.

- ☐ **F6.4.9.1.** Interfaz Marimo, versión mínima: tabla de loci de repetición
  y tabla de alelos por individuo y copia, con equivalencia
  interfaz → TOML → configuración efectiva.
  - Depende de: F6.4.8.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; supuesto **S-17**
    corregido por la contradicción **C-6**: sin este control, el control de
    tartamudeo de la entrega 2 no permite producir ni una lectura tartamudeada.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`, `ui/gps_app.py`;
    prueba nueva en `tests_ui/`.
  - Se comprueba con: ida y vuelta de la interfaz al TOML y del TOML a la
    configuración efectiva, al estilo de
    `tests_ui/test_f5_6_1_scenario_ui.py`; ausencia de controles para lo que el
    corte no admite —pureza objetivo, generación aleatoria de loci, modelo de
    mutación escalonado—; `MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo
    check --strict ui/gps_app.py` sin avisos. **Debe fallar contra:** un
    control de pureza objetivo presente; una ida y vuelta que pierda las
    interrupciones declaradas.

> **Revisión del hito de la entrega 1.** Comprobar el modelo del locus, la
> realización por el extremo derecho, el alelo por copia, el umbral, la verdad
> v3, el VCF y la interfaz, con sus evidencias. Ejecutar la batería completa de
> puertas de la sección 4.2 con
> `./scripts/check-independent-validators.sh evidence/repeats-f6-4`. Confirmar
> que `sha256sum` de los artefactos de una ejecución **sin** loci coincide con
> la línea base de F6.4.0. No continuar hasta confirmar el resultado.

---

## 6. Entrega 2 — Tartamudeo de base, sin amplificar (F4.8)

Resultado: una molécula que contiene un tracto declarado puede mostrar un
número de unidades distinto del que el individuo tiene, por deslizamiento de la
polimerasa durante la amplificación en puente sobre la celda de flujo y durante
la síntesis de la secuenciación. **La verdad del individuo no cambia.**

**Cómo se nombra, y por qué no contradice el contrato de la fase.** Una
biblioteca «sin PCR» evita la amplificación en la preparación, pero no evita
toda polimerasa: la formación de agrupaciones en la celda es amplificación en
puente, y la secuenciación por síntesis es copia con polimerasa. El
deslizamiento ocurre en ambas, y por eso los microsatélites tartamudean también
en datos de genoma completo sin PCR. El requisito nuevo debe decir
explícitamente que **este artefacto no es una preparación amplificada**, que no
introduce duplicados y que el perfil sigue registrándose como sin PCR
(supuesto **S-12**).

- ☐ **F4.8.1.** Decisión técnica nueva y catálogo de ámbitos aleatorios
  **`f4-technical-v4`**, que extiende el vigente `f4-technical-v3` con la etapa
  `stutter` y los componentes
  `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`, con
  `decision` en `bernoulli` y `step`, indexada por **identificador de locus** y
  nunca por la posición del tracto dentro de la molécula.
  - Depende de: cierre de la entrega 1.
  - Procede de: REQ-DAT-15, DTC-037, DTC-081, DTC-086, DTC-088; supuestos
    **S-11** y **S-12**. Fila nueva propuesta en
    `design/decisiones-tecnicas-0.1.md`: **DTC-090** (provisional).
  - Cambia: `src/randomness.rs`.
  - **Aclaración que la decisión debe llevar por escrito:** DTC-067 obliga a
    que los **recuentos** aleatorios salgan de `poisson-count-v1`. Aquí no hay
    ningún recuento de sucesos sobre oportunidades: hay un ensayo binario y una
    elección entre alternativas. Los consumidores que corresponden ya existen
    —`bernoulli-v1`, `uniform-index-v1` (DTC-088) y `categorical-f64-v1`
    (DTC-072)— y **DTC-067 no los prohíbe aquí**. Decirlo evita que alguien lea
    la regla al revés.
  - Se comprueba con: vectores congelados de derivación de subsemilla para una
    terna evento/biblioteca/locus; el catálogo vigente `f4-technical-v3` sigue
    produciendo exactamente las mismas subsemillas que antes, comparadas con la
    línea base. **Debe fallar contra:** una extensión que renumere o reordene
    las etapas existentes —los vectores del catálogo vigente lo cazan.

- ☐ **F4.8.2.** Modelo `stutter-v1`, puro: probabilidad
  `p = min(p_max, max(0, a + b·(u − u₀)))` por longitud de motivo, donde `u` es
  el número de unidades completas **contenidas en la molécula**; y sorteo
  categórico sobre un conjunto declarado de pasos en unidades enteras con sus
  pesos. **Sin reintentos**: un paso solo entra en el sorteo si deja al menos
  una unidad completa en el tracto y deja la molécula no más corta que la
  lectura. Si ningún paso es viable, no hay tartamudeo y **no se consume
  aleatoriedad**.
  - Depende de: F4.8.1.
  - Procede de: REQ-DAT-77 propuesto, DTC-072, DTC-088; supuesto **S-14**;
    pendientes **P-1** y **P-3**.
  - Cambia: módulo nuevo `src/stutter.rs`.
  - Se comprueba con: **(a)** la función de probabilidad, llamada
    directamente, es estrictamente mayor para un tracto de 30 unidades que para
    uno de 6 con el mismo motivo cuando `b > 0`, y **igual** cuando `b = 0`,
    que es el valor de partida de P-1; **(b)** con la misma longitud **en pares
    de base** —un mononucleótido de 30 unidades y un tetranucleótido de 7 más
    dos bases— la probabilidad es mayor para el mononucleótido; **(c)** el
    filtrado previo excluye los pasos inviables y un generador instrumentado
    confirma que **no se consume** aleatoriedad cuando no queda ningún paso
    viable; **(d)** con peso 1 a la contracción y 0 a la expansión, todas las
    moléculas que tartamudean salen más cortas, y un segundo caso con pesos
    repartidos confirma que el categórico sí se ejercita. **Debe fallar
    contra:** un modelo con probabilidad constante que ignore la longitud del
    tracto (prueba a); un modelo que dependa solo de los pares de base del
    tracto y no de la longitud del motivo (prueba b); una implementación que
    sortee y descarte en vez de filtrar antes (prueba c, que cuenta el consumo,
    no el resultado); un signo invertido (prueba d). *Cuidado:* el caso con
    peso 1 no consume aleatoriedad en el categórico, de modo que por sí solo no
    discrimina un categórico roto; de ahí el segundo caso.
  - **P-3 y P-1, dichos en el paso:** el conjunto de pasos y sus pesos son
    **obligatorios** cuando el tartamudeo está activo, sin valor
    predeterminado, porque la calibración establece el sesgo hacia la
    contracción pero no publica los pesos. El coeficiente `b` de la dependencia
    con la longitud del tracto vale **0** de partida y el manifiesto lo publica
    como no calibrado.

- ☐ **F4.8.3.** Tasas por longitud de motivo tomadas de la calibración: la
  columna **sin PCR** de `design/calibracion-tartamudeo-0.1.md`, que con `b = 0`
  es directamente el valor de `a` por longitud de motivo. La procedencia
  —Gymrek (2016), CC BY 4.0, **humano, Illumina**— viaja con los valores en la
  configuración efectiva y en el manifiesto.
  - Depende de: F4.8.2.
  - Procede de: REQ-CONF-131/132 propuestos; contradicciones **C-1** y **C-2**;
    calibración §3.
  - Cambia: `src/stutter.rs`, `src/technical_profile.rs`,
    `schemas/config-general-v4.schema.json`,
    `schemas/effective-config-v6.schema.json`.
  - **Valores de partida** (fracción de lecturas afectadas, columna sin PCR):
    motivo de 1 pb, 0,027; 2 pb, 0,011; 3 pb, 0,0029; 4 pb, 0,0016; 5 pb,
    0,0015; 6 pb, 0,0010, **este último marcado como el menos fiable de la
    tabla**, porque el propio autor advierte que la diferencia de los
    hexanucleótidos procede probablemente del ruido de analizar pocos loci.
  - Se comprueba con: la configuración efectiva publica los seis valores, la
    cita de la fuente y la etiqueta de especie y tecnología; una configuración
    con tartamudeo activo y **sin** pesos de paso declarados se rechaza
    nombrando el parámetro; una tasa fuera de `[0,1]`, un paso expresado en
    bases en vez de en unidades y un modelo distinto del admitido se rechazan
    antes de simular. **Debe fallar contra:** unos valores predeterminados
    embebidos en el código sin la etiqueta de procedencia; una configuración
    sin pesos que se acepte con un reparto implícito.

- ☐ **F4.8.4.** `stutter::transform_paired_end`, con la misma forma que los dos
  transformadores que ya existen, aplicado **después de la quimera y antes de
  calidad y errores**. Una única decisión por molécula y tracto, aplicada a la
  molécula completa antes de volver a emitir R1 y R2.
  - Depende de: F4.8.3.
  - Procede de: REQ-DAT-74/75/76 propuestos; supuesto **S-13**.
  - Cambia: `src/stutter.rs`, `src/execution.rs` (el orden del recorrido),
    `src/wgs_sampling.rs` si el punto de inserción lo exige.
  - **Por qué ahí y no antes.** `src/chimera.rs` vuelve a materializar el
    prefijo del fragmento desde la copia haplotípica y lo compara **byte a
    byte** con el prefijo del fragmento que recibe; si no coinciden devuelve
    `ChimeraError::PairSourceMismatch`. Un fragmento ya modificado por el
    tartamudeo rompería esa comprobación cada vez que la molécula resultara
    elegida para formar quimera. Después de la quimera, la pareja todavía lleva
    su fragmento completo, de modo que la molécula sigue disponible para
    reescribirla y volver a emitir las dos lecturas.
  - Se comprueba con: **(a) coherencia del par** sobre un fixture
    **saturado**: copia haplotípica lo bastante corta para forzar el inserto
    mínimo del perfil, 250 pb, de modo que R1 y R2 —de 150 pb— compartan las 50
    bases centrales, y un tracto de 30 pb colocado **íntegramente dentro de ese
    solapamiento**; se cuentan las unidades en R1 y en el complemento inverso
    de R2 y se exige que sean iguales. **(b) quimera**: un corpus con tasa de
    quimera y tasa de tartamudeo altas a la vez no produce ni un
    `ChimeraError::PairSourceMismatch`. **Debe fallar contra:** un tartamudeo
    aplicado **por lectura** (prueba a); una implementación que modifique el
    fragmento **antes** de la quimera (prueba b). *Cuidado, y va en la
    docstring:* un fixture en el que solo R1 cubra el tracto **no discrimina**
    —la aserción se cumple sola—, y con la media de inserto de 350 pb la
    mayoría de las parejas están en ese caso; de ahí el fixture saturado. Y sin
    tasa de quimera alta, la prueba (b) no ejercita el camino que dice cubrir.

- ☐ **F4.8.5.** La molécula cambia de longitud en `k·(longitud del motivo)` y
  con ella el tamaño de inserto observado. No se compensa materializando
  secuencia flanqueante ni recortando. Las dos lecturas conservan la longitud
  del perfil.
  - Depende de: F4.8.4.
  - Procede de: REQ-DAT-78 propuesto; supuesto **S-15**.
  - Cambia: `src/stutter.rs`, `src/paired_end.rs`.
  - Se comprueba con: una molécula que pierde una unidad de un motivo de 4 pb
    produce un fragmento exactamente 4 pb más corto, observado **sobre la
    molécula**, no sobre las lecturas; las dos lecturas siguen midiendo 150 pb;
    una contracción que dejaría la molécula más corta que la lectura queda
    excluida del sorteo por el filtrado de F4.8.2 y no produce ningún error
    duro. **Debe fallar contra:** una implementación que recorte o rellene para
    conservar la longitud sorteada. *Cuidado:* observar las lecturas no
    discrimina, porque siguen midiendo lo mismo pase lo que pase; la aserción
    tiene que mirar la molécula.

- ☐ **F4.8.6.** Compatibilidad con lo ya publicado, **acotada**: una ejecución
  de un escenario que **no declara loci de repetición** publica exactamente los
  mismos bytes que antes de que el tartamudeo existiera. Un escenario que **sí**
  los declara tartamudea con las tasas de la calibración, y esa es la corrección
  C-1.
  - Depende de: F4.8.5.
  - Procede de: REQ-DAT-80 propuesto **con la acotación de C-1**, DTC-081;
    supuesto **S-16**.
  - Cambia: nada; es una prueba de sistema y su evidencia.
  - Se comprueba con: **(a)** `sha256sum` de los dos FASTQ, la verdad, el VCF,
    la configuración efectiva y los dos manifiestos de un escenario **sin**
    loci, comparados con la línea base de F6.4.0; **(b)** un generador
    instrumentado confirma **cero consumo** en la etapa `stutter` en ese
    escenario; **(c)** el mismo escenario **con** un locus declarado y su tasa
    de calibración produce FASTQ **distintos** y verdad, VCF y manifiestos de
    verdad **idénticos**; **(d)** dos ejecuciones idénticas salvo en la tasa de
    tartamudeo difieren solo en los FASTQ y en los recuentos técnicos
    agregados. **Debe fallar contra:** una implementación que consuma
    aleatoriedad con la tasa a cero o que desplace el ordinal de plantilla
    (pruebas a y b); una implementación que aplique el tartamudeo a la **copia
    haplotípica** en vez de a la molécula, es decir, que cambie el genotipo
    (prueba c); una lectura de la congelación como absoluta, que dejaría el
    caso con loci en cero tartamudeo (prueba c, que exige que los FASTQ
    **difieran**). *Cuidado:* si en la prueba (c) los FASTQ salieran iguales, la
    comparación de verdad pasaría por el motivo equivocado; por eso la
    diferencia de los FASTQ es una aserción, no una observación.

- ☐ **F4.8.7.** Configuración efectiva y manifiesto del tartamudeo: activación,
  parámetros efectivos por longitud de motivo y por locus, versión del
  algoritmo, versión del catálogo de ámbitos aleatorios y, como **recuentos
  agregados**, moléculas evaluadas, moléculas que tartamudearon y distribución
  de los pasos aplicados por locus. Nada por lectura ni por pareja.
  - Depende de: F4.8.6.
  - Procede de: REQ-MAN-93 propuesto, REQ-DAT-56 —que prohíbe conservar
    procedencia por lectura—; supuesto **S-14**; pendientes **P-1** y **P-3**.
  - Cambia: `src/event_manifest.rs`, `schemas/manifest-event-v3.schema.json`
    nuevo.
  - Se comprueba con: validación externa del manifiesto con
    `python3 -m jsonschema`; la suma de la distribución de pasos por locus es
    igual al número de moléculas que tartamudearon; el manifiesto marca los
    parámetros con `calibration: "gymrek-2016-human-illumina"` y la dependencia
    con la longitud del tracto como `not_calibrated`; ningún campo nombra una
    lectura o una pareja concreta (`grep` negativo con control positivo). **Debe
    fallar contra:** un manifiesto que publique el número **esperado** de
    moléculas tartamudeadas en vez del realizado; una procedencia por pareja
    filtrada.

- ☐ **F4.8.8.** Validador independiente del tartamudeo: reproduce las
  decisiones desde la semilla de la ejecución **sin reutilizar el código
  productor**, y comprueba sobre los FASTQ finales que las dos lecturas de cada
  pareja que cubren un mismo tracto muestran el mismo número de unidades.
  - Depende de: F4.8.7.
  - Procede de: REQ-VAL-57 propuesto, arquitectura §11.
  - Cambia: validador nuevo bajo `tests/`, más su entrada en
    `scripts/check-independent-validators.sh` si procede.
  - Se comprueba con: el validador coincide con los recuentos agregados del
    manifiesto sobre el fixture saturado de F4.8.4;
    `./scripts/check-independent-validators.sh evidence/stutter-f4-8` da **5
    PASS**. **Debe fallar contra:** un validador que lea los bytes publicados
    por el productor en vez de recomponerlos desde la semilla —se detecta
    mutando un recuento del manifiesto y exigiendo que el validador lo rechace.

- ☐ **F4.8.9.** Interfaz Marimo del tartamudeo: activación, parámetros por
  longitud de motivo, conjunto de pasos y pesos, y sustituciones por locus, con
  equivalencia interfaz → TOML → configuración efectiva.
  - Depende de: F4.8.8, F6.4.9.1.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; supuesto **S-17**
    corregido por **C-6**.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`; prueba nueva en
    `tests_ui/`.
  - Se comprueba con: ida y vuelta completa; ausencia de control para un modelo
    distinto del admitido y para un paso expresado en bases;
    `marimo check --strict` sin avisos; y **la comprobación que C-6 motiva**:
    partiendo solo de la interfaz, se puede construir un escenario completo con
    un locus y una tasa que produzca lecturas tartamudeadas. **Debe fallar
    contra:** una interfaz con parámetros de tartamudeo y sin forma de declarar
    un locus, que es exactamente lo que la ficha D-17 recomendaba y este plan
    corrige.

> **Revisión del hito de la entrega 2.** Comprobar el punto de inserción en el
> recorrido, la coherencia del par sobre el fixture saturado, la congelación
> acotada, la verdad intacta, el validador independiente y la interfaz.
> Ejecutar la batería completa con
> `./scripts/check-independent-validators.sh evidence/stutter-f4-8`. Aplicar el
> cambio de contrato de la sección 9 en su parte de tartamudeo y **repetir la
> evidencia de cierre de la fase de realismo técnico**. No continuar hasta
> confirmar el resultado.

---

## 7. Entrega 3 — Preparación amplificada por PCR (F4.9)

Resultado: un modo de preparación `pcr`, alternativo al único que hoy existe,
en el que el recorrido técnico deja de producir una pareja por fragmento y pasa
a producir **familias de amplificación**: grupos de parejas que descienden de
la misma molécula original. De esa estructura salen los duplicados, la
sobrerrepresentación de unas moléculas frente a otras y el sesgo por
composición de bases. **Ningún paso de este hito cita un locus de repetición.**

**Cómo se acota lo inacotable.** No se simula la población de moléculas, que
crece como `2^ciclos`. Se muestrea **hacia atrás** desde las parejas que se van
a publicar: como el número de parejas `N` de una biblioteca ya es determinista
desde la cobertura, se generan familias una tras otra hasta cubrirlas, y cada
familia cuesta un número de operaciones acotado por constantes de contrato. La
cota resultante es lineal en `N`.

- ☐ **F4.9.1.** Perfil hermano `resources/f4/gps-hiseqx-pcr-pe150-v1.json`,
  idéntico al vigente en lectura, inserto, calidad y error, con
  `"preparation": "pcr"` y un bloque de amplificación con sus dominios; y
  variante `Pcr` en `LibraryPreparation`.
  - Depende de: cierre de la entrega 2 (por el orden aritmético de §3.1; no hay
    dependencia técnica).
  - Procede de: REQ-CONF-82, REQ-CONF-89 reescrito; supuesto **S-25**.
  - Cambia: recurso nuevo; `src/technical_profile.rs`.
  - Se comprueba con: el recurso vigente y su SHA-256 **no cambian**, comparados
    byte a byte contra la línea base; el perfil nuevo se carga y declara su
    preparación. **Debe fallar contra:** una edición del recurso vigente para
    meter las dos preparaciones en uno solo, que rompería los `const` de los
    tres esquemas publicados y el hash congelado.

- ☐ **F4.9.2.** Constantes de contrato en una sola constante de Rust cada una
  —30 ciclos, 6 ciclos explícitos, 256 de tamaño máximo de familia, 16 intentos
  de aceptación, 4 clases de eficiencia— y un **juez** que las compare con los
  máximos de los esquemas y de la interfaz.
  - Depende de: F4.9.1.
  - Procede de: precedente del juez de techos de la sección 11.4 del plan de la
    fase de población fundadora.
  - Cambia: `src/technical_profile.rs`, los esquemas nuevos, `ui/`.
  - Se comprueba con: una prueba que compare cada constante con su máximo en
    cada esquema y en la interfaz, **con control positivo de que las
    comparaciones ocurrieron**. **Debe fallar contra:** un máximo cambiado en
    un solo esquema. *Motivo:* una copia decimal silenciosa es exactamente lo
    que sobrevivió a un cambio de contrato la vez anterior.

- ☐ **F4.9.3.** Configuración del modo amplificado y sus rechazos: tabla
  `[technical_event.pcr]` admitida **solo** con `library_preparation = "pcr"` y
  obligatoria en ese caso; ciclos, eficiencia y parejas por molécula **sin
  valor predeterminado**; ciclos explícitos, tamaño máximo de familia, intentos,
  curva de composición y clases de eficiencia opcionales dentro de su dominio.
  - Depende de: F4.9.2.
  - Procede de: REQ-CONF-134/135/136 propuestos; supuestos **S-19**, **S-22**,
    **S-26**.
  - Cambia: `schemas/config-general-v4.schema.json`,
    `schemas/effective-config-v6.schema.json`, `src/configuration.rs`,
    `src/adapters.rs`, `src/application.rs`.
  - Se comprueba con: la batería de rechazos, cada uno citando parámetro y
    valor — preparación sin PCR con tabla de amplificación presente; `pcr` sin
    ciclos; 31 ciclos; eficiencia 0; parejas por molécula negativas; curva con
    eficiencia 0 o con nudos decrecientes; clases cuyos pesos no suman 1 dentro
    de la tolerancia de DTC-049. **Debe fallar contra:** un esquema con
    `additionalProperties: true` en la tabla nueva; un dominio abierto en
    cualquiera de los siete parámetros.

- ☐ **F4.9.4.** Distribución del número de copias tras los primeros ciclos, por
  programación dinámica exacta sobre los estados `1..=2^(ciclos explícitos)`, y
  su versión **sesgada por tamaño** —`p*(z) ∝ z · p(z)`—, que es la que se
  sortea, porque una familia observada es más probable cuantas más copias tenga
  su molécula.
  - Depende de: F4.9.3.
  - Procede de: diseño con PCR §5.5; Sun (1995) y Weiss y von Haeseler (1995)
    como origen del proceso de ramificación.
  - Cambia: módulo nuevo `src/pcr_amplification.rs`.
  - Se comprueba con: con 3 ciclos explícitos y eficiencia 0,5, la distribución
    de copias en `1..8` coincide con una **enumeración a mano** de los ocho
    caminos, y la sesgada por tamaño tiene media `E[Z²]/E[Z]`. **Debe fallar
    contra:** una programación dinámica que olvide el término «no se copia»;
    una que no sesgue por tamaño. *Cuidado:* la enumeración se escribe a mano;
    reusar la propia programación dinámica como oráculo sería circular.

- ☐ **F4.9.5.** Tamaño de familia: Poisson de media `λ` **condicionada a al
  menos una pareja** y **pinzada** en el máximo declarado, construida por
  multiplicación sucesiva y normalizada con la suma compensada de DTC-049, y
  sorteada con **un solo** `categorical-f64-v1`. La configuración efectiva
  publica la fracción esperada de duplicados de una biblioteca homogénea,
  `1 − (1 − e^-m)/m`, como orientación.
  - Depende de: F4.9.4.
  - Procede de: DTC-049, DTC-072; REQ-CONF-137 propuesto; supuesto **S-19**.
  - **Por qué no `poisson-count-v1`:** su muestreo por rechazo no acota el
    número de sorteos y su cola no está pinzada, y el proyecto ya fijó en
    DTC-089 y en el commit `a1f8104` que un recuento se **pinza** en vez de
    rechazarse, «porque el rechazo no acota el número de sorteos y el
    presupuesto de operaciones sí». La decisión técnica nueva debe decirlo
    expresamente para que nadie lea DTC-067 como una prohibición.
  - Cambia: `src/pcr_amplification.rs`.
  - Se comprueba con: para `λ = 0,3` y máximo 4, los pesos coinciden con los
    calculados a mano y el último recoge **toda** la masa de la cola,
    `P(K ≥ 4 | K ≥ 1)`; con `λ = 0` no se consume aleatoriedad, contado con un
    generador instrumentado; la fracción esperada publicada vale 1 % con
    `m = 0,02`, 4,8 % con `m = 0,1` y 21 % con `m = 0,5`. **Debe fallar
    contra:** unos pesos normalizados sobre `1..máximo` **sin pinzar la cola**
    —la comparación con los pesos a mano lo caza—; un `λ = 0` que consuma.

- ☐ **F4.9.6.** Aceptación por composición de bases **al elegir la molécula**:
  fracción GC del fragmento sobre A/C/G/T (las `N` no cuentan en numerador ni
  en denominador; un fragmento sin A/C/G/T recibe factor 1), factor
  `G = ((1 + e·g(GC)) / (1 + e))^c`, y un `bernoulli-v1(G)` que acepta. Si
  rechaza, se repite el sorteo del fragmento con un componente `attempt`
  añadido a la clave; **agotado el tope, se acepta el último intento y se
  cuenta**.
  - Depende de: F4.9.5.
  - Procede de: REQ-CONF-94/95 reescritos, REQ-DAT-52 reescrito; supuesto
    **S-20**; Benjamini y Speed (2012), que muestran que es el contenido del
    **fragmento completo**, no el de la lectura, lo que importa.
  - **Por qué al elegir y no al repartir copias.** Con pocas parejas por
    molécula —el caso normal— casi todas las familias tienen tamaño 1 tanto en
    regiones ricas como pobres en GC, y meter la composición en el tamaño de
    familia no cambia casi nada. La infrarrepresentación real de una región
    viene de las moléculas que **no llegan a observarse**, y eso solo puede
    expresarse al elegir cuál da lugar a la familia.
  - Cambia: `src/pcr_amplification.rs`, `src/wgs_sampling.rs`.
  - Se comprueba con: **(a)** referencia con una mitad al 30 % de GC y otra al
    70 %, curva que da factor 1 a la primera y 0,2 a la segunda, 10 ciclos y
    `m = 0,01`: la fracción de familias en la segunda mitad cae al orden de
    0,009, muy por debajo de 0,5; **(b)** con una curva que rechaza casi
    siempre y un tope de 2 intentos, el recuento de familias aceptadas al
    agotar el tope es alto, y con 16 intentos es menor. **Debe fallar contra:**
    una implementación que aplique la composición solo al tamaño de familia
    —con `m = 0,01` daría aproximadamente 0,5 (prueba a)—; un rechazo **sin**
    tope, que no termina (prueba b, ejecutada con límite de tiempo externo).
    *Cuidado:* con `m` grande la variante defectuosa también produce sesgo y no
    se distingue; el fixture **exige** `m` pequeño, y su docstring lo dice.

- ☐ **F4.9.7.** Muestreo por familias integrado en `WgsPairedEndSampler`,
  conservando su contrato externo: sigue siendo un iterador que produce
  exactamente `N` parejas con ordinales `0..N-1`. Clase de eficiencia, copias
  tras los primeros ciclos, molécula aceptada, tamaño de familia, y después una
  pareja por miembro con **identidad técnica propia**. La quimera se sigue
  decidiendo **por pareja** sobre la molécula sin modificar, de modo que su
  comparación byte a byte se conserva intacta.
  - Depende de: F4.9.6.
  - Procede de: REQ-DAT-18, REQ-DAT-47 reescrito, REQ-DAT-81 propuesto,
    DTC-087; supuestos **S-22** y **S-23**.
  - Cambia: `src/wgs_sampling.rs`, `src/pcr_amplification.rs`,
    `src/randomness.rs` (etapa `amplification` del catálogo, sin las decisiones
    de linaje, que llegan en la entrega 4).
  - Se comprueba con: **(a)** con errores de secuenciación desactivados,
    quimera a cero y sin loci, el número de parejas cuyas secuencias R1 y R2
    coinciden con otra, **contado por un oráculo que solo lee el FASTQ**, es
    igual al recuento de duplicados del manifiesto; y con `m = 0,5` ese número
    es **mayor que cero**; **(b)** las parejas de una familia no aparecen
    adyacentes en el fichero publicado, porque la escritura ordena por
    identificador de instancia, como en un FASTQ real; **(c)** añadir una
    biblioteca no cambia las familias de las demás. **Debe fallar contra:** un
    manifiesto que publique el recuento **esperado** en vez del realizado; una
    implementación que sortee el tamaño de familia pero emita una sola pareja
    (prueba a, cuya exigencia de «mayor que cero» evita que pase por el motivo
    equivocado); un ámbito indexado por orden de biblioteca (prueba c).
    *Cuidado:* con los errores de secuenciación activados los duplicados dejan
    de ser idénticos y la prueba (a) no discrimina; el fixture los desactiva y
    la docstring lo dice.

- ☐ **F4.9.8.** Recuentos agregados del manifiesto: familias, parejas, parejas
  duplicadas, fracción realizada de duplicados, histograma de tamaños de
  familia, familias truncadas al final de la biblioteca, intentos de
  composición totales y agotados; y la **cobertura de moléculas distintas**,
  que es la que un analista ve tras eliminar duplicados. Nada por pareja.
  - Depende de: F4.9.7.
  - Procede de: REQ-MAN-45 reescrito, REQ-WGS-11, REQ-DAT-56; supuesto **S-28**.
  - Cambia: `src/event_manifest.rs`, `schemas/manifest-event-v3.schema.json`.
  - Se comprueba con: la suma del histograma por su índice es igual al número
    de parejas; la fracción realizada es `1 − familias/N`; la cobertura de
    moléculas distintas es la básica multiplicada por `familias/N`; validación
    externa del manifiesto. **Debe fallar contra:** un histograma que cuente
    familias en vez de parejas en alguna celda; una cobertura de moléculas
    distintas calculada desde la fracción esperada en vez de la realizada.

- ☐ **F4.9.9.** Cota de recursos **publicada, no impuesta**: operaciones y
  memoria derivadas de la configuración y de la geometría conocida antes de
  abrir el primer ámbito aleatorio, con aritmética entera, publicadas en el
  manifiesto como valores **calculados sin ejecutar** y no como medidas.
  - Depende de: F4.9.8.
  - Procede de: REQ-MAN-94 propuesto; supuesto **S-27**; diseño con PCR §6.
  - **Qué se publica.** Por encima del coste que el recorrido ya tiene,
    `ops ≤ N × [22 211 + 16·P + 60·T] + N × X_pair + 98 304`, donde `P` es el
    número de copias haplotípicas, `T` una cota del número de tractos que caben
    en un inserto y `X_pair` el coste por pareja que este trabajo no modifica;
    y una memoria añadida **independiente de `N` y del genoma**. Se declara
    además que el recorrido regenera cada biblioteca **dos veces**, de modo que
    el coste efectivo es el doble.
  - **Qué no se acota, dicho expresamente.** El número de palabras de generador
    de los sorteos por rechazo que el recorrido ya usa —la normal truncada del
    inserto y el índice uniforme— no tiene cota superior en el peor caso, y
    DTC-088 y el commit `a1f8104` ya lo reconocen. Este trabajo no añade ningún
    sorteo por rechazo nuevo: todas sus distribuciones son categóricos de una
    palabra o ensayos de Bernoulli. Sí multiplica por el tope de intentos los
    sorteos de inserto de una familia.
  - Cambia: `src/pcr_amplification.rs`, `src/event_manifest.rs`.
  - Se comprueba con: un contador de unidades instrumentado, solo en pruebas,
    sobre un fixture con intentos, tamaño máximo de familia, ciclos explícitos
    y número de tractos **en sus máximos de contrato** y una curva que rechace
    casi siempre: las unidades reales son menores o iguales que la cota, con
    control positivo de que el contador contó. **Debe fallar contra:** un
    mutante que suprima el término de los intentos de la fórmula. *Cuidado:*
    con curva neutra el tope de intentos **no se ejercita** y el mutante
    sobrevive; sin el fixture saturado la prueba no discrimina, y la docstring
    lo dice.

- ☐ **F4.9.10.** Reducción al caso sin amplificar, y **rechazo temporal** de la
  combinación de preparación amplificada con loci de repetición declarados.
  - Depende de: F4.9.9.
  - Procede de: REQ-DAT-84 propuesto; supuesto **S-16**; restricción aritmética
    de §3.1.
  - Cambia: `src/configuration.rs`, `src/diagnostics.rs`.
  - **El rechazo, y por qué.** Mientras la entrega 4 no exista, una biblioteca
    amplificada sobre un escenario con loci declarados tartamudearía **solo con
    la parte de base**, y por tanto por debajo de lo medido. Antes que publicar
    una cifra que se sabe baja, la configuración **rechaza** esa combinación con
    un diagnóstico que nombra el paso pendiente. El rechazo se levanta en
    F4.10.6.
  - Se comprueba con: **(a)** con preparación amplificada, parejas por molécula
    a cero, curva neutra y tasas a cero, los FASTQ son **byte a byte** iguales
    a los de la preparación sin PCR, y el manifiesto difiere solo en los campos
    declarativos; **(b)** un generador instrumentado confirma **cero consumo**
    en la etapa de amplificación con esos parámetros; **(c)** preparación
    amplificada más un locus declarado se rechaza citando el paso F4.10.
    **Debe fallar contra:** una implementación que consuma aleatoriedad con
    parámetros neutros o que desplace los ordinales de plantilla (pruebas a y
    b); una que acepte en silencio la combinación y publique una fracción de
    tartamudeo entre un 14 % y un 28 % por debajo de la medida (prueba c).

- ☐ **F4.9.11.** Validador independiente de la amplificación: rederiva desde la
  semilla y la configuración efectiva las familias, las copias de los primeros
  ciclos y los tamaños de familia, **sin llamar al productor**, y coincide con
  los recuentos del manifiesto.
  - Depende de: F4.9.10.
  - Procede de: REQ-VAL-58 propuesto, arquitectura §11.
  - Cambia: validador nuevo bajo `tests/`.
  - Se comprueba con: coincidencia con el manifiesto sobre dos fixtures;
    `./scripts/check-independent-validators.sh evidence/pcr-f4-9` da **5 PASS**.
    **Debe fallar contra:** un oráculo que lea el histograma publicado en vez
    de rederivarlo —se detecta mutando el histograma del manifiesto y exigiendo
    que el validador lo rechace. *Nota:* la independencia algebraica de la
    Poisson pinzada la da la prueba de pesos a mano de F4.9.5, no este oráculo.

- ☐ **F4.9.12.** Interfaz Marimo de la preparación amplificada: selector de
  preparación, ciclos, eficiencia, parejas por molécula, ciclos explícitos,
  tamaño máximo de familia, intentos, curva de composición y clases de
  eficiencia; y la **fracción esperada de duplicados** mostrada al usuario
  según las parejas por molécula que elija.
  - Depende de: F4.9.11.
  - Procede de: REQ-UI-01/03/04/05, hoja de ruta §7; REQ-CONF-137 propuesto.
  - Cambia: `ui/config_editor.py`, `ui/ui_panels.py`; prueba nueva en
    `tests_ui/`.
  - Se comprueba con: ida y vuelta completa; la tabla de amplificación no
    aparece con preparación sin PCR; los límites de los controles coinciden con
    las constantes de contrato, comprobado por el juez de F4.9.2;
    `marimo check --strict` sin avisos. **Debe fallar contra:** un control cuyo
    máximo no coincida con la constante; una tabla de amplificación editable
    con preparación sin PCR.

> **Revisión del hito de la entrega 3.** Comprobar familias, duplicados
> contados por bytes, composición al elegir la molécula, tope de intentos, cota
> publicada, reducción al caso sin PCR y el rechazo temporal. Ejecutar la
> batería completa con
> `./scripts/check-independent-validators.sh evidence/pcr-f4-9`. Aplicar el
> cambio de contrato de la sección 9 en su parte de preparación amplificada y
> **repetir la evidencia de cierre de la fase de realismo técnico con
> preparación sin PCR**, que debe dar los mismos bytes. No continuar hasta
> confirmar el resultado.

---

## 8. Entrega 4 — Tartamudeo clonal heredable por linaje (F4.10)

Resultado: con preparación amplificada, el deslizamiento puede ocurrir en
**cualquier ciclo**, y una copia que resbala transmite su número de unidades a
**toda su descendencia**. Un resbalón en el primer ciclo alcanza a media
familia; uno en el último, a una sola copia. El resultado observable es que una
familia de duplicados puede mostrar dos o tres alelos distintos en proporciones
que imitan a un heterocigoto, cuando procede de una sola molécula.

- ☐ **F4.10.1.** Extensión del catálogo de ámbitos aleatorios con las
  decisiones de linaje: progenitor, ciclo de divergencia, resbalón por ciclo y
  tracto, y paso; indexadas por **identificador de locus**, nunca por posición.
  - Depende de: cierre de las entregas 2 y 3.
  - Procede de: REQ-DAT-85 propuesto, DTC-081, DTC-086; supuesto **S-11**.
  - Cambia: `src/randomness.rs`.
  - Se comprueba con: vectores congelados; las etapas ya existentes producen
    exactamente las mismas subsemillas que antes; **añadir un locus en otro
    cromosoma no cambia los resbalones de los loci previos**. **Debe fallar
    contra:** un ámbito indexado por la posición del tracto dentro de la
    molécula —la prueba del locus añadido en otro cromosoma es la única que lo
    caza.

- ☐ **F4.10.2.** Genealogía mínima de la familia: la primera pareja tiene un
  linaje propio completo de un paso por ciclo; la pareja `i` elige un
  **progenitor uniforme** entre las anteriores y un **ciclo de divergencia**
  geométrico con `q = 1/(1+e)`, hereda los resbalones del progenitor en los
  ciclos anteriores a la divergencia y sortea los posteriores.
  - Depende de: F4.10.1.
  - Procede de: REQ-DAT-82 propuesto; supuesto **S-21**; diseño con PCR §5.7.
  - **Declarado como aproximación.** Es una aproximación del coalescente del
    proceso de ramificación, no el coalescente exacto; no se afirma tener una
    forma cerrada de este último. La forma geométrica refleja que dos
    descendientes de una amplificación eficiente divergen casi siempre pronto.
  - Cambia: módulo `src/pcr_amplification.rs` o hermano.
  - Se comprueba con: **fixture saturado** — un locus dentro de la molécula,
    `m = 8` para forzar familias grandes, divergencia forzada al último ciclo
    para todas las parejas y tasa de resbalón 1 en el primer ciclo y 0 en los
    demás (sustitución por ciclo, solo en prueba): **todas** las parejas de
    cada familia muestran el **mismo** número de unidades, distinto del
    verdadero. Segunda variante con divergencia forzada a cero: las parejas son
    **independientes**. **Debe fallar contra:** un tartamudeo sorteado por
    pareja o por lectura; una herencia que copie los resbalones **posteriores**
    a la divergencia en vez de los anteriores. *Cuidado, en la docstring:* sin
    forzar la divergencia y sin tasa 1, algunas discordancias dentro de la
    familia son legítimas y la aserción no discrimina.

- ☐ **F4.10.3.** Derivación de la tasa de resbalón por ciclo desde la
  calibración: `r_pcr = 1 − (1 − r_con)/(1 − r_sin)` y
  `p = 1 − (1 − r_pcr)^(1/c)`, con los tres insumos y la fórmula publicados en
  la configuración efectiva, y sustitución directa admitida. Si `r_con ≤ r_sin`,
  la tasa vale cero.
  - Depende de: F4.10.2.
  - Procede de: REQ-CONF-138 propuesto; supuesto **S-24**; calibración §3.2.
  - Cambia: `src/pcr_amplification.rs`, `src/technical_profile.rs`.
  - Se comprueba con: con diez ciclos, las tasas por ciclo valen 0,0158 para
    motivo de 1 pb; 0,00276 para 2 pb; 0,000815 para 3 pb; 0,000532 para 4 pb;
    0,000593 para 5 pb y 0,00131 para 6 pb, comparadas con el cálculo a mano a
    tres cifras; la configuración efectiva publica `r_con`, `r_sin`, el número
    de ciclos y la fórmula; una sustitución directa del usuario prevalece y se
    marca como suya. **Debe fallar contra:** una implementación que **sume** la
    columna sin PCR a la columna con PCR, que es el error contrario y cuenta el
    mismo fenómeno dos veces; una que use `r_con` directamente como tasa por
    ciclo.

- ☐ **F4.10.4.** Aplicación en el recorrido: sobre la molécula ya quimérica si
  lo es, **primero** el paso neto heredado del linaje de PCR y **después** la
  decisión de base por pareja de la entrega 2. Los tractos se procesan en orden
  canónico de inicio sobre la molécula. Las guardas de viabilidad se aplican al
  paso **neto acumulado**, excluyendo antes de sortear y sin reintentar.
  - Depende de: F4.10.3.
  - Procede de: REQ-DAT-76 y REQ-DAT-82 propuestos; supuestos **S-13**,
    **S-14**, **S-15**.
  - Cambia: `src/stutter.rs`, `src/pcr_amplification.rs`, `src/execution.rs`.
  - Se comprueba con: el orden se observa sobre un fixture con paso de PCR
    fijado y tasa de base 1 con paso fijado, cuyo resultado es la **suma** de
    los dos; un paso neto que dejaría el tracto sin ninguna unidad completa
    queda excluido y no consume; un corpus con tasa de quimera alta no produce
    `ChimeraError::PairSourceMismatch`. **Debe fallar contra:** un orden
    invertido que aplique la base antes que la de PCR y luego recalcule la
    viabilidad sobre el tracto equivocado; una guarda aplicada a cada paso por
    separado en vez de al neto.

- ☐ **F4.10.5.** **Comprobación aritmética de la composición**, que es la razón
  de ser del orden de este plan: con las dos partes activas, la fracción de
  lecturas afectadas por tartamudeo que un oráculo cuenta sobre los FASTQ
  reproduce la fracción **medida** con PCR, `1 − (1−r_sin)(1−r_pcr) = r_con`,
  dentro de la tolerancia muestral declarada.
  - Depende de: F4.10.4.
  - Procede de: calibración §3.2; restricción aritmética de §3.1.
  - Cambia: prueba de sistema nueva bajo `tests/`.
  - Se comprueba con: un corpus con un homopolímero declarado y un
    dinucleótido, número de moléculas suficiente para que el intervalo de
    confianza sea estrecho, y semilla congelada: la fracción observada con
    preparación amplificada cae dentro de la tolerancia de 0,17 para el
    homopolímero y 0,038 para el dinucleótido; con preparación sin PCR, dentro
    de 0,027 y 0,011. **Debe fallar contra:** una entrega de la preparación
    amplificada **sin** la parte de base, que daría 0,147 y 0,0273 —un 14 % y
    un 28 % por debajo—; una implementación que sume las dos columnas, que
    daría por encima. *Cuidado:* la tolerancia se calcula desde el número de
    moléculas del corpus y se escribe en la docstring; una tolerancia holgada
    dejaría pasar las dos variantes defectuosas y la prueba no discriminaría.

- ☐ **F4.10.6.** Cierre: recuentos de resbalones y de resbalones heredados por
  locus en el manifiesto; validador independiente que rederiva los pasos netos
  de linaje desde la semilla; control en la interfaz de la tabla de tartamudeo
  con PCR; y **levantar el rechazo temporal** de F4.9.10.
  - Depende de: F4.10.5.
  - Procede de: REQ-MAN-45 y REQ-VAL-58 propuestos, hoja de ruta §7.
  - Cambia: `src/event_manifest.rs`, `ui/`, validador bajo `tests/`,
    `src/configuration.rs`.
  - Se comprueba con: los recuentos de resbalones heredados son mayores que
    cero en el fixture de familias grandes y cero con familias de tamaño 1;
    el validador coincide sin llamar al productor;
    `./scripts/check-independent-validators.sh evidence/clonal-f4-10` da **5
    PASS**; la combinación antes rechazada ahora se acepta y produce la
    fracción de F4.10.5; ida y vuelta de la interfaz. **Debe fallar contra:**
    un recuento de resbalones heredados que cuente todos los resbalones; un
    levantamiento del rechazo sin la prueba aritmética de F4.10.5 verde.

> **Revisión del hito de la entrega 4.** Comprobar la genealogía sobre el
> fixture saturado, la derivación de tasas, el orden de aplicación, la
> composición aritmética y el levantamiento del rechazo. Ejecutar la batería
> completa con `./scripts/check-independent-validators.sh
> evidence/clonal-f4-10`. Cerrar el cambio de contrato de la sección 9. No
> continuar hasta confirmar el resultado.

---

## 9. Cambio de contrato posterior al cierre de la fase de realismo técnico

Con la forma de la sección 11 del plan de la fase de población fundadora. **No
es un añadido silencioso y no debe presentarse como tal.**

### 9.1. Qué **no** se reabre

Los hitos F4.0 a F4.5 y sus evidencias de cierre no se reabren. Ninguna
decisión de la fase se revoca. Lo que cambia es un contrato posterior y
acotado, por decisión expresa del usuario del 19 de septiembre de 2026:
simular, además de la biblioteca sin amplificar que ya existe, una **biblioteca
amplificada por PCR**, y reconocer el **tartamudeo** como artefacto del perfil
vigente.

### 9.2. Qué cambia, y por qué no es una contradicción

- **El tartamudeo sin PCR no contradice la preparación sin PCR.** Los cinco
  requisitos que declaran la preparación exclusivamente sin PCR hablan de la
  **preparación de biblioteca**. El deslizamiento de la polimerasa ocurre
  además en la amplificación en puente sobre la celda y en la síntesis de la
  secuenciación, que existen en el perfil vigente. El requisito nuevo debe
  decirlo, para que nadie lo lea como una grieta.
- **La preparación amplificada sí reabre una ampliación aplazada.** DTF-056 se
  **divide**: salen de ella la preparación amplificada, los duplicados de PCR,
  el sesgo de composición **configurable**, la sobrerrepresentación estocástica
  y el tartamudeo heredable; **permanecen** el perfil de sesgo GC calibrado con
  datos reales, los duplicados ópticos y de agrupación, y la combinación con
  excepciones regionales. DTF-055, el marcado de duplicados en BAM y CRAM, no
  cambia.

### 9.3. Dónde vive cada literal que hay que tocar

Una copia silenciosa de un literal es lo que sobrevive a un cambio de contrato.
Los sitios conocidos, comprobados en `11d26e5`:

| Sitio | Papel | Qué dice hoy |
|---|---|---|
| `src/technical_profile.rs`, líneas 28–30 | canónico | `LibraryPreparation` con una sola variante, `PcrFree` |
| `src/application.rs`, líneas 776–786 | constantes | duplicados de PCR no aplicables y ópticos a cero |
| `src/event_manifest.rs`, líneas 1516–1531 | contrato | exige `pcr.applicable == false`, `gc.mode == "neutral"`, preparación sin PCR |
| `schemas/effective-config-v4.schema.json`, líneas 312–313 | contrato | `const` del identificador del perfil, de su ruta, de su SHA-256, de la preparación, de `pcr.applicable` y de `gc.mode` |
| `schemas/manifest-event-v2.schema.json`, líneas 115–118 y 157 | contrato | los mismos `const` |
| `schemas/config-general-v2.schema.json`, `$defs/technical_event` | contrato | objeto cerrado de once claves; ninguna nombra la PCR, de modo que hoy una clave de PCR se rechaza **como campo desconocido**, no como capacidad aplazada |
| `ui/config_editor.py`, `ui/ui_panels.py` | contrato | controles del evento técnico |
| `resources/f4/gps-hiseqx-pcrfree-pe150-v1.json` | recurso | perfil único, con su hash congelado |

**Regla:** los esquemas publicados **no se editan**. Las versiones nuevas
—general v4, efectivo v6, manifiesto de evento v3— llevan la capacidad, y los
`const` de las versiones heredadas se comparan byte a byte en cada hito. Un
fichero general v2 sigue publicando efectivo v4 y manifiesto de evento v2.

### 9.4. Evidencia que hay que repetir

Al cerrar la entrega 2 y otra vez al cerrar la entrega 3, se repite la
evidencia de cierre de la fase de realismo técnico **con preparación sin PCR y
sin loci de repetición declarados**, y debe dar los mismos bytes que la línea
base de F6.4.0. Esa repetición es la que sostiene la afirmación de que el
contrato se ha ampliado sin romperse.

### 9.5. Dónde se registra

Sección **16 nueva** del plan de la fase de realismo técnico, titulada «Cambio
de contrato posterior al cierre — tartamudeo y preparación amplificada
(fecha)», que recoja: qué no se reabre, qué cambia y por decisión de quién,
dónde vive cada literal, las pruebas nuevas con sus variantes defectuosas, los
fixtures que se reabren y la repetición de la evidencia de cierre. Y una nota
fechada en la sección de la fase de realismo técnico de la hoja de ruta.

---

## 10. Correcciones que este plan **propone** a otros documentos

Ninguna está aplicada. Este plan no toca ningún otro fichero.

| Documento | Dónde | Qué | Origen |
|---|---|---|---|
| `especificaciones-0.1.md` | §4.3, glosario de la repetición en tándem | Interrupciones declaradas, pureza realizada; rechazo de la pureza como entrada | S-06 |
| `especificaciones-0.1.md` | §4.3, mismo párrafo | Las unidades se ganan o se pierden por el **extremo derecho** del tracto | S-07 |
| `especificaciones-0.1.md` | §4.4, frase de orden normativo, y REQ-DAT-17 | Insertar el tartamudeo entre la quimera y calidad y errores, y las familias de amplificación antes de la quimera | S-13, diseño con PCR §8 |
| `especificaciones-0.1.md` | §4.4, tras REQ-DAT-73 | Añadir REQ-DAT-74 a 80 (tartamudeo) y REQ-DAT-81 a 85 (amplificación) | Expediente §10.2, diseño §9.2 |
| `especificaciones-0.1.md` | **REQ-DAT-80** | **Reescribir**: la identidad byte a byte vale para los escenarios **que no declaran loci de repetición** | **C-1** |
| `especificaciones-0.1.md` | §4.1, tras REQ-CONF-130 | Añadir REQ-CONF-131 a 133 (tartamudeo) y 134 a 138 (amplificación) | Expediente §10.2, diseño §9.2 |
| `especificaciones-0.1.md` | **REQ-CONF-132** | **Reescribir**: lo desactivado por omisión es el tartamudeo **heredable de PCR**; el de base se activa al declarar loci | **C-1** |
| `especificaciones-0.1.md` | REQ-CONF-88/89/94/95, REQ-TEC-07/08, REQ-DAT-46/47/52/53, REQ-MAN-45/48/60 | Las redacciones de la sección 9.1 del diseño con PCR | Diseño §9.1 |
| `especificaciones-0.1.md` | §4.7, §4.6, §7 | Añadir REQ-MAN-93/94, REQ-VAL-57/58, ACE-118 y ACE-119; reescribir ACE-49/51/68 | Ambos documentos |
| `decisiones-tecnicas-0.1.md` | Tabla | **DTC-090** (provisional): `stutter-v1`, catálogo `f4-technical-v4`, y la aclaración de que DTC-067 no prohíbe usar `bernoulli-v1` ni `categorical-f64-v1` aquí | S-11, S-14 |
| `decisiones-tecnicas-0.1.md` | Tabla | **DTC-091** (provisional): `pcr-amplification-v1`, la programación dinámica de los primeros ciclos, la Poisson pinzada como categórico de una palabra, el rechazo con tope, la genealogía, las constantes de contrato y su juez, y la derivación de la tasa por ciclo | S-19 a S-28 |
| `decisiones-tecnicas-0.1.md` | Tabla | Decisión sobre la realización de interrupciones y la pureza como salida | S-06 |
| `plan-implementacion-f6-0.1.md` | §8, F6.4.3 | Corregir: **ya no depende** del primer paso del hito de variantes estructurales | S-02, C-4 |
| `plan-implementacion-f6-0.1.md` | §8, F6.4.2 | Corregir: **no depende** de derivar una referencia estructural | S-03, C-4 |
| `plan-implementacion-f6-0.1.md` | §3 | DP-F6-03, DP-F6-06 y DP-F6-07 quedan cerradas **en su parte de repeticiones** por S-08, S-09 y S-06; DP-F6-04 anota que el catálogo técnico nuevo lleva las etapas `stutter` y `amplification` | S-06/08/09/11 |
| `plan-implementacion-f6-0.1.md` | §8 | Añadir un paso F6.4.10 «Tartamudeo heredable con PCR», que remite a la entrega 4 de este plan | Diseño §9.5 |
| `plan-implementacion-f4-0.1.md` | Sección 16 nueva | El cambio de contrato de la sección 9 de este plan | S-18 |
| `hoja-de-ruta-0.1.md` | Fases de realismo técnico y de variación compleja | Registrar el adelanto de los pasos de repeticiones y la reapertura acotada del contrato técnico, con fecha | S-18 |
| `registro-ampliaciones-futuras-0.1.md` | DTF-056 | Dividir según §9.2; añadir **DTF-059**: errores de polimerasa heredables, quimeras de PCR heredadas, sesgo de longitud y genealogía exacta | Diseño §9.4 |
| `calibracion-tartamudeo-0.1.md` | §2 y §5 | Anotar que la columna con PCR aplica a la preparación amplificada y se convierte en tasa por ciclo con la fórmula de §3.2; y **resolver C-2**, dejando dicho que el mecanismo adoptado es el de la ficha D-14 y que la forma de tres parámetros se conserva como referencia del sesgo, no como estructura de datos | **C-2** |
| `expediente-decisiones-repeticiones-tandem-0.1.md` | D-14 | Anotar que su premisa «no existe ninguna calibración en el repositorio» **ya no es cierta**, y que por tanto su recomendación de valores cero queda sustituida | **C-2** |
| `expediente-decisiones-repeticiones-tandem-0.1.md` | D-16 y D-17 | Anotar las correcciones **C-1** y **C-6** | C-1, C-6 |
| `expediente-decisiones-repeticiones-tandem-0.1.md` | Parte B, encabezado | Identificar **cuáles** son las tres confirmaciones: hoy solo D-10 lleva la marca | Aviso de recuento, §2 |
| `trazabilidad-requisitos-0.1.md` | Matriz | Filas nuevas para los requisitos nuevos; las filas que hoy apuntan a DTF-056 pasan a la sección 16 de la fase de realismo técnico | Ambos documentos |

---

## 11. Matriz de cierre

| Entrega | Resultado construido | Evidencia | Riesgo protegido |
|---|---|---|---|
| **1** (F6.4 adelantado) | Loci declarados, alelos por copia que difieren entre individuos, verdad v3, VCF 4.5 con campos de repetición, manifiesto e interfaz | Diploide 5/7 con secuencias, verdad y VCF comprobados; flancos byte a byte iguales con interrupción; umbral de 50 pb; hashes sin loci iguales a la línea base; validadores independientes 5 PASS | Que el alelo pierda el número de unidades; que la realización recoloque las interrupciones; que una ejecución sin repeticiones cambie de bytes |
| **2** (F4.8) | Tartamudeo de base después de la quimera, por molécula y no por lectura, con tasas calibradas | Fixture saturado con el tracto dentro del solapamiento; corpus con quimera alta sin `PairSourceMismatch`; verdad idéntica y FASTQ distintos; consumo cero sin loci | Que el par salga incoherente; que el tartamudeo se aplique antes de la quimera y reviente la comparación byte a byte; que cambie la verdad |
| **3** (F4.9) | Familias de amplificación, duplicados, composición al elegir la molécula, cota publicada | Duplicados contados por un oráculo que solo lee el FASTQ; fracción de familias en la mitad rica en GC; tope de intentos; contador de unidades bajo la cota; identidad byte a byte con parámetros neutros | Que los duplicados se cuenten esperados en vez de realizados; que la composición se aplique al tamaño de familia y no a la elección; que el rechazo no termine |
| **4** (F4.10) | Tartamudeo heredable por linaje, y la fracción compuesta que reproduce lo medido | Fixture con divergencia forzada; tasas por ciclo contra el cálculo a mano; fracción observada dentro de la tolerancia de 0,17 y 0,038 | Que el tartamudeo clonal se sortee por pareja; que se sumen las dos columnas y se cuente dos veces; que la biblioteca amplificada subestime entre un 14 % y un 28 % |

---

## 12. Lo que queda fuera a propósito, con destino

| Queda fuera | Destino |
|---|---|
| Duplicaciones, deleciones estructurales, inversiones y translocaciones | Hitos de referencias estructurales y de variantes estructurales de la fase de variación compleja |
| Referencias estructurales derivadas, su mapa y su linaje | Mismo sitio |
| Ploidía generalizada por secuencia | Hito de ploidía generalizada |
| **Modelo escalonado simple de mutación**, que hace que las repeticiones **muten entre generaciones** | Paso F6.4.4. El usuario pidió que **varíen entre individuos**, que se consigue declarando alelos distintos; son cosas diferentes |
| Mutaciones internas de pureza, interrupciones o motivo | Paso F6.4.5. REQ-CONF-63 las exige desactivadas por omisión, y desactivadas no cuestan nada |
| Generación aleatoria de loci por regiones | Paso F6.4.6. Con loci declarados no se consume aleatoriedad ninguna |
| Frecuencias alélicas por locus y población, que es lo que el usuario acabará usando | Segunda mitad de S-04, inmediatamente después de este trabajo; es la salida del límite C-5 |
| Fichero complementario del VCF y catálogo de códigos de pérdida | Hito de variantes complejas. No hace falta: los campos del estándar 4.5 no pierden detalle |
| Repeticiones cuyo alelo se separa 50 pb o más de la referencia | Hito de variantes complejas, con representación simbólica |
| Conjuntos y matrices por referencia de mapeo, BAM, CRAM y mapeo real | Hito de conjuntos y matrices por referencia |
| Fundadores de la fase de población inicial con alelos de repetición | El corte usa la vía directa |
| **Errores de polimerasa heredables** por familia | DTF-059 nuevo. Kebschull y Zador (2015) midieron que apenas afectan a la distribución, lo que respalda aplazarlos |
| **Quimeras formadas durante la PCR y heredadas** por la familia | DTF-059. La misma fuente las encuentra raras y limitadas a números de copia bajos |
| Sesgo de longitud de la amplificación y genealogía exacta del proceso de ramificación | DTF-059 |
| Duplicados ópticos y de agrupación, y marcado de duplicados en BAM/CRAM | DTF-056 y DTF-055, sin cambios |
| Perfil de sesgo GC **calibrado con datos reales** | DTF-056, que conserva esa parte |
| **Mosaicismo somático** modelado aparte | No se modela: se queda dentro de la tasa sin PCR, porque la medición no permite separarlo. Si algún día se modelara, pertenecería a la **verdad** del individuo y no a las lecturas, y por tanto no cabe en una fase cuyo contrato prohíbe que los eventos técnicos alteren la verdad |
| Imposición de la cota de recursos como rechazo previo | Fase de presupuesto de recursos. Aquí se **publica**, no se impone (S-27) |
| Optimizar el coste del tope de intentos y el factor dos de las dos pasadas del recorrido | Fase de rendimiento. Se declaran y no se optimizan |
| La deuda de `tests/test_process_adapter.py`, que aborta bajo la puerta canónica desde el commit `6139fe9` | Sigue registrada en el plan de la fase de población fundadora, sin destino asignado. Este plan no la toca ni la cita como puerta |
