# Plan de implementación — indels pequeños con frecuencia poblacional 0.1

## 1. Estado, autoridad y alcance

- **Fecha:** 27 de octubre de 2026.
- **Estado:** **propuesta, no ejecutable todavía.** Antes de escribir código
  hay que aprobar los cambios normativos de la sección 4 y cerrar las
  decisiones pendientes de la sección 3. Ningún paso de la sección 6 está
  hecho.
- **Qué es.** Una **reapertura acotada de F2.2** —la matriz poblacional de
  variantes y su generador—, por decisión del responsable del 27/10/2026, con
  la misma forma que las reaperturas de F1 (sección 14 de
  `plan-implementacion-f1-0.1.md`) y de F4 (sección 16 de
  `plan-implementacion-f4-0.1.md`), y con el registro de bytes de la sección
  14 de `plan-implementacion-repeticiones-pcr-0.1.md`. Toca además, por
  decisión expresa, los contratos de la población fundadora sintética (F5),
  cerrada.
- **Numeración.** Los pasos heredan la coordenada F2.2 y continúan después del
  último existente, F2.2.12 (`plan-implementacion-f2-0.1.md` §7). Los pasos de
  la entrega que toca F5 también son F2.2.x, porque la decisión 4 del
  responsable sitúa todo el trabajo en F2.2; cada uno cita además el hito F5
  cuyo contrato modifica.
- **Prohibiciones heredadas.** Este plan no autoriza *commit* ni publicación.
  No revierte trabajo ajeno. No se etiqueta una puerta con un comando que no se
  haya ejecutado. Las copias de seguridad de cualquier mutación van fuera de
  git; nunca `git checkout --`. Nada de bucles de sondeo.

### 1.1. Qué pide el responsable

Hoy la matriz poblacional (`[population_snvs]`, recurso TSV+JSON de DTC-054)
solo admite cambios de una base. Los indels existen escritos a mano por
individuo (`[[direct_events]]`) o en un generador aleatorio privado por
individuo, sin exponer y aplazado al subpaso F6.5.a. ACE-10 dice que los SNV
iniciales salen de la matriz y los indels de una cantidad media o densidad.

El responsable quiere **indels pequeños con frecuencia en la población,
compartidos entre individuos como los SNV**. Sus decisiones del 27/10/2026,
que este plan toma como fijas:

1. **Matriz suministrada y generada** (y reutilizada, como hoy).
2. **Un indel puede solaparse con SNV en las mismas posiciones**; los solapes
   no se rechazan.
3. **Los indels entran también en la población fundadora de F5**, que hoy por
   contrato solo tiene SNP bialélicos (REQ-EVO-42, REQ-EVO-74, DTC-089).
4. **Es una ampliación de F2.2**, no un hito de F6.
5. **Tamaño máximo:** menos de 50 pb de secuencia nueva, con el umbral
   configurable de REQ-CONF-55; un indel de la matriz tiene como mucho
   umbral − 1 bases (49 por omisión).

### 1.2. Qué **no** hace este plan

| Fuera | Destino |
|---|---|
| Hacer configurable el umbral de REQ-CONF-55 | F6.2, donde lo sitúa `trazabilidad-requisitos-0.1.md` (fila REQ-CONF-55, «Planificada»). Aquí se lee de un único sitio que hoy vale 50 (`src/indel_application.rs:23`) |
| Indels aleatorios privados por individuo (cantidad media o densidad) | F6.5.a, sin cambio |
| Alelos complejos en la matriz (`REF=ACG`, `ALT=T`, MNP) | Fuera; solo SNV, inserción pura anclada y deleción pura anclada |
| Preferencia por homopolímeros o repeticiones en el generador | Ampliación futura nueva (DTF), ver DP-IND-07 |
| Conversión de VCF poblacionales a la matriz | DTF-038, sin cambio |
| Solapes dentro de los bloques de F5 | Ver DP-IND-11 |
| Indels de la matriz combinados con inserciones y deleciones directas en la misma secuencia | Ver DP-IND-03; queda rechazado y su levantamiento, sin destino propuesto |
| Indels de matriz junto a un tracto insertado en referencia suministrada | Ya rechazado hoy para matriz suministrada o reutilizada (`src/application.rs:2694-2711`); no cambia |

## 2. Lo que el código hace hoy, comprobado

Cada afirmación cita el sitio leído. Lo que no se ha ejecutado se dice.

- **La matriz solo admite una base en `ref` y `alt`.** `parse_snv_base`
  (`src/population_snvs.rs:1970`) rechaza cualquier otra cosa. La cabecera del
  TSV tiene ocho columnas fijas (`src/population_snvs.rs:28`).
- **La multialelia suministrada ya existe.** Varias filas con la misma
  `contig`, `position` y `ref` forman un locus; la suma ALT por población se
  normaliza con DTC-049 (`src/population_snvs.rs:695-990`).
- **El genotipo se sortea por locus `(contig, position, ref)`, con
  Hardy-Weinberg o `F`, y es multialélico de serie.** `canonical_loci` y
  `sample_initial_genotype_compositions` (`src/initial_genotypes.rs:448-668`)
  pasan un vector de frecuencias de longitud arbitraria a
  `hwe_composition` o `f_diploid_genotype`. El ámbito de azar lleva `ref` y la
  lista de ALT unida por comas (`initial_genotype`, catálogo de ámbitos §3.1 de
  `decisiones-tecnicas-0.1.md`).
- **El identificador de variante ya admite cadenas.** `derive_matrix_variant_id`
  (`src/identity.rs:1388`) convierte `ref` y `alt` en texto y los codifica con
  longitud; un SNV conserva su identificador si `ref` y `alt` pasan a ser
  cadenas de una base.
- **Los SNV de la matriz se aplican antes que los eventos directos.**
  `integrate_initial_individual` (`src/initial_individual_integration.rs:215`)
  monta la lista «SNV de la matriz, después eventos directos» y llama al
  compositor. Cada alelo de la matriz es un evento **por copia**, con
  identificador derivado de individuo, población, variante y fase
  (`derive_matrix_event_id`, línea 583).
- **Los cambios de base se anclan a la referencia; las inserciones y
  deleciones directas, a las coordenadas vigentes.** Lo dice
  `src/execution.rs:2752-2755` y lo confirma `with_snp`
  (`src/segments.rs:545`), que busca la posición en los segmentos de fuente.
- **La verdad ya convive con una copia cuyo locus SNP borró una deleción.** El
  caso SNP omite esa copia (`src/truth_projection.rs:771`) y la deleción tiene
  su caso estructural propio.
- **La verdad exige un caso SNP por fila de la matriz**, aunque sea todo-REF
  (`schemas/truth-v4.md`, «El conjunto de casos»).
- **El VCF ya funde las variantes cuyo `REF` comparte una base, contando el
  ancla**, de forma transitiva (commit `6bef50e`).
- **Un locus de microsatélite rechaza cualquier fila de la matriz dentro de su
  tracto**, en la frontera de configuración (`admit_repeat_loci_against_small_variants`,
  `src/execution.rs:2718-2790`), solo por la posición de la fila.
- **El generador de SNV no excluye los tractos declarados.**
  `select_generated_loci` y `candidate_weight`
  (`src/population_snvs.rs:3265-3349`) solo miran A/C/G/T y las regiones. Ver
  la contradicción C-2.
- **Las versiones publicadas se eligen así:** documento general v2 → efectiva
  v4, v4 → efectiva v6 (`DeclaredDocumentVersion`, `src/application.rs:105-130`);
  manifiesto general v4, o v6 si hay loci de repetición
  (`src/results.rs:68-85` y `1373-1381`); la verdad, v4 o v5 por contenido.
  F5 elige su catálogo de azar por contenido: v2 solo si el perfil declara la
  distribución de H (`src/f5_publication.rs:136-158`).
- **Las estadísticas por locus de la configuración efectiva v6 solo admiten
  una base.** `initial_population_statistics_locus.ref` e
  `initial_population_statistics_allele.allele` son el enumerado A/C/G/T
  (`schemas/effective-config-v6.schema.json:476` y `539`).
- **La verdad de F5 solo tiene casos SNP.** `genotype.case_kind` es la
  constante `snp` y `f5_locus.ref/alt` son bases (`schemas/truth-v2.schema.json`).
- **F5 coloca sitios y luego sortea su ALT.** Las posiciones salen de los gaps
  del perfil (`src/f5_layout.rs:1728-1790`) y el ALT, de un ámbito propio por
  sitio (`src/f5_layout.rs:1916`). F5 no publica VCF ni lecturas.
- **La interfaz escribe el documento general v4 solo si hay loci de
  repetición; si no, v2** (`ui/config_editor.py:2438`). El perfil de F5 no se
  edita en Marimo: se enlaza por ruta y huella (`tests_ui/test_f5_6_1_scenario_ui.py`).
- **Huellas congeladas que ya existen:** el recorrido técnico sin
  amplificación (`tests/pcr_frozen_identity.rs`), ocho huellas de F5 en
  `tests/f5_publication.rs`, y un vector de dos SHA-256 en `tests/f5_cli_run.rs:242`.
  No he comprobado si alguna congela una matriz **generada** ni un escenario
  con microsatélites generados.

## 3. Decisiones que el responsable no ha fijado

> **Decisiones del responsable, 28 de octubre de 2026.** Cerradas sobre las
> recomendaciones de este apartado:
>
> - **DP-IND-01 — aceptada.** Un SNV dentro de un *indel* es el mismo locus:
>   las filas solapadas forman un locus compuesto y cada copia lleva uno solo
>   de sus alelos.
> - **DP-IND-03 — aceptada.** En una misma secuencia, los *indels* escritos a
>   mano y los de la matriz son incompatibles: o unos u otros. Los SNP escritos
>   a mano se permiten fuera de los locus con *indel*.
> - **DP-IND-07 — aceptada, con su razón biológica.** El generador no da
>   preferencia a homopolímeros ni repeticiones porque **esa variación
>   pertenece al modelo de microsatélites** (alelos por deslizamiento,
>   DTC-094 a 097), no a los *indels* de la matriz; por eso DP-IND-13 rechaza
>   un *indel* que toque un tracto declarado. Limitación que se declara: solo
>   varían los microsatélites declarados; mientras no existan el
>   descubrimiento de microsatélites de un genoma suministrado (F6.4.a) ni su
>   generación al azar (F6.4.6), la variación en repeticiones de una
>   simulación queda por debajo de la real.
> - **DP-IND-11 — aceptada, con aviso obligatorio.** En F5 no hay solapes
>   dentro de un bloque; los solapes solo existen por la vía de la matriz. La
>   limitación —ningún SNP dentro de una deleción en la población fundadora, y
>   deleciones acortadas o convertidas en SNV cuando los sitios están muy
>   juntos— **debe avisarse**: en la documentación, en la configuración
>   efectiva y en el manifiesto de toda ejecución F5 con *indels*, junto con
>   el recuento de sitios que quedaron como SNV por no caber la deleción.
>
> - **DP-IND-13 — sustituida.** Un *indel* de la matriz que toque un
>   microsatélite **no se rechaza**: forma un solo locus con él y pasa a ser un
>   alelo más de ese locus; cada copia lleva uno solo de sus alelos (una
>   longitud del microsatélite o el *indel*). Es la misma regla que DP-IND-01.
>   El generador de la matriz **evita** los microsatélites, porque su
>   variación de longitud ya la produce el deslizamiento y se contaría dos
>   veces. Los pasos que dependían de la recomendación anterior deben
>   reescribirse.
> - **Orden de ejecución.** Este plan se ejecuta **después** del plan de
>   detección de microsatélites (`plan-implementacion-deteccion-microsatelites-0.1.md`):
>   los microsatélites de una simulación serán siempre los detectados sobre la
>   referencia final, y la regla anterior necesita ese catálogo.
>
> - **TAREA PENDIENTE heredada del plan de detección (28/10/2026), que no
>   puede perderse:** al construir los locus compuestos, **un SNV dentro de un
>   microsatélite forma también un solo locus con él**, como el *indel*. Esto
>   sustituye la regla provisional del plan de detección (DP-SSR-11), que
>   retira el tracto y publica la causa.
> - **DP-IND-09 y DP-IND-14 — simplificadas.** La etiqueta `v0.2.0-test` es
>   una versión alfa, no una publicación (decisión del responsable del
>   28/10/2026): los esquemas se retocan en su versión actual, sin versiones
>   nuevas. La huella contra el commit base queda solo como comprobación de
>   cambios involuntarios, no como obligación de compatibilidad.
>
> Siguen pendientes las demás DP-IND, sobre las que el agente raíz recomienda
> aceptar la recomendación de cada una.
>
> **Nota del 28/10/2026 sobre números y versiones, escrita al aplicar la
> entrega 0 del plan de detección.** (1) El plan de detección **ya tomó**
> DTC-103 y DTC-104 en `decisiones-tecnicas-0.1.md`, y REQ-REF-38 a 40,
> REQ-CONF-142/143, REQ-MAN-98, REQ-VAL-64, REQ-UI-17, ACE-121 y DTF-064 a
> DTF-067 en sus documentos: los DTC de este plan pasan a **DTC-105 a
> DTC-107** (no hay conflicto con los REQ-IND-69 a 71). (2) Por la regla de
> la etapa alfa (nota de DTC-036), el plan de detección **no reserva** ninguna
> versión de esquema: el general v5 y la efectiva y el manifiesto v7 que este
> plan propone quedan libres, y por la misma regla de la etapa
> alfa este plan puede retocar los esquemas en su versión actual. (3) La tarea obligatoria de los loci compuestos con
> SNV está registrada como DTF-066 y asignada a este plan. (4) Este plan
> hereda del de detección la exclusión de los tractos detectados en su
> generador y la consulta por intervalo del catálogo (REQ-CONF-143).

Cada una con opciones, consecuencias y una recomendación. Se identifican como
DP-IND-nn (decisión pendiente). **Ninguna se da por tomada en los pasos**: los
pasos están escritos contra la recomendación y se marcan con «(según
DP-IND-nn)» donde dependen de ella.

### DP-IND-01. Cómo se reparten los alelos cuando se solapan

**El problema.** Una deleción de `ACG` anclada en la base 100 y un SNV en la
101 no pueden estar en la misma copia: la base 101 ya no existe. En una
población real eso se resuelve porque los dos alelos van en haplotipos
distintos.

| Opción | Qué hace | Consecuencias |
|---|---|---|
| **A1. Locus compuesto por solape del `REF` con ancla** | Dos filas se solapan si sus intervalos `REF` —el del VCF, con la base de anclaje— comparten alguna base. Las filas que se solapan, de forma transitiva, forman **un locus compuesto multialélico**: cada copia lleva **uno solo** de sus alelos o la referencia | Reutiliza el sorteo multialélico que ya existe (HWE para P = 1–10 y `F` diploide con REQ-EVO-48). La frecuencia pedida es exactamente la esperada; HWE y `F` son exactos por construcción. Exige que la suma de alelos del locus no pase de 1 por población, la misma regla que hoy para un locus multialélico (REQ-IND-60). Coincide con la regla con que el VCF funde registros (commit `6bef50e`) y con la de DTC-064 de que ningún `REF` pise a otro. Un SNV en la base de anclaje de una deleción **no** puede ir en la misma copia que la deleción, aunque biológicamente podría |
| A2. Locus compuesto por solape de las bases **cambiadas** | Como A1, pero el ancla no cuenta: un SNV en la base 100 y la deleción anclada en 100 son independientes | Más realista en el ancla; las dos reglas —la de la matriz y la del VCF— dejan de coincidir y hay que explicar por qué. Casos de borde más finos (inserción en la junta de una deleción) |
| B. Sorteo independiente por fila y resolución por copia | Cada fila se sortea por su cuenta; si una copia recibe dos alelos incompatibles, uno gana (como DTC-080 hace con las deleciones directas) | La frecuencia realizada del perdedor queda por debajo de la pedida de forma sistemática, no por azar; HWE y `F` dejan de cumplirse en los loci afectados; dos inserciones en el mismo punto no tienen orden definido |

**Recomendación: A1.** Es la única que no fabrica copias imposibles ni
distorsiona las frecuencias, reutiliza lo que hay y se enuncia en una frase
que ya rige en el VCF. Un locus compuesto **es** un locus multialélico; un
locus formado por un solo SNV queda exactamente como hoy, y esa es la
condición de los bytes idénticos (sección 5).

**Lo que hay que fijar con A1**, y que la recomendación propone:

- el **locus compuesto** se identifica por la primera base de su tramo; su
  `REF` es la referencia del tramo completo, y cada alelo se expresa sobre ese
  tramo para el ámbito de azar `initial_genotype` y para `initial_phase_locus`;
- el **orden de los alelos** es el orden canónico de las filas (posición,
  `ref`, `alt`), que para un locus de SNV coincide con el orden A/C/G/T de hoy;
- las **frecuencias**: la suma por población de los alelos del locus se valida
  y normaliza con DTC-049 igual que hoy; la diferencia es que el grupo pasa a
  ser el locus compuesto en vez de `(población, contig, posición, ref)`;
- un **límite de tamaño**, ver DP-IND-12.

### DP-IND-02. Cómo se aplican a las copias los indels de la matriz

| Opción | Consecuencias |
|---|---|
| **Como eventos por copia, igual que los SNV de la matriz, aplicados antes que los eventos directos y en orden decreciente de posición de referencia** | En orden decreciente, cada indel encuentra sus coordenadas vigentes iguales a las de la referencia, porque nada a su izquierda ha cambiado todavía de longitud. No hace falta una segunda aritmética de coordenadas, que S-32 descartó expresamente. La verdad ya sabe publicar inserciones y deleciones por evento (`schemas/truth-v4.md`) |
| Como un evento compartido por todas las copias portadoras del individuo | Un homocigoto tendría un caso con dosis 2 en vez de dos casos con dosis 1. El alelo estructural de la verdad enlaza `piece_links` de «la copia portadora», en singular; no he comprobado que admita varias portadoras de una inserción |

**Recomendación: eventos por copia, en orden decreciente.** Los SNV de la
matriz siguen aplicándose como hoy, en el mismo orden, para no mover bytes; los
indels de la matriz van después de ellos y antes de los eventos directos.
Consecuencia en la verdad: un individuo homocigoto para una deleción de la
matriz publica **dos casos *indel***, uno por copia, cada uno con su evento.
Es lo que ya ocurre con dos deleciones directas iguales en dos fases.

**Lo que hay que comprobar antes de fijarlo** (paso F2.2.14): que el contrato
`truth.jsonl` v4 acepta esos eventos sin cambiar de esquema —solo el juez del
conjunto de casos cambia, porque hoy exige un caso SNP por cada fila y una
fila *indel* no lo tiene—. Si la comprobación dice lo contrario, se para y se
informa: la verdad necesitaría versión nueva y la entrega 2 crece.

### DP-IND-03. Indels de la matriz y eventos directos en la misma secuencia

**El problema, que ya costó una ronda entera.** Las inserciones y deleciones
directas se declaran en coordenadas **vigentes** de su copia. Si los indels
de la matriz se aplican antes, esas coordenadas pasan a depender de lo que el
sorteo puso en cada copia, que el usuario no conoce al escribir el TOML.

| Opción | Consecuencias |
|---|---|
| **Rechazar, antes de sortear, las inserciones y deleciones directas en una secuencia donde la matriz activa tenga alguna fila *indel*; y rechazar un SNP directo dentro del `REF` de un locus compuesto que contenga un *indel*** | Determinista: se decide con lo que el usuario escribió y la matriz, sin semilla. El significado de los TOML actuales no cambia. Se pierde la combinación de las dos fuentes en la misma secuencia |
| Aplicar los indels de la matriz **después** de los directos, traduciendo su coordenada de referencia a la vigente de cada copia | Mantiene la combinación, pero introduce la traducción referencia → vigente y una política para cada colisión (deleción directa que borra el locus, inserción directa dentro de una deleción de la matriz) que depende del sorteo |
| Redefinir las coordenadas de los eventos directos como «vigentes sin contar la matriz» | Segunda aritmética de coordenadas; la descartó S-32 |

**Recomendación: la primera.** Es la única sin rechazos dependientes de la
semilla y sin aritmética nueva. Los SNP directos siguen permitidos fuera de
esos loci porque se anclan a la referencia. DTC-080 (SNV de la matriz borrado
por una deleción directa) sigue rigiendo en las secuencias sin indels de
matriz.

### DP-IND-04. Representación y normalización de una fila *indel*

**Recomendación**, que no tiene alternativa razonable salvo en el segundo
punto:

- **Forma VCF anclada.** `position` es la base de anclaje, en base 1 y en
  coordenadas de la **referencia** de la matriz. Deleción: `ref` = ancla +
  bases borradas, `alt` = ancla. Inserción: `ref` = ancla, `alt` = ancla +
  bases nuevas. Cualquier otra forma (complejo, MNP, ancla a la derecha) se
  rechaza.
- **Matriz suministrada: se exige ya normalizada** —alineada a la izquierda y
  mínima—, y una fila que no lo esté se rechaza con un diagnóstico que
  **muestra su forma normalizada**. La alternativa, normalizar al preparar,
  podría fundir en silencio dos filas distintas de la entrada en la misma
  variante. Las herramientas habituales (`bcftools norm -m-`) ya producen esa
  forma.
- **Criterio de normalización, el mismo en todo el proyecto:** una deleción
  de `D` tras el ancla `a` está alineada si `a ≠ último(D)`; una inserción de
  `I` tras el ancla `a`, si `último(I) ≠ a`.
- **Bases.** `ref` y `alt` solo A/C/G/T; ni `N` ni otra letra, ni en la
  suministrada ni en la generada. Es más estricto que REQ-IND-58, que tolera
  `N` dentro de una deleción, y es deliberado: el ancla y la comparación de
  normalización necesitan bases concretas.
- **Bordes**, los de REQ-IND-40/41: la inserción necesita ancla y no puede ir
  tras la última base; la deleción deja al menos una base a cada lado.
- **Longitud**: secuencia nueva entre 1 y umbral − 1. Ver DP-IND-05.
- **Coordenadas.** Toda fila está en coordenadas de la referencia sobre la que
  se preparó la matriz. Las coordenadas vigentes solo aparecen al aplicar el
  evento a una copia, y por DP-IND-02 coinciden con las de referencia en ese
  instante.

### DP-IND-05. Relación con el umbral configurable

| Opción | Consecuencias |
|---|---|
| **Cada ejecución valida las filas contra su umbral efectivo; la matriz no guarda umbral** | Hoy el umbral es la constante 50. Cuando F6.2 lo haga configurable, una matriz reutilizada con un umbral menor se **rechaza** si tiene una fila igual o mayor —no se reclasifica como estructural—. No hay campo nuevo |
| Registrar en el JSON de la matriz el umbral con que se preparó y exigir igualdad | Una matriz válida dejaría de reutilizarse solo porque cambió un número que no la afecta |
| Fijar 49 como tope duro de la matriz | Contradice REQ-CONF-55, que pide coherencia con el umbral efectivo (REQ-IND-31) |

**Recomendación: la primera**, leyendo el umbral del mismo sitio que ya usa
`SmallIndelLengths` para no crear otra fuente. El parámetro de longitud
máxima del generador (DP-IND-06) se valida contra ese mismo umbral.

**Observación, no comprobada a fondo:** no encontré en `src/application.rs` ni
en `src/small_event_composition.rs` una comprobación del umbral para las
inserciones y deleciones **directas**. Si no la hay, REQ-IND-31 ya no se
cumple hoy para ellas. No es objeto de este plan; se anota como sospecha para
el registro de deudas.

### DP-IND-06. Parámetros del generador

**Principio:** lo mínimo útil, y que una configuración sin la tabla nueva
genere exactamente lo de hoy.

Propuesta de tabla opcional `[population_snvs.generate.indels]`:

| Parámetro | Propuesta | Por omisión | Fuente o justificación |
|---|---|---|---|
| Cantidad | El mismo campo que la matriz usa según su `count_mode`: `locus_count` (exacto) o `mean_loci` / `mean_density_per_base` (Poisson). Se hereda `count_mode` del padre | **Sin valor por omisión**: si la tabla existe, el usuario declara la cantidad, como para los SNV | Referencia documental, no valor por omisión: 1000 Genomas fase 3 publicó 84,7 millones de SNP y 3,6 millones de indels cortos ([Nature 2015](https://doi.org/10.1038/nature15393); comprobado en el resumen), unos 4 indels por cada 100 SNV, **infraestimado** porque la sensibilidad para indels de ese proyecto era menor |
| Distribución espacial | Se hereda la del padre (`uniform` o `weighted_regions`, con sus regiones y pesos) | — | Una sola declaración espacial; separar las dos es ampliación posterior |
| `insertion_fraction` | Probabilidad de que un indel sea inserción | 0,5 | **Sin preferencia, y es una elección, no una medida.** Ver la nota sobre fuentes, abajo |
| Longitud | Geométrica truncada en `[1, max_length]` con parámetro `length_geometric_p` | `p = 0,5` | **Forma convencional, no estimación.** Reproduce el rasgo cualitativo pedido —la masa se concentra en 1–2 pb (P(1)=0,5, P(2)=0,25)—. No he podido verificar una cifra publicada de la proporción de indels de 1 pb |
| `max_length` | Longitud máxima de secuencia nueva | umbral − 1 (49) | REQ-CONF-55 |
| Contenido de la inserción | A/C/G/T uniformes por base; la última, entre las tres distintas del ancla (DP-IND-07) | — | REQ-IND-59 ya exige solo A/C/G/T en inserciones aleatorias |
| Espectro de frecuencias | **El mismo de cada población que usan los SNV** | — | Separarlo no aporta nada que el usuario haya pedido |

**Nota sobre fuentes.** Montgomery y col. (2013) analizan 1,6 millones de
indels de 179 genomas humanos y afirman, en el resumen, que entre el 43 % y el
48 % caen en el 4,03 % del genoma y que el deslizamiento de la polimerasa
explica más de tres cuartas partes ([Genome Research 23:749](https://doi.org/10.1101/gr.148718.112)).
Las cifras de proporción de 1 pb y de inserciones frente a deleciones están en
el texto completo, que **no he leído**. Antes de aprobar el DTC hay que leerlo
y decidir si sustituye los dos valores por omisión anteriores. Mientras tanto,
el plan los marca como convencionales, igual que el plan de microsatélites
hizo con su fracción de duplicados (§2.3, «Por qué no se fija un valor por
defecto»).

**Recomendación:** la tabla anterior. Alternativa más estricta: sin valores
por omisión para `insertion_fraction` ni `length_geometric_p`, obligando a
declararlos. Es más honesta y menos cómoda; la decisión es del responsable.

### DP-IND-07. Normalización en el generador y preferencia por homopolímeros

| Opción | Consecuencias |
|---|---|
| **Normalizado por construcción en su ancla.** Se eligen anclas sin reemplazo; la inserción sortea su última base entre las tres distintas del ancla; la deleción sortea su longitud **condicionada** a `ref[a] ≠ ref[a+k]`, renormalizando los pesos de las longitudes posibles. Si ninguna longitud lo cumple, el indel es inserción | Cada fila sale ya alineada; como hay un indel por ancla y la forma alineada es única, **no puede haber dos filas iguales** y el modo exacto sigue siendo exacto. Sin desplazamiento ni reintentos. **Sin preferencia por homopolímeros**: cada tramo repetido aporta un ancla, no una por base. Sesgo leve en la última base de la inserción, declarado |
| Sortear libremente y alinear a la izquierda después | Los homopolímeros largos atraen indels en proporción a su longitud, que se parece más a lo real. Pero dos sorteos pueden acabar en la misma variante: hay que fundirlos (el recuento exacto deja de serlo) o reintentar (REQ-IND-56 lo prohíbe) |
| Preferencia explícita por repeticiones | Modelo nuevo de mutación por contexto; sobreingeniería para este encargo |

**Recomendación: la primera.** La infrarrepresentación de los indels en
repeticiones queda **declarada** en la configuración efectiva como
limitación, con la cita de Montgomery, y se abre una ampliación futura (DTF)
para la preferencia por contexto.

### DP-IND-08. Solapes en la matriz generada

DTF-039 aplaza la «generación aleatoria de loci multialélicos» porque exige
un espectro conjunto cuyas frecuencias no sumen más de 1. Si el generador
permite que un indel caiga sobre un SNV, crea justo eso.

| Opción | Consecuencias |
|---|---|
| **Permitir solapes.** SNV e indels se eligen con selecciones independientes; cuando caen solapados, el locus compuesto reparte sus frecuencias **por el resto**: recorriendo sus alelos en orden canónico, `af_efectiva_i = af_sorteada_i × (1 − Σ af_efectiva_j, j<i)` | Los solapes aparecen con la frecuencia natural de dos densidades independientes, como pidió el responsable. La suma nunca pasa de 1, no se reintenta ni falla. El primer alelo no se toca; los demás bajan en proporción a lo que ya ocupan los anteriores, que con espectros dominados por frecuencias bajas es casi nada. Se publican la frecuencia sorteada y la efectiva. **Reabre DTF-039 solo para este caso** |
| Excluir los solapes en la generación | DTF-039 intacta; los solapes solo llegan por matriz suministrada. Exige excluir del sorteo de indels los tramos tocados por SNV, que es más código que permitirlos |
| Permitir solapes y fallar si una suma pasa de 1 | Rechazo dependiente de la semilla |

**Recomendación: permitir solapes con reparto por el resto.** Un locus de SNV
solo no se toca, de modo que una matriz generada sin indels da los mismos
bytes. Hay que anotar en DTF-039 que esta parte sale de ella.

### DP-IND-09. Versiones que suben, y qué las decide

| Opción | Consecuencias |
|---|---|
| **Versiones nuevas, promovidas por contenido.** JSON de la matriz v2 si tiene alguna fila *indel*; configuración efectiva v7 y manifiesto general v7 si la matriz activa la tiene; documento general v5 solo para declarar `generate.indels` | Una ejecución sin indels no cambia nada. Es el criterio de la verdad v4/v5, del manifiesto v4/v6 y del catálogo de F5. Coste: otra rama en la escalera de versiones y su juez |
| Excepción expresa, como DT-041: relajar en las v6 publicadas los enumerados A/C/G/T a cadenas A/C/G/T | No mueve bytes y es menos código. Pero cambia el significado de un contrato publicado («los alelos son bases») sin que su versión lo diga, y DT-041 se concedió porque la alternativa movía bytes; aquí no los mueve |

**Recomendación: versiones nuevas, promovidas por contenido.** Detalle:

- **efectiva v7** es superconjunto de v6 (cadenas en `ref`, `allele` y en la
  auditoría de normalización; tramo del locus compuesto; bloque de auditoría
  del reparto por el resto; limitaciones declaradas del generador). Se publica
  con un documento general v2, v4 o v5 siempre que la matriz activa tenga un
  *indel*;
- **manifiesto general v7**, superconjunto de v6, con la misma regla;
- **general v5** = v4 + `[population_snvs.generate.indels]`. Un documento v2 o
  v4 que la declare se rechaza;
- la **verdad** no cambia de versión si se confirma DP-IND-02;
- F5: ver DP-IND-10.

Una consecuencia que hay que aceptar: con la matriz **suministrada**, la
versión publicada depende del contenido de un fichero y no del TOML. Ya pasa
con los loci de repetición.

### DP-IND-10. F5: dónde se declaran los indels y qué se versiona

| Opción | Consecuencias |
|---|---|
| **En el perfil**, metadata v3 con un objeto cerrado `indel_sites` (fracción de sitios *indel*, `insertion_fraction`, `length_geometric_p`, `max_length`), obligatorio en v3; `haplotype_count_distribution` opcional en v3 | El perfil ya contiene la forma estadística (bloques, gaps, espectro, H). No necesita control en Marimo porque el perfil no se edita allí. Rompe la regla de «la versión fija el conjunto cerrado de claves» para H en v3; la protección que daba —que un perfil antiguo no adquiera una opción por descuido— se mantiene, porque v1 y v2 no cambian |
| En el escenario F5 (v2) | Exige control en Marimo, paridad y un escenario v2 |

**Recomendación: en el perfil v3.** Lo que se versiona, solo si el perfil
declara `indel_sites`: catálogo de azar `f5-rng-catalog-v3` (superconjunto de
v2), contrato de sitios `f5-02-sites-v2` (`ref`/`alt` como cadenas), verdad F5
nueva —la siguiente libre es `"6"`, igual que v2 con un caso *indel*—,
configuración efectiva y manifiesto de F5 en versiones nuevas superconjunto de
v5 (numeración a fijar en el DTC: la 6 y la 7 ya las usa el recorrido
general), y el algoritmo de admisión de recursos v9, porque las cotas de bytes
publicados de v8 suponen alelos de una base. Sin `indel_sites`, todo como hoy.

### DP-IND-11. F5: cómo entran los indels en los bloques

**Recomendación**, con sus piezas:

- **Sin solapes dentro de un bloque.** La matriz de haplotipos de F5 es
  binaria sitio a sitio; dos sitios solapados exigirían prohibir haplotipos
  con los dos alelos, y eso cambia la búsqueda de DTC-089. Con este modelo los
  solapes viven en la matriz (decisión 2) y no en F5. **Hay que confirmarlo**:
  la decisión 2 del responsable habla de la matriz, y la recomendación lo
  interpreta así.
- **Las posiciones no cambian.** El layout y los gaps se sortean exactamente
  como hoy; después, cada sitio decide su clase con un sorteo propio
  (`Bernoulli(fracción)`), en un ámbito nuevo, por sitio.
- **Deleción que cabe.** Su longitud se sortea truncada a lo que cabe sin tocar
  el `REF` del sitio siguiente ni salirse del bloque (`k ≤ s' − s − 1`; en el
  último sitio, `k ≤ fin − s`) y condicionada a la normalización de
  DP-IND-07. Si no cabe ninguna longitud, el sitio es SNV y se cuenta. Con un
  perfil de gaps de 1 no habrá deleciones, y se publica.
- **Inserción**: ancla en el sitio, contenido como en DP-IND-07.
- **Marginales**: el mismo `snp_spectrum.tsv` para los dos tipos de sitio.
- **`r²(d)`**: sin cambio; la distancia es entre anclas y un sitio *indel* es
  bialélico como un SNP.
- **Oráculo**: reimplementa la clase, la longitud y el contenido; una
  transcripción, con la misma salvedad de independencia que DTC-089 ya
  declara para H.

### DP-IND-12. Cota del locus compuesto

**El riesgo.** Las estadísticas por locus enumeran todos los genotipos
esperados (`expected_genotypes`, `src/initial_individual_materialization.rs`)
y los publican en la configuración efectiva. Con `n` alelos y ploidía `P` son
`C(n+P−1, P)`: con 4 alelos y P = 10, 286; con 20 alelos y P = 10, unos 20
millones de filas en un JSON. Hoy no pasa porque un locus de SNV tiene como
mucho 4 alelos; un locus compuesto encadenado no tiene ese techo.

**Recomendación:** rechazar antes de sortear un locus compuesto cuya
enumeración supere una cota fija de genotipos, con un diagnóstico que nombre
el locus y la cuenta. El valor de la cota se fija en el DTC; el índice de
alelo es `u8`, así que 255 alelos es además un techo duro.

### DP-IND-13. Microsatélites: indel dentro o junto a un tracto declarado

**Hoy:** una fila de la matriz se rechaza si su **posición** cae dentro del
tracto (`src/execution.rs:2772-2789`).

**Recomendación:** para una fila *indel*, se rechaza si su intervalo `REF`
**toca** el tracto o la base inmediatamente anterior, que es la base en que se
ancla una inserción alineada en el borde izquierdo de una repetición. Esa
inserción no entra en el intervalo del tracto, pero cambia su número aparente
de unidades, y la anotación del alelo de repetición sería falsa, que es lo que
S-32 prohíbe. **Las filas SNV conservan la regla de hoy**, para no cambiar el
comportamiento de ningún escenario existente. El generador de indels excluye
de sus candidatos esos mismos tramos, para que el rechazo nunca dependa de la
semilla.

### DP-IND-14. Catálogo de ámbitos de azar

Los ámbitos nuevos del generador (recuento, selección de anclas, clase y
longitud, contenido de inserción) van en `subseed-scopes-v11`, superconjunto
de v10. **Pero la versión del catálogo sí se publica**: `slippage_json` escribe
`rng_scope_catalog_version` = `SUBSEED_CATALOG_VERSION`
(`src/results.rs:970`, `src/repeat_slippage.rs:100`), y el esquema efectivo v6
lo exige (`schemas/effective-config-v6.schema.json:3055`). Subir la constante a
`"11"` cambiaría los bytes de todo escenario con microsatélites generados.

**Recomendación:** fijar la constante del deslizamiento al literal `"10"` —la
versión en que nació su ámbito— en el mismo paso que se añade v11, con una
huella congelada que lo demuestre. Ver la contradicción C-1.

## 4. Cambios previos en los documentos normativos

Redactados como **propuesta**. Se aprueban antes de la entrega 1. El texto
definitivo lo escribe el paso F2.2.13.

### 4.1. Especificación (`especificaciones-0.1.md`)

| Elemento | Hoy | Propuesta |
|---|---|---|
| **ACE-10** | «Los SNV heredables iniciales proceden de una matriz…; *Indels*, variantes estructurales y mutaciones *de novo* pueden indicar cantidad media o densidad» | «Los SNV **y los indels pequeños** heredables iniciales **compartidos por la población** proceden de una matriz poblacional suministrada o generada; su aleatoriedad está en las frecuencias y genotipos. Los indels **privados de un individuo**, las variantes estructurales y las mutaciones *de novo* pueden indicar cantidad media o densidad…» (resto igual) |
| **REQ-IND-04, REQ-IND-09** | Matriz solo para SNV | Añadir «y los indels pequeños poblacionales»; la cantidad media o densidad sigue para los indels privados |
| **REQ-IND-60** | `ALT` A/C/G/T distinto de `REF` A/C/G/T | Una fila es un SNV o un indel pequeño anclado y normalizado (DP-IND-04); secuencia nueva entre 1 y umbral − 1; A/C/G/T; **las filas cuyos `REF` comparten una base forman un locus compuesto** y su suma ALT por población cumple la regla actual |
| **REQ-IND-61** | «F2 genera solo loci bialélicos» | Genera SNV bialélicos y, opcionalmente, indels; los solapes entre clases forman loci compuestos con reparto por el resto (DP-IND-08) |
| **REQ-IND-62** | Muestreo por locus | Precisar que el locus es el compuesto y que cada copia recibe un único alelo suyo |
| **REQ-IND-63** | Subsemillas listadas | Añadir las del generador de indels |
| **REQ-IND-64** | Exclusión matriz SNV / SNP por cantidad | Añadir que los indels de matriz y los privados por individuo de F6.5.a se decidirán al planificar F6.5.a |
| **REQ-IND-65** | Estadísticas por locus | Por locus compuesto, con alelos como cadenas |
| **Nuevo REQ-IND-69** | — | Restricción con eventos directos (DP-IND-03) |
| **Nuevo REQ-IND-70** | — | Microsatélites (DP-IND-13) |
| **Nuevo REQ-IND-71** | — | Cota del locus compuesto (DP-IND-12) |
| **REQ-CONF-128/129** | «matriz de frecuencias SNV» | «matriz de variantes pequeñas»; parámetros de `generate.indels` |
| **REQ-VAL-55, REQ-MAN-91, REQ-UI-14, ACE-116** | SNV | Añadir los indels, el locus compuesto, la normalización y el reparto por el resto |
| **Texto de §4.3.1** sobre el TSV, el JSON v1 y `[population_snvs]` | Gramática SNV | Gramática de la fila *indel*, JSON v2 condicional, tabla `generate.indels` |
| **REQ-EVO-42** | «Los bloques F5 contienen exclusivamente SNP bialélicos» | «…SNP bialélicos y, si el perfil lo declara, indels pequeños bialélicos no solapados…» |
| **REQ-EVO-74** | «F5 no admite … indels» | Quitar «indels»; «cada bloque contiene al menos dos sitios bialélicos» |
| **REQ-EVO-75** | «marginales SNP» | «marginales por sitio» |
| **REQ-EVO-77** | Siete TSV, verdad v2 | Precisar las versiones condicionales de DP-IND-10 |

El nombre `population_snvs` (tabla, ficheros, identificadores de artefacto) se
**conserva**: cambiarlo movería bytes de toda ejecución. El texto normativo
dirá una vez que «la matriz SNV» contiene también indels pequeños.

### 4.2. Decisiones técnicas (`decisiones-tecnicas-0.1.md`)

- **DTC-103 (nueva). Indels en la matriz poblacional:** representación,
  normalización, locus compuesto, reparto de genotipos, aplicación por copia en
  orden decreciente, restricción con eventos directos, cota, microsatélites,
  versiones condicionales y garantía de bytes. Enmienda DTC-052 y DTC-054.
- **DTC-104 (nueva). Generador de indels de la matriz:** parámetros, valores
  por omisión con su fuente o su etiqueta de convencionales, normalización por
  construcción, reparto por el resto, ámbitos de azar, limitación declarada.
  Enmienda DTC-071 (filas v11) y DTC-072 (consumidores).
- **DTC-105 (nueva). Indels en F5:** perfil v3, clase por sitio, longitud que
  cabe, contratos versionados, admisión v9, oráculo. Enmienda DTC-089, que hoy
  dice que F5 «rechaza … *indels*».
- **DTC-061**: el conjunto de casos de la verdad (filas *indel* sin caso SNP).
- **DTC-064/065**: ninguna regla cambia; se añade el corpus de indels de matriz
  a la puerta externa.
- **DTC-101**: la batería de parejas finales gana escenarios (paso F2.2.21).
- **Nota en §14.1 del plan de microsatélites** (contradicción C-1).

### 4.3. Otros documentos

- **Hoja de ruta:** nota de «reapertura acotada» en F2 y en F5, con la forma
  de las de F1 y F4, apuntando a la sección de registro de este plan.
- **Registro de ampliaciones futuras:** DTF-039 (sale el caso del solape),
  DTF-040 (sale la parte de indels de F5), DTF nueva para la preferencia por
  contexto repetido.
- **Trazabilidad:** filas de los REQ nuevos y modificados.
- **Registro de deudas:** C-2 y la sospecha de DP-IND-05.
- `schemas/truth-v4.md`, `schemas/f5-tsv-v1.md`: secciones de lo que cambia.

## 5. Garantía de compatibilidad

**Promesa:** una ejecución que no use indels de matriz —ni filas *indel* en la
matriz activa, ni `generate.indels`, ni `indel_sites` en el perfil de F5—
publica **exactamente los mismos bytes** que el commit base `d101f8f`.

**Cómo se ancla, y por qué no basta comparar dos ejecuciones del binario
nuevo:** esa comparación es circular. Se congelan huellas **antes de tocar el
código** (paso F2.2.13) con el binario de `d101f8f`, compilado en un árbol
extraído con `git archive`, sobre estos escenarios:

1. matriz suministrada con un locus bialélico, uno multialélico y uno todo-REF,
   con VCF activado;
2. matriz generada con espectro y regiones ponderadas, modo `exact`;
3. el mismo, modo `stochastic`, con `F` diploide;
4. un escenario con un microsatélite **generado** (el literal del catálogo,
   DP-IND-14);
5. F5 con perfil v1 y F5 con perfil v2 (las ocho huellas de
   `tests/f5_publication.rs` ya existen; se comprueba que cubren los dos);
6. el escenario técnico de `tests/pcr_frozen_identity.rs`, que ya existe.

Cada huella cubre **todos** los ficheros publicados: TSV y JSON de la matriz,
configuración efectiva, manifiesto, verdad, VCF e índice. La prueba se
verifica fallando: se muta el orden de las filas de la matriz en una copia de
trabajo y la huella del escenario 1 tiene que morir.

## 6. Entregas y pasos

Orden: primero los contratos y la línea base; después la ruta suministrada
—que no necesita azar nuevo y permite probar el reparto con verdad conocida—;
después el generador; F5 al final porque no depende del resto y es lo más
grande. Cada entrega termina en pausa.

### Entrega 0 — Documentos y línea base

- ☐ **F2.2.13.** Documentos normativos aprobados y huellas de la línea base
  congeladas.
  - Depende de: cierre de las decisiones de la sección 3.
  - Procede de: decisiones 1–5 del responsable; sección 4.
  - Cambia: `design/especificaciones-0.1.md`, `design/decisiones-tecnicas-0.1.md`,
    `design/hoja-de-ruta-0.1.md`, `design/registro-ampliaciones-futuras-0.1.md`,
    `design/trazabilidad-requisitos-0.1.md`, `design/registro-deudas-tecnicas-0.1.md`;
    pruebas nuevas de huella congelada (sección 5) con su variante defectuosa.
    Incluye la lectura del texto completo de Montgomery 2013 que pide DP-IND-06.
  - Se comprueba con: `cargo test --workspace --all-targets --locked` en verde
    con las huellas nuevas; la mutación del orden de filas las hace fallar
    (SHA del fichero mutado antes y después anotados).

> Revisión de la entrega 0: el responsable aprueba el texto normativo. Sin
> esto no se escribe código.

### Entrega 1 — La matriz admite indels (suministrada y reutilizada)

- ☐ **F2.2.14.** Comprobación previa: la verdad v4 representa los indels de
  matriz sin cambiar de esquema.
  - Depende de: F2.2.13.
  - Procede de: DP-IND-02.
  - Cambia: nada del producto; una prueba de exploración que construye a mano
    una copia con un evento *indel* de identificador derivado de la matriz y la
    pasa por el lector y el juez de la verdad.
  - Se comprueba con: la prueba pasa con el esquema actual. **Si no pasa, se
    para y se informa**: DP-IND-02 y la entrega 2 se replantean.
- ☐ **F2.2.15.** Analizador y preparación de filas *indel*.
  - Depende de: F2.2.13.
  - Procede de: REQ-IND-60 (enmendado), DTC-103, DP-IND-04/05.
  - Cambia: `src/population_snvs.rs` (gramática de alelo, normalización con
    diagnóstico que propone la forma correcta, bordes, umbral, `N`,
    identificador desde cadenas, orden canónico por bytes).
  - Se comprueba con: pruebas unitarias de cada rechazo; **una matriz solo SNV
    produce el mismo TSV y el mismo identificador de variante que antes**
    (huella del escenario 1).
- ☐ **F2.2.16.** Locus compuesto, suma y normalización de frecuencias.
  - Depende de: F2.2.15.
  - Procede de: DP-IND-01 (A1), DP-IND-12, DTC-049.
  - Cambia: `src/population_snvs.rs` (agrupación transitiva por `REF` con
    ancla, auditoría de normalización por locus compuesto, cota del locus).
  - Se comprueba con: fixture **saturado** que contiene a la vez: SNV dentro de
    una deleción, inserción en el ancla de la deleción, una cadena de tres
    filas que solo se unen por transitividad, y dos filas **adyacentes sin
    base común**, que tienen que quedar separadas. Un agrupador sin
    transitividad, o que cuente la adyacencia como solape, falla en un caso
    distinto cada uno.
- ☐ **F2.2.17.** JSON de la matriz v2 condicional, reutilización y su juez.
  - Depende de: F2.2.16.
  - Procede de: DP-IND-09.
  - Cambia: `src/population_snvs.rs` (metadata v2 solo con *indels*;
    `validate_reusable_population_snvs*` acepta v1 y v2 y exige la coherencia
    versión↔contenido en las dos direcciones), esquema o descripción del JSON
    v2, `scripts/guard-coverage-baseline.tsv` (fila de
    `src/population_snvs.rs`, las dos columnas).
  - Se comprueba con: un JSON v1 con filas *indel* y un v2 sin ellas se
    rechazan; `./scripts/check-guard-coverage.sh <dir-nuevo>` sin regresión.
- ☐ **F2.2.18.** Microsatélites y eventos directos: rechazos deterministas.
  - Depende de: F2.2.16.
  - Procede de: DP-IND-03, DP-IND-13; REQ-IND-69/70 nuevos.
  - Cambia: `src/execution.rs` (`admit_repeat_loci_against_small_variants`),
    admisión de `[[direct_events]]` contra la matriz.
  - Se comprueba con: el mismo escenario se rechaza con **dos semillas**
    distintas (prueba de que no depende del sorteo); una fila SNV en la base
    anterior al tracto **sigue admitiéndose** (comportamiento conservado); una
    inserción anclada en esa base se rechaza.

> Revisión de la entrega 1: la matriz con indels se prepara, se reutiliza y
> se rechaza donde debe; las huellas sin indels siguen en su sitio.

### Entrega 2 — Genotipos, copias, verdad, VCF y lecturas

- ☐ **F2.2.19.** Sorteo del locus compuesto y aplicación a las copias.
  - Depende de: F2.2.14, F2.2.18.
  - Procede de: REQ-IND-62 (enmendado), DP-IND-01/02.
  - Cambia: `src/initial_genotypes.rs` (locus compuesto, componentes del
    ámbito con el tramo), `src/initial_individual_integration.rs` y
    `src/small_event_composition.rs` (eventos *indel* de matriz por copia,
    después de los SNV de matriz y en orden decreciente de posición).
  - Se comprueba con: (a) ninguna copia lleva dos alelos de un locus: fixture
    saturado con dos alelos de frecuencia 0,5 cada uno, ploidía 10 y muchos
    individuos, donde un sorteo independiente produciría copias dobles casi
    con certeza; (b) **orden de aplicación**: varias filas *indel* fijadas
    (`alt_af = 1`) en la misma secuencia, una inserción a la izquierda de una
    deleción, y la secuencia de cada copia comparada con una reconstrucción
    **escrita en la prueba** desde las filas; con orden creciente, la
    deleción cae desplazada y la prueba muere; (c) las probabilidades
    genotípicas esperadas publicadas se comparan con las fórmulas de
    REQ-EVO-47/48 calculadas en la prueba, no con el productor.
- ☐ **F2.2.20.** Estadísticas por locus compuesto, configuración efectiva v7 y
  manifiesto v7.
  - Depende de: F2.2.19.
  - Procede de: REQ-IND-65, REQ-MAN-91, DP-IND-09.
  - Cambia: `src/initial_individual_materialization.rs` (observación del
    alelo sobre el tramo), `src/results.rs` (selección de versión por
    contenido y su juez), `schemas/effective-config-v7.schema.json`,
    `schemas/manifest-general-v7.schema.json`, baseline de la puerta de
    cobertura para `src/results.rs`.
  - Se comprueba con: el juez rechaza v7 sin *indels* y v6 con ellos; las
    huellas sin *indels* no se mueven; `check-guard-coverage.sh` sin
    regresión.
- ☐ **F2.2.21.** Juez de la verdad: filas *indel* y conjunto de casos.
  - Depende de: F2.2.19.
  - Procede de: DTC-061 (enmendado), `schemas/truth-v4.md`.
  - Cambia: `src/truth_validation.rs` (una fila *indel* no tiene caso SNP; un
    caso SNP en una posición que comparte con una fila *indel* toma solo los
    alternativos de las filas SNV), `schemas/truth-v4.md`, baseline de los
    módulos `src/truth_validation.rs` y `src/truth_schema.rs`.
  - Se comprueba con: mutación que quita la exención «fila *indel* sin caso
    SNP» → la prueba muere; mutación que deja entrar la fila *indel* en los
    alelos del caso SNP → otra prueba muere.
- ☐ **F2.2.22.** VCF y juez externo.
  - Depende de: F2.2.21.
  - Procede de: DTC-064/065/066.
  - Cambia: corpus de `src/vcf_semantic_validation.rs` y de las pruebas con
    `bcftools` 1.24 y `tabix` 1.24.
  - Se comprueba con: (a) un *indel* de matriz aislado publica en el VCF
    exactamente su fila (`POS`, `REF`, `ALT`) —oráculo independiente, porque
    la matriz se escribe en la prueba—; (b) `bcftools norm` no altera el
    fichero; (c) el locus compuesto sale en un registro fundido y el juez
    reconstruye cada copia desde el `GT`. Comprobar antes que `bcftools` y
    `tabix` del camino de búsqueda son 1.24.
- ☐ **F2.2.23.** Batería de parejas finales contra la verdad.
  - Depende de: F2.2.22.
  - Procede de: DTC-101.
  - Cambia: `tests/final_pairs_against_truth.rs` y, si hace falta,
    `tests/final_pairs_oracle/`.
  - Se comprueba con: dos escenarios nuevos —diploide con *indels* de matriz
    suministrada y errores de secuenciación; tetraploide con un locus
    compuesto saturado (SNV dentro de una deleción en copias distintas)—. Cada
    uno exige que el oráculo **haya visto** una pareja que cruza un *indel* de
    matriz; sin esa exigencia, un escenario cuyas lecturas no tocaran el
    *indel* pasaría sin probar nada.

> Revisión de la entrega 2: revisión adversarial única de la ruta
> suministrada (regla de corte de `CLAUDE.md` §4).

### Entrega 3 — Generador de indels y Marimo

- ☐ **F2.2.24.** Ámbitos de azar y catálogo v11.
  - Depende de: F2.2.13.
  - Procede de: DP-IND-14, DTC-104.
  - Cambia: `catalogs/subseed-scopes-v11.tsv`, `catalogs/README.md`,
    `src/randomness.rs`, `src/repeat_slippage.rs` (literal `"10"`).
  - Se comprueba con: vectores de cada fila nueva con oráculo BLAKE3
    independiente, como las anteriores; la huella del escenario 4 no se mueve;
    variante defectuosa: dejar la constante ligada a la versión del catálogo →
    la huella muere.
- ☐ **F2.2.25.** Generador normalizado por construcción.
  - Depende de: F2.2.16, F2.2.24.
  - Procede de: REQ-IND-61 (enmendado), DP-IND-06/07/13.
  - Cambia: `src/population_snvs.rs` (selección de anclas sin reemplazo con
    la distribución espacial heredada y excluyendo los tramos de DP-IND-13;
    clase, longitud condicionada y contenido).
  - Se comprueba con: una referencia **hecha de homopolímeros** (donde casi
    cualquier deleción libre estaría desalineada) y un normalizador escrito en
    la prueba que comprueba cada fila; quitar la condición `ref[a] ≠ ref[a+k]`
    → la prueba muere. Recuento exacto conservado; sin *indels* declarados,
    huellas 2 y 3 intactas.
- ☐ **F2.2.26.** Reparto por el resto y su auditoría.
  - Depende de: F2.2.25.
  - Procede de: DP-IND-08.
  - Cambia: `src/population_snvs.rs`, auditoría publicada en efectiva v7.
  - Se comprueba con: fixture **saturado** con espectro uniforme en
    `[0,6; 0,99)` y densidades altas para forzar solapes; sin reparto, la suma
    pasa de 1 y la preparación rechaza → la prueba muere. Comprobar también
    que existe al menos un locus compuesto en la salida, o la prueba no prueba
    nada.
- ☐ **F2.2.27.** Documento general v5 y configuración efectiva.
  - Depende de: F2.2.25.
  - Procede de: DP-IND-09, REQ-CONF-129.
  - Cambia: `schemas/config-general-v5.schema.json`, `src/configuration.rs`,
    `src/application.rs`.
  - Se comprueba con: v2/v4 con `generate.indels` se rechazan; v5 sin ella
    funciona; validador de esquemas externo en Python.
- ☐ **F2.2.28.** Paridad Marimo.
  - Depende de: F2.2.27.
  - Procede de: hoja de ruta §2.9 y §7; REQ-UI-14.
  - Cambia: `ui/gps_app.py`, `ui/config_editor.py` (controles de
    `generate.indels`; escribe v5 solo si se usan; vista previa de la matriz
    suministrada con filas *indel* y sumas por locus compuesto),
    `tests_ui/test_ui_contract.py`, `tests_ui/test_population_snvs_editor.py`.
  - Se comprueba con: `uv run --frozen python -m unittest discover -s tests_ui`;
    equivalencia UI→TOML→configuración efectiva; marimo `check --strict`; un
    TOML de la UI **sin** indels sigue siendo v2 o v4.

> Revisión de la entrega 3.

### Entrega 4 — Indels en la población fundadora (F5)

Toca los contratos de F5.1, F5.2, F5.4 y F5.5, cerrados. Todo lo nuevo se
activa solo con un perfil v3.

- ☐ **F2.2.29.** Perfil v3 y catálogo F5 v3.
  - Depende de: F2.2.13.
  - Procede de: DP-IND-10, DTC-105.
  - Cambia: `src/f5_profile.rs`, `schemas/profile-metadata-f5-v3.schema.json`,
    `catalogs/f5-rng-catalog-v3.tsv`, `src/randomness.rs`,
    `src/f5_publication.rs` (`f5_rng_catalog_for`).
  - Se comprueba con: v1/v2 con `indel_sites` y v3 sin ella se rechazan; las
    ocho huellas de `tests/f5_publication.rs` no se mueven.
- ☐ **F2.2.30.** Clase de sitio y forma del *indel*.
  - Depende de: F2.2.29.
  - Procede de: DP-IND-11.
  - Cambia: `src/f5_layout.rs` (después de colocar los sitios),
    `src/f5_marginals.rs` si hace falta.
  - Se comprueba con: perfil saturado (fracción 1, gaps pequeños) donde la
    truncación «cabe antes del sitio siguiente» se ejerce; quitarla → dos
    sitios solapados y la prueba muere. Se publica cuántos sitios quedaron SNV
    por no caber.
- ☐ **F2.2.31.** Fundadores, sitios v2 y verdad F5 nueva.
  - Depende de: F2.2.30.
  - Procede de: REQ-EVO-77 (enmendado), DP-IND-10.
  - Cambia: `src/f5_founders.rs`, `src/f5_artifacts.rs`,
    `schemas/f5-tsv-v1.md` (o su sucesor), esquema de la verdad F5 nueva,
    `src/truth_schema.rs`.
  - Se comprueba con: reconstrucción independiente de cada copia fundadora
    desde los TSV publicados y la referencia.
- ☐ **F2.2.32.** Admisión de recursos v9.
  - Depende de: F2.2.31.
  - Procede de: REQ-EVO-77, DTC-089 (admisión).
  - Cambia: `src/f5_admission.rs`, `src/f5_preflight.rs`.
  - Se comprueba con: la cota de bytes publicados no queda por debajo de lo
    realizado en el fixture saturado con `max_length` máximo; con v8 sobre ese
    fixture, el juez encuentra lo realizado por encima de lo acotado y la
    prueba muere.
- ☐ **F2.2.33.** Oráculo F5, verificación del paquete y auditoría de cierre.
  - Depende de: F2.2.32.
  - Procede de: DTC-089 (oráculo), F5.4/F5.5.
  - Cambia: `src/f5_oracle.rs`, `src/f5_publication.rs`, `src/f5_closure.rs`,
    baseline de la puerta de cobertura de los tres módulos.
  - Se comprueba con: mutación en el productor de la longitud de la deleción →
    el oráculo la rechaza; `check-guard-coverage.sh` sin regresión.

> Revisión de la entrega 4: revisión adversarial única de F5.

### Entrega 5 — Cierre

- ☐ **F2.2.34.** Comparación de bytes y registro de la reapertura.
  - Depende de: F2.2.23, F2.2.28, F2.2.33.
  - Procede de: sección 5; forma de §14 del plan de microsatélites.
  - Cambia: sección 13 nueva en `plan-implementacion-f2-0.1.md` y sección
    nueva en `plan-implementacion-f5-0.1.md` («qué no se reabre», «qué
    cambia y por decisión de quién», «qué bytes se mueven, medido», «dónde
    vive cada literal», «pruebas y su variante defectuosa»).
  - Se comprueba con: binario de `d101f8f` y binario final, en árboles
    separados, sobre los seis escenarios de la sección 5: **cero diferencias**.
- ☐ **F2.2.35.** Puertas completas.
  - Depende de: F2.2.34.
  - Procede de: `CLAUDE.md` §7.
  - Se comprueba con: las nueve puertas, con `uv sync --frozen` hecho y
    `bcftools`/`tabix` 1.24 en el camino; las dos puertas de script en
    directorios nuevos.

> Revisión de cierre: el agente raíz verifica las puertas y las regresiones
> mutando y restaurando por su cuenta. Commit solo con autorización.

## 7. Riesgos

| Paso | Forma del defecto recurrente que más amenaza | Cómo se evita |
|---|---|---|
| F2.2.13 | **Oráculo circular**: «sin *indels*, mismos bytes» comparando dos ejecuciones del binario nuevo | Huellas tomadas con el binario base, en árbol aparte |
| F2.2.16 | **Fixture degenerado**: filas que no se solapan no distinguen un agrupador correcto de uno sin transitividad | Fixture saturado con cadena transitiva y adyacencia sin base común |
| F2.2.18 | **Aserción implicada**: un rechazo que se prueba con una sola semilla podría depender de ella | Dos semillas, mismo diagnóstico |
| F2.2.19 | **Fixture degenerado**: con frecuencias bajas, dos alelos de un locus casi nunca caen en la misma copia aunque el sorteo fuera independiente | Frecuencias 0,5/0,5, ploidía 10 |
| F2.2.19 | **Confusión de coordenadas** referencia/vigentes, que ya costó una ronda | Varias filas fijadas en la misma copia y reconstrucción escrita en la prueba |
| F2.2.20 | Esquema nuevo que el juez nunca recorre | Mutar el selector de versión y ver morir la prueba |
| F2.2.21 | **Guarda descubierta**: el juez acepta el conjunto porque otra comprobación caza antes | Una mutación por regla; si no puede ejercitarse, decirlo en el docstring |
| F2.2.22 | Normalización que el juez externo da por buena por ser la misma que la del productor | Registro esperado escrito a mano desde la fila |
| F2.2.23 | **Escenario que no toca lo que dice cubrir** | Exigir que el oráculo vea la transformación |
| F2.2.24 | Cambio de catálogo que mueve bytes publicados (C-1) | Huella del escenario 4 |
| F2.2.25 | **Fixture degenerado**: en una referencia aleatoria casi todas las deleciones ya están alineadas | Referencia de homopolímeros |
| F2.2.26 | La rama del reparto nunca se ejecuta con espectros realistas | Espectro alto y comprobación de que hay locus compuesto |
| F2.2.30 | Sitios solapados que la búsqueda acepta porque nada los mira | Fracción 1 con gaps pequeños |
| F2.2.32 | Cota de bytes calculada con alelos de una base | Fixture con `max_length` máximo |
| Todas | Puerta medida con el entorno equivocado (`python3`, `bcftools`) | Comprobar versiones antes de dar una cifra |

Un riesgo más, que no es una prueba: **la expectativa del usuario sobre el
realismo**. El generador no pone más indels en repeticiones (DP-IND-07); un
conjunto usado para evaluar detectores de indels en homopolímeros saldría
optimista. Se publica como limitación en la configuración efectiva.

## 8. Tamaño relativo

Escala: S pequeño, M mediano, L grande, XL muy grande. Comparado con las
entregas del plan de microsatélites, donde la entrega 1 sería L.

| Entrega | Tamaño | Por qué |
|---|---|---|
| 0 — documentos y línea base | M | Mucho texto normativo y seis escenarios congelados |
| 1 — matriz suministrada | M | Analizador, agrupación y juez; el sorteo ya es multialélico |
| 2 — genotipos, verdad, VCF, lecturas | L | Toca cinco módulos y dos esquemas; si F2.2.14 falla, XL |
| 3 — generador y Marimo | M–L | Generador nuevo pequeño; la paridad de UI y la versión v5 pesan |
| 4 — F5 | XL | Perfil, catálogo, sitios, verdad nueva, admisión y oráculo de casi 8.000 líneas |
| 5 — cierre | S | Medición y registro |

## 9. Contradicciones encontradas

- **C-1. El catálogo de ámbitos sí se publica.** La §14.1 de
  `plan-implementacion-repeticiones-pcr-0.1.md` dice que la versión del
  catálogo de ámbitos «no se publica en ningún documento». Pero
  `slippage_json` (`src/results.rs:959-978`) escribe
  `rng_scope_catalog_version` con `SUBSEED_CATALOG_VERSION` en la efectiva y el
  manifiesto v6 de todo escenario con microsatélites generados, y el esquema lo
  exige. No es solo una frase inexacta: subir el catálogo por cualquier motivo
  cambia esos bytes. DP-IND-14 propone cómo evitarlo aquí; la corrección del
  texto es aparte.
- **C-2. Rechazo dependiente de la semilla con matriz generada y
  microsatélites.** S-32 justifica juzgar la matriz en la frontera de
  configuración porque «ahí no hay sorteo». Para una matriz **generada** sí lo
  hay: el generador de SNV no excluye los tractos
  (`src/population_snvs.rs:3265-3349`) y la admisión rechaza cualquier fila
  dentro de uno (`src/execution.rs:2772-2789`), de modo que el mismo escenario
  puede aceptarse o rechazarse según la semilla. **Leído en el código, no
  reproducido ejecutando.** Arreglarlo para los SNV movería los bytes de las
  matrices generadas con microsatélites, por eso este plan no lo toca: el
  generador de **indels** sí excluye los tractos desde el principio. Se
  propone como deuda.
- **C-3. DTF-039 frente a la decisión 2.** Generar solapes es generar loci
  multialélicos, que DTF-039 aplaza. Resuelto como decisión pendiente
  (DP-IND-08), no por este plan.
- **C-4. REQ-EVO-42, REQ-EVO-74, DTC-089 y DTF-040 frente a la decisión 3.**
  Los cuatro dicen que F5 no tiene indels. Se resuelve con los cambios de la
  sección 4, que necesitan aprobación.

## 10. Matriz de cierre

| Entrega | Resultado | Evidencia | Riesgo protegido |
|---|---|---|---|
| 0 | Normas aprobadas y línea base congelada | Huellas del binario base que mueren al mutar | Cambio silencioso de bytes |
| 1 | Matriz con indels suministrada y reutilizada | Fixture saturado de solapes; rechazos con dos semillas | Agrupación errónea; rechazos dependientes del azar |
| 2 | Individuos, verdad, VCF y lecturas con indels de matriz | Reconstrucción escrita en la prueba; `bcftools` 1.24; batería de parejas | Copias imposibles; coordenadas desplazadas; VCF o lecturas incoherentes con la verdad |
| 3 | Generador y Marimo | Referencia de homopolímeros; reparto forzado; equivalencia UI→TOML | Filas sin normalizar; sumas mayores que 1; opción sin control |
| 4 | F5 con sitios *indel* | Huellas F5 intactas; oráculo; cota de bytes | Solapes en bloques; cotas falsas |
| 5 | Reapertura registrada | Cero diferencias medidas contra `d101f8f` | Ampliación presentada como añadido silencioso |
