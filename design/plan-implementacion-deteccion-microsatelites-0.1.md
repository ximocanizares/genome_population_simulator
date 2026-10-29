# Plan de implementación — detección de microsatélites sobre la referencia final 0.1

## 1. Estado, autoridad y alcance

- **Fecha:** 28 de octubre de 2026.
- **Estado:** **propuesta, no ejecutable todavía.** Antes de escribir código
  hay que aprobar los cambios normativos de la sección 4 y cerrar las
  decisiones pendientes de la sección 3 que bloquean la entrega 0. La entrega 1
  es una **medición** y termina en una **pausa de decisión**: lo que venga
  después se dimensiona con sus cifras. Ningún paso de la sección 6 está hecho.
- **Qué es.** El desarrollo del hito **F6.4.a** de `hoja-de-ruta-0.1.md`
  («localizar los microsatélites de un genoma suministrado»), que la hoja deja
  sin desarrollar y con cuatro preguntas abiertas: el criterio de
  descubrimiento, qué pasa con un locus declarado que el barrido no encuentra,
  el coste de barrer un genoma completo y el régimen del tartamudeo para los
  descubiertos. Las decisiones del responsable del 27 y 28 de octubre de 2026
  (sección 1.1) responden a parte de ellas y **amplían el alcance**: la
  detección corre sobre toda referencia final, suministrada o simulada.
- **Numeración.** Los pasos heredan la coordenada del hito: F6.4.a.1,
  F6.4.a.2… Las decisiones pendientes se nombran DP-SSR-nn («SSR», del inglés
  *simple sequence repeat*, microsatélite).
- **Orden respecto a otros planes.** Este plan se ejecuta **antes** que
  `plan-implementacion-indels-poblacionales-0.1.md`, por decisión del
  responsable. La sección 5 fija lo que ese plan recibirá de este.
- **Fuentes:** `especificaciones-0.1.md`, `decisiones-tecnicas-0.1.md`
  (DTC-089, DTC-094 a DTC-098, DTC-101, DTC-102), `hoja-de-ruta-0.1.md` (F6,
  §2 y §7), `plan-implementacion-repeticiones-pcr-0.1.md` (decisiones S-03,
  S-06, S-11, S-26, S-30, S-31, S-32),
  `expediente-decisiones-repeticiones-tandem-0.1.md`,
  `plan-implementacion-f6-0.1.md` (§8, paso F6.4.6),
  `registro-ampliaciones-futuras-0.1.md` (DTF-062), y el detector del proyecto
  `pop_var_caller` en `/home/ximo/devel/join_per_sample_vcfs`.
- **Prohibiciones heredadas.** Este plan no autoriza *commit* ni publicación.
  No revierte trabajo ajeno. No se etiqueta una puerta con un comando que no se
  haya ejecutado. Las copias de seguridad de cualquier mutación van fuera de
  git; nunca `git checkout --`. Nada de bucles de sondeo.

### 1.1. Lo que decidió el responsable (27–28/10/2026), fijo para este plan

1. **Tres etapas, siempre en este orden:** construir la referencia (FASTA
   suministrado, reducido si se pide, o secuencia al azar; a las dos se les
   pueden plantar o insertar tractos); **detectar los microsatélites sobre la
   referencia final**; variar (SNV, indels y alelos de microsatélites se
   definen sobre esa referencia). **No puede existir un microsatélite variable
   que no esté en la secuencia de referencia que usa la simulación.**
2. **Inserción de tractos:** sí con genoma suministrado sin matriz o con matriz
   generada; no con matriz suministrada o reutilizada; no con cambios escritos
   a mano; no con reducción en la misma ejecución —se reduce en una, se
   publica el FASTA reducido y se suministra en otra, y el diagnóstico debe
   indicar ese camino—; sí en genoma al azar, plantando tractos.
3. **Un tracto insertado o un locus declarado a mano debe coincidir con uno
   detectado**, o se rechaza con diagnóstico claro. Declarar a mano sirve
   **solo para ajustar parámetros** de un locus concreto (deslizamiento o
   tartamudeo).
4. **Un indel que toque un microsatélite forma un solo locus con él.** Lo
   implementa el plan de indels; este deja el catálogo en forma consumible.
5. **Criterios de detección:** motivos de 1 a 6 bases; número mínimo de copias
   por longitud de motivo, configurable, por omisión `[8, 6, 6, 6, 5, 4]`;
   pureza con la puntuación del propio detector y con la medida de pureza que
   ya usa el simulador.
6. **Parámetros de un locus detectado:** el deslizamiento lo declara el usuario
   **una vez para todos**, sin valor por omisión oculto; un locus declarado a
   mano puede fijar el suyo. El tartamudeo sigue apagado salvo que el usuario
   lo encienda, y encendido usa las tasas calibradas por longitud de motivo.
7. **La primera entrega es una medición** sobre el tomate reducido: cuántos
   loci salen y cuánto ocupan y tardan verdad, VCF y memoria. Con esas cifras
   se decidirá si hace falta limitar la detección a regiones. **El plan deja
   esa decisión como pausa.**
8. **Reutilizar el detector** de `pop_var_caller` copiando su núcleo con
   atribución, sin depender del crate.
9. **Generar microsatélites al azar** en genomas simulados es posterior.

### 1.2. Qué **no** hace este plan

| Fuera | Destino |
|---|---|
| Generar loci de repetición al azar en un genoma simulado | Paso F6.4.6 de `plan-implementacion-f6-0.1.md`; encaja como se describe en la sección 5.2 |
| Indels de la matriz que forman locus con un microsatélite | `plan-implementacion-indels-poblacionales-0.1.md` (DP-IND-13 sustituida) |
| Motivos de más de 6 bases detectados | DTF nueva (sección 4.3) |
| Descomponer automáticamente un tracto compuesto en bloques | DTF nueva; S-30 dice que la estructura en bloques se declara y no se deriva |
| Limitar la detección a regiones | Se decide en la pausa de la entrega 1 (DP-SSR-18) |
| Detección en la población fundadora F5 | Fuera: F5 no tiene microsatélites (REQ-EVO-74; ningún módulo `src/f5_*.rs` usa loci de repetición, comprobado con `grep`) |
| Cambiar el modelo de deslizamiento, el de tartamudeo o sus calibraciones | Sin cambio (DTC-091 a DTC-097, DTC-102) |
| Mutaciones internas de interrupciones | F6.4.5, sin cambio |

---

## 2. Lo que el código hace hoy, comprobado

Cada afirmación cita lo leído o lo ejecutado. Lo que no se ha comprobado se
dice.

### 2.1. El simulador

- **Todo locus lo declara el usuario.** `DeclaredRepeatLocus`
  (`src/configuration.rs:245-287`) lleva identificador, secuencia, tramo
  (`start`/`end`) o punto de inserción (`insert_after`), motivos,
  `reference_units`, `max_units`, mínimos propios y, excluyentes, generador de
  deslizamiento o alelos con frecuencias.
- **El identificador lo pone el usuario** y no hay clase automática (decisión
  S-11; `src/identity.rs:127-129`). Admite de 1 a 64 caracteres ASCII,
  alfanuméricos más `.`, `_` y `-`, empezando por alfanumérico
  (`src/identity.rs:147-156`).
- **La guarda de admisión comprueba que el locus existe** donde se declara:
  pureza sobre bases y tirada más larga en unidades, contra los mínimos del
  locus o los generales (`admit_declared_repeat_loci`,
  `src/execution.rs:2607-2675`), con referencia suministrada y con simulada.
  Los generales por omisión son pureza 0,8 y tirada 3
  (`src/application.rs:656-657`).
- **La medida de pureza del simulador** es `observe_tract`
  (`src/tandem_repeat.rs:863-916`): descomposición voraz de izquierda a
  derecha en copias exactas de los motivos, en la **fase en que se declaró el
  motivo**; una `N` no encaja con ningún motivo.
- **El motivo admite solo A, C, G y T en mayúsculas** (`RepeatMotif::new`,
  `src/tandem_repeat.rs:111-124`). El lector del FASTA pasa todo a mayúsculas
  y admite solo A, C, G, T y N (`src/reference_fasta.rs:297-298`): no quedan
  minúsculas ni códigos IUPAC en la referencia.
- **Tope por omisión de 200 unidades** por locus
  (`DEFAULT_REPEAT_MAX_UNITS`, `src/application.rs:74`); el generador de
  deslizamiento rechaza un tramo de referencia que ya lo rebase (DTC-095,
  REQ-CONF-140).
- **Dos loci que se solapan se rechazan**, comparando cada locus con todos los
  ya resueltos (`src/application.rs:4919-4931`). Es un recorrido cuadrático;
  la sección 2.3 lo mide.
- **El generador de alelos por deslizamiento ya corre en el sitio donde iría la
  detección:** justo después de preparar la referencia y antes de crear nada
  publicado, igual en `validate` y en `run` (`src/execution.rs:143` y `167`,
  `EffectiveConfig::with_generated_repeat_alleles` en
  `src/application.rs:1614`). Devuelve una configuración efectiva nueva o
  `None`, y con `None` se sigue usando la original: es lo que preserva los
  bytes de quien no lo usa.
- **Con simulación, la referencia entera está en memoria.** Si hay
  poblaciones, matriz, cambios directos, evento técnico, VCF o loci
  (`requires_materialized_bases`, `src/execution.rs:2569-2576`), se cargan
  todas las bases. Solo una ejecución que únicamente reduce valida el FASTA
  sin retenerlas (`src/execution.rs:2829-2836`).
- **La inserción ya se rechaza casi como pide la decisión 2**
  (`effective_tract_insertion`, `src/application.rs:2631-2722`): con
  referencia simulada (el tracto se coloca en su tramo), con cambios directos,
  con matriz suministrada o reutilizada, y con reducción. Lo único que cambia
  es el **texto** del rechazo por reducción, que hoy dice
  «la reducción se haría sobre la suministrada, sin el tracto»
  (`src/application.rs:2714-2720`) y no indica el camino en dos ejecuciones.
  Lo prueba `tests/reference_tract_insertion_run.rs`.
- **Una modificación pequeña dentro de un tracto rechaza el escenario**
  (decisión S-32): cambios de base y filas de la matriz en la frontera de
  configuración (`admit_repeat_loci_against_small_variants`,
  `src/execution.rs:2718-2790`), recorriendo cada locus contra cada evento y
  cada fila; inserciones y deleciones en la cadena de segmentos.
- **El generador de la matriz de SNV no excluye los tractos**
  (`select_generated_loci` y `candidate_weight`,
  `src/population_snvs.rs:3265-3349`): solo mira A/C/G/T y las regiones con
  peso. Leído en el código, no reproducido. Es la contradicción C-2 del plan de
  indels, y con detección se vuelve casi segura (sección 2.3).
- **El tartamudeo recorre todos los loci por cada tramo de molécula**
  (`place_tracts_in_segments`, `src/stutter.rs:1467`). Con decenas de miles de
  loci el coste sería lecturas × loci. Leído, no medido.
- **El tartamudeo se configura en el evento técnico**, no en el locus
  (`DeclaredStutter`, `src/configuration.rs:470-477`): interruptor, tasas por
  longitud de motivo, pesos del paso y tasas por locus (`locus_rates`).
- **Versiones publicadas:** documento general v2 → efectiva v4; general v4
  (el que admite loci) → efectiva v6 (`DeclaredDocumentVersion`,
  `src/application.rs:105-130`).
- **El documento general tiene un techo físico de 65 536 bytes**, también fuera
  de F5: el lector acotado `read_bounded_f5_toml` lo aplica a todo documento
  (`src/application.rs:1885-1889`, constante en la línea 60). Comprobado
  ejecutando: un general con 500 loci declarados (91 351 bytes) se rechaza con
  `GPS-INPUT-001`. Unos 300 loci declarados caben. Esto acota cuántos ajustes
  a mano puede declarar un escenario.
- **La arquitectura no describe los microsatélites.** `arquitectura-0.1.md`
  los menciona una sola vez de pasada (línea 442); ningún componente de
  `src/tandem_repeat.rs`, `src/repeat_*.rs` ni `src/stutter.rs` tiene sección.
  Ver la contradicción C-2.

### 2.2. El detector de `pop_var_caller`

- **Dónde está y qué licencia tiene.** `src/tandem_repeat.rs` del proyecto
  `pop_var_caller` (2123 líneas), MIT, «Copyright (c) 2026 Jose Blanca &
  Joaquin Cañizares» (`LICENSE`). Los autores coinciden con los de este
  proyecto.
- **Nota de sala limpia** (líneas 26-28): derivado solo de la idea publicada
  por Benson (1999) —comparar con una copia desplazada un periodo y hacer
  crecer las tiradas de coincidencias— y de código propio; **no** del código
  AGPL-v3 de `TRF-mod`, que no se leyó. Copiarlo no arrastra obligación AGPL.
- **El núcleo** es `find_tandem_repeats(seq, PeriodRange, &ScanParams)`
  (línea 483): para cada periodo `p` puntúa +2 donde `seq[j] == seq[j-p]` y −7
  donde no (una base no A/C/G/T nunca coincide), toma los segmentos de
  puntuación máxima de Ruzzo y Tompa (1999), emite solo motivos **primitivos**
  (`AA` se reduce a `A`) y descarta la re-detección de un tracto con un
  periodo múltiplo cuando el periodo menor cubre la mitad o más. Devuelve
  intervalos semiabiertos en base cero con su periodo y su puntuación, **sin
  el motivo**: el motivo es `seq[start..start+p]`. Por omisión: periodos 1 a
  6, recompensa 2, penalización 7 (tolera una pureza de alrededor de 0,78) y
  dos copias mínimas.
- **Qué arrastra.** El fichero importa `crate::fasta` (línea 42), que solo usan
  el error de lectura `ScanError` y el iterador `RegionScanner`; usa
  `thiserror` en dos enumerados de error (líneas 89-98 y 273-281), y las
  pruebas del recorrido por ventanas usan `tempfile` y el lector de FASTA del
  otro proyecto (líneas 1981-1998). **El simulador no tiene `thiserror`** en
  sus dependencias (`Cargo.toml`). No hay `unsafe` en el fichero (`grep -c`
  devuelve 0).
- **Comprobado compilando.** Extraje las líneas 44-617 a un crate aparte, quité
  `thiserror`, `ScanError` y la parte de regiones, y **compila solo con la
  biblioteca estándar** y `unsafe_code = "forbid"`, con rustc 1.94 (el
  simulador fija 1.97.1 en `rust-toolchain.toml`).
- **Lo reutilizable además del núcleo.** `WindowCursor`, `WindowPlan`,
  `scan_window` y los métodos de `ScannedWindow` (líneas 741-1003) son puros:
  recorren una secuencia por ventanas con margen y atribuyen cada intervalo a
  la ventana donde empieza. `scan_windowed` también es puro (recibe una
  clausura que trae las bases). `RegionScanner` depende del lector del otro
  proyecto y no se copia.
- **Las pruebas del núcleo** (líneas 1140-1480) incluyen una prueba de
  propiedad de Ruzzo–Tompa contra fuerza bruta, casos de sustitución, indel,
  `N`, minúsculas, motivo primitivo, homopolímero y re-detección. Todas usan
  solo el núcleo.
- **Su catálogo aplica más criterios**, en
  `src/region_typing/segment_criteria.rs` y `src/repeat_catalog/criteria.rs`:
  suelo de copias por periodo (`MinCopies`, por omisión `[8, 6, 6, 6, 5, 4]`,
  líneas 412-455, medido en tomate), descarte de motivos «compuestos»
  (`is_compound`, umbral 0,8, línea 1213), descarte de racimos de repeticiones
  a menos de 15 pb (`DEFAULT_BUNDLE_THRESHOLD`, línea 301), recorte de los
  extremos a dos copias limpias (`minimal_trim`, línea 1239), pureza por fase
  cero (`recompute_purity`, línea 1288), flanco mínimo hasta el borde del
  contig, y tope de satélite (500 pb en el catálogo, 1000 en el detector).
  Están pensados para **genotipar lecturas**: un racimo sin flanco limpio no
  se puede anclar. La sección 3 decide cuáles convienen a un simulador.
- **Tiene un oráculo externo pequeño:** `tests/data/tandem_repeat/synthetic_ref.fa`
  y `golden.ssr_catalog.bed.gz` (27 líneas), producido con `trf-mod` 4.10.0
  sobre una referencia sintética generada por
  `generate_synthetic_ref.py`. El catálogo dorado está **después** del filtro
  de GangSTR (periodos 2 a 6, suelos de catálogo, racimos), así que solo sirve
  para el núcleo, no para nuestros criterios.

### 2.3. Sondeo previo, hecho al planificar (no es la entrega 1)

Lo hice para saber si el plan tenía el tamaño correcto. **No sustituye a la
medición de la entrega 1**: usa una copia del detector fuera del repositorio y
una aproximación de los criterios, y mide el formato de publicación de hoy, no
uno decidido. Scripts y salidas en
`/tmp/claude-1000/-home-ximo-devel-genome-population-simulator/953bb011-7c8c-40a4-9795-e622407914a7/scratchpad/plan-ssr/`.

**Referencia:** `results/tomato-sl4-retain-10pct-canonical/resources/references/tomato-sl4-retain-10pct.fa`
(79 230 297 bytes, 78 252 003 bases, 13 secuencias, SHA-256
`f5a8d13a2fb91e0d62157733645b7cf5727ffe7d892451ab994b9f08d3bfea0a`). Se
construyó con la configuración de `inputs/config/general.toml` de ese mismo
directorio: `reduction_percent = 90` y `region_size_bases = 10000` sobre
`/home/ximo/genomes/s_lycopersicum/4.00/S_lycopersicum_chromosomes.4.00.fa`
(795 562 190 bytes). El directorio `results/` no está en git.

**Cuántos tramos salen** (detector con sus valores por omisión; suelo de copias
sobre la longitud detectada; pureza con la descomposición voraz del simulador
en la fase del inicio del tramo, mínimo 0,8, y tirada mínima 3):

| Etapa | Periodo 1 | 2 | 3 | 4 | 5 | 6 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| Intervalos del detector, sin filtrar | | | | | | | 20 288 980 |
| Con el suelo `[8,6,6,6,5,4]` | 53 328 | 6 567 | 1 737 | 190 | 80 | 565 | 62 467 |
| Y pureza ≥ 0,8 y tirada ≥ 3 | 44 564 | 4 586 | 902 | 82 | 26 | 73 | **50 233** |

De esos 50 233, dos superan 200 unidades, 25 pares de vecinos se solapan, y
juntos cubren 562 200 bases (el 0,72 % de la referencia). El barrido tarda
4,6 s y la ejecución entera alcanza 227 MB de memoria residente, sobre una
referencia de 78 MB en memoria. Para el tomate completo, `pop_var_caller`
publica 462 701 loci con estos mismos suelos más su flanco de 30 pb y su tope
de 100 pb (`doc/devel/ng/reports/joint_fit_estimator_2026-08-12.md:64`); no lo
he comprobado.

**Cuánto cuesta publicarlos hoy.** Declaré muestras repartidas de esos loci con
deslizamiento (`mean_events = 1`) sobre 10 individuos diploides, con VCF de
texto, y ejecuté la biblioteca del simulador (`admit_toml_header`,
`parse_declared_config`, `build_effective_config`,
`execute_validated_config`, todas públicas) para saltar el techo de 64 KiB del
lector:

| Loci | Configuración | Preflight | Ejecución | Total | Memoria | Efectiva | Manifiesto | Verdad | VCF |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 300 (binario) | — | — | — | 10,2 s | 163 MB | 1,4 MB | 1,5 MB | 4,8 MB | 0,2 MB |
| 2 000 | 0,9 s | — | — | 28,6 s | 528 MB | 9,3 MB | 9,8 MB | 31 MB | 1,4 MB |
| 4 000 | 3,8 s | 2,1 s | 66,0 s | 71,8 s | — | — | — | — | — |
| 8 000 | 13,7 s | — | — | 290,5 s | 1 826 MB | 37 MB | 39 MB | 124 MB | 5,5 MB |

Tres conclusiones, que la entrega 1 tiene que confirmar o refutar:

1. **Los tamaños crecen en línea recta**: unos 4,6 kB por locus en la
   configuración efectiva, 4,9 kB en el manifiesto, 1,5 kB por locus e
   individuo en la verdad y 0,7 kB por locus en el VCF de 10 muestras. Para
   50 000 loci y 10 individuos serían unos 230 MB de efectiva, 245 MB de
   manifiesto y 760 MB de verdad; para 100 individuos, del orden de 7,5 GB de
   verdad. Es extrapolación.
2. **El tiempo crece de forma cuadrática**: cuatro veces más loci cuestan diez
   veces más. La resolución de la configuración es una causa
   (`src/application.rs:4921`, 0,9 → 3,8 → 13,7 s), pero **la mayor parte está
   en la publicación o en su validación**: el preflight, que prepara la
   referencia, genera los abanicos y materializa los individuos, tarda 2,1 s
   con 4 000 loci, y la ejecución, 66 s. No he localizado dónde.
3. **La memoria también crece más que los ficheros**: 1,8 GB con 8 000 loci.

Consecuencia para el plan: con los criterios decididos, el catálogo del tomate
reducido no se puede publicar con el código de hoy en un tiempo razonable. Por
eso la entrega 1 mide también dónde está el coste, y la pausa que la cierra
decide a la vez las regiones (decisión 7) y la forma de la publicación
(DP-SSR-13).

---

## 3. Decisiones que el responsable no ha fijado

> **Pausa de la entrega 1, 29 de octubre de 2026:** DP-SSR-02 a 08, 13, 17 y
> 18 quedan decididas según las propuestas del apartado 11.7, que recoge la
> decisión.

> **Decisiones del responsable, 28 de octubre de 2026.**
>
> - **DP-SSR-11 — aceptada como paso provisional, con tarea obligatoria
>   anotada.** Un SNV suministrado o escrito a mano dentro de un tracto
>   detectado retira ese tracto (queda sin variación) y se publica la causa;
>   el generador de la matriz excluye los tractos. **TAREA PENDIENTE, que no
>   puede perderse:** cuando el plan de *indels*
>   (`plan-implementacion-indels-poblacionales-0.1.md`) construya los locus
>   compuestos, la misma regla se extiende a los SNV: **un SNV dentro de un
>   microsatélite forma un solo locus con él** y pasa a ser un alelo más,
>   sustituyendo la retirada provisional. Ese plan lo lleva anotado igual.
> - **DP-SSR-09 — aceptada en lo pendiente:** un locus declarado a mano
>   **conserva sus alelos escritos por el usuario**.
> - **Motivos de más de 6 pb (contradicción con REQ-CONF-60) — resuelta:** un
>   locus declarado con motivo de 7 o más bases es una **excepción** a la regla
>   de coincidir con uno detectado, porque el detector solo busca periodos de 1
>   a 6.
> - **DP-SSR-15 — resuelta por el responsable en sentido contrario a la
>   recomendación:** la etiqueta `v0.2.0-test` es una versión de prueba
>   (alfa), **no una publicación**. Los esquemas y formatos se **retocan en
>   su versión actual** en vez de crear versiones nuevas; desaparece la
>   colisión de números de versión con el plan de *indels*. Hay que reflejar
>   esta regla en DTC-036 (versionado) como parte de la entrega 0.
> - **DP-SSR-01 — aceptada por el responsable el 28/10/2026 en la forma
>   propuesta por el agente raíz** (párrafo siguiente), que sustituye a la
>   opción (a) de la tabla. Antecedente: el responsable preguntó si
>   la activación explícita servía solo para conservar resultados o pruebas
>   anteriores. Propuesta del agente raíz del 28/10: la activación es la
>   **presencia de la sección de microsatélites** (que hace falta de todos
>   modos para declarar el deslizamiento común) y hay **un solo modo**: con la
>   sección, se detecta y los declarados solo ajustan; sin ella, no hay
>   microsatélites. Se retira el modo actual de loci declarados sin
>   detección, y las pruebas que lo usan se reescriben. La comparación de
>   bytes contra el commit base se conserva solo como comprobación barata de
>   cambios involuntarios en escenarios sin microsatélites, no como
>   obligación de compatibilidad.

Cada una tiene opciones, consecuencias y una recomendación. **Ninguna se da por
tomada en los pasos**: donde un paso depende de una, lo dice con «(según
DP-SSR-nn)». Las que bloquean la entrega 0 son DP-SSR-01 y DP-SSR-15. Las que
dependen de números (DP-SSR-02 a 08, 13, 17 y 18) se deciden en la pausa de la
entrega 1, con las cifras por etapa delante.

### DP-SSR-01. Activación y compatibilidad

**Hoy:** no hay detección. Los escenarios existentes declaran sus loci a mano,
y algunos son válidos aunque el tracto quede por debajo del suelo de copias:
por ejemplo `STR-A`, un `(ACG)3` de nueve bases sobre referencia simulada,
publicado sin rechazo (`tests/simulated_reference_repeat_tract.rs:1103-1125`),
cuando el suelo por omisión de un trinucleótido es 6.

| Opción | Qué pasa con lo existente | Coste |
|---|---|---|
| **(a) Explícita**: una tabla nueva, `[reference.repeat_detection]`, la activa | Sin la tabla, **mismos bytes** que hoy y mismo régimen: los loci declarados siguen siendo la única fuente | Las decisiones 1 y 3 solo rigen con la tabla. Hay que documentar dos modos |
| (b) Siempre activa | Toda ejecución con referencia gana loci; los declarados que no cumplan el suelo se rechazan; cambian los bytes de casi todo escenario | Reescribir huellas, fixtures y la batería de parejas; romper escenarios de usuarios |
| (c) Activa por versión: un documento general nuevo la lleva siempre; los v2 y v4, nunca | Lo existente no cambia; lo nuevo cumple 1 y 3 siempre | Un escenario nuevo que no quiera microsatélites publica igualmente decenas de miles: con las cifras de la sección 2.3 es inviable |

**Recomendación: (a).** Es la única que conserva los bytes y no impone el
coste del catálogo a quien no lo pide. Se protege con una **huella congelada
contra el commit base** (paso F6.4.a.2). Consecuencia que hay que decir en la
especificación: sin la tabla, un locus declarado no se contrasta con el
detector, como hoy; con la tabla, rigen las decisiones 1 y 3 al pie de la
letra.

### DP-SSR-02. Qué tramo detectado cuenta como microsatélite

La decisión 5 fija los motivos (1 a 6), el suelo y que la pureza se juzga «con
la puntuación del detector y con la medida del simulador». Quedan tres
detalles.

- **Sobre qué se cuentan las copias.** (a) La longitud del tramo dividida por
  el motivo, como el prefiltro del otro proyecto; (b) las **copias completas**
  de la descomposición voraz del simulador (bases cubiertas entre longitud del
  motivo). Con (a) un tramo con interrupciones pasa el suelo con copias que no
  tiene. **Recomendación: (b)**, que es la misma medida con la que el
  simulador anota el alelo; el otro proyecto llegó a lo mismo al contar sobre
  el tramo recortado (`segment_criteria.rs`, comentario de `finish_locus`).
- **La puntuación del detector** actúa ya al formar el segmento (2 y 7 fijos).
  **Recomendación:** no exponerla como parámetro; publicarla como constante
  del algoritmo. La pureza mínima y la tirada mínima son los **mínimos
  generales que ya existen** (`[reference.repeat_admission]`, 0,8 y 3), sin
  una segunda perilla de pureza.
- **El suelo**, configurable como tabla de seis valores por longitud de
  motivo, con `[8, 6, 6, 6, 5, 4]` por omisión, publicado siempre en la
  configuración efectiva aunque no se declare.

### DP-SSR-03. Motivo y rotación

El detector no da el motivo; se lee en el inicio del tramo, de modo que un
`(CA)n` puede salir como `CA` o como `AC` según dónde empiece.

| Opción | Consecuencia |
|---|---|
| **(a) El motivo en la fase del inicio del tramo** | Es el que hace que `observe_tract` mida la pureza correcta: la fase importa, porque la descomposición busca copias exactas desde la izquierda. Dos loci `(CA)n` pueden publicarse uno como `CA` y otro como `AC` |
| (b) Rotación canónica (la menor en orden alfabético) y fronteras sin tocar | `observe_tract` mediría con una fase desplazada: pierde la primera copia y baja la pureza. Habría que ajustar el inicio del tramo, que es cambiar fronteras |
| (c) (a) más un campo con la rotación canónica | Útil para agrupar; nadie lo consume hoy |

**Recomendación: (a).** Nada del simulador agrupa por motivo (las tasas del
tartamudeo van por **longitud** de motivo). El emparejamiento con un locus
declarado (DP-SSR-09) compara rotaciones, no cadenas.

### DP-SSR-04. Fronteras

El tramo de Ruzzo–Tompa puede incluir una copia parcial al final (un
`ACACACA` antes de una `G`) y no tiene por qué coincidir con lo que un usuario
declararía.

| Opción | Consecuencia |
|---|---|
| **(a) Las del detector, sin recortar** | El modelo del simulador admite unidades fraccionarias y la copia parcial baja la pureza, que se mide igual que hoy. Ningún paso nuevo de recorte |
| (b) Recortar a copias completas en los dos extremos | Fronteras más «limpias»; cambia lo que el detector certifica y es un algoritmo más que probar |
| (c) Recorte de GangSTR (`minimal_trim`) | Exige dos copias limpias en cada extremo y descarta lo que no las tiene; es política de genotipado |

**Recomendación: (a).** «¿No sería más fácil?»: sí, y el modelo ya lo admite.

### DP-SSR-05. Tramos solapados de periodos distintos

El detector conserva un tramo corto que vive dentro de otro largo (la `AAA`
dentro de un `(GAAA)n`); el simulador rechaza dos loci solapados
(`src/application.rs:4919-4931`). En el sondeo quedaron 25 pares solapados.

| Opción | Consecuencia |
|---|---|
| **(a) Gana uno, por regla fija**: mayor puntuación del detector; empate, más bases; empate, periodo menor; empate, el de la izquierda. Los demás se descartan y se cuentan | Determinista; no depende del orden de entrada; la puntuación es la medida del detector de longitud y pureza a la vez |
| (b) Se descartan todos los solapados | Pierde el locus bueno con el malo |
| (c) Se funden en un locus compuesto | Exige derivar bloques, que S-30 prohíbe |

**Recomendación: (a).** Dos tramos que solo se tocan (uno acaba donde empieza
el otro) no se solapan y se conservan, como hoy con los declarados.

### DP-SSR-06. Tractos compuestos e impuros

Un `(TCTA)n(TCTG)n` puede salir como un solo intervalo de periodo 4 con pureza
baja, o como dos. El simulador modela compuestos solo si se declaran (S-30).

**Recomendación:** la detección produce **solo loci de un motivo**. Un tramo
que no alcance la pureza mínima con su motivo se descarta y se cuenta. Un
compuesto se declara a mano con `motifs`, y coincide con lo detectado según
DP-SSR-09. No se trae el descarte de motivos «compuestos» de GangSTR
(`is_compound`): el detector ya emite solo motivos primitivos, y lo que esa
regla añade (un motivo que es casi una repetición de otro) lo resuelve la
pureza. Tampoco se trae el descarte de racimos: es una condición para anclar
lecturas al genotipar, no una propiedad del genoma; un simulador debe hacer
variar también los microsatélites que estén juntos.

### DP-SSR-07. Tope de unidades y satélites

Un tramo de más de `max_units` unidades haría fallar la ejecución entera,
porque el generador de deslizamiento lo rechaza (DTC-095).

**Recomendación:** un tope **único para los detectados**, declarable en la
tabla de detección, que por omisión es el de siempre (200 unidades,
`DEFAULT_REPEAT_MAX_UNITS`). Un tramo que lo rebase se descarta como satélite
y se cuenta. Un locus declarado a mano conserva su propio `max_units`.

### DP-SSR-08. `N` y bordes de secuencia

El detector nunca hace coincidir una `N`, de modo que una `N` corta el tramo.
Queda el borde: un tramo que empieza en la primera base o acaba en la última.

| Opción | Consecuencia |
|---|---|
| **(a) Al menos una base a cada lado** | La misma regla que la inserción (REQ-IND-40, `insert_after` de 1 a longitud − 1); el VCF tiene base de anclaje |
| (b) Un flanco mínimo de 15 pb, como el catálogo del otro proyecto | Es una condición de anclaje de lecturas, no del genoma |

**Recomendación: (a).** Cuenta también como borde una `N` inmediatamente
contigua: el tracto sigue siendo real, pero sin base de anclaje conocida. **No
lo he comprobado** contra el VCF de un locus pegado a una `N`; la entrega 3 lo
prueba.

### DP-SSR-09. Cuándo un locus declarado o insertado «coincide» con uno detectado

La decisión 3 exige coincidencia y deja la regla por fijar. Un usuario no sabe
de antemano las fronteras que dará el detector: `(CA)24` declarado en
1001-1048 puede detectarse en 1001-1049 si la base 1049 es una `C`.

**Recomendación:**

1. **Regla de coincidencia.** Un locus declarado coincide si su intervalo (en
   coordenadas finales, ya traducido si hay inserción) **se solapa** con al
   menos un tramo detectado de su misma secuencia cuyo motivo sea **una
   rotación** de uno de los motivos declarados. Sin ninguno, se rechaza
   nombrando el locus y lo detectado más cercano. La comparación se hace con
   los **tramos que pasan los criterios** de DP-SSR-02 a 08: si el declarado
   no cumple el suelo, el diagnóstico dice qué criterio no alcanza.
2. **Qué fronteras mandan.** Las **del locus declarado**, y el declarado
   **sustituye** a todos los detectados que solapa. Así un escenario que hoy
   declara un locus con alelos escritos a mano sigue teniendo el tracto que
   describen esos alelos, y un compuesto declarado ocupa su tramo entero.
3. **Qué puede declarar.** La decisión 3 dice «solo para ajustar parámetros».
   Con detección, un locus declarado admite su **deslizamiento**, su
   `max_units` y sus mínimos; si no declara deslizamiento, **hereda el común**.
   El tartamudeo por locus ya se declara en el evento técnico (`locus_rates`).
   Queda abierta una pregunta que **no** puedo resolver: si un locus declarado
   con detección puede seguir escribiendo sus alelos con frecuencias. Opciones:
   (i) sí, como hoy, porque es otro modo de fijar su variabilidad; (ii) no,
   lectura literal de la decisión 3. Recomiendo (i), porque (ii) deja sin vía
   de expresar un alelo concreto, que es la salida que S-32 ofrece al usuario.
4. **Tractos insertados y plantados.** Igual: tras insertar o plantar, el
   tracto tiene que detectarse con esa regla. Un `(ACG)` plantado en ocho
   bases (`ACGACGAC`, pureza 0,75) ya hoy no pasa la guarda; con detección se
   rechaza con el mismo motivo.
5. **Motivos de más de 6 bases.** REQ-CONF-60 los admite y el detector no los
   ve, de modo que con detección **no podrían coincidir nunca**. Es la
   contradicción C-1. Recomendación: rechazarlos con detección activada y
   diagnóstico que lo explique, y registrar la ampliación (DTF).

### DP-SSR-10. Identificadores de los loci detectados

El ámbito de azar del deslizamiento se indexa por identificador del locus, y
DTC-096 exige que añadir un locus no cambie el abanico de los demás.

| Opción | Consecuencia |
|---|---|
| **(a) Por coordenadas y motivo**: `SSR-<índice de secuencia>-<inicio>-<fin>-<motivo>` en base uno, p. ej. `SSR-2-120345-120360-CA` | Legible; depende solo del propio tramo, así que cumple DTC-096; cabe en 64 caracteres (como mucho 3+10+10+6 cifras y letras más separadores); no estrena clase de identificador |
| (b) Con el nombre de la secuencia | Más legible, pero un nombre de 64 caracteres no cabe |
| (c) Derivado por resumen criptográfico, como `VAR_…` | Clase nueva de identificador y versión nueva del catálogo de identificadores, que se publica |
| (d) Por orden (`SSR-000001`) | Insertar un locus renumera los siguientes y cambia su azar: incumple DTC-096 |

**Recomendación: (a).** El índice es la posición de la secuencia en la
referencia, que ya forma parte de su identidad. Un declarado conserva el
identificador que le puso el usuario; la unicidad se sigue comprobando sobre
el conjunto final, y un declarado que imite el patrón y choque con un
detectado se rechaza como identificador repetido.

### DP-SSR-11. Variantes pequeñas dentro de un tracto detectado

S-32 rechaza todo escenario con un cambio de base o una fila de la matriz
dentro de un tracto. Con unos 50 000 tractos que cubren el 0,72 % de la
referencia:

- **Matriz generada de 1000 SNV:** la probabilidad de que ninguna caiga en un
  tracto es de alrededor de `e^-7,2 ≈ 0,0007`. Sin cambio, casi todo escenario
  se rechazaría según la semilla.
- **Matriz suministrada desde un VCF real:** tendrá SNV dentro de
  microsatélites con toda seguridad.
- **Cambios directos:** los SNP en coordenadas de referencia se comparan; los
  indels, en la cadena de segmentos.

| Opción | Consecuencia |
|---|---|
| (a) Mantener el rechazo | La detección es inservible con matrices reales y frágil con las generadas |
| **(b) Matriz generada: el generador excluye los tractos detectados** de sus candidatos; **suministrada o cambio directo: el tracto detectado que los contiene se retira del catálogo**, y se cuenta y se publica cuál y por qué | La variante del usuario manda sobre lo detectado; el rechazo no depende de la semilla; S-32 sigue intacta para los declarados |
| (c) La variante pasa a ser alelo del locus | Es lo que el plan de indels hará con los indels; para los SNV es un modelo nuevo de locus |

**Recomendación: (b)**, solo con detección activada, de modo que los bytes de
lo existente no cambian. El generador de indels del plan posterior hereda la
misma exclusión. Un locus **declarado** que contenga una variante sigue
rechazándose como hoy (S-32).

### DP-SSR-12. Detección y reducción en la misma ejecución

Con reducción declarada, la simulación corre sobre la **suministrada**, no
sobre la reducida (DTC-098: «la reducción es un subproducto de preparación y
nunca se simula sobre ella»). Detectar en esa ejecución barrería el genoma
entero: unos 460 000 loci en tomate.

| Opción | Consecuencia |
|---|---|
| (a) Permitirla y detectar sobre la suministrada | Coherente con la regla de las tres etapas; sorprende a quien esperaba el catálogo del reducido |
| **(b) Rechazarla**, con el mismo diagnóstico de dos pasos que la decisión 2 pide para la inserción | Simétrico con la decisión 2; evita una ejecución enorme por accidente |

**Recomendación: (b).**

### DP-SSR-13. Qué se publica y con qué forma

| Opción | Consecuencia |
|---|---|
| (a) Cada detectado, como un locus más en `repeat_loci` de la configuración efectiva y del manifiesto, con origen `detected` | Ningún consumidor cambia; con los costes medidos (unos 9,5 kB por locus entre los dos documentos, sin contar los alelos de cada individuo), 50 000 loci son del orden de 475 MB de JSON y el juez tendría que cargarlos |
| (b) Un **catálogo aparte**: TSV con una fila por locus (identificador, secuencia, inicio, fin, motivo, puntuación, origen `detected`/`declared`) más su metadato JSON con criterios, versión del algoritmo, recuentos por motivo de descarte y huella, como el par TSV+JSON de la matriz (DTC-054); la configuración efectiva lleva los criterios, el deslizamiento común (con sus tres limitaciones una sola vez) y el enlace con huella; los abanicos de alelos, en un TSV propio | Los documentos JSON dejan de crecer con los loci; hay que escribir esquema y juez nuevos |
| (c) (a) por debajo de un umbral de loci y (b) por encima | Dos formas de publicar lo mismo: dos jueces, dos caminos que se separan |

**Recomendación:** (b), **pendiente de confirmar en la pausa de la entrega 1**,
porque depende de cuánto pese cada parte medida y de si se limitan regiones.
La verdad y el VCF no cambian de formato en ninguna opción.

### DP-SSR-14. El juez independiente y la puerta de cobertura

El juez que relee lo publicado **no puede llamar al detector**: compararía el
catálogo consigo mismo (el espejo del defecto número uno).

**Recomendación:**

- **Lo que el juez comprueba sin detector**, contra el FASTA final publicado y
  con aritmética propia (sin llamar a `observe_tract` ni a la conversión):
  cada locus está dentro de su secuencia con una base a cada lado; su motivo
  es primitivo y sale de las bases del inicio; sus copias completas alcanzan
  el suelo de su longitud de motivo; su pureza y su tirada alcanzan los
  mínimos; no rebasa el tope; no se solapa con otro; su identificador sale de
  sus coordenadas por la fórmula; cada declarado se solapa con un tracto de un
  motivo rotación del suyo; los recuentos por motivo de descarte suman lo que
  el metadato dice; y cada retirado por una variante (DP-SSR-11) contiene de
  verdad esa variante.
- **Lo que el juez no puede comprobar: la exhaustividad** —que no falte ningún
  tramo que cumpla—. Exige un detector. Se dice así en su documentación, y se
  cubre en pruebas (sección 7) con verdad plantada y con un buscador de
  tiradas exactas por fuerza bruta.
- **Módulo:** uno nuevo y pequeño, `src/repeat_catalog_validation.rs`, con
  cada rechazo por `rejected(..)` (`src/guard_trace.rs`). Entra en la puerta
  de cobertura tocando los cuatro sitios que exige `CLAUDE.md`: el array
  `judges` de `scripts/check-guard-coverage.sh`, una fila en
  `scripts/guard-coverage-baseline.tsv`, la lista de
  `tests/guard_constructors.rs` y el párrafo de `CLAUDE.md`. La alternativa
  es meterlo en `src/results.rs`, ya vigilado; es menos trabajo hoy y deja la
  cobertura del catálogo mezclada con la de seis documentos.

### DP-SSR-15. Versiones de esquema

S-26 permitió completar el general v4, la efectiva v6 y el manifiesto v6
«mientras no se publiquen fuera». El 27/10/2026 se subió a GitHub la etiqueta
`v0.2.0-test` sobre `866bd67`, un ejecutable de pruebas que ya escribe esas
versiones. El plan de indels, por su parte, reserva el general v5 y la
efectiva y el manifiesto v7 para su propio contenido.

| Opción | Consecuencia |
|---|---|
| (a) Ampliar v4/v6 con campos opcionales, que sin detección no aparecen | Lo más barato; mismos bytes sin detección (precedente: DTC-098). Solo vale si `v0.2.0-test` no cuenta como publicación |
| **(b) Versiones nuevas promovidas por contenido**: general v5 solo si declara `[reference.repeat_detection]`; efectiva v7 y manifiesto v7 solo si la ejecución detecta | Mismos bytes sin detección; mismo criterio que la verdad v4/v5 y que el plan de indels; otra rama en la escalera de versiones y su juez |

**Recomendación: (b)**, y que el responsable diga si `v0.2.0-test` cuenta como
publicación. Como este plan va primero, **toma el general v5 y la efectiva y
el manifiesto v7**, y el plan de indels debe ampliarlos (si siguen sin
publicarse) o abrir los siguientes. Es la contradicción C-4.

### DP-SSR-16. Paridad Marimo

La regla del proyecto (`hoja-de-ruta-0.1.md` §2.9 y §7) exige control y
equivalencia UI→TOML→configuración efectiva para toda opción.

**Recomendación:** un panel «Detectar microsatélites» junto al editor de loci
de `ui/config_editor.py`, con el interruptor, los seis suelos (mostrando los
de omisión), el tope y el deslizamiento común **sin valores rellenados**. El
editor de loci existente sigue sirviendo para los ajustes. Con detección, el
editor escribe el general v5 (según DP-SSR-15); sin ella, lo de hoy
(`ui/config_editor.py:2413`). Las pruebas, en `tests_ui/` con `unittest`.

### DP-SSR-17. Escala: cómo se recorre el genoma y dónde está el coste

**Memoria del barrido.** Con simulación la referencia ya está entera en memoria
(sección 2.1), así que el detector no necesita lector propio.

| Opción | Consecuencia |
|---|---|
| **(a) Secuencia a secuencia sobre las bases ya cargadas**, aplicando los criterios a los intervalos de cada una antes de pasar a la siguiente | Lo más sencillo. Pico añadido: una copia canónica de la secuencia más larga y sus intervalos sin filtrar (en el sondeo, 227 MB en total para 78 MB de referencia) |
| (b) Por ventanas con margen (`WindowCursor` y `scan_window` copiados) | Pico acotado por la ventana; hay que probar que el resultado no depende del tamaño de ventana, como hace el otro proyecto |

**Recomendación: (a)**, y (b) solo si la medición lo exige para el genoma
completo.

**El coste cuadrático** de publicar miles de loci (sección 2.3) no lo crea este
plan, pero lo destapa. **Recomendación:** la entrega 1 localiza cada tramo
cuadrático; la pausa decide cuáles se arreglan, en una entrega propia y antes
de conectar la detección, cada arreglo con la huella congelada intacta.

### DP-SSR-18. Limitar la detección a regiones

**No se decide en este plan.** Es la pausa de la decisión 7: se toma con las
cifras de la entrega 1 delante. Si se limita, la tabla de detección ganaría
una lista de regiones en coordenadas finales, y fuera de ellas no habría loci
detectados.

---

## 4. Cambios normativos previos (propuesta)

Se aprueban antes de la entrega 0. El texto definitivo lo escribe el paso
F6.4.a.1. Los números nuevos chocan con los que **propone** el plan de indels,
que tampoco está aprobado (contradicción C-5); como este va primero, los toma.

### 4.1. Especificación (`especificaciones-0.1.md`)

| Elemento | Hoy | Propuesta |
|---|---|---|
| **Nuevo REQ-REF-38** | — | Con la detección activada, el sistema localiza los microsatélites de **la referencia final de la simulación** —suministrada, con tractos insertados o simulada con tractos plantados—, después de construirla y antes de toda variación; ningún locus variable existe fuera de ella. Motivos de 1 a 6 pb |
| **Nuevo REQ-REF-39** | — | Criterios de un tramo detectado (DP-SSR-02 a 08), con el suelo de copias por longitud de motivo `[8, 6, 6, 6, 5, 4]` por omisión, configurable, y descarte contado de lo que no los cumple |
| **Nuevo REQ-REF-40** | — | Con detección, todo locus declarado o insertado coincide con un tramo detectado (DP-SSR-09) o se rechaza, nombrando el locus y el criterio que no alcanza |
| **Nuevo REQ-CONF-142** | — | La tabla `[reference.repeat_detection]` y su deslizamiento común, obligatorio y sin valores por omisión |
| **Nuevo REQ-CONF-143** | — | Variantes pequeñas y tractos detectados (DP-SSR-11) |
| **REQ-CONF-56** | Loci concretos y generados al azar | Añadir los **detectados** como tercera procedencia |
| **REQ-CONF-60** | Motivo de cualquier longitud | Precisar que la detección solo ve 1 a 6 pb y qué pasa con un declarado más largo (DP-SSR-09.5) |
| **REQ-CONF-139** | Un locus sin alelos ni generador se rechaza | Con detección, un locus declarado sin uno ni otro hereda el deslizamiento común |
| **REQ-CONF-132** | El tartamudeo de base está desactivado por omisión | Precisar que, encendido, alcanza a los loci detectados con las tasas calibradas por longitud de motivo |
| **REQ-REF-37** | Rechazo de inserción junto a reducción | El diagnóstico indica el camino de dos ejecuciones; lo mismo para la detección (DP-SSR-12) |
| **Nuevo REQ-MAN-98** | — | Publicación del catálogo, de sus criterios y de los recuentos de descarte (DP-SSR-13) |
| **Nuevo REQ-VAL-64** | — | Juez independiente del catálogo y su límite declarado: no prueba la exhaustividad (DP-SSR-14) |
| **Nuevo REQ-UI-17** | — | Control Marimo de la detección (DP-SSR-16) |
| **Nuevo ACE-119** | — | Aceptación de la detección: bytes intactos sin ella, verdad plantada encontrada, coincidencia, rechazos, juez, batería de parejas |

### 4.2. Decisiones técnicas (`decisiones-tecnicas-0.1.md`)

- **DTC-103 (nueva). Detección de microsatélites:** el núcleo copiado y su
  procedencia (MIT, autores, nota de sala limpia, commit de origen), la
  puntuación fija 2/7, los criterios, la conversión de tramo en locus, los
  identificadores, el momento en que corre (después de preparar la referencia,
  antes del generador de deslizamiento, igual en `validate` y `run`), la
  garantía de bytes sin detección y la versión del algoritmo
  (`ssr-detection-v1`).
- **DTC-104 (nueva). Publicación y juez del catálogo** (según DP-SSR-13 y 14).
- **DTC-095**: el deslizamiento común de los detectados usa el mismo generador;
  sus tres limitaciones se publican **una vez** si se elige DP-SSR-13 (b).
- **DTC-096**: el identificador del detectado depende solo de su tramo, y por
  eso se conserva la propiedad de que añadir un locus no mueve a los demás.
- **DTC-098**: detección después de la inserción; texto del rechazo con
  reducción.
- **DTC-101**: la batería de parejas finales gana escenarios (paso F6.4.a.17).

### 4.3. Registro de ampliaciones futuras

- **DTF-064 (nueva):** detectar motivos de más de 6 pb.
- **DTF-065 (nueva):** descomponer automáticamente tractos compuestos en
  bloques.
- **DTF-066 (nueva, si se elige DP-SSR-11 b):** un SNV dentro de un
  microsatélite como alelo del locus, en vez de retirar el tracto.
- **DTF-067 (nueva, si se elige DP-SSR-12 b):** detectar en la misma ejecución
  que reduce.
- **DTF-062:** añadir que la inserción junto a reducción tiene camino
  documentado en dos ejecuciones.

### 4.4. Hoja de ruta

Desarrollar F6.4.a: enlazar este plan, recoger las decisiones del 27–28/10 que
contestan sus cuatro preguntas, y decir que **el alcance pasa a toda
referencia final**, no solo a la suministrada (contradicción C-6).

### 4.5. Arquitectura

Una sección nueva para el modelo de microsatélites, que hoy no tiene ninguna
(contradicción C-2): el dominio (`src/tandem_repeat.rs`), el deslizamiento
(`src/repeat_slippage.rs`), el sorteo (`src/repeat_genotypes.rs`), el
tartamudeo (`src/stutter.rs`, `src/clonal_stutter.rs`) y el componente nuevo
de detección con su contrato de entrada (bases de la referencia final y
criterios) y de salida (catálogo ordenado, sección 5.1). Este plan no puede
escribir la arquitectura del trabajo ya hecho; propone que la escriba el paso
F6.4.a.1 como requisito de entrada, porque la skill de planificación exige
arquitectura antes de implementar.

### 4.6. Otros

- Trazabilidad: filas de los REQ nuevos y modificados.
- Registro de deudas: el coste cuadrático de publicar muchos loci (sección
  2.3), con lo que la entrega 1 localice.
- Plan de indels: nota de que las versiones y los números DTC cambian (C-4,
  C-5) y de que hereda la exclusión de tractos en su generador.

---

## 5. Lo que este plan entrega a los planes posteriores

### 5.1. Al plan de indels

El plan de indels necesita saber, para cada indel de la matriz, si **toca** un
microsatélite, y entonces hacerlo alelo de ese locus. Este plan le deja:

- **Un solo catálogo**: la lista de loci de la configuración efectiva
  (`EffectiveConfig::repeat_loci()`), detectados y declarados juntos, sin
  solapes, **ordenada por secuencia (en el orden de la referencia) y por
  inicio**, en coordenadas de la referencia final. Cada locus dice su origen
  (`detected` o `declared`).
- **Una consulta por intervalo**: qué loci se solapan con o tocan
  `[inicio, fin)` de una secuencia, por búsqueda binaria sobre ese orden. La
  usan aquí la exclusión del generador de SNV, la retirada por variante
  suministrada (DP-SSR-11) y el juez; el plan de indels la reutiliza para su
  regla de «tocar», que define ese plan (por ejemplo, si cuenta la base
  anterior al tracto, donde se ancla una inserción alineada a la izquierda).
- **La exclusión en el generador** como función reutilizable: el generador de
  indels excluye los mismos tramos que el de SNV.
- Sin detección, el catálogo es la lista de declarados de hoy: la regla del
  plan de indels funciona igual en los dos modos.

### 5.2. A la generación de microsatélites al azar (F6.4.6)

Encaja en la etapa 1 sin tocar la detección: el generador de la referencia
simulada **planta** tractos al azar como regiones declaradas, antes del
relleno (lo que ya hace con los declarados, S-03), y la etapa 2 los detecta
como a cualquier otro. Este plan **no** hace: elegir cuántos, dónde, con qué
motivo ni con qué ámbito de azar; ni exigir que un plantado al azar coincida
con uno detectado (con DP-SSR-09 lo hará solo, porque el plantado es un
declarado).

---

## 6. Entregas y pasos

Las puertas de `CLAUDE.md` §7 se ejecutan al cerrar cada paso que toque
código; en cada paso solo se nombran las comprobaciones propias. **Toda
regresión se verifica fallando contra la variante defectuosa que se nombra,
mutando una copia fuera del árbol real** (memoria del proyecto «mutaciones de
verificación sobre copia»), con SHA antes y después.

### Entrega 0 — normas y línea base

- ✅ **F6.4.a.1.** Textos normativos de la sección 4 escritos y aprobados,
  incluida la sección de arquitectura.
  - Depende de: aprobación de la sección 4 y de DP-SSR-01 y DP-SSR-15.
  - Procede de: sección 4.
  - Cambia: `design/especificaciones-0.1.md`, `design/decisiones-tecnicas-0.1.md`,
    `design/registro-ampliaciones-futuras-0.1.md`, `design/hoja-de-ruta-0.1.md`,
    `design/arquitectura-0.1.md`, `design/trazabilidad-requisitos-0.1.md`.
  - Se comprueba con: revisión del responsable; cada REQ nuevo tiene fila de
    trazabilidad; cada «fuera» de la sección 1.2 tiene destino.
- ✅ **F6.4.a.2.** Huella congelada de lo existente.
  - Depende de: ninguno.
  - Procede de: DP-SSR-01; forma de `tests/pcr_frozen_identity.rs`.
  - Cambia: `tests/repeat_detection_frozen_identity.rs` (nuevo).
  - Se comprueba con: seis escenarios pequeños ejecutados con el **binario
    del commit base** (hoy `d101f8f`), con las SHA-256 de **todos** los
    ficheros publicados escritas como literales: sin loci; loci declarados
    con alelos sobre referencia suministrada; compuesto en simulada con
    `reference_units`; inserción con matriz generada; deslizamiento con
    tartamudeo; PCR con tartamudeo clonal. Uno de ellos lleva un locus por
    debajo del suelo de copias, para que la huella proteja también el régimen
    manual. **Debe fallar contra:** consumir una palabra de azar más al
    colocar el tracto en la referencia simulada, y añadir un campo vacío a
    `repeat_loci` en la efectiva. Se comprueba que la mutación mutó. (El texto
    de un diagnóstico no se publica en ningún fichero y la huella no lo ve:
    lo cubren sus propias pruebas.)

> **Revisión de la entrega 0:** normas aprobadas, huella congelada que muere
> al mutar. No continuar hasta confirmarlo.
>
> **Cerrada el 28/10/2026.** El responsable aprobó los textos y dos retoques:
> el deslizamiento común sigue la regla de REQ-CONF-139 (modelo y media
> obligatorios; sesgo simétrico y paso de una unidad si faltan), y un SNP
> escrito a mano dentro de un microsatélite declarado sigue rechazándose hasta
> DTF-066, que lo resuelve junto con el SNV. La huella
> (`tests/repeat_detection_frozen_identity.rs`, cuatro escenarios sin
> microsatélites) la verificó también el agente raíz sobre una copia aparte:
> un salto de línea de más al final de los JSON publicados hace fallar los
> cuatro.

### Entrega 1 — medición

- ☐ **F6.4.a.3.** Núcleo del detector copiado, sin conectar.
  - Depende de: F6.4.a.1.
  - Procede de: decisión 8; DTC-103; sección 2.2.
  - Cambia: `src/repeat_detection.rs` (nuevo; el nombre `tandem_repeat` ya lo
    usa el dominio del simulador) y `src/lib.rs`. Se copian `PeriodRange`,
    `ScanParams`, `RepeatInterval`, `CANONICAL`, `is_primitive_motif`,
    Ruzzo–Tompa y `find_tandem_repeats`, con sus pruebas de las líneas
    1140-1480 del original. Se quitan `thiserror` (errores con `Display`
    escrito a mano), `crate::fasta`, `ScanError`, la parte de regiones y la
    función muerta `bases_match`. Cabecera con la atribución: proyecto,
    licencia MIT, autores, commit de origen (`950909f9` hoy) y la nota de sala
    limpia literal. Comentarios en castellano; los del original se traducen.
  - Se comprueba con: las pruebas copiadas pasan; `cargo clippy -D warnings`
    con las reglas del simulador; `grep unsafe` vacío. **Debe fallar contra:**
    la búsqueda binaria que el original documenta como errónea en Ruzzo–Tompa
    (la prueba de fuerza bruta la caza) y contra quitar el filtro de motivo
    primitivo (la prueba del homopolímero).
- ☐ **F6.4.a.4.** Conversión de tramos en candidatos, como función pura.
  - Depende de: F6.4.a.3.
  - Procede de: DP-SSR-02 a 08, con la recomendación de cada una de forma
    **provisional** y con los criterios como parámetros.
  - Cambia: `src/repeat_detection.rs`. Recibe las bases de una secuencia y los
    criterios, y devuelve los candidatos y un recuento **por etapa** y por
    longitud de motivo (suelo, pureza, tirada, tope, borde, solape).
  - Se comprueba con: verdad plantada (sección 7, R-1); cada etapa con un
    fixture **saturado** en el que solo esa regla distingue (por ejemplo,
    pureza: dos tramos idénticos en copias y tirada, uno con 0,79 y otro con
    0,81). **Debe fallar contra:** contar copias por longitud y no por copias
    completas; empate de solape resuelto por orden de entrada.
- ☐ **F6.4.a.5.** Medición sobre el tomate reducido.
  - Depende de: F6.4.a.4.
  - Procede de: decisión 7; DP-SSR-13, 17 y 18.
  - Cambia: `tests/repeat_detection_measurement.rs` (nuevo, `#[ignore]`, lee la
    ruta del FASTA de una variable de entorno y no corre en la batería). Y una
    sección nueva, la 11, en este plan con las cifras.
  - Mide, y deja escrito con su comando:
    1. los candidatos por etapa y por longitud de motivo con el suelo por
       omisión y con dos suelos alternativos (más bajo y más alto), para
       decidir DP-SSR-02 a 08 con números;
    2. tiempo y memoria del barrido (`/usr/bin/time -v`);
    3. tiempo, memoria y bytes de **publicar** ese catálogo con el código de
       hoy, declarando los candidatos por la biblioteca —saltando el techo de
       64 KiB del lector, como en el sondeo de la sección 2.3— con 10 y 100
       individuos, con y sin evento técnico con tartamudeo, por fases
       (configuración, preflight, ejecución);
    4. **dónde** está cada coste cuadrático, con `heaptrack` (instalado) o con
       tiempos por fase si no basta; el sondeo ya señala la resolución de la
       configuración y la publicación;
    5. el peso de cada parte de lo publicado por locus, para decidir
       DP-SSR-13.
  - Se comprueba con: la medición se repite dos veces con el mismo resultado de
    recuentos y bytes; las cifras se miden con el binario optimizado
    (`--release`), nunca con el de depuración. Ningún fichero de `src/` cambia
    de comportamiento: la huella de F6.4.a.2 sigue intacta.

> **Pausa de decisión de la entrega 1 (decisión 7 del responsable).** Con las
> cifras delante, el responsable decide: si se limita a regiones (DP-SSR-18);
> los criterios definitivos (DP-SSR-02 a 08); la forma de la publicación
> (DP-SSR-13); y qué costes cuadráticos se arreglan (DP-SSR-17). Las entregas
> 2 a 6 se reescriben con lo decidido. No continuar sin esa decisión escrita.

### Entrega 2 — escala

Reescrita el 29/10/2026 con la decisión de la pausa (apartado 11.7,
DP-SSR-17): **se arreglan los cuatro costes cuadráticos** del apartado 11.5 —el
VCF y su oráculo, la validación física final, el tartamudeo y la lectura del
TOML— antes de conectar la detección, más la comprobación de solapes, que es
barata. El punto 5 de 11.5 (los recorridos localizados solo por lectura, sin
peso medido) **no** se toca. Cada arreglo es de **rendimiento**: no puede
cambiar ni un byte de lo publicado.

**Condición de entrada de todos los pasos: una huella con microsatélites.**
La de F6.4.a.2 solo cubre escenarios **sin** microsatélites, así que no ve
nada de lo que esta entrega toca. Antes de cambiar `src/` se congela
`tests/repeat_scale_frozen_identity.rs`, con los literales obtenidos con el
binario del commit `deecb79` compilado aparte: (1) quince loci en dos
secuencias, dos parejas pegadas, alelos que varían entre individuos, VCF,
tartamudeo de base y quimeras entre secuencias (DTC-102); (2) lo mismo con PCR
y tartamudeo clonal; (3) los veintiún escenarios de la revisión de la fusión
del VCF (`tests/fixtures/vcf-fusion/review`). El juez semántico del VCF tiene
que aceptar lo publicado en (1) y (2). **Debe fallar contra:** un byte de más
en un registro publicado (se comprueba que la mutación mutó).

Los números de paso siguen a F6.4.a.5 y, para no chocar con la entrega 3, el
del VCF y los que la medición añadió llevan sufijo.

- ☐ **F6.4.a.6.** Comprobación de solapes de loci sin comparar cada locus con
  todos los anteriores.
  - Depende de: la huella de entrada.
  - Procede de: DP-SSR-17; sección 2.3.
  - Cambia: `src/application.rs` (`effective_repeat_loci`): los intervalos
    admitidos, por secuencia y por inicio, en un mapa ordenado; el nuevo se
    compara solo con el último que empieza antes que él y con los que empiezan
    dentro, y se nombra el **primero en orden de declaración** de los que
    pisa.
  - Se comprueba con: `tests/repeat_scale_regressions.rs`: tres loci donde el
    tercero solapa al primero y no al segundo dan el mismo diagnóstico de hoy
    (nombra el primero), y uno que pisa a dos nombra el declarado antes, no el
    primero por posición; la huella intacta.
  - **Debe fallar contra:** comparar solo con el vecino declarado
    inmediatamente antes, sin ordenar.
- ☐ **F6.4.a.7.** Índice de tractos por copia para el tartamudeo, de base y
  clonal.
  - Depende de: la huella de entrada.
  - Procede de: DP-SSR-17; apartado 11.5, punto 4.
  - Cambia: `src/stutter.rs` (`StutterTractIndex`: las piezas de cada copia se
    recorren una vez por individua y los tractos quedan ordenados por inicio
    con su fin acumulado; `place_tracts` y `transform_paired_end` reciben el
    índice en vez de la individua) y `src/wgs_sampling.rs` (el muestreador lo
    construye al conectar el plan; el linaje clonal usa el mismo).
  - Se comprueba con: la prueba unitaria
    `the_tract_index_keeps_a_tract_that_starts_before_the_segment_and_ends_inside`;
    la huella intacta, con quimeras entre dos secuencias (DTC-102) y
    tartamudeo de base y clonal.
  - **Debe fallar contra:** un índice que busque el primer tracto por su
    inicio y pierda el que empieza antes del tramo y acaba dentro.
- ☐ **F6.4.a.8a.** VCF: asignar los casos a sus registros y repartir las
  bases de tracto por locus sin recorrer todos los registros ni todos los
  loci.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 1, **corregido por la medición**: el
    coste de la etapa del VCF no estaba en la proyección sino en el
    agrupador (apartado 11.9).
  - Cambia: `src/vcf_grouping.rs` (`finish_records`: por copia, los tramos
    fuente de sus `ALT` ordenados por inicio con su fin acumulado) y
    `src/vcf_publication.rs` (`tract_bases_by_locus`: búsqueda por rango en el
    mapa de loci, conservando el primero en su orden).
  - Se comprueba con: la prueba unitaria
    `carrier_sources_keep_a_long_source_span_that_starts_before_the_locus`; la
    huella intacta.
  - **Debe fallar contra:** buscar el primer tramo candidato por su inicio, que
    pierde el tramo largo que empieza antes del locus.
- ☐ **F6.4.a.8b.** Proyección del VCF y su oráculo semántico: los tractos
  realizados de una secuencia, recogidos una vez.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 1 (`src/vcf_projection.rs:1532, 1743,
    1800`; `src/vcf_semantic_validation.rs:726, 1449, 1566`).
  - Cambia: `src/vcf_projection.rs` (`RealizedTracts`, que sustituye a los
    recorridos de todos los tractos de todas las copias en el registro
    fundido, en `tracts_in`, en la anotación y en el límite de alineación;
    `CopyPath::pieces_near`, que acota con dos búsquedas binarias las piezas
    que mira cada ventana; y `normalization_traversal_closures`, que mira
    solo los componentes que empiezan dentro del tramo atravesado) y
    `src/vcf_semantic_validation.rs`
    (`OracleTracts`, escrito aparte, sin el índice del productor).
  - Se comprueba con: las pruebas unitarias
    `realized_tracts_and_piece_windows_keep_what_starts_before_the_window` y
    `the_oracle_tract_index_keeps_a_tract_that_starts_before_the_window`; la
    huella intacta (los veintiún escenarios de la revisión de la fusión
    ejercitan cada rama: tracto cruzado, anotación por posición, derivada,
    límite de alineación); el juez acepta lo publicado.
  - **Debe fallar contra:** un índice que pierda el tracto que empieza antes de
    la ventana y acaba dentro, en el productor y en el oráculo.
- ☐ **F6.4.a.8c.** Lectura del TOML: la línea de cada valor por búsqueda
  binaria entre los saltos de línea, guardados una vez.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 3.
  - Cambia: `src/configuration.rs` (`AdmittedToml` guarda las posiciones de
    los saltos; `AdmittedToml::location`) y `src/adapters.rs` (la función que
    sitúa cada valor la usa). Los diagnósticos de la cabecera, que se emiten
    una sola vez, siguen con `SourceLocation::from_span`.
  - Se comprueba con: la prueba unitaria
    `the_indexed_location_matches_counting_from_the_start_at_every_offset`
    (todos los desplazamientos de un documento con `CRLF`, caracteres de
    varios bytes y desplazamientos más allá del final) y
    `a_type_error_is_located_at_the_same_line_and_column_after_a_large_document`
    (al final de un documento de 4 000 líneas, en la primera línea y con
    `CRLF`).
  - **Debe fallar contra:** contar el salto de línea que cae justo en el
    desplazamiento (`<=` en vez de `<`).
- ☐ **F6.4.a.8d.** Validación física final: el locus de cada caso y de cada
  registro, por su intervalo.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 2 (`src/results.rs:714` y `:4595`). Las
    comparaciones de árboles JSON enteros (`:4205`, `:4243`) no se tocan: la
    medición no muestra que dominen (apartado 11.9) y el catálogo aparte de la
    entrega 4 dejará de hacerlos crecer con los loci detectados.
  - Cambia: `src/results.rs` (`truth_repeat_tally`: los loci por intervalo y
    las copias por secuencia en mapas; `recount_vcf_tract_bases`: búsqueda
    por rango en el conjunto de loci, conservando el primero en su orden).
  - Se comprueba con: la huella intacta; con la mutación de abajo la
    ejecución misma falla, porque el juez no encuentra el locus de un
    registro.
  - **Debe fallar contra:** un rango que empiece después de la posición del
    registro y pierda el locus que empieza en su primera base.
- ☐ **F6.4.a.8e.** Manifiesto: los sorteos de alelos agrupados por locus una
  vez.
  - Depende de: la huella de entrada.
  - Procede de: la medición de esta entrega (apartado 11.9): la etapa
    «configuración efectiva y manifiesto», que 11.5 midió como «algo más que
    lineal», es cuadrática por `repeat_loci_manifest_json`, que recorría todos
    los sorteos por cada locus. Entra por la regla de la versión anterior de
    esta entrega: «los tramos cuadráticos de la publicación que localice la
    medición, uno por paso».
  - Cambia: `src/results.rs` (`repeat_loci_manifest_json`).
  - Se comprueba con: la huella intacta (publica `drawn_copies` y
    `realized_copies` de cada alelo).
  - **Debe fallar contra:** perder el primer sorteo de cada locus al agrupar.

> **Revisión de la entrega 2:** cada etapa arreglada crece en línea recta al
> doblar los loci, con la huella intacta (apartado 11.9).
>
> **Estado, 29/10/2026:** los siete pasos están implementados, con las
> puertas en verde y cada «debe fallar contra» comprobado mutando una copia
> aparte (apartado 11.9.3). Falta la revisión adversarial; por eso siguen sin
> marcar.

### Entrega 3 — detección en el producto

- ☐ **F6.4.a.9.** Tabla `[reference.repeat_detection]` y su resolución.
  - Depende de: pausa de la entrega 1; F6.4.a.1.
  - Procede de: REQ-CONF-142; decisiones 5 y 6; DP-SSR-01, 02, 07, 15.
  - Cambia: `src/configuration.rs`, `src/adapters.rs`, `src/application.rs`,
    `schemas/config-general-v5.schema.json` (según DP-SSR-15).
  - Se comprueba con: sin la tabla, la huella intacta; con ella, el
    deslizamiento común es obligatorio y su ausencia se rechaza nombrando el
    campo y con un ejemplo en el mensaje (como los pesos del tartamudeo); los
    suelos por omisión aparecen en la efectiva aunque no se declaren; un suelo
    0 o una tabla de cinco valores se rechazan. **Debe fallar contra:** un
    valor por omisión del deslizamiento (la decisión 6 lo prohíbe).
- ☐ **F6.4.a.10.** La detección corre en su sitio.
  - Depende de: F6.4.a.4, F6.4.a.9.
  - Procede de: decisión 1; REQ-REF-38; DTC-103.
  - Cambia: `src/application.rs` (`EffectiveConfig::with_detected_repeat_loci`,
    hermana de `with_generated_repeat_alleles`) y `src/execution.rs:143` y
    `167`: detectar antes de generar abanicos, en `validate` y en `run`; y
    `requires_materialized_bases` cuenta la tabla de detección.
  - Se comprueba con: `validate` y `run` ven el mismo catálogo; una referencia
    con tractos plantados en posiciones conocidas da exactamente esos loci;
    con inserción, el tracto insertado se detecta en coordenadas **finales**;
    con referencia simulada, el plantado. **Debe fallar contra:** detectar
    sobre la referencia de partida en un escenario con inserción (el locus
    sale desplazado la longitud del tracto anterior).
- ☐ **F6.4.a.11.** Coincidencia de declarados e insertados.
  - Depende de: F6.4.a.10.
  - Procede de: decisión 3; REQ-REF-40; DP-SSR-09.
  - Cambia: `src/application.rs`.
  - Se comprueba con: un declarado con fronteras distintas de las detectadas
    pero solapadas coincide y manda sus fronteras; uno de motivo `AC` sobre un
    tracto detectado `CA` coincide; uno sin tracto se rechaza nombrando el
    locus y lo más cercano; uno de 7 pb se rechaza con el diagnóstico de
    DP-SSR-09.5; un compuesto declarado sustituye a los dos tramos que solapa.
    **Debe fallar contra:** emparejar por igualdad de cadena del motivo;
    emparejar por igualdad exacta de fronteras.
- ☐ **F6.4.a.12.** Variantes pequeñas frente a tractos detectados.
  - Depende de: F6.4.a.10.
  - Procede de: REQ-CONF-143; DP-SSR-11; S-32.
  - Cambia: `src/population_snvs.rs` (`select_generated_loci` recibe los
    intervalos excluidos), `src/execution.rs`
    (`admit_repeat_loci_against_small_variants` y la retirada), y la consulta
    por intervalo de la sección 5.1.
  - Se comprueba con: matriz generada de 1000 SNV sobre una referencia con un
    0,7 % de tractos, aceptada con **dos semillas** y ninguna fila dentro; fila
    suministrada dentro de un detectado: el tracto se retira y se publica la
    causa; la misma fila dentro de un **declarado**: rechazo como hoy. **Debe
    fallar contra:** un generador que excluya solo el interior y no la base
    inicial del tracto; una exclusión que cambie el consumo de azar sin
    detección (la huella lo caza).
- ☐ **F6.4.a.13.** Diagnósticos de dos pasos.
  - Depende de: F6.4.a.9.
  - Procede de: decisión 2; REQ-REF-37; DP-SSR-12.
  - Cambia: `src/application.rs:2714-2720` y el rechazo nuevo de detección con
    reducción.
  - Se comprueba con: `tests/reference_tract_insertion_run.rs` actualizado
    exige que el mensaje nombre las dos ejecuciones y el campo `fasta` de la
    segunda. **Debe fallar contra:** el texto de hoy.

> **Revisión de la entrega 3:** la detección se ejecuta en las tres clases de
> referencia, coincide o rechaza los declarados, y no deja rechazos que
> dependan de la semilla. Revisión adversarial de la entrega.

### Entrega 4 — publicación y juez

- ☐ **F6.4.a.14.** Publicación del catálogo.
  - Depende de: F6.4.a.10; DP-SSR-13 decidida.
  - Procede de: REQ-MAN-98; DTC-104.
  - Cambia: `src/results.rs`, `schemas/` (efectiva v7, manifiesto v7 y, con
    DP-SSR-13 b, el esquema del catálogo), `src/execution.rs`.
  - Se comprueba con: los esquemas nuevos validados por el validador de
    esquemas en Python que ya lanzan las pruebas, y la puerta
    `check-independent-validators.sh` sobre un directorio nuevo; recuentos de
    descarte que suman los candidatos; sin detección, la huella intacta.
- ☐ **F6.4.a.15.** Juez independiente del catálogo.
  - Depende de: F6.4.a.14.
  - Procede de: REQ-VAL-64; DP-SSR-14.
  - Cambia: `src/repeat_catalog_validation.rs` (nuevo), y los cuatro sitios de
    la puerta de cobertura.
  - Se comprueba con: cada comprobación de la lista de DP-SSR-14 tiene una
    prueba que **la dispara** con un catálogo publicado alterado a mano (y no
    reproduciendo la conversión); `check-guard-coverage.sh` con el juez dentro.
    **Debe fallar contra:** un juez que recalcule la pureza llamando a
    `observe_tract`. Prueba de independencia: con `observe_tract` mutado en
    una copia del árbol para que sobrestime la pureza, el productor deja pasar
    un locus de pureza real 0,79 y el juez **lo rechaza**.

> **Revisión de la entrega 4:** lo publicado se valida sin el detector.

### Entrega 5 — Marimo

- ☐ **F6.4.a.16.** Panel de detección.
  - Depende de: F6.4.a.9.
  - Procede de: REQ-UI-17; DP-SSR-16; `hoja-de-ruta-0.1.md` §2.9.
  - Cambia: `ui/config_editor.py`, `ui/gps_app.py`, `tests_ui/test_repeat_detection_ui.py`
    (nuevo).
  - Se comprueba con: ida y vuelta UI→TOML→efectiva para cada campo;
    `marimo check --strict`; sin la tabla, el TOML escrito es el de hoy. **Debe
    fallar contra:** un deslizamiento rellenado por el panel.

### Entrega 6 — batería de parejas finales y cierre

- ☐ **F6.4.a.17.** Escenarios nuevos en la batería de DTC-101.
  - Depende de: F6.4.a.14.
  - Procede de: DTC-101; REQ-VAL-03.
  - Cambia: `tests/final_pairs_against_truth.rs`.
  - Escenarios: (1) suministrada con tractos plantados de periodos 1, 2 y 4,
    detección, alelos generados para todos y tartamudeo de base encendido;
    (2) inserción detectada con matriz generada; (3) simulada con plantados y
    PCR con tartamudeo clonal; (4) un declarado que ajusta su deslizamiento y
    su tasa de tartamudeo junto a detectados.
  - Se comprueba con: el oráculo lee los loci de la **verdad publicada**, no del
    catálogo ni del detector; cada escenario exige `repeat_bases_sampled > 0`
    en un locus **detectado** (no solo en el declarado), para que no pueda
    pasar sin ejercitar la detección.
- ☐ **F6.4.a.18.** Puerta externa del VCF sobre un escenario con detección.
  - Depende de: F6.4.a.17.
  - Procede de: DTC-066.
  - Se comprueba con: `scripts/validate_vcf_external.sh` con `bcftools` y
    `tabix` 1.24, incluido un locus pegado a una `N` (DP-SSR-08).
- ☐ **F6.4.a.19.** Demostración y registro.
  - Depende de: F6.4.a.18.
  - Procede de: decisión 7; `hoja-de-ruta-0.1.md` §7.
  - Cambia: sección de cierre de este plan; trazabilidad.
  - Se comprueba con: la ejecución sobre el tomate reducido con lo decidido en
    la pausa, con tiempo, memoria y bytes medidos y escritos con su comando.

> **Revisión de cierre:** revisión adversarial independiente sobre una copia
> extraída con `git archive`, con su propio directorio de compilación
> (`CLAUDE.md` §2).

---

## 7. Riesgos y la forma del defecto número uno en cada paso

| Riesgo | Forma concreta del defecto número uno | Cómo se evita |
|---|---|---|
| **R-1. Oráculo del detector que llama al detector** (F6.4.a.4, F6.4.a.15) | Un oráculo que construye lo esperado con `find_tandem_repeats` o con la conversión: compara la salida consigo misma | Tres oráculos ajenos: (1) **verdad plantada**: fondo aleatorio sin repeticiones de más de tres copias, con tractos en posiciones escritas en la prueba; (2) un **buscador de tiradas exactas por fuerza bruta**, de diez líneas, para tractos puros aislados; (3) el **catálogo dorado de `trf-mod`** del otro proyecto, copiado con su origen, solo para el núcleo: cada locus dorado lo cubre un intervalo del mismo periodo |
| R-2. Pruebas copiadas que no discriminan aquí | Las pruebas del otro proyecto pasan porque prueban su contexto; alguna aserción puede quedar implicada por otra | Mutar el núcleo copiado en dos sitios (F6.4.a.3) y ver morir una prueba distinta en cada uno |
| R-3. Fixture ordinario donde la regla no se ve | El criterio de solape, de pureza o de tope sin fixture que lo aísle: otro criterio descarta antes | Fixture **saturado** por regla (F6.4.a.4); el recuento por etapa demuestra en qué etapa cayó cada tramo |
| R-4. Huella que no puede fallar | Huella calculada en la propia prueba con el binario nuevo | Literales escritos con el binario del commit base (F6.4.a.2), y mutación comprobada |
| R-5. Rechazo dependiente de la semilla | Una prueba con una sola semilla pasa por suerte | Dos semillas en F6.4.a.12; la exclusión en el generador hace imposible la colisión |
| R-6. Coordenadas de partida frente a finales | Detectar en la de partida con inserción: el catálogo sale desplazado y nadie lo nota si el tracto insertado va después de todos los loci | El fixture de F6.4.a.10 inserta **antes** de un locus detectado |
| R-7. Exhaustividad no demostrable | Un tracto puro de muchas copias absorbido en un segmento impuro que falla la pureza desaparece sin aviso | Se dice en la documentación del juez; prueba de propiedad con tractos puros **aislados**; y el caso no aislado se declara en el docstring como no garantizado, sin fingir que discrimina |
| R-8. Batería de parejas que no toca lo detectado | Los escenarios pasan ejercitando solo el locus declarado | `repeat_bases_sampled > 0` en un detectado (F6.4.a.17) |
| R-9. Juez que llama al productor | El juez usa `observe_tract` o la conversión | Prueba de independencia de F6.4.a.15 |
| R-10. Medición con el entorno equivocado | Cifras con el binario de depuración (137 s frente a 3 s en la segunda prueba de uso del 13/10/2026) o con un bucle de sondeo que consume un núcleo | `--release`; esperar con la notificación del trabajo, nunca con bucles |

---

## 8. Tamaño relativo de cada entrega

| Entrega | Tamaño | Por qué |
|---|---|---|
| 0 — normas y línea base | M | Muchos documentos y una sección de arquitectura que no existe; la huella es pequeña |
| 1 — medición | M | Copia del núcleo con pruebas (unas 900 líneas con pruebas), conversión pura y arnés de medición |
| 2 — escala | S a L | Depende de lo que localice la medición; hoy se conoce uno seguro y dos probables |
| 3 — detección en el producto | L | Configuración, conexión, coincidencia, variantes y diagnósticos |
| 4 — publicación y juez | L | Esquemas, escritura, juez nuevo y puerta de cobertura |
| 5 — Marimo | S | Un panel y sus pruebas |
| 6 — batería y cierre | M | Cuatro escenarios, puerta externa y demostración |

---

## 9. Contradicciones encontradas

Pendientes; este plan no las resuelve.

- **C-1. REQ-CONF-60 frente a la decisión 3.** REQ-CONF-60 dice que la
  longitud del motivo no se limita a 1–6 pb; el detector solo ve de 1 a 6, y
  la decisión 3 exige que todo declarado coincida con un detectado. Con
  detección, un declarado de 7 pb o más no puede coincidir nunca. Propuesta en
  DP-SSR-09.5.
- **C-2. La arquitectura no describe los microsatélites.** La skill de
  planificación exige arquitectura antes del plan; `arquitectura-0.1.md` no
  tiene sección para ningún módulo de repeticiones. Propuesta en la sección
  4.5.
- **C-3. Decisión 3 frente a REQ-CONF-139 y a lo que se puede declarar.**
  REQ-CONF-139 rechaza un locus sin alelos ni generador; con detección, un
  declarado de ajuste puede heredar el deslizamiento común. Y la decisión 3
  («solo para ajustar parámetros») no dice si un declarado conserva sus
  alelos escritos a mano (DP-SSR-09.3).
- **C-4. Versiones de esquema.** S-26 permite completar versiones mientras no
  se publiquen fuera, y `v0.2.0-test` está en GitHub desde el 27/10/2026. Y el
  plan de indels reservó el general v5 y la efectiva y el manifiesto v7, que
  este plan, por ir primero, toma (DP-SSR-15).
- **C-5. Numeración de los documentos normativos.** El plan de indels propone
  DTC-103 a 105 y REQ-IND-69 a 71; este propone DTC-103 y 104. Ninguno está
  aprobado; el que se apruebe después se renumera.
- **C-6. Alcance de F6.4.a.** La hoja de ruta lo titula «de un genoma
  suministrado»; la decisión 1 lo extiende a toda referencia final, simulada
  incluida.
- **C-7. Sin matriz no hay individuos.** La decisión 2 admite la inserción
  «con genoma suministrado sin matriz»; al ejecutar un escenario con
  población y sin `[population_snvs]` el simulador responde «la matriz SNV
  falta después de la validación efectiva» (`GPS-CONF-008`, comprobado). Sin
  matriz solo cabe una ejecución sin individuos. No es un fallo, pero la frase
  de la decisión puede llevar a esperar otra cosa.
- **C-8. S-32 frente a la detección.** S-32 se justificó con loci declarados a
  mano, pocos. Con detección, su rechazo sería casi seguro con una matriz
  generada y seguro con una suministrada (DP-SSR-11). Es la contradicción C-2
  del plan de indels, agravada.

---

## 10. Matriz de cierre

| Entrega | Resultado | Evidencia | Riesgo protegido |
|---|---|---|---|
| 0 | Normas aprobadas; huella congelada | Literales del commit base que mueren al mutar | Cambio silencioso de bytes en lo existente |
| 1 | Núcleo copiado; cifras del tomate | Pruebas copiadas que mueren al mutar; medición repetida dos veces | Detector roto al copiarlo; decidir sin números |
| 2 | Publicación que escala | Tiempos lineales medidos; huella intacta | Ejecuciones inviables con el catálogo real |
| 3 | Detección conectada | Verdad plantada encontrada en las tres referencias; dos semillas | Microsatélite fuera de la referencia final; rechazos por azar |
| 4 | Catálogo publicado y juzgado | Juez sin detector, con cada comprobación disparada | Catálogo falso que nadie juzga; espejo |
| 5 | Panel | Ida y vuelta UI→TOML→efectiva | Opción sin control |
| 6 | Parejas finales y VCF | Batería con detectados ejercitados; `bcftools` 1.24 | Lecturas o VCF incoherentes con la verdad |

---

## 11. Medición de la entrega 1

Medida el 28/10/2026 por el implementador de la entrega 1. **Son cifras, no
decisiones**: las propuestas del apartado 11.7 las decide el responsable en la
pausa.

### 11.1. Entorno y forma de medir

- **Código:** `43d0ecd` (la entrega 0, que solo añade documentos y la huella
  a `d101f8f`) más los pasos F6.4.a.3 y F6.4.a.4 sin *commit*; `rustc 1.97.1 (8bab26f4f 2026-07-14)`.
  Ningún fichero de `src/` cambia de comportamiento: `src/repeat_detection.rs`
  es nuevo y no lo llama nada del producto.
- **Binario:** el de pruebas optimizado (`--release`) de
  `tests/repeat_detection_measurement.rs`, copiado fuera de `target/` para que
  ninguna compilación lo pisara. Recuentos y barrido con el de SHA-256
  `c7e468d3…`; publicación con el `7268cb7e…`, que solo difiere en el texto
  del evento técnico (faltaban dos pasos de tartamudeo).
- **FASTA:** `results/tomato-sl4-retain-10pct-canonical/resources/references/tomato-sl4-retain-10pct.fa`,
  SHA-256 `f5a8d13a2fb91e0d62157733645b7cf5727ffe7d892451ab994b9f08d3bfea0a`,
  13 secuencias, 78 252 003 bases.
- **Máquina:** WSL2, 14 núcleos, 27 GB de RAM. Nada más corría durante las
  mediciones (comprobado con `ps`), sin bucles de espera; cada medición se
  lanzó en primer plano o en segundo plano esperando su notificación.
- **Repetición:** todo se ejecutó **dos veces**. Los recuentos de las dos
  pasadas son idénticos línea a línea, y en las diez publicaciones las SHA-256
  de **todos** los ficheros publicados coinciden entre las dos pasadas.
- **Comandos.** `B` es el binario y `FA` el FASTA:

  ```text
  cargo test --release --locked --test repeat_detection_measurement --no-run
  GPS_SSR_FASTA=$FA /usr/bin/time -v $B --ignored --exact measure_load_only --nocapture
  GPS_SSR_FASTA=$FA /usr/bin/time -v $B --ignored --exact measure_scan --nocapture
  GPS_SSR_FASTA=$FA /usr/bin/time -v $B --ignored --exact measure_detection_stages --nocapture
  GPS_SSR_FASTA=$FA GPS_SSR_DIR=<dir> GPS_SSR_LOCI=<n> GPS_SSR_INDIVIDUALS=<i> GPS_SSR_STUTTER=<0|1> \
      $B --ignored --exact prepare_publication --nocapture
  GPS_SSR_TOML=<dir>/general.toml GPS_SSR_DEST=<nuevo> /usr/bin/time -v \
      $B --ignored --exact measure_publication --nocapture
  ```

  Para reproducir cualquier cifra basta con estos comandos y las funciones de
  `tests/repeat_detection_measurement.rs`, cuya cabecera explica cada
  variable. Las publicaciones se lanzaron una tras otra, cada una dos veces,
  con `prepare_publication` una vez y `measure_publication` dos sobre el mismo
  documento. Las salidas brutas y los guiones auxiliares vivían en un
  directorio temporal que ya no existe; no se conservan.

### 11.2. Punto 1: candidatos por etapa y longitud de motivo

Criterios provisionales de DP-SSR-02 a 08 (copias completas, pureza ≥ 0,8,
tirada ≥ 3, tope 200 unidades, una base A/C/G/T a cada lado, desempate de
solapes por puntuación, bases, periodo e inicio). Cada tramo se cuenta en la
**primera** etapa que no supera. `measure_detection_stages`:

| Suelo | Etapa | p1 | p2 | p3 | p4 | p5 | p6 | Total |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| todos | Intervalos del detector | 15 813 346 | 2 781 558 | 1 168 533 | 337 030 | 132 229 | 56 284 | 20 288 980 |
| todos | Motivo con base no ACGT | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| **`[8,6,6,6,5,4]`** | Suelo de copias | 15 766 649 | 2 776 236 | 1 167 318 | 336 932 | 132 200 | 56 138 | 20 235 473 |
| | Pureza | 2 133 | 736 | 313 | 29 | 6 | 71 | 3 288 |
| | Tirada | 0 | 0 | 0 | 0 | 0 | 2 | 2 |
| | Tope | 0 | 0 | 2 | 0 | 0 | 0 | 2 |
| | Borde | 1 | 0 | 1 | 0 | 0 | 0 | 2 |
| | Solape | 7 | 7 | 9 | 1 | 0 | 0 | 24 |
| | **Aceptados** | 44 556 | 4 579 | 890 | 68 | 23 | 73 | **50 189** |
| `[6,4,4,4,3,3]` (bajo) | Aceptados | 217 412 | 34 675 | 5 042 | 901 | 1 965 | 638 | 260 633 |
| | Pureza / tirada / solape | | | | | | | 7 882 / 56 / 652 |
| `[10,8,8,8,7,6]` (alto) | Aceptados | 10 343 | 2 885 | 398 | 12 | 3 | 5 | 13 646 |
| | Pureza / tirada / solape | | | | | | | 1 418 / 0 / 5 |

Bases cubiertas por los aceptados: 559 260 (0,715 % del genoma) con el suelo
por omisión; 2 017 050 (2,58 %) con el bajo; 238 774 (0,31 %) con el alto. Con
el de omisión, entre 582 y 769 candidatos por Mb en cada cromosoma (299 en el
*scaffold* `ch00`). El sondeo de la sección 2.3 contó 50 233: la diferencia
de 44 la explican los solapes (24) y los bordes (2), que el sondeo no
resolvía, los 16 tramos que el sondeo admitía contando copias por longitud y
aquí no pasan por copias completas, y los 2 del tope, que el sondeo no
descartaba.

Cifras auxiliares, con el suelo por omisión:

- **DP-SSR-02.** 8 960 intervalos alcanzan el suelo contando copias por
  **longitud** y no por **copias completas**; solo **16** de ellos (13 de
  periodo 4 y 3 de periodo 5) superarían además pureza y tirada. Elegir una u
  otra forma de contar cambia el catálogo en 16 loci de 50 189.
- **DP-SSR-04.** 2 873 aceptados (5,7 %) acaban en una copia parcial
  (ninguno de periodo 1). Pureza media de los aceptados: 0,9955. Longitud en
  pb: mediana 9, percentil 90 15, percentil 99 55, máxima 459.
- **DP-SSR-07.** Sin tope saldrían 50 190 (uno más: de los dos de periodo 3
  que el tope descarta, uno perdería además un solape). Con más de 50
  unidades, 50; de 100, 6; de 200, 1; ninguno de 500. Máximo: 422 unidades.
- **DP-SSR-08.** Ningún aceptado contiene una `N`; solo 2 se descartan por
  borde. Un solo par de aceptados se toca sin solaparse.

### 11.3. Punto 2: tiempo y memoria del barrido

`/usr/bin/time -v` sobre `measure_load_only` y `measure_scan`, dos pasadas:

| Medida | Pasada 1 | Pasada 2 |
|---|---:|---:|
| Leer el FASTA (memoria residente máxima) | 156 924 kB | 157 052 kB |
| Leer y barrer (memoria residente máxima) | 229 540 kB | 229 752 kB |
| Barrido de las 13 secuencias (tiempo propio) | 5,27 s | 5,44 s |
| Ejecución entera (pared) | 5,39 s | 5,55 s |

El barrido secuencia a secuencia (DP-SSR-17, opción (a)) añade unos **73 MB**
sobre lo que ya ocupa la referencia cargada, y tarda unos **5,3 s**.

### 11.4. Punto 3: publicar el catálogo con el código de hoy

`prepare_publication` toma uno de cada `k` de los 50 189 candidatos, en orden
de referencia, y los declara a mano con deslizamiento `stepwise`,
`mean_events = 1`, una población diploide HWE, VCF de texto y una matriz de
una fila fuera de todo tracto; con tartamudeo, añade un evento técnico sobre
los diez primeros individuos a cobertura 0,1 con el tartamudeo de base
encendido. `measure_publication` lo publica por la biblioteca. «Ejecución»
incluye un segundo preflight interno; lo propio de publicar es ejecución menos
preflight. Pasada 1 / pasada 2:

| Loci | Indiv. | Tartamudeo | Lectura TOML (s) | Resolución (s) | Preflight (s) | Ejecución (s) | Pared | Memoria máx. |
|---:|---:|:---:|---:|---:|---:|---:|---:|---:|
| 1 000 | 10 | no | 0,20 / 0,21 | 0,00 / 0,00 | 0,99 / 1,03 | 15,7 / 15,8 | 0:17 / 0:15 | 317 MB |
| 2 000 | 10 | no | 0,75 / 0,82 | 0,01 / 0,01 | 1,13 / 1,30 | 26,8 / 27,2 | 0:28 / 0:27 | 539 MB |
| 4 000 | 10 | no | 3,15 / 3,18 | 0,03 / 0,04 | 1,89 / 2,02 | 60,8 / 60,7 | 1:03 / 1:03 | 967 MB |
| 8 000 | 10 | no | 12,5 / 13,1 | 0,15 / 0,17 | 4,41 / 4,40 | 244 / 243 | 4:08 / 4:08 | 1 848 MB |
| 1 000 | 100 | no | 0,21 / 0,20 | 0,00 / 0,00 | 7,25 / 7,52 | 147 / 146 | 2:27 / 2:27 | 835 MB |
| 2 000 | 100 | no | 0,83 / 0,78 | 0,01 / 0,01 | 9,52 / 9,81 | 282 / 282 | 4:39 / 4:39 | 1 505 MB |
| 4 000 | 100 | no | 3,44 / 4,16 | 0,03 / 0,05 | 16,5 / 20,0 | 1 081 / 1 165 | 17:36 / 18:55 | 2 887 MB |
| 1 000 | 10 | sí | 0,26 / 0,25 | 0,00 / 0,00 | 1,20 / 1,26 | 111 / 111 | 1:48 / 1:48 | 317 MB |
| 2 000 | 10 | sí | 0,96 / 0,93 | 0,01 / 0,01 | 1,45 / 1,41 | 162 / 168 | 2:44 / 2:36 | 539 MB |
| 4 000 | 10 | sí | 4,35 / 3,55 | 0,03 / 0,03 | 2,41 / 2,31 | 353 / 354 | 5:44 / 5:44 | 968 MB |

«Lectura TOML» es `admit_toml_header` más `parse_declared_config`;
«Resolución», `build_effective_config`. La memoria de las dos pasadas difiere
en menos de 1 MB. Una anomalía del reloj: el tiempo interno del proceso supera
en un 5 % a la pared de `/usr/bin/time` en las ejecuciones largas (261 s frente
a 248 s con 8 000 loci); en WSL los dos relojes no siempre concuerdan. Las
tendencias no cambian con uno u otro.

**El catálogo entero (50 189 loci) no se ha publicado**, porque es inviable
con el código de hoy: extrapolando el último tramo, ya cuadrático (de 4 000 a
8 000 loci la ejecución se multiplica por 4,0), serían del orden de **2 a 3
horas** con 10 individuos y de **dos días** con 100; y la memoria, que crece en
línea recta (unos 0,22 MB por locus con 10 individuos y 0,68 MB con 100), rondaría
**11 GB** y **34 GB**: lo segundo no cabe en esta máquina. Es extrapolación.

### 11.5. Punto 4: dónde está cada coste cuadrático

Tres instrumentos. **(1) Tiempos por fase** (tabla anterior). **(2) Línea de
tiempo de la publicación**, leyendo la hora de modificación de cada fichero
publicado con 2 000, 4 000 y 8 000 loci y 10 individuos (una pasada, con
`find -printf '%T@'` sobre el directorio publicado por `measure_publication`),
que dice cuándo termina cada etapa. **(3) `heaptrack`** con 1 000 y 2 000 loci
(`heaptrack` y `heaptrack_print` sobre `measure_publication`). No hay `perf` ni `gdb` en la máquina, así que **no hay
perfil de CPU por función**: lo que sigue combina la etapa medida con la
lectura del código, y se dice cuál es cuál.

Duración de cada etapa de la ejecución (segundos de pared):

| Etapa (fichero que la cierra) | 2 000 | 4 000 | 8 000 | Crecimiento |
|---|---:|---:|---:|---|
| Preflight interno y copia de entradas y referencia | 1,4 | 2,3 | ≈ 4 | — |
| Verdad (`truth.jsonl`) | 6,2 | 11,0 | 19,4 | lineal |
| **VCF (`truth.vcf`)** | **3,9** | **17,8** | **132,0** | **×4,5 y ×7,4** |
| Configuración efectiva y manifiesto | 0,8 | 2,3 | 11,3 | algo más que lineal |
| **Validación física final** (tras el manifiesto) | **12,0** | **23,5** | **62,6** | ×2,0 y ×2,7 |

Lo localizado:

1. **La proyección del VCF y su oráculo** (medido por etapa; función por
   lectura). Por cada registro con tracto, `tracts_in`
   (`src/vcf_projection.rs:1743`) y la construcción de `realized`
   (`src/vcf_projection.rs:1532`) recorren **todos los tractos de todas las
   copias** de la secuencia, y `annotate_tract_alleles`
   (`src/vcf_projection.rs:1800`) busca la copia y su tracto recorriendo
   listas enteras. El oráculo semántico del VCF repite el patrón
   (`src/vcf_semantic_validation.rs:726`, `:1449` y `:1566`). Coste: registros
   × copias × loci de la secuencia. Es **el mayor**: 132 de 245 s con 8 000
   loci.
2. **La validación física final** (medido por etapa; función por lectura).
   `truth_repeat_tally` busca el locus de cada caso de la verdad recorriendo
   todos los loci (`src/results.rs:4595`), y el recuento de bases del VCF hace
   lo mismo por cada registro (`src/results.rs:714`). Además se reconstruyen y
   comparan enteros, como árboles JSON, los loci de la efectiva y del
   manifiesto (`src/results.rs:4205` y `:4243`).
3. **La lectura del TOML** (medido: 0,20 → 0,75 → 3,15 → 12,5 s; función por
   lectura). `SourceLocation::from_offset` (`src/configuration.rs:113`) cuenta
   los saltos de línea **desde el principio del fichero** para cada valor
   leído: valores × bytes. **Corrige la sección 2.3**, que atribuía el coste
   de la configuración a la comprobación de solapes de
   `src/application.rs:4921`: esa, la de F6.4.a.6, solo cuesta 0,15 s con
   8 000 loci.
4. **El tartamudeo** (medido: el evento con tartamudeo añade 96, 135 y 292 s
   con 1 000, 2 000 y 4 000 loci; función por lectura). Por cada tramo de
   molécula, `place_tracts_in_segments` (`src/stutter.rs:1467`) recorre todos
   los loci de la secuencia, y para cada uno `locate_tract_in_copy`
   (`src/stutter.rs:1526`) recorre los fragmentos de la copia hasta el suyo:
   moléculas × loci², no solo moléculas × loci como decía la sección 2.1.
5. **Por lectura, sin peso medido apreciable hasta 8 000 loci:** aplicar los
   alelos (`apply_repeat_alleles`, `src/repeat_genotypes.rs:295`) busca el
   tramo intacto recorriendo los fragmentos de la copia y copia la lista
   entera de fragmentos por cada locus (`with_current_repeat_allele`,
   `src/segments.rs:830`); la verdad busca el alelo de cada caso recorriendo
   los fragmentos (`src/truth_projection.rs:1177`); y la admisión contra
   variantes pequeñas compara cada locus con cada fila de la matriz
   (`src/execution.rs:2725`), lineal en loci pero multiplicada por las filas.

**Memoria (`heaptrack`).** Las llamadas a asignación crecen en línea recta (49
millones con 1 000 loci, 98 millones con 2 000), y también el pico (263 MB y
447 MB de montón): **la memoria no tiene tramo cuadrático**. El pico lo
dominan los árboles JSON en memoria de la configuración efectiva y del
manifiesto (`serde_json::Value` en `BTreeMap`, unos 180 MB de los 447) y sus
clones al compararlos (`clone_subtree`, 51 MB), además de los 78 MB de la
referencia.

### 11.6. Punto 5: peso de cada parte publicada, por locus

Pendiente de la recta entre 1 000 y 8 000 loci (10 individuos) y entre 1 000 y
4 000 (100 individuos); los bytes son idénticos en las dos pasadas:

| Fichero | 10 individuos | 100 individuos | Para 50 189 loci (10 / 100) |
|---|---:|---:|---:|
| Configuración efectiva (`effective-config.json`) | 4 634 B/locus | 6 036 B/locus | 233 MB / 303 MB |
| Manifiesto (`manifest.json`) | 4 888 B/locus | 6 595 B/locus | 245 MB / 331 MB |
| **Verdad (`truth.jsonl`)** | 15 463 B/locus | 154 544 B/locus | **776 MB / 7,8 GB** |
| VCF de texto (`truth.vcf`) | 682 B/locus | 5 686 B/locus | 34 MB / 285 MB |
| Documento general copiado | 180 B/locus | 180 B/locus | 9 MB |
| Manifiesto del evento técnico | 226 B/locus | — | 11 MB |
| FASTQ (0,1×, 10 individuos) | 0 | — | no depende de los loci |

La verdad crece con locus × individuo, unos **1,55 kB por par**, igual con 10
que con 100 individuos; el VCF, unos 57 B por par con 100 individuos y 68 B
con 10. La configuración efectiva y
el manifiesto crecen con el locus y algo con los individuos (el abanico de
alelos tiene más linajes). Un catálogo aparte en TSV ocuparía del orden de
50 B por locus. **La verdad no es objeto de DP-SSR-13** (su formato no
cambia), y con 100 individuos es ella, no la configuración, lo que manda en
bytes.

### 11.6 bis. Hallazgos y límites de la medición (revisión adversarial)

- **Resuelto el 29/10/2026 (apartado 11.8): un microsatélite pegado a un
  homopolímero de su misma base se perdía entero.** `A×6 (GA)×10` no daba
  ningún intervalo de periodo 2: el segmento de Ruzzo y Tompa empieza dentro
  del poli-A, el motivo leído en su inicio es `AA`, que no es primitivo, y el
  filtro de motivos del núcleo copiado lo descartaba entero. La revisión
  estimó unos **126 tramos** así en el tomate reducido. El responsable decidió
  arreglarlo en el núcleo; el arreglo y sus cifras están en el apartado 11.8.
- **«Una `N` corta el tramo» (DP-SSR-08) es falso.** `(CA)×20 N (CA)×20` da
  **un solo** candidato, con la `N` dentro (pureza 80/81). En el tomate
  reducido no ocurre (cero candidatos con `N` interna, apartado 11.2), pero la
  entrega 3 debe probar la `N` interna y decidir qué hace con ella.
- **No se midió 100 individuos con tartamudeo**: la escalera con tartamudeo
  es solo de 10 individuos (1 000, 2 000 y 4 000 loci).
- **El tercer oráculo del riesgo R-1, el catálogo dorado de `trf-mod`** del
  otro proyecto, **no se copió**: queda para el paso del juez independiente.
  La entrega 1 usa los otros dos (verdad plantada y buscador de tiradas
  exactas por fuerza bruta).

### 11.7. Propuestas para la pausa

> **Decisión del responsable, 29 de octubre de 2026: se aceptan todas las
> propuestas de este apartado tal como están escritas.** En concreto:
> regiones opcionales con el genoma entero por omisión (DP-SSR-18); los cuatro
> costes cuadráticos se arreglan en una entrega propia antes de conectar la
> detección (DP-SSR-17); catálogo aparte en TSV con su metadato JSON
> (DP-SSR-13, opción (b)); y los criterios de DP-SSR-02 a 08 como se proponen
> abajo. Además, el defecto del núcleo copiado que perdía un microsatélite
> pegado a una tira de su misma base **se corrige** en el núcleo (apartado
> 11.8), y el inicio recortado sale de la tira para que la tira y el
> microsatélite solo se toquen.

- **DP-SSR-02:** copias completas (opción (b)), como recomendaba el plan. La
  medida dice que la elección apenas cambia el catálogo (16 loci de 50 189), así
  que conviene la que coincide con cómo anota el simulador. Suelo
  `[8, 6, 6, 6, 5, 4]` por omisión: da 0,7 % del genoma, y el bajo lo
  quintuplica.
- **DP-SSR-03:** motivo en la fase del inicio del tramo (opción (a)). Ningún
  motivo detectado llevó una base no ACGT.
- **DP-SSR-04:** fronteras del detector sin recortar (opción (a)); el 5,7 % con
  copia parcial final lo absorbe el modelo de unidades fraccionarias.
- **DP-SSR-05:** la regla fija propuesta (opción (a)). Solo hay 24 solapes:
  la regla importa por determinismo, no por volumen.
- **DP-SSR-06:** solo loci de un motivo, como se recomendaba; la pureza
  descarta 3 288 tramos, que es el techo de lo que una descomposición en
  bloques podría recuperar.
- **DP-SSR-07:** tope común de 200 unidades, declarable. Solo descarta 2
  tramos; el mayor tiene 422 unidades.
- **DP-SSR-08:** una base A/C/G/T a cada lado (opción (a)); solo 2 descartes
  y ningún tracto con `N` dentro.
- **DP-SSR-13:** catálogo aparte (opción (b)). Con 50 189 loci, la efectiva y
  el manifiesto sumarían 480–630 MB de JSON que hoy se cargan enteros en
  memoria como árboles y se comparan enteros en la validación; en TSV serían
  unos 2,5 MB.
- **DP-SSR-17:** barrido secuencia a secuencia (opción (a)): 5,3 s y 73 MB
  sobre la referencia cargada. Arreglar en la entrega 2, en este orden de
  peso medido: el VCF y su oráculo, la validación física final, el tartamudeo
  y la lectura del TOML. La comprobación de solapes de F6.4.a.6 es la menor de
  todas (0,15 s con 8 000 loci); se puede hacer de paso, pero no resuelve nada.
- **DP-SSR-18:** ofrecer una **lista opcional de regiones**, con el genoma
  entero por omisión. Aun con todos los tramos cuadráticos arreglados, lo
  lineal del catálogo entero con 100 individuos (unos 7,8 GB de verdad y unos
  34 GB de memoria al ritmo medido) no cabe en esta máquina; con 10
  individuos sí cabría. Si la meta es 100 individuos sobre el genoma entero,
  hace falta además bajar la memoria por locus, que hoy dominan los árboles
  JSON.

### 11.8. Arreglo del inicio en fase primitiva y nueva medición del barrido (29/10/2026)

**Qué cambia en el núcleo copiado** (decisión del responsable; única
desviación del algoritmo original, documentada en la cabecera de
`src/repeat_detection.rs`). Si el motivo leído al inicio de un segmento no es
primitivo, el original descartaba el intervalo entero; ahora el inicio
**avanza base a base hasta la primera fase cuyo motivo de periodo `p` sea
primitivo**, y solo se descarta si no queda ninguna antes de que el tramo baje
del mínimo de copias del detector. La puntuación del tramo recortado se
**recalcula** (la del segmento menos la del prefijo quitado), de modo que
siempre es la suma de las comparaciones del propio tramo, que es lo que usa el
desempate de solapes. Un tramo cuyo motivo inicial ya era primitivo sale
idéntico; uno puro de motivo no primitivo (alias de un homopolímero o de un
dinucleótido) no tiene ninguna fase primitiva y se sigue descartando. Ejemplos,
con sus pruebas: `A×6 (GA)×10 C` da `[5, 26)` con motivo `AG`; `A×8 (GA)×10`
da el poli-A `[0, 8)` y `[7, 28)` con motivo `AG`; `T×12 (CTT)×7` da el poli-T
`[0, 12)` y `[10, 33)` con motivo `TTC`; `A×4 (GA)×10`, que ya funcionaba, no
cambia.

**Nueva medición** (mismos comandos del apartado 11.1, binario de pruebas
optimizado de SHA-256 `71c54c58…`, dos pasadas con recuentos idénticos línea a
línea). Antes → ahora, suelo por omisión:

| Etapa | p1 | p2 | p3 | p4 | p5 | p6 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| Intervalos del detector | 15 813 346 → 15 813 346 | 2 781 558 → 2 784 397 | 1 168 533 → 1 170 281 | 337 030 → 338 635 | 132 229 → 132 490 | 56 284 → 56 606 | 20 288 980 → **20 295 755** |
| Suelo de copias | 15 766 649 → 15 766 649 | 2 776 236 → 2 778 949 | 1 167 318 → 1 169 053 | 336 932 → 338 484 | 132 200 → 132 461 | 56 138 → 56 445 | 20 235 473 → 20 242 041 |
| Pureza | 2 133 → 2 133 | 736 → 752 | 313 → 316 | 29 → 59 | 6 → 6 | 71 → 80 | 3 288 → 3 346 |
| Tirada | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 2 | 2 → 2 |
| Tope | 0 → 0 | 0 → 0 | 2 → 2 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 2 |
| Borde | 1 → 1 | 0 → 0 | 1 → 1 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 2 |
| Solape | 7 → 60 | 7 → 27 | 9 → 9 | 1 → 10 | 0 → 0 | 0 → 1 | 24 → **107** |
| **Aceptados** | 44 556 → **44 503** | 4 579 → **4 669** | 890 → **900** | 68 → **82** | 23 → 23 | 73 → **78** | 50 189 → **50 255** |

Con los otros suelos: el bajo pasa de 260 633 a 260 663 aceptados (solapes de
652 a 1 457) y el alto de 13 646 a 13 682 (solapes de 5 a 24). Bases en
candidatos con el suelo por omisión: 559 260 → 561 327 (0,715 % → 0,717 %).
Cifras auxiliares: copia parcial final 2 873 → 2 937; pares que se tocan 1 →
2; pureza media 0,9955 → 0,9954; DP-SSR-02, 8 960 → 9 083 intervalos que pasan
por longitud y no por copias completas, y siguen siendo **16** los que
superarían además pureza y tirada; DP-SSR-07 y la longitud de los tramos, sin
cambio. Barrido: 5,3 s → **6,2 y 6,6 s** (pasadas 1 y 2), memoria máxima
229 540 kB → 229 268 y 229 780 kB.

**Lo que hay que saber de estas cifras, sin decidirlo aquí.** Se ganan 119
microsatélites de periodo 2 a 6 (90, 10, 14, 0 y 5), pero se **pierden 53
homopolímeros** que antes se aceptaban: los descartes por solape de periodo 1
suben exactamente en 53 (de 7 a 60). La causa, comprobada en las pruebas con
`A×8 (GA)×10`: el tramo recortado empieza en la **última base del
homopolímero** (`[7, 28)`, fase `AG`), los dos se solapan en esa base y la
regla de solapes de DP-SSR-05 deja solo el de más puntuación. El saldo es
+66. Avanzar el
inicio hasta la primera fase primitiva **que no empiece dentro del
homopolímero** (`[8, 28)`, fase `GA`, que solo toca al poli-A) conservaría los
dos, pero es una regla de fronteras (DP-SSR-04) y de solapes (DP-SSR-05), así
que la decide el responsable, no el implementador.

**Tercera versión: el inicio recortado sale de la tira (decisión del
responsable, 29/10/2026).** Un intervalo cuyo inicio se recortó no puede
empezar dentro de la tira de una sola base que lo precede: avanza hasta la
primera fase primitiva **cuya base sea distinta de la anterior**. Los
intervalos no recortados no cambian. `A×8 (GA)×10` da el poli-A `[0, 8)` y
`[8, 28)` con motivo `GA`, que solo se tocan; `T×12 (CTT)×7` da `[0, 12)` y
`[12, 33)` con motivo `CTT`; `A×6 (GA)×10 C` da `[6, 26)`. Medido con el
binario de SHA-256 `5bc796da…`, dos pasadas idénticas línea a línea. Original →
primera fase primitiva → tercera versión, suelo por omisión:

| Etapa | p1 | p2 | p3 | p4 | p5 | p6 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| Intervalos del detector | 15 813 346 (igual en las tres) | 2 781 558 → 2 784 397 → 2 784 387 | 1 168 533 → 1 170 281 → 1 169 475 | 337 030 → 338 635 → 338 019 | 132 229 → 132 490 → 132 284 | 56 284 → 56 606 → 56 536 | 20 288 980 → 20 295 755 → **20 294 047** |
| Pureza | 2 133 (igual) | 736 → 752 → 744 | 313 → 316 → 315 | 29 → 59 → 59 | 6 (igual) | 71 → 80 → 80 | 3 288 → 3 346 → 3 337 |
| Solape | 7 → 60 → **7** | 7 → 27 → 20 | 9 (igual) | 1 → 10 → 10 | 0 (igual) | 0 → 1 → 1 | 24 → 107 → **47** |
| **Aceptados** | 44 556 → 44 503 → **44 556** | 4 579 → 4 669 → **4 651** | 890 → 900 → **895** | 68 → 82 → **81** | 23 (igual) | 73 → 78 → **78** | 50 189 → 50 255 → **50 284** |

Tirada, tope y borde no cambian (2, 2 y 2). Suelo bajo: 260 633 → 260 663 →
261 122 aceptados; suelo alto: 13 646 → 13 682 → 13 690. Bases en candidatos:
559 260 → 561 327 → 561 431. Pares de candidatos que se tocan: 1 → 2 → **48**.
Barrido: 5,3 s → 6,2–6,6 s → 7,1–8,1 s (pasadas 1 y 2); memoria máxima, 230 MB
en las tres.

**Se recuperan las 53 tiras, pero de los 119 microsatélites ganados se quedan
95.** Comparando, locus a locus, los candidatos de la segunda y la tercera
versión: los 53 homopolímeros vuelven (los descartes por solape de periodo 1
bajan de 60 a 7, los del original); 83 microsatélites siguen, con el inicio
desplazado fuera de la tira; y **32 microsatélites de la segunda versión
desaparecen** (26 de periodo 2, 5 de periodo 3 y 1 de periodo 4). Todos
alcanzaban su suelo de copias **solo contando bases de la tira**. Ejemplos del
tomate: `…TTTTTTTTT ATATATATATA…`, que en la segunda versión era `(TA)×6`
empezando en la última `T` y ahora es `(AT)×5 A`, por debajo del suelo de 6;
y `…AAAAAAAAA GAAGAAGAAGAAGAAGA…`, que era `(AAG)×6` empezando en las dos
últimas `A` y ahora es `(GAA)×5 GA`, por debajo de 6. En el segundo tipo (el
motivo empieza por la base de la tira) no hay forma de saber dónde acaba la
tira y empieza el microsatélite. **Decisión del responsable, 29/10/2026: se
acepta esa pérdida.** Esos 32 solo alcanzaban el suelo contando dos veces una
base ambigua, el núcleo original tampoco los encontraba, y recuperarlos
exigiría contar la base compartida para el suelo pero no para las
coordenadas, una regla más compleja que no compensa.

### 11.9. Entrega 2: los costes cuadráticos arreglados, y lo que queda (29/10/2026)

Medido el 29/10/2026 por el implementador de la entrega 2. Mismos comandos
que el apartado 11.1 (`prepare_publication` y `measure_publication`), mismo
FASTA y misma máquina, con el binario de pruebas optimizado (`--release`)
copiado fuera de `target/`.

**Una salvedad sobre «antes».** Los documentos se prepararon de nuevo con el
detector de hoy (50 284 candidatos, apartado 11.8), y no con el de 11.4
(50 189): con el mismo número de loci, los loci elegidos difieren un poco. Por
eso «antes» no se toma de 11.4: se midió otra vez, **sobre los mismos
documentos**, con el binario de `deecb79` sin cambios (SHA-256 `5bc796da…`,
el mismo de 11.8).

#### 11.9.1. Dónde estaba de verdad cada coste

La medición de 11.5 localizaba las etapas por la hora de los ficheros y las
funciones **por lectura**. Para esta entrega se midió el tiempo propio de cada
función sospechosa en una copia instrumentada aparte (cronómetros
acumulados por función; no hay `perf` en la máquina), con 4 000 y 8 000 loci
y 10 individuos. Lo que encontró, y que corrige 11.5:

Segundos de pared acumulados por función, copia instrumentada de `deecb79`
(una pasada; la máquina tenía otros trabajos en marcha, así que las cifras
absolutas salen infladas y solo vale la proporción entre 4 000 y 8 000):

| Función (fichero) | 4 000 | 8 000 | Crecimiento |
|---|---:|---:|---|
| **Etapa del VCF entera** (`publish_truth_vcf_within`) | 30,3 | 218,4 | ×7,2 |
| — asignar cada caso de la verdad a sus registros (`finish_records`, `src/vcf_grouping.rs`) | **27,5** | **212,7** | ×7,7 |
| — la proyección entera (`project_view_from_validated_truth`, `src/vcf_projection.rs`) | 0,82 | 2,88 | ×3,5 |
| —— de ella, los tractos de todas las copias por registro fundido (`:1532`) | 0,48 | 2,00 | ×4,2 |
| —— el límite de alineación (`left_alignment_floor`) | 0,10 | 0,31 | ×3,2 |
| — repartir las bases de tracto por locus (`tract_bases_by_locus`, `src/vcf_publication.rs`) | 0,08 | 0,42 | ×5,1 |
| **Validación física final** (`validate_published_resources`, dos llamadas) | 27,6 | 82,0 | ×3,0 |
| — releer y resolver el TOML copiado (`validate_effective_targets_against_configuration`) | 7,6 | 30,4 | ×4,0 |
| — el recuento de la verdad por locus (`truth_repeat_tally`, `src/results.rs:4595`) | 1,3 | 9,0 | ×6,8 |
| — rehacer las individuas iniciales (`validate_physical_initial_population_statistics`) | 3,8 | 9,6 | ×2,5 |
| — el recuento de bases del VCF (`recount_vcf_tract_bases`, `src/results.rs:714`) | 0,10 | 0,48 | ×5,0 |
| **Lectura del TOML** (`parse_declared_config`) | 3,1 | 14,8 | ×4,8 |

Y con tartamudeo, 2 000 loci: situar los tractos en cada molécula
(`place_tracts_in_segments`) sumaba 44,5 s de los 147 s de la ejecución.

Lo que esto corrige de 11.5:

1. **El coste de la etapa del VCF no estaba en la proyección** que 11.5
   señalaba por lectura (`src/vcf_projection.rs:1532, 1743, 1800`), sino en el
   agrupador: por cada portadora de cada caso de la verdad, `finish_records`
   recorría **todos los registros** del VCF (casos × registros). La proyección
   también era cuadrática, pero pesaba 2,9 s de 218 con 8 000 loci. Se
   arreglaron las dos (F6.4.a.8a y F6.4.a.8b). Es la misma etapa y la misma
   forma de coste que la decisión de la pausa manda arreglar; solo estaba en
   otra función.
2. **La mayor parte de la validación física cuadrática era la lectura del
   TOML**, que el juez repite sobre la copia publicada: arreglar
   `SourceLocation` (F6.4.a.8c) la arregla también ahí. Los dos recorridos de
   `src/results.rs` que 11.5 nombraba también lo eran (F6.4.a.8d). Las
   comparaciones de árboles JSON (`:4205`, `:4243`) no aparecen como coste
   propio: van dentro de `validate_published_repeat_loci`, que creció ×2,1
   de 4 000 a 8 000 loci tras los arreglos; no se tocaron.
3. **La etapa «configuración efectiva y manifiesto»**, que 11.5 midió como
   «algo más que lineal» (0,8 → 2,3 → 11,3 s), es cuadrática:
   `repeat_loci_manifest_json` recorría todos los sorteos de alelos por cada
   locus. Se arregló (F6.4.a.8e).
4. **Rehacer las individuas iniciales crece ×2,5 al doblar loci**, en el
   preflight y en la validación. Es el punto 5 de 11.5 (`apply_repeat_alleles`
   y `with_current_repeat_allele`), que la decisión de la pausa deja fuera:
   **no se tocó** (apartado 11.9.6).


#### 11.9.2. Qué cambió en cada paso

El detalle está en cada paso de la entrega 2 (sección 6). En todos la idea
es la misma y sencilla: ordenar una vez por inicio, guardar el fin máximo
acumulado y buscar con dos búsquedas binarias, o agrupar en un mapa, en vez de
recorrer la lista entera por cada consulta. El fin **acumulado** es lo que
impide perder el tramo que empieza antes de la ventana y acaba dentro, que es
la forma del defecto que cada prueba de regresión caza. Ninguna estructura
genérica nueva: cada módulo tiene su índice de veinte o treinta líneas, y el
del juez del VCF está escrito aparte del del productor.

| Paso | Ficheros | Qué recorría antes | Qué consulta ahora |
|---|---|---|---|
| F6.4.a.6 | `src/application.rs` | cada locus contra todos los anteriores | el último que empieza antes y los que empiezan dentro |
| F6.4.a.7 | `src/stutter.rs`, `src/wgs_sampling.rs` | por tramo de molécula, todos los loci y las piezas de la copia hasta cada uno | los tractos de la copia, situados una vez por individua |
| F6.4.a.8a | `src/vcf_grouping.rs`, `src/vcf_publication.rs` | por portadora, todos los registros; por registro, todos los loci | por copia, sus tramos fuente; los loci por rango |
| F6.4.a.8b | `src/vcf_projection.rs`, `src/vcf_semantic_validation.rs` | por registro, componente o portadora, todos los tractos de todas las copias; por ventana, todas las piezas de la copia; por alelo, todos los componentes | los tractos de la secuencia por intervalo; las piezas de la ventana; los componentes del tramo atravesado |
| F6.4.a.8c | `src/configuration.rs`, `src/adapters.rs` | por valor, los saltos de línea desde el principio | los saltos guardados, por búsqueda binaria |
| F6.4.a.8d | `src/results.rs` | por caso y por registro, todos los loci; por locus, todos los registros | los loci por intervalo; las copias por secuencia |
| F6.4.a.8e | `src/results.rs` | por locus, todos los sorteos | los sorteos agrupados por locus |


#### 11.9.3. Huella y mutaciones

**La huella** (`tests/repeat_scale_frozen_identity.rs`, tres pruebas, unos
20 s en depuración) se obtuvo con el binario de `deecb79` compilado en un árbol
extraído con `git archive`, con su propio directorio de compilación; dos
ejecuciones dieron los mismos literales. Con todos los arreglos sigue intacta,
igual que las dos huellas anteriores (`tests/repeat_detection_frozen_identity.rs`
y `tests/pcr_frozen_identity.rs`).

**Mutaciones**, todas sobre una copia extraída con `git archive` más los
ficheros cambiados, con su propio `CARGO_TARGET_DIR`, la copia de seguridad
fuera de git y el SHA-256 del fichero antes, mutado y restaurado (doce
primeras cifras). Las tres pruebas de la huella se abrevian «huella: revisión»
(`the_vcf_fusion_review_scenarios_keep_their_frozen_bytes`), «huella: base» y
«huella: clonal».

| Mutación | Paso | Fichero | SHA antes / mutado / restaurado | Lo que muere |
|---|---|---|---|---|
| Un salto de línea de más al final de los JSON publicados | huella | `src/results.rs` | `7ca2a1aacf06` / `26b15d12d5f7` / `7ca2a1aacf06` | las tres pruebas de la huella |
| Comparar solo con el locus declarado justo antes | F6.4.a.6 | `src/application.rs` | `5e291aed2a10` / `5d9e1aaa131c` / `5e291aed2a10` | las dos pruebas de solapes |
| Índice del tartamudeo que busca el primer tracto por su inicio | F6.4.a.7 | `src/stutter.rs` | `8b1e0ae07559` / `35e7d0f87c93` / `8b1e0ae07559` | la prueba unitaria; huella: base y clonal |
| Índice de los tramos fuente que busca por inicio | F6.4.a.8a | `src/vcf_grouping.rs` | `1eb55c90ccde` / `736967596f34` / `1eb55c90ccde` | la prueba unitaria; huella: revisión |
| Reparto de bases por locus que empieza el rango después de la posición | F6.4.a.8a | `src/vcf_publication.rs` | `b02223885d1c` / `893561e38e01` / `b02223885d1c` | las tres de la huella (la ejecución falla) |
| Tractos realizados de la proyección buscados por inicio | F6.4.a.8b | `src/vcf_projection.rs` | `ddc8b3cefa33` / `5cf500f29625` / `ddc8b3cefa33` | **solo** la prueba unitaria |
| Piezas de una ventana buscadas por inicio | F6.4.a.8b | `src/vcf_projection.rs` | `ddc8b3cefa33` / `8c2c1f36a04d` / `ddc8b3cefa33` | la prueba unitaria; las tres de la huella |
| Componentes atravesados: se salta el que empieza en el inicio del tramo | F6.4.a.8b | `src/vcf_projection.rs` | `a0bf2d57baf3` / `91de2dae05ba` / `a0bf2d57baf3` | huella: revisión; dos de `tests/vcf_overlap_fusion.rs` |
| Tractos del juez semántico buscados por inicio | F6.4.a.8b | `src/vcf_semantic_validation.rs` | `ec9c84948442` / `8f63db48957e` / `ec9c84948442` | **solo** la prueba unitaria |
| La línea cuenta el salto que cae en el desplazamiento (`<=`) | F6.4.a.8c | `src/configuration.rs` | `ab627416eca2` / `8085c754a11f` / `ab627416eca2` | **solo** la prueba unitaria |
| Recuento de bases del VCF que empieza el rango después de la posición | F6.4.a.8d | `src/results.rs` | `7ca2a1aacf06` / `9c1c6597e9cb` / `7ca2a1aacf06` | las tres de la huella (la ejecución falla) |
| El manifiesto pierde el primer sorteo | F6.4.a.8e | `src/results.rs` | `7ca2a1aacf06` / `f12fa86df4a6` / `7ca2a1aacf06` | las tres de la huella |

Lo que esto dice, sin adornos: tres arreglos —el índice de tractos de la
proyección, el del juez semántico y la posición de los errores del TOML—
**solo** los protege su prueba unitaria. En los dos primeros la huella no
puede verlo con los escenarios de hoy: ninguno tiene un tracto realizado que
empiece antes de un registro fundido y lo pise (`tracts_in` lo rechazaría como
invariante rota), y un juez que pierde un tracto es más permisivo, no más
estricto, con un VCF correcto. La posición de un error del TOML no se publica
en ningún fichero. La mutación de los componentes atravesados se
verificó aparte, después de las demás, porque esa función se indexó al final.


#### 11.9.4. Tiempos antes y después

Segundos de pared por etapa, como en 11.4 y 11.5. «Antes» es el binario de
`deecb79` (una pasada); «después», el de este árbol (`rustc 1.97.1`, binario de
pruebas optimizado de SHA-256 `8ce0291e…`), **dos pasadas** separadas por una
barra. «Preflight» es el de la biblioteca; «preflight interno» es el que
`execute` repite por dentro más la copia de entradas y referencia. La lectura
del TOML es `admit_toml_header` más `parse_declared_config`. La memoria
máxima de las dos pasadas difiere en menos de 1 MiB. Los bytes publicados no
se compararon entre pasadas (cada directorio se borró al medirlo, para no
llenar el disco): que no cambian lo vigila la huella. La máquina no tenía
nada más en marcha (comprobado con `ps` al lanzar); la campaña entera corrió
en segundo plano y se esperó a que terminara.

| Loci | Indiv. | Tart. | Binario | Lectura TOML | Preflight | Preflight interno | Verdad | VCF | Evento | Efectiva y manifiesto | Validación final | Total | Memoria máx. (MiB) |
|---:|---:|:---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 2 000 | 10 | no | antes | 0,80 | 1,2 | 1,5 | 7,1 | 6,5 | — | 0,8 | 13,1 | 31,0 | 525 |
| 2 000 | 10 | no | después | 0,01 / 0,01 | 1,5 / 1,2 | 1,4 / 1,6 | 6,9 / 7,3 | 1,2 / 1,4 | — | 0,4 / 0,4 | 9,2 / 9,3 | 20,7 / 21,1 | 526 |
| 4 000 | 10 | no | antes | 3,09 | 1,9 | 2,3 | 12,8 | 27,5 | — | 2,7 | 26,0 | 76,3 | 946 |
| 4 000 | 10 | no | después | 0,03 / 0,03 | 1,9 / 1,9 | 2,1 / 2,2 | 10,9 / 11,9 | 1,7 / 2,0 | — | 0,9 / 1,0 | 18,4 / 16,8 | 35,9 / 35,8 | 946 |
| 8 000 | 10 | no | antes | 13,49 | 4,7 | 5,7 | 23,2 | 207,3 | — | 13,6 | 68,1 | 336,3 | 1 805 |
| 8 000 | 10 | no | después | 0,06 / 0,07 | 4,7 / 4,5 | 4,8 / 4,7 | 21,0 / 22,1 | 2,5 / 2,6 | — | 1,7 / 1,9 | 35,3 / 34,4 | 70,0 / 70,3 | 1 801 |
| 16 000 | 10 | no | después | 0,15 / 0,14 | 17,6 / 15,6 | 16,6 / 16,8 | 48,8 / 48,1 | 4,2 / 4,3 | — | 4,2 / 4,2 | 84,4 / 86,5 | 175,9 / 175,6 | 3 512 |
| 50 284 | 10 | no | después | 0,61 | 201,4 | 205,6 | 290,1 | 21,7 | — | 44,5 | 673,3 | 1437,2 | 10 801 |
| 2 000 | 10 | sí | antes | 0,93 | 1,3 | 2,3 | 7,2 | 6,9 | 139,9 | 1,1 | 15,5 | 175,2 | 525 |
| 2 000 | 10 | sí | después | 0,02 / 0,02 | 1,3 / 1,3 | 1,6 / 1,6 | 6,2 / 6,7 | 1,4 / 1,5 | 74,7 / 82,0 | 0,4 / 0,5 | 8,5 / 8,0 | 94,1 / 101,6 | 526 |
| 4 000 | 10 | sí | antes | 4,84 | 2,7 | 3,1 | 17,0 | 38,8 | 345,7 | 3,4 | 26,3 | 441,9 | 946 |
| 4 000 | 10 | sí | después | 0,03 / 0,04 | 1,9 / 2,0 | 2,2 / 2,4 | 10,6 / 12,0 | 1,7 / 2,8 | 79,7 / 82,1 | 0,9 / 0,9 | 15,8 / 14,9 | 112,7 / 117,1 | 947 |
| 8 000 | 10 | sí | después | 0,06 / 0,06 | 4,8 / 4,6 | 4,9 / 4,9 | 22,4 / 22,6 | 2,7 / 2,6 | 93,4 / 93,7 | 1,8 / 2,9 | 30,5 / 33,3 | 160,4 / 164,5 | 1 806 |
| 1 000 | 100 | no | después | 0,01 / 0,01 | 7,0 / 7,4 | 9,8 / 7,7 | 44,9 / 42,9 | 14,8 / 14,6 | — | 0,6 / 0,5 | 54,2 / 52,5 | 131,3 / 125,7 | 808 |
| 2 000 | 100 | no | después | 0,02 / 0,02 | 9,7 / 9,3 | 11,4 / 10,0 | 69,7 / 73,2 | 18,0 / 19,7 | — | 1,1 / 1,2 | 71,7 / 72,8 | 181,5 / 186,2 | 1 468 |
| 4 000 | 100 | no | después | 0,04 / 0,03 | 16,4 / 16,5 | 18,3 / 17,9 | 123,7 / 134,4 | 29,0 / 30,5 | — | 2,1 / 1,9 | 122,6 / 130,7 | 312,0 / 332,0 | 2 818 |

La ejecución de 50 284 loci —**el catálogo entero del tomate reducido**— se
hizo una vez: 24 min y 10 801 MiB (11,1 GB). Con el código de antes no se
intentó; 11.4 la estimaba en 2 a 3 horas.

Con 100 individuos «antes» no se repitió (4 000 loci tardaban 18 min): las
cifras de 11.4, con un juego de loci algo distinto, eran 147, 282 y 1 081–1 165 s
para 1 000, 2 000 y 4 000 loci; ahora, 128, 184 y 322 s.

**Crecimiento de cada etapa arreglada al doblar los loci** (después, media de
las dos pasadas, 10 individuos sin tartamudeo, de 2 000 a 4 000, a 8 000 y a
16 000; y de 16 000 a 50 284, que es ×3,14):

| Etapa | ×2 | ×2 | ×2 | ×3,14 | Lineal |
|---|---:|---:|---:|---:|---|
| Lectura del TOML | ×3,0 | ×2,2 | ×2,2 | ×4,2 | sí (céntimas de segundo) |
| Resolución (solapes) | ×2,0 | ×2,0 | ×2,4 | ×6,5 | sí hasta 16 000; 0,06 s con 50 284 |
| VCF | ×1,4 | ×1,4 | ×1,7 | ×5,1 | sí hasta 16 000; algo más con 50 284 |
| Evento con tartamudeo (10 individuos, 0,1×) | ×1,03 | ×1,16 | — | — | no depende de los loci |
| Efectiva y manifiesto | ×2,4 | ×1,9 | ×2,3 | ×10,6 | sí hasta 16 000; **no** con 50 284 |
| Validación final | ×1,9 | ×2,0 | ×2,5 | ×7,9 | **no**: dentro va el punto 5 (11.9.6) |
| Verdad | ×1,6 | ×1,9 | ×2,2 | ×6,0 | **no** a partir de 16 000 |
| Preflight | ×1,4 | ×2,4 | ×3,6 | ×12,1 | **no**: es el punto 5 (11.9.6) |

Lo que esta tabla dice, sin adornos:

- **Las cuatro etapas que la pausa mandó arreglar lo están** hasta 16 000 loci:
  la lectura del TOML pasa de 13,5 s a 0,06 s con 8 000 loci; el VCF, de
  207 s a 2,5 s; el evento con tartamudeo deja de depender de los loci (de
  346 s a 81 s con 4 000); la validación final, de 68 s a 35 s, y lo que en
  ella queda superlineal es lo que la pausa dejó fuera.
- **Con el catálogo entero aparecen tres crecimientos más que lineales que no
  estaban en 11.5**: la verdad (×6,0 para ×3,14 loci), la etapa de la
  efectiva y el manifiesto (×10,6) y, en menor medida, el VCF (×5,1). De la
  verdad, la lectura del código apunta al punto 5 de 11.5
  (`src/truth_projection.rs:1177`, que busca el alelo de cada caso
  recorriendo las piezas de la copia), pero **no se midió por función**. De
  la efectiva y el manifiesto no sé la causa: hasta 16 000 loci crece en
  línea recta, el arreglo de los sorteos (F6.4.a.8e) está y la etapa escribe
  y relee 480 MB de JSON; puede ser memoria (11 GB residentes) u otro
  recorrido. Pesan 290 s y 44 s de 1 437.
- **Lo que domina el catálogo entero es el punto 5**: rehacer las individuas
  iniciales se hace en el preflight, en el preflight interno y en el juez de
  las estadísticas de la población inicial. Medido por función en una copia
  instrumentada con 8 000 y 16 000 loci, `preflight_initial_individuals`
  pasa de 12,8 a 58,9 s (dos llamadas, ×4,6) y el juez de 12,5 a 39,7 s
  (×3,2); el resto de la validación final, de 35 a 62 s (×1,8). Con 50 284
  loci el preflight solo ya son 201 s, y sumando las tres veces es del orden
  de 800 s de los 1 437.


#### 11.9.5. Extrapolación al catálogo del tomate reducido

Con los arreglos de esta entrega, y sin tocar el punto 5:

| Catálogo de 50 284 loci | Antes (11.4, extrapolado) | Ahora |
|---|---|---|
| 10 individuos, sin tartamudeo | 2 a 3 horas; unos 11 GB | **medido: 1 437 s (24 min), 10 801 MiB** |
| 10 individuos, con evento técnico a 0,1× | — | unos 1 530 s: lo anterior más los 80–95 s del evento, que no dependen de los loci (extrapolado) |
| 100 individuos, sin tartamudeo | unos dos días; unos 34 GB | **unas 3,6 a 3,8 horas** y unos **33 800 MiB** (extrapolado) |

Cómo se extrapola, para que se pueda rehacer:

- **Tiempo con 100 individuos**, dos cuentas que coinciden: (a) un ajuste
  `a + b·n + c·n²` por mínimos cuadrados a la media de las dos pasadas con
  1 000, 2 000 y 4 000 loci da 13 743 s; (b) la razón entre 100 y 10
  individuos con 4 000 loci (322 / 35,9 = 9,0) aplicada a los 1 437 s
  medidos da 12 907 s. El término cuadrático es el punto 5, que crece además
  con los individuos. El mismo ajuste con 10 individuos (2 000 a 16 000
  loci) predecía 1 138 s y se midieron 1 437: **la extrapolación se queda
  corta**, porque con el catálogo entero aparecen los crecimientos del
  apartado anterior; la de 100 individuos puede quedarse corta igual.
- **Memoria**: en línea recta, como en 11.5 (0,213 MiB por locus con 10
  individuos, entre 2 000 y 16 000; 0,670 con 100, entre 1 000 y 4 000). La
  de 10 individuos predecía 10 824 MiB y se midieron 10 801. **Con 100
  individuos no cabe en esta máquina** (27 GB): eso no lo arregla esta
  entrega, sino el catálogo aparte de la entrega 4 y, si hace falta, bajar
  la memoria por locus (DP-SSR-18 en 11.7).


#### 11.9.6. Lo que queda sin arreglar

1. **Rehacer las individuas iniciales** (punto 5 de 11.5:
   `apply_repeat_alleles` busca el tramo intacto recorriendo las piezas de la
   copia y `with_current_repeat_allele` copia la lista entera de piezas por
   cada locus). Se hace tres veces por ejecución —preflight, preflight interno
   y juez de las estadísticas de la población inicial— y crece ×3 a ×4,6 al
   doblar loci. La decisión de la pausa lo deja fuera y **no se tocó**; con
   el catálogo entero es del orden de 800 s de los 1 437 medidos: lo que más
   pesa de lo que no es lineal.
2. **En el juez semántico del VCF**, que solo corre en las pruebas, quedan
   dos recorridos cuadráticos por copia: al normalizar, cada alelo se compara
   con **todos** los componentes de la secuencia para decidir qué funde
   (`normalize_oracle_components`), y cada ventana recorre todas las piezas
   de la copia (`OraclePath::window_from`). Se arregló lo que 11.5 nombraba
   (`:726`, `:1449`, `:1566`) y el límite de alineación; estos dos no, porque
   el juez no corre en la ejecución del usuario y la condición de fundido mira
   la posición del alelo de **otro** componente, que la alineación puede
   llevar lejos de su componente: acotarlo por rango exige una regla nueva
   que el juez no tiene, no un índice. Con las pruebas de hoy no pesa.
3. **`copy_is_reference_at`** (`src/vcf_grouping.rs`) recorre todas las
   ventanas iguales a la referencia por cada caso **sin** `ALT`. Es cuadrático
   en teoría, pero en las mediciones no se llamó ni una vez.
4. **La memoria** sigue creciendo en línea recta, 0,21 MiB por locus con
   10 individuos: los árboles JSON de la configuración efectiva y del
   manifiesto (apartado 11.5). No era objeto de esta entrega; el catálogo
   aparte de la entrega 4 la baja.
5. **El evento técnico** cuesta unos 80 a 95 s con 10 individuos a 0,1× y no
   depende del número de loci: es el muestreo de las lecturas del genoma
   entero.

