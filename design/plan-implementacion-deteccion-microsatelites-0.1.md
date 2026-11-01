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

**Reescrita el 01/11/2026** (decisión del responsable, opción A; ver el
apartado 11.12). La regla original, opción (a), recorría los tramos por
prioridad y descartaba **entero** cualquier aspirante que pisara a uno ya
aceptado, aunque compartiera una sola base. Eso perdía microsatélites reales
pegados a otro: un microsatélite suele acabar con una copia parcial del motivo
que coincide con el principio de lo que viene detrás, y el buscador entrega
los dos tramos pisándose esa copia (`(CA)×10 A×12`; el `(TCTG)` central de
`(TCTA)×8 (TCTG)×8 (TCTA)×8`). Es DT-061, parte (a).

| Opción | Consecuencia |
|---|---|
| (a) Gana uno, por regla fija; los demás se descartan enteros (la original) | Determinista; pierde el microsatélite que solo pisa la copia parcial de la frontera |
| **(a′) La misma prioridad, pero el aspirante que pisa una copia de su motivo o menos por cada lado se recorta y se vuelve a juzgar** | Conserva los dos microsatélites pegados; el resultado sigue sin depender del orden de entrada; el recortado pasa por los mismos criterios que cualquier tramo |
| (b) Se descartan todos los solapados | Pierde el locus bueno con el malo |
| (c) Se funden en un locus compuesto | Exige derivar bloques, que S-30 prohíbe |

**Regla vigente: (a′).** Se recorren los candidatos por prioridad: mayor
puntuación del detector; empate, más bases; empate, periodo menor; empate, el
de la izquierda. De cada aspirante:

1. **Sin solape** con los ya aceptados: se acepta.
2. **Solape de una copia de su propio motivo o menos por cada lado** (solape ≤
   periodo del aspirante), sin ningún aceptado dentro de él ni que lo cubra
   entero: **se recorta** a la parte que no solapa y ese tramo **vuelve a
   pasar por todos los criterios** (`classify_interval`: suelo de copias,
   pureza, tirada, tope y base A/C/G/T a cada lado). Entra solo si los supera.
3. **Cualquier otro caso** (más de una copia de solape por algún lado, un
   aceptado dentro del aspirante o que lo cubre): se descarta entero, como
   antes.

Decisiones de detalle, con su porqué:

- **Puntuación del recortado:** se recalcula con las comparaciones de su
  propio tramo (+2 y −7), como el recorte del inicio del buscador, de modo que
  todo candidato lleva la suma de sus propias comparaciones. No decide nada
  en la resolución, porque el recortado se acepta en el acto y no vuelve a
  competir.
- **Motivo del recortado:** se lee al inicio del tramo recortado (DP-SSR-03),
  así que tras recortar el inicio sale una **rotación**: el `(TCTG)` central
  del compuesto sale como `GTCT`, y el identificador lleva ese motivo
  (`SSR-1-2736-2767-GTCT`). La coincidencia de un declarado (DP-SSR-09)
  compara rotaciones y no se ve afectada.
- **Dónde se cuenta:** un recortado que entra se cuenta como aceptado en su
  longitud de motivo; **uno que cae** (porque ya no llega al suelo, a la
  pureza, etc.) se cuenta en **solape**, no en la etapa del criterio que ya no
  supera: superó todos por sí solo y lo pierde el solape, de modo que las
  columnas de suelo, pureza y demás siguen contando solo lo que cae por sí
  solo y el balance detectados = descartados + aceptados se mantiene.
- **Regiones:** se aplican **antes** de los solapes (DP-SSR-18). El recortado
  es un subconjunto del aspirante, que ya cabía en su región, así que no
  necesita volver a comprobarlas; y un tramo que no cabe en su región no
  recorta a nadie.
- **Motivo no primitivo tras recortar:** comprobado por fuerza bruta
  (1 800 000 secuencias de dos a cuatro bloques con mutaciones, 665 000
  recortados aceptados, con pureza mínima de 0 a 0,5): ninguno. Con la pureza
  por omisión es inalcanzable en la práctica, y no se añade guarda.

**El umbral es «una copia o menos», no «menos de una copia».** Con «menos de
una copia» una tira de una sola base (periodo 1) no podría recortarse nunca,
porque un solape de una base es ya una copia, y `(CA)×10 A×12` seguiría
perdiendo su tira de `A`. El responsable fijó «una copia o menos» el
01/11/2026. Consecuencia: `(CA)×10 G (CAG)×7`, donde el `(CAG)×8` pisa dos
bases al `(CA)×10` (una copia), da ahora los dos, el `(CA)×10` recortado a
`(CA)×9`.

**Recomendación original (a), sustituida por (a′).** Dos tramos que solo se
tocan (uno acaba donde empieza el otro) no se solapan y se conservan, como hoy
con los declarados.

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
peso medido) **no** se tocaba; **el 29/10/2026, con las mediciones de 11.9
delante, el responsable decide arreglarlo también** (pasos F6.4.a.8f a
F6.4.a.8j, apartado 11.10). Cada arreglo es de **rendimiento**: no puede
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
- ☐ **F6.4.a.8f.** Alelos de repetición: comprobar y aplicar los de cada copia
  en un lote, sin recorrer ni copiar sus fragmentos por cada locus.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 5, **medido** en 11.10: con 50 284 loci
    es el coste mayor (546 s de 973 en la copia instrumentada), y se paga
    cuatro veces por ejecución —los dos preflights y las dos llamadas del
    juez de las estadísticas de la población inicial—.
  - Cambia: `src/segments.rs` (`RepeatAlleleBatch`: los fragmentos de fuente
    de la copia ordenados por inicio una vez, el de cada locus por búsqueda
    binaria, los cambios de base del fragmento también por búsqueda binaria,
    y la cadena final construida en una sola pasada; si la cadena no admite
    el índice —fragmentos desordenados o que se pisan, que ninguna operación
    del módulo produce— o dos alelos encolados se pisan, hace lo de antes,
    alelo a alelo) y `src/repeat_genotypes.rs` (`apply_repeat_alleles`: la
    comprobación de cada locus y copia en el orden de siempre, contra la
    cadena original, y un lote por copia aplicado al final).
  - Se comprueba con: las pruebas unitarias
    `the_repeat_allele_batch_matches_applying_one_allele_at_a_time` (misma
    cadena y misma respuesta a «¿llega intacto?» que uno a uno, con loci que
    empiezan justo donde empieza un fragmento y un cambio de base en la
    primera base del intervalo) y
    `the_repeat_allele_batch_falls_back_to_one_at_a_time_when_it_cannot_index`;
    las dos huellas intactas.
  - **Debe fallar contra:** un índice que busque el fragmento con `<` en vez
    de `<=` y pierda el que empieza justo en el inicio del locus.
- ☐ **F6.4.a.8g.** Verdad: el alelo de cada caso y la huella de cada
  procedencia, sin recorrer listas enteras.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 5 (`src/truth_projection.rs:1177`) y,
    **corregido por la medición** (11.10), la comprobación de huellas: el
    crecimiento ×6,0 de la verdad no estaba en el alelo de cada caso sino en
    `require_digest`, que recorría por cada fragmento y cada evento de cada
    copia la lista de procedencias, que tiene una por locus aplicado.
  - Cambia: `src/truth_projection.rs` (`project_repeat_genotypes`: una tabla
    por copia de sus alelos por intervalo, rehecha al cambiar de secuencia y
    con el primer fragmento de cada intervalo; `project_truth` y los que
    proyectan cada copia: las procedencias con huella en un conjunto).
  - Se comprueba con: las dos huellas intactas (dos secuencias con loci, de
    modo que la tabla tiene que rehacerse).
  - **Debe fallar contra:** reutilizar la tabla de la secuencia anterior al
    cambiar de secuencia.
- ☐ **F6.4.a.8h.** Admisión de los loci contra las variantes pequeñas por
  índices de posición.
  - Depende de: la huella de entrada.
  - Procede de: apartado 11.5, punto 5 (`src/execution.rs:2725`). Con la
    matriz de una fila de la medición no pesa; con una matriz grande crece
    con loci × filas.
  - Cambia: `src/execution.rs` (`admit_repeat_loci_against_small_variants`:
    los cambios de base directos y las filas de la matriz en mapas por
    secuencia y posición, con el primero en su orden; los identificadores de
    los cambios directos en otro; por cada locus se nombra el primero en su
    orden de los que caen dentro, y entre un cambio llamado como el locus y
    uno que lo invade, el declarado antes).
  - Se comprueba con: `tests/repeat_scale_regressions.rs`,
    `the_small_variant_admission_names_the_same_invader_as_comparing_everything`
    (primera y última base de un locus, las mismas coordenadas en la otra
    secuencia, dos cambios en el mismo locus, el cambio llamado como el locus
    en los dos órdenes, cambio contra fila y locus contra locus), que pasa
    igual con el código de antes; `tests/repeat_locus_invasion.rs`.
  - **Debe fallar contra:** un rango que empiece en la base siguiente al
    inicio del locus y deje pasar la fila que cae en la primera.
- ☐ **F6.4.a.8i.** Configuración efectiva: las bases del VCF por locus en un
  mapa.
  - Depende de: la huella de entrada.
  - Procede de: la medición (11.10): el crecimiento ×10,6 de la etapa
    «configuración efectiva y manifiesto» que 11.9 dejó sin explicar es
    `publish_measured_vcf_tract_bases`, que buscaba las bases de cada locus
    recorriendo la lista de todos (16 s con 50 284 loci).
  - Cambia: `src/results.rs` (`publish_measured_vcf_tract_bases`).
  - Se comprueba con: la huella intacta (publica `vcf_tract_bases_bound` de
    cada locus en la configuración efectiva).
  - **Debe fallar contra:** perder la primera entrada de la lista al hacer el
    mapa.
- ☐ **F6.4.a.8j.** Juez de la verdad: el alelo de cada caso de repetición sin
  recorrer el historial entero de cada copia ni todas las copias.
  - Depende de: la huella de entrada.
  - Procede de: la medición de los pasos anteriores (11.10): arreglados
    F6.4.a.8f a F6.4.a.8i, lo único que seguía creciendo con el cuadrado era
    `validate_repeat_genotype` (`src/truth_validation.rs`), que por cada caso
    de repetición recorría todas las copias de todos los individuos y, en
    cada una de las suyas, el historial entero, con una entrada por locus
    aplicado (×9,9 al pasar de 16 000 a 50 284 loci; 31 s). Es la misma
    clase de coste; el juez corre en cada ejecución.
  - Cambia: `src/truth_validation.rs` (`CopySummary` guarda las posiciones de
    los alelos de repetición del historial por intervalo, en su orden;
    `validate_repeat_genotype` toma el tramo del mapa de copias de ese
    individuo y esa secuencia y busca en la tabla, conservando el primero que
    cumple). No añade ni quita comprobaciones: la puerta de cobertura de los
    jueces no cambia.
  - Se comprueba con: las dos huellas intactas;
    `tests/truth_validation_guards.rs` y `tests/truth_validation.rs`.
  - **Debe fallar contra:** un tramo de copias que acabe antes de la última
    copia del individuo.

> **Revisión de la entrega 2:** cada etapa arreglada crece en línea recta al
> doblar los loci, con la huella intacta (apartado 11.9).
>
> **Estado, 29/10/2026:** los siete pasos están implementados, con las
> puertas en verde y cada «debe fallar contra» comprobado mutando una copia
> aparte (apartado 11.9.3). Falta la revisión adversarial; por eso siguen sin
> marcar.
>
> **Estado, 29/10/2026 (noche):** los cinco pasos añadidos, F6.4.a.8f a
> F6.4.a.8j, están implementados, con las puertas en verde, las huellas
> intactas y cada «debe fallar contra» comprobado mutando una copia aparte
> (apartado 11.10). Sin revisión adversarial y sin *commit*.

### Entrega 3 — detección en el producto

Reescrita el 29/10/2026 con las decisiones del responsable de la sección 3 y
del apartado 11.7. Lo que cambia respecto de la versión anterior:

- **Modo único (DP-SSR-01).** La detección se activa por la presencia de la
  sección `[reference.repeat_detection]`, que declara el deslizamiento común.
  Sin la sección **no hay microsatélites**: un locus declarado, un tracto
  insertado o uno plantado sin ella se rechazan con un diagnóstico que remite
  a la sección. Se retira el régimen de loci declarados sin detección.
- **Tres etapas siempre:** construir la referencia → detectar sobre la
  referencia **final** → variar.
- **Criterios de 11.7:** copias completas con suelo `[8, 6, 6, 6, 5, 4]`
  configurable; pureza 0,8 y tirada 3, los de la admisión
  (`[reference.repeat_admission]`); tope de 200 unidades configurable; una base
  A/C/G/T a cada lado; la regla de solapes de DP-SSR-05; solo loci de un motivo.
- **Regiones (DP-SSR-18)**, paso nuevo F6.4.a.10b.
- **Esquemas en su versión actual** (nota alfa de DTC-036): se retocan in situ
  el general v4, la efectiva v6 y el manifiesto v6. **No** se crean versiones.
- **Publicación mínima:** el catálogo aparte es de la entrega 4. Aquí cada
  locus de `repeat_loci` (efectiva y manifiesto) gana su `origin`
  (`detected`/`declared`), y la efectiva una sección `repeat_detection` con
  los criterios, el deslizamiento común, las regiones y los tramos retirados.
  Nada más.

Los números F6.4.a.9 a F6.4.a.13 se conservan; lo añadido lleva sufijo.

- ☑ **F6.4.a.9.** Sección `[reference.repeat_detection]` y su resolución.
  - Depende de: pausa de la entrega 1; F6.4.a.1.
  - Procede de: REQ-CONF-142; REQ-CONF-139; decisiones 5 y 6; DP-SSR-01, 02,
    07, 15 (en la forma del 28/10: retocar in situ).
  - Cambia: `src/configuration.rs` (`DeclaredRepeatDetection`),
    `src/adapters.rs` (la lectura; un locus sin alelos ni generador se admite
    con la sección, porque hereda el común), `src/application.rs`
    (`EffectiveRepeatDetection`; un locus declarado sin generador ni alelos
    recibe el deslizamiento común; el rechazo de loci sin la sección),
    `schemas/config-general-v4.schema.json` **en su versión actual**.
  - Forma: `min_copies` (seis enteros positivos, opcional), `max_units`
    (entero positivo, opcional, 200 por omisión), la tabla obligatoria
    `[reference.repeat_detection.slippage]` con `model` y `mean_events`
    obligatorios y `expansion_bias` y `step_units` como en REQ-CONF-139, y la
    lista opcional `[[reference.repeat_detection.regions]]` del paso F6.4.a.10b.
  - Se comprueba con: `tests/repeat_detection_config.rs`: sin la sección, la
    huella de F6.4.a.2 intacta y un locus declarado rechazado nombrando la
    sección; con ella, la falta de `slippage`, de `model` o de `mean_events` se
    rechaza nombrando el campo y con un ejemplo en el mensaje; los suelos por
    omisión aparecen en la efectiva aunque no se declaren; un suelo 0, una
    tabla de cinco valores y un tope 0 se rechazan; un declarado sin alelos ni
    generador hereda el común. **Debe fallar contra:** un valor por omisión de
    `mean_events` (la decisión 6 lo prohíbe).
- ☑ **F6.4.a.10.** La detección corre en su sitio.
  - Depende de: F6.4.a.4, F6.4.a.9.
  - Procede de: decisión 1; REQ-REF-38; DTC-103; DP-SSR-10.
  - Cambia: `src/application.rs` (`EffectiveConfig::with_detected_repeat_loci`,
    hermana de `with_generated_repeat_alleles`: barre cada secuencia de la
    referencia final, convierte cada candidato en locus con origen
    `detected`, identificador `SSR-<índice de secuencia en base uno>-<inicio>-<fin>-<motivo>`
    en base uno, el deslizamiento común, el tope de la sección y los mínimos
    generales; la lista final queda **ordenada por secuencia y por inicio**),
    `src/execution.rs` (detectar antes de preparar la matriz y de generar
    abanicos, en `validate` y en `run`; `requires_materialized_bases` cuenta la
    sección) y `src/results.rs` (la reconstrucción del juez físico, que rehace
    la configuración desde el TOML copiado, detecta igual sobre la referencia
    publicada; el juez independiente sin detector es de la entrega 4). Las
    comprobaciones del tartamudeo y de la PCR contra los loci (REQ-CONF-132,
    REQ-CONF-138) se aplican también a los detectados.
  - Se comprueba con: `tests/repeat_detection_run.rs`: `validate` y `run` ven
    el mismo catálogo; una referencia suministrada con tractos plantados en
    posiciones **escritas en la prueba** da exactamente esos loci (riesgo R-1:
    el oráculo no llama al detector), con **dos semillas** (R-5); con
    inserción, el tracto insertado **antes** de un tracto existente se detecta
    en coordenadas finales y el existente sale desplazado la longitud del
    insertado (R-6); con referencia simulada, el plantado; la `N` suelta
    dentro de un tracto largo (apartado 11.6 bis) da un solo locus, con la `N`
    contada como discrepancia de pureza, y la ejecución publica su VCF.
    **Debe fallar contra:** detectar sobre la referencia de partida en un
    escenario con inserción.
- ☑ **F6.4.a.10b.** Regiones de detección (DP-SSR-18, añadido el 29/10/2026).
  - Depende de: F6.4.a.10.
  - Procede de: DP-SSR-18 y apartado 11.7.
  - Cambia: `src/adapters.rs`, `src/application.rs` (cada región en base uno,
    cerrada, en coordenadas **finales**, dentro de su secuencia; dos regiones
    de la misma secuencia no se solapan; una lista vacía se rechaza) y
    `src/repeat_detection.rs` (un tramo se conserva solo si cabe **entero** en
    una región; el descarte se cuenta en una etapa propia, `outside_regions`,
    que se aplica **antes** de los solapes, de modo que lo de fuera no puede
    quitarle sitio a lo de dentro).
  - Se comprueba con: dos tractos plantados, una región que contiene uno: sale
    solo ese; una región que corta un tracto por la mitad no lo detecta; sin
    lista, los dos; una región con coordenadas de partida en un escenario con
    inserción se rechaza o no contiene el tracto (las coordenadas son
    finales). **Debe fallar contra:** conservar un tramo que solo **solapa** la
    región.
- ☑ **F6.4.a.11.** Coincidencia de declarados, insertados y plantados.
  - Depende de: F6.4.a.10.
  - Procede de: decisión 3; REQ-REF-40; DP-SSR-09 (aceptado: conserva sus
    alelos escritos a mano).
  - Cambia: `src/application.rs`. Un declarado coincide si su intervalo final
    **se solapa** con un locus detectado (después de criterios, regiones y
    solapes) de la misma secuencia cuyo motivo sea una **rotación** de uno de
    los suyos; mandan sus fronteras y **sustituye** a los detectados que
    solapa; conserva su identificador y sus alelos. Un declarado cuyos motivos
    tienen **todos** 7 pb o más queda exento (REQ-CONF-60) y también sustituye
    a lo que solapa. Sin coincidencia se rechaza nombrando el locus, el
    detectado más cercano de su secuencia y la **etapa** en que cayó el tramo
    de su motivo que lo solapaba (suelo, pureza, tirada, tope, borde, región o
    solape). Un declarado cuyo identificador repite el de un detectado que
    queda se rechaza como identificador repetido.
  - Se comprueba con: `tests/repeat_detection_run.rs`: fronteras
    distintas pero solapadas coinciden y mandan las declaradas; `AC` sobre un
    detectado `CA` coincide; uno sin tracto se rechaza nombrando el locus y lo
    más cercano; uno por debajo del suelo nombra la etapa `min_copies`; uno de
    7 pb pasa sin detector; un compuesto declarado sustituye a los dos tramos
    que solapa. **Debe fallar contra:** emparejar por igualdad de cadena del
    motivo; emparejar por igualdad exacta de fronteras.
- ☑ **F6.4.a.12.** Variantes pequeñas frente a tractos detectados.
  - Depende de: F6.4.a.10.
  - Procede de: REQ-CONF-143; DP-SSR-11 (provisional, con la tarea de DTF-066);
    S-32.
  - Cambia: `src/population_snvs.rs` (el generador recibe los intervalos
    excluidos: una posición dentro de un locus del catálogo tiene peso cero y
    no consume azar; sin intervalos, el recorrido es el de siempre),
    `src/application.rs` (`EffectiveConfig::with_retired_repeat_loci`: un
    cambio de base directo o una fila de la matriz —suministrada o
    reutilizada— dentro de un locus **detectado** lo retira del catálogo, y la
    efectiva publica cuál y por qué en `repeat_detection.retired_loci`),
    `src/execution.rs` (orden: detectar → preparar la matriz excluyendo el
    catálogo → retirar → generar abanicos → resto del preflight; un SNV dentro
    de un **declarado** se sigue rechazando en
    `admit_repeat_loci_against_small_variants`) y la reconstrucción del juez
    físico en `src/results.rs` en el mismo orden.
  - Se comprueba con: `tests/repeat_detection_run.rs`: matriz generada de
    1000 SNV sobre una referencia con tractos plantados que cubren alrededor
    del 0,7 %, aceptada con **dos semillas** y ninguna fila dentro de un locus;
    fila suministrada dentro de un detectado: el tracto se retira y la causa
    se publica; la misma fila dentro de un **declarado**: rechazo como hoy.
    **Debe fallar contra:** un generador que excluya solo el interior y no la
    primera base del tracto; una exclusión que cambie el consumo de azar sin
    detección (la huella de F6.4.a.2 lo caza).
- ☑ **F6.4.a.13.** Diagnósticos de dos pasos.
  - Depende de: F6.4.a.9.
  - Procede de: decisión 2; REQ-REF-37; DP-SSR-12.
  - Cambia: `src/application.rs` (el texto del rechazo de la inserción junto a
    la reducción, y el rechazo nuevo de la sección de detección junto a una
    reducción en la misma ejecución).
  - Se comprueba con: `tests/reference_tract_insertion_run.rs` actualizado y
    `tests/repeat_detection_config.rs`: el mensaje nombra las dos ejecuciones y
    el campo `fasta` de la segunda. **Debe fallar contra:** el texto de hoy.
- ☑ **F6.4.a.13b.** Reescritura de las pruebas del régimen retirado y de la
  huella de escala.
  - Depende de: F6.4.a.9 a F6.4.a.12.
  - Procede de: DP-SSR-01 («las pruebas que lo usan se reescriben»).
  - Cambia: las pruebas del inventario de abajo, que declaran loci sin la
    sección; cada una gana la sección con **regiones que cubren exactamente
    sus loci declarados** (y, si un locus suyo queda por debajo del suelo por
    omisión, un `min_copies` propio del escenario), de modo que lo detectado
    coincide con lo declarado y la prueba sigue comprobando lo mismo.
    `tests/repeat_scale_frozen_identity.rs` se reescribe así y **no** se
    recongela a ciegas: los literales siguen siendo los obtenidos de
    `dcae406`; solo pueden cambiar la configuración efectiva, el manifiesto y
    la copia del TOML, en los campos nuevos, y lo que dependa de la huella del
    TOML (el identificador de ejecución). La prueba escribe qué cambió y por
    qué.
  - Se comprueba con: la huella de F6.4.a.2 sin tocar un literal; la de escala
    con la demostración anterior; las puertas.
  - **Inventario** (29/10/2026, `grep -l repeat_loci tests/`): 
    `clonal_composition`, `clonal_independent_validator`, `clonal_lineage`,
    `event_manifest_guards`, `final_pairs_against_truth` y
    `final_pairs_oracle/mod.rs`, `manifest_minimal`, `pcr_admission`,
    `publication_process`, `published_resource_validation`,
    `reference_tract_insertion_run`, `repeat_admission_guard`,
    `repeat_genotype_draw`, `repeat_locus_contract`, `repeat_locus_invasion`,
    `repeat_manifest`, `repeat_scale_frozen_identity`,
    `repeat_scale_regressions`, `repeat_slippage`, `repeat_tract_budget`,
    `repeat_truth_v3`, `repeat_vcf_projection`, `results_effective_guards`,
    `results_manifest_guards`, `results_population_snvs_guards`,
    `results_repeat_guards`, `simulated_reference_repeat_tract`,
    `stutter_admission`, `stutter_independent_validator`, `stutter_pipeline`,
    `stutter_tract_length`, `truth_guard_scenario/mod.rs`,
    `vcf_null_effect_case`, `vcf_overlap_fusion`, y los TOML de
    `tests/fixtures/vcf-fusion/review/*` (quince de veintiuno),
    `tests/fixtures/vcf-fusion/tract-without-indels` y
    `tests/fixtures/vcf-null-effect/direct-cancel`. Las pruebas de la interfaz
    (`tests_ui/test_repeat_loci_ui.py`, `tests_ui/test_stutter_ui.py`) **no**
    se tocan aquí: son de la entrega 5.

> **Revisión de la entrega 3:** la detección se ejecuta en las tres clases de
> referencia, coincide o rechaza los declarados, y no deja rechazos que
> dependan de la semilla. Revisión adversarial de la entrega.
>
> **Estado, 29/10/2026: detenida en F6.4.a.13b, pendiente de una decisión.**
> F6.4.a.9 a F6.4.a.13 están implementados con sus pruebas
> (`tests/repeat_detection_config.rs`, `tests/repeat_detection_run.rs`), y cada
> «debe fallar contra» se comprobó mutando una copia aparte. La reescritura de
> las pruebas del régimen retirado se detuvo porque **el detector funde en un
> solo tramo dos tractos contiguos que comparten periodo**, y entonces un
> locus declarado no encuentra tramo con el que coincidir: `(GA)×12` pegado a
> `C×12` sale como un tramo de periodo 2 de pureza 24/36, que se descarta, y
> la `C×12` aparte (escenario 1 de `tests/repeat_scale_frozen_identity.rs`);
> `G×60` pegado a `A×15` sale como un tramo de periodo 1 con motivo `G`, y la
> `A×15` desaparece (la evidencia de DT-009 y `tests/clonal_composition.rs`);
> `G×50 (GT)×25 (CA)×50 G×100` da un solo tramo de periodo 2 que no alcanza la
> pureza, y el `(CA)×50` desaparece (`tests/repeat_slippage.rs` y otras); un
> compuesto declarado `(TCTA)(TCTG)(TCTA)` sale como un solo tramo de periodo
> 4 que no alcanza la pureza y no puede coincidir. Lo segundo deja sin
> respuesta en las normas a los compuestos declarados (DP-SSR-06 y DP-SSR-09
> dicen que se declaran a mano y coinciden con lo detectado). La huella de
> escala no puede reproducirse sin cambiar su referencia, y no se ha
> recongelado.
>
> **Estado, 30/10/2026: se reanuda con el núcleo de `pop_var_caller`
> `fa215c23`** (apartado 11.11), que ya no funde tramos pegados de motivos
> distintos. Un compuesto declarado `(TCTA)(TCTG)(TCTA)` coincide con los dos
> `(TCTA)` detectados (se solapa con ellos y su motivo es una rotación de uno
> de los suyos) y los sustituye, aunque el `(TCTG)` central no se detecte. Las
> pérdidas que quedan en la frontera son en parte de nuestra regla de solapes
> y en parte del buscador original (DT-061, corregida el 31/10/2026). Decisiones menores del responsable del 30/10/2026, las cuatro
> según la recomendación:
>
> 1. **Sin individuos no se generan alelos** de los loci detectados.
>    Confirmada por el responsable el 30/10/2026: los detectados se dejan
>    fuera del catálogo; un declarado sin alelos en la misma situación se
>    sigue rechazando, y la diferencia se acepta.
> 2. **Con PCR, el tartamudeo solo se exige si hay al menos un locus
>    detectado**; sin ninguno, la PCR se admite sin tartamudeo.
> 3. **Una `N` suelta dentro de un tracto cuenta como fallo de pureza**, no
>    corta el tracto.
> 4. **Un tramo solo se detecta si cabe entero en una región**; uno que solo
>    la solapa se descarta en la etapa `outside_regions`.
>
> **Estado, 30/10/2026 (noche): entrega 3 terminada, pendiente de la revisión
> adversarial.** Sin *commit*.
>
> - **Las cuatro decisiones menores.** La 3 (la `N` cuenta como fallo de
>   pureza) y la 4 (un tramo solo cuenta si cabe entero en una región) ya
>   estaban implementadas y probadas; se comprobó que sus pruebas mueren con la
>   variante contraria. La 1 y la 2 se implementaron: sin copias haplotípicas
>   los detectados **no entran en el catálogo** (no hay alelos que generar ni
>   otros que darles) y la ejecución se publica; la PCR sin tartamudeo de base
>   solo se rechaza al resolver la configuración si hay loci declarados, y
>   después de detectar y de retirar si el catálogo final tiene algún locus.
>   El tartamudeo clonal se conecta con la sección al resolver la
>   configuración y **se desconecta si el catálogo final se queda vacío**
>   (ronda 2 de la revisión, abajo).
> - **F6.4.a.13b.** Las pruebas del régimen retirado ganan la sección con
>   `tests/support/detection_section.rs`: una región por locus con seis bases
>   de margen y suelos que son la mitad de las copias del locus más corto de
>   cada longitud de motivo (inalcanzables en las longitudes que nadie
>   declara). Se comprobó ejecutando la batería entera con un registro de cada
>   locus detectado que queda en el catálogo que, fuera de las pruebas de la
>   detección, ninguna gana un locus que no declare. Cuatro grupos de pruebas
>   necesitaron algo más, y se dice en cada una: las que rodeaban su tracto de
>   tiras largas de `G` usan un fondo de `G` y `T` (el detector no ve un
>   microsatélite corto entre dos tiras largas de una misma base; anotado en
>   DT-061); las que declaran un locus que es parte de una repetición más
>   larga llevan su propia región; las que suponían el orden del TOML usan el
>   del catálogo (por secuencia y por inicio); y tres cambian lo que comprueban
>   porque el régimen lo cambia —una región sin repetición ya se rechaza por
>   la coincidencia (REQ-REF-40) y no por la guarda de admisión, las variantes
>   generadas ya no pueden caer dentro de un tracto (REQ-CONF-143), y el VCF de
>   `61b9177` se compara sin la línea del identificador de ejecución—.
> - **La huella de escala** (`tests/repeat_scale_frozen_identity.rs`) tiene
>   sus literales: la tabla normalizada sale de los árboles que publicó el
>   binario de `dcae406`, compilado aparte, y la verdad, el VCF y los FASTQ de
>   los 23 escenarios son idénticos a los de `dcae406` tras la normalización.
>   Solo cambian la efectiva (`repeat_detection` y el `origin` de cada locus),
>   el manifiesto (el `origin`, el identificador de ejecución y las huellas y
>   tamaños de lo que depende del TOML) y la copia del TOML.
> - **Las huellas de la evidencia de DT-009** en
>   `tests/technical_catalog_contract.rs` y `tests/clonal_homopolymer_fold.rs`
>   cambian con la misma justificación, comparada también contra `dcae406`.
> - **Juez físico.** La comparación de «lo declarado de cada locus» con el
>   TOML deja de ir por identificador: compara las posiciones declaradas como
>   lista ordenada, como antes comparaba la lista en el orden del TOML, y la
>   comparación exacta de los loci ata cada declaración a su locus. Con la
>   versión por identificador, quitarle el identificador a un locus se
>   rechazaba ahí y dejaba sin alcanzar la guarda del juez del abanico que lo
>   comprueba.
> - Quitada la marca «provisional» de REQ-REF-39, DTC-103 y DTC-096 en lo que
>   la pausa de la entrega 1 y las decisiones del 30/10 resolvieron.
>   Actualizados los tres ejemplos de `uso-con-agentes` y su guía.
> - **Puertas** (con `python3` 3.12.7 del entorno del proyecto y `bcftools` y
>   `tabix` 1.24): formato y `clippy` limpios; batería 2208 pasan, 0 fallan,
>   14 ignoradas; pruebas de documentación 2; compilación optimizada; cinco
>   PASS de los validadores independientes; cobertura de los jueces dentro de
>   su techo, con `src/results.rs` en 725 comprobaciones y 153 sin disparar
>   (una nueva, que se dispara; nota en la línea base); Marimo limpio; 178
>   pruebas de Python con los 7 fallos esperados (los dos antiguos y los cinco
>   de la entrega 5).
>
> **Ronda 2 de la revisión adversarial (30/10/2026).** Un hallazgo grave y dos
> medios, aplicados:
>
> - **Grave: PCR con el catálogo final vacío.** `validate` aceptaba y `run`
>   moría con GPS-FMT-006, porque el tartamudeo clonal, conectado al leer la
>   sección, dejaba el evento con el catálogo de azar v6 y el manifiesto general
>   en la versión 4. Pasaba con regiones que no contienen nada y con una fila
>   de la matriz que retira el único detectado. Ahora un paso final sobre el
>   catálogo (`with_final_repeat_catalog`, en `validate`, en `run` y en la
>   reconstrucción del juez) aplica la regla de la PCR **después de retirar** y
>   desconecta el tartamudeo clonal si no queda ningún locus; el juez compara
>   el evento técnico con el catálogo publicado. La prueba ejecuta `run` en
>   los dos casos y exige `f4-technical-v5` y los FASTQ y el manifiesto de
>   evento idénticos a los de la misma ejecución sin la sección.
> - **Medio: efectiva sin su sección.** El juez físico leía la referencia sin
>   bases si la efectiva no declaraba la sección, aunque el TOML copiado sí, y
>   la reconstrucción moría con un defecto interno. Ahora mira también la
>   sección del TOML, y el rechazo es el de la sección.
> - **Medio: el compuesto que sustituye a dos tramos.** El fixture cambia a
>   un `(CA)×10` y un `(GAT)×7` separados por dos bases de flanco, con un
>   control que demuestra que sin el declarado salen las dos piezas.
> - Pruebas nuevas para el identificador repetido, las regiones en
>   coordenadas finales con inserción y la retirada por un cambio de base
>   escrito a mano; `validate_and_run_see_the_same_catalog` pasa a llamarse
>   `validate_and_run_accept_and_reject_the_same_scenarios`, porque no compara
>   catálogos. Quitada la marca «provisional» de REQ-CONF-142 y REQ-REF-40 en
>   lo decidido.
> - Cada arreglo tiene una prueba que muere con la variante contraria,
>   comprobado mutando una copia aparte. Puertas tras la ronda: formato y
>   `clippy` limpios; batería 2213 pasan, 0 fallan, 14 ignoradas; pruebas de
>   documentación 2; compilación optimizada; cinco PASS de los validadores
>   independientes; cobertura de los jueces sin cambios (`src/results.rs` en
>   725 y 153); Marimo limpio; 178 pruebas de Python con los mismos 7 fallos
>   esperados.

### Entrega 4 — publicación y juez

- ☑ **F6.4.a.14.** Publicación del catálogo.
  - Depende de: F6.4.a.10; DP-SSR-13 decidida.
  - Procede de: REQ-MAN-98; DTC-104.
  - Cambia: `src/results.rs`, `schemas/` (la efectiva v6 y el manifiesto v6
    **retocados in situ** por la nota alfa de DTC-036 —el texto anterior
    decía «efectiva v7, manifiesto v7», superado por esa nota y por la forma
    en que la entrega 3 ya los retocó— y, con DP-SSR-13 b, el esquema del
    catálogo, nuevo en su versión `"1"`), `src/execution.rs`.
  - Se comprueba con: los esquemas nuevos validados por el validador de
    esquemas en Python que ya lanzan las pruebas, y la puerta
    `check-independent-validators.sh` sobre un directorio nuevo; recuentos de
    descarte que suman los candidatos; sin detección, la huella intacta.
- ☑ **F6.4.a.15.** Juez independiente del catálogo.
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

#### Diseño concreto de la entrega 4 (01/11/2026), antes de tocar el código

**Estado:** escrito antes de tocar el código y detenido en la decisión D4-1,
que el responsable resolvió el 01/11/2026 (abajo). Lo demás de este apartado
está fijado por las normas o es mecánico; lo decidió el implementador y el
responsable lo aceptó tal cual el 01/11/2026.

**Lo que hay hoy (comprobado en `45af26b`).** Cada locus detectado viaja
**entero** en `repeat_loci` de la configuración efectiva y del manifiesto
general, como un declarado: identificador, intervalo, motivo, tope, mínimos,
auditoría de normalización, abanico de alelos con la anotación de cada uno
(unidades, pureza, bloques, interrupciones) y el bloque del deslizamiento con
las **tres limitaciones palabra por palabra repetidas en cada locus**
(`src/results.rs`, `repeat_loci_json`, `slippage_json`,
`repeat_loci_manifest_json`). El manifiesto añade por alelo las copias y la
frecuencia realizadas. La efectiva lleva además `repeat_detection` con la
versión del algoritmo, los suelos, el tope, el deslizamiento común, las
regiones y `retired_loci`. Los recuentos de descarte por etapa que calcula el
barrido (`DetectionCounts`) **no se guardan** ni se publican, y la puntuación
de cada tramo se pierde al convertirlo en locus. Esa forma es la que la
entrega 3 llamó «publicación mínima» a la espera de esta entrega, y es la que
el apartado 11.7 midió en unos 9,5 kB por locus entre los dos documentos (de
ellos, alrededor de 1 kB por documento son el bloque del deslizamiento con las
limitaciones).

**Decidido por las normas o mecánico (lo decide el implementador):**

1. **Ficheros publicados**, con la forma del par de la matriz de SNV
   (DTC-054): `repeat_catalog.tsv` (artefacto `repeat-catalog-tsv`, clase
   `repeat_catalog_tsv`) y `repeat_catalog.metadata.json` (artefacto
   `repeat-catalog-metadata`, clase `repeat_catalog_metadata`), solo con la
   sección de microsatélites. Sin ella no se escribe nada nuevo y todo byte
   publicado queda idéntico.
2. **Columnas del TSV**, con cabecera y en este orden:
   `locus_id`, `sequence_id`, `start`, `end`, `motif`, `score`, `origin`.
   Coordenadas en **base uno y cerradas**, las de los identificadores
   `SSR-…`, las regiones y `retired_loci`. Una fila por locus del **catálogo
   final** (después de retirar), en su orden: por secuencia en el orden de la
   referencia y por inicio. `motif` es el motivo; un declarado compuesto lleva
   sus motivos separados por coma, como ya hace `retired_loci`. `score` es la
   puntuación del detector para un detectado (la de su propio tramo, +2 por
   base igual a la de un periodo antes y −7 por distinta o no ACGT, que el juez
   puede rehacer con aritmética propia); un declarado no tiene puntuación del
   detector y lleva `.` (sus fronteras son las del usuario).
3. **Metadato JSON**, `document_kind: "repeat_catalog_metadata"`,
   `schema_version: "1"`, `coordinate_system: "base1_closed"`: la referencia
   final (identificador y huella del FASTA final publicado), la tabla
   (`path` y `sha256`), la versión del algoritmo (`ssr-detection-v1`), los
   criterios completos (suelos, pureza y tirada mínimas, tope, regiones), los
   **recuentos del barrido por etapa y longitud de motivo** (intervalos del
   detector; descartados en cada una de las ocho etapas; aceptados), y la
   **contabilidad del paso de aceptados a catálogo** por longitud de motivo:
   aceptados = sustituidos por un declarado + retirados por una variante +
   fuera por falta de copias haplotípicas (decisión 1 del 30/10) + publicados
   como detectados. Por longitud de motivo, intervalos = descartes de las ocho
   etapas + aceptados (lo que REQ-MAN-98 llama «sumar lo que el barrido
   encontró»), y el juez lo exige. Los **tramos retirados** con su causa van
   aquí, con la posición de la variante en un campo propio para que el juez
   compruebe que cae dentro sin leer la frase.
4. **Enlace con huella desde la efectiva**: `repeat_detection.catalog` con
   `metadata_path` y `metadata_sha256`; el metadato lleva la huella de la
   tabla. Los dos ficheros figuran además en la lista de artefactos del
   manifiesto, como la matriz.
5. **Esquemas**: `schemas/repeat-catalog-metadata-v1.schema.json` (nuevo) y la
   efectiva v6 retocada in situ para el enlace; el TSV se juzga por el juez,
   como el de la matriz. Validados por el validador de Python que lanzan las
   pruebas.
6. **Juez** `src/repeat_catalog_validation.rs`, con la lista de REQ-VAL-64 y
   DP-SSR-14, contra el FASTA final publicado y con aritmética propia: no
   importa nada de `src/repeat_detection.rs`, de la anotación de tractos ni de
   la conversión. Se encadena en `validate_published_resources_v3`
   (`src/results.rs`) después de rehacer la referencia final y antes de atar
   los loci al TOML. **No entra en `check-independent-validators.sh`**: ese
   script certifica herramientas externas (samtools, sha256sum, tar,
   jsonschema, Docker) sobre fixtures, no jueces propios.

**D4-1, resuelta el 01/11/2026 por decisión del responsable: opción (C).**
**Qué sale de `repeat_loci` y dónde van los abanicos de los detectados.** El apartado 11.7 aceptó DP-SSR-13 (b) «tal como está escrita»,
y su texto dice que la efectiva lleva los criterios, el deslizamiento común
«con sus tres limitaciones una sola vez» y el enlace, que los abanicos van en
un TSV propio y que los documentos JSON dejan de crecer con los loci. Pero la
forma concreta no está en ninguna norma, y lo que hoy dicen REQ-MAN-95 («la
configuración efectiva y el manifiesto general deben publicar, por alelo, el
origen…»), REQ-VAL-60 y REQ-VAL-62 («en la configuración efectiva, los loci
—con sus alelos— deben coincidir exactamente…») supone los abanicos dentro de
los dos documentos. Opciones:

| Opción | Consecuencia |
|---|---|
| (A) **Nada sale de `repeat_loci`**: el catálogo se publica además, como índice juzgado contra el FASTA, y `retired_loci` pasa de la efectiva al metadato | Ningún consumidor cambia y la entrega queda en su tamaño; los documentos JSON siguen creciendo unos 9,5 kB por locus (unos 540 MB con los 56 736 loci del tomate reducido), que es lo que 11.7 juzgó inviable |
| (B) **Los detectados salen de los dos documentos**: identidad en el catálogo, abanicos y lo realizado por alelo en un TSV de alelos, contabilidad del deslizamiento por locus en una tabla, limitaciones una vez en `repeat_detection`; los declarados se quedan en `repeat_loci` | Los JSON dejan de crecer con los detectados; hay que reescribir REQ-MAN-95, REQ-MAN-96, REQ-VAL-60 y REQ-VAL-62 y las seis funciones de `src/results.rs` que leen `repeat_loci` de lo publicado (contadas con `grep`), y decidir las columnas del TSV de alelos; la entrega pasa a ser dos o tres veces mayor |
| (C) (A) ahora y (B) en una entrega propia, después del juez | El juez del catálogo y la invariante llegan ya; la reestructuración se decide con su norma escrita; durante un tiempo el catálogo duplica lo que dice `repeat_loci` |
| (D) (A) más quitar de cada detectado lo repetido: el bloque fijo del deslizamiento y las limitaciones, una vez en `repeat_detection` | Ahorra alrededor del 20 % por locus sin mover los abanicos; cambia la forma de cada detectado en los dos documentos y el juez de los loci |

**Recomendación del implementador: (C).** Es lo más pequeño que cumple las
normas tal como están hoy, entrega ya el juez independiente, y deja la parte
que cambia lo que ve un consumidor (y reescribe cuatro requisitos) para cuando
el responsable fije su forma.

**Decisión del responsable (01/11/2026): (C).** En esta entrega se hace (A):
la lista de loci de la configuración efectiva y del manifiesto general **no
cambia**, y el catálogo aparte (tabla, metadato y enlace con huella en
`repeat_detection.catalog`) se publica además, con el juez independiente.
`retired_loci` pasa de la efectiva al metadato. La parte (B) es la entrega 4
bis, inmediatamente después de esta.

> **Estado, 01/11/2026: entrega 4 implementada, pendiente de la revisión
> adversarial.** Sin *commit*.
>
> - **Qué se publica.** Con la sección, `repeat_catalog.tsv` y
>   `repeat_catalog.metadata.json` en la raíz del paquete, como en el diseño de
>   arriba; la efectiva gana `repeat_detection.catalog` y pierde
>   `repeat_detection.retired_loci`, que pasa al metadato con la posición y el
>   origen de la variante. Para publicarlo, el barrido guarda sus recuentos,
>   los tramos que sustituye un declarado y los que quedan fuera por falta de
>   copias haplotípicas (`DetectionScan` en `src/application.rs`), cada
>   detectado guarda su puntuación y cada retirado su variante. Esquemas: el
>   nuevo `repeat-catalog-metadata-v1.schema.json` y la efectiva v6 retocada
>   in situ.
> - **El juez** (`src/repeat_catalog_validation.rs`, 57 comprobaciones; 59 tras la revisión) corre
>   dentro del juez de los recursos publicados, después de atar los loci de la
>   efectiva al TOML y antes de una comprobación nueva de `src/results.rs` que
>   rehace la tabla y el metadato con el productor
>   (`validate_reconstructed_repeat_catalog`): es la que ata los descartes por
>   etapa y la exhaustividad, que el juez independiente no puede ver. El juez
>   compara además la tabla con `repeat_loci` de la efectiva, lo único que caza
>   un escritor de la tabla que perdiera una fila (la reconstrucción usa el
>   mismo escritor). No entra en `check-independent-validators.sh` (ver el
>   diseño).
> - **Pruebas** (`tests/repeat_catalog_validation.rs`): una referencia con
>   tractos plantados en posiciones escritas en la prueba, incluidos los que
>   cada criterio descarta y el tracto `GATGATGATGAC × 6` de pureza real 0,792
>   (57/72); el catálogo publicado se compara con una tabla escrita a mano; los
>   recuentos se suman con aritmética propia; el metadato y la efectiva pasan
>   su esquema por el validador de Python (y una instancia alterada de cada uno
>   no). Las 57 comprobaciones del juez las dispara un caso que altera a mano
>   la tabla, el metadato, el enlace, el manifiesto o la efectiva publicados y
>   exige el motivo exacto; ningún esperado sale del detector.
> - **Prueba de independencia**, en una copia extraída con `git archive` de
>   `45af26b` más los ficheros de esta entrega, con su propio
>   `CARGO_TARGET_DIR`. `observe_tract` (`src/tandem_repeat.rs`, SHA-256
>   `dbb2bb97…`) mutado para sumar un 5 % de las bases a las cubiertas
>   (`433ce4f5…`): el productor publica `SSR-2-2000-2071-TGA`, de pureza real
>   57/72, y **la ejecución muere en el juez** con «su pureza 57/72 no alcanza
>   0.8» (las once pruebas que publican la plantilla fallan con ese motivo).
>   Con, además, la pureza del juez calculada con `observe_tract`
>   (`src/repeat_catalog_validation.rs`, `30530647…` → `bde181e7…`), la
>   ejecución se publica con el locus dentro y mueren la comparación con la
>   tabla escrita a mano y la guarda léxica: es el «debe fallar contra» de
>   F6.4.a.15. Restaurados los dos ficheros, con su SHA de antes.
> - **Mutaciones** (misma copia; cada una con SHA antes, mutado y restaurado
>   igual al de antes): sin la comprobación de pureza, con el suelo rebajado en
>   una copia, con el solape medido como estricto, sin la suma de los
>   intervalos del detector, sin la de la puntuación, sin comparar el origen con
>   la efectiva y aceptando cualquier variante de la matriz: cada una mata la
>   prueba de su comprobación. Publicar el catálogo sin la sección mata las
>   cuatro huellas de `tests/repeat_detection_frozen_identity.rs` y la prueba
>   de «sin la sección»; no contar los sustituidos hace que la publicación muera
>   en el juez; quitar la comprobación de la reconstrucción mata
>   `the_reconstruction_ties_what_the_judge_cannot_see`.
> - **Huellas.** Sin la sección no cambia un byte
>   (`tests/repeat_detection_frozen_identity.rs` intacta). Con ella cambian a
>   propósito la efectiva y el manifiesto, y aparecen los dos ficheros del
>   catálogo: comparado campo a campo contra el binario de `45af26b` compilado
>   aparte (dos escenarios de la revisión de la fusión del VCF y los dos de la
>   evidencia de DT-009 con loci), no cambia nada más. La tabla normalizada de
>   `tests/repeat_scale_frozen_identity.rs` conserva los literales de `dcae406`
>   (la normalización deja fuera el catálogo); su tabla de documentos y las dos
>   huellas de `tests/technical_catalog_contract.rs` se volvieron a medir, con
>   la justificación escrita en cada prueba.
> - **Puerta de cobertura**: el juez entra con 57 comprobaciones (59 tras la revisión) y techo 0;
>   `src/results.rs` pasa de 725 a 726 (la reconstrucción del catálogo) y sigue
>   en 153 sin disparar.
> - **Puertas** (con `python3` 3.12.7 del entorno del proyecto y `bcftools` y
>   `tabix` 1.24): formato y `clippy` limpios; batería 2237 pasan, 0 fallan,
>   14 ignoradas; pruebas de documentación 2; compilación optimizada; cinco
>   PASS de los validadores independientes; cobertura de los jueces dentro de
>   su techo en los quince; Marimo limpio; 178 pruebas de Python con los 7
>   fallos esperados (los dos antiguos y los cinco de la entrega 5).

> **Revisión adversarial de la entrega 4 (01/11/2026).** Ningún defecto del
> programa; dos hallazgos medios de pruebas que no podían fallar, aplicados:
>
> - **Medio: la variante de un retirado no estaba fijada a su posición.** Los
>   casos solo usaban identificadores inexistentes. Ahora: el identificador
>   real de otra fila de la matriz puesto en la posición del retirado; un
>   cambio directo con ese identificador en otra posición, en otra secuencia
>   y de otra clase; y el borde (primera y última base dentro, la anterior
>   fuera). Mueren la búsqueda en la matriz que ignora la posición, la de los
>   cambios directos que ignora secuencia, posición o clase, y `<`→`<=` o
>   `>`→`>=` en el borde.
> - **Medio: «una base A/C/G/T a cada lado» solo se probaba en el borde
>   izquierdo de la secuencia.** La referencia de las pruebas gana un `A×12`
>   con una `N` delante, otro con una `N` detrás y otro que acaba en la última
>   base; mueren quitar las dos condiciones de base y quitar el límite derecho
>   (que hace entrar en pánico al juez).
> - **Menores aplicados**: dos aceptados que comparten una sola base (la
>   referencia gana `(CA)×8` pegado a once `A`, retiradas por una fila de la
>   matriz); cada campo de la cabecera del metadato y la huella de la tabla
>   dentro de él; **el tamaño de los dos ficheros frente al manifiesto**, como
>   en la matriz de SNV (dos comprobaciones nuevas); el comentario de la
>   función que ata los loci al TOML, que había quedado encima de la nueva;
>   `path` y `schema` del rechazo nombran el fichero que falla (el texto fijo
>   del código sigue hablando de un documento JSON: cambiarlo exigiría un
>   código nuevo); los dos recorridos sobre los declarados ya no recorren las
>   filas; y cada caso publica su propio paquete, que se borra al acabar, en
>   vez de dejar una plantilla en el directorio temporal.
> - **Anotado sin tocar**: el juez completo no compara el tamaño del FASTA
>   publicado con el manifiesto (laguna previa a esta entrega).
> - Cada prueba nueva se vio morir contra su mutación en la copia aparte, con
>   SHA antes, mutado y restaurado iguales. Ninguna huella congelada cambia
>   en esta ronda.
> - Puertas tras la ronda: formato y `clippy` limpios; batería 2240 pasan, 0
>   fallan, 14 ignoradas; pruebas de documentación 2; compilación optimizada;
>   cinco PASS de los validadores independientes; cobertura de los jueces
>   dentro de su techo (el del catálogo, 59 y 0; `src/results.rs`, 726 y 153);
>   Marimo limpio; 178 pruebas de Python con los 7 fallos esperados.

> **Revisión de la entrega 4:** lo publicado se valida sin el detector.

### Entrega 4 bis — los detectados salen de los documentos JSON

Añadida el 01/11/2026 por la decisión D4-1 (opción C) del responsable. **No
implementada.** Va inmediatamente después de la entrega 4.

- ☑ **F6.4.a.15b.** Cambio normativo previo.
  - Depende de: entrega 4.
  - Procede de: DP-SSR-13 (b), aceptada en el apartado 11.7; D4-1.
  - Cambia: `especificaciones-0.1.md` (REQ-MAN-95, REQ-MAN-96, REQ-VAL-60 y
    REQ-VAL-62, que hoy suponen los abanicos dentro de la efectiva y del
    manifiesto), DTC-104 y, si hace falta, la arquitectura.
  - Fija: qué queda de cada detectado en los dos documentos (nada, según
    DP-SSR-13 b); las columnas de la tabla de alelos (abanico, frecuencias
    declarada y efectiva, origen por alelo, copias y frecuencia realizadas, y
    si la anotación de cada alelo se publica o se deriva de la secuencia);
    dónde va la contabilidad del deslizamiento por locus (linajes y pasos); y
    las tres limitaciones una sola vez en `repeat_detection`.
  - Se comprueba con: revisión del responsable antes de tocar el código.
  - **Estado, 01/11/2026: hecho.** El responsable aprobó la propuesta de abajo
    y decidió A1 a A3 según la recomendación. Aplicada ese día a la
    especificación (REQ-MAN-29, REQ-MAN-95, REQ-MAN-96, REQ-MAN-98, REQ-VAL-60,
    REQ-VAL-62 y REQ-VAL-64), a las decisiones técnicas (DTC-035, DTC-045 en su
    ampliación del paquete fallido, DTC-095, DTC-097 y DTC-104) y al apartado
    16.4 de la arquitectura, con nota fechada donde un texto quedó superado.
    Sin código.
- ☑ **F6.4.a.15c.** Publicación y jueces.
  - Depende de: F6.4.a.15b.
  - Procede de: la propuesta normativa aprobada de abajo (F6.4.a.15b).
  - Cambia: `src/results.rs`, en la escritura (`repeat_loci_json`,
    `repeat_loci_manifest_json`, `repeat_output_budget_json`,
    `repeat_detection_json` y el escritor nuevo de `repeat_alleles.tsv`) y en
    los jueces que leen `repeat_loci` de lo publicado y cambian: el del
    abanico (`validate_published_slippage_fans`), el de los loci frente al
    TOML (`validate_published_repeat_loci`), el de la contabilidad de tracto
    (`judge_tract_output_accounting`), el de las cotas frente al TOML
    (`validate_published_tract_bounds`), el de la efectiva frente al
    manifiesto (`validate_effective_config_v3_against_manifest_bytes`), el de
    los objetivos frente al TOML
    (`validate_effective_targets_against_configuration`), la forma de la
    sección del manifiesto (`validate_repeat_loci_section`) y el orquestador
    (`validate_published_resources_v3`); `src/execution.rs` (orden de
    escritura y paquete fallido); `src/repeat_catalog_validation.rs` (la tabla
    de alelos entra en su alcance); `schemas/` (efectiva y manifiesto
    retocados in situ por la nota alfa de DTC-036). **La tabla de alelos no
    lleva esquema**: la juzga el código, como la del catálogo y la de la
    matriz. No cambian, aunque la lista anterior los nombraba, el juez de la
    referencia simulada (`published_simulated_tracts`) ni el de la referencia
    final (`judge_inserted_reference`), que ya saltan los detectados, ni el de
    los nombres (`validate_published_reference_names`), que solo lee la
    referencia de la sección.
  - Se comprueba con: sin la sección de detección, todo byte publicado
    idéntico; con ella, los documentos JSON no crecen con los detectados
    (medido con dos tamaños de catálogo); cada comprobación nueva disparada
    con lo publicado alterado a mano. **Debe fallar contra:** un abanico de un
    detectado que se quede en la efectiva; un escenario con solo detectados y
    tartamudeo clonal que el juez de los objetivos tome por catálogo vacío
    (hoy lo decide mirando si `repeat_loci` de la efectiva está vacío, y con
    solo detectados lo estaría siempre); un juez que atribuya las bases del
    VCF solo a los declarados; un escritor que pierda la misma fila en el
    catálogo y en la tabla de alelos (lo caza una tabla escrita a mano en la
    prueba).

> **Estado, 01/11/2026: F6.4.a.15c implementado, pendiente de la revisión
> adversarial.** Sin *commit*.
>
> - **Qué se publica.** Con la sección, `repeat_loci` y
>   `repeat_output_budget.per_locus` de la efectiva y `repeat_loci.loci` del
>   manifiesto llevan solo los declarados (la lista del manifiesto puede quedar
>   vacía); `repeat_detection` gana, en `slippage`, las versiones, el ámbito,
>   la versión de su catálogo, las tres limitaciones y los linajes,
>   `frequency_normalization` con sus tres constantes y el enlace
>   `alleles_table`; `repeat_alleles.tsv` (artefacto `repeat-alleles-tsv`,
>   clase `repeat_alleles_tsv`, en el mínimo) lleva las trece columnas de la
>   propuesta, solo con la cabecera si no queda ningún detectado. Se escribe
>   justo después de la verdad y el VCF y antes del evento técnico y de la
>   efectiva, y entra en la unidad genotípica con el par del catálogo (que
>   hasta hoy no estaba en ella). Las frecuencias se escriben con el decimal
>   más corto de ida y vuelta de Rust, sin exponente. Comprobado con un guion
>   en Python sobre 250 loci: cada fila de la tabla es, campo a campo, lo que la
>   efectiva y el manifiesto de `afda48b` publicaban de ese alelo.
> - **Jueces.** El del catálogo (`src/repeat_catalog_validation.rs`) gana el
>   enlace, la huella y el tamaño de la tabla de alelos, su lector estricto y
>   la lista de la propuesta; las filas declaradas se comparan con
>   `repeat_loci` y las detectadas con la tabla de alelos. Nuevo en
>   `src/results.rs`: el juez del abanico de los detectados
>   (`validate_published_detected_fans`), que corre justo después del de los
>   declarados y antes de atar nada al TOML, con el mismo recorrido escrito
>   aparte del productor; y la atadura de la tabla al TOML
>   (`validate_published_repeat_alleles`, copias contra la verdad y la tabla
>   rehecha byte a byte con las copias publicadas y las bases del VCF
>   recontadas), detrás de recontar la verdad y el VCF. La contabilidad de
>   tracto atribuye cada registro de repetición al primer locus del catálogo
>   y exige la parte VCF de cada fila de cada detectado. El lector del
>   manifiesto admite la lista de loci vacía si el manifiesto lista el
>   catálogo o es de un paquete fallido; el juez de los dos documentos admite
>   una efectiva sin loci frente a un manifiesto con ellos si la efectiva
>   tiene la sección.
> - **Defecto previo encontrado y arreglado: el evento técnico frente al
>   TOML.** El juez de los objetivos comparaba `technical_event` con la
>   configuración recién leída del TOML, que aún no ha detectado nada: con
>   solo loci detectados y tartamudeo, el plan del tartamudeo salía apagado en
>   lo esperado y **la ejecución moría en su propio juez**, también en
>   `afda48b` (reproducido en una copia de `afda48b` con PCR, tartamudeo de
>   base y un solo detectado: «technical_event efectivo no coincide
>   exactamente con el TOML resuelto»). Es lo que el riesgo 1 de la propuesta
>   veía a medias: el criterio del catálogo vacío no era lo único mal. La
>   comparación pasa al juez completo, contra la configuración del **catálogo
>   final** (`validate_published_technical_event_block`, justo después de
>   `with_final_repeat_catalog`), y el criterio del catálogo vacío desaparece:
>   lo decide el propio catálogo final.
> - **Detalle de la propuesta que no se sostenía (apartado 7):** «secuencias
>   de A, C, G y T». Un tramo detectado puede llevar una `N` (cuenta como base
>   que no encaja en la pureza) y sus alelos la heredan; lo hizo ver
>   `a_lone_n_inside_a_long_tract_counts_as_a_purity_mismatch`. El juez admite
>   A, C, G, T y N, como el lector de los alelos del manifiesto. Ninguna norma
>   decía lo contrario: REQ-VAL-64 solo pide «su forma».
> - **Pesos medidos** con los binarios de `afda48b` y de esta entrega en
>   `--release`, sobre una referencia sintética con 250 y 1 000 microsatélites
>   plantados (el detector halla 250 y 1 005), diploides HWE, VCF de texto:
>
>   | Fichero | 250 loci, 10 ind. | 1 005 loci, 10 ind. | 250, 100 ind. | 1 005, 100 ind. |
>   |---|---:|---:|---:|---:|
>   | Efectiva, `afda48b` | 1 203 723 B | 4 823 707 B | 1 570 032 B | 6 285 341 B |
>   | Efectiva, ahora | 7 412 B | 7 415 B | 8 779 B | 8 782 B |
>   | Manifiesto, `afda48b` | 1 264 160 B | 5 060 387 B | 1 705 584 B | 6 824 285 B |
>   | Manifiesto, ahora | 8 488 B | 8 495 B | 9 858 B | 9 865 B |
>   | `repeat_alleles.tsv` | 127 770 B | 512 911 B | 207 529 B | 830 621 B |
>
>   Los dos documentos cambian 3 y 7 bytes entre 250 y 1 005 loci (los dígitos
>   de los totales); la tabla de alelos crece 510 B por locus con 10
>   individuos y 825 B con 100, como predijo el apartado 9. La prueba
>   `the_json_documents_do_not_grow_with_the_detected_loci` lo exige con 12 y
>   48 plantados.
> - **Pruebas** (`tests/repeat_catalog_validation.rs`, `tests/repeat_detection_run.rs`
>   y una de `src/execution.rs`): cada comprobación nueva de los dos jueces la
>   dispara un caso que altera a mano la tabla de alelos, el catálogo, la
>   efectiva o el manifiesto y exige el motivo exacto; los detectados se fijan
>   con una lista escrita a mano; con PCR, tartamudeo clonal y solo detectados
>   la ejecución se publica; el paquete fallido conserva o pierde los tres
>   ficheros con la unidad.
> - **Mutaciones**, en una copia extraída con `git archive` de `afda48b` más
>   los ficheros de esta entrega, con su propio directorio de compilación, SHA
>   antes, mutado y restaurado (igual al de antes) en cada una: veinticuatro
>   mueren con la prueba prevista — los detectados en la efectiva, en el
>   manifiesto o en `per_locus`; la misma fila perdida en el catálogo y en la
>   tabla de alelos con el metadato rehecho (el juez y la reconstrucción la
>   aceptan; la matan las dos tablas escritas a mano); el catálogo vacío
>   decidido por `repeat_loci` de la efectiva; la atribución del VCF solo a
>   los declarados; sin la parte VCF de cada fila; el abanico sin frecuencia,
>   con la contabilidad solo de la primera fila o sin cotejar los linajes; la
>   tabla sin copias frente a la verdad o sin comparación byte a byte; el juez
>   del catálogo sin suma de copias, sin frecuencia realizada, sin el orden de
>   los loci, contando los declarados contra todas las filas o sin la huella
>   del enlace; la tabla de alelos o el catálogo fuera de la unidad
>   genotípica; el manifiesto fallido exigiendo un locus; la lista vacía
>   siempre admitida; la efectiva sin loci siempre admitida; y **el generador
>   con el espectro invertido** (cada alelo con la frecuencia del simétrico),
>   que hace morir la publicación con solo detectados en el juez del abanico
>   de los detectados («alelo 0: publicada 0.75, rehecha 0.25 = 1/4»); con ese
>   juez quitado, la misma publicación pasa: es la única guarda de eso.
> - **Una mutación sobrevive, y es un hueco previo:** invertir solo la
>   frecuencia **efectiva** de un abanico generado (la declarada intacta). Ningún
>   juez independiente compara la efectiva con la declarada, ni para los
>   detectados ni para los declarados con generador: la ata solo la
>   reconstrucción con el productor. Anotado como DT-064 en el registro de
>   deudas (no lo toca esta entrega).
> - **Huellas.** Sin la sección no cambia un byte (las cuatro de
>   `tests/repeat_detection_frozen_identity.rs` intactas; los cinco escenarios
>   de la revisión de la fusión del VCF sin la sección y las ejecuciones
>   amplificada sin loci y sin evento de `tests/technical_catalog_contract.rs`,
>   idénticos al binario de `afda48b`). Con ella cambian a propósito la
>   efectiva y el manifiesto y aparece la tabla de alelos: comparado campo a
>   campo contra el binario de `afda48b` compilado aparte (los dos escenarios
>   de muchos loci, los dieciséis de la revisión con la sección y los dos de
>   la evidencia de DT-009 con loci), no cambia nada más. La tabla normalizada
>   de `tests/repeat_scale_frozen_identity.rs` conserva los literales de
>   `dcae406`; su tabla de documentos y las dos huellas de
>   `tests/technical_catalog_contract.rs` se volvieron a medir, con la
>   justificación escrita en cada prueba.
> - **Puerta de cobertura**: el juez del catálogo pasa de 59 a 78
>   comprobaciones con techo 0; `src/results.rs` de 726 a 739 y de 153 a 154
>   sin disparar (la derivación de la subsemilla del abanico de un detectado,
>   inalcanzable en la práctica como su gemela de los declarados).
> - **Puertas** (`python3` 3.12.7 del entorno del proyecto, `bcftools` y
>   `tabix` 1.24): formato y `clippy` limpios; batería 2254 pasan, 0 fallan,
>   14 ignoradas; pruebas de documentación 2; compilación optimizada; cinco
>   PASS de los validadores independientes; cobertura dentro de su techo en
>   los quince jueces; Marimo limpio; 178 pruebas de Python con los 7 fallos
>   esperados.

#### Propuesta normativa de la entrega 4 bis (F6.4.a.15b, 01/11/2026)

**Estado: aprobada por el responsable el 01/11/2026, con A1 a A3 decididas
según la recomendación, y aplicada ese día** a `especificaciones-0.1.md`,
`decisiones-tecnicas-0.1.md` y `arquitectura-0.1.md` (ver el estado de
F6.4.a.15b). El código no se ha tocado. Lo que sigue se conserva como se
escribió; el apartado 8 lleva el texto que se aplicó.

**Aviso previo: el cambio toca más normas que las cuatro que nombra el paso.**
Además de REQ-MAN-95, REQ-MAN-96, REQ-VAL-60, REQ-VAL-62 y DTC-104, sacar los
detectados de los dos documentos choca con dos normas generales que no son de
microsatélites: **DTC-035** («cada manifiesto contiene directamente todo lo
que REQ-MAN obliga a registrar […]; V1 no fragmenta manifiestos en JSONL ni
**tablas auxiliares**») y **REQ-MAN-29** («el manifiesto debe registrar los
loci y parámetros configurados […] y conservar los valores realizados»). La
dirección ya la fijó el responsable al aceptar DP-SSR-13 (b) el 29/10/2026
(«los abanicos de alelos, en un TSV propio»; «los documentos JSON dejan de
crecer con los loci»), pero el texto de esas dos normas no se cambió. Aquí se
propone una nota fechada para cada una y se deja como decisión abierta A1. Por
consecuencia mecánica cambian también REQ-MAN-98 y REQ-VAL-64 (el catálogo y
su juez, que hoy comparan la tabla con `repeat_loci`), DTC-097 (la forma
publicada de cada alelo y del bloque `slippage`), una frase de DTC-095 y la
ampliación de DTC-045 sobre el paquete fallido, que nombra `repeat_loci`.

##### 1. Lo que hay hoy, comprobado

Comprobado en `afda48b` leyendo `src/results.rs` y publicando dos escenarios
propios (apartado 9 de esta propuesta). Cada **detectado** publica en la
configuración efectiva, dentro de `repeat_loci`: `id`, `origin`,
`sequence_id`, `start`, `end`, `motifs`, `motif_lengths`, `max_units`
(siempre el tope común de la sección), `admission` (siempre los mínimos
comunes de `repeat_admission`), `frequency_normalization` (siete campos),
`alleles` (trece campos por alelo: secuencia, origen, las dos frecuencias y
nueve de anotación) y `slippage` (catorce campos, de ellos 743 bytes de las
tres limitaciones). Además, una entrada en `repeat_output_budget.per_locus`
(identificador, cota de la verdad y bases del VCF). El manifiesto general
repite el locus sin la normalización y con `drawn_copies` y, por alelo,
`realized_copies` y `realized_frequency`. Medido: unos 4,2 kB por locus en
cada documento, más 0,1 kB de la cota por locus.

##### 2. Qué queda de cada detectado en la efectiva y en el manifiesto: nada

- `repeat_loci` de la efectiva lleva **solo los declarados**, con su forma de
  hoy; sin declarados es `[]`.
- `repeat_output_budget.per_locus` lleva solo los declarados. Los totales
  (`truth_tract_bases_bound`, `vcf_tract_bases_bound`,
  `total_tract_bases_bound`) siguen cubriendo **todos** los loci.
- `repeat_loci.loci` del manifiesto lleva solo los declarados y **puede quedar
  vacía**. La sección sigue existiendo mientras el catálogo final tenga algún
  locus, detectado o declarado, porque lleva la referencia, los mínimos y la
  contabilidad de tracto total.
- **Los declarados se quedan como están**, con su bloque entero (también las
  limitaciones si generan sus alelos). No hay norma que pida moverlos y son
  pocos; pueden ser compuestos, llevar alelos escritos a mano con
  interrupciones, una posición declarada o insertada y un generador propio, y
  meterlos en la tabla obligaría a columnas que un detectado nunca usa.
  Consecuencia para quien lea los alelos: los de los declarados en los JSON,
  los de los detectados en la tabla; `repeat_catalog.tsv` sigue siendo el
  índice de todos, con su origen.

##### 3. Lo que se publica una sola vez en `repeat_detection` de la efectiva

- `slippage`: a lo de hoy (`model`, `step_units`, `mean_events`,
  `expansion_bias`) se añaden `algorithm`, `count_algorithm`,
  `bias_algorithm`, `rng_scope`, `rng_scope_catalog_version`, `limitations`
  (las tres, palabra por palabra) y `lineages`, que es común a todos los
  detectados porque es el número de copias haplotípicas.
- `frequency_normalization`: las tres constantes (`tolerance`,
  `sum_algorithm`, `canonical_order`). Las cuatro sumas por locus **no** se
  publican para un detectado: salen de la columna de frecuencias con el
  algoritmo publicado, y en un abanico generado son sumas de `k/L`.
- `alleles_table`: el enlace con huella a la tabla nueva, `{"path":
  "repeat_alleles.tsv", "sha256": …}`.
- `catalog`, `max_units`, `min_copies` y `regions`, como hoy; los mínimos de
  pureza y tirada siguen en `repeat_admission`.

##### 4. La tabla nueva: `repeat_alleles.tsv`

Artefacto `repeat-alleles-tsv`, clase `repeat_alleles_tsv`, en la raíz del
paquete junto al catálogo. Se escribe **siempre que exista la sección** (solo
la cabecera si no queda ningún detectado), igual que el catálogo; sin la
sección no se escribe y todo byte publicado queda idéntico. Una fila por alelo
de cada **detectado**, en el orden del catálogo y, dentro de cada locus, por
`allele_index`. Columnas, con cabecera y en este orden:

| Columna | Qué es | Dónde está hoy |
|---|---|---|
| `locus_id` | Identificador `SSR-…` del detectado; une con `repeat_catalog.tsv` | `id` |
| `allele_index` | Posición del alelo en el orden canónico del abanico, desde 0 | orden de `alleles` |
| `sequence` | Bases del alelo | `alleles[].sequence` |
| `origin` | `generated` en toda fila: un detectado no admite alelos escritos a mano. Se conserva por alelo porque lo pide REQ-MAN-95 | `alleles[].origin` |
| `declared_frequency` | Proporción de linajes que terminan en el alelo | efectiva |
| `effective_frequency` | La anterior tras la normalización de DTC-049 | efectiva |
| `realized_copies` | Copias haplotípicas que el sorteo dio al alelo | manifiesto |
| `realized_frequency` | `realized_copies` entre las copias sorteadas del locus; `.` si no hay ninguna | manifiesto (`null`) |
| `drawn_events` | Pasos sorteados del locus | `slippage` |
| `applied_events` | Pasos aplicados del locus | `slippage` |
| `unapplied_events` | Pasos no aplicados al detenerse en un borde | `slippage` |
| `truth_tract_bases_bound` | Cota de bases de tracto del locus en la verdad | `repeat_output_budget.per_locus` |
| `vcf_tract_bases_bound` | Bases de tracto que el VCF publicado lleva en el locus, medidas al publicarlo (DTC-097) | `repeat_output_budget.per_locus` |

Las cinco últimas son **del locus** y se repiten en cada fila suya (decisión
abierta A3). Enteros en decimal; frecuencias como decimal de ida y vuelta, que
el juez lee como `f64` y compara bit a bit; `.` solo en `realized_frequency`;
sin comillas ni campos vacíos.

**Lo que no se publica porque se deriva**, y quién lo deriva:

- `drawn_copies`: la suma de `realized_copies` del locus.
- La **anotación de cada alelo** (unidades, pureza, tirada más larga, bloques e
  interrupciones): sale de `sequence` y del motivo del catálogo. Las unidades
  son bases entre longitud del motivo; lo demás, de la descomposición que la
  arquitectura (apartado 16.6) declara salida derivada y nunca entrada.
  Decisión abierta A2.
- Del locus: `max_units`, `admission` y la normalización (comunes, en la
  efectiva), `motifs`, `motif_lengths`, el intervalo y el origen (en el
  catálogo).

##### 5. Enlace, manifiesto y orden de escritura

- **Enlace:** `repeat_detection.alleles_table` de la efectiva, con ruta y
  huella. **Sin metadato JSON propio**: el contexto de la tabla (parámetros,
  versiones, ámbito, linajes, limitaciones) ya está en `repeat_detection`. La
  matriz de SNV lleva metadato (DTC-054) porque es un recurso reutilizable
  fuera de una ejecución; esta tabla no lo es.
- **Manifiesto:** la tabla figura en su lista de artefactos con ruta, clase,
  tamaño y huella, como el par del catálogo.
- **Orden:** se escribe **después de la unidad genotípica** (la verdad y, si se
  pidió, el VCF), porque las bases de tracto del VCF por locus solo se conocen
  al publicarlo, y **antes de la efectiva**, que lleva su huella. Las copias
  realizadas ya se conocen desde el preflight.

##### 6. Paquete fallido

El par del catálogo y la tabla de alelos **entran en la unidad genotípica**:
se conservan si y solo si se conserva la unidad, porque la verdad nombra los
loci detectados por su identificador y sin ellos no se explica. Sin unidad no
hay tabla, y el `repeat_loci` del manifiesto fallido registra solo los
declarados, con cero copias como hoy; los detectados se rehacen desde el TOML
copiado y la referencia. Hoy, sin unidad, el manifiesto fallido los registra
todos con cero copias: es un cambio visible, y por eso se propone texto para la
ampliación de DTC-045 (apartado 8).

##### 7. Jueces de la tabla

**Independientes, sin llamar al productor:**

- **El juez del abanico** (REQ-VAL-60; `validate_published_slippage_fans`, en
  `src/results.rs`): por cada detectado lee su fila del catálogo (secuencia,
  inicio, fin, motivo), el deslizamiento común, los linajes y el tope de
  `repeat_detection`, la referencia publicada y la semilla publicada; rehace
  el abanico con su propio recorrido, como hoy, y compara **fila a fila**: el
  número de filas del locus, `sequence` y `declared_frequency` bit a bit, y
  los tres recuentos de pasos **de cada fila** contra los rehechos, lo que
  además obliga a que sean iguales en todas las filas del locus.
- **El juez del catálogo** (`src/repeat_catalog_validation.rs`), ampliado con
  la forma y la aritmética de la tabla: cabecera exacta; trece campos por fila;
  los `locus_id` son exactamente las filas detectadas del catálogo, en su orden,
  sin que falte ni sobre ninguna; `allele_index` de 0 a k−1 sin huecos;
  `origin` siempre `generated`; secuencias de A, C, G y T, distintas dentro del
  locus; la suma de `realized_copies` del locus es el número de copias
  haplotípicas publicado; `realized_frequency` es copias entre esa suma, bit a
  bit, y `.` solo si la suma es cero; aplicados más no aplicados suman los
  sorteados; la huella y el tamaño coinciden con el enlace de la efectiva y con
  el manifiesto. La comprobación de hoy «la tabla dice lo mismo que
  `repeat_loci` de la efectiva» pasa a ser: las filas declaradas del catálogo
  son los loci de `repeat_loci`, y las detectadas, los loci de la tabla de
  alelos, en el mismo orden.
- **Lo realizado contra los ficheros** (REQ-VAL-62, en `src/results.rs`):
  `realized_copies` contra lo que se recuenta en los bytes de la verdad, y
  `vcf_tract_bases_bound` de cada locus contra lo que se recuenta en los bytes
  del VCF, atribuyendo cada registro de repetición al primer locus que lo
  contiene **entre todos los del catálogo** (leídos de `repeat_catalog.tsv`,
  no de `repeat_loci` del manifiesto, que ya no los tiene).

**Atadura al TOML** (la reconstrucción con el productor que REQ-VAL-62 ya
admite, porque la independencia la pone REQ-VAL-60): la tabla se rehace desde
el TOML copiado, con las columnas medidas sustituidas por lo recontado, y se
compara byte a byte.

**Lo que no ve**, igual que hoy: las primitivas de azar las comparte con el
productor (con sus propios vectores congelados), y la exhaustividad del
catálogo no se puede juzgar sin un detector.

**Para F6.4.a.15c, tres sitios donde el cambio puede romper en silencio**
(cada uno con su «debe fallar contra»):

1. `validate_effective_targets_against_configuration` decide que el catálogo
   está vacío —y que el tartamudeo clonal se desconecta— mirando si
   `repeat_loci` de la efectiva está vacío. Con solo detectados lo estaría
   siempre. Debe fallar contra: un escenario con solo detectados y tartamudeo
   clonal que se juzgue como catálogo vacío.
2. La atribución del VCF por locus en `judge_tract_output_accounting` usa los
   intervalos de `repeat_loci.loci` del manifiesto. Debe fallar contra: un
   juez que atribuya solo a los declarados.
3. La tabla y el catálogo salen del mismo recorrido del productor, y la
   reconstrucción usa los mismos escritores. Debe fallar contra: un escritor
   que pierda la misma fila en las dos tablas (lo caza una tabla escrita a mano
   en la prueba, como en la entrega 4).

##### 8. Texto normativo propuesto (aplicado el 01/11/2026)

**REQ-MAN-95.** *Hoy:* «La configuración efectiva y el manifiesto general
deben publicar, **por alelo y no por locus**, el **origen** de cada alelo de
repetición: declarado por el usuario o generado por deslizamiento. El origen
del alelo y el del locus deben concordar.»
*Propuesto:* «Se debe publicar, **por alelo y no por locus**, el **origen** de
cada alelo de repetición: declarado por el usuario o generado por
deslizamiento. Los alelos de un locus **declarado** se publican en la
configuración efectiva y en el manifiesto general; los de un locus
**detectado**, en la tabla de alelos de REQ-MAN-98, y no en esos dos
documentos. El origen del alelo y el del locus deben concordar.»

**REQ-MAN-96.** *Hoy:* «Un locus generado debe publicar además sus
parámetros, la versión de su algoritmo, las de los algoritmos de azar que
consume, el nombre de su ámbito de azar y la versión del catálogo de ámbitos,
las tres limitaciones declaradas del modelo palabra por palabra, y la
contabilidad de la generación: linajes, pasos sorteados, pasos aplicados y
pasos no aplicados al detenerse en un borde. Los aplicados más los no
aplicados deben sumar exactamente los sorteados.»
*Propuesto:* el mismo texto, y al final: «Para los loci **detectados**, lo que
es común a todos —los parámetros del deslizamiento común, las versiones, el
ámbito y la versión de su catálogo, los linajes y las tres limitaciones
palabra por palabra— se publica **una sola vez**, en la sección de
microsatélites de la configuración efectiva, y la contabilidad de cada uno
—pasos sorteados, aplicados y no aplicados— en la tabla de alelos de
REQ-MAN-98. La suma se exige locus a locus.»

**REQ-MAN-98.** *Hoy* termina en «El catálogo viaja **aparte** […]. La verdad y
el VCF no cambian de formato.» *Propuesto:* antes de la última frase, añadir:
«**Los alelos de los loci detectados** viajan también aparte, en una tabla con
una fila por alelo, enlazada con huella desde la configuración efectiva y
listada en el manifiesto general: su secuencia, su origen, sus frecuencias
declarada y efectiva, sus copias y su frecuencia realizadas y, por locus, la
contabilidad del deslizamiento y sus bases de tracto. La configuración efectiva
y el manifiesto general no publican nada propio de cada locus detectado, de
modo que no crecen con su número; los loci declarados siguen publicándose
enteros en ellos.»

**REQ-VAL-60.** *Hoy:* «Un juez independiente debe **rehacer el abanico
publicado** —desde la referencia publicada, los parámetros publicados y la
semilla publicada— y contrastarlo con el escrito, alelo a alelo. No basta con
reconstruir la configuración llamando al productor: eso compara el abanico
consigo mismo y no puede fallar.» *Propuesto:* el mismo texto, y tras «alelo a
alelo»: «El de un locus detectado lo rehace desde su fila del catálogo y el
deslizamiento común de la sección de microsatélites publicados, y lo contrasta
fila a fila con la tabla de alelos, incluida en cada fila la contabilidad del
locus.»

**REQ-VAL-62.** *Hoy* (fragmentos que cambian): «En la configuración efectiva,
los loci —con sus alelos, declarados o generados— y los mínimos de admisión
deben coincidir **exactamente** con los que resuelve el TOML copiado sobre la
referencia publicada. En el manifiesto general, la sección de loci debe
coincidir exactamente con la que resuelve el TOML **salvo lo medido** […]; y
su parte VCF, en total y por locus, debe ser la que lleva el VCF publicado […].
Las copias realizadas de cada alelo deben ser las que lleva la verdad
publicada […].»
*Propuesto:* «En la configuración efectiva, los loci declarados —con sus
alelos, declarados o generados— y los mínimos de admisión, y en la tabla de
alelos los alelos de los loci detectados, deben coincidir **exactamente** con
los que resuelve el TOML copiado sobre la referencia publicada, salvo, en la
tabla, las columnas medidas: copias y frecuencia realizadas y bases de tracto
del VCF. En el manifiesto general, la sección de loci —la de los declarados—
debe coincidir exactamente con la que resuelve el TOML **salvo lo medido** […];
y su parte VCF, en total y por locus —los declarados en la configuración
efectiva y los detectados en la tabla de alelos—, debe ser la que lleva el VCF
publicado, atribuyendo cada registro de repetición al primer locus del
catálogo que lo contiene […]. Las copias realizadas de cada alelo —en el
manifiesto las de los declarados y en la tabla las de los detectados— deben
ser las que lleva la verdad publicada […].» El resto, sin cambios.

**REQ-VAL-64.** *Hoy* termina en «[…] y que la tabla dice lo mismo que
`repeat_loci` de la configuración efectiva (`src/repeat_catalog_validation.rs`,
01/11/2026).» *Propuesto:* «[…] que las filas declaradas de la tabla son los
loci de `repeat_loci` de la configuración efectiva y las detectadas, los loci
de la tabla de alelos, en el mismo orden. De la tabla de alelos exige además
su forma, que las copias realizadas de cada locus suman las copias
haplotípicas publicadas, que cada frecuencia realizada sale de sus copias y
que los pasos aplicados y no aplicados suman los sorteados.»

**REQ-MAN-29**, nota fechada nueva: «Con la sección de microsatélites, los loci
**detectados** se registran en el catálogo y en la tabla de alelos, enlazados
desde la configuración efectiva y listados con su huella en el manifiesto
general (REQ-MAN-98), y no dentro del manifiesto (decisión del responsable del
29/10/2026; DTC-035).»

**DTC-035**, nota fechada nueva: «Excepción a "V1 no fragmenta manifiestos en
JSONL ni tablas auxiliares": lo que crece con los loci **detectados** sobre el
genoma —el catálogo de microsatélites y los alelos de los detectados— viaja en
tablas propias, enlazadas con huella desde la configuración efectiva y listadas
con tamaño y huella en el manifiesto general (REQ-MAN-98, DTC-104). Decidido
por el responsable en la pausa de la entrega 1 del plan de detección
(29/10/2026, DP-SSR-13 b): con 56 736 loci los dos documentos pasaban de medio
gigabyte.»

**DTC-104.** *Hoy:* «**Por ahora los loci siguen además en `repeat_loci`** de la
efectiva y del manifiesto, con sus abanicos: sacar los detectados de esos dos
documentos y llevar sus abanicos a una tabla propia es la entrega 4 bis del
plan, con su cambio normativo previo (REQ-MAN-95, REQ-MAN-96, REQ-VAL-60 y
REQ-VAL-62).» *Propuesto, en su lugar:* «**Los alelos de los detectados** van
en `repeat_alleles.tsv`, una fila por alelo en el orden del catálogo y del
abanico, con las columnas `locus_id`, `allele_index`, `sequence`, `origin`,
`declared_frequency`, `effective_frequency`, `realized_copies`,
`realized_frequency`, `drawn_events`, `applied_events`, `unapplied_events`,
`truth_tract_bases_bound` y `vcf_tract_bases_bound` (las cinco últimas, del
locus, repetidas en cada fila suya); se escribe siempre con la sección,
después de la unidad genotípica y antes de la efectiva, que la enlaza con
huella en `repeat_detection.alleles_table`, y figura en el manifiesto general.
La anotación de cada alelo no se publica: se deriva de su secuencia y del
motivo del catálogo. Lo común del deslizamiento —parámetros, versiones, ámbito,
catálogo de ámbitos, linajes y las tres limitaciones— va una vez en
`repeat_detection.slippage`, y las constantes de la normalización en
`repeat_detection.frequency_normalization`. `repeat_loci` y
`repeat_output_budget.per_locus` de la efectiva, y `repeat_loci.loci` del
manifiesto, llevan solo los declarados. En un paquete fallido, el par del
catálogo y la tabla de alelos forman parte de la unidad genotípica.» Y en el
párrafo del juez, añadir la ampliación de su lista (REQ-VAL-64) y que el
abanico de cada detectado lo rehace el juez de REQ-VAL-60 desde la tabla.

**DTC-097**, nota fechada nueva: «Para los loci detectados (DTC-104), los
alelos con su origen, sus frecuencias y lo realizado, la contabilidad del
bloque `slippage` y la parte de la cota de cada locus viajan en la tabla de
alelos, y lo fijo del bloque —modelo, tamaño de paso, media, sesgo, versiones,
ámbito, versión del catálogo de ámbitos y limitaciones— y los linajes, una vez
en `repeat_detection.slippage`. Lo que esta decisión exige de cada locus,
incluida la parte VCF locus a locus, vale igual en la tabla.»

**DTC-095**, ampliación del 28/10/2026: «las tres limitaciones declaradas, que
hoy viajan en cada locus generado, se publicarían **una sola vez** si el
catálogo se publica aparte (DTC-104)» pasa a «[…] se publican **una sola vez**
para los detectados, en la sección de microsatélites de la efectiva (DTC-104);
un declarado que genera sus alelos las sigue llevando en su bloque».

**Ampliación de DTC-045** (qué conserva un paquete fallido), en «la unidad
genotípica»: añadir «con la sección de microsatélites, también el par del
catálogo y la tabla de alelos», y cambiar «Con loci de repetición,
`repeat_loci` publica los sorteos […]» por «Con loci de repetición, `repeat_loci`
publica los de los **declarados**, y la tabla de alelos conservada los de los
detectados, con los sorteos y las bases de tracto que la verdad y el VCF
conservados llevan de verdad; sin unidad, `repeat_loci` registra los declarados
con cero copias sorteadas y cifras realizadas a cero».

**Arquitectura, apartado 16.4**, viñeta «Lo publicado». *Hoy:* «los loci con sus
alelos y su origen por alelo (REQ-MAN-95, REQ-MAN-96) en la configuración
efectiva y el manifiesto general; […]». *Propuesto:* «los loci declarados con
sus alelos y su origen por alelo (REQ-MAN-95, REQ-MAN-96) en la configuración
efectiva y el manifiesto general; el catálogo de todos los loci y los alelos
de los detectados en tablas aparte enlazadas con huella (REQ-MAN-98), de modo
que los documentos JSON no crecen con los detectados; […]».

##### 9. Pesos: medidos sobre lo publicado y extrapolados

**Cómo se midió.** Binario de `afda48b` compilado con `--release` en un
directorio de compilación propio. Dos referencias sintéticas de una secuencia
con 250 y 1 000 microsatélites plantados (motivos de 1 a 6 pb al azar, entre 2
y 12 copias por encima del suelo, uno cada 400 bases sobre un fondo sin
tiradas de periodo 1 ni 2); el detector halló 251 y 1 006 loci. Una población
diploide HWE de 10 y de 100 individuos, deslizamiento común con media 1, VCF de
texto y una matriz de una fila fuera de todo tracto, como en el apartado 11.4.
La forma propuesta se obtuvo **transformando los JSON publicados** con un
guion propio (quitar los detectados, añadir lo común una vez, el enlace y la
entrada del manifiesto) y escribiendo la tabla con sus columnas desde los
mismos datos; la reescritura de los JSON reproduce el tamaño del fichero
publicado con un error del 0,02 %. Las pendientes son la diferencia entre los
dos tamaños dividida por los 755 loci de diferencia.

| Fichero | 10 individuos | 100 individuos | Con 56 736 detectados (10 / 100) |
|---|---:|---:|---:|
| Efectiva hoy | 5 006 B/locus | 6 558 B/locus | 284 MB / 372 MB |
| Manifiesto hoy | 5 256 B/locus | 7 116 B/locus | 298 MB / 404 MB |
| **Efectiva propuesta** | 7,4 kB fijos | 8,8 kB fijos | **7,4 kB / 8,8 kB** |
| **Manifiesto propuesto** | 8,4 kB fijos | 9,8 kB fijos | **8,4 kB / 9,8 kB** |
| Tabla de alelos, solo secuencia (A2 a) | 526 B/locus | 848 B/locus | 30 MB / 48 MB |
| … con las unidades (A2 b) | 558 B/locus | 896 B/locus | 32 MB / 51 MB |
| … con toda la anotación (A2 c) | 657 B/locus | 1 048 B/locus | 37 MB / 59 MB |
| `repeat_catalog.tsv` (sin cambio) | 60 B/locus | 60 B/locus | 3,4 MB |

Los dos documentos propuestos cambian **3 y 7 bytes** entre 251 y 1 006 loci
(los dígitos de los totales): ya no crecen con los detectados; lo que queda
crece con los individuos (la lista de identificadores) y con los declarados.
Hoy los dos suman unos 10,3 kB por locus con 10 individuos, algo más que los
9,5 kB del apartado 11.6 porque los tractos sintéticos son algo más largos; con
los pesos de 11.6 la tabla rondaría el 5 % de lo que hoy ocupan los dos JSON.
La tabla la domina la secuencia de cada alelo (43 bases de media aquí; 4,6
alelos por locus con 10 individuos y 7,0 con 100). Con la forma propuesta, la
efectiva cae además por debajo de los 2 MB que la interfaz Marimo acepta para
sus resúmenes; hoy los supera a partir de unos 400 loci.

**Fuera de alcance, medido en 11.6:** el manifiesto del evento técnico también
crece con los loci (226 B por locus con tartamudeo, unos 13 MB con 56 736:
parámetros y pasos por locus); la verdad, con locus por individuo (DP-SSR-13:
su formato no cambia).

##### 10. Consumidores que hoy leen `repeat_loci` de lo publicado

Buscado con `grep` en `src/`, `tests/`, `ui/`, `tests_ui/` y `scripts/`.

- **Cambian** (`src/results.rs`): el juez del abanico
  (`validate_published_slippage_fans`, no estaba en la lista de seis); el de
  los loci frente al TOML (`validate_published_repeat_loci`); el de la
  contabilidad de tracto (`judge_tract_output_accounting`); el de las cotas
  frente al TOML (`validate_published_tract_bounds`, no estaba en la lista);
  el de la efectiva frente al manifiesto
  (`validate_effective_config_v3_against_manifest_bytes`: hoy exige «efectiva
  con loci si y solo si manifiesto v6», y con solo detectados la efectiva no
  tendría loci); el de los objetivos frente al TOML
  (`validate_effective_targets_against_configuration`, por el catálogo vacío
  del apartado 7); la forma de la sección del manifiesto
  (`validate_repeat_loci_section`, que hoy rechaza una lista de loci vacía);
  y el orquestador `validate_published_resources_v3`. En
  `src/repeat_catalog_validation.rs`, la comparación con `repeat_loci`.
  Escritura: `repeat_loci_json`, `repeat_loci_manifest_json`,
  `repeat_output_budget_json`, `repeat_detection_json` y el escritor nuevo, en
  `src/results.rs`; el orden de escritura y el paquete fallido, en
  `src/execution.rs`.
- **No cambian, aunque la lista de seis los nombraba**: el de la referencia
  simulada (`published_simulated_tracts`) y el de la referencia final
  (`judge_inserted_reference`) ya saltan los detectados; el de los nombres
  (`validate_published_reference_names`) solo lee la referencia de la sección,
  que se queda.
- **La verdad, el VCF y el evento técnico** no dependen de que los abanicos
  estén en la efectiva: sus jueces leen la verdad y el VCF publicados o rehacen
  desde el TOML, y ninguno lee `repeat_loci` publicado.
- **Pruebas**: la batería de parejas finales (`tests/final_pairs_oracle/mod.rs`
  lee los intervalos de `repeat_loci` de la efectiva y debe leer los de los
  detectados del catálogo), `tests/repeat_detection_run.rs` (su ayudante
  `catalog`), `tests/repeat_catalog_validation.rs`, las cuatro huellas con la
  sección de `tests/repeat_detection_frozen_identity.rs` y las de
  `tests/technical_catalog_contract.rs` (cambian a propósito). Las demás
  pruebas que leen `repeat_loci` usan el ayudante que deja el catálogo igual a
  los declarados (`tests/support/detection_section.rs`), y no cambian.
- **Interfaz Marimo: no se rompe nada.** El editor solo escribe
  `[[reference.repeat_loci]]` en el TOML de entrada, que no cambia, y los
  resúmenes de resultados (`ui/process_adapter.py`) no leen `repeat_loci`;
  como solo abren documentos de hasta 2 MB, con la forma propuesta vuelven a
  cargarse en ejecuciones con muchos loci. La entrega no añade ninguna opción
  al TOML, así que no hay control nuevo que exigir (DP-SSR-16).
- **Esquemas**: la efectiva v6 y el manifiesto v6 se retocan in situ por la
  nota alfa de DTC-036 (`repeat_detection` gana campos; la lista de loci del
  manifiesto admite estar vacía con la sección). La tabla, como la del catálogo
  y la de la matriz, la juzga el código, sin esquema.
- **Puerta de cobertura**: cambian las cifras del juez del catálogo y de
  `src/results.rs` en `scripts/guard-coverage-baseline.tsv`; no entra ningún
  módulo nuevo.

##### 11. Decisiones (abiertas al escribir la propuesta; decididas el 01/11/2026)

**El responsable decidió las tres el 01/11/2026 según la recomendación:
A1 (a), A2 (a) y A3 (a).**

**A1. ¿Se amplía el cambio normativo a DTC-035 y REQ-MAN-29?**

| Opción | Consecuencia |
|---|---|
| **(a) Sí, con las notas fechadas del apartado 8** | La norma general recoge la excepción que el 29/10 ya se aceptó; los dos JSON quedan en unos 8 kB |
| (b) No: el manifiesto conserva los alelos y lo realizado de los detectados | Solo adelgaza la efectiva; el manifiesto sigue creciendo unos 5,3 kB por locus (300 MB con 56 736), que es lo que 11.7 juzgó inviable |

**Recomendación: (a). Decidido (a) por el responsable el 01/11/2026.**

**A2. ¿Se publica la anotación de cada alelo de un detectado?**

| Opción | Consecuencia |
|---|---|
| **(a) No: solo la secuencia** | Lo más pequeño; la arquitectura ya declara la anotación salida derivada; quien quiera las unidades divide la longitud entre la del motivo |
| (b) Solo las unidades | Un 6 % más de tabla; ahorra una división a quien lea; una columna más que juzgar |
| (c) Toda, como hoy (unidades, pureza, tirada, bloques, interrupciones) | Un 25 % más; el juez tendría que rehacer la descomposición con aritmética propia, sin la del productor |

**Recomendación: (a). Decidido (a) por el responsable el 01/11/2026.**

**A3. ¿Dónde va la contabilidad de cada locus (pasos y bases de tracto)?**

| Opción | Consecuencia |
|---|---|
| **(a) Columnas de la tabla de alelos, repetidas en cada fila del locus** | Un solo fichero nuevo; unos 25 bytes por fila; el juez del abanico compara cada fila con lo rehecho y eso ya obliga a que coincidan |
| (b) Una tabla aparte con una fila por detectado | Sin repetición, pero un fichero, un enlace y una huella más |
| (c) Columnas nuevas en `repeat_catalog.tsv`, vacías para los declarados | Mezcla el resultado de la generación con un índice que se juzga contra el FASTA, y obliga a escribir el catálogo después del VCF |

**Recomendación: (a). Decidido (a) por el responsable el 01/11/2026.**

**Decidido aquí por ser mecánico o venir fijado** (el responsable puede
reabrir cualquiera): los declarados se quedan en `repeat_loci`; nombre,
artefacto, clase y sitio de la tabla; enlace `alleles_table` sin metadato
propio; tabla siempre con la sección; orden y formato de las columnas;
`drawn_copies` y las sumas de la normalización no se publican por detectado;
linajes y constantes de la normalización una vez en `repeat_detection`; la
sección de loci del manifiesto existe mientras el catálogo tenga algún locus y
su lista puede quedar vacía; el par del catálogo y la tabla, dentro de la unidad
genotípica del paquete fallido; el manifiesto del evento técnico no se toca en
esta entrega.

> **Revisión de la entrega 4 bis:** los JSON dejan de crecer con los loci y
> los jueces siguen sin llamar al productor.

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
- **C-9. DTC-035 y REQ-MAN-29 frente a DP-SSR-13 (b).** DTC-035 dice que un
  manifiesto contiene directamente todo lo que debe registrar y que no se
  fragmenta en tablas auxiliares; REQ-MAN-29, que el manifiesto registra los
  loci y lo realizado. La entrega 4 bis saca los detectados a una tabla, como
  aceptó el responsable el 29/10/2026. **Resuelta el 01/11/2026:** el
  responsable decidió A1 (a), y las dos normas llevan desde ese día una nota
  fechada que recoge la excepción (propuesta normativa de la entrega 4 bis).

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

Los apartados 11.2 a 11.10 se midieron con el núcleo copiado de `pop_var_caller`
en el commit `950909f`. Desde el 30/10/2026 el núcleo procede de `fa215c23`,
que ya no funde dos repeticiones pegadas de motivos distintos; sus cifras están
en el apartado 11.11.

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



### 11.10. Entrega 2, segunda parte: el punto 5 y los tres crecimientos sin explicar (29/10/2026)

Medido el 29/10/2026 por el implementador de los pasos F6.4.a.8f a
F6.4.a.8j, que el responsable añadió a la entrega 2 ese mismo día. Mismos
comandos que el apartado 11.1 (`measure_publication`), mismo FASTA, misma
máquina y **los mismos documentos** que 11.9 (preparados con el detector de
hoy; el de 16 000 loci y 100 individuos se preparó aparte con los mismos
16 000 loci). `rustc 1.97.1`.

#### 11.10.1. Dónde estaba cada coste, medido

Cronómetros acumulados por función, en una copia instrumentada aparte de
`ed30469` (el código de 11.9), con 8 000, 16 000 y 50 284 loci y 10
individuos; una pasada cada una, con la máquina en reposo. Segundos de
pared; entre paréntesis, el número de llamadas con 50 284 loci. El
crecimiento esperado si fuera lineal es ×2 de 8 000 a 16 000 y ×3,14 de
16 000 a 50 284.

| Función (fichero) | 8 000 | 16 000 | 50 284 | ×2 | ×3,14 |
|---|---:|---:|---:|---:|---:|
| **Aplicar los alelos** (`apply_repeat_alleles`, `src/repeat_genotypes.rs`; 40) | 13,8 | 54,5 | 545,8 | ×3,96 | **×10,0** |
| — copiar la cadena por cada locus (`with_current_repeat_allele`, `src/segments.rs`; 2,2 millones) | 7,3 | 30,5 | 319,5 | ×4,2 | ×10,5 |
| — comprobar que el tramo llega intacto (`intact_repeat_interval`; 6,2 millones) | 1,1 | 4,4 | 42,4 | ×4,0 | ×9,7 |
| de ella, en los dos preflights (`preflight_initial_individuals`; 2) | 8,1 | 28,8 | 273,6 | ×3,6 | ×9,5 |
| y en el juez de la población inicial (`validate_physical_initial_population_statistics`; 2) | 8,1 | 28,4 | 275,2 | ×3,5 | ×9,7 |
| **Proyectar cada copia en la verdad** (`project_copy`, `src/truth_projection.rs`; 260) | 2,0 | 8,4 | 89,3 | ×4,3 | **×10,6** |
| Buscar el alelo de cada caso (`project_repeat_genotypes`; 1) | 0,16 | 0,41 | 2,9 | ×2,7 | ×7,0 |
| **Bases del VCF por locus en la efectiva** (`publish_measured_vcf_tract_bases`, `src/results.rs`; 3) | 0,33 | 1,4 | 16,1 | ×4,1 | **×11,8** |
| Admisión contra variantes pequeñas (`admit_repeat_loci_against_small_variants`; 2) | 0,00 | 0,00 | 0,00 | — | — |
| Escribir y validar la verdad (`write_truth_artifact`; 1) | 15,9 | 32,3 | 100,9 | ×2,0 | ×3,1 |
| Ejecución entera | 62,9 | 145,3 | 834,6 | ×2,3 | ×5,7 |

Lo que esto dice de los tres crecimientos que 11.9 dejó sin explicar:

1. **Preflight ×12,1**: es el punto 5, como 11.9 suponía. Pero el coste no
   estaba en recorrer las piezas para buscar el tramo intacto (42 s), sino
   en **copiar la lista entera de fragmentos** de la copia por cada locus
   aplicado (320 s). Los unos 190 s de `apply_repeat_alleles` que quedan
   fuera de las funciones cronometradas no se midieron por separado; por
   lectura, son sobre todo soltar la lista anterior de fragmentos cada vez
   que se sustituye la copia, que también es proporcional a su longitud.
2. **Verdad ×6,0**: **no** era la búsqueda del alelo de cada caso que 11.5
   señalaba por lectura (`src/truth_projection.rs:1177`, 2,9 s), sino
   `require_digest`: al proyectar cada fragmento y cada evento de cada copia
   se comprobaba que su procedencia tuviera huella **recorriendo la lista de
   procedencias**, que tiene una entrada por locus aplicado (89 s). Se
   arreglaron las dos (F6.4.a.8g).
3. **Efectiva y manifiesto ×10,6**: `publish_measured_vcf_tract_bases`, que
   para escribir en la configuración efectiva las bases del VCF de cada locus
   las buscaba en la lista de todos (16 s). No es memoria ni la escritura de
   los 480 MB de JSON: el resto de la etapa crece ×3,1 a ×3,5.
4. **La admisión contra variantes pequeñas no pesa** con la matriz de una
   fila de la medición (menos de 5 ms). Se arregló porque es de la misma
   clase y crece con loci × filas (F6.4.a.8h), no por lo medido.

Con esos cuatro arreglos, una segunda copia instrumentada (del árbol ya
arreglado) mostró un quinto coste de la misma clase, más pequeño:
`validate_repeat_genotype` (`src/truth_validation.rs`), el juez de cada caso
de repetición, crecía ×9,9 de 16 000 a 50 284 loci (3,2 → 31,3 s) porque
recorría todas las copias de todos los individuos y el historial entero de
cada una. Se arregló también (F6.4.a.8j).

#### 11.10.2. Qué cambió en cada paso

| Paso | Ficheros | Qué recorría antes | Qué consulta ahora |
|---|---|---|---|
| F6.4.a.8f | `src/segments.rs`, `src/repeat_genotypes.rs` | por locus y copia, los fragmentos de la copia hasta el suyo, y copiaba la lista entera | los fragmentos de fuente de la copia por búsqueda binaria; la cadena final se construye una vez por copia |
| F6.4.a.8g | `src/truth_projection.rs` | por caso, los fragmentos de cada copia; por fragmento y evento, todas las procedencias | una tabla de alelos por copia y secuencia; las procedencias en un conjunto |
| F6.4.a.8h | `src/execution.rs` | por locus, todos los cambios directos y todas las filas de la matriz | los que caen en su intervalo, en mapas por secuencia y posición |
| F6.4.a.8i | `src/results.rs` | por locus, la lista de bases del VCF de todos | un mapa por locus |
| F6.4.a.8j | `src/truth_validation.rs` | por caso, todas las copias de todos los individuos y el historial entero de cada una | las copias de ese individuo y esa secuencia; los alelos del historial con ese intervalo |

En todos se conserva **el primero** que encontraba el recorrido de antes:
el primer fragmento que contiene el intervalo, el primer alelo del
intervalo, el primer cambio o fila en su orden, la primera entrada de cada
locus, el primer evento del historial que cumple. En F6.4.a.8f la
comprobación de cada locus sigue haciéndose en el mismo orden (locus a
locus, copia a copia), de modo que el primer rechazo es el mismo; y el lote
renuncia al índice y aplica los alelos uno a uno si la cadena no lo admite
(ver la documentación de `RepeatAlleleBatch`). F6.4.a.8j no añade ni quita
ninguna comprobación del juez: la puerta de cobertura sale con las mismas
cifras (121 comprobaciones y 32 sin disparar en `src/truth_validation.rs`) y
la línea base no se tocó.

#### 11.10.3. Huella y mutaciones

**La huella.** Las dos congeladas siguen pasando sin tocar un literal:
`tests/repeat_scale_frozen_identity.rs` (tres pruebas) y
`tests/repeat_detection_frozen_identity.rs` (cuatro). Cubren los caminos
tocados así: los escenarios 1 y 2 de la de escala llevan quince loci en dos
secuencias **y una matriz de tres filas** entre ellos (los cambios de base
de la matriz caen dentro de fragmentos que los alelos parten, que es lo que
el lote reparte), VCF y configuración efectiva con las bases por locus; los
veintiún escenarios de la revisión de la fusión del VCF llevan indels y
tractos, de modo que hay copias con varios fragmentos de fuente. **No hacía
falta un escenario nuevo**: la matriz con microsatélites ya está; lo que la
huella no puede ver de la admisión son sus **rechazos**, que no se publican,
y para eso está la regresión de F6.4.a.8h. Esa regresión se ejecutó también
contra el código de `ed30469` (extraído con `git archive`, con su propio
directorio de compilación) y **pasa igual**: el diagnóstico es el de antes.

**Mutaciones**, todas sobre una copia aparte del árbol (extraída con
`git ls-files -co --exclude-standard`, sin `.vscode/`), con su propio
`CARGO_TARGET_DIR`, la copia de seguridad fuera de git y el SHA-256 del
fichero antes, mutado y restaurado (doce primeras cifras). Las de
`src/truth_projection.rs`, `src/execution.rs` y `src/results.rs` se hicieron
antes de renumerar en los comentarios los pasos F6.4.a.8g a F6.4.a.8i, por
eso su SHA difiere del del árbol; el código es el mismo. «Huella: revisión»,
«base» y «clonal» son las tres pruebas de la huella de escala; «huella F6»
son las cuatro de `tests/repeat_detection_frozen_identity.rs`. Cuando una
batería de pruebas muere, `cargo test` no ejecuta las siguientes: «solo»
quiere decir que las demás que se ejecutaron pasaron.

| Mutación | Paso | Fichero | SHA antes / mutado / restaurado | Lo que muere |
|---|---|---|---|---|
| El índice busca el fragmento con `<` y pierde el que empieza en el inicio del locus | F6.4.a.8f | `src/segments.rs` | `63bbfef94552` / `b277841bdd93` / `63bbfef94552` | la prueba unitaria; huella: revisión |
| La búsqueda de cambios de base con `<=` pierde el de la primera base | F6.4.a.8f | `src/segments.rs` | `63bbfef94552` / `1e810714e32c` / `63bbfef94552` | **solo** la prueba unitaria |
| El lote pierde los cambios de base de los restos anteriores a cada alelo | F6.4.a.8f | `src/segments.rs` | `63bbfef94552` / `2674527fadba` / `63bbfef94552` | la prueba unitaria; huella: revisión, base y clonal |
| El lote escribe los eventos en orden de posición y no de aplicación | F6.4.a.8f | `src/segments.rs` | `63bbfef94552` / `cb33c153a1b4` / `63bbfef94552` | **solo** la prueba unitaria |
| La tabla de alelos de la verdad no se rehace al cambiar de secuencia | F6.4.a.8g | `src/truth_projection.rs` | `85a1d8aaf8ce` / `0a277d040b86` / `85a1d8aaf8ce` | huella: base y clonal |
| El conjunto de procedencias pierde la primera | F6.4.a.8g | `src/truth_projection.rs` | `85a1d8aaf8ce` / `5debf179a5c7` / `85a1d8aaf8ce` | huella F6 (las cuatro) |
| El rango de la admisión empieza en la base siguiente al inicio del locus | F6.4.a.8h | `src/execution.rs` | `2f0dc7d315f1` / `205b6fb76e47` / `2f0dc7d315f1` | **solo** la regresión de F6.4.a.8h |
| La admisión nombra el primero por posición y no el primero declarado | F6.4.a.8h | `src/execution.rs` | `2f0dc7d315f1` / `6bb386d7ca12` / `2f0dc7d315f1` | **solo** la regresión de F6.4.a.8h |
| El mapa de bases del VCF pierde la primera entrada | F6.4.a.8i | `src/results.rs` | `cd82b57f2b7a` / `0af9761921b9` / `cd82b57f2b7a` | huella: revisión, base y clonal |
| El tramo de copias del juez acaba antes de la segunda copia | F6.4.a.8j | `src/truth_validation.rs` | `7c99bc18f6ae` / `d8a95d114642` / `7c99bc18f6ae` | huella: revisión, base y clonal; 27 pruebas de `tests/truth_validation.rs` y `tests/truth_validation_guards.rs` |
| El juez busca el alelo por un intervalo que acaba una base después | F6.4.a.8j | `src/truth_validation.rs` | `7c99bc18f6ae` / `2977ef681b1a` / `7c99bc18f6ae` | las mismas |

Lo que esto dice, sin adornos: la huella **no** ve que el lote pierda un
cambio de base en la primera base de un locus (ningún escenario lo tiene:
la admisión lo rechaza antes) ni el orden de los eventos (en los escenarios
los loci se declaran en orden de posición); los dos los caza la prueba
unitaria. Los rechazos de la admisión no se publican: solo los caza su
regresión. Que el conjunto de procedencias o el mapa de bases del VCF
conserven **la primera** entrada repetida no se puede comprobar: los
productores no repiten ni procedencias (se deduplican antes) ni loci.

#### 11.10.4. Tiempos antes y después

Segundos de pared por etapa, como en 11.9. «Antes (11.9)» es la columna
«después» de 11.9 (dos pasadas, media; la de 50 284 loci, una). «Antes
(`ed30469`)» es el mismo código medido otra vez ahora, una pasada, para
tener un antes en las mismas condiciones. «Después» es el binario de pruebas
optimizado de este árbol (SHA-256 `07c9bb90a4a6…`), **dos pasadas**
separadas por una barra. La máquina no tenía nada más en marcha (comprobado
con `ps` antes de cada ejecución); la campaña corrió en segundo plano y se
esperó a su notificación.

| Loci | Indiv. | Binario | Lectura TOML | Preflight | Preflight interno | Verdad | VCF | Efectiva y manifiesto | Validación final | Total | Memoria máx. (MiB) |
|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 000 | 10 | antes (11.9) | 0,07 | 4,6 | 4,8 | 21,6 | 2,6 | 1,8 | 34,9 | 70,2 | 1 801 |
| 8 000 | 10 | después | 0,07 / 0,08 | 1,1 / 1,5 | 1,5 / 2,2 | 22,1 / 26,2 | 2,7 / 2,6 | 1,9 / 1,9 | 33,5 / 31,7 | 62,9 / 66,1 | 1 806 |
| 16 000 | 10 | antes (11.9) | 0,15 | 16,6 | 16,7 | 48,5 | 4,3 | 4,2 | 85,5 | 175,8 | 3 512 |
| 16 000 | 10 | antes (`ed30469`) | 0,17 | 22,3 | 20,7 | 55,5 | 5,4 | 5,8 | 96,0 | 205,8 | 3 499 |
| 16 000 | 10 | después | 0,16 / 0,16 | 1,9 / 1,9 | 2,7 / 2,3 | 43,4 / 42,6 | 4,2 / 4,4 | 4,0 / 4,4 | 62,6 / 58,2 | 119,0 / 114,0 | 3 511 |
| 50 284 | 10 | antes (11.9) | 0,61 | 201,4 | 205,6 | 290,1 | 21,7 | 44,5 | 673,3 | 1 437,2 | 10 801 |
| 50 284 | 10 | antes (`ed30469`) | 0,53 | 178,3 | 190,7 | 243,8 | 21,3 | 25,4 | 530,1 | 1 190,2 | 11 294 |
| 50 284 | 10 | después | 0,58 / 0,55 | 4,8 / 3,6 | 4,7 / 4,1 | 131,5 / 131,1 | 12,0 / 12,1 | 15,1 / 12,0 | 211,3 / 186,8 | 380,2 / 350,3 | 10 795 |
| 4 000 | 100 | antes (11.9) | 0,04 | 16,5 | 18,1 | 129,1 | 29,8 | 2,0 | 126,7 | 322,0 | 2 818 |
| 4 000 | 100 | después | 0,04 / 0,04 | 11,4 / 10,3 | 10,8 / 10,3 | 124,6 / 139,5 | 19,8 / 21,9 | 2,3 / 2,3 | 96,2 / 107,1 | 265,1 / 291,5 | 2 869 |
| 16 000 | 100 | después | 0,16 / 0,17 | 14,2 / 15,3 | 14,6 / 15,4 | 424,6 / 460,5 | 47,4 / 49,2 | 9,2 / 9,9 | 307,2 / 319,0 | 817,3 / 869,4 | 11 022 |

El catálogo entero del tomate reducido con 10 individuos pasa de 1 190 a
1 437 s (según cuál de los dos «antes») a **350–380 s**, unos seis minutos,
con la misma memoria. Las dos mediciones del «antes» difieren entre sí un
17 % con 16 000 loci y un 21 % con 50 284, en direcciones opuestas: la
máquina no da más precisión que esa entre sesiones, y la de 50 284 de 11.9
fue una sola ejecución. Ninguna de las conclusiones de abajo depende de
cuál se tome.

**Crecimiento de cada etapa de 16 000 a 50 284 loci** (×3,14; media de las
dos pasadas, 10 individuos). El criterio de la entrega es ≈ ×3,1:

| Etapa | Antes (11.9) | Después | Lineal |
|---|---:|---:|---|
| Preflight | ×12,1 | ×2,2 | sí (lo que queda son costes fijos: 4 s) |
| Preflight interno | ×12,3 | ×1,8 | sí (ídem) |
| Verdad | ×6,0 | ×3,05 | **sí** |
| VCF | ×5,1 | ×2,8 | sí |
| Efectiva y manifiesto | ×10,6 | ×3,2 | **sí** |
| Validación final | ×7,9 | ×3,3 | casi: ver 11.10.6 |
| Total | ×8,2 | ×3,1 | sí |

Y con 100 individuos, de 4 000 a 16 000 loci (×4): verdad ×3,4, VCF ×2,3,
validación final ×3,1, total ×3,0.

#### 11.10.5. Extrapolación al catálogo del tomate reducido

| Catálogo de 50 284 loci | 11.9 | Ahora |
|---|---|---|
| 10 individuos, sin tartamudeo | medido: 1 437 s, 10 801 MiB | **medido: 350–380 s, 10 795 MiB** |
| 10 individuos, con evento técnico a 0,1× | unos 1 530 s | unos 430–475 s (lo anterior más los 80–95 s del evento, que no dependen de los loci; extrapolado) |
| 100 individuos, sin tartamudeo | unas 3,6 a 3,8 horas; unos 33 800 MiB | **unos 41 a 44 minutos**; unos **34 300 MiB** (extrapolado) |

Cómo se extrapola, para que se pueda rehacer:

- **Tiempo con 100 individuos**, dos cuentas: (a) en línea recta por los
  loci desde las medias medidas con 4 000 y 16 000 loci (278,3 y 843,4 s:
  0,0471 s por locus) da 2 458 s; (b) el crecimiento medido con 10
  individuos de 16 000 a 50 284 loci (×3,135) aplicado a los 843,4 s da
  2 644 s. Ya no hay término cuadrático que ajustar: 11.9 avisaba de que
  su ajuste se quedaba corto porque lo había; ahora el crecimiento con 10
  individuos es el de los loci.
- **Memoria**: en línea recta, 0,679 MiB por locus con 100 individuos (de
  2 869 a 11 022 MiB entre 4 000 y 16 000 loci). **Con 100 individuos el
  catálogo entero sigue sin caber en esta máquina** (27 GB): la memoria no
  era objeto de estos pasos. La bajan el catálogo aparte de la entrega 4 y,
  si hace falta, DP-SSR-18 (11.7).

#### 11.10.6. Lo que queda

1. **La validación final crece ×3,3** para ×3,14 loci. Medido por función en
   la copia instrumentada del árbol arreglado, lo que crece algo más que los
   loci es de otra clase o pesa poco: `json_document_bytes` ×6,2 (0,3 →
   1,7 s; árboles JSON de 480 MB), `repeat_loci_manifest_json` ×4,1 (6 →
   25 s en cinco llamadas; recorre los sorteos agrupados, sin búsquedas: es
   la construcción de los árboles JSON) y la lectura de la verdad
   publicada, ×3,4. Memoria y árboles JSON: no se tocaron.
2. **Dos recorridos de la misma clase en el juez de la verdad, pequeños, sin
   arreglar**: `expected_repeat_alleles` (`src/truth_validation.rs`, el
   juez del conjunto de casos) busca el alelo de cada caso recorriendo el
   historial de cada copia (0,4 → 2,6 s de 16 000 a 50 284 loci, ×6,9), y la
   réplica del historial (`replay_event`) busca el fragmento de cada alelo
   de repetición recorriendo los fragmentos y lo sustituye desplazando la
   lista (0,7 → 3,9 s, ×5,6). Entre los dos son 6,5 s de 374 con el catálogo
   entero; el primero se arregla con la tabla de F6.4.a.8j, y el segundo
   exige rehacer la réplica por lotes en un juez que es independiente a
   propósito. Se dejaron para no reabrir la campaña de medición por menos
   del 2 %.
3. **Sin medir, de la misma clase**: el juez de cada caso SNP
   (`validate_snp_genotype`) recorre todas las copias de todos los
   individuos y, en cada una, sus tramos de fuente, que crecen con los loci
   de la secuencia. Con la matriz de una fila de la medición no pesa; con
   una matriz grande y muchos loci, sí.
4. **Lo de 11.9.6 que sigue igual**: los dos recorridos cuadráticos del
   juez semántico del VCF, `copy_is_reference_at` y la memoria.

### 11.11. Núcleo nuevo de `pop_var_caller` (`fa215c23`): tramos pegados de motivos distintos (30/10/2026)

**Qué cambia.** El núcleo de `950909f` fundía en un solo intervalo dos
repeticiones pegadas de motivos distintos, con el motivo de la primera, y el
tramo fundido no pasaba la pureza: `A×12 T×12`, `(CA)×10 (GT)×10` y
`(TCTA)×8 (TCTG)×8 (TCTA)×8` no daban ningún candidato; `(GA)×12 C×12` daba
solo la tira de `C`, y `G×60 A×15` un solo candidato `G` de 75 bases que se
tragaba la tira de `A`. El núcleo de `fa215c23` parte cada segmento donde
cambia el motivo, vuelve a unir una interrupción corta cuando el motivo de
fuera reaparece y emite el tramo unido junto con sus trozos; la cadena de
`detect_repeat_candidates` juzga cada lectura por separado y la regla de
solapes (DP-SSR-05) se queda con una. Se conserva la desviación propia
decidida el 29/10 (apartado 11.8, tercera versión), aplicada ahora también al
trozo nacido de un corte: sin ella, el trozo empieza en la última base de la
tira previa (`A×8 (GA)×10` daría `[7, 28)`). Detalle y pruebas, en la cabecera
de `src/repeat_detection.rs`.

**Medición.** Mismos comandos y FASTA que el apartado 11.1; `rustc 1.97.1`.
Binario de pruebas optimizado de SHA-256 `85da6f0a…`, del árbol con
`src/repeat_detection.rs` de SHA-256 `de1c504a…` (la entrega 3, a medias, sin
*commit*). Dos pasadas con recuentos idénticos línea a línea. El «antes» es
la tercera versión del apartado 11.8; se volvió a medir en la misma sesión con
un binario del núcleo viejo compilado en una copia aparte (SHA-256
`0d2d8de2…`), y sus recuentos coinciden con 11.8 cifra a cifra. Suelo por
omisión `[8,6,6,6,5,4]`, antes → después:

| Etapa | p1 | p2 | p3 | p4 | p5 | p6 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| Intervalos del detector | 15 813 346 → 15 822 760 | 2 784 387 → 2 789 346 | 1 169 475 → 1 178 087 | 338 019 → 341 924 | 132 284 → 132 757 | 56 536 → 57 244 | 20 294 047 → **20 322 118** |
| Suelo de copias | 15 766 649 → 15 773 834 | 2 778 972 → 2 783 034 | 1 168 253 → 1 175 382 | 337 869 → 341 711 | 132 255 → 132 727 | 56 375 → 57 022 | 20 240 373 → 20 263 710 |
| Pureza | 2 133 → **27** | 744 → 117 | 315 → 826 | 59 → 49 | 6 → 5 | 80 → 86 | 3 337 → **1 110** |
| Tirada | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 3 | 2 → 3 |
| Tope | 0 → 0 | 0 → 0 | 2 → 10 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 10 |
| Borde | 1 → 1 | 0 → 0 | 1 → 3 | 0 → 0 | 0 → 0 | 0 → 0 | 2 → 4 |
| Fuera de regiones | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 |
| Solape | 7 → 306 | 20 → 552 | 9 → 200 | 10 → 28 | 0 → 0 | 1 → 9 | 47 → **1 095** |
| **Aceptados** | 44 556 → **48 592** | 4 651 → **5 643** | 895 → **1 666** | 81 → **136** | 23 → 25 | 78 → **124** | 50 284 → **56 186** |

Los intervalos del detector suben porque un tramo unido sale con sus trozos
(cada lectura se cuenta), y los solapes suben por lo mismo: el trozo que
pierde frente a su tramo unido, o al revés, se cuenta en «Solape». Con los
otros suelos: el bajo pasa de 261 122 a **274 286** aceptados (pureza 7 986 →
1 744; solapes 710 → 3 341) y el alto de 13 690 a **16 122** (pureza 1 445 →
919; solapes 11 → 515). Bases en candidatos con el suelo por omisión:
561 431 → **673 012** (0,717 % → 0,860 %).

Cifras auxiliares, suelo por omisión, antes → después: copia parcial final
2 921 → 4 065; pares de candidatos que se tocan 48 → **1 012**; pureza media
0,9955 → 0,9937; longitud en pb, mediana 9 → 9, percentil 90 15 → 17,
percentil 99 55 → 62, máxima 459 → 691. DP-SSR-02: 9 073 → 512 intervalos
pasan por longitud y no por copias completas, y 16 → 26 superarían además
pureza y tirada. DP-SSR-07: sin tope saldrían 56 175 (11 menos que con él:
los 10 tramos largos que el tope descarta ganarían solapes a varios trozos);
con más de 50 unidades, 50 → 102; de 100, 6 → 23; de 200, 1 → 7; de 500,
0 → 2; de 1 000, 0 → 1; máximo 422 → 1 262 unidades.

Barrido (`measure_scan`, misma sesión): núcleo viejo 12,0 y 12,1 s, núcleo
nuevo **7,8 y 6,9 s**; memoria residente máxima 229 952 kB → **181 724 y
181 704 kB**. El núcleo viejo midió 7,1–8,1 s el 29/10 (apartado 11.8) en la
misma máquina; la diferencia con los 12 s de hoy es de la máquina, no del
código (mismos recuentos), así que la comparación que vale es la de la misma
sesión. La mejora viene de la revisión de rendimiento del original
(`cdce162c`): la pila de Ruzzo y Tompa emite los segmentos en cuanto son
definitivos, y los segmentos demasiado cortos no se recorren.

**Cuántos tramos pegados recupera.** Comparando locus a locus los candidatos
de los dos núcleos (suelo por omisión): 50 125 comunes, **6 061 solo en el
nuevo** (4 110 de periodo 1, 1 062 de 2, 785 de 3, 56 de 4, 2 de 5 y 46 de 6)
y **159 solo en el viejo**. De los 6 061 nuevos, **1 920 tocan a otro
candidato** (1 904 homopolímeros, 10 de periodo 2 y 6 de periodo 3): son los
tramos pegados que el núcleo viejo fundía. 5 897 caen donde el viejo no tenía
candidato, 103 junto a o dentro de un candidato viejo de otro motivo y 61 son
la misma repetición que un candidato viejo con otras fronteras.

**Lo que se pierde, sin decidirlo aquí.** De los 159 candidatos que solo da el
viejo:

- **54** siguen ahí con otras fronteras (la misma repetición, por ejemplo un
  `(TA)n` que ahora acaba 6 bases antes porque su cola de otro motivo se
  cortó).
- **100** se pierden por solape con un candidato nuevo de otro motivo: 58
  solapan **1 base**, casi siempre una tira de una sola base que empieza en la
  última base de un microsatélite (`…TATATATATATA AAAAAAAA`: el `(TA)n`
  nuevo acaba en la primera `A` de la tira y la regla de solapes deja uno). Es
  la imagen especular, en el **final** del tramo, del caso que el responsable
  decidió el 29/10 para el **inicio**; esa decisión no cubre el final. Los
  otros 42 solapan de 2 a 5 bases o con dos trozos, casi siempre la copia
  parcial que un trozo comparte con el motivo siguiente (`(AC)n` contra un
  `(CACG)n` pegado; el compuesto `(TCTA)×8 (TCTG)×8 (TCTA)×8` de las pruebas
  pierde así su `(TCTG)n`, que solapa 3 bases con cada vecino).
- **5** desaparecen sin nada nuevo encima. Son repeticiones con una
  interrupción que desplaza la fase (`ATATATATA AAA TATATATATAT`): el viejo las
  aceptaba enteras, con pureza ≥ 0,8 por la descomposición en copias
  completas; el nuevo las corta en la interrupción y no las vuelve a unir
  (comprobado en una copia aparte sobre cuatro de los cinco: salen dos trozos
  y ningún tramo unido), porque la unión mide la pureza en una sola fase, la
  del trozo de fuera, y la interrupción la desplaza, como documenta el
  original; los dos trozos no alcanzan el suelo de copias.

El saldo con el suelo por omisión es **+5 902** candidatos.

**Decisión del responsable (30/10/2026), corregida el 31/10/2026.** El
30/10 se decidió esperar a una versión nueva de `pop_var_caller` para todas
las pérdidas de la frontera. El 31/10 una sonda demostró que el buscador
emite bien los dos tramos en `(CA)×10 A×12` y los tres trozos del compuesto,
y que es **nuestra** regla de solapes (DP-SSR-05) la que descarta uno entero:
esas pérdidas (los 100 por solape) son nuestras, y el arreglo propuesto
(recortar en vez de descartar cuando el solape es menor que una copia del
motivo) queda pendiente de decisión. Solo el microsatélite corto entre dos
tiras largas de una misma base y los 5 de desplazamiento de fase son del
buscador y esperan su nueva versión. Detalle en el registro de deudas
técnicas (DT-061). Mientras tanto, la entrega 3 sigue con el núcleo de
`fa215c23` y congela el catálogo que este da.

*Nota del 01/11/2026: el arreglo propuesto para la regla de solapes se decidió
y se aplicó; ver el apartado 11.12.*

### 11.12. La regla de solapes recorta en vez de descartar (01/11/2026)

**Qué cambia.** DP-SSR-05 se reescribe (apartado 3): un aspirante que pisa a los
aceptados **una copia de su propio motivo o menos por cada lado** (solape ≤
periodo) se recorta a la parte que no solapa y vuelve a pasar por todos los
criterios; con más de una copia por algún lado, o con un aceptado dentro o
alrededor, se descarta entero como antes. Decisión del responsable del
01/11/2026 (opción A: con «menos de una copia» una tira de una sola base no
se habría recortado nunca, porque un solape de una base ya es una copia). El
núcleo copiado del buscador no cambia; solo cambian `resolve_overlaps` y el
recuento de `detect_repeat_candidates_in_regions`. Es el cierre de la parte
(a) de DT-061.

**Medición.** Mismos FASTA (SHA-256 `f5a8d13a…`, 13 secuencias, 78 252 003
bases) y mismos comandos que el apartado 11.1, `rustc 1.97.1`. «Antes» es el
árbol de la entrega 3 (`src/repeat_detection.rs` de SHA-256 `de1c504a…`,
binario de medida `ecab10be…`), que reproduce **cifra a cifra** el apartado
11.11; «después», el árbol con la regla nueva (`254ca582…`, binario
`86fdb761…`). Los dos binarios se compilaron en copias aparte, sin tocar el
árbol de trabajo. **Dos pasadas de `measure_detection_stages` por árbol, con
recuentos idénticos línea a línea** (se compararon sin los tiempos).
Suelo por omisión `[8,6,6,6,5,4]`, antes → después; solo cambian las dos
últimas filas:

| Etapa | p1 | p2 | p3 | p4 | p5 | p6 | Total |
|---|---:|---:|---:|---:|---:|---:|---:|
| Intervalos del detector | 15 822 760 | 2 789 346 | 1 178 087 | 341 924 | 132 757 | 57 244 | 20 322 118 |
| Suelo de copias | 15 773 834 | 2 783 034 | 1 175 382 | 341 711 | 132 727 | 57 022 | 20 263 710 |
| Pureza | 27 | 117 | 826 | 49 | 5 | 86 | 1 110 |
| Tirada | 0 | 0 | 0 | 0 | 0 | 3 | 3 |
| Tope | 0 | 0 | 10 | 0 | 0 | 0 | 10 |
| Borde | 1 | 0 | 3 | 0 | 0 | 0 | 4 |
| Fuera de regiones | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| Solape | 306 → **269** | 552 → **209** | 200 → **59** | 28 → **7** | 0 → 0 | 9 → **1** | 1 095 → **545** |
| **Aceptados** | 48 592 → **48 629** | 5 643 → **5 986** | 1 666 → **1 807** | 136 → **157** | 25 → 25 | 124 → **132** | 56 186 → **56 736** |

Los intervalos del detector y las etapas de cada tramo por sí solo no cambian
(la regla actúa después). El saldo es **+550** candidatos: cada recortado que
entra sale de «solape», y los demás descartes por solape también bajan porque
un recortado deja sitio al siguiente. Con los otros suelos: el bajo pasa de
274 286 a **275 626** aceptados (solapes 3 341 → 2 001) y el alto, de 16 122 a
**16 402** (solapes 515 → 235). Bases en candidatos con el suelo por omisión:
673 012 → **684 434** (0,860 % → 0,875 %). Pares de candidatos que se tocan:
1 012 → **1 616**. Entre los candidatos nuevos no queda ningún par solapado
(comprobado sobre el volcado de todo el genoma).

**Cuántos de los 100 se recuperan.** Los 100 candidatos de la sección 11.11
(los que el núcleo viejo `deecb79` aceptaba y se perdían por solape con uno
nuevo de otro motivo) se identifican comparando los volcados de candidatos de
los tres árboles: un candidato viejo «sigue ahí» si hay uno actual de su
motivo, rotado o no, que lo solape; da 54 que siguen, **100 perdidos por
solape** y 5 desaparecidos, como en 11.11. Con la regla nueva:

- **36 se recuperan** (34 de periodo 1, 1 de periodo 2 y 1 de periodo 3);
- **24 caen aun así**: todas son tiras de una sola base de **exactamente 8
  bases** (el suelo) que pisan una base; al recortarlas quedan en 7 y **el
  recortado ya no llega al suelo**. Se cuenta como solape (decisión de este
  apartado) y no en el suelo;
- **40 siguen descartados** porque el solape es de más de una copia del que
  pierde: 12 tiras de una base que pisan 2 bases, una cubierta entera, 22
  dinucleótidos que pisan 3 bases, 4 que pisan 5 y un trinucleótido que pisa 5
  (`(AC)n` contra un `(CACG)n` pegado y casos así). Es lo que la regla
  aprobada deja fuera;
- los 5 desaparecidos (repeticiones con una interrupción que desplaza la
  fase) siguen sin aparecer: son del buscador, parte (b) de DT-061.

De las 58 pérdidas por **una sola base**, **34 se recuperan** y 24 caen por el
suelo tras recortar. En total hay 583 candidatos que antes no estaban y 33
que antes estaban y ya no figuran con esas coordenadas (22 de ellos están
ahora recortados por un extremo; los otros 11 no se han clasificado).

**Barrido** (`measure_scan`, cinco pasadas por árbol, alternadas, en la misma
sesión y sin nada más en marcha): `scan_seconds` antes 3,96–5,43 s (mediana
5,02), después 5,21–5,71 s (mediana 5,41). La dispersión entre pasadas del
**mismo** binario (1,5 s) es mayor que la diferencia entre árboles, así que
no se afirma ninguna diferencia de tiempo: el recorte solo corre en los aspirantes que pisan
a uno ya aceptado, una fracción mínima de los candidatos (los descartes por
solape de antes eran 1 095 de 56 000). Memoria residente
máxima: antes
181 072–181 768 kB, después 181 604–181 772 kB (cargar el FASTA solo, 156 996
kB): **no cambia**.

**Lo que no cambia.** Las huellas congeladas (la de escala y las demás) no
cambian: ninguna de las pruebas que las comparan falla. El catálogo detectado
del tomate sí cambia (56 186 → 56 736 loci), pero ninguna cifra congelada
depende de él.

**Comprobación de que ninguna prueba gana un locus que no declara.** Con un
registro de cada locus detectado que entra en el catálogo, la batería entera
(2 213 pruebas antes, 2 223 después) en el árbol de antes y en el de después da
**el mismo conjunto de loci detectados, prueba a prueba**, salvo los de las dos
pruebas nuevas; y todos los loci detectados de la batería salen de
`tests/repeat_detection_run.rs` (ningún otro fichero de pruebas ni ningún
subproceso registra uno). Las secciones de detección que fabrica
`tests/support/detection_section.rs` (seis bases de margen) no ganan nada con
el recorte.
