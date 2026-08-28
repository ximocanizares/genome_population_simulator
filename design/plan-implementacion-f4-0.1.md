# Plan de implementación de F4 0.1

## 1. Cabecera y estado

- **Título:** realismo técnico WGS: perfiles, pares, degradaciones y modos de
  alineamiento.
- **Versión:** 0.1.
- **Fecha:** 28 de agosto de 2026.
- **Estado:** borrador de planificación; no autoriza implementación mientras
  permanezcan abiertas las decisiones bloqueantes enumeradas en §4.
- **Ámbito:** F4.1–F4.7 de la hoja de ruta vigente, construido sobre el recorrido
  mínimo F3.1/F3.5 y sus artefactos `.fastq.gz`.
- **Fuentes:** especificación §4.4, §4.6 y §4.9; arquitectura §3, §5.5, §6–§7,
  §8–§10 y §12; hoja de ruta F4; trazabilidad vigente; DTC-037–045,
  DTC-049, DTC-078 y DTC-079; DTPD-104–108.

F4 amplía el único recorrido técnico real de F3 sin crear otra vía biológica:
la simulación WGS lee instantáneas faseadas, y artefactos, contaminación,
errores y alineamientos nunca modifican la verdad genotípica, `truth.jsonl` ni
el VCF de verdad. La base reutilizable actual es `src/wgs_sampling.rs`,
`src/fastq_output.rs`, `src/fastq_validation.rs`,
`src/technical_identity.rs` y `src/event_manifest.rs`. El FASTQ F3 se escribe
por flujo BGZF/gzip, se valida descomprimido y publica tamaño y SHA-256 del
archivo comprimido; F4 conserva esas propiedades.

En este primer alcance de F4, una ejecución técnica consume individuos de una
sola población y declara un único `[technical_event]`. No se introduce
`[[technical_events]]`, ni se materializan varias poblaciones o eventos
técnicos en una misma ejecución. La extensión multievento y multipoblación se
recuperará solo desde DTF-050.

## 2. Alcance, exclusiones y orden

**Incluye:** perfil versionado incorporado de Illumina PE;
longitud, cobertura, error y calidad; pares completos; quimeras, contaminación
biológica y demultiplexado;
resumen técnico agregado sin secuencias duplicadas y FASTQ finales comprimidos.

El perfil predeterminado de F4 es
`gps-hiseqx-pcrfree-pe150-v1`, definido literalmente por
`resources/f4/gps-hiseqx-pcrfree-pe150-v1.json` (SHA-256
`ae8d598a95f7a9826b3a8c563a023d745f85031651c1e8c03ee9f54310162d3b`). Es una
aproximación propia y sintética, no una calibración empírica: R1 recorre
Q36→Q32→Q28 y R2 Q35→Q30→Q25 en los tres tramos declarados del recurso. ART
`HSXn` es una referencia comparativa externa y no se incorpora ni redistribuye.
El inserto usa `normal_truncated_box_muller_v1` con media 350 pb, desviación
20 pb y límites inclusivos 250–450 pb. Longitud, inserto y perfil son
configurables. Cuando ocurre un error, el perfil inicial
elige uniformemente una de las otras tres bases; una matriz empírica queda en
DTF-053. Sus curvas admiten una longitud común a R1/R2 de 1–150 pb: una menor
usa su prefijo y una mayor se rechaza sin extrapolar; longitudes asimétricas
entre mates quedan fuera de F4. El inserto mínimo debe ser al menos esa longitud
común, pues F4 no modela adaptadores. Por defecto, los duplicados PCR no aplican,
los ópticos, quimeras, sesgo GC e *index hopping* tienen tasa cero; la
contaminación usa rango predeterminado `0–0 %`; los adaptadores no aparecen
porque los insertos mínimos superan la
longitud de lectura y la capacidad se aplaza en DTF-054. Los duplicados y el
sesgo GC no neutro quedan aplazados en DTF-056. F4 exige UDI explícitos de
ocho bases `ACGT` para `i7` e `i5`, únicos por componente dentro de la tanda, y
no genera UDI inesperados; si alguno alcanza el demultiplexado se descarta
defensivamente y su generación explícita queda en DTF-025. Permanecen configurables únicamente
los mecanismos definidos en este plan: contaminación, quimeras e *index
hopping*; DTPD-106 conserva el perfil inicial y sus ajustes aprobados.

**Excluye:** estimar, calibrar o reutilizar perfiles empíricos; tecnologías de
lectura larga, incluido PacBio HiFi, Oxford Nanopore y PacBio CLR; protocolos
dirigidos, carriles, varias bibliotecas/tandas dentro de un evento, un segundo
FASTQ limpio, BAM/CRAM, mapeo real y el modo directo de alineamientos sintéticos
(DTF-051/055), perfiles o
degradaciones no especificados, automatización de callers, escenarios/sesiones
paralelos y endurecimiento contra atacante local activo. F5–F9 conservan sus
destinos de referencias complejas, poblaciones, paralelismo, presupuestos y
distribución.

El orden normativo de F4 es cobertura → fragmento → quimera opcional →
calidad/error por mate → contaminación biológica → *index hopping* y
demultiplexado → FASTQ/resumen. DTC-081 fija `f4-technical-v1`: cada decisión
reinicia su ChaCha12 desde una clave estable, con `event_id`, `library_id`,
`stage` y plantilla o instancia pre-salto antes de los componentes específicos.
`bernoulli-v1`, `uniform-index-v1` y Fisher–Yates parcial sobre candidatos
canónicamente ordenados son los únicos sorteos discretos. Probabilidades 0/1
no sortean y los espacios vacíos o incapaces se rechazan antes de simular. El
*index hopping* usa una instantánea post-contaminación, evalúa cada instancia
una sola vez y aplica transferencias simultáneas; no hay cadenas ni dependencia
del orden de bibliotecas. La distribución cuantitativa del perfil permanece en
DTPD-106. DTC-079 permite una solución local proporcionada: entradas externas
se validan, Rust conserva la autoridad semántica y no se crean protocolos,
almacenes o defensas no necesarios.

## 3. Requisitos previos reutilizables

- F3 está cerrado y aporta el único evento explícito, biblioteca y grupo de
  lectura únicos, FASTQ final por flujo, `event_manifest`, publicación
  no-clobber, hashes y validación independiente de origen, cobertura y `N`.
- F2 debe continuar aportando referencia tipada, estado haplotípico faseado,
  materialización temporal, identidades, semilla raíz y subsemillas aisladas.
- Todo campo F4 admitido por la configuración general tendrá, dentro del mismo
  hito que lo introduce, su recorrido local Marimo → TOML → validación Rust y
  una prueba de paridad. Marimo solo autoriza y muestra la configuración: no
  incorpora lógica biológica, sesión ni servicio.
- Todos los cambios usan Rust 1.97.1, Cargo bloqueado y las puertas existentes:

      cargo fmt --all -- --check
      cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
      cargo test --workspace --all-targets --all-features --locked
      cargo test --workspace --doc --all-features --locked
      cargo doc --workspace --no-deps --all-features --locked
      git diff --check

- Las pruebas focales se añadirán en `tests/f4_*.rs`; un productor no se valida
  solo reproduciendo su propia lógica. Las pruebas sintéticas pequeñas usarán
  oráculos independientes y los resúmenes agregados no retendrán procedencia
  oráculo independiente.

## 4. Decisiones bloqueantes y pausa de diseño

Este documento no resuelve las decisiones siguientes. Son necesarias antes del
primer paso ejecutable que dependa de ellas.

| Bloqueo | Bloquea | Lo que debe decidir el documento propietario |
|---|---|---|
| DTPD-106 | — | Resuelta: perfil inicial propio `gps-hiseqx-pcrfree-pe150-v1`, calibrado externamente frente a ART HSXn sin incorporar ni redistribuir sus tablas, y ajustes de longitud, inserto, errores, calidad, quimeras e índices aprobados. Duplicados y GC no neutro quedan en DTF-056; adaptadores, en DTF-054; anomalías PE, en DTF-055. |
| DTPD-108 | — | Resuelta: manifiesto/resumen agregado sin fichero ni origen individual de lectura; conserva semilla, perfil, coberturas, contaminación/aportes, artefactos, descartes y hashes/tamaños. |
| DTPD-104 | DTF-055 | Adaptador SAM/BAM/CRAM, versión, cabeceras/grupos de lectura, orden e índices. |
| DTPD-105 | DTF-055 | Descubrimiento e invocación reproducible de minimap2/alternativas, compatibilidad de versión, *presets*, parámetros, diagnósticos y tratamiento de fallos. |

> **PAUSA-DISEÑO-F4.** Cerrada: DTPD-106 y DTPD-108 están resueltas. BAM/CRAM y
> mapeo no pertenecen a F4 y quedan en DTF-055. La implementación no podrá
> reabrir tasas, distribuciones, formatos ni dependencias desde el código.

F4.7 es independiente de esos bloqueos: su rechazo explícito puede ejecutarse
cuando F3 esté cerrado, sin anticipar un formato de perfil empírico.

## 4. F4.0 — Catálogo RNG técnico

- ☐ F4.0.1. Publicar `f4-technical-v1` exactamente como DTC-081 y registrar su
  versión en la configuración efectiva y el manifiesto de evento.
  - Depende de: PAUSA-DISEÑO-F4.
  - Procede de: DTC-037/038/081, REQ-DAT-15, REQ-CONF-113 y REQ-REP-22/23.
  - Cambia: proveedor de aleatoriedad y contratos de reproducibilidad; no
    introduce perfiles, salidas ni estados persistentes por lectura.
  - Se comprueba con: `cargo test --test f4_rng_catalog --locked`; vectores
    literales de subsemillas/resultados, ramas 0/1, rechazo de cardinalidad
    imposible, permutación de bibliotecas y ejecución paralela.

- ☐ F4.0.2. Probar la regeneración de una pareja desde su clave y que cobertura,
  quimera, error, contaminación y *index hopping* no consumen ni alteran el
  ámbito de otra decisión.
  - Depende de: F4.0.1.
  - Procede de: DTC-081 y REQ-DAT-15/70–71.
  - Cambia: sólo fixtures y comparadores de prueba; conserva el catálogo como
    contrato único.
  - Se comprueba con: `cargo test --test f4_rng_isolation --locked`; un vector
    manual compara regeneración, reordenación, hilos y una transferencia
    simultánea sin cadenas.

## 5. F4.1 — Perfiles tecnológicos y pares

**Resultado:** el único evento selecciona un perfil Illumina efectivo; las
bibliotecas de la única población seleccionada son compatibles y una unidad PE
conserva dos *mates* con una identidad de pareja única.

- ☐ F4.1.1. Convertir el contrato fijo `illumina-se-minimal-v1` de F3 en la
  frontera del único perfil F4 Illumina PE versionado, sin asignar valores
  cuantitativos no aprobados.
  - Depende de: PAUSA-DISEÑO-F4, F4.0.2 y DTPD-106 aprobada.
  - Procede de: REQ-TEC-01/05/09–10, REQ-DAT-12/14/16/18, REQ-CONF-66/82/101–102,
    REQ-MAN-54/60/78 y arquitectura §5.5/§6.
  - Cambia: `src/application.rs`, `src/configuration.rs`,
    `src/technical_identity.rs`, `src/event_manifest.rs` y contratos TOML/JSON
    vigentes; actualiza sin compatibilidad el TOML general a v2, la configuración
    efectiva y el manifiesto general a v4 y el manifiesto de evento a v2; verifica
    el recurso y SHA-256 fijados por DTPD-106 y registra ID,
    versión, recurso, digest, preparación y parámetros/sustituciones con
    procedencia. Rechaza Illumina SE, tecnologías de lectura larga, perfiles
    desconocidos y cualquier intento de declarar más de un evento o población
    en F4.
  - Se comprueba con: `cargo test --test f4_profiles --locked`; oráculo de
    admisión comprueba perfiles, versiones, exclusividad y hechos del manifiesto;
    puertas comunes.

- ☐ F4.1.2. Generalizar las identidades técnicas para plantilla, pareja y mate
  sin renombrar ni colisionar con las plantillas SE F3.
  - Depende de: F4.1.1.
  - Procede de: REQ-DAT-18/68/71, REQ-RG-03–04 y DTC-037–039.
  - Cambia: `src/identity.rs`, `src/technical_identity.rs` y pruebas; la pareja
    tiene nombre único dentro del evento y la cabecera FASTQ es exactamente
    `@GPS:<event_id>:<library_id>:<technical_pair_id> <mate>:N:0:<index>`;
    `<mate>` vale `1`/`2` y `<index>` vale siempre `i7+i5`: ambos componentes
    son `ACGT` de longitud ocho, se escriben en orientación de secuenciación y
    no se complementa `i5`; `0` se rechaza en F4. No incluye campos de
    instrumento, carril, *tile*, coordenadas ni BCL. Ámbitos
    RNG de plantilla, fragmento, mate y artefacto se documentan y no consumen el
    de cobertura.
  - Se comprueba con: `cargo test --test f4_technical_identity --locked`; vector
    literal comprueba repetición, aislamiento al reordenar bibliotecas,
    unicidad y ausencia de pares huérfanos; puertas comunes.

- ☐ F4.1.3. Implementar el muestreador de fragmentos y la emisión PE sólo con
  orientación `FR` y distribución de inserto ya aprobadas.
  - Depende de: F4.1.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-12/15/18/65–66, REQ-WGS-01–11 y REQ-VAL-07.
  - Cambia: `src/wgs_sampling.rs`, nuevo módulo de pares y
    `src/fastq_output.rs`; R1 es directo desde el extremo izquierdo y R2 el
    complemento inverso desde el derecho; ambos mates proceden del mismo
    fragmento y mantienen extremos y `N` literal/Phred 0. F4 no admite un
    objetivo numérico: F4.2 deriva el número entero de parejas desde cobertura.
    Cada fragmento se sortea directamente entre
    inicios válidos que permiten alojar el inserto completo, sin recortar ni
    reintentar extremos inválidos. La longitud configurable es única y común a
    R1/R2, entre 1 y 150 pb, y no puede superar el inserto mínimo. Un
    solapamiento entre mates se conserva como
    dos lecturas completas, sin fusión.
  - Se comprueba con: `cargo test --test f4_paired_end --locked`; verdad manual
    de fragmentos cortos verifica orientación, inserto, extremos, calidad e
    imposibilidad de un único mate; puertas comunes.

> **PAUSA-F4.1.** Revisar perfiles, semántica de pareja, IDs y manifiestos antes
> de aplicar errores o artefactos. No continuar si el perfil aún codifica una
> distribución o un *preset* no aprobado.

## 6. F4.2 — Cobertura, longitud, error y calidad

**Resultado:** cada biblioteca aplica el perfil y sus sustituciones efectivas
con muestreo WGS básico separado de la evidencia final.

- ☐ F4.2.1. Extender el plan WGS para longitudes variables y objetivos Illumina
  PE sin alterar el denominador haploide acordado.
  - Depende de: F4.1.3 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-14–15/24, REQ-WGS-01–11/14 y REQ-CONF-72–74/83–84.
  - Cambia: `src/wgs_sampling.rs`, configuración efectiva e informe de
    validación; conserva selección uniforme básica, reparto por copia/ploidía,
    `N` elegible y cobertura básica previa a artefactos. Una pareja contribuye
    con las bases emitidas por ambos mates, también cuando se solapan. Un objetivo
    no entero se redondea a la pareja más próxima, con empate hacia arriba. Las
    secuencias menores que el inserto mínimo se excluyen y se informan a cero.
    F4 admite únicamente una cobertura objetivo haploide común por evento;
    rechaza objetivos por número de parejas y sustituciones por individuo,
    secuencia o región.
  - Se comprueba con: `cargo test --test f4_wgs_profiles --locked`; oráculo
    independiente comprueba denominador, objetivos alternativos, pares enteros,
    cobertura por copia y `N`; puertas comunes.

- ☐ F4.2.2. Aplicar el modelo de error y calidad aprobado como transformación
  transitoria de una lectura verdadera, sin resolver ni crear `N`.
  - Depende de: F4.2.1 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-37–39/72–73, REQ-MAN-42 y REQ-VAL-01/06/09–13.
  - Cambia: nuevo módulo de perfil/error y `src/fastq_output.rs`; acumula solo
    los contadores necesarios para el resumen de DTPD-108, sin conservar por
    instancia diferencias verdad/observación. Tras un error, sortea
    uniformemente una de las otras tres bases. R1 y R2 consumen sorteos
    independientes de calidad/error aun cuando se solapen. Cuando haya
    desplazamiento de calidad y escala de error, aplica el primero, deriva la
    probabilidad y aplica después la escala; la probabilidad final gobierna a
    la vez el sorteo de error y la calidad emitida. La curva por ciclo se
    interpola y redondea exactamente como el recurso DTPD-106; el inserto usa
    su `normal_truncated_box_muller_v1`, sin caché de variates. Rechaza probabilidades
    finales inválidas, sin saturarlas silenciosamente. Estos ajustes son
    globales por evento y pueden diferir entre R1 y R2; no admite selectores
    por individuo, secuencia o región. Si se omiten, materializa y publica sus
    valores neutros: desplazamiento `0 Phred` y escala `1` para cada mate.
    Codifica las calidades como FASTQ Sanger (`Phred+33`), redondea al entero
    más próximo y rechaza un resultado fuera de `0..=93`, sin saturarlo.
  - Se comprueba con: `cargo test --test f4_error_quality --locked`; vector
    externo al productor comprueba sustituciones, posición/contexto si procede,
    relación calidad-error, `N` y flujos RNG aislados; puertas comunes.

- ☐ F4.2.3. Publicar e informar las tres coberturas de F4, sin convertir la
  cobertura FASTQ final en un objetivo nuevo.
  - Depende de: F4.2.1.
  - Procede de: REQ-DAT-24, REQ-WGS-11, REQ-MAN-37/61 y REQ-VAL-04.
  - Cambia: `src/fastq_validation.rs`, `src/event_manifest.rs` y resultados de
    validación; separa objetivo, WGS básica y FASTQ final propia/contaminante.
  - Se comprueba con: `cargo test --test f4_coverage_reporting --locked`; caso
    controlado demuestra que las tres medidas no se sustituyen entre sí;
    puertas comunes.

> **PAUSA-F4.2.** Confirmar mediante vectores que se conserva literalmente el
> recurso propio y su digest, y que las diferencias aleatorias se informan, no
> se convierten en tolerancias de aceptación. No se declara calibración empírica.

## 7. F4.3 — Artefactos WGS

- ☐ F4.3.1. Fijar y validar la preparación PCR-free del perfil inicial, sin
  generar duplicados ni admitir su configuración.
  - Depende de: F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-TEC-06–08, REQ-CONF-88–90 y REQ-DAT-46–49.
  - Cambia: configuración efectiva y validación; declara PCR no aplicable y
    ópticos a cero, y rechaza sus parámetros. La preparación amplificada y sus
    duplicados quedan en DTF-056.
  - Se comprueba con: `cargo test --test f4_profile_rejections --locked`;
    fixtures rechazan perfiles/tasas de duplicación y confirman los valores
    efectivos neutros; puertas comunes.

- ☐ F4.3.2. Declarar GC neutro y rechazar su configuración no nula, junto con
  adaptadores para el perfil inicial.
  - Depende de: F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-50–53, REQ-CONF-91–95, REQ-WGS-08/14 y REQ-MAN-47–48.
  - Cambia: configuración efectiva y validación; el perfil inicial usa GC
    neutro, tiene insertos mayores que las lecturas y rechaza adaptadores. El
    perfil GC calibrado queda en DTF-056 y los fragmentos cortos, en DTF-054.
  - Se comprueba con: `cargo test --test f4_profile_rejections --locked`;
    fixtures rechazan curva/intensidad GC y adaptadores, y confirman valores
    efectivos neutros; puertas comunes.

- ☐ F4.3.3. Crear quimeras desde fragmentos verdaderos, sin usar la
  contaminación como sustituto semántico.
  - Depende de: F4.1.3, F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-60–64, REQ-CONF-99, REQ-MAN-51 y REQ-VAL-07.
  - Cambia: módulo de fragmentos/artefactos y datos transitorios de origen; los
    fragmentos quiméricos son del mismo individuo/biblioteca y la pareja sigue
    siendo una unidad completa. Para inserto `L`, el primer segmento es el
    prefijo de la plantilla hasta un punto `b` uniforme en `1..L-1`; el segundo
    aporta `L-b` bases y se elige uniformemente de la lista canónica de todos
    los intervalos válidos del mismo individuo/biblioteca, ordenada por copia,
    secuencia e inicio. En la misma copia se excluyen intervalos que solapen o
    toquen el primero; ambos se leen en orientación directa y la PE normal se
    aplica al inserto compuesto. Si una pareja Bernoulli seleccionada no tiene
    segunda fuente válida, se falla antes de publicar. El perfil inicial fija
    frecuencia cero. Una tasa no nula se sortea Bernoulli e
    independientemente por pareja limpia antes de contaminación y
    demultiplexado, sin añadir ni retirar parejas ni publicar procedencia
    individual; informa el recuento realizado sin forzarlo al esperado.
  - Se comprueba con: `cargo test --test f4_chimeras_pairs --locked`; verdad
    manual comprueba punto, lista canónica, exclusiones por solape/contacto,
    orientación, inserto, fallo sin candidata y ausencia de procedencia
    individual; puertas comunes.

> **PAUSA-F4.3.** Revisar que cada artefacto tiene parámetro/perfil y que no
> altera la verdad biológica. Solo puede aumentar los contadores agregados de
> DTPD-108, sin retener todas las lecturas en memoria.

## 8. F4.4 — Contaminación biológica

- ☐ F4.4.1. Sortear uniformemente por receptora el porcentaje total dentro de un
  rango y repartir sus unidades al azar entre los otros individuos de su misma
  población. Por defecto recibe toda la población seleccionada, con subconjunto
  explícito opcional, sin crear otra biblioteca ni mezclar disposiciones.
  - Depende de: PAUSA-DISEÑO-F4, F4.0.2, F4.1.3, F4.2.2, F4.3.3 y
    PAUSA-F4.3.
  - Procede de: REQ-CONF-75–77, REQ-CONT-01–09/15–16 y REQ-DAT-25–28/69/71.
  - Cambia: configuración, planificador de asignación y manifiesto;
    resuelve receptoras (toda la población por omisión), sortea el objetivo de
    cada una, lo convierte desde su porcentaje y parejas limpias al entero más
    próximo (empate hacia arriba), preflighta donantes y capacidad conjunta; si
    falta, rechaza antes de escribir FASTQ. Registra subsemillas y aportes
    realizados. Selecciona sin reemplazo para cada receptora sobre el conjunto
    de parejas limpias de todos sus donantes elegibles, siempre desde sus
    lecturas propias previas a contaminación; no admite matriz dirigida ni FASTQ
    externo en F4. El resultado es un plan temporal de asignaciones por claves,
    sin secuencias, calidades ni FASTQ limpio.
  - Se comprueba con: `cargo test --test f4_biological_contamination --locked`;
    oráculo de asignación entero verifica rango por receptora, 0/100 %, reparto
    aleatorio reproducible, aislamiento poblacional, PE completo y ausencia de
    reemplazo de una unidad de origen dentro de la misma receptora, sin cadenas
    ni dependencia del orden de receptoras;
    puertas comunes.

- ☐ F4.4.2. Materializar sustitución biológica y sus recuentos agregados sin
  duplicar secuencias/calidades ni conservar vínculos por lectura.
  - Depende de: F4.4.1.
  - Procede de: REQ-CONT-01–05/13–16, REQ-DAT-25–27/34–36/69 y REQ-NOF-08.
  - Cambia: módulos de instancias, FASTQ final y manifiesto; la receptora recibe
    una copia técnica nueva y la donante no pierde la unidad. Retira al azar, sin
    reemplazo, una pareja propia por cada incorporada. Solo se publican recuentos
    agregados por donante y receptora. La copia conserva sus bases y calidades
    realizadas; no vuelve a consumir el RNG de error/calidad. El emisor la
    regenera bajo demanda desde su clave y la descarta tras escribirla.
  - Se comprueba con: `cargo test --test f4_contamination_provenance --locked`;
    lector independiente comprueba recuentos, ausencia de pares huérfanos y que
    no se publica procedencia individual; puertas comunes.

> **PAUSA-F4.4.** Confirmar que contaminación biológica copia a la receptora y
> demultiplexado técnico transferirá, no copiará, la misma instancia.

## 9. F4.5 — Demultiplexado técnico y FASTQ final

- ☐ F4.5.1. Validar UDI obligatorios y simular *index hopping* únicamente entre
  bibliotecas UDI de la tanda, con descarte defensivo de combinaciones
  inesperadas.
  - Depende de: PAUSA-DISEÑO-F4, F4.0.2, F4.4.2 y PAUSA-F4.4.
  - Procede de: REQ-CONF-78–81, REQ-DAT-29–32/70–71 y REQ-CONT-10–11.
  - Cambia: configuración de tanda, plan de demultiplexado y validación de
    colisiones; cada biblioteca declara `i7` e `i5` de exactamente ocho bases
    `ACGT`, ambos únicos dentro de la tanda y escritos como `i7+i5` sin invertir
    ni complementar `i5`. F4 rechaza ausencia de UDI, marcador `0`, índices
    simples o combinatorios, longitudes/alfabetos inválidos y duplicados. UDI
    permite reconocer y filtrar, sin asumir que elimina el *index hopping*. El
    perfil inicial fija *index hopping* a cero. F4 no genera combinaciones UDI
    inesperadas; si alguna alcanzara esta frontera se descarta defensivamente, y
    su modelo explícito queda en DTF-025.
    Una tasa no nula se sortea Bernoulli e independientemente por pareja de la
    tanda; cada seleccionada se transfiere uniformemente a otra biblioteca y
    exige al menos dos. Opera después de contaminación biológica, por lo que
    puede transferir una pareja ya contaminante, e informa el recuento realizado
    sin forzarlo al esperado. La transferencia no duplica contenido y la salida
    recibe el identificador técnico local de la biblioteca final.
  - Se comprueba con: `cargo test --test f4_demultiplexing --locked`; tabla de
    casos cubre gramática 8+8, orientación literal de `i5`, unicidad de ambos
    componentes, rechazo de `0`, traslado al UDI declarado de destino, descarte
    defensivo de UDI inesperado, conservación de tanda y colisión; puertas
    comunes.

- ☐ F4.5.2. Publicar una pareja final PE por individuo/evento y su resumen
  agregado en manifiesto.
  - Depende de: F4.5.1, F4.2.3 y PAUSA-F4.4.
  - Procede de: REQ-DAT-17/33–36/40, REQ-RG-02/04/06, REQ-SAL-08/09,
    REQ-MAN-40–41 y REQ-NOF-08.
  - Cambia: `src/fastq_output.rs`, `src/fastq_validation.rs`,
    `src/event_manifest.rs` y ensamblado. Mantiene gzip/BGZF, validación
    descomprimida, hash del comprimido y publica
    `<individual_id>_R1.fastq.gz` y `<individual_id>_R2.fastq.gz`; no guarda
    versión limpia ni procedencia individual.
  - Se comprueba con: `cargo test --test f4_final_fastq_provenance --locked`;
    validador independiente comprueba FASTQ, pareja, hashes y resúmenes, y que
    no se publica procedencia individual; puertas comunes.

- ☐ F4.5.3. Validar por flujo cada pareja final contra la reconstrucción
  independiente de su clave técnica transitoria antes de destruirla.
  - Depende de: F4.5.2.
  - Procede de: REQ-VAL-01–03/09–11, REQ-DAT-25–36/69–71 y DTPD-108.
  - Cambia: `src/fastq_validation.rs` y el orquestador provisional; el
    reconstructor consulta verdad biológica, plan compacto, perfil y catálogo
    RNG, pero no al escritor FASTQ. Compara literalmente R1/R2, calidad,
    cabecera, mate, donante o retirada biológica y biblioteca tras *hopping*.
    En éxito solo conserva el resultado agregado; en discrepancia detiene antes
    de publicar y conserva el diagnóstico sin un fichero de procedencia.
  - Se comprueba con: `cargo test --test f4_technical_replay_validation --locked`;
    fixture pequeño con haplotipos, calidades y bibliotecas distinguibles, más
    mutaciones negativas independientes de base, calidad, donante, retirada y
    salto, demuestra que ninguna supera la validación.

> **PAUSA-F4.5.** Revisar coherencia transversal FASTQ–resumen agregado–manifiesto,
> la reconstrucción independiente de cada pareja final y que los descartes
> reducen cobertura final sin regeneración.

## 10. F4.6 — Aplazado

Los BAM/CRAM y el mapeo real quedan fuera de F4 en DTF-055. La revisión de F4
termina tras F4.5: comprueba que toda degradación afecta únicamente al FASTQ y
que el resumen agregado no contiene procedencia individual.

## 11. F4.7 — Calibración retirada de v1

- ☐ F4.7.1. Rechazar antes de simular solicitudes de calibración, FASTQ o
  referencia destinados a ella y perfiles empíricos reutilizables.
  - Depende de: cierre F3; no depende de DTPD-106/108.
  - Procede de: REQ-CONF-85–87, REQ-DAT-41–45, REQ-VAL-49, REQ-MAN-43–44,
    REQ-SAL-10 y hoja de ruta F4.7.
  - Cambia: admisión TOML, diagnóstico y contratos de configuración; no crea
    artefactos, perfil ni entrada de manifiesto para una capacidad retirada.
  - Se comprueba con: `cargo test --test f4_calibration_rejection --locked`;
    fixtures cubren cada campo retirado, error previo y ausencia de salida;
    puertas comunes.

> **PAUSA-F4.7.** Confirmar que el rechazo no bloquea perfiles incorporados de
> F4 ni reintroduce estimación, contrato de entrada o reutilización empírica.

## 12. Matriz de evidencia y riesgos

| Hito | Resultado construido | Evidencia | Riesgo protegido |
|---|---|---|---|
| F4.1 | Perfiles, identidades y pares | vectores de perfil/ID, verdad manual de fragmentos | pares huérfanos, perfil o tecnología implícitos |
| F4.2 | Muestreo y errores cuantitativos | oráculo de cobertura, longitud, calidad y `N` | confundir cobertura básica/final o resolver `N` |
| F4.3 | Rechazos de artefactos y quimeras | rechazos de perfil y oráculo de quimeras | aceptar un modelo no declarado o alterar verdad |
| F4.4 | Contaminación biológica | oráculo de asignación conjunta y agregados | copia/retirada equivocada, exceso o par roto |
| F4.5 | Demultiplexado y FASTQ definitivo | lector independiente FASTQ+agregados+manifiesto | pérdida de instancia, secuencias duplicadas o salida limpia extra |
| F4.6 | Aplazado en DTF-055 | — | introducir BAM/CRAM o mapeo en F4 |
| F4.7 | Rechazo de calibración retirada | fixtures de admisión negativa | introducir perfil empírico o salida no admitida |

## 13. Recursos y trabajo posterior

F4 conservará producción/validación por flujo y no mantendrá una colección de
lecturas, parejas o secuencias sólo para trazabilidad. Antes de cada artefacto
grande se medirán memoria residente, tamaño FASTQ/procedencia y tiempo en un
escenario pequeño reproducible, pero los presupuestos, advertencias de recursos
y la evaluación a escala permanecen en F9.4 (REQ-NOF-01–07). La procedencia
obligatoria y los agregados de descartes cuentan desde F4 en cualquier medición;
no se reduce detalle para que un escenario continúe.

F5 debe reutilizar este único recorrido con verdad/referencias complejas; F6
debe demostrar que sus poblaciones alimentan el mismo recorrido; F9 cierra
paralelismo, endurecimiento, paquetes, formatos/validadores consolidados y
presupuestos. Ninguno justifica adelantar su implementación a F4.
