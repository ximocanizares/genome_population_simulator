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
`gps-hiseqx-pcrfree-pe150-v1`: aproximación propia a Illumina HiSeq X con
biblioteca WGS sin PCR, pares de 150 pb e inserto medio de 350 pb. Sus curvas
R1/R2 son un recurso propio versionado; ART `HSXn` se usa exclusivamente como
referencia externa de calibración y no se incorpora ni redistribuye. El inserto sigue una normal
truncada con desviación estándar de 20 pb y límites 250–450 pb. Longitud,
inserto y perfil son configurables. Cuando ocurre un error, el perfil inicial
elige uniformemente una de las otras tres bases; una matriz empírica queda en
DTF-053. Sus curvas admiten una longitud común a R1/R2 de 1–150 pb: una menor
usa su prefijo y una mayor se rechaza sin extrapolar; longitudes asimétricas
entre mates quedan fuera de F4. El inserto mínimo debe ser al menos esa longitud
común, pues F4 no modela adaptadores. Por defecto, los duplicados PCR no aplican,
los ópticos, quimeras, sesgo GC e *index hopping* tienen tasa cero; la
contaminación usa rango predeterminado `0–0 %`; los adaptadores no aparecen
porque los insertos mínimos superan la
longitud de lectura y la capacidad se aplaza en DTF-054. Los duplicados y el
sesgo GC no neutro quedan aplazados en DTF-056. F4 no genera UDI
inesperados; si alguno alcanza el demultiplexado se descarta defensivamente y
su generación explícita queda en DTF-025. Permanecen configurables únicamente
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
demultiplexado → FASTQ/resumen. Ninguna etapa que use una pareja, una instancia
transferida o una lectura quimérica se apoya en identificadores derivados del
orden físico o de un RNG compartido. El catálogo RNG F4 extenderá, sin modificar
las filas previas de DTC-037, ámbitos locales para: plantilla/fragmento y
quimera (`event_id`, `library_id`, `template_ordinal`); calidad/error
(`event_id`, `library_id`, `template_ordinal`, mate, ciclo); porcentaje,
selección de donante y retirada propia de contaminación (receptora y pareja
candidata); e *index hopping* (biblioteca actual y pareja técnica). DTC-079 permite una solución local proporcionada: entradas
externas se validan, Rust conserva la autoridad semántica y no se crean
protocolos, almacenes o defensas no necesarios.

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

## 5. F4.1 — Perfiles tecnológicos y pares

**Resultado:** el único evento selecciona un perfil Illumina efectivo; las
bibliotecas de la única población seleccionada son compatibles y una unidad PE
conserva dos *mates* con una identidad de pareja única.

- ☐ F4.1.1. Convertir el contrato fijo `illumina-se-minimal-v1` de F3 en una
  frontera de perfiles versionados, sin asignar valores cuantitativos no
  aprobados.
  - Depende de: PAUSA-DISEÑO-F4 y DTPD-106 aprobada.
  - Procede de: REQ-TEC-01/05/09–10, REQ-DAT-12/14/16/18, REQ-CONF-66/82/101–102,
    REQ-MAN-54/60/78 y arquitectura §5.5/§6.
  - Cambia: `src/application.rs`, `src/configuration.rs`,
    `src/technical_identity.rs`, `src/event_manifest.rs` y contratos TOML/JSON
    vigentes; registra ID, versión, preparación y parámetros/sustituciones con
    procedencia. Rechaza tecnologías de lectura larga, perfiles desconocidos y
    cualquier intento de declarar más de un evento o población en F4.
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
    `<mate>` vale `1`/`2` y, con UDI, `<index>` vale `i7+i5` (sin UDI, `0`). No
    incluye campos de instrumento, carril, *tile*, coordenadas ni BCL. Ámbitos
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
    fragmento y mantienen extremos y `N` literal/Phred 0. El objetivo numérico PE se
    expresa en parejas completas. Cada fragmento se sortea directamente entre
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
    la vez el sorteo de error y la calidad emitida. Rechaza probabilidades
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

> **PAUSA-F4.2.** Confirmar que los perfiles cuantitativos se validan contra sus
> datos y que las diferencias aleatorias se informan, no se convierten en
> tolerancias de aceptación.

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
    siendo una unidad completa. Cada quimera une exactamente dos segmentos no
    contiguos con punto de unión uniforme dentro del inserto. El perfil inicial
    fija frecuencia cero. Una tasa no nula se sortea Bernoulli e
    independientemente por pareja limpia antes de contaminación y
    demultiplexado, sin añadir ni retirar parejas ni publicar procedencia
    individual; informa el recuento realizado sin forzarlo al esperado.
  - Se comprueba con: `cargo test --test f4_chimeras_pairs --locked`; verdad
    manual comprueba uniones, orígenes transitorios, orientación, inserto y
    ausencia de procedencia individual; puertas comunes.

> **PAUSA-F4.3.** Revisar que cada artefacto tiene parámetro/perfil y que no
> altera la verdad biológica. Solo puede aumentar los contadores agregados de
> DTPD-108, sin retener todas las lecturas en memoria.

## 8. F4.4 — Contaminación biológica

- ☐ F4.4.1. Sortear uniformemente por receptora el porcentaje total dentro de un
  rango y repartir sus unidades al azar entre los otros individuos de su misma
  población. Por defecto recibe toda la población seleccionada, con subconjunto
  explícito opcional, sin crear otra biblioteca ni mezclar disposiciones.
  - Depende de: PAUSA-DISEÑO-F4 y F4.1.3.
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

- ☐ F4.5.1. Admitir índices simples, duales combinatorios y UDI, junto con
  tasas, procedencia y descarte defensivo de combinaciones inesperadas.
  - Depende de: PAUSA-DISEÑO-F4 y F4.4.2.
  - Procede de: REQ-CONF-78–81, REQ-DAT-29–32/70–71 y REQ-CONT-10–11.
  - Cambia: configuración de tanda, plan de demultiplexado y validación de
    colisiones; UDI permite reconocer y filtrar, sin asumir que elimina el
    *index hopping*. El perfil inicial usa UDI y fija *index hopping* a cero.
    Los UDI son configurables por biblioteca; sin ellos, la cabecera usa el
    marcador neutro literal `0`. F4 no genera combinaciones UDI inesperadas;
    si alguna alcanzara esta frontera se descarta defensivamente, y su modelo
    explícito queda en DTF-025.
    Una tasa no nula se sortea Bernoulli e independientemente por pareja de la
    tanda; cada seleccionada se transfiere uniformemente a otra biblioteca y
    exige al menos dos. Opera después de contaminación biológica, por lo que
    puede transferir una pareja ya contaminante, e informa el recuento realizado
    sin forzarlo al esperado. La transferencia no duplica contenido y la salida
    recibe el identificador técnico local de la biblioteca final.
  - Se comprueba con: `cargo test --test f4_demultiplexing --locked`; tabla de
    casos cubre traslado, descarte defensivo de UDI inesperado, conservación de
    tanda y colisión; puertas comunes.

- ☐ F4.5.2. Publicar una pareja final PE por individuo/evento y su resumen
  agregado en manifiesto.
  - Depende de: F4.5.1.
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

> **PAUSA-F4.5.** Revisar coherencia transversal FASTQ–resumen agregado–manifiesto y
> que los descartes reducen cobertura final sin regeneración.

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
| F4.3 | Rechazos de artefactos y quimeras | rechazos de perfil y oráculo de quimeras | aceptar modelo no calibrado o alterar verdad |
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
