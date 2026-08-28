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

## 2. Alcance, exclusiones y orden

**Incluye:** perfiles versionados incorporados de Illumina SE/PE y PacBio HiFi;
longitud, cobertura, error y calidad; pares completos; duplicados, adaptadores,
sesgo GC, quimeras y anomalías PE; contaminación biológica; demultiplexado;
procedencia durable mínima sin secuencias duplicadas; modo directo; y, cuando
las decisiones de alineamiento estén cerradas, BAM/CRAM e integración con el
mapeador.

**Excluye:** estimar, calibrar o reutilizar perfiles empíricos; Oxford Nanopore,
PacBio CLR, protocolos dirigidos, carriles, varias bibliotecas/tandas dentro de
un evento, un segundo FASTQ limpio, un segundo BAM/CRAM limpio, perfiles o
degradaciones no especificados, automatización de callers, escenarios/sesiones
paralelos y endurecimiento contra atacante local activo. F5–F9 conservan sus
destinos de referencias complejas, poblaciones, paralelismo, presupuestos y
distribución.

El orden protege primero contratos y verdades de origen, después el muestreo y
los perfiles, luego degradaciones y procedencia, y finalmente las vistas de
alineamiento. Ninguna etapa que use una pareja, una instancia transferida o una
lectura quimérica se apoya en identificadores derivados del orden físico o de un
RNG compartido. DTC-079 permite una solución local proporcionada: entradas
externas se validan, Rust conserva la autoridad semántica y no se crean
protocolos, almacenes o defensas no necesarios.

## 3. Requisitos previos reutilizables

- F3 debe cerrar su pausa de revisión y aportar evento explícito, biblioteca y
  grupo de lectura únicos, FASTQ final por flujo, `event_manifest`, publicación
  no-clobber, hashes y validación independiente de origen, cobertura y `N`.
- F2 debe continuar aportando referencia tipada, estado haplotípico faseado,
  materialización temporal, identidades, semilla raíz y subsemillas aisladas.
- Todos los cambios usan Rust 1.97.1, Cargo bloqueado y las puertas existentes:

      cargo fmt --all -- --check
      cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
      cargo test --workspace --all-targets --all-features --locked
      cargo test --workspace --doc --all-features --locked
      cargo doc --workspace --no-deps --all-features --locked
      git diff --check

- Las pruebas focales se añadirán en `tests/f4_*.rs`; un productor no se valida
  solo reproduciendo su propia lógica. Los ficheros de procedencia se recorrerán
  por flujo y las pruebas sintéticas pequeñas conservarán verdad manual u
  oráculo independiente.

## 4. Decisiones bloqueantes y pausa de diseño

Este documento no resuelve las decisiones siguientes. Son necesarias antes del
primer paso ejecutable que dependa de ellas.

| Bloqueo | Bloquea | Lo que debe decidir el documento propietario |
|---|---|---|
| DTPD-106 | F4.1–F4.3 | Catálogo de perfiles incorporados, fuentes/datos de validación, distribuciones y valores cuantitativos de longitud, inserto, errores, calidad, duplicados, adaptadores, GC, quimeras y anomalías PE; sustituciones permitidas y versiones. |
| DTPD-108 | F4.4–F4.6 y artefactos finales de F4.3 | Esquema y formato de procedencia por lectura o pareja, enlaces de origen, diferencias verdad/observación, agregados de descartes y límites de volumen. |
| DTPD-104 | BAM/CRAM de F4.6 | Adaptador SAM/BAM/CRAM, versión, cabeceras/grupos de lectura, orden e índices. |
| DTPD-105 | mapeo real de F4.1/F4.6 | Descubrimiento e invocación reproducible de minimap2/alternativas, compatibilidad de versión, *presets*, parámetros, diagnósticos y tratamiento de fallos. |

> **PAUSA-DISEÑO-F4.** No iniciar F4.1–F4.6 hasta aprobar DTPD-106 y DTPD-108;
> no activar BAM/CRAM o mapeo hasta aprobar además DTPD-104 y DTPD-105. La
> resolución debe actualizar especificación/arquitectura/decisiones y
> trazabilidad antes de desmarcar los pasos afectados. Esta pausa no autoriza a
> elegir biología, tasas, distribuciones, formatos ni dependencias desde el
> código.

F4.7 es independiente de esos bloqueos: su rechazo explícito puede ejecutarse
cuando F3 esté cerrado, sin anticipar un formato de perfil empírico.

## 5. F4.1 — Perfiles tecnológicos y pares

**Resultado:** un evento selecciona una sola tecnología y perfil efectivo; las
bibliotecas de ese evento son compatibles y una unidad PE conserva dos *mates*
con una identidad de pareja única.

- ☐ F4.1.1. Convertir el contrato fijo `illumina-se-minimal-v1` de F3 en una
  frontera de perfiles versionados, sin asignar valores cuantitativos no
  aprobados.
  - Depende de: PAUSA-DISEÑO-F4 y DTPD-106 aprobada.
  - Procede de: REQ-TEC-01/05/09–10, REQ-DAT-12/14/16/18, REQ-CONF-66/82/101–102,
    REQ-MAN-54/60/78 y arquitectura §5.5/§6.
  - Cambia: `src/application.rs`, `src/configuration.rs`,
    `src/technical_identity.rs`, `src/event_manifest.rs` y contratos TOML/JSON
    vigentes; registra ID, versión, preparación y parámetros/sustituciones con
    procedencia. Rechaza tecnología mezclada, ONT, CLR y perfiles desconocidos.
  - Se comprueba con: `cargo test --test f4_profiles --locked`; oráculo de
    admisión comprueba perfiles, versiones, exclusividad y hechos del manifiesto;
    puertas comunes.

- ☐ F4.1.2. Generalizar las identidades técnicas para plantilla, pareja y mate
  sin renombrar ni colisionar con las plantillas SE F3.
  - Depende de: F4.1.1.
  - Procede de: REQ-DAT-18/68/71, REQ-RG-03–04 y DTC-037–039.
  - Cambia: `src/identity.rs`, `src/technical_identity.rs` y pruebas; la pareja
    tiene nombre único dentro del evento y cada miembro usa solo el indicador de
    mate admitido por FASTQ/SAM. Ámbitos RNG de plantilla, fragmento, mate y
    artefacto se documentan y no consumen el de cobertura.
  - Se comprueba con: `cargo test --test f4_technical_identity --locked`; vector
    literal comprueba repetición, aislamiento al reordenar bibliotecas,
    unicidad y ausencia de pares huérfanos; puertas comunes.

- ☐ F4.1.3. Implementar el muestreador de fragmentos y la emisión PE sólo con
  orientación y distribución de inserto ya aprobadas.
  - Depende de: F4.1.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-12/15/18/65–66, REQ-WGS-01–11 y REQ-VAL-07.
  - Cambia: `src/wgs_sampling.rs`, nuevo módulo de pares y
    `src/fastq_output.rs`; ambos mates proceden del mismo fragmento, mantienen
    orientación, extremos y `N` literal/Phred 0. El objetivo numérico PE se
    expresa en parejas completas.
  - Se comprueba con: `cargo test --test f4_paired_end --locked`; verdad manual
    de fragmentos cortos verifica orientación, inserto, extremos, calidad e
    imposibilidad de un único mate; puertas comunes.

> **PAUSA-F4.1.** Revisar perfiles, semántica de pareja, IDs y manifiestos antes
> de aplicar errores o artefactos. No continuar si el perfil aún codifica una
> distribución o un *preset* no aprobado.

## 6. F4.2 — Cobertura, longitud, error y calidad

**Resultado:** cada biblioteca aplica el perfil y sus sustituciones efectivas
con muestreo WGS básico separado de la evidencia final.

- ☐ F4.2.1. Extender el plan WGS para longitudes variables y objetivos SE, PE y
  HiFi sin alterar el denominador haploide acordado.
  - Depende de: F4.1.3 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-13–15/24, REQ-WGS-01–11/14 y REQ-CONF-72–74/83–84.
  - Cambia: `src/wgs_sampling.rs`, configuración efectiva e informe de
    validación; conserva selección uniforme básica, reparto por copia/ploidía,
    `N` elegible y cobertura básica previa a artefactos.
  - Se comprueba con: `cargo test --test f4_wgs_profiles --locked`; oráculo
    independiente comprueba denominador, objetivos alternativos, pares enteros,
    cobertura por copia y `N`; puertas comunes.

- ☐ F4.2.2. Aplicar el modelo de error y calidad aprobado como transformación
  transitoria de una lectura verdadera, sin resolver ni crear `N`.
  - Depende de: F4.2.1 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-37–40/72–73, REQ-MAN-42 y REQ-VAL-01/06/09–13.
  - Cambia: nuevo módulo de perfil/error y `src/fastq_output.rs`; conserva por
    instancia la diferencia verdad/observación para el productor de procedencia
    que DTPD-108 defina, sin serializar aún un formato supuesto.
  - Se comprueba con: `cargo test --test f4_error_quality --locked`; vector
    externo al productor comprueba sustituciones, posición/contexto si procede,
    relación calidad-error, `N` y flujos RNG aislados; puertas comunes.

- ☐ F4.2.3. Publicar e informar las cuatro coberturas distintas, sin convertir
  la cobertura FASTQ final ni la alineada en objetivos nuevos.
  - Depende de: F4.2.1; la alineada requiere DTPD-104/105.
  - Procede de: REQ-DAT-24, REQ-WGS-11, REQ-MAN-37/61 y REQ-VAL-04.
  - Cambia: `src/fastq_validation.rs`, `src/event_manifest.rs` y resultados de
    validación; separa objetivo, WGS básica, FASTQ final propia/contaminante y
    alineada por referencia cuando exista.
  - Se comprueba con: `cargo test --test f4_coverage_reporting --locked`; caso
    controlado demuestra que las cuatro medidas no se sustituyen entre sí;
    puertas comunes.

> **PAUSA-F4.2.** Confirmar que los perfiles cuantitativos se validan contra sus
> datos y que las diferencias aleatorias se informan, no se convierten en
> tolerancias de aceptación.

## 7. F4.3 — Artefactos WGS

- ☐ F4.3.1. Modelar preparación sin PCR por omisión y preparación amplificada
  sólo según los mecanismos y valores ya aprobados.
  - Depende de: F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-TEC-06–08, REQ-CONF-88–90 y REQ-DAT-46–49.
  - Cambia: configuración efectiva y módulos de duplicación; PCR y ópticos se
    distinguen, una duplicación PE afecta a la pareja completa y el estado
    marcado/sin marcar no se finge hasta que DTPD-104 cierre SAM/BAM/CRAM.
  - Se comprueba con: `cargo test --test f4_duplicates --locked`; verdad manual
    comprueba tipos, origen, recuentos y ausencia de huérfanos; puertas comunes.

- ☐ F4.3.2. Añadir adaptadores sin recortar y sesgo GC después del muestreo WGS
  básico, con la excepción explícita de ventana solo `N`.
  - Depende de: F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-50–53, REQ-CONF-91–95, REQ-WGS-08/14 y REQ-MAN-47–48.
  - Cambia: módulo de artefactos y validación; calcula GC como `(G+C)/(A+C+G+T)`
    excluyendo `N`, registra la excepción y no reabre dominancia de configuración.
  - Se comprueba con: `cargo test --test f4_adapters_gc --locked`; oráculo de
    ventanas cubre GC 0/1, mezcla con `N`, todo `N`, extremos e informe por
    intervalo; puertas comunes.

- ☐ F4.3.3. Crear quimeras y anomalías PE desde fragmentos verdaderos, sin usar
  la contaminación como sustituto semántico.
  - Depende de: F4.1.3, F4.2.2 y DTPD-106 aprobada.
  - Procede de: REQ-DAT-60–67, REQ-CONF-99, REQ-MAN-51–52 y REQ-VAL-07.
  - Cambia: módulo de fragmentos/artefactos y datos transitorios de origen; los
    fragmentos quiméricos son del mismo individuo/biblioteca y la pareja anómala
    sigue siendo una unidad completa.
  - Se comprueba con: `cargo test --test f4_chimeras_pairs --locked`; verdad
    manual comprueba uniones, orígenes, orientación, inserto y mates no
    alineados; puertas comunes.

> **PAUSA-F4.3.** Revisar que cada artefacto tiene parámetro/perfil/origen y que
> no altera la verdad biológica. La procedencia durable queda pendiente de
> DTPD-108 y no se sustituye por retener todas las lecturas en memoria.

## 8. F4.4 — Contaminación biológica

- ☐ F4.4.1. Admitir y resolver fuentes compatibles, matriz dirigida o reparto
  automático sin crear otra biblioteca ni mezclar disposiciones de lectura.
  - Depende de: PAUSA-DISEÑO-F4, F4.1.3 y DTPD-108 aprobada.
  - Procede de: REQ-CONF-75–77, REQ-CONT-01–09/15–16 y REQ-DAT-25–28/69/71.
  - Cambia: configuración, planificador de asignación conjunta y manifiesto;
    convierte unidades tras conocer el número limpio, preflighta disponibilidad
    y registra empate/subsemilla. FASTQ externo se limita al modo realista y se
    valida/copia antes de generar.
  - Se comprueba con: `cargo test --test f4_biological_contamination --locked`;
    oráculo de asignación entero verifica mínimos globales, matriz asimétrica,
    0/100 %, PE completo, externo incompatible y ausencia de reemplazo doble;
    puertas comunes.

- ☐ F4.4.2. Materializar sustitución biológica y sus vínculos de instancia sin
  duplicar secuencias/calidades en la procedencia.
  - Depende de: F4.4.1 y DTPD-108 aprobada.
  - Procede de: REQ-CONT-01–05/13–16, REQ-DAT-25–27/34–36/69 y REQ-NOF-08.
  - Cambia: módulos de instancias y procedencia, FASTQ final y manifiesto; la
    receptora recibe ID técnico nuevo enlazado al origen estable y la donante no
    pierde la unidad.
  - Se comprueba con: `cargo test --test f4_contamination_provenance --locked`;
    lector independiente reconstruye origen/mecanismo de FASTQ+procedencia y
    detecta secuencias duplicadas o un par huérfano; puertas comunes.

> **PAUSA-F4.4.** Confirmar que contaminación biológica copia a la receptora y
> demultiplexado técnico transferirá, no copiará, la misma instancia.

## 9. F4.5 — Demultiplexado técnico y FASTQ final

- ☐ F4.5.1. Admitir índices simples, duales combinatorios y UDI, junto con
  tasas, procedencia y política de combinaciones inesperadas.
  - Depende de: PAUSA-DISEÑO-F4, F4.4.2 y DTPD-108 aprobada.
  - Procede de: REQ-CONF-78–81, REQ-DAT-29–32/70–71 y REQ-CONT-10–11.
  - Cambia: configuración de tanda, plan de demultiplexado y validación de
    colisiones; UDI permite reconocer y filtrar, sin asumir que elimina el
    *index hopping*.
  - Se comprueba con: `cargo test --test f4_demultiplexing --locked`; tabla de
    casos cubre traslado, no asignada, filtrada, UDI inesperado, conservación de
    tanda y colisión; puertas comunes.

- ☐ F4.5.2. Publicar un único FASTQ final SE o una pareja final PE por
  individuo/evento y el fichero de procedencia definido, con agregados para
  descartadas/no asignadas.
  - Depende de: F4.5.1 y DTPD-108 aprobada.
  - Procede de: REQ-DAT-17/33–36/40, REQ-RG-02/04/06, REQ-SAL-08/09,
    REQ-MAN-40–41 y REQ-NOF-08.
  - Cambia: `src/fastq_output.rs`, `src/fastq_validation.rs`,
    `src/event_manifest.rs`, ensamblado y nuevo escritor de procedencia por
    flujo. Mantiene gzip/BGZF, validación descomprimida, hash del comprimido y
    no guarda versión limpia.
  - Se comprueba con: `cargo test --test f4_final_fastq_provenance --locked`;
    validador independiente comprueba FASTQ, pareja, hashes, IDs, vínculos y que
    los descartes sólo tienen agregados; puertas comunes.

> **PAUSA-F4.5.** Revisar coherencia transversal FASTQ–procedencia–manifiesto y
> que los descartes reducen cobertura final sin regeneración.

## 10. F4.6 — Modo directo y alineamientos

- ☐ F4.6.1. Admitir modo realista/directo y categorías sintéticas, sin permitir
  FASTQ externo o FASTQ obligatorio en directo.
  - Depende de: PAUSA-DISEÑO-F4, DTPD-108 aprobada y contrato de salida de
    DTPD-104 aprobado.
  - Procede de: REQ-CONF-67–68/97–98/113, REQ-DAT-19–20/57–59/62/67,
    REQ-CONT-12–13 y REQ-MAN-34/50.
  - Cambia: admisión, configuración efectiva, categorías de alineamiento y
    manifiesto; conserva origen/haplotipo y categoría aun cuando posición
    observada difiera.
  - Se comprueba con: `cargo test --test f4_direct_admission --locked`; rechaza
    proporciones inválidas y FASTQ externo, y comprueba trazabilidad del modo;
    puertas comunes.

- ☐ F4.6.2. Generar registros directos internamente coherentes y validar sus
  pares, quimeras, CIGAR, banderas, RNEXT/PNEXT, TLEN y MAPQ contra la categoría
  acordada.
  - Depende de: F4.6.1 y DTPD-104 aprobada.
  - Procede de: REQ-DAT-57–59/62/67, REQ-VAL-07–08 y arquitectura §6.
  - Cambia: nuevo adaptador de alineamiento directo, validador e informe;
    no invoca mapeador ni crea FASTQ.
  - Se comprueba con: `cargo test --test f4_direct_alignment --locked`; oráculo
    SAM/BAM/CRAM independiente valida una muestra de cada categoría y pareja;
    puertas comunes.

- ☐ F4.6.3. Conectar el modo realista al adaptador de mapeador elegido y medir
  cobertura/alineamiento observado sin sobrescribir MAPQ real.
  - Depende de: F4.5.2, DTPD-104 y DTPD-105 aprobadas.
  - Procede de: REQ-DAT-17/21–23/40/54–56/61/66, REQ-VAL-08 y REQ-MAN-34/36.
  - Cambia: adaptador de proceso, conjuntos por referencia, validación y
    manifiesto; usa exactamente FASTQ final, registra herramienta/versión/
    parámetros y conserva no mapeadas.
  - Se comprueba con: `cargo test --test f4_realistic_mapping --locked` y el
    validador externo aprobado por DTPD-104/105; fixture de baja mapeabilidad
    verifica que MAPQ procede del mapeador; puertas comunes.

- ☐ F4.6.4. Comparar funcionalmente modos directo y realista sin exigir
  identidad de decisiones o bytes.
  - Depende de: F4.6.2–F4.6.3.
  - Procede de: REQ-CONF-113, REQ-CONT-13, REQ-VAL-01–13 y REQ-MAN-34.
  - Cambia: informe estructurado y prueba de concordancia; compara objetivos,
    categorías/procedencia y medidas declaradas, preservando que el realista es
    resultado del mapeador y el directo sintético.
  - Se comprueba con: `cargo test --test f4_mode_concordance --locked`; caso
    pequeño con misma semilla/configuración demuestra equivalencia funcional y
    diferencias permitidas explícitas; puertas comunes.

> **PAUSA-F4.6 / revisión de F4.** Aprobar perfiles y artefactos antes de F5.
> Comprobar que toda degradación tiene causa trazada, afecta únicamente evidencia
> técnica, y que el modo directo se etiqueta siempre como sintético.

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
| F4.3 | Artefactos WGS | casos independientes de duplicado, GC, adaptador, quimera y PE | degradación sin origen o alteración de verdad |
| F4.4 | Contaminación biológica | oráculo de asignación conjunta y procedencia | copia/retirada equivocada, exceso o par roto |
| F4.5 | Demultiplexado y FASTQ definitivo | lector independiente FASTQ+procedencia+manifiesto | pérdida de instancia, secuencias duplicadas o salida limpia extra |
| F4.6 | Directo y mapeo real | validador de alineamiento y concordancia funcional | presentar sintético como mapeado o forzar MAPQ |
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
