# Plan de implementación de F2 0.1

## 1. Cabecera y estado

- **Título:** plan detallado de F2, verdad genómica mínima.
- **Versión del plan:** 0.1.
- **Fecha:** 10 de agosto de 2026.
- **Estado:** borrador técnico condicionado; F1.4 está aprobada, pero no está
  autorizado para ejecución hasta que F1.5 y la revisión de cierre de F1 estén
  terminadas y este plan sea aprobado explícitamente.
- **Ámbito temporal:** únicamente F2.1–F2.5.
- **Fuentes normativas y de diseño:**
  [especificaciones-0.1.md](./especificaciones-0.1.md),
  [arquitectura-0.1.md](./arquitectura-0.1.md),
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md),
  [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md),
  [trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md),
  [registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md)
  y [plan-implementacion-f1-0.1.md](./plan-implementacion-f1-0.1.md).
- **Decisiones confirmadas aplicables:** DTC-001, DTC-003, DTC-004, DTC-006,
  DTC-010, DTC-011, DTC-013, DTC-020 a DTC-022, DTC-024, DTC-031 a
  DTC-049.
- **Decisiones pendientes bloqueantes:** DTPD-101 para F2.1/F2.2, DTPD-102
  para F2.1/F2.4 y DTPD-103 para F2.5. También falta documentar y aprobar la
  concreción técnica de Poisson indicada en §4.2. Para los *indels* aleatorios
  de F2 ya se acordó una longitud exacta configurable por separado para
  inserciones y deleciones, con valor predeterminado de 1 pb; falta trasladar
  este acuerdo al registro canónico de decisiones.

Todas las casillas permanecen abiertas. El documento ordena trabajo futuro: no
acredita que F1 haya cerrado, no resuelve decisiones pendientes, no selecciona
crates, versiones, formatos ni herramientas externas y no autoriza commits.

## 2. Resultado y alcance

F2 termina con un recorrido directo pequeño pero completo que, partiendo de una
referencia de F1, construye individuos de ploidía uniforme, conserva copias
haplotípicas separadas, aplica SNP e inserciones/deleciones pequeñas inequívocas
y publica una verdad estructurada canónica junto con exactamente un VCF de verdad
multimuestra por conjunto. La fase, la dosis, la procedencia, las coordenadas y
los alelos deben reconstruirse sin consultar lecturas, alineamientos ni el propio
VCF como fuente de verdad.

### 2.1. Incluido

- Una referencia pequeña y una ploidía uniforme configurable entre 1 y 10 para
  todas sus secuencias en el recorrido de F2; la variación de ploidía por
  secuencia sigue reservada a F5.3.
- Estado haplotípico por individuo, secuencia y copia, con fase, dosis,
  procedencia y materialización exacta, sin duplicar secuencias completas por
  defecto cuando la decisión DTPD-102 elija otra representación.
- SNP e inserciones/deleciones pequeñas inequívocas, explícitas y aleatorias,
  sobre toda una secuencia o regiones declaradas; solo el número aleatorio de
  variantes usa Poisson. La longitud de los *indels* aleatorios de F2 es un
  valor exacto configurable por clase, predeterminado a 1 pb y limitado al
  intervalo desde 1 hasta `umbral estructural efectivo - 1`.
- Coordenadas externas basadas en 1 y cerradas, inserción anclada en la base de
  referencia precedente y solo interna, y deleción interna que deja al menos
  una base a cada lado.
- Conservación literal de `N`: un SNP nunca tiene origen ni destino `N`, una
  inserción aleatoria usa solo `A/C/G/T`, y una inserción explícita o copiada
  puede conservar `N`.
- Configuración declarada, validada y efectiva del recorrido directo, reutilizando
  el resolvedor de reglas, IDs tipados, subsemillas locales, integridad,
  manifiestos y publicación que cierre F1.
- Verdad estructurada obligatoria e independiente del VCF, y un VCF de verdad
  mínimo con GT faseado, dosis correcta, orden de referencia, alelos derivados
  de la verdad y validación externa independiente.

### 2.2. Excluido

- Lecturas, FASTQ, BAM/CRAM, mapeo y cualquier error técnico; comienzan en F3.
- Variantes estructurales, repeticiones, referencias derivadas, ploidía distinta
  por secuencia y la cobertura VCF/complemento completa de casos complejos; F5.
- Frecuencias poblacionales, poblaciones, evolución, recombinación, selección,
  familias y cruces; F6–F8.
- Distribuciones no Poisson para recuentos, distribuciones probabilísticas de
  longitud de *indels*, gVCF con posiciones invariantes y cualquier resolución
  silenciosa de `N`. Las distribuciones configurables de longitud corresponden
  a REQ-MUT-08 en F7.2; F2 solo admite la longitud exacta acordada.
- DTF-003, DTF-008, DTF-013, DTF-014, DTF-023, DTF-033, DTF-034 y DTF-036,
  además de las demás ampliaciones post-v1 del registro vigente.
- Afirmar el cierre total de requisitos asignados a fases posteriores. F2 puede
  dejar evidencia parcial reutilizable para, por ejemplo, REQ-IND-58,
  REQ-IND-59, REQ-VER-01, REQ-VER-24 a REQ-VER-26, REQ-FMT-07 y REQ-SAL-02,
  pero no los marcará satisfechos en su matriz de cierre.

## 3. Orden, reutilización y fronteras

1. Cerrar F1 y congelar su línea base antes de modificar contratos de dominio,
   configuración, referencia, aleatoriedad o publicación.
2. Resolver DTPD-101 y DTPD-102 con prototipos y oráculos independientes antes
   de implementar coordenadas o estado haplotípico; resolver DTPD-103 antes de
   implementar VCF.
3. Mantener un único servicio de coordenadas. Las fronteras traducen la base 1
   cerrada configurada; ningún adaptador ni algoritmo realiza conversiones
   implícitas.
4. Extender `DeclaredConfig` y `EffectiveConfig`, el catálogo `ParameterKey`,
   `ResolutionContext` y los IDs tipados existentes. No crear un segundo parser,
   resolvedor, catálogo de semillas, mecanismo de hash o publicador.
5. Crear primero verdad de dominio y materialización; derivar después el fichero
   estructurado y, por último, el VCF. Ninguna salida técnica se convierte en
   estado canónico.
6. Aplicar eventos explícitos en orden declarado y después clases aleatorias en
   el orden normativo aprobado. Los ámbitos aleatorios se registran antes de
   que exista un consumidor y no dependen del orden físico de ejecución.
7. Introducir cada dependencia mediante DTC-043 y cada cambio de algoritmo,
   catálogo o versión como cambio aislado que repite vectores, reproducibilidad
   y validadores.
8. Terminar cada hito con una pausa obligatoria. Un hito consumidor no comienza
   mientras la revisión o las evidencias del productor estén pendientes.

### 3.1. Línea base observada, no garantía de cierre

El árbol actual ya contiene el parser TOML, `DeclaredConfig`, `EffectiveConfig`,
el resolvedor de reglas con `ParameterKey::VariantRate`, selectores espaciales,
diagnósticos, RNG ChaCha12 local y configuración de Rust 1.97.1/edición 2024.
F1.4 está aprobada y F1.5 continúa pendiente; por tanto, los IDs tipados,
SHA-256, publicación y manifiestos de F1.4 pueden reutilizarse, mientras el
adaptador FASTA solo puede consumirse después de la revisión y cierre de F1.5.
`src/domain.rs` aún no constituye un modelo
genómico. En el manifiesto actual no hay una dependencia aprobada para Poisson,
serialización de verdad ni VCF; este plan no presupone ninguna.

## 4. Prerrequisitos, decisiones y comprobaciones comunes

### 4.1. Puerta de entrada a F2

Antes de A1 deben constar: F1.4 y F1.5 terminadas; revisión de cierre de F1
aprobada; pruebas y evidencias F1 verdes; hashes/esquemas/catálogos F1 congelados;
árbol de trabajo atribuido; y aprobación explícita de este plan. Si F1 cambia
un contrato que F2 pretendía reutilizar, se actualiza primero este plan y su
trazabilidad.

### 4.2. Decisiones que detienen solo a sus consumidores

- DTPD-101 debe elegir la representación interna y demostrar bordes,
  conversiones y composición de segmentos sin pérdidas.
- DTPD-102 debe elegir por separado el estado haplotípico interno y el fichero
  obligatorio de verdad, demostrando reconstrucción, memoria, fase, dosis y
  procedencia.
- DTPD-103 debe fijar versión VCF, normalización, campos `INFO`, representación
  complementaria y validadores independientes, incluidos BGZF e índice.
- La documentación no concreta todavía la implementación Poisson. B1 debe
  crear una entrada de decisión conforme a §9 de
  `decisiones-tecnicas-0.1.md`, evaluar opciones, fijar el algoritmo y aprobarla.
  Esta puerta debe cerrarse antes de B2.
  Debe conservar Poisson como única distribución del número de variantes y
  especificar espacios vacíos, límites y consumo aleatorio. La misma entrada
  debe registrar la regla ya acordada para F2: inserciones y deleciones
  aleatorias reciben longitudes exactas configurables e independientes, con
  valor predeterminado de 1 pb y dominio inequívoco
  `1 ≤ longitud < umbral estructural efectivo`; no se
  sortea la longitud ni se introduce una distribución probabilística hasta
  F7.2.

Que una puerta quede bloqueada no impide preparar matrices, casos literales o
validadores que no dependan de su resultado. Sí impide escribir el contrato de
producción que la consuma.

### 4.3. Puertas de calidad futuras

Estos comandos son criterios de ejecución futura y no se presentan como
ejecutados por este documento:

```bash
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --format-version 1
cargo tree --locked -e features
```

Cada validador externo elegido debe registrarse con versión y comando exactos.
Los oráculos de aceptación no pueden llamar al productor para calcular el
resultado esperado ni leer el VCF para reconstruir la verdad con la que se lo
compara.

## 5. F2.1 — Estado haplotípico mínimo

**Resultado del hito:** coordenadas y estado propios representan cada copia
haplotípica, fase, dosis, procedencia y `N`, y permiten materializar secuencias
exactas desde la referencia cerrada en F1.

- ☐ A1. Congelar y auditar la línea base que F2 consumirá.
  - Depende de: cierre explícito de F1.4, F1.5 y F1, y aprobación de este plan.
  - Procede de: dependencia F1→F2 de la hoja de ruta; DTC-043, DTC-045 y DTC-048.
  - Cambia: evidencia de preflight de F2 con commit base, toolchain, lock, esquemas, catálogos, comandos verdes y contratos F1 reutilizados; no cambia producción.
  - Se comprueba con: revisión del acta de cierre F1, `git status --short`, `rustc --version`, `cargo metadata --locked --format-version 1` y las puertas de §4.3 sobre el commit base.

- ☐ A2. Resolver y aprobar DTPD-101 sin filtrar tipos de dependencias a dominio.
  - Depende de: A1.
  - Procede de: DTPD-101; DTC-010; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: registro de decisiones y prototipo desechable que compara las representaciones candidatas en límites, inserciones ancladas, deleciones internas y composición de ediciones.
  - Se comprueba con: tabla de opciones y costes, vectores literales de primera/última base y rangos vacíos, demostración de ida/vuelta en fronteras, revisión independiente y aprobación fechada.

- ☐ A3. Resolver y aprobar DTPD-102 para estado interno y fichero de verdad como contratos distintos.
  - Depende de: A1 y A2 para interpretar coordenadas del prototipo.
  - Procede de: DTPD-102; arquitectura §§6.2, 9 y 10; REQ-IND-28, REQ-IND-29 y REQ-VER-31.
  - Cambia: registro de decisiones y prototipos comparables con copias, SNP, indels, `N`, materialización, fase, dosis y procedencia; documenta formato y versionado elegidos sin confundirlos con la estructura interna.
  - Se comprueba con: reconstrucción por un comprobador independiente, medida reproducible de memoria en el caso pequeño acordado, casos poliploides 1–10 y aprobación explícita conforme a §9 del registro técnico.

- ☐ A4. Evaluar cada dependencia necesaria para las decisiones A2–A3 antes de incorporarla.
  - Depende de: A2–A3; se omite si las opciones aprobadas no necesitan dependencia nueva.
  - Procede de: DTC-043, DTC-047 y DTC-048.
  - Cambia: evaluación individual con propósito, REQ, versión exacta, licencia, mantenimiento, seguridad, características, árbol, código nativo, scripts y coste de distribución; después, solo las dependencias aprobadas y `Cargo.lock`.
  - Se comprueba con: revisión de la evaluación, `cargo metadata --locked --format-version 1`, `cargo tree --locked -e features` y repetición de §4.3.

- ☐ A5. Implementar tipos propios de coordenada y el único servicio de conversión.
  - Depende de: A2 y A4 cuando corresponda.
  - Procede de: DTC-010, DTC-013 y DTC-039; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: dominio y servicio de coordenadas según DTPD-101, con constructores validados para posición, intervalo, ancla de inserción y deleción interna; no expone tipos de crates.
  - Se comprueba con: pruebas de tabla escritas desde la semántica externa base 1 cerrada, propiedades de ida/vuelta y fallos explícitos en cero, fuera de rango, desbordamiento y límites terminales.

- ☐ A6. Definir bases y vistas de referencia que preserven `N` literalmente.
  - Depende de: A1, A5 y contrato FASTA cerrado en F1.5.
  - Procede de: DTC-032; REQ-REF-14 y REQ-IND-28.
  - Cambia: tipos de dominio para base y acceso inmutable por secuencia/coordenada, reutilizando la referencia F1 sin reinterpretar ni sustituir `N`.
  - Se comprueba con: referencia sintética que contiene `N` en principio, interior y final, comparación byte a byte con la referencia F1 y casos de error de secuencia o coordenada inexistente.

- ☐ A7. Modelar individuo, secuencia y copia haplotípica con ploidía uniforme explícita.
  - Depende de: A3, A5–A6 y catálogo de IDs cerrado en F1.4.
  - Procede de: DTC-003 y DTC-039; REQ-IND-28, REQ-IND-29 y REQ-VER-27; alcance deliberado de F2.
  - Cambia: tipos de estado que contienen exactamente una copia por índice de fase entre 1 y la ploidía 1–10, sin habilitar ploidía distinta por secuencia.
  - Se comprueba con: constructores y pruebas para ploidías 1, 2 y 10, rechazo de 0 y 11, unicidad de IDs/índices y cardinalidad exacta por secuencia.

- ☐ A8. Representar ediciones y procedencia sin duplicación innecesaria de la referencia.
  - Depende de: A3 y A5–A7.
  - Procede de: DTPD-102; DTC-003, DTC-031 y DTC-039; REQ-IND-28 y REQ-IND-29.
  - Cambia: estado aprobado de segmentos/eventos por copia con alelo de referencia, alelo verdadero, origen explícito o aleatorio y vínculo a la referencia; no añade clases fuera de F2.
  - Se comprueba con: casos literales de copia intacta, SNP, inserción y deleción, inspección de que cada evento conserva identidad/procedencia y medida que demuestra que una copia intacta no clona toda la secuencia.

- ☐ A9. Materializar una copia y calcular fase y dosis desde el estado canónico.
  - Depende de: A7–A8.
  - Procede de: DTC-003 y DTC-031; REQ-IND-28, REQ-IND-29 y REQ-VER-27.
  - Cambia: operaciones de materialización determinista y consulta genotípica por locus; fase y dosis se derivan de copias ordenadas, no de una salida VCF.
  - Se comprueba con: oráculo manual que aplica ediciones a cadenas cortas, casos haploide/diploide/poliploide y propiedad de que dosis equivale al número exacto de copias portadoras.

- ☐ A10. Completar la evidencia y trazabilidad de F2.1.
  - Depende de: A5–A9 y puertas de calidad verdes.
  - Procede de: F2.1; REQ-REF-14, REQ-IND-28, REQ-IND-29 y REQ-VER-27; ACE-03, ACE-10, ACE-11, ACE-28 y ACE-66.
  - Cambia: pruebas de integración, matriz de evidencia F2.1, casos de aceptación y documentación de límites frente a F5.3.
  - Se comprueba con: §4.3, recuento automático de cuatro REQ asignados sin huérfanos y revisión cruzada de materialización, fase, dosis, procedencia y conservación de `N`.

> **Pausa obligatoria F2.1.** Revisar y aprobar DTPD-101/DTPD-102, tipos de
> coordenada, representación haplotípica, consumo de memoria, oráculos y los
> cuatro REQ. No iniciar B6 ni F2.3/F2.4 si la representación no reconstruye
> exactamente cada copia o si una conversión de coordenadas queda implícita.

## 6. F2.2 — Variantes pequeñas explícitas y aleatorias

**Resultado del hito:** el motor aplica SNP e indels pequeños inequívocos con
coordenadas, orden, regiones y conteos correctos, y rechaza cualquier caso que
pretenda convertir `N` en una base conocida sin declaración explícita válida.

- ☐ B1. Registrar y aprobar la concreción técnica de generación aleatoria mínima de F2 antes de diseñar sus consumidores.
  - Depende de: pausa F2.1 aprobada, A2 y A3.
  - Procede de: REQ-CONF-23 a REQ-CONF-25 y REQ-IND-16 a REQ-IND-18; DTC-022, DTC-037, DTC-043 y DTC-049; §4.2 de este plan.
  - Cambia: nueva entrada con identificador asignado por el registro técnico que fija implementación Poisson, selección espacial uniforme sobre candidatos válidos, longitud exacta configurable por clase con predeterminado de 1 pb, contenido `A/C/G/T` de inserciones, límites, espacios vacíos y consumo aleatorio; evalúa conforme a DTC-043 cualquier dependencia candidata y aplaza las distribuciones probabilísticas de longitud a F7.2.
  - Se comprueba con: opciones y dependencias evaluadas, vectores u oráculo independiente para Poisson, prueba estadística con umbral fijado antes de ejecutarla, casos de espacio vacío, longitudes 1/2/`umbral - 1` y rechazo de 0 o `>= umbral`, más aprobación conforme a §9 de decisiones técnicas.

- ☐ B2. Especificar la extensión de autoría TOML mínima de variantes F2 antes de codificarla.
  - Depende de: B1.
  - Procede de: REQ-CONF-08, REQ-CONF-20, REQ-CONF-21, REQ-CONF-22, REQ-CONF-24, REQ-CONF-25 y REQ-CONF-28; DTC-034 y DTC-049.
  - Cambia: esquema documental de configuración con forma inequívoca para eventos SNP/inserción/deleción, ámbito completo o región, objetivo individual/copia, modo explícito o aleatorio, tasa, longitud exacta aleatoria independiente para inserción/deleción y orden; la longitud se omite para obtener 1 pb y otras clases continúan rechazadas o reservadas.
  - Se comprueba con: ejemplos TOML válidos e inválidos, revisión contra coordenadas base 1 cerradas, longitudes omitida/1/2/`umbral - 1` y rechazo de 0 o `>= umbral`; no existen campos que permitan distribuciones no Poisson, distribuciones probabilísticas de longitud o ploidía por secuencia.

- ☐ B3. Extender configuración declarada y parser conservando ubicaciones de origen.
  - Depende de: B2.
  - Procede de: DTC-034, DTC-036 y DTC-044; REQ-CONF-08, REQ-CONF-20 a REQ-CONF-22 y REQ-DIR-03.
  - Cambia: tipos `DeclaredConfig`, adaptador TOML y diagnósticos para los campos aprobados, con orden de eventos y literales originales; no construye estado genómico.
  - Se comprueba con: pruebas de parsing para cada forma, claves desconocidas/duplicadas, clase/versión erróneas y snapshots estructurados de ruta, campo y ubicación.

- ☐ B4. Extender catálogo de parámetros y resolución efectiva solo donde lo exijan los modelos F2.
  - Depende de: B2–B3.
  - Procede de: DTC-006, DTC-033 y DTC-049; REQ-CONF-20 a REQ-CONF-22, REQ-CONF-24, REQ-CONF-28, REQ-IND-09, REQ-IND-10, REQ-IND-13 y REQ-IND-15.
  - Cambia: `ParameterKey`, metadatos de familia/ejes/dominio, solicitudes de resolución y configuración efectiva para distinguir los parámetros realmente necesarios por SNP/inserción/deleción, incluidas las dos longitudes exactas, sin duplicar el resolvedor.
  - Se comprueba con: pruebas de catálogo completas, dominancia por individuo y espacio, cero explícito que sustituye a una tasa heredada, longitudes de inserción y deleción resueltas por separado, ausencia de fallback ambiguo y rechazo de ejes no admitidos.

- ☐ B5. Ampliar el catálogo normativo de subsemillas antes de generar variantes.
  - Depende de: B1 y catálogo de aleatoriedad cerrado en F1.3.
  - Procede de: DTC-011, DTC-037 y DTC-038; REQ-CONF-26 y REQ-CONF-27.
  - Cambia: ámbitos/versiones y vectores BLAKE3 para individuo, secuencia, copia y clase de variante según la decisión aprobada, sin reutilizar etiquetas ad hoc.
  - Se comprueba con: vectores literales independientes, unicidad de claves lógicas, repetición idéntica con la misma semilla y cambio localizado al variar una sola identidad.

- ☐ B6. Construir y contar oportunidades elegibles completas y regionales con densidad cero válida.
  - Depende de: A5–A6, B4 y B5.
  - Procede de: REQ-CONF-21, REQ-CONF-23 y REQ-CONF-28; REQ-IND-13, REQ-IND-14, REQ-IND-19, REQ-IND-38, REQ-IND-40 y REQ-IND-41; REQ-DIR-03.
  - Cambia: enumeración de oportunidades válidas por secuencia, clase y región después de resolver reglas, excluyendo `N`, bordes imposibles y objetivos inválidos; expone su cardinalidad sin seleccionar todavía posiciones ni sortear recuentos.
  - Se comprueba con: enumeración independiente de oportunidades en secuencias cortas, pruebas de primera/última coordenada, región de una base, región vacía tras filtros y sustitución local de tasa por individuo.

- ☐ B7. Implementar conteos Poisson sobre el número ya calculado de oportunidades válidas, sin sorteo cuando la tasa es cero.
  - Depende de: B1, B4–B6.
  - Procede de: DTC-022 y DTC-037; REQ-CONF-23 a REQ-CONF-25 y REQ-IND-16 a REQ-IND-18.
  - Cambia: servicio de conteo que recibe tasa efectiva y cardinalidad elegible calculada por B6, registra esperado/realizado y usa exclusivamente la implementación aprobada; cero o espacio vacío producen cero sin consumir RNG.
  - Se comprueba con: vectores Poisson independientes, semillas fijas, límites grandes/pequeños, media/varianza sobre repeticiones con criterio predeclarado, igualdad `media = tasa × oportunidades válidas` y prueba de consumo nulo para tasa o espacio cero.

- ☐ B8. Validar y aplicar SNP explícitos en el orden declarado.
  - Depende de: A8–A9 y B3–B4.
  - Procede de: REQ-IND-05, REQ-IND-29 y REQ-REC-24; DTC-020 y DTC-032; alcance F2.2 de la hoja de ruta.
  - Cambia: evento SNP explícito con referencia esperada, alelo alternativo `A/C/G/T`, destino de copia y error contextual; rechaza origen o destino `N`.
  - Se comprueba con: casos literales de transición/transversión, desacuerdo de alelo, copia ausente, coordenada límite y SNP sobre `N`, verificando que todo rechazo deja el estado sin cambio.

- ☐ B9. Generar SNP aleatorios solo sobre candidatos y alelos válidos.
  - Depende de: B5–B7.
  - Procede de: REQ-IND-04, REQ-IND-17 a REQ-IND-19 y REQ-IND-28; DTC-032 y DTC-037; alcance F2.2 de la hoja de ruta.
  - Cambia: generador de posición y ALT que excluye origen `N`, nunca produce destino `N` y conserva el orden lógico independiente del planificador.
  - Se comprueba con: referencia exhaustiva corta con `A/C/G/T/N`, oráculo de candidatos, vectores con semilla fija y prueba estadística predeclarada de uniformidad sin exigir una secuencia concreta fuera de los vectores normativos.

- ☐ B10. Validar y aplicar inserciones pequeñas explícitas y aleatorias.
  - Depende de: A5, A8–A9, B1, B3 y B5–B7.
  - Procede de: REQ-IND-05, REQ-IND-29, REQ-IND-40 y REQ-REC-24; DTC-032.
  - Cambia: evento de inserción interna anclado en la base precedente; la forma explícita conserva literalmente `N` y la aleatoria genera únicamente `A/C/G/T` con la longitud exacta efectiva de inserción.
  - Se comprueba con: casos manuales de ancla y secuencia resultante, rechazo de inserción terminal imposible, inserción explícita con `N` conservada, longitudes aleatorias exactas 1/2/`umbral - 1` y barrido que no contiene `N`.

- ☐ B11. Validar y aplicar deleciones pequeñas solo cuando quedan bases a ambos lados.
  - Depende de: A5, A8–A9, B1, B3–B7.
  - Procede de: REQ-IND-05, REQ-IND-29, REQ-IND-41 y REQ-REC-24.
  - Cambia: evento de deleción con intervalo y alelo de referencia explícitos; la forma aleatoria usa la longitud exacta efectiva de deleción; se rechazan longitud cero, desacuerdo, solapamiento inválido y cualquier deleción terminal.
  - Se comprueba con: oráculo de cadenas cortas, longitudes aleatorias exactas 1/2/`umbral - 1`, casos que intentan borrar primera/última/toda la secuencia, verificación de bases flanqueantes y ausencia de mutación tras error.

- ☐ B12. Componer eventos sobre el estado vigente con orden y conflictos explícitos.
  - Depende de: B8–B11.
  - Procede de: DTC-020 y DTC-031; REQ-REC-24 a REQ-REC-27.
  - Cambia: aplicación secuencial de explícitos y después aleatorios según orden de clase aprobado, con coordenada de referencia trazable, alelo vigente comprobado y diagnóstico si el objetivo ya no existe o el conflicto no es representable.
  - Se comprueba con: pares SNP→deleción, deleción→SNP, inserciones vecinas y eventos repetidos cuyo resultado esperado se calcule manualmente; permutar el orden debe cambiar solo los casos semánticamente ordenados.

- ☐ B13. Cerrar diagnósticos y validación anticipada de entradas F2.2.
  - Depende de: B3–B12.
  - Procede de: DTC-044; REQ-CONF-25, REQ-CONF-28, REQ-IND-40, REQ-IND-41, REQ-REC-27 y REQ-DIR-03.
  - Cambia: catálogo estable de errores de variante, agregado de fallos independientes y validación antes de publicar, incluidos tasa/distribución no admitida, coordenada, alelo, `N`, región, copia y conflicto.
  - Se comprueba con: tabla diagnóstico→caso, texto y JSON equivalentes, orden determinista de errores y prueba de que una configuración inválida no crea artefactos finales.

- ☐ B14. Completar evidencia y trazabilidad de F2.2 con oráculos ajenos al productor.
  - Depende de: B5–B13 y puertas de calidad verdes.
  - Procede de: los 22 REQ asignados a F2.2 en la matriz de §13; ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109 y ACE-111.
  - Cambia: pruebas unitarias, de propiedades, estadísticas y de integración; corpus sintético versionado; matriz de evidencia F2.2 con esperado calculado fuera del motor.
  - Se comprueba con: §4.3, recuento automático de 22 REQ sin huérfanos, caso negativo que rechaza una frecuencia alélica poblacional como valor literal de un individuo para REQ-IND-11, revisión independiente de los vectores y repetición exacta de verdad lógica con igual semilla. REQ-CONF-08, REQ-CONF-21, REQ-IND-02 y REQ-IND-03 no se marcan cerrados hasta resolver RAF2-004.

> **Pausa obligatoria F2.2.** Revisar generación Poisson, espacios elegibles,
> orden, bordes de indels y todos los casos con `N`. No iniciar F2.3 si un
> oráculo se deriva del motor, si una tasa cero consume RNG o si un rechazo deja
> una edición parcial.

## 7. F2.3 — Creación directa de individuos

**Resultado del hito:** la vía directa construye uno o varios individuos desde
la referencia y su configuración efectiva, aplica cambios por copia y permite
reconstruir exactamente cada haplotipo sin crear población ni escenario.

- ☐ C1. Extender la autoría de individuos directos sin anticipar estructuras poblacionales.
  - Depende de: pausa F2.2 aprobada y B1–B3.
  - Procede de: REQ-CONF-20, REQ-CONF-22, REQ-DIR-01 y REQ-DIR-02; DTC-034.
  - Cambia: esquema documental de `individuals` para ruta `direct`, referencia, ploidía uniforme, eventos/cambios por copia y ajustes aleatorios permitidos; no añade población, parentesco ni ploidía por secuencia.
  - Se comprueba con: ejemplos de cero, uno y varios individuos, revisión contra el alcance F2 y errores para escenario presente, referencia desconocida, copia fuera de ploidía y campos futuros.

- ☐ C2. Parsear, validar y resolver individuos directos con IDs tipados y procedencia.
  - Depende de: C1 y contratos de identidad/configuración cerrados en F1.
  - Procede de: DTC-033, DTC-039 y DTC-044; REQ-CONF-20, REQ-CONF-22, REQ-DIR-02 y evidencia compartida de REQ-IND-25/REQ-MAN-10.
  - Cambia: configuración declarada y efectiva de individuos, catálogo de identidades y solicitudes de resolución por individuo, conservando valor, regla ganadora, selectores y origen.
  - Se comprueba con: IDs repetidos/desconocidos, orden de declaración, sustitución individual y prueba de que dos individuos no comparten accidentalmente una identidad o valor efectivo mutable.

- ☐ C3. Construir el estado inicial faseado de cada individuo desde la referencia.
  - Depende de: A6–A9 y C2.
  - Procede de: DTC-003 y DTC-006; REQ-DIR-01, REQ-DIR-02, REQ-IND-06, REQ-IND-28 y REQ-IND-29.
  - Cambia: caso de uso de creación directa que instancia copias intactas por secuencia/ploidía y enlaza referencia, individuo y procedencia antes de aplicar variantes.
  - Se comprueba con: individuo sin variantes cuyas copias materializadas igualan byte a byte la referencia, cardinalidad de copias 1–10 y ausencia de escenarios o poblaciones implícitas.

- ☐ C4. Aplicar todos los eventos explícitos de un individuo en orden declarado.
  - Depende de: B8, B10–B13 y C3.
  - Procede de: DTC-020; REQ-IND-05, REQ-REC-24, REQ-DIR-01 y REQ-DIR-02.
  - Cambia: orquestación directa de eventos explícitos por individuo/copia, con transacción lógica que conserva el estado previo si la validación global falla.
  - Se comprueba con: configuración corta cuyo resultado por copia se calcula manualmente, permutación controlada del orden y fallo tardío que no publica ni devuelve un individuo parcialmente válido.

- ☐ C5. Aplicar variantes aleatorias con ámbitos locales independientes por individuo, secuencia, copia y clase.
  - Depende de: B5–B9, B10–B12 y C3–C4.
  - Procede de: DTC-011, DTC-037 y DTC-038; REQ-CONF-22 y REQ-IND-15.
  - Cambia: orquestación de conteos y candidatos aleatorios después de los explícitos, sin compartir RNG mutable entre identidades lógicas.
  - Se comprueba con: añadir o quitar un individuo no altera los restantes, desactivar una clase no altera las demás, ejecutar en distinto orden físico conserva la verdad y los vectores normativos permanecen estables.

- ☐ C6. Materializar y validar cada secuencia haplotípica final contra su historial.
  - Depende de: C4–C5 y A9.
  - Procede de: DTC-003 y DTC-031; REQ-IND-15 y REQ-DIR-01.
  - Cambia: resultado de creación directa con acceso a secuencia final, eventos, fase, dosis y procedencia, sin convertir la secuencia materializada en una segunda verdad divergente.
  - Se comprueba con: dos reconstructores —aplicación de eventos y ensamblado del estado aprobado— que producen los mismos bytes, longitudes esperadas y `N` conservadas en cada copia.

- ☐ C7. Integrar múltiples individuos y sustituciones locales sin acoplamiento.
  - Depende de: C2–C6.
  - Procede de: DTC-006 y DTC-038; REQ-CONF-22, REQ-IND-15 y REQ-DIR-01.
  - Cambia: recorrido directo completo para un conjunto ordenado de individuos con valores efectivos y verdad separados por identidad.
  - Se comprueba con: caso de al menos tres individuos donde solo uno sustituye una tasa/región, comparación de los otros dos con una ejecución sin ese individuo y revisión de orden estable por ID normativo.

- ☐ C8. Completar evidencia y trazabilidad de F2.3.
  - Depende de: C2–C7 y puertas de calidad verdes.
  - Procede de: REQ-CONF-20, REQ-CONF-22, REQ-IND-15, REQ-DIR-01 y REQ-DIR-02; evidencia compartida de REQ-IND-25 y REQ-MAN-10; ACE-01, ACE-10, ACE-12, ACE-22 y ACE-25.
  - Cambia: pruebas de integración de la vía directa, matriz de evidencia F2.3 y escenario sintético reproducible con verdad manual por haplotipo; registra por separado la evidencia parcial de los REQ compartidos.
  - Se comprueba con: §4.3, recuento automático de cinco REQ de cierre y dos REQ compartidos sin huérfanos, asignación real de ID estable cuando se omite y registro de valores globales/específicos/efectivos, además de revisar que el recorrido no crea escenarios, poblaciones, lecturas ni formatos técnicos.

> **Pausa obligatoria F2.3.** Reconstruir manualmente todos los haplotipos del
> escenario directo de aceptación y revisar aislamiento de subsemillas. No
> iniciar F2.4 si la secuencia final no se deriva de la misma verdad canónica o
> si añadir un individuo altera a otro.

## 8. F2.4 — Verdad estructurada independiente del VCF

**Resultado del hito:** una ejecución directa válida publica un fichero de
verdad estructurada versionado que basta para reconstruir referencia, individuos,
copias, fase, dosis, eventos y procedencia aunque no se produzca ningún VCF.

- ☐ D1. Formalizar el esquema externo aprobado por DTPD-102 y su contrato de compatibilidad.
  - Depende de: pausa F2.3 aprobada y A3.
  - Procede de: DTPD-102; DTC-031, DTC-035 y DTC-036; REQ-VER-31 y REQ-SAL-06.
  - Cambia: esquema versionado del fichero de verdad con campos obligatorios, tipos, orden normativo, coordenadas, unidades y política ante versión/claves desconocidas; no reutiliza el esquema del estado interno.
  - Se comprueba con: ejemplos válidos e inválidos, recorrido campo→REQ y validación mediante una implementación de esquema independiente elegida y registrada.

- ☐ D2. Implementar una proyección pura desde la verdad canónica al modelo externo.
  - Depende de: C6–C7 y D1.
  - Procede de: DTC-003, DTC-004 y DTC-031; REQ-VER-31 y evidencia compartida de REQ-IND-26.
  - Cambia: transformación que incluye referencia/huella, IDs, secuencias, copias, ploidía, fase, dosis, eventos, alelos, coordenadas, `N`, origen y subsemilla o procedencia necesaria, sin consultar VCF.
  - Se comprueba con: prueba exhaustiva campo a campo sobre el escenario sintético y propiedad de que desactivar la salida VCF deja idéntica la proyección de verdad.

- ☐ D3. Escribir la verdad de forma determinista y acotada por el mecanismo de publicación F1.
  - Depende de: D1–D2 y cierre de publicación/integridad F1.4.
  - Procede de: DTC-035, DTC-040 y DTC-045; REQ-SAL-06.
  - Cambia: adaptador del formato aprobado, artefacto temporal, hash SHA-256 y promoción atómica/no-clobber; el orden físico estable no concede semántica adicional.
  - Se comprueba con: dos escrituras lógicas iguales con bytes/hash iguales cuando el contrato lo exija, fallo inyectado sin destino completo y `sha256sum` comparado con el manifiesto.

- ☐ D4. Construir un lector/comprobador independiente del productor para la aceptación.
  - Depende de: D1; no depende de D2–D3 para calcular resultados esperados.
  - Procede de: DTC-047; REQ-VER-31 y REQ-SAL-06; ACE-34 y ACE-62.
  - Cambia: herramienta o implementación de prueba separada que valida esquema, relaciones y reconstruye haplotipos desde el fichero y la referencia; registra versión/comando si es externa.
  - Se comprueba con: corpus escrito manualmente, mutaciones deliberadas de cada invariante y comparación byte a byte contra secuencias esperadas precomputadas, no contra una nueva llamada al productor.

- ☐ D5. Demostrar que la verdad permanece disponible sin ninguna salida VCF.
  - Depende de: D2–D4.
  - Procede de: DTC-004 y DTC-031; REQ-VER-31 y REQ-SAL-06.
  - Cambia: opción efectiva de salida y recorrido que publica verdad/manifiesto aun con VCF desactivado, sin cambiar el genotipo ni eliminar detalles necesarios.
  - Se comprueba con: ejecuciones gemelas VCF activado/desactivado que tienen la misma verdad lógica y reconstrucción, y ausencia verificada de dependencia desde D2 hacia cualquier módulo VCF.

- ☐ D6. Integrar verdad, configuración efectiva, referencia y manifiestos mediante identidades e integridad.
  - Depende de: D3 y contratos cerrados de F1.4/F1.5.
  - Procede de: DTC-039, DTC-040 y DTC-045; REQ-VER-31, REQ-SAL-06 y evidencia compartida de REQ-IND-26.
  - Cambia: manifiestos de ejecución/evento y catálogo de artefactos para enlazar fichero de verdad, esquema, referencia exacta, configuración efectiva, hashes y estado de publicación.
  - Se comprueba con: recorrido bidireccional de cada ID/ruta/hash, `sha256sum`, rechazo de referencia o artefacto intercambiado y revisión de que no existen rutas huérfanas.

- ☐ D7. Cubrir corrupción, truncado, versión y pérdida de campos sin degradación silenciosa.
  - Depende de: D3–D6.
  - Procede de: DTC-031, DTC-036, DTC-044 y DTC-045; REQ-VER-31 y REQ-SAL-06.
  - Cambia: validación y diagnósticos de lectura/publicación para bytes corruptos, JSON incompleto o equivalente del formato elegido, versión desconocida, duplicados y relaciones rotas.
  - Se comprueba con: mutación sistemática del corpus, códigos/ubicaciones estables y prueba de que ningún fallo se acepta como verdad parcial ni deja una ejecución marcada completa.

- ☐ D8. Completar evidencia y trazabilidad de F2.4.
  - Depende de: D1–D7 y puertas de calidad verdes.
  - Procede de: REQ-VER-31 y REQ-SAL-06; evidencia compartida de REQ-IND-26; ACE-14, ACE-34 y ACE-62.
  - Cambia: corpus dorado pequeño, pruebas de integración, informe del validador independiente y matriz de evidencia F2.4 que separa cierre y evidencia compartida.
  - Se comprueba con: §4.3, recuento automático de dos REQ de cierre y un REQ compartido sin huérfanos, validador de esquema independiente, coherencia del mismo ID entre estado y fichero, y reconstrucción exacta de todas las copias del corpus.

> **Pausa obligatoria F2.4.** Entregar el fichero y la referencia a una revisión
> que no use el productor y exigir reconstrucción exacta. No iniciar F2.5 si el
> VCF es necesario para interpretar la verdad, si falta fase/dosis/procedencia o
> si una corrupción puede parecer una ejecución completa.

## 9. F2.5 — VCF de verdad mínimo

**Resultado del hito:** cada conjunto publica exactamente un VCF de verdad
multimuestra, ordenado por referencia, derivado de la verdad estructurada y
validado de forma independiente, con GT faseado y sin campos observacionales
inventados.

- ☐ E1. Resolver y aprobar DTPD-103 con los casos mínimos y futuros explícitamente separados.
  - Depende de: pausa F2.4 aprobada.
  - Procede de: DTPD-103; DTC-021; REQ-FMT-05, REQ-FMT-06, REQ-VER-06 y REQ-VER-34.
  - Cambia: registro técnico con versión VCF, normalización, campos `INFO`, BGZF/índice, representación complementaria para pérdidas de F2, herramientas de validación y límites que quedan en F5.5.
  - Se comprueba con: evaluación de capacidad haploide/poliploide, SNP/indels/multialélicos/`N`, prototipos interoperables, reconstrucción sin invención y aprobación conforme a §9 de decisiones técnicas.

- ☐ E2. Evaluar e incorporar solo dependencias y validadores aprobados para VCF.
  - Depende de: E1.
  - Procede de: DTC-043, DTC-047 y DTC-048.
  - Cambia: evaluaciones con versión exacta y, tras aprobarlas, manifiesto/lock y registro de versión/comando de cada validador; no selecciona herramientas por costumbre.
  - Se comprueba con: revisión DTC-043, `cargo metadata --locked --format-version 1`, `cargo tree --locked -e features`, §4.3 y ejecución del prototipo con el validador independiente registrado.

- ☐ E3. Definir la proyección VCF desde verdad final y referencia exacta.
  - Depende de: E1–E2 y D2.
  - Procede de: DTC-003, DTC-004 y DTC-021; REQ-VER-06, REQ-VER-15, REQ-VER-16, REQ-VER-17 y REQ-VER-30.
  - Cambia: modelo intermedio propio con contigs, muestras, posiciones, REF, ALT y genotipos derivados únicamente del estado canónico; no parsea secuencias de lecturas ni infiere observaciones.
  - Se comprueba con: fixtures construidos a mano, prueba de dependencia de módulos y comparación de cada registro con referencia/verdad mediante un oráculo distinto.

- ☐ E4. Normalizar y ordenar SNP e indels según la decisión aprobada.
  - Depende de: E1 y E3.
  - Procede de: REQ-FMT-05, REQ-FMT-06, REQ-VER-15 y REQ-VER-16; DTC-013 y DTC-021.
  - Cambia: representación REF/ALT mínima para SNP, inserción y deleción, conversión exclusiva en el servicio de coordenadas y orden por catálogo de referencia, posición y desempate aprobado.
  - Se comprueba con: casos manuales en primera/interior/última zona admitida, repetición local, ancla de indel y comparación contra el normalizador independiente fijado por E1.

- ☐ E5. Agrupar alelos por sitio y calcular GT faseado/dosis para todas las muestras.
  - Depende de: E3–E4.
  - Procede de: DTC-003 y DTC-021; REQ-VER-02, REQ-VER-12, REQ-VER-17, REQ-VER-30, REQ-VER-32, REQ-VER-33 y REQ-VER-34.
  - Cambia: agrupación multialélica con ALT fijos incluidos, índices de alelo estables, una columna por individuo y GT con exactamente una entrada por copia en orden de fase.
  - Se comprueba con: loci bialélicos/multialélicos, dosis 0..P y ploidías 1, 2 y 10, oráculo de conteo por copia y permutación de muestras que no cambia alelos ni coordenadas.

- ☐ E6. Omitir invariantes y alelos eliminados sin perder la verdad canónica.
  - Depende de: E3–E5.
  - Procede de: REQ-REC-26, REQ-VER-30, REQ-VER-36, REQ-VER-37, REQ-VER-39, REQ-VER-40 y REQ-VER-41; DTC-021.
  - Cambia: filtros de proyección que excluyen posiciones invariantes, variantes ausentes del conjunto y SNP borrados por una deleción posterior, manteniendo el historial completo en el fichero estructurado.
  - Se comprueba con: verdad manual que contiene alelo fijo, evento revertido o eliminado y sitio invariante; el VCF contiene solo estados finales representables y el fichero D conserva todo el historial exigido.

- ☐ E7. Aplicar la regla de `N` y emitir el complemento mínimo aprobado cuando el VCF pierda detalle.
  - Depende de: E1, E3–E6.
  - Procede de: DTC-021 y DTC-032; REQ-VER-24, REQ-VER-34, REQ-VER-47 y REQ-MAN-27; cierre F2 de la hoja de ruta.
  - Cambia: rechazo de cualquier SNP con `REF=N` y conservación de inserciones explícitas con `N` en verdad; solo crea complemento si DTPD-103 identifica una pérdida real del alcance F2. Si ningún caso F2 pierde detalle, no se inventa una entrada: se conserva la prueba negativa y RAF2-007 impide aprobar el cierre hasta corregir la evidencia de la hoja de ruta.
  - Se comprueba con: casos `N` de referencia, inserción explícita y aleatoria, ausencia de complemento para todo caso plenamente representable y prueba negativa de que ninguna base desconocida se sustituye por `A/C/G/T`; una reconstrucción conjunta VCF+complemento solo se exige si DTPD-103 documenta antes un caso F2 literal con pérdida real.

- ☐ E8. Escribir exactamente un VCF multimuestra por conjunto con cabecera y campos mínimos correctos.
  - Depende de: E2–E7.
  - Procede de: DTC-021; REQ-FMT-05, REQ-FMT-06, REQ-VER-06, REQ-VER-43, REQ-VER-44 y REQ-MAN-57.
  - Cambia: adaptador VCF con versión/cabecera aprobadas, `QUAL=.`, `FILTER=PASS`, GT faseado y ausencia de `DP`, `AD`, `GQ`, `PL` u otros campos observacionales no justificados.
  - Se comprueba con: parser independiente de cabecera/registros, búsqueda negativa de campos prohibidos, cardinalidad de un fichero por conjunto y comparación exacta del orden de muestras/contigs.

- ☐ E9. Publicar BGZF e índice por defecto y texto solo cuando se solicite explícitamente.
  - Depende de: E1–E2, E8 y publicación F1.4 cerrada.
  - Procede de: REQ-SAL-12 y REQ-MAN-58 como fuentes directas de BGZF, índice y modo texto; REQ-SAL-21 a REQ-SAL-23 y evidencia compartida de REQ-MAN-83 para los estados y la publicación; DTC-040 y DTC-045.
  - Cambia: artefactos temporales, compresión/índice aprobados, hashes y manifiestos con modo efectivo; todos se promueven como unidad coherente o la ejecución falla.
  - Se comprueba con: validador externo que abre BGZF e índice, consultas por región, modo texto sin índice según contrato, `sha256sum` y fallos inyectados sin conjunto parcialmente completo.

- ☐ E10. Ejecutar validación sintáctica independiente y validación semántica contra la verdad.
  - Depende de: E7–E9.
  - Procede de: DTC-047; REQ-VAL-28, REQ-VER-02, REQ-VER-12, REQ-VER-43, REQ-VER-44 y REQ-VER-47.
  - Cambia: pruebas/CI con validador VCF externo registrado y comprobador semántico separado que lee referencia, verdad y VCF/complemento, pero no llama a la proyección productora.
  - Se comprueba con: comando externo exacto sobre texto y BGZF/índice, mutaciones deliberadas de REF/ALT/GT/fase/dosis/orden y detección de cada discordancia por el comprobador independiente.

- ☐ E11. Verificar el recorrido completo configuración→verdad→VCF y sus invariantes cruzadas.
  - Depende de: E3–E10.
  - Procede de: F2.5 y cierre F2; DTC-003, DTC-004, DTC-021 y DTC-031.
  - Cambia: escenario de aceptación con varias muestras, ploidías del alcance, SNP, inserción/deleción, sitio multialélico, alelo fijo y `N`, más inventario de artefactos y métricas pequeñas de tiempo/memoria.
  - Se comprueba con: reconstrucción independiente de todos los haplotipos, validación VCF externa, contraste fase/dosis/alelos/coordenadas y repetición con igual semilla que conserva verdad lógica e integridad esperada.

- ☐ E12. Completar evidencia y trazabilidad de F2.5 y preparar la revisión de cierre F2.
  - Depende de: E1–E11 y puertas de calidad verdes.
  - Procede de: los 25 REQ de cierre asignados a F2.5 en la matriz de §13 y evidencia compartida de REQ-MAN-83; ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100 y ACE-103.
  - Cambia: matriz de evidencia F2.5, informe de validación independiente, inventario de limitaciones y dossier de cierre; separa la evidencia compartida y no marca requisitos de fases posteriores como completos.
  - Se comprueba con: §4.3, recuento automático de 25 REQ de cierre y un REQ compartido sin huérfanos, comandos registrados de todos los validadores y revisión cruzada de los 58 REQ de cierre más los cuatro REQ compartidos de §13.

> **Pausa obligatoria F2.5 y cierre de F2.** No autorizar F3 hasta que una
> revisión independiente reconstruya todos los haplotipos, confirme
> coordenadas/alelos/fase/dosis, valide VCF/BGZF/índice y demuestre los cuatro
> casos de `N`: SNP rechazado, inserción aleatoria sin `N`, inserción explícita
> con `N` conservada y VCF/complemento sin resolución inventada. Aprobar por
> separado el cierre de F2; completar E12 no lo concede automáticamente.

## 10. Riesgos genómicos silenciosos y oráculos obligatorios

| Riesgo silencioso | Prevención en producción | Oráculo independiente obligatorio | Evidencia de fallo que debe detectarse |
|---|---|---|---|
| Desfase base 0/base 1 o extremo abierto/cerrado | Tipos y único servicio de A5 | Vectores literales sobre cadenas etiquetadas por posición | SNP desplazado, indel con ancla incorrecta o acceso a posición 0 |
| Longitud de *indel* mal resuelta o convertida en distribución implícita | Parámetros exactos separados B1–B4 y constructores B10–B11 | Casos 1/2/`umbral - 1` y RNG instrumentado | Longitud distinta, consumo RNG para longitud o clasificación estructural incorrecta |
| Inserción/deleción terminal inválida | Constructores validados B10–B11 | Enumeración manual de todos los indels de una cadena corta | Pérdida del ancla o secuencia sin base a uno de los lados |
| Sustitución silenciosa de `N` | Dominio A6 y reglas B8–B10/E7 | Corpus con `N` en borde/interior y barrido de salidas | SNP `REF=N`, ALT `N` o inserción aleatoria con `N` |
| Fase o dosis alterada | Copias ordenadas A7/A9 y GT E5 | Conteo manual por copia, no parsing circular del VCF | GT sin P alelos, separador no faseado o dosis distinta |
| Evento perdido por solapamiento | Composición B12 e historial D2 | Aplicador de cadenas corto escrito desde fixtures | SNP superviviente tras ser borrado o conflicto aceptado |
| Tasa cero que consume RNG | Ruta explícita B6 | RNG instrumentado y comparación de ámbitos vecinos | Cambios en otras clases o individuos al activar cero |
| Poisson incorrecto o confundido con otra distribución | Decisión B1 y servicio B7 | Vectores o implementación ajena y test estadístico prefijado | Media/varianza fuera del umbral o resultado no reproducible |
| Verdad y secuencia materializada divergen | Una verdad canónica A8–A9/C6 | Segundo reconstructor de D4 | Bytes o longitud distintos para la misma copia |
| VCF convertido en fuente de verdad | Dirección de dependencias D2/E3 | Prueba sin VCF y auditoría de módulos | No se puede reconstruir o cambia la verdad al desactivar VCF |
| REF/ALT/normalización incorrectos | Proyección E3–E4 | Normalizador/validador independiente aprobado | Registro válido sintácticamente pero no representa la verdad |
| Pérdida de detalle al agrupar multialélicos | Agrupación E5 y complemento E7 | Tabla manual sitio×copia y reconstructor combinado | ALT o copia portadora desaparece, o se inventa una base |
| Artefactos mezclados entre ejecuciones | IDs/hashes/publicación F1 reutilizados | `sha256sum`, referencias cruzadas y fallos inyectados | Índice, VCF, verdad o referencia de ejecuciones distintas |

Una prueba que serializa y vuelve a leer con el mismo código solo demuestra
consistencia interna. No satisface por sí sola ninguna fila de esta tabla.

## 11. Cambios pequeños, revisión y política de bloqueo

- Cada paso se implementa y revisa de forma aislada; una revisión puede dividirlo
  aún más, pero no fusionarlo con un paso dependiente para ocultar su evidencia.
- Si A2, A3, B1 o E1 no obtienen aprobación, se detienen sus consumidores y se
  actualiza este plan con el bloqueo. Se puede continuar únicamente con corpus,
  matrices o validadores que no presupongan la opción elegida.
- Si una decisión aprobada contradice especificación, arquitectura o DTC
  confirmada, no se adapta el requisito al código: se pausa y se eleva la
  contradicción para decisión explícita.
- Las mediciones pequeñas observan regresiones y consumo, pero no prometen los
  objetivos de escalabilidad de F9 ni justifican optimizaciones especulativas.
- Ninguna pausa se supera solo porque compila. Requiere revisión de contratos,
  pruebas negativas, oráculo independiente y matriz REQ/ACE actualizada.

## 12. Validación mecánica del documento

Mientras el fichero sea nuevo, `git diff --no-index` termina con código 1 por
existir contenido; es satisfactorio si `--check` no imprime diagnósticos. Antes
de aprobar el plan se ejecutarán, como mínimo:

```bash
git diff --no-index --check /dev/null design/plan-implementacion-f2-0.1.md
rg -n '^- ☐ [A-E][0-9]+\.' design/plan-implementacion-f2-0.1.md
rg -n '^> \*\*Pausa obligatoria F2\.[1-5]' design/plan-implementacion-f2-0.1.md
rg -n 'DTPD-101|DTPD-102|DTPD-103|F1\.4|F1\.5' design/plan-implementacion-f2-0.1.md
```

Además, un comprobador mecánico debe verificar que cada casilla contiene una y
solo una línea `Depende de`, `Procede de`, `Cambia` y `Se comprueba con`; que la
matriz siguiente contiene una sola vez los 58 REQ cuyo cierre está asignado
exclusivamente a F2; y que la lista separada de evidencia compartida contiene
los otros cuatro REQ cuyo campo de hito incluye F2 junto con fases posteriores.
La unión debe coincidir exactamente con las 62 filas de la trazabilidad cuyo
hito contiene F2, sin confundir evidencia parcial con cierre total.

## 13. Matriz de trazabilidad REQ/ACE/DTC

La lista de REQ de cierre exclusivo es exhaustiva para la asignación F2 vigente:
4 + 22 + 5 + 2 + 25 = **58**. Se añaden por separado **cuatro REQ compartidos**
que exigen evidencia parcial en F2 y cierre posterior. Los ACE son la unión de
los criterios enlazados por esas filas de trazabilidad. Una mención no acredita
ejecución ni cierre. REQ-CONF-08, REQ-CONF-21, REQ-IND-02 y REQ-IND-03 siguen
señalados como cierres exclusivos solo porque así consta en la matriz maestra;
RAF2-004 impide marcarlos satisfechos hasta revisar esa asignación frente al
alcance deliberado de F2.

| Hito | REQ asignados | ACE cubiertos | Decisiones/puertas principales | Pasos de cierre |
|---|---|---|---|---|
| F2.1 | REQ-REF-14; REQ-IND-28; REQ-IND-29; REQ-VER-27 | ACE-03, ACE-10, ACE-11, ACE-28, ACE-66 | DTC-003, DTC-010, DTC-013, DTC-031, DTC-032, DTC-039; DTPD-101 y DTPD-102 | A2–A10 |
| F2.2 | REQ-CONF-08; REQ-CONF-21; REQ-CONF-24; REQ-CONF-25; REQ-CONF-28; REQ-IND-02; REQ-IND-03; REQ-IND-05; REQ-IND-06; REQ-IND-10; REQ-IND-11; REQ-IND-13; REQ-IND-14; REQ-IND-17; REQ-IND-18; REQ-IND-19; REQ-IND-21; REQ-IND-38; REQ-IND-40; REQ-IND-41; REQ-DIR-03; REQ-REC-24 | ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109, ACE-111 | DTC-006, DTC-011, DTC-020, DTC-022, DTC-032, DTC-033, DTC-037, DTC-038, DTC-044, DTC-049; DTPD-101 y decisión B1 aún sin identificador | B1–B14 |
| F2.3 | REQ-CONF-20; REQ-CONF-22; REQ-IND-15; REQ-DIR-01; REQ-DIR-02 | ACE-01, ACE-10, ACE-12, ACE-111 | DTC-003, DTC-006, DTC-011, DTC-020, DTC-037, DTC-038, DTC-039 | C1–C8 |
| F2.4 | REQ-VER-31; REQ-SAL-06 | ACE-34, ACE-62 | DTC-003, DTC-004, DTC-031, DTC-035, DTC-036, DTC-040, DTC-045, DTC-047; DTPD-102 | D1–D8 |
| F2.5 | REQ-FMT-05; REQ-FMT-06; REQ-REC-26; REQ-VER-02; REQ-VER-06; REQ-VER-12; REQ-VER-15; REQ-VER-16; REQ-VER-17; REQ-VER-30; REQ-VER-32; REQ-VER-33; REQ-VER-34; REQ-VER-36; REQ-VER-37; REQ-VER-39; REQ-VER-40; REQ-VER-41; REQ-VER-43; REQ-VER-44; REQ-VER-47; REQ-VAL-28; REQ-MAN-27; REQ-MAN-57; REQ-SAL-21 | ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100, ACE-103 | DTC-003, DTC-004, DTC-013, DTC-021, DTC-032, DTC-040, DTC-043, DTC-045, DTC-047; DTPD-103 | E1–E12 |

### 13.1. Evidencia compartida obligatoria en F2

Estas cuatro filas tienen un destino compuesto en la matriz maestra. F2 debe
producir la evidencia concreta indicada, pero no declarará el REQ completamente
satisfecho hasta que cierren también sus fases posteriores:

| REQ compartido | Evidencia parcial exigida en F2 | Pasos |
|---|---|---|
| REQ-IND-25 (F2.3/F6.1) | Un creador directo real asigna un ID estable cuando falta, comprueba unicidad y no reutilización dentro del recorrido disponible | C2, C8 |
| REQ-IND-26 (F2.4/F3.5/F7.7/F8.5) | El mismo ID aparece sin traducción en estado genotípico, verdad estructurada y metadatos F2; observaciones y linaje quedan para sus fases | D2, D6, D8 |
| REQ-MAN-10 (F2.3/F5.5/F9.7) | El manifiesto distingue valores globales, específicos y efectivos consumidos por individuos/regiones/clases de F2 | C2, C8 |
| REQ-MAN-83 (F2.5/F3.5) | El primer conjunto real y sus artefactos VCF/índice manifiestan estados y diagnósticos coherentes con ejecución y evento | E9, E12 |

### 13.2. Evidencia adyacente no asignada a F2

Los casos de aceptación de F2 deben conservar evidencia útil sobre `N`,
inserciones aleatorias, VCF y salida directa, pero la matriz no incorporará como
cerrados REQ-IND-58, REQ-IND-59, REQ-VER-01, REQ-VER-24 a REQ-VER-26,
REQ-FMT-07 ni REQ-SAL-02 porque la trazabilidad vigente los asigna a hitos
posteriores con alcance mayor. Cualquier discrepancia entre esa asignación y la
evidencia de cierre pedida por la hoja de ruta se documentará como cobertura
parcial, no se resolverá alterando silenciosamente la fase.

## 14. Criterio de salida del plan

El plan queda listo para aprobación cuando su validación mecánica no detecta
errores, los 58 REQ de cierre aparecen exactamente una vez en §13, los cuatro
REQ compartidos aparecen exactamente una vez en §13.1, cada hito tiene su pausa
y cada paso conserva las cuatro líneas obligatorias. La ejecución seguirá
condicionada al cierre de F1 y a las decisiones A2, A3, B1 y E1. F2 solo podrá
cerrarse cuando E12, la pausa final y la aprobación independiente de cierre
consten por separado.
