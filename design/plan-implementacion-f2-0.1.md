# Plan de implementación de F2 0.1

## 1. Cabecera y estado

- **Título:** plan detallado de F2, verdad genómica mínima.
- **Versión del plan:** 0.1.
- **Fecha:** 13 de agosto de 2026.
- **Estado:** borrador técnico condicionado; F1.4 está aprobada, pero no está
  autorizado para ejecución hasta que F1.5, F1.6 y la revisión de cierre de F1
  estén terminadas y este plan sea aprobado explícitamente.
- **Ámbito temporal:** únicamente F2.0–F2.5.
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
  DTC-049 y DTC-051 a DTC-057.
- **Decisiones pendientes bloqueantes:** DTPD-102 bloquea
  F2.1/F2.4 y DTPD-103, F2.5. También falta documentar y aprobar la
  concreción técnica de Poisson indicada en §4.2. Para los *indels* aleatorios
  de F2 ya se acordó una longitud exacta configurable por separado para
  inserciones y deleciones, con valor predeterminado de 1 pb; falta trasladar
  este acuerdo al registro canónico de decisiones.
  DTPD-117 bloquea F2.0 hasta evaluar y fijar Marimo, su entorno, empaquetado y
  contrato versionado con el binario. DTC-054 resuelve DTPD-118 y fija para
  F2.2/F2.3 los formatos de matriz/espectro, el algoritmo reproducible, los
  perfiles versionados y las excepciones limitadas. DTC-055 resuelve DTPD-119 y
  deja preparada la inicialización fundadora desde paneles para F6.4; no se
  consume dentro de F2. DTC-055 ya cierra para F6.4 el constructor común y su
  variante uniforme; F7.1 lo reutilizará después para la recombinación meiótica
  y DTPD-111 solo resolverá mapas y semántica meiótica.
- **Pendientes del informe adversarial que siguen abiertos:** RAF2-001,
  RAF2-004 y RAF2-007. Esta ampliación no los da por resueltos; DTC-052 aporta
  la semántica que faltaba para reevaluar RAF2-004, pero su cierre requiere
  evidencia y revisión formal de la asignación maestra.

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
  sobre toda una secuencia o regiones declaradas. Los SNV heredables iniciales
  proceden de una matriz común a la población, suministrada o generada desde
  número/densidad, espacio, espectro de sustituciones y espectro alélico; no se
  generan por densidades independientes por individuo. Solo un recuento pedido
  como media o densidad usa Poisson: ni frecuencias, alelos ni genotipos se
  muestrean mediante Poisson. La longitud de los *indels* aleatorios de F2 es un
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
- Aplicación Marimo base con nueve pestañas, cobertura completa de F1/F1.6,
  importación/exportación del conjunto TOML e invocación segura del binario; durante F2 se
  amplía para individuos, verdad estructurada y VCF.

### 2.2. Excluido

- Lecturas, FASTQ, BAM/CRAM, mapeo y cualquier error técnico; comienzan en F3.
- Variantes estructurales, repeticiones, referencias derivadas, ploidía distinta
  por secuencia y la cobertura VCF/complemento completa de casos complejos; F5.
- Dinámica de poblaciones, evolución, recombinación, selección, familias y
  cruces; F6–F8. F2 solo usa la población de origen necesaria para interpretar
  la matriz de frecuencias y crear individuos directos, sin simularla.
- Paneles de haplotipos, desequilibrio de ligamiento y creación fundadora
  intacta o recombinante desde ellos; DTC-055 destina ambas alternativas a
  F6.4. F7.1 reutiliza el motor común en la reproducción, pero no inicializa el
  conjunto fundador.
- Distribuciones no Poisson para recuentos, distribuciones probabilísticas de
  longitud de *indels*, gVCF con posiciones invariantes y cualquier resolución
  silenciosa de `N`. Las distribuciones configurables de longitud corresponden
  a REQ-MUT-08 en F7.2; F2 solo admite la longitud exacta acordada.
- DTF-003, DTF-008, DTF-013, DTF-014, DTF-023, DTF-033, DTF-034, DTF-036 y
  DTF-037,
  además de las demás ampliaciones post-v1 del registro vigente.
- Afirmar el cierre total de requisitos asignados a fases posteriores. F2 puede
  dejar evidencia parcial reutilizable para, por ejemplo, REQ-IND-58,
  REQ-IND-59, REQ-VER-01, REQ-VER-24 a REQ-VER-26, REQ-FMT-07 y REQ-SAL-02,
  pero no los marcará satisfechos en su matriz de cierre.
- Reimplementar en Marimo validación, dominancia, biología, generación aleatoria
  o publicación; la aplicación solo configura e invoca al binario.

## 3. Orden, reutilización y fronteras

1. Cerrar F1.6 y F1, implementando DTC-056 y el servicio común de coordenadas, y
   congelar su línea base antes de modificar contratos de dominio,
   configuración, referencia, aleatoriedad o publicación.
2. Resolver DTPD-117 y cerrar F2.0 antes de F2.1. Después resolver DTPD-102 con
   prototipos y oráculos independientes antes de implementar estado haplotípico;
   resolver DTPD-103 antes de implementar VCF.
3. Reutilizar y extender el único servicio de coordenadas decidido e iniciado
   en F1.6. Las fronteras traducen la base 1 cerrada configurada; ningún
   adaptador ni algoritmo realiza conversiones implícitas ni abre otra
   representación interna.
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

Antes de M1 deben constar: F1.4, F1.5 y F1.6 terminadas; revisión de cierre de F1
aprobada; pruebas y evidencias F1 verdes; hashes/esquemas/catálogos F1 congelados;
árbol de trabajo atribuido; y aprobación explícita de este plan. Si F1 cambia
un contrato que F2 pretendía reutilizar, se actualiza primero este plan y su
trazabilidad.

### 4.2. Decisiones que detienen solo a sus consumidores

- DTPD-101 está resuelta mediante DTC-056. F1.6 implementa representación común,
  tipos de frontera y vectores de intervalos origen↔derivada; F2.1/F2.2 deben
  auditarlos y extender el mismo servicio para segmentos, anclas y composición
  sin reabrir la decisión.
- DTPD-102 debe elegir por separado el estado haplotípico interno y el fichero
  obligatorio de verdad, demostrando reconstrucción, memoria, fase, dosis y
  procedencia.
- DTPD-103 debe fijar versión VCF, normalización, campos `INFO`, representación
  complementaria y validadores independientes, incluidos BGZF e índice.
- DTPD-117 debe evaluar y fijar versión compatible de Marimo, entorno,
  dependencias, licencia, seguridad, empaquetado Ubuntu/WSL2, contrato de
  esquemas/catálogos y adaptador de proceso sin shell. Ninguna versión se
  presupone en este plan.
- DTC-054 resuelve DTPD-118 y debe consumirse sin reinterpretación: solo admite
  `alt_af`; separa el espectro, la matriz larga canónica y sus metadatos; fija
  generación bialélica, poblaciones independientes, selección sin reemplazo,
  perfiles explícitos, ámbitos RNG, preparación/reutilización y excepciones
  limitadas. B1 debe congelar su versión documental antes de implementar sus
  consumidores. La capacidad de panel resuelta por DTC-055 permanece fuera de
  este alcance.
- DTC-052 fija la separación funcional posterior a esa puerta: F2.2 produce o
  admite la matriz poblacional y F2.3 muestrea genotipos HWE o `F` diploide.
  Poisson solo puede obtener un recuento de loci/eventos desde una media o
  densidad; no gobierna frecuencias, alelos, genotipos ni fase.
- La documentación no concreta todavía la implementación Poisson. B2 debe
  crear una entrada de decisión conforme a §9 de
  `decisiones-tecnicas-0.1.md`, evaluar opciones, fijar el algoritmo y aprobarla.
  Esta puerta debe cerrarse antes de sus consumidores de recuentos.
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

## 5. F2.0 — Aplicación Marimo base

**Resultado del hito:** una aplicación local versionada permite configurar,
validar y lanzar todas las capacidades cerradas en F1/F1.6 mediante el mismo
TOML y binario, sin duplicar semántica y con las nueve pestañas preparadas para
las ampliaciones F2–F9.

- ☐ M1. Resolver DTPD-117 y aprobar la dependencia y distribución de Marimo.
  - Depende de: cierre explícito de F1.6/F1 y aprobación de este plan.
  - Procede de: DTPD-117, DTC-043, DTC-048 y DTC-051; REQ-UI-01, REQ-UI-12 y REQ-UI-13.
  - Cambia: evaluación con versión exacta, entorno, lock/pinning, licencia, mantenimiento, seguridad, empaquetado local Ubuntu/WSL2 y protocolo de proceso; solo después incorpora los archivos/dependencias aprobados de la app.
  - Se comprueba con: revisión DTC-043 equivalente para el entorno, instalación reproducible aislada, inventario bloqueado, arranque local y aprobación fechada de DTPD-117 sin introducir Python en el núcleo Rust.

- ☐ M2. Publicar el contrato versionado aplicación↔esquema↔catálogos↔binario.
  - Depende de: M1 y esquemas/catálogos cerrados en F1.6.
  - Procede de: REQ-UI-01, REQ-UI-09 y REQ-UI-12; DTC-036 y DTC-051.
  - Cambia: metadatos o comando de capacidades aprobado que declara versiones compatibles de esquema TOML general/escenario, campos, diagnósticos y protocolo `validate`/`run`, sin exponer tipos internos.
  - Se comprueba con: conjuntos compatible, anterior, posterior, ausente y corrupto; toda incompatibilidad de cualquier documento bloquea importación/validación/ejecución antes de interpretar parcialmente opciones.

- ☐ M3. Crear el armazón accesible de las nueve pestañas y su visibilidad por vía/capacidad.
  - Depende de: M1–M2.
  - Procede de: REQ-UI-02 y REQ-UI-13; DTC-051.
  - Cambia: aplicación con `General/Ejecución`, `Genoma`, `Individuos`, `Poblaciones`, `Evolución`, `Familias`, `Eventos técnicos`, `Salidas` y `Revisar y ejecutar`; las no implementadas quedan deshabilitadas/ocultas con causa explícita.
  - Se comprueba con: prueba de presencia y orden estable, rutas direct/evolution/families, navegación por teclado y comprobación de que una pestaña no disponible no pierde su destino trazable.

- ☐ M4. Mantener un único estado de formulario tipado e invalidar toda validación tras una edición.
  - Depende de: M2–M3.
  - Procede de: REQ-UI-03, REQ-UI-07 y REQ-UI-09; arquitectura §§6.1 y 7.1.
  - Cambia: modelo de autoría multidocumento separado de `EffectiveConfig`, con un general y cero o varios escenarios identificados/ordenados/enlazados; digest/revisión del conjunto validado y transición explícita editado→no validado→validado→ejecutable.
  - Se comprueba con: editar, añadir, quitar, reordenar o renombrar cualquier documento/enlace después de validar, verificar invalidación inmediata y demostrar que la app nunca calcula regla ganadora, valor efectivo o verdad biológica.

- ☐ M5. Cubrir con controles tipados todas las opciones escalares F1/F1.6 de General y Genoma.
  - Depende de: M2–M4.
  - Procede de: REQ-UI-01, REQ-UI-03 y REQ-CONF-124 a REQ-CONF-127; DTC-050/DTC-051.
  - Cambia: controles de clase/versión, vía, seed, referencia simulada/suministrada/reducida, IDs, clases, tamaños, `region_size_bases`, unión de objetivo y opciones de salida existentes, con unidad normativa en bases y ayuda en kb.
  - Se comprueba con: matriz opción F1/F1.6→pestaña/control/campo TOML, límites y exclusividad de reducción, conversión visual kb↔bases sin pérdida y ausencia de opción implementada solo en TOML.

- ☐ M6. Implementar tablas editables para secuencias, eventos, reglas y demás listas F1.
  - Depende de: M3–M5.
  - Procede de: REQ-UI-04, REQ-UI-05 y DTC-051.
  - Cambia: tablas con añadir/eliminar/reordenar, IDs y orden visibles y editores de selectores/ámbitos permitidos según catálogo; no concede precedencia por orden accidental de celdas.
  - Se comprueba con: alta/baja/reordenación, duplicados, región/locus, selector no admitido y comparación del TOML resultante con fixtures manuales equivalentes.

- ☐ M7. Importar y exportar un conjunto TOML autoconsistente.
  - Depende de: M2 y M4–M6.
  - Procede de: REQ-UI-06 y REQ-UI-12; DTC-034, DTC-036 y DTC-051.
  - Cambia: adaptador de importación al estado de formulario multidocumento, previsualización individual y exportación/guardado conjunto de exactamente un general y cero o varios escenarios, conservando IDs, orden y enlaces/rutas relativas. F2.0 activa el caso general F1/F1.6 y prepara escenarios sin implementar sus controles biológicos.
  - Se comprueba con: general sin escenarios y fixtures estructurales con varios escenarios, rutas relativas y nombres con espacios; faltante, duplicado, ajeno y enlace roto se rechazan; el ciclo importar→previsualizar→exportar→`validate` produce la misma `EffectiveConfig` o diagnóstico que el conjunto manual.

- ☐ M8. Invocar `validate` sin shell y asociar diagnósticos a campos sin reinterpretarlos.
  - Depende de: M1–M2 y M7.
  - Procede de: REQ-UI-07, REQ-UI-09 y REQ-VAL-54; DTC-044 y DTC-051.
  - Cambia: adaptador de proceso con argumentos separados para el general y sus escenarios enlazados, captura acotada de stdout/stderr/código y presentación por documento/pestaña/tabla/fila/campo usando el diagnóstico JSON del binario.
  - Se comprueba con: proceso falso que registra cada argumento, rutas relativas con espacios/metacaracteres, binario real válido/inválido, códigos 0–5/70 y prueba de que ningún texto se ejecuta como shell.

- ☐ M9. Seleccionar el destino de forma segura y lanzar `run` solo por acción explícita.
  - Depende de: M4, M8 y publicación F1.4.
  - Procede de: REQ-UI-08 a REQ-UI-10; DTC-045 y DTC-051.
  - Cambia: selector/revisión de destino, confirmación visible de no-clobber y botón de ejecución habilitado solo para el digest validado vigente; abrir/importar/editar nunca ejecuta.
  - Se comprueba con: destino vacío/existente/válido, doble clic, edición posterior, recarga de app y proceso falso; solo una acción explícita válida lanza una lista de argumentos y el binario conserva el rechazo no-clobber.

- ☐ M10. Presentar estado, logs, diagnósticos y resultados sin cargar artefactos grandes.
  - Depende de: M8–M9 y manifiestos F1.4/F1.6.
  - Procede de: REQ-UI-11 y DTC-051.
  - Cambia: vista de proceso/estado/código, logs acotados o paginados y resumen desde manifiestos con enlaces a artefactos; nunca lee completos FASTA/FASTQ/BAM/CRAM/VCF para resumir.
  - Se comprueba con: proceso largo/fallido, log grande, manifiesto válido/roto y artefacto señuelo grande instrumentado para demostrar que solo se leen metadatos o fragmentos solicitados.

- ☐ M11. Cerrar accesibilidad y funcionamiento local Ubuntu/WSL2.
  - Depende de: M3–M10.
  - Procede de: REQ-UI-13; DTC-042, DTC-051 y DTPD-117 aprobada.
  - Cambia: etiquetas, ayuda/unidades, foco, teclado, estados textuales además de color y guía de arranque local sobre sistema de ficheros Linux; no promete Windows nativo.
  - Se comprueba con: recorrido completo solo teclado, inspección de etiquetas/foco/contraste acordada, Ubuntu 24.04 y WSL2/ext4, y rechazo/documentación de `/mnt/c` conforme al entorno admitido.

- ☐ M12. Completar equivalencia y trazabilidad de F2.0.
  - Depende de: M1–M11.
  - Procede de: REQ-UI-01 a REQ-UI-13, REQ-VAL-54 y ACE-114/115.
  - Cambia: 14 casos `TEST-REQ-*`, matriz opción→control, pruebas app→conjunto TOML→binario y expediente DTPD-117; no marca controles de escenarios evolutivos/familias ni otras capacidades futuras como implementados.
  - Se comprueba con: par UI/manual para cada capacidad F1/F1.6 y para la estructura general+0..N escenarios que produce igual `EffectiveConfig` o diagnóstico, puertas propias de la app y recuento automático de 14 REQ sin huérfanos.

> **Pausa obligatoria F2.0.** No iniciar F2.1 hasta aprobar DTPD-117, la
> equivalencia de todas las opciones F1/F1.6, la invocación sin shell, la
> autoría/exportación autoconsistente del conjunto TOML, la invalidación tras
> cambiar cualquier documento/enlace, el destino no-clobber y la revisión de
> accesibilidad. RAF2-001, RAF2-004 y RAF2-007 continúan abiertos.

## 6. F2.1 — Estado haplotípico mínimo

**Resultado del hito:** coordenadas y estado propios representan cada copia
haplotípica, fase, dosis, procedencia y `N`, y permiten materializar secuencias
exactas desde la referencia cerrada en F1.

- ☐ A1. Congelar y auditar la línea base que F2 consumirá.
  - Depende de: cierre explícito de F1.4, F1.5 y F1, y aprobación de este plan.
  - Procede de: dependencia F1→F2 de la hoja de ruta; DTC-043, DTC-045 y DTC-048.
  - Cambia: evidencia de preflight de F2 con commit base, toolchain, lock, esquemas, catálogos, comandos verdes y contratos F1 reutilizados; no cambia producción.
  - Se comprueba con: revisión del acta de cierre F1, `git status --short`, `rustc --version`, `cargo metadata --locked --format-version 1` y las puertas de §4.3 sobre el commit base.

- ☐ A2. Auditar y congelar para F2 la DTPD-101 resuelta en F1.6.
  - Depende de: A1.
  - Procede de: DTPD-101; DTC-010; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: expediente de consumo F2 que enlaza la decisión y los tipos/servicio comunes de F1.6, completa vectores de inserciones ancladas, deleciones internas y composición y prohíbe otra representación; no reabre opciones ni cambia DTPD-101.
  - Se comprueba con: huella/versión coincidentes con F1.6, vectores literales de primera/última base y rangos vacíos, ida/vuelta en fronteras y revisión independiente que detectaría cualquier conversión o tipo local nuevo.

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

- ☐ A5. Extender el servicio común de coordenadas de F1.6 a los casos de F2.
  - Depende de: A2 y A4 cuando corresponda.
  - Procede de: DTC-010, DTC-013 y DTC-039; REQ-CONF-28, REQ-IND-38 y REQ-IND-39.
  - Cambia: los tipos y el servicio ya usados por el mapa de reducción ganan constructores validados para posición, intervalo, ancla de inserción y deleción interna; no expone tipos de crates ni crea un segundo módulo/convenio.
  - Se comprueba con: pruebas de tabla compartidas con F1.6 desde la semántica externa base 1 cerrada, propiedades de ida/vuelta y fallos explícitos en cero, fuera de rango, desbordamiento y límites terminales.

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

> **Pausa obligatoria F2.1.** Revisar el consumo sin reapertura de DTPD-101 y
> aprobar DTPD-102, la extensión de tipos de coordenada, representación
> haplotípica, consumo de memoria, oráculos y los
> cuatro REQ. No iniciar B27 ni F2.3/F2.4 si la representación no reconstruye
> exactamente cada copia o si una conversión de coordenadas queda implícita.

## 7. F2.2 — Variantes pequeñas y matriz de frecuencias SNV

**Resultado del hito:** el motor admite o genera una matriz poblacional de
frecuencias SNV reproducible y validada, aplica SNP explícitos e *indels*
pequeños inequívocos y nunca usa Poisson para muestrear frecuencias, alelos o
genotipos.

- ☐ B1. Congelar como entrada de implementación DTC-054, que resuelve DTPD-118.
  - Depende de: pausa F2.1 aprobada, A2 y A3.
  - Procede de: DTC-052, DTC-054 y DTPD-118 resuelta; REQ-CONF-128, REQ-CONF-129, REQ-IND-60 y REQ-IND-61.
  - Cambia: expediente de consumo que enlaza la versión aprobada y enumera formato, algoritmos, perfiles, límites, ámbitos aleatorios y excepciones; declara fuera de F2 VCF→TSV, generación multialélica, estructura/evolución poblacional y la capacidad de panel de DTC-055.
  - Se comprueba con: revisión requisito→cláusula DTC-054→paso que no deja comportamiento aprobado sin consumidor ni reabre DTPD-118 al implementar.

- ☐ B2. Registrar y aprobar la concreción Poisson y las longitudes exactas de *indels* de F2.
  - Depende de: pausa F2.1 aprobada y A2–A3; es independiente de la elección de formato B1.
  - Procede de: REQ-CONF-23 a REQ-CONF-25 y REQ-IND-16 a REQ-IND-18; DTC-022, DTC-037, DTC-043 y DTC-049; §4.2.
  - Cambia: decisión técnica que fija implementación Poisson, límites, espacio vacío y consumo RNG, selección espacial y longitudes exactas por clase con predeterminado 1 pb; declara expresamente que Poisson solo produce recuentos desde medias/densidades.
  - Se comprueba con: oráculo Poisson independiente, prueba estadística con umbral previo, casos cero/vacío y longitudes 1/2/`umbral - 1`, rechazo de 0 o `>= umbral` y aprobación formal.

- ☐ B3. Publicar los esquemas versionados de espectro y matriz SNV.
  - Depende de: B1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-CONF-129 y REQ-IND-60.
  - Cambia: esquema del espectro TSV por intervalos y filas `fixed`; entrada suministrada mínima `contig/position/ref/alt/population_id/alt_af` con `source_id` opcional; pareja canónica `population_snvs.tsv` + `population_snvs.metadata.json`, sin entrada VCF.
  - Se comprueba con: fixtures manuales válidos e inválidos validados por una implementación de esquema separada, rechazo de columnas de frecuencia distintas de `alt_af` y prueba negativa de que `run` no acepta VCF ni TSV sin metadatos canónicos.

- ☐ B4. Fijar el contrato de metadatos, integridad y referencia de la pareja canónica.
  - Depende de: B3 y contratos SHA-256 de F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: JSON con versión, referencia esperada, SHA-256 del FASTA, coordenadas `1-based-closed`, poblaciones, procedencia sintética/empírica y SHA-256 del TSV; muestra/fuente/filtros desconocidos se declaran `unknown`, y configuración/manifiesto enlazan el SHA-256 del JSON sin digest circular.
  - Se comprueba con: pareja escrita a mano, `sha256sum`, intercambio independiente de TSV/JSON/FASTA, contig literal `1` frente a `chr1` y rechazo previo a individuos ante cualquier discordancia.

- ☐ B5. Especificar la autoría TOML excluyente de matriz preparada, reutilizada o generada y de los demás eventos F2.
  - Depende de: B2–B4.
  - Procede de: REQ-CONF-08, REQ-CONF-21, REQ-CONF-24, REQ-CONF-25, REQ-CONF-28, REQ-CONF-128 y REQ-CONF-129; DTC-034 y DTC-052.
  - Cambia: esquema documental para fuente canónica o generación, referencia/poblaciones, `exact_loci`/`mean_loci`/`density_per_mb`, espacio, sustituciones y espectro explícito, más SNP e *indels*; elimina densidad SNP inicial por individuo.
  - Se comprueba con: ejemplos de las tres rutas y rechazo de fuentes simultáneas, perfil alélico ausente, más de un modo de recuento, densidad SNP individual, VCF y campos futuros.

- ☐ B6. Extender configuración declarada y parser conservando literales y ubicaciones.
  - Depende de: B3–B5.
  - Procede de: DTC-034, DTC-036, DTC-044 y DTC-052; REQ-CONF-128, REQ-CONF-129 y REQ-DIR-03.
  - Cambia: `DeclaredConfig`, adaptador TOML y diagnósticos para matriz, metadatos, espectros, regiones y eventos pequeños, sin construir verdad, calcular digests ni normalizar silenciosamente.
  - Se comprueba con: parsing de cada forma, claves desconocidas/duplicadas, clase/versión/ruta errónea y snapshots de documento, campo, literal y ubicación.

- ☐ B7. Extender catálogo y resolución efectiva con parámetros de matriz e *indels*.
  - Depende de: B5–B6.
  - Procede de: DTC-006, DTC-033, DTC-049 y DTC-052; REQ-IND-10, REQ-IND-13, REQ-IND-14 y REQ-CONF-129.
  - Cambia: `ParameterKey`, metadatos, dominios y solicitudes para recuentos excluyentes, espacio, sustituciones, perfil, modo `exact/stochastic` y longitudes exactas de inserción/deleción, sin duplicar el resolvedor.
  - Se comprueba con: catálogo completo, dominios/unidades, cero explícito, precedencia espacial, valor efectivo y rechazo de ejes no admitidos o fallback empírico implícito.

- ☐ B8. Versionar los ámbitos aleatorios independientes aprobados para la matriz.
  - Depende de: B1, B7 y catálogo RNG cerrado en F1.3.
  - Procede de: DTC-011, DTC-037, DTC-038 y DTC-052; REQ-IND-61 y REQ-VAL-55.
  - Cambia: ámbitos `matrix/locus_selection`, `matrix/alternate_allele/<variant_id>` y `matrix/frequency/<population_id>/<variant_id>`, separados de *indels*, genotipo y fase, con versión y vectores BLAKE3.
  - Se comprueba con: vectores literales, claves únicas y pruebas de que añadir/reordenar una población no cambia loci, ALT ni frecuencias existentes y cambiar un perfil solo afecta a esa población.

- ☐ B9. Enumerar bases elegibles antes de determinar o seleccionar loci.
  - Depende de: A5–A6 y B7–B8.
  - Procede de: REQ-CONF-21, REQ-CONF-28, REQ-IND-13, REQ-IND-14, REQ-IND-19 y REQ-DIR-03.
  - Cambia: servicio que cuenta posiciones `A/C/G/T` elegibles tras inclusiones, exclusiones, posiciones protegidas/ocupadas y pesos cero; excluye `N` y no sortea todavía.
  - Se comprueba con: enumeración ajena al productor en secuencias cortas, contigs/regiones literales, espacio vacío y denominador de bases elegibles por cromosoma.

- ☐ B10. Implementar `uniform` y `weighted_regions` sin solapamientos.
  - Depende de: B3, B7 y B9.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-19 y REQ-VAL-55.
  - Cambia: perfil incorporado único `uniform` con peso 1; tabla propia `1-based-closed` con pesos finitos no negativos, `outside_weight` explícito, regiones no solapadas y probabilidad por base proporcional al peso.
  - Se comprueba con: tabla manual por longitud×peso, casos `outside_weight` 0/1, peso 0/0.5/1/2, frontera contigua, solapamiento rechazado y requisito de al menos una base con peso positivo.

- ☐ B11. Determinar el número de loci candidatos sin confundirlo con loci activos.
  - Depende de: B2 y B7–B10.
  - Procede de: REQ-CONF-129, REQ-IND-17, REQ-IND-18 y REQ-IND-61; DTC-022 y DTC-052.
  - Cambia: `exact_loci` fija candidatos; `mean_loci` y `density_per_mb × bases_elegibles/1_000_000` fijan una media y sortean Poisson; los tres modos son excluyentes y el recuento se conserva aunque después se eliminen ALT globalmente ausentes.
  - Se comprueba con: exactos, vectores Poisson para media/densidad, denominador manual, tasa cero, espacio vacío/insuficiente y prueba de que cambiar frecuencias no vuelve a calcular candidatos.

- ☐ B12. Validar espectros alélicos por intervalos no solapados.
  - Depende de: B3, B6–B7 y DTC-049.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-61 y REQ-VAL-55.
  - Cambia: lector `f64` de decimal/científica para `[min_af,max_af)`, último límite 1 y filas `fixed`; admite huecos, rechaza solapamiento, porcentajes, NaN/infinito y pesos inválidos; `log_uniform` exige mínimo positivo y solo `fixed` genera exactamente 0/1.
  - Se comprueba con: corpus de fronteras 0/1, intervalos contiguos/solapados/con huecos, pesos no normalizados, todo cero, `1e-3`, porcentaje, NaN/infinito y conservación del literal y valor efectivo sin redondeo.

- ☐ B13. Implementar los repartos `stochastic` y `exact` por población.
  - Depende de: B8, B11–B12.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-61 y REQ-VAL-55.
  - Cambia: sorteo ponderado de fila en modo estocástico predeterminado y mayores restos en modo exacto, con parte entera primero y desempate por orden estable del espectro; cada población aplica su reparto independientemente.
  - Se comprueba con: ejemplos manuales, incluido 11 loci con pesos 60/25/10/5→7/3/1/0, suma exacta por población, empate literal y prueba estadística prefijada del modo estocástico.

- ☐ B14. Incorporar solo los perfiles alélicos empíricos 1000 Genomes seleccionables de forma explícita.
  - Depende de: B1, B3, B12–B13 y evaluación DTC-043 de datos/licencia/herramienta.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-61 y ACE-116.
  - Cambia: recurso `1000g_phase3_v5_sfs_v1` global/AFR/AMR/EAS/EUR/SAS derivado de autosomas, PASS y SNV bialélicos; separa `fixed 0`, `min_positive–0.001`, 0.001–0.005, 0.005–0.01, 0.01–0.05, 0.05–0.10, 0.10–0.50, 0.50–0.90, 0.90–0.99, 0.99–1 y `fixed 1`, con recuentos, fuente, filtros, receta y SHA-256; nunca es predeterminado.
  - Se comprueba con: regeneración independiente desde la fuente versionada, comparación de recuentos/digests, selección nominal obligatoria y configuraciones sin perfil, con perfil propio y con perfil empírico explícito.

- ☐ B15. Resolver el espectro dirigido de las 12 sustituciones.
  - Depende de: B3 y B7.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129 y REQ-IND-61.
  - Cambia: `uniform_snv12_v1` predeterminado, tabla propia, atajo `ti_tv_ratio=R` convertido a peso `2R` para la transición y 1 para cada transversión, y perfil explícito `1000g_phase3_v5_snv12_v1`; la configuración efectiva conserva siempre los 12 pesos.
  - Se comprueba con: las tres ALT por cada REF, R=0/1/2, rechazo de R negativo/infinito, razón esperada independiente, perfil empírico reproducible y prueba de que este último nunca se carga por defecto.

- ☐ B16. Crear IDs canónicos y orden independiente del orden de entrada.
  - Depende de: B3–B4 y servicios de identidad F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60, REQ-IND-61 y REQ-VAL-55.
  - Cambia: `site_id` por referencia/contig/posición/REF y `variant_id` por referencia/contig/posición/REF/ALT; conserva `source_id` como alias y ordena por catálogo FASTA, posición, REF, ALT y población.
  - Se comprueba con: permutación de filas/poblaciones, mismo SNV en referencias distintas, sitio multialélico y comparación literal de IDs y bytes ordenados.

- ☐ B17. Preparar una matriz suministrada sin modificar el original.
  - Depende de: B3–B4, B9 y B16.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-IND-60, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: orden `prepare-snv-matrix --reference ... --input ... --output ...` que calcula SHA-256, IDs y orden, comprueba contigs/coordenadas/REF y publica de forma segura la pareja canónica nueva.
  - Se comprueba con: fixture mínimo, original inalterado, salida reconstruible, referencia intercambiada, fallo sin pareja parcial y equivalencia con la acción Marimo «Vincular y verificar referencia».

- ☐ B18. Validar filas completas y poblaciones de una matriz suministrada.
  - Depende de: B3–B4 y B16–B17.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-IND-60 y REQ-VAL-55.
  - Cambia: exige una fila por `ALT × población`, `alt_af` explícito en `[0,1]`, población declarada y claves únicas; rechaza vacío/NA, MAF/DAF, AC/AN, inferencias y alias de población, pero admite y valida poblaciones no seleccionadas.
  - Se comprueba con: corpus que omite una combinación, usa cero explícito, añade población válida no usada, solicita una inexistente e intenta cada semántica de frecuencia prohibida.

- ☐ B19. Validar loci multialélicos suministrados con suma estable y tolerancia aprobada.
  - Depende de: B18 y DTC-049.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60 y REQ-VAL-55.
  - Cambia: agrupa por sitio/población, calcula suma ALT con Neumaier en orden canónico y `ref_af=1-suma`; acepta suma ≤1, normaliza solo un exceso ≤`1e-9` y lo registra, y rechaza exceso mayor sin recortar valores individuales.
  - Se comprueba con: sumas manuales 0.5, 1, `1+5e-10` y 1.01, orden permutado y lector independiente que confirma frecuencias/normalización y diagnóstico.

- ☐ B20. Normalizar alelos globalmente ausentes sin reponer candidatos.
  - Depende de: B18–B19.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-60, REQ-IND-61, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: elimina cada ALT con `alt_af=0` en todas las poblaciones y el sitio si pierde todos sus ALT; conserva `alt_af=1`, cuenta candidatos/ALT descartados/loci activos por separado y nunca selecciona sustitutos.
  - Se comprueba con: bialélico global 0/1, multialélico con un ALT global 0 y otro activo, matriz suministrada/generada y recuentos antes/después calculados fuera del normalizador.

- ☐ B21. Implementar autorización explícita para una matriz activa vacía.
  - Depende de: B20 y diagnósticos F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-VAL-55 y REQ-MAN-91.
  - Cambia: `validate` avisa; `run` se detiene salvo `--allow-empty-snv-matrix`; la autorización no relaja otros errores y queda en configuración efectiva/manifiesto.
  - Se comprueba con: matriz vacía inicial y tras normalización, ejecución con/sin flag, SNP explícitos que aún se aplican y prueba negativa de que el flag no acepta digest, REF o frecuencia inválidos.

- ☐ B22. Seleccionar las posiciones candidatas sin reemplazo.
  - Depende de: B8–B11.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-03, REQ-IND-19 y REQ-IND-61.
  - Cambia: selector ponderado sobre las bases elegibles que obtiene exactamente el recuento candidato ya fijado, sin duplicados, y falla si no hay suficientes posiciones; no promete prefijo al cambiar N.
  - Se comprueba con: referencia exhaustiva `A/C/G/T/N`, uniformidad/ponderación con umbral previo, N=0/1/máximo/exceso, semilla fija y comparación con un enumerador independiente.

- ☐ B23. Asignar un único ALT generado por locus candidato.
  - Depende de: B8, B15–B16 y B22.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129 y REQ-IND-61.
  - Cambia: generador bialélico que lee REF del FASTA y sortea una de sus tres ALT mediante los 12 pesos; no genera multialelia ni contexto trinucleotídico.
  - Se comprueba con: cada REF A/C/G/T, pesos degenerados, perfiles uniforme/propio/empírico y prueba de que jamás produce REF, `N`, segundo ALT o cambio fuera del perfil.

- ☐ B24. Asignar frecuencias independientes a todas las poblaciones.
  - Depende de: B8, B12–B13, B16 y B23.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-61 y REQ-VAL-55.
  - Cambia: para cada `población × variant_id` elige fila y valor `alt_af` con su perfil, admite 0/1, conserva fila/subsemilla y no introduce población ancestral, FST, parentesco, migración, correlación o LD.
  - Se comprueba con: dos poblaciones y perfiles distintos, cambio aislado de una, frecuencias 0/intermedias/1, orden permutado y prueba de independencia sin exigir distribuciones correlacionadas.

- ☐ B25. Publicar la matriz generada mediante el mismo contrato canónico.
  - Depende de: B4, B16 y B20–B24.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-129, REQ-IND-61 y REQ-MAN-91.
  - Cambia: escritor determinista de TSV+JSON con origen generado, parámetros/perfiles efectivos, candidatos/descartes/activos, subsemillas, fuente/versión y digests, usando publicación segura.
  - Se comprueba con: dos generaciones de igual referencia/configuración/semilla con mismos bytes, hashes externos, fallo inyectado sin pareja parcial y lector independiente que reconstruye todos los metadatos.

- ☐ B26. Ofrecer generación independiente y reutilización sin reescribir la matriz.
  - Depende de: B17 y B25.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-CONF-128, REQ-CONF-129, REQ-IND-60 y REQ-IND-61.
  - Cambia: orden `generate-snv-matrix config.toml --output ...` y consumo posterior de la pareja canónica; el recurso conserva bytes/digest/origen y el manifiesto consumidor registra su reutilización.
  - Se comprueba con: generación autónoma e integrada equivalentes byte a byte, segunda ejecución con individuos distintos que no modifica `mtime`/hash del recurso y rechazo de pareja parcial o mezclada.

- ☐ B27. Validar y aplicar SNP explícitos como restricciones trazables, no como otra fuente aleatoria inicial.
  - Depende de: A8–A9 y B6–B7.
  - Procede de: REQ-IND-05, REQ-IND-11 y REQ-REC-24; DTC-020, DTC-032 y DTC-052.
  - Cambia: evento SNP explícito con referencia esperada, ALT `A/C/G/T`, copia y procedencia; prepara su contraste posterior con la matriz y rechaza origen/destino `N`.
  - Se comprueba con: transición/transversión, desacuerdo, copia ausente, borde y `N`, más caso negativo que no acepta una frecuencia poblacional como valor literal individual.

- ☐ B28. Validar y aplicar inserciones pequeñas explícitas y aleatorias.
  - Depende de: A5, A8–A9, B2, B6–B11.
  - Procede de: REQ-IND-05, REQ-IND-40 y REQ-REC-24; DTC-032.
  - Cambia: inserción interna anclada en la base precedente; la explícita conserva `N` literal y la aleatoria genera solo `A/C/G/T` con longitud exacta efectiva.
  - Se comprueba con: oráculo manual, terminal imposible, `N` explícita, longitudes 1/2/`umbral - 1` y barrido aleatorio sin `N`.

- ☐ B29. Validar y aplicar deleciones pequeñas solo con bases flanqueantes.
  - Depende de: A5, A8–A9, B2 y B6–B11.
  - Procede de: REQ-IND-05, REQ-IND-41 y REQ-REC-24.
  - Cambia: deleción con intervalo/REF explícitos o longitud exacta aleatoria; rechaza cero, desacuerdo, solapamiento inválido y cualquier deleción terminal.
  - Se comprueba con: cadenas cortas, longitudes 1/2/`umbral - 1`, intentos sobre primera/última/toda la secuencia y estado intacto tras error.

- ☐ B30. Componer eventos pequeños sobre el estado vigente con conflictos explícitos.
  - Depende de: B27–B29.
  - Procede de: DTC-020 y DTC-031; REQ-REC-24 y REQ-IND-21.
  - Cambia: aplicación secuencial de explícitos y demás clases aleatorias en orden normativo, con coordenada de referencia, alelo vigente y diagnóstico de objetivo inexistente; la matriz todavía no se muestrea.
  - Se comprueba con: SNP→deleción, deleción→SNP, inserciones vecinas y repetidas con resultado manual y prueba de ausencia de edición parcial.

- ☐ B31. Cerrar diagnósticos y preflight de las entradas F2.2.
  - Depende de: B6–B30.
  - Procede de: DTC-044; REQ-CONF-25, REQ-CONF-28, REQ-IND-40, REQ-IND-41 y REQ-VAL-55.
  - Cambia: catálogo estable y agregado para formato/digest, referencia, población, coordenada, alelo, espectro/intervalo, región/peso, recuentos, `N`, *indel*, matriz vacía y conflicto.
  - Se comprueba con: tabla diagnóstico→caso, texto/JSON equivalentes, orden determinista, flags limitados y configuración inválida que no crea artefactos finales.

- ☐ B32. Construir un validador independiente de los formatos, digests y normalización.
  - Depende de: B3–B4; no depende de B17–B26 para calcular el esperado.
  - Procede de: REQ-VAL-55 y ACE-116; DTC-047 y DTC-054 (resuelve DTPD-118).
  - Cambia: corpus escrito a mano y lector separado para TSV+JSON, referencia, filas completas, IDs, orden, multialelia, sumas Neumaier, tolerancia y eliminación de ALT cero.
  - Se comprueba con: mutación deliberada de cada campo/digest/invariante y reconstrucción de resultados literales sin llamar al parser, normalizador ni escritor de producción.

- ☐ B33. Construir oráculos independientes para selección, espectros y perfiles.
  - Depende de: B8–B15; no depende de B22–B25 para producir el esperado.
  - Procede de: REQ-VAL-55 y ACE-116; DTC-047 y DTC-054 (resuelve DTPD-118).
  - Cambia: casos sintéticos y cálculo separado de elegibilidad ponderada, Poisson, mayores restos, muestreo de intervalos, ti/tv y resúmenes 1000 Genomes.
  - Se comprueba con: vectores literales, semillas fijas, pruebas estadísticas con umbrales previos y mutaciones del productor que el oráculo detecta sin reutilizar sus funciones.

- ☐ B34. Integrar fuente, matriz efectiva, digest y procedencia sin crear individuos.
  - Depende de: B17–B26 y B31–B33.
  - Procede de: REQ-CONF-128, REQ-CONF-129, REQ-IND-02, REQ-IND-06, REQ-IND-60 y REQ-IND-61; DTC-040 y DTC-052.
  - Cambia: resultado F2.2 inmutable con referencia, poblaciones disponibles/seleccionadas, matriz preparada/generada/reutilizada, perfiles, candidatos/descartes/activos, subsemillas, versiones y SHA-256 para F2.3.
  - Se comprueba con: las tres rutas, referencias cruzadas, `sha256sum`, igualdad con semilla fija, poblaciones adicionales válidas y ausencia demostrada de genotipos o SNP individuales sorteados.

- ☐ B35. Completar evidencia y trazabilidad de F2.2 con oráculos ajenos al productor.
  - Depende de: B8–B34 y puertas de calidad verdes.
  - Procede de: los 26 REQ asignados a F2.2 en §15; ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109, ACE-110, ACE-111 y ACE-116.
  - Cambia: pruebas unitarias, de propiedades, estadísticas e integración y matriz de evidencia de los 22 REQ previos más CONF-128/129 e IND-60/61.
  - Se comprueba con: §4.3, recuento automático de 26 REQ sin huérfanos, revisión de B32/B33 y repetición exacta; DTC-052 aporta semántica a RAF2-004, pero no lo cierra sin revisar formalmente la matriz maestra.

> **Pausa obligatoria F2.2.** Confirmar el consumo sin reinterpretación de
> DTC-054 y el cierre de DTPD-118, la pareja TSV+JSON y sus digests, perfiles solo explícitos,
> candidatos/descartes/activos, poblaciones independientes, validadores,
> Poisson limitado a recuentos, espacios elegibles y casos con `N`. No iniciar
> F2.3 si un oráculo se deriva del productor, una densidad SNP individual sigue
> activa, un flag relaja otra validación o un rechazo deja cambios.

## 8. F2.3 — Individuos directos desde una población de origen

**Resultado del hito:** la vía directa muestrea una matriz común mediante HWE o
`F` diploide, conserva fase y valores esperados/realizados, integra los demás
cambios y reconstruye cada haplotipo sin simular dinámica poblacional.

- ☐ C1. Extender la autoría de individuos, población de origen y modelo genotípico.
  - Depende de: pausa F2.2 aprobada y B1–B6.
  - Procede de: REQ-CONF-20, REQ-CONF-22, REQ-CONF-130, REQ-DIR-01 y REQ-DIR-02; DTC-034 y DTC-052.
  - Cambia: esquema `individuals` para ruta `direct`, exactamente una población de origen disponible en la matriz por individuo, asignada explícitamente o mediante el reparto configurado de individuos entre poblaciones, ploidía uniforme, HWE predeterminado o `F` diploide; no admite mezcla/admixture ni evolución.
  - Se comprueba con: ejemplos explícitos y por reparto, HWE y `F`, cero/uno/varios individuos, población desconocida o múltiple, suma/cantidad de reparto inválida, F fuera de `[0,1]` y F en ploidía distinta de 2.

- ☐ C2. Parsear, validar y resolver individuos con IDs tipados y procedencia completa.
  - Depende de: C1 y contratos F1 de identidad/configuración.
  - Procede de: DTC-033, DTC-039 y DTC-044; REQ-CONF-20, REQ-CONF-22, REQ-CONF-130 y evidencia compartida de REQ-IND-25/REQ-MAN-10.
  - Cambia: configuración declarada/efectiva, identidades, población, modelo y solicitudes de resolución, conservando valor, regla, selectores y origen.
  - Se comprueba con: IDs repetidos/desconocidos, orden, sustitución individual, HWE omitido/predeterminado y aislamiento de valores efectivos.

- ☐ C3. Versionar subsemillas de genotipo y fase separadas de la generación de matriz.
  - Depende de: B8, C1–C2.
  - Procede de: REQ-IND-63 y REQ-VAL-55; DTC-011, DTC-037, DTC-038 y DTC-052.
  - Cambia: ámbitos `individual/genotype/<population_id>/<individual_id>/<variant_id>` e `individual/phase/<individual_id>/<variant_id>` que consumen la matriz inmutable sin compartir RNG mutable.
  - Se comprueba con: vectores literales, mismo resultado con distinto orden físico y prueba de estabilidad de matriz e individuos existentes al añadir otro ID.

- ☐ C4. Construir el estado inicial faseado de cada individuo desde la referencia.
  - Depende de: A6–A9 y C2.
  - Procede de: DTC-003 y DTC-006; REQ-DIR-01, REQ-DIR-02, REQ-IND-28 y REQ-IND-29.
  - Cambia: copias intactas ordenadas por secuencia/ploidía y enlaces a referencia, población, individuo y procedencia antes de asignar alelos.
  - Se comprueba con: copias 1–10 iguales byte a byte a la referencia, cardinalidad/IDs exactos y ausencia de escenarios o poblaciones simuladas.

- ☐ C5. Muestrear alelos HWE por locus y copia conservando fase.
  - Depende de: B34 y C3–C4.
  - Procede de: REQ-CONF-130 y REQ-IND-62; DTC-052.
  - Cambia: muestreo independiente de la composición genotípica desde frecuencias REF/ALT por población y locus; `alt_af=0/1` es determinista y ningún recuento se corrige para acercar frecuencias realizadas al objetivo.
  - Se comprueba con: corpus bialélico/multialélico y ploidías 1/2/10, frecuencias 0/intermedias/1, oráculo de probabilidades/vectores y prueba de ausencia de Poisson o autocorrección.

- ☐ C6. Muestrear genotipos diploides con `F` usando las fórmulas normativas.
  - Depende de: B34 y C3–C4.
  - Procede de: REQ-CONF-130, REQ-IND-62 y REQ-EVO-47/48; DTC-052.
  - Cambia: distribución bialélica/multialélica diploide para `F` efectivo 0–1, con asignación faseada a dos copias; no generaliza `F` a otras ploidías.
  - Se comprueba con: probabilidades calculadas fuera del motor para F=0/intermedio/1, suma/normalización, casos multialélicos y rechazo haploide/poliploide.

- ☐ C7. Asignar fase independiente después de sortear la composición genotípica.
  - Depende de: C3 y C5–C6.
  - Procede de: REQ-IND-62, REQ-IND-63 y REQ-VAL-55; DTC-052 y DTC-054 (resuelve DTPD-118).
  - Cambia: permutación reproducible de alelos entre las copias para heterocigotos diploides, multialélicos y poliploides; la fase es conocida pero no representa LD ni haplotipos poblacionales.
  - Se comprueba con: genotipos 0/1, 1/1, 1/2 y poliploides, vectores de fase separados y prueba de que cambiar solo fase no vuelve a sortear genotipo ni correlaciona loci.

- ☐ C8. Contrastar SNP explícitos con matriz y modelo antes de modificar copias.
  - Depende de: B27, C2 y C5–C7.
  - Procede de: REQ-IND-64 y REQ-DIR-02; DTC-020 y DTC-052.
  - Cambia: preflight transaccional que acepta restricciones compatibles y diagnostica probabilidad cero, ALT ajeno al sitio, dosis/genotipo contradictorio, población ausente y densidad SNP individual; no modifica la matriz ni sus frecuencias.
  - Se comprueba con: tabla matriz/modelo/restricción para `alt_af=0/1`, F=1, ALT nuevo, dosis imposible y contradicción tardía sin individuo parcial.

- ☐ C9. Aplicar solo la excepción `--allow-incompatible-genotypes` aprobada.
  - Depende de: C8 y catálogo de diagnósticos F1.
  - Procede de: DTC-054 (resuelve DTPD-118); REQ-IND-64, REQ-VAL-55 y REQ-MAN-91.
  - Cambia: sin flag, toda incompatibilidad de C8 bloquea; con flag, conserva el genotipo explícito y registra individuo/locus/alelos/motivo, sin relajar referencia, coordenadas, ploidía/dosis imposible, formato, sumas ni no-clobber.
  - Se comprueba con: cada incompatibilidad permitida con/sin flag y matriz de errores no anulables que siguen fallando; el manifiesto y las frecuencias realizadas incluyen el genotipo forzado.

- ☐ C10. Integrar *indels* y los demás cambios directos después de fijar genotipos SNV.
  - Depende de: B28–B31, C4 y C8–C9.
  - Procede de: REQ-IND-15, REQ-DIR-01, REQ-DIR-02 y REQ-REC-24; DTC-020 y DTC-031.
  - Cambia: orquestación ordenada por individuo/copia sobre el estado faseado, con transacción lógica y RNG locales para clases no SNV iniciales.
  - Se comprueba con: caso manual matriz+SNP restringido+inserción+deleción, permutación controlada y fallo tardío sin publicación parcial.

- ☐ C11. Materializar y validar cada secuencia haplotípica final contra su historial.
  - Depende de: C5–C10 y A9.
  - Procede de: DTC-003 y DTC-031; REQ-IND-15 y REQ-DIR-01.
  - Cambia: resultado con secuencia final, eventos, fase, dosis y procedencia sin convertir bytes materializados en una segunda verdad.
  - Se comprueba con: aplicación de eventos y ensamblado independiente que producen los mismos bytes/longitudes y conservan `N` por copia.

- ☐ C12. Integrar múltiples individuos sin acoplar matriz, genotipos ni cambios locales.
  - Depende de: C2–C11.
  - Procede de: REQ-IND-63, REQ-CONF-22 y REQ-DIR-01; DTC-006 y DTC-038.
  - Cambia: recorrido ordenado de varios individuos con matriz común inmutable y valores/verdad separados por identidad.
  - Se comprueba con: al menos tres individuos, adición/eliminación/reordenación, ejecución física distinta y comparación exacta de matriz e IDs previos.

- ☐ C13. Calcular objetivos y realizados sin autocorregir el muestreo finito.
  - Depende de: C5–C12.
  - Procede de: REQ-IND-65 y REQ-VAL-55; DTC-049 y DTC-052.
  - Cambia: agregados canónicos por población/locus de `alt_af` objetivo, copias totales/ALT, `alt_af` realizada y genotipos esperados/realizados, incluyendo individuos forzados.
  - Se comprueba con: cálculo manual sobre corpus pequeño, orden/suma canónicos, desviación finita informada, extremos deterministas y prueba de que ningún genotipo se modifica para acercarla al objetivo.

- ☐ C14. Integrar procedencia de matriz y muestreo en configuración efectiva y manifiesto.
  - Depende de: B34 y C2–C13.
  - Procede de: REQ-MAN-91, REQ-IND-65 y evidencia compartida de REQ-MAN-10; DTC-040 y DTC-052.
  - Cambia: manifiesto con modo/ruta/SHA-256 de matriz y JSON, referencia, poblaciones disponibles/seleccionadas, parámetros/perfiles/fuente/versión, semillas, HWE/F, autorizaciones y valores objetivo/realizados.
  - Se comprueba con: inspección campo a campo, `sha256sum`, enlaces resolubles y ejecuciones preparada/generada/reutilizada/forzada que no omiten procedencia declarada o efectiva.

- ☐ C15. Validar estadísticamente HWE/`F` con un corpus y cálculos independientes.
  - Depende de: C5–C7 y C13; no usa C11 para calcular el esperado.
  - Procede de: REQ-VAL-55 y ACE-116; DTC-047 y DTC-052.
  - Cambia: pruebas de aceptación con probabilidades y umbrales predeclarados para HWE/`F`, fase y aislamiento de subsemillas, separadas de la validez formal B32.
  - Se comprueba con: muestras reproducibles de tamaños pequeño/grande, comparación de esperado/realizado sin exigir igualdad exacta, límites F=0/1 y fallo deliberado detectado por el oráculo.

- ☐ C16. Completar evidencia y trazabilidad de F2.3.
  - Depende de: C2–C15 y puertas de calidad verdes.
  - Procede de: los 13 REQ exclusivos y dos compartidos asignados a F2.3 en §15/§15.1; ACE-01, ACE-10, ACE-12, ACE-22, ACE-23, ACE-28, ACE-104, ACE-111, ACE-114 y ACE-116.
  - Cambia: matriz de evidencia y escenario sintético reproducible con matriz, genotipos, fase, haplotipos, métricas y manifiesto calculados independientemente.
  - Se comprueba con: §4.3, recuento automático de 13 exclusivos y dos compartidos sin huérfanos, ID estable omitido y ausencia de dinámica poblacional, lecturas o formatos técnicos.

- ☐ C17. Activar matriz, población y creación directa completa en Marimo.
  - Depende de: pausa F2.0 y C1–C16.
  - Procede de: REQ-UI-01 a REQ-UI-07, REQ-UI-13, REQ-UI-14 y DTC-051/DTC-052.
  - Cambia: pestaña `Individuos` para preparar/vincular/verificar referencia, generar/previsualizar/reutilizar matriz, elegir población/perfiles/modos/regiones/HWE/F y confirmar matriz vacía o genotipos incompatibles; invoca al binario sin reimplementar reglas.
  - Se comprueba con: matriz campo→control, SHA/REF y resúmenes por cromosoma/sustitución/población/intervalo, accesibilidad e igualdad UI→TOML→`EffectiveConfig`→artefactos; confirmaciones producen solo los flags específicos y una edición invalida la ejecución.

> **Pausa obligatoria F2.3.** Reconstruir manualmente todos los haplotipos,
> aprobar HWE/`F`, fase, métricas, manifiesto y aislamiento de subsemillas. No
> iniciar F2.4 si Poisson aparece en muestreo de alelos, añadir un individuo
> altera a otro, una excepción limitada oculta otra validación o C17 deja una
> opción F2.3 solo en TOML.

## 9. F2.4 — Verdad estructurada independiente del VCF

**Resultado del hito:** una ejecución directa válida publica un fichero de
verdad estructurada versionado que basta para reconstruir referencia, individuos,
copias, fase, dosis, eventos y procedencia aunque no se produzca ningún VCF.

- ☐ D1. Formalizar el esquema externo aprobado por DTPD-102 y su contrato de compatibilidad.
  - Depende de: pausa F2.3 aprobada y A3.
  - Procede de: DTPD-102; DTC-031, DTC-035 y DTC-036; REQ-VER-31 y REQ-SAL-06.
  - Cambia: esquema versionado del fichero de verdad con campos obligatorios, tipos, orden normativo, coordenadas, unidades y política ante versión/claves desconocidas; no reutiliza el esquema del estado interno.
  - Se comprueba con: ejemplos válidos e inválidos, recorrido campo→REQ y validación mediante una implementación de esquema independiente elegida y registrada.

- ☐ D2. Implementar una proyección pura desde la verdad canónica al modelo externo.
  - Depende de: C7–C11 y D1.
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

- ☐ D9. Añadir controles y resumen enlazado de verdad estructurada a Marimo.
  - Depende de: pausa F2.0 y D1–D8.
  - Procede de: REQ-UI-01, REQ-UI-07 a REQ-UI-13 y DTC-051; F2.4 de la hoja de ruta.
  - Cambia: controles de activación/formato ya aprobado en `Salidas`, revisión del TOML y resumen desde manifiesto en `Revisar y ejecutar`, sin cargar el fichero de verdad completo.
  - Se comprueba con: configuración UI/manual equivalente con VCF activado/desactivado, validación invalidada tras cambio, manifiesto grande simulado y prueba de que la UI solo lee metadatos/enlaces.

> **Pausa obligatoria F2.4.** Entregar el fichero y la referencia a una revisión
> que no use el productor y exigir reconstrucción exacta. No iniciar F2.5 si el
> VCF es necesario para interpretar la verdad, si falta fase/dosis/procedencia o
> si una corrupción puede parecer una ejecución completa. D9 debe demostrar que
> la app resume y enlaza la verdad sin cargarla ni interpretarla.

## 10. F2.5 — VCF de verdad mínimo

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
  - Procede de: los 25 REQ de cierre asignados a F2.5 en la matriz de §15 y evidencia compartida de REQ-MAN-83; ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100 y ACE-103.
  - Cambia: matriz de evidencia F2.5, informe de validación independiente, inventario de limitaciones y dossier de cierre; separa la evidencia compartida y no marca requisitos de fases posteriores como completos.
  - Se comprueba con: §4.3, recuento automático de 25 REQ de cierre y un REQ compartido sin huérfanos, comandos registrados de todos los validadores y revisión cruzada de los 84 REQ de cierre más los cuatro REQ compartidos de §15.

- ☐ E13. Incorporar todas las opciones VCF F2.5 a Marimo y cerrar la paridad F2.
  - Depende de: pausa F2.0 y E1–E12.
  - Procede de: REQ-UI-01, REQ-UI-03, REQ-UI-07 a REQ-UI-13 y DTC-051; F2.5 de la hoja de ruta.
  - Cambia: controles `Salidas` para VCF BGZF/indexado predeterminado o texto, revisión de artefactos esperados, validación/lanzamiento explícitos y resumen de manifiesto/VCF sin cargar registros completos.
  - Se comprueba con: pares UI/manual para ambos modos con igual `EffectiveConfig`, destino existente, edición tras validar, ejecución real mínima y matriz de todas las opciones configurables F2 sin huecos Marimo.

> **Pausa obligatoria F2.5 y cierre de F2.** No autorizar F3 hasta que una
> revisión independiente reconstruya todos los haplotipos, confirme
> coordenadas/alelos/fase/dosis, valide VCF/BGZF/índice y demuestre los cuatro
> casos de `N`: SNP rechazado, inserción aleatoria sin `N`, inserción explícita
> con `N` conservada y VCF/complemento sin resolución inventada. Aprobar por
> separado el cierre de F2; completar E12/E13 no lo concede automáticamente y
> ninguna opción configurable F2 puede quedar solo en TOML.

## 11. Dependencias posteriores confirmadas, fuera de la ejecución de F2

Esta sección no añade hitos ni casillas a F2. Conserva el destino explícito de
DTC-053 sin anticipar la planificación completa de F6 o F7:

| Destino | Contrato que deberá consumir | Puertas previas | Evidencia mínima futura |
|---|---|---|---|
| F6.4, panel fundador intacto o recombinante | REQ-EVO-74 a REQ-EVO-77; REQ-VAL-56, REQ-MAN-92 y REQ-UI-15 | Consumir íntegramente DTC-055: paquete/TSV, referencia única, mínimos, ejes/modelos, productor directo, constructor común uniforme, puntos x\|x+1, F, IDs DTC-039, RNG DTC-037, origen/proceso y salidas canónicas | Probar catálogo/ejes y rechazos de combinaciones, incluida distribución explícita solo en intacto; ejecutar Bernoulli(F), productores y fuentes solo del panel; probar referencia/secuencias, K, puntos no SNP/invisibles; verificar C, ID, `n_h/C`, recuentos por origen y los tres esquemas/órdenes/claves TSV; reconstruir con oráculo independiente y corrupción deliberada; publicar/enlazar las tres tablas; probar Marimo, validación vigente y paridad TOML |
| F7.1, recombinación meiótica | REQ-REC y requisitos reproductivos aplicables | F6.4 cerrado y DTPD-101/DTPD-111 resueltas | Conectar los mapas y la semántica meiótica decididos por DTPD-111 al constructor común ya implementado en F6.4 y demostrar que no reinterpreta la inicialización fundadora como evento histórico ni crea una segunda implementación |

DTPD-119 está resuelta mediante DTC-055. Quedan decididos el paquete y su preparación, el
alcance SNP, bloques e identidades locales, tablas y vectores densos/únicos, los
modelos y modos del muestreo intacto, la asignación a copias, el uso de la matriz
para loci exteriores y el rechazo de cualquier locus `alt_af` dentro de un
bloque. La alternativa recombinante se ha trasladado a la formación del conjunto
fundador F6.4 y ya tiene productor directo, puntos, F, IDs, RNG, tres salidas y
Marimo definidos, incluida la sintaxis TSV, la semántica de origen/proceso y el
rechazo de distribuciones genotípicas explícitas en recombinante.
Ningún paso B/C
puede aceptar un panel,
inferir haplotipos desde la matriz marginal ni implementar recombinación local.
Marimo solo expondrá REQ-UI-15 cuando cierren las capacidades correspondientes.

## 12. Riesgos genómicos silenciosos y oráculos obligatorios

| Riesgo silencioso | Prevención en producción | Oráculo independiente obligatorio | Evidencia de fallo que debe detectarse |
|---|---|---|---|
| Desfase base 0/base 1 o extremo abierto/cerrado | Tipos/servicio DTPD-101 iniciados en F1.6 y extendidos en A5 | Vectores compartidos de mapa y cadenas etiquetadas por posición | Mapa, SNP o indel desplazado, ancla incorrecta o acceso a posición 0 |
| Longitud de *indel* mal resuelta o convertida en distribución implícita | Parámetros exactos B2/B7 y constructores B28–B29 | Casos 1/2/`umbral - 1` y RNG instrumentado | Longitud distinta, consumo RNG para longitud o clasificación estructural incorrecta |
| Inserción/deleción terminal inválida | Constructores validados B28–B29 | Enumeración manual de todos los indels de una cadena corta | Pérdida del ancla o secuencia sin base a uno de los lados |
| Sustitución silenciosa de `N` | Dominio A6 y reglas B9/B17/B22–B29/E7 | Corpus con `N` en borde/interior y barrido de salidas | SNP `REF=N`, ALT `N` o inserción aleatoria con `N` |
| Fase o dosis alterada | Copias ordenadas A7/A9 y GT E5 | Conteo manual por copia, no parsing circular del VCF | GT sin P alelos, separador no faseado o dosis distinta |
| Evento perdido por solapamiento | Composición B30 e historial D2 | Aplicador de cadenas corto escrito desde fixtures | SNP superviviente tras ser borrado o conflicto aceptado |
| Tasa cero que consume RNG | Rutas explícitas B11/B20/B24 | RNG instrumentado y comparación de ámbitos vecinos | Cambios en otras clases o individuos al activar cero |
| Poisson incorrecto o aplicado a alelos/genotipos | Decisión B2 y servicio B11/C5–C6 | Vectores o implementación ajena y test estadístico prefijado | Media/varianza fuera del umbral, Poisson fuera del recuento o resultado no reproducible |
| Matriz o espectro incompatibles con referencia/población | Esquemas/validación DTC-054 y B3–B4/B12/B17–B21 | Corpus/lector B32 calculado fuera del productor | REF/ALT, suma, duplicado, población, perfil o digest incorrectos aceptados |
| Recuento candidato confundido con loci activos | Separación B11/B20/B25 | Tabla manual candidatos→ALT cero→activos | Reposición de un locus globalmente ausente o `exact_loci` interpretado como activos |
| Distribución o perfil implícito/incorrecto | Perfiles B10/B12–B15 y selección B22–B24 | Oráculos B33 y regeneración 1000 Genomes | Solapamiento aceptado, mayor resto incorrecto o perfil empírico cargado por defecto |
| Identidad o azar cambia por orden/adición | IDs B16 y ámbitos B8/C3 | Permutaciones y vectores literales | Reordenar filas/poblaciones o añadir individuo altera resultados existentes |
| Excepción limitada relaja otra validación | Flags B21/C9 y preflight B31 | Matriz independiente error×flag | Digest, REF, frecuencia, ploidía/dosis o no-clobber inválidos aceptados |
| HWE/`F`, fase o frecuencia realizada incorrectos | Muestreo C5–C9 y agregados C13 | Probabilidades/corpus independiente C15 | Genotipo fuera del modelo, fase cambiada o muestra finita autocorregida |
| Verdad y secuencia materializada divergen | Una verdad canónica A8–A9/C11 | Segundo reconstructor de D4 | Bytes o longitud distintos para la misma copia |
| VCF convertido en fuente de verdad | Dirección de dependencias D2/E3 | Prueba sin VCF y auditoría de módulos | No se puede reconstruir o cambia la verdad al desactivar VCF |
| REF/ALT/normalización incorrectos | Proyección E3–E4 | Normalizador/validador independiente aprobado | Registro válido sintácticamente pero no representa la verdad |
| Pérdida de detalle al agrupar multialélicos | Agrupación E5 y complemento E7 | Tabla manual sitio×copia y reconstructor combinado | ALT o copia portadora desaparece, o se inventa una base |
| Artefactos mezclados entre ejecuciones | IDs/hashes/publicación F1 reutilizados | `sha256sum`, referencias cruzadas y fallos inyectados | Índice, VCF, verdad o referencia de ejecuciones distintas |
| UI y conjunto TOML producen semánticas distintas, rompen enlaces o ejecutan una revisión obsoleta | Modelo/digest multidocumento M4/M7, binario único M8/M9 y pasos C17/D9/E13 | Conjunto manual gemelo, escenarios faltante/duplicado y proceso falso que registra argumentos | `EffectiveConfig` distinta, ruta/orden perdido, ejecución automática o validación antigua reutilizada |

Una prueba que serializa y vuelve a leer con el mismo código solo demuestra
consistencia interna. No satisface por sí sola ninguna fila de esta tabla.

## 13. Cambios pequeños, revisión y política de bloqueo

- Cada paso se implementa y revisa de forma aislada; una revisión puede dividirlo
  aún más, pero no fusionarlo con un paso dependiente para ocultar su evidencia.
- Si M1, A2, A3, B2 o E1 no obtienen aprobación, o B1 no confirma sin
  contradicción DTC-054 y el cierre de DTPD-118, se detienen sus consumidores y se
  actualiza este plan con el bloqueo. Se puede continuar únicamente con corpus,
  matrices o validadores que no presupongan la opción elegida.
- Si una decisión aprobada contradice especificación, arquitectura o DTC
  confirmada, no se adapta el requisito al código: se pausa y se eleva la
  contradicción para decisión explícita.
- Las mediciones pequeñas observan regresiones y consumo, pero no prometen los
  objetivos de escalabilidad de F9 ni justifican optimizaciones especulativas.
- Ninguna pausa se supera solo porque compila. Requiere revisión de contratos,
  pruebas negativas, oráculo independiente y matriz REQ/ACE actualizada.

## 14. Validación mecánica del documento

Antes de aprobar el plan se ejecutarán, como mínimo:

```bash
git diff --check -- design/plan-implementacion-f2-0.1.md
rg -n '^- ☐ [A-EM][0-9]+\.' design/plan-implementacion-f2-0.1.md
rg -n '^> \*\*Pausa obligatoria F2\.[0-5]' design/plan-implementacion-f2-0.1.md
rg -n 'DTC-052|DTC-053|DTC-054|DTC-055|DTPD-101|DTPD-102|DTPD-103|DTPD-117|DTPD-118|DTPD-119|F1\.5|F1\.6' design/plan-implementacion-f2-0.1.md
rg -n 'RAF2-001|RAF2-004|RAF2-007' design/plan-implementacion-f2-0.1.md
```

Además, un comprobador mecánico debe verificar que cada casilla contiene una y
solo una línea `Depende de`, `Procede de`, `Cambia` y `Se comprueba con`; que la
matriz siguiente contiene una sola vez los 84 REQ cuyo cierre está asignado
exclusivamente a F2; y que la lista separada de evidencia compartida contiene
los otros cuatro REQ cuyo campo de hito incluye F2 junto con fases posteriores.
La unión debe coincidir exactamente con las 88 filas de la trazabilidad cuyo
hito contiene F2, sin confundir evidencia parcial con cierre total.

## 15. Matriz de trazabilidad REQ/ACE/DTC

La lista de REQ de cierre exclusivo es exhaustiva para la asignación F2 vigente:
14 + 4 + 26 + 13 + 2 + 25 = **84**. Se añaden por separado **cuatro REQ compartidos**
que exigen evidencia parcial en F2 y cierre posterior. Los ACE son la unión de
los criterios enlazados por esas filas de trazabilidad. Una mención no acredita
ejecución ni cierre. REQ-CONF-08, REQ-CONF-21, REQ-IND-02 y REQ-IND-03 siguen
señalados como cierres exclusivos solo porque así consta en la matriz maestra;
RAF2-004 impide marcarlos satisfechos hasta revisar esa asignación frente al
alcance deliberado de F2. La auditoría de la matriz maestra asigna los doce REQ
nuevos como cuatro a F2.2 y ocho a F2.3; no se adopta el reparto 5/7 porque
contradiría esas filas sin una revisión explícita de trazabilidad.

| Hito | REQ asignados | ACE cubiertos | Decisiones/puertas principales | Pasos de cierre |
|---|---|---|---|---|
| F2.0 | REQ-UI-01; REQ-UI-02; REQ-UI-03; REQ-UI-04; REQ-UI-05; REQ-UI-06; REQ-UI-07; REQ-UI-08; REQ-UI-09; REQ-UI-10; REQ-UI-11; REQ-UI-12; REQ-UI-13; REQ-VAL-54 | ACE-114, ACE-115 | DTC-034, DTC-036, DTC-042, DTC-043, DTC-044, DTC-045, DTC-048, DTC-051; DTPD-117 | M1–M12 |
| F2.1 | REQ-REF-14; REQ-IND-28; REQ-IND-29; REQ-VER-27 | ACE-03, ACE-10, ACE-11, ACE-28, ACE-66 | DTC-003, DTC-010, DTC-013, DTC-031, DTC-032, DTC-039 y DTC-056; DTPD-102 | A2–A10 |
| F2.2 | REQ-CONF-08; REQ-CONF-21; REQ-CONF-24; REQ-CONF-25; REQ-CONF-28; REQ-CONF-128; REQ-CONF-129; REQ-IND-02; REQ-IND-03; REQ-IND-05; REQ-IND-06; REQ-IND-10; REQ-IND-11; REQ-IND-13; REQ-IND-14; REQ-IND-17; REQ-IND-18; REQ-IND-19; REQ-IND-21; REQ-IND-38; REQ-IND-40; REQ-IND-41; REQ-IND-60; REQ-IND-61; REQ-DIR-03; REQ-REC-24 | ACE-01, ACE-09, ACE-10, ACE-21, ACE-22, ACE-28, ACE-85, ACE-95, ACE-103, ACE-104, ACE-109, ACE-110, ACE-111, ACE-116 | DTC-006, DTC-011, DTC-020, DTC-022, DTC-032, DTC-033, DTC-037, DTC-038, DTC-044, DTC-049, DTC-052, DTC-054 y DTC-056; decisión B2 aún sin identificador; DTPD-118 resuelta | B1–B35 |
| F2.3 | REQ-CONF-20; REQ-CONF-22; REQ-CONF-130; REQ-IND-15; REQ-IND-62; REQ-IND-63; REQ-IND-64; REQ-IND-65; REQ-DIR-01; REQ-DIR-02; REQ-UI-14; REQ-VAL-55; REQ-MAN-91 | ACE-01, ACE-10, ACE-12, ACE-22, ACE-23, ACE-28, ACE-104, ACE-111, ACE-114, ACE-116 | DTC-003, DTC-006, DTC-011, DTC-020, DTC-037, DTC-038, DTC-039, DTC-049, DTC-051, DTC-052 y DTC-054; DTPD-118 resuelta | C1–C17 |
| F2.4 | REQ-VER-31; REQ-SAL-06 | ACE-34, ACE-62 | DTC-003, DTC-004, DTC-031, DTC-035, DTC-036, DTC-040, DTC-045, DTC-047, DTC-051; DTPD-102 | D1–D9 |
| F2.5 | REQ-FMT-05; REQ-FMT-06; REQ-REC-26; REQ-VER-02; REQ-VER-06; REQ-VER-12; REQ-VER-15; REQ-VER-16; REQ-VER-17; REQ-VER-30; REQ-VER-32; REQ-VER-33; REQ-VER-34; REQ-VER-36; REQ-VER-37; REQ-VER-39; REQ-VER-40; REQ-VER-41; REQ-VER-43; REQ-VER-44; REQ-VER-47; REQ-VAL-28; REQ-MAN-27; REQ-MAN-57; REQ-SAL-21 | ACE-19, ACE-25, ACE-28, ACE-29, ACE-30, ACE-60, ACE-61, ACE-62, ACE-63, ACE-64, ACE-66, ACE-67, ACE-93, ACE-94, ACE-100, ACE-103 | DTC-003, DTC-004, DTC-013, DTC-021, DTC-032, DTC-040, DTC-043, DTC-045, DTC-047, DTC-051; DTPD-103 | E1–E13 |

### 15.1. Evidencia compartida obligatoria en F2

Estas cuatro filas tienen un destino compuesto en la matriz maestra. F2 debe
producir la evidencia concreta indicada, pero no declarará el REQ completamente
satisfecho hasta que cierren también sus fases posteriores:

| REQ compartido | Evidencia parcial exigida en F2 | Pasos |
|---|---|---|
| REQ-IND-25 (F2.3/F6.1) | Un creador directo real asigna un ID estable cuando falta, comprueba unicidad y no reutilización dentro del recorrido disponible | C2, C16 |
| REQ-IND-26 (F2.4/F3.5/F7.7/F8.5) | El mismo ID aparece sin traducción en estado genotípico, verdad estructurada y metadatos F2; observaciones y linaje quedan para sus fases | D2, D6, D8 |
| REQ-MAN-10 (F2.3/F5.5/F9.7) | El manifiesto distingue valores globales, específicos y efectivos consumidos por individuos/regiones/clases de F2 | C2, C14, C16 |
| REQ-MAN-83 (F2.5/F3.5) | El primer conjunto real y sus artefactos VCF/índice manifiestan estados y diagnósticos coherentes con ejecución y evento | E9, E12 |

### 15.2. Evidencia adyacente no asignada a F2

Los casos de aceptación de F2 deben conservar evidencia útil sobre `N`,
inserciones aleatorias, VCF y salida directa, pero la matriz no incorporará como
cerrados REQ-IND-58, REQ-IND-59, REQ-VER-01, REQ-VER-24 a REQ-VER-26,
REQ-FMT-07 ni REQ-SAL-02 porque la trazabilidad vigente los asigna a hitos
posteriores con alcance mayor. Cualquier discrepancia entre esa asignación y la
evidencia de cierre pedida por la hoja de ruta se documentará como cobertura
parcial, no se resolverá alterando silenciosamente la fase.

## 16. Criterio de salida del plan

El plan queda listo para aprobación cuando su validación mecánica no detecta
errores, los 84 REQ de cierre aparecen exactamente una vez en §15, los cuatro
REQ compartidos aparecen exactamente una vez en §15.1, cada hito tiene su pausa
y cada paso conserva las cuatro líneas obligatorias. La ejecución seguirá
condicionada al cierre de F1 y a las decisiones M1, A2, A3, B2 y E1, además de
la auditoría de consumo B1. F2 solo podrá
cerrarse cuando E12/E13, la pausa final y la aprobación independiente de cierre
consten por separado.
