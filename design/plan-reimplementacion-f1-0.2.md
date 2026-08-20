# Plan de simplificación y reimplementación de F1 0.2

## 1. Estado y autoridad

- **Estado:** listo para ejecución.
- **Base recuperable:** commit `66ad8b2`.
- **Sustituye para la ejecución:** el recorrido de cierre de
  `plan-implementacion-f1-0.1.md`, que se conserva como historial.
- **Fuentes:** especificación 0.1, arquitectura 0.1, DTC-078, hoja de ruta 0.1 y
  `auditoria-sistematica-cierre-f1-0.1.md`.

Este plan no rehace automáticamente todo el código. Conserva una pieza cuando
cumple el contrato mínimo con menos riesgo que reescribirla; la simplifica o la
reemplaza cuando anticipa capacidades posteriores o mezcla autoridades.

## 2. Resultado y límites

F1 termina con un binario Rust que:

- admite y valida el TOML general necesario para una referencia inicial o
  reducida;
- conserva configuración declarada, validada y efectiva sin dominios
  biológicos sintéticos;
- usa semilla raíz y subsemillas reproducibles para generación y reducción;
- genera o copia una referencia FASTA y puede derivar una referencia reducida;
- publica `effective-config.json`, FASTA, mapa de reducción cuando corresponda
  y un `manifest.json` general con tamaños y SHA-256;
- distingue rechazo, fallo controlado e interrupción mediante provisional,
  `INCOMPLETE` y renombrado atómico no-clobber.

Quedan fuera de F1: población/subpoblación, fase/generación, locus/posición como
jerarquías de reglas, eventos técnicos, manifiestos de evento, normalización
poblacional, `resolution-evidence`, defensa exhaustiva frente a atacante local
activo y revalidación combinatoria del paquete completo. Sus destinos son
F2/F3/F4/F6/F7/F9 conforme a DTC-078.

## 3. Criterios de ejecución

1. Retirar primero obligaciones y rutas sin consumidor; después simplificar el
   código que las sostenía.
2. Mantener interfaces útiles para F2 solo cuando tengan un consumidor real y
   una prueba directa.
3. Probar cada productor con un oráculo pequeño que no replique toda su lógica.
4. No usar el hash para convertir una salida derivada en autoridad semántica.
5. Detenerse en cada pausa; una suite verde no sustituye la revisión del
   alcance.

## 4. F1.1 — Base mínima recuperable

- ✅ F1.1.1. Congelar el inventario de componentes que se conservan, simplifican,
  trasladan o eliminan.
  - Depende de: ninguno.
  - Procede de: DTC-078 y auditoría sistemática §§10–13.
  - Cambia: acta de ejecución F1 0.2 y mapa de módulos/pruebas.
  - Se comprueba con: correspondencia completa entre clasificación A/B/C/D y
    rutas reales; ningún módulo queda sin destino.

- ✅ F1.1.2. Inventariar el grafo objetivo y fijar qué dependencia se retirará
  después de eliminar su último consumidor.
  - Depende de: F1.1.1.
  - Procede de: DTC-043, DTC-048 y DTC-078.
  - Cambia: acta de dependencias; no retira todavía crates que el checkpoint
    necesita para compilar.
  - Se comprueba con: `cargo metadata --locked --format-version 1`,
    `cargo tree --locked -e features` y correspondencia consumidor→crate.

- ✅ F1.1.3. Mantener toolchain, formato, Clippy y una prueba de humo del CLI.
  - Depende de: F1.1.2.
  - Procede de: DTC-042, DTC-047 y DTC-048.
  - Cambia: solo infraestructura o pruebas que hayan quedado acopladas al
    alcance retirado.
  - Se comprueba con: `cargo fmt --all -- --check`, Clippy con `-D warnings`,
    prueba de ayuda/versión y construcción release `--locked`.

> Pausa F1.1: confirmar que el árbol compila, que cada dependencia y familia de
> pruebas trasladada tiene hito propietario y que se ha identificado su último
> consumidor. La ausencia física se exige tras F1.4.5, no en esta pausa.

## 5. F1.2 — Configuración efectiva limitada a consumidores reales

- ✅ F1.2.1. Fijar la superficie TOML exacta que consumen referencia inicial y
  referencia reducida.
  - Depende de: F1.1.3.
  - Procede de: REQ-CONF-01/04/05/14/26 y REQ-CONF-124/125/127;
    DTC-034 y DTC-078. La parte F1 de REQ-CONF-126 es evidencia parcial.
  - Cambia: esquema de autoría, `DeclaredConfig` y DTO de adaptadores. La raíz
    activa contiene exactamente `config_kind`, `schema_version`, `seed`,
    `reference` y `reduced_references`; vía, ploidía, reglas, normalización
    genérica y entidades de fases posteriores se rechazan como desconocidas.
  - Se comprueba con: fixtures válidos y negativos de versión, campo
    desconocido, duplicado, exclusividad y ubicación del diagnóstico.

- ✅ F1.2.2. Construir `EffectiveConfig` únicamente desde valores y selectores
  que F1 aplica.
  - Depende de: F1.2.1.
  - Procede de: DTC-034, DTC-036 y DTC-078.
  - Cambia: `src/configuration.rs`, `src/rules.rs`, `src/application.rs` y tipos
    asociados; conserva únicamente la composición efectiva A/C/G/T, clases y
    longitudes de secuencia y objetivos de reducción. Para un FASTA suministrado,
    la longitud TOML es una aserción que F1.5 contrasta con los bytes, no una
    autoridad independiente.
  - Se comprueba con: referencia simulada/suministrada/reducida, permutación de
    entrada, composición con tolerancia fija `1e-9` y rechazo de cualquier
    selector o parámetro de fases posteriores.

- ✅ F1.2.3. Retirar `ResolutionDomain`, `resolution-evidence` y la resolución de
  jerarquías sin datos reales.
  - Depende de: F1.2.2.
  - Procede de: DTC-076/DTC-077 superadas y DTC-078.
  - Cambia: módulos, esquemas, catálogos de artefactos, escritores, lectores y
    tests dedicados exclusivamente a esa evidencia.
  - Se comprueba con: búsqueda sin referencias runtime residuales y prueba de
    que el CLI mínimo no crea `inputs/resolution/evidence.json`.

- ✅ F1.2.4. Retirar del runtime F1 la composición y semántica de escenarios
  poblacionales.
  - Depende de: F1.2.1–F1.2.3.
  - Procede de: DTC-078; REQ-CONF-09/10 se cierran en F8.1 y la composición
    multidocumento se incorpora en F2.0 cuando tenga consumidor.
  - Cambia: F1 conserva los esquemas publicados para compatibilidad documental,
    pero su runtime no carga escenarios ni crea eventos, fases o poblaciones.
  - Se comprueba con: el recorrido F1 no abre TOML de escenario ni crea unidades
    o manifiestos de evento; una petición fuera de alcance se diagnostica sin
    publicación y no se ignora silenciosamente.

> Pausa F1.2: inspeccionar `EffectiveConfig` y demostrar que todo campo tiene un
> consumidor F1 o una obligación documental explícita.

## 6. F1.3 — Aleatoriedad reproducible mínima

- ✅ F1.3.1. Conservar la semilla raíz de 256 bits y la derivación BLAKE3 estable.
  - Depende de: F1.2.2.
  - Procede de: REQ-CONF-26, REQ-REP-23 y DTC-037/DTC-038; F1 aporta
    evidencia parcial para REQ-CONF-27 y REQ-REP-22, que cierran en F9.7.
  - Cambia: `src/randomness.rs` y registro efectivo de raíz/origen.
  - Se comprueba con: vectores congelados, semilla suministrada/generada y
    rechazo de representación inválida.

- ✅ F1.3.2. Activar en runtime solo los ámbitos `reference_generation` y
  `reference_reduction` sin reescribir catálogos publicados.
  - Depende de: F1.3.1.
  - Procede de: DTC-037, DTC-050, DTC-057 y DTC-078.
  - Cambia: consumidores y catálogo activo; v1–v3 permanecen byte a byte
    compatibles para fases posteriores.
  - Se comprueba con: vectores de ambos ámbitos y prueba de que ningún evento o
    población inexistente consume RNG en F1.

- ✅ F1.3.3. Demostrar repetibilidad lógica y aislamiento entre generación y
  reducción.
  - Depende de: F1.3.2.
  - Procede de: REQ-REP-23 y ACE-23; REQ-REP-22 recibe solo evidencia parcial.
  - Cambia: vectores y pruebas de consumo de los flujos, no el contrato.
  - Se comprueba con: mismas entradas de derivación producen las mismas
    subsemillas y alterar un componente de un ámbito no cambia el otro. La
    repetición integral de FASTA y mapa se ejecuta después, en F1.6.

> Pausa F1.3: revisar vectores y consumo RNG real; no aceptar consumidores
> preparatorios de F2–F9 como evidencia F1.

## 7. F1.4 — Manifiesto general y publicación proporcional

- ✅ F1.4.1. Definir un manifiesto general mínimo y su configuración efectiva.
  - Depende de: F1.2.4 y F1.3.3.
  - Procede de: REQ-MAN-02/08/32, DTC-035/DTC-036 y DTC-078. El
    manifiesto general aporta evidencia parcial a REQ-MAN-01/05/69/72, que
    cierran en sus fases propietarias.
  - Cambia: esquemas y tipos de resultados; elimina eventos, conjuntos,
    `resolution-evidence` y enlaces duplicados sin autoridad independiente.
  - Se comprueba con: JSON Schema externo en pruebas para positivo completo,
    positivo fallido, truncado, campo desconocido y enlace ausente.

- ✅ F1.4.2. Mantener catálogo simple de artefactos con ruta canónica, tamaño y
  SHA-256.
  - Depende de: F1.4.1.
  - Procede de: REQ-MAN-08, REQ-SAL-03/22 y DTC-040/DTC-078.
  - Cambia: `src/results.rs` y `src/integrity.rs`; cada fichero aparece una sola
    vez como artefacto y los enlaces semánticos lo referencian por ID.
  - Se comprueba con: `sha256sum`, cambio de un byte, truncado, ID/ruta duplicada
    y todos los enlaces resolubles.

- ✅ F1.4.3. Simplificar el provisional y conservar publicación atómica
  no-clobber.
  - Depende de: F1.4.2.
  - Procede de: REQ-SAL-22/26, DTC-045 y DTC-078.
  - Cambia: `src/publication.rs` y orquestación; mantiene directorio privado,
    `INCOMPLETE`, mismo sistema de ficheros y rename no-clobber, pero retira
    capacidades, inventario de inodos, hardlinks y hooks TOCTOU de producción.
  - Se comprueba con: destino ausente/preexistente, carrera de dos procesos,
    fallo antes de publicar y SIGKILL que deja `INCOMPLETE`.

- ✅ F1.4.4. Cerrar rechazo, fallo controlado y éxito sin unidades ficticias.
  - Depende de: F1.4.1–F1.4.3.
  - Procede de: REQ-SAL-22/26 y DTC-044/DTC-045/DTC-078. Los estados
    subordinados de eventos y conjuntos se introducen con sus consumidores.
  - Cambia: `src/execution.rs`, diagnósticos y pruebas de proceso.
  - Se comprueba con: códigos de salida, ausencia de destino en rechazo,
    manifiesto general `failed` tras fallo controlado y `complete` solo al éxito.

- ✅ F1.4.5. Retirar dependencias que hayan quedado sin consumidor runtime.
  - Depende de: F1.4.1–F1.4.4 y F1.1.2.
  - Procede de: DTC-043, DTC-048 y DTC-078.
  - Cambia: `Cargo.toml`, `Cargo.lock` y evaluación de dependencias; retira el
    motor JSON Schema runtime si ya solo pertenece a pruebas externas.
  - Se comprueba con: `cargo metadata --locked --format-version 1`,
    `cargo tree --locked -e features`, `cargo audit --no-fetch` y búsqueda de
    imports/usos del crate retirado.

> ✅ Pausa F1.4 aprobada: se inspeccionaron tres árboles reales —complete,
> failed e INCOMPLETE— y se confirmó que no contienen eventos, conjuntos ni
> defensas trasladadas. La séptima revisión no encontró defectos funcionales;
> su única deriva de trazabilidad se corrigió con una prueba CLI que demuestra
> copia TOML byte a byte, enlace, clase, ruta, tamaño y SHA-256 externo. F1.5 no
> se ha iniciado.

## 8. F1.5 — Referencia inicial

- ☐ F1.5.1. Conservar generación simulada y admisión/copia exacta de FASTA.
  - Depende de: F1.3.3 y F1.4.4.
  - Procede de: REQ-REF-01 a REQ-REF-07, REQ-REF-12 a
    REQ-REF-18, REQ-REF-22 y DTC-041. F1 aporta evidencia parcial al alfabeto
    universal de REQ-REF-23, que cierra en F9.7.
  - Cambia: adaptador FASTA y caso de uso, reutilizando código ya correcto.
  - Se comprueba con: fixtures de cabecera, duplicado, A/C/G/T/N, copia `cmp`,
    `samtools faidx` en pruebas y origen inalterado.

- ☐ F1.5.2. Enlazar la referencia autoritativa al manifiesto mínimo.
  - Depende de: F1.4.2 y F1.5.1.
  - Procede de: REQ-MAN-06/86 y REQ-SAL-01. REQ-MAN-72 y
    REQ-SAL-18 se cierran con eventos reales en F3.5.
  - Cambia: `EffectiveConfig` y `manifest.json`; conserva ID, origen, cabecera,
    descripción, longitud, N, ruta, tamaño y SHA-256.
  - Se comprueba con: lectura independiente del FASTA y comparación con ambos
    documentos sin revalidar relaciones futuras.

- ☐ F1.5.3. Repetir los escenarios pequeños de referencia y medir sin convertir
  el benchmark en contrato.
  - Depende de: F1.5.2.
  - Procede de: DTC-046/DTC-047 y REQ-NOF-01 a REQ-NOF-04.
  - Cambia: evidencias A/B/C y scripts de validación necesarios.
  - Se comprueba con: suite, release portable, FASTA/SHA/JSON externos y medidas
    documentadas como observación.

> Pausa F1.5: revisar una referencia simulada y otra suministrada completas,
> reproducibles y sin manifiestos de evento.

## 9. F1.6 — Referencia reducida

- ☐ F1.6.1. Reutilizar el cálculo decimal, coordenadas y selector de reducción
  tras aislarlos del grafo retirado.
  - Depende de: F1.5.3.
  - Procede de: DTC-050, DTC-056 y DTC-057.
  - Cambia: `src/exact_decimal.rs`, `src/coordinates.rs` y
    `src/reference_reduction.rs` solo donde dependan de alcance eliminado.
  - Se comprueba con: vectores exactos, bordes de coordenadas, preflight,
    uniformidad predeclarada y tamaño exacto.

- ☐ F1.6.2. Materializar FASTA y mapa de reducción sin
  `resolution-evidence`.
  - Depende de: F1.6.1.
  - Procede de: REQ-REF-24 a REQ-REF-36 y REQ-VAL-52/53.
  - Cambia: `src/reference_reduction_io.rs`, esquema del mapa y orquestación.
  - Se comprueba con: reconstructor independiente origen+mapa, `cmp`,
    `sha256sum`, `samtools faidx`, N, múltiples secuencias y residual.

- ☐ F1.6.3. Publicar la reducción como una unidad simple y coherente.
  - Depende de: F1.4.4 y F1.6.2.
  - Procede de: REQ-MAN-89/90, REQ-SAL-27 y DTC-078.
  - Cambia: `EffectiveConfig`, manifiesto general y catálogo simple para FASTA
    fuente, FASTA derivado y mapa.
  - Se comprueba con: fallo individual antes de publicación, ningún resultado
    parcial presentado como completo y enlaces/tamaños/SHA exactos.

- ☐ F1.6.4. Sustituir el expediente 83/83 por la asignación estricta vigente
  de 54 REQ F1 —4/1/6/22/21— y pruebas
  ejecutables honestas.
  - Depende de: F1.1.1–F1.6.3.
  - Procede de: DTC-078 y trazabilidad vigente.
  - Cambia: plan, matriz, evidencia y punto de reanudación; conserva la evidencia
    histórica sin usarla como cierre de código retirado.
  - Se comprueba con: recuento automático 4+1+6+22+21=54 por IDs exactos sin duplicados, cada
    TEST apunta a una función existente y ningún REQ trasladado figura como
    superado por F1.

- ☐ F1.6.5. Ejecutar cierre completo y dos revisiones independientes limpias.
  - Depende de: F1.6.4.
  - Procede de: estrategia acordada tras la auditoría sistemática.
  - Cambia: informe de revisión y acta de pausa; no cambia comportamiento.
  - Se comprueba con: fmt, Clippy, tests all-targets/all-features, doctest, doc,
    release, metadata/tree, audit y validadores TOML/JSON/FASTA/SHA; después
    dos revisiones consecutivas sin hallazgos graves ni bloqueantes.

> Pausa F1.6 y cierre de F1: no comenzar F2 hasta aprobar explícitamente el
> alcance mínimo, las evidencias vigentes y las dos revisiones limpias.

## 10. Matriz de salida

| Hito | Resultado | Evidencia principal | Riesgo protegido |
| --- | --- | --- | --- |
| F1.1 | Base y dependencias mínimas | Cargo/Clippy/test/release | Arrastrar infraestructura sin consumidor |
| F1.2 | Configuración efectiva real | Fixtures TOML y campos consumidos | Biología sintética y autoridad circular |
| F1.3 | Dos flujos RNG reales | Vectores y aislamiento | Repetición dependiente del orden accidental |
| F1.4 | Manifiesto y publicación mínimos | complete/failed/INCOMPLETE | Salida parcial o sobrescrita |
| F1.5 | Referencia inicial íntegra | `cmp`, `samtools`, SHA | Reescritura o identidad falsa del FASTA |
| F1.6 | Reducción reconstruible | origen+mapa, tamaño exacto | Coordenadas, sesgo o publicación parcial |
