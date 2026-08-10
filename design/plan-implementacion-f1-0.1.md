# Plan de implementación de F1 0.1

## 1. Cabecera y estado

- **Título:** plan detallado de F1, núcleo reproducible.
- **Versión del plan:** 0.1.
- **Fecha:** 8 de agosto de 2026.
- **Estado:** aprobado para ejecución; F0 está cerrada y ninguna casilla acredita todavía trabajo ejecutado.
- **Ámbito temporal:** únicamente F1.1–F1.5.
- **Fuentes normativas y de diseño:**
  [especificaciones-0.1.md](./especificaciones-0.1.md),
  [arquitectura-0.1.md](./arquitectura-0.1.md),
  [decisiones-tecnicas-0.1.md](./decisiones-tecnicas-0.1.md),
  [hoja-de-ruta-0.1.md](./hoja-de-ruta-0.1.md),
  [trazabilidad-requisitos-0.1.md](./trazabilidad-requisitos-0.1.md) y
  [registro-ampliaciones-futuras-0.1.md](./registro-ampliaciones-futuras-0.1.md).
- **Decisiones confirmadas aplicables:** DTC-034 a DTC-049. Este plan no
  convierte una propuesta DTP, una DTPD posterior ni una candidata de
  dependencia en decisión aprobada.

Este documento ordena trabajo ya acordado. F0 quedó aprobada el 8 de agosto de
2026 y su ejecución comienza en A1. El plan no autoriza a ampliar el alcance, a
marcar F1 como cerrada ni a hacer commits; deben respetarse sus pausas de
revisión obligatorias.

## 2. Resultado y alcance

F1 termina con una aplicación Rust mínima que admite una configuración TOML,
obtiene una configuración efectiva, prepara una ejecución identificada, usa una
semilla raíz y flujos aleatorios reconstruibles, conserva una referencia FASTA
simulada o suministrada y publica resultados con manifiestos, integridad y
estado coherentes. Dos ejecuciones bajo las mismas condiciones deben producir
el mismo resultado lógico; una entrada inválida no puede dejar un destino final
marcado como completo.

### 2.1. Incluido

- Proyecto Cargo, política de dependencias, calidad automática y construcción
  portable con Rust 1.97.1, edición 2024 y `Cargo.lock` versionado.
- Fronteras de dominio, aplicación y adaptadores suficientes para F1, sin
  filtrar tipos de crates a los contratos propios.
- Configuración general y envoltura de escenario TOML con clase y
  `schema_version = "1"`, copia exacta, `DeclaredConfig` y `EffectiveConfig`.
- Motor general de reglas base, selectores y dominancia parcial. Los requisitos
  asignados a F1 que pertenecen a capacidades biológicas posteriores se prueban
  aquí sobre el catálogo y el resolvedor, no mediante simuladores anticipados.
- Catálogo público de diagnósticos F1 y sus dos representaciones, texto humano
  y un único documento JSON por invocación.
- Catálogo normativo de ámbitos y vectores de subsemillas; catálogo normativo de
  clases, causas, prefijos y vectores de identificadores. Ambos preceden a todo
  consumidor.
- Identidad de entradas y artefactos mediante SHA-256, estados de ejecución y
  publicación Linux atómica no-clobber.
- `effective-config.json`, `manifest.json` general y manifiestos de los eventos
  explícitos, con esquemas JSON Schema Draft 2020-12.
- Referencia inicial pequeña simulada o FASTA suministrado sin comprimir,
  validación interna, copia exacta cuando proceda e índices solo si una operación
  de F1 los necesita.
- Validación independiente pertinente a F1 y escenarios sintéticos A, B y C de
  DTC-046, sin convertir sus medidas en objetivos de producto.

### 2.2. Excluido

- Individuos y haplotipos, variantes, VCF, verdad genotípica, lecturas,
  BAM/CRAM, mapeo y realismo técnico, que empiezan en F2–F4.
- Referencias derivadas, variantes estructurales, repeticiones y múltiples
  referencias complejas, que corresponden a F5.
- Construcción efectiva de poblaciones, evolución, recombinación, selección,
  familias y cruces, que corresponden a F6–F8. F1 solo deja preparado el
  resolvedor común que sus configuraciones necesitarán.
- Optimización, paralelismo, presupuestos de recursos, exportación
  autocontenida y publicación de una versión, que corresponden a F9. Los
  contrastes mínimos de `.tar.gz` y OCI exigidos por DTC-047 son pruebas de la
  base de distribución, no la implementación anticipada de F9.6/F9.7.
- Todas las ampliaciones DTF-001 a DTF-036: entre otras, migración de esquemas,
  reproducibilidad entre versiones, reanudación, ejecución distribuida,
  auditoría externa en runtime, calibración empírica y niveles configurables de
  trazabilidad.
- Elección de crates por conveniencia. Solo se incorpora cada dependencia tras
  su puerta de evaluación; `noodles-fasta` 0.66.0 superó la puerta E1.

## 3. Criterios de orden

1. Fijar primero el toolchain, las comprobaciones y las puertas de dependencias.
2. Crear los catálogos normativos de diagnósticos, subsemillas e IDs antes que
   cualquiera de sus consumidores.
3. Implementar reglas e invariantes en el dominio antes de TOML, JSON, FASTA o
   sistema de ficheros.
4. Separar siempre configuración declarada, validada y efectiva; ningún
   adaptador puede conceder precedencia por el orden o la jerarquía física.
5. Validar completamente antes de crear el plan de ejecución; publicar solo
   artefactos cerrados, validados y catalogados.
6. Probar cada contrato con una verdad independiente: vectores literales,
   casos manuales, `sha256sum`, `samtools faidx` y una implementación JSON
   Schema distinta del productor.
7. Integrar en cambios pequeños. Cada hito termina con una pausa obligatoria y
   no se inicia su consumidor mientras falle la evidencia de su productor.
8. Ejecutar habitualmente en Ubuntu 24.04 sobre WSL2 dentro del sistema de
   ficheros Linux/ext4; reservar Ubuntu 24.04 nativo y OCI para la referencia de
   publicación. Nunca activar `target-cpu=native`.

## 4. Prerrequisitos y puertas

### 4.1. Condiciones previas

- F0 debe aprobar la arquitectura y este plan.
- El trabajo puede realizarse en Ubuntu 24.04 nativo o en Ubuntu 24.04 sobre
  WSL2; cuando use WSL2 debe trabajar sobre el sistema de ficheros Linux/ext4,
  no en `/mnt/c`. Las evidencias registran plataforma y sistema de ficheros.
- Toda nueva dependencia debe tener propósito, REQ atendidos, versión exacta,
  licencia, mantenimiento, seguridad, estabilidad, rendimiento, características
  activadas, árbol final, código nativo, scripts de construcción y coste de
  distribución revisados conforme a DTC-043.
- La selección de una implementación independiente de JSON Schema y de una
  herramienta de inspección/arranque OCI debe quedar registrada con versión y
  comando exactos antes de usarse en aceptación. El plan no inventa esas
  herramientas pendientes.
- Una actualización posterior de Rust, del lock, de un crate, de un catálogo o
  de un algoritmo aleatorio exige un cambio aislado y repetir vectores,
  validadores, reproducibilidad y escenarios de rendimiento.

### 4.2. Comandos disponibles y comandos posteriores

Ahora, antes de crear Cargo, solo son comprobables la documentación y el árbol
de trabajo, por ejemplo:

```bash
git diff --no-index --check /dev/null design/plan-implementacion-f1-0.1.md
rg -n '^- ☐ [A-E][0-9]+\.' design/plan-implementacion-f1-0.1.md
rg -n 'DTC-049|noodles-fasta|Rust 1\.97\.1|edition = "2024"' design/plan-implementacion-f1-0.1.md
```

Mientras el plan sea un fichero nuevo, la primera orden termina con código 1
porque compara dos contenidos distintos; la comprobación es satisfactoria si no
imprime diagnósticos de espacios o líneas inválidas.

Los siguientes comandos son criterios de pasos futuros. No se presentarán como
ejecutados hasta que A1 haya creado el proyecto y A2 haya producido el lock:

```bash
rustc --version
cargo --version
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --format-version 1
cargo tree --locked -e features
```

Los validadores `samtools faidx`, `sha256sum`, listado/extracción estándar de
`.tar.gz`, JSON Schema independiente e inspección/arranque OCI pertenecen a
pruebas, CI y aceptación, no al recorrido normal ni a la instalación de la
persona usuaria.

## 5. F1.1 — Base del proyecto

**Resultado del hito:** proyecto Rust reproducible y portable con puertas de
calidad, dependencias y validadores, además de los catálogos que deben preceder
a la implementación funcional.

- ✅ A1. Crear la raíz Cargo y las fronteras mínimas de dominio, aplicación y adaptadores.
  - Depende de: aprobación de F0 y de este plan.
  - Procede de: F1.1; arquitectura §§5–7; DTC-043 y DTC-048.
  - Cambia: `Cargo.toml`, `rust-toolchain.toml` y estructura mínima bajo `src/`; fija Rust 1.97.1, perfil `minimal`, `rustfmt`, `clippy`, target `x86_64-unknown-linux-gnu` y edición 2024, sin dependencias no evaluadas.
  - Se comprueba con: `rustc --version`, `cargo --version` y revisión literal de `rust-toolchain.toml` y `Cargo.toml`; no debe aparecer `target-cpu=native`.

- ✅ A2. Generar y versionar la resolución reproducible de dependencias de la aplicación.
  - Depende de: A1 y aceptación individual de cualquier dependencia inicial.
  - Procede de: DTC-043 y DTC-048.
  - Cambia: `Cargo.lock` y el registro de evaluación de cada dependencia realmente incorporada; no añade capacidades anticipadas.
  - Se comprueba con: `cargo metadata --locked --format-version 1` y `cargo tree --locked -e features`; ambos deben corresponder al lock versionado.

- ✅ A3. Establecer las puertas automáticas de formato, análisis estático, pruebas y construcción portable.
  - Depende de: A1–A2.
  - Procede de: F1.1; DTC-042, DTC-043 y DTC-048.
  - Cambia: configuración de automatización y documentación de comandos, separando pruebas rápidas, integración, validadores independientes y observaciones manuales.
  - Se comprueba con: `cargo fmt --all -- --check`, `cargo clippy --workspace --all-targets --locked -- -D warnings`, `cargo test --workspace --all-targets --locked`, `cargo test --workspace --doc --locked` y `cargo build --release --locked --target x86_64-unknown-linux-gnu`.

- ✅ A4. Crear el catálogo y el esquema públicos versionados de diagnósticos necesarios en F1 antes de implementar sus emisores.
  - Depende de: A1.
  - Procede de: DTC-036 y DTC-044; REQ-CONF/REF/VAL/MAN/SAL asignados a F1.
  - Cambia: artefacto normativo de códigos `GPS-{CONF|INPUT|FMT|DOMAIN|RES|TOOL|IO|INTERNAL}-NNN`, significado inmutable, contexto exigido y salida 0/1/2/3/4/5/70; además define el esquema Draft 2020-12 del documento único con `diagnostic_kind: "gps_diagnostics"`, `schema_version: "1"` y la lista `diagnostics`.
  - Se comprueba con: prueba de catálogo que rechace códigos repetidos, categorías desconocidas y significados vacíos; fixtures humano/JSON que conserven el mismo diagnóstico y revisión que asigne al menos un código a cada rechazo previsto por B, C, D y E.

- ✅ A5. Crear el catálogo normativo de ámbitos de subsemilla y sus vectores literales antes del generador.
  - Depende de: A1.
  - Procede de: DTC-037 y DTC-038; REQ-CONF-111 a REQ-CONF-115 y REQ-REP-22/23.
  - Cambia: artefacto versionado con ámbitos, nombres ASCII de campos, orden, valores UTF-8, versión de derivación, material binario esperado y subsemilla esperada; cubre al menos población inicial, generación/proceso biológico, muestreo y evento técnico sin implementar esos consumidores posteriores.
  - Se comprueba con: revisión byte a byte de vectores publicados que incluya longitudes y enteros `u32` big-endian y demuestre que no usa ruta, reloj, hilo, posición accidental ni orden de iteración.

- ✅ A6. Crear el catálogo normativo de IDs y sus vectores literales antes del registro de identidades.
  - Depende de: A1 y A5 para reutilizar únicamente la codificación canónica confirmada.
  - Procede de: DTC-039; REQ-CONF-38/101, REQ-IND-25/26 y REQ-EVO-63.
  - Cambia: artefacto versionado de clases, prefijos ASCII, ámbitos, causas mínimas, orden de campos, material binario, primeros 128 bits y Base32 RFC 4648 mayúscula sin `=`; incluye identidad de ejecución causada por raíz y digest de configuración.
  - Se comprueba con: vectores literales revisados independientemente y casos de separación entre `ReferenceId`, `IndividualId`, `EventId` e ID de ejecución; un cambio sin nueva versión debe fallar la prueba.

- ✅ A7. Aprobar las herramientas independientes y el formato de sus evidencias antes de introducirlas en CI.
  - Depende de: A3.
  - Procede de: DTC-042, DTC-043 y DTC-047.
  - Cambia: registro versionado con versión y comando literal de `samtools`, `sha256sum`, listado/extracción de `.tar.gz`, implementación JSON Schema ajena al productor y herramienta de inspección/arranque OCI; fija por digest la base y receta OCI cuando se creen.
  - Se comprueba con: acta de evaluación completa y una prueba deliberadamente corrupta por clase; no puede cerrarse A7 con nombres, versiones o comandos pendientes.

### Pausa obligatoria de revisión de F1.1

No comenzar consumidores de configuración, RNG o IDs hasta verificar A4–A6.
La revisión debe conservar la salida de las cinco órdenes Cargo de A3, el árbol
de características, el contrato Ubuntu/WSL2 y las evaluaciones de dependencias.
Si falla una puerta, se corrige F1.1 y se repite la pausa.

**Evidencia ejecutada en la pausa de F1.1 (8 de agosto de 2026):** las cinco
órdenes Cargo de A3, `cargo metadata --locked --format-version 1`,
`cargo tree --locked -e features`, las siete pruebas de catálogos y
`scripts/check-independent-validators.sh` terminaron correctamente en Ubuntu
24.04 sobre WSL2/ext4. El grafo Cargo no contiene dependencias externas. Las
evaluaciones están en `evaluaciones-dependencias-f1.md` y
`evaluacion-validadores-independientes-f1.md`; la salida negativa de los cinco
validadores se conserva en `../evidence/f1.1/validators-corrupt-20260808.md`.
La pausa queda alcanzada y no autoriza comenzar F1.2 sin revisión del usuario.

**Aprobación de la pausa (9 de agosto de 2026):** F1.1 queda aprobada tras la
revisión independiente, la corrección de los dos hallazgos y la repetición de
las puertas Cargo y de los cinco validadores. El commit `a95e614` conserva su
implementación y evidencias; queda autorizado comenzar F1.2.

## 6. F1.2 — Configuración validada

**Resultado del hito:** admisión TOML y motor de resolución comprobados antes de
simular, con configuración efectiva inmutable y procedencia de la regla
ganadora. Los consumidores biológicos posteriores no se implementan aquí.

- ✅ B1. Fijar ejemplos mínimos y esquemas de autoría para configuración general y envolturas de escenario.
  - Depende de: A1, A4 y aprobación de los campos derivados de los REQ asignados a F1.
  - Procede de: DTC-034 y DTC-036; REQ-CONF-01/02/03/04/07/09/10/101.
  - Cambia: esquemas y fixtures TOML con `config_kind`, `schema_version = "1"`, vía, referencia simulada/suministrada, escenarios como listas identificadas y cero/uno/varios eventos explícitos.
  - Se comprueba con: revisión contra ejemplos válidos e inválidos que cubran clase o versión ausente, desconocida o intercambiada sin interpretación parcial.

- ✅ B2. Evaluar y aprobar las dependencias exactas para TOML y configuración antes de añadirlas al lock.
  - Depende de: A2 y B1.
  - Procede de: DTC-034, DTC-043 y DTC-048.
  - Cambia: registro de evaluación, versiones y características mínimas; solo después, `Cargo.toml` y `Cargo.lock`.
  - Se comprueba con: `cargo tree --locked -e features` más fixtures que demuestren detección localizable de campos y claves duplicados y desconocidos; una biblioteca que no conserve la evidencia necesaria se rechaza.

- ✅ B3. Implementar la frontera de versión y clase antes del parseo semántico.
  - Depende de: B1–B2 y catálogo A4.
  - Procede de: DTC-034, DTC-036 y DTC-044.
  - Cambia: adaptador TOML de admisión que conserva bytes exactos, identifica ruta y ubicación y rechaza clase/versión inválida antes de construir el dominio.
  - Se comprueba con: `cargo test --locked` sobre versión ausente, `"0"`, `"2"`, no textual, clase cruzada y TOML malformado; cada caso obtiene código, contexto y salida 2.

- ✅ B4. Construir `DeclaredConfig` sin conceder autoridad a la jerarquía u orden físico del TOML.
  - Depende de: B3.
  - Procede de: DTC-034; REQ-CONF-01/09/10/13/14.
  - Cambia: tipos declarados con valor, fuente, ruta y ubicación; fases, eventos y cruces conservan IDs explícitos y orden solo donde sea normativo.
  - Se comprueba con: `cargo test --locked` comparando TOML reordenados semánticamente equivalentes y casos con IDs ausentes o duplicados; la copia original debe coincidir byte a byte.

- ✅ B5. Materializar el catálogo normativo de familias de parámetros y selectores permitidos.
  - Depende de: B4 y A6; no depende de consumidores de F2–F8.
  - Procede de: REQ-CONF-29/117. La auditoría de cierre conserva aquí la infraestructura parcial de REQ-CONF-17/18/20/22/43/83/84/116/120/121/122, REQ-IND-13/15/21 y REQ-REP-02, pero reasigna su cierre al hito que prueba todos sus consumidores.
  - Cambia: tipos de eje y metadatos de cada familia para individuo preexistente, población/subpoblación, fase/generación, secuencia/región/locus, grupo/bloque y evento técnico solo donde esté permitido.
  - Se comprueba con: tabla exhaustiva de aceptación/rechazo independiente del resolvedor; incluye evento técnico aplicado a parámetro biológico, individuo futuro, ámbito desconocido y extremos de la jerarquía espacial.

- ✅ B6. Implementar reglas base y coincidencia por todos los selectores.
  - Depende de: B5.
  - Procede de: REQ-CONF-29/117. REQ-CONF-16/17/20/22/24/116/121 conserva aquí infraestructura reutilizable, pero se traza al hito donde puede demostrarse íntegramente.
  - Cambia: tipos de regla, regla base obligatoria donde proceda, predicado de coincidencia y valor aplicable sin prioridades entre ejes distintos.
  - Se comprueba con: `cargo test --locked` mediante una tabla manual de casos global, población, fase, secuencia, región y locus; cada resultado se compara con el ganador esperado literal.

- ✅ B7. Implementar la dominancia parcial y exigir una única máxima.
  - Depende de: B6.
  - Procede de: REQ-CONF-29/31/117/118 y REQ-VAL-51; la base parcial de REQ-CONF-116 no constituye su cierre completo.
  - Cambia: comparación parcial entre reglas y resolución que devuelve una máxima única o un conflicto de máximas incomparables, sin desempate por orden.
  - Se comprueba con: `cargo test --locked` sobre un retículo pequeño enumerado manualmente; permutar las reglas no cambia el ganador y dos máximas incomparables producen el diagnóstico previsto antes de simular.

- ✅ B8. Rechazar duplicados exactos, solapes espaciales y selectores o referencias inválidos.
  - Depende de: B5–B7.
  - Procede de: REQ-CONF-30/38/65/119 y REQ-VAL-51. REQ-CONF-18/120/122 conserva validación parcial aquí y queda cerrado solo tras probar sus consumidores reales en F9.7.
  - Cambia: validadores previos de unicidad, pertenencia a catálogos, intervalos, referencias de selector y solape para el mismo parámetro y selectores no espaciales.
  - Se comprueba con: `cargo test --locked` sobre duplicado exacto con valores iguales/diferentes, bordes de regiones cerradas, solape parcial/total, adyacencia válida, ID inexistente y referencia de selector de categoría incorrecta.

- ✅ B9. Validar coherencia mínima de referencia, vía, escenarios, eventos y parámetros declarados sin ejecutar biología futura.
  - Depende de: B4–B8.
  - Procede de: REQ-CONF-04/05/07/09/10/65/101 y REQ-REF-22. REQ-CONF-32 y REQ-CONF-43 se reasignan respectivamente a F6.2 y F7.6, donde existe el comportamiento consumidor.
  - Cambia: validación cruzada y cálculo comprobado del tamaño total; clasifica cada secuencia como cromosoma o `scaffold`, exige escenario para las vías que lo requieren y conserva cero eventos sin crear uno implícito.
  - Se comprueba con: `cargo test --locked` sobre suma conocida y desbordamiento, vía/escenario incompatibles, evento duplicado/implícito, ploidía efectiva ausente y formas de frecuencia incompatibles; no se crea directorio final.

- ✅ B10. Implementar la validación y normalización numérica reproducible fijada por DTC-049.
  - Depende de: B6–B9 y DTC-049.
  - Procede de: DTC-049; REQ-CONF-106/107 y ACE-105; afecta a REQ-MAN-81. La infraestructura general de REQ-VAL-40 queda aquí, pero su cierre se reasigna a F6.2 para incluir la tolerancia 0,5/C del modo controlado.
  - Cambia: valores y tolerancia como `f64` finitos; probabilidades y frecuencias en `[0,1]`; tolerancia absoluta configurable en `[0, 1e-6]` y `1e-9` predeterminada; suma Neumaier en orden canónico estable por identidad y categoría, sin reducción paralela; cálculo de suma y delta, borde inclusivo `abs(suma - 1) <= tolerancia`, normalización y recomputación de suma y delta o rechazo previo.
  - Se comprueba con: vectores literales justo dentro, en el borde y justo fuera; orden de entrada permutado; valores/tolerancias no finitos o fuera de rango; suma cero/no positiva/no normalizable; bits `f64` esperados tras parsear los decimales TOML; los esperados no se calculan con el normalizador probado.

- ✅ B11. Producir una `EffectiveConfig` inmutable con valores, regla ganadora, selectores y procedencia.
  - Depende de: B7–B10.
  - Procede de: DTC-034 y DTC-035; REQ-CONF-123 y REQ-MAN-05/10/88.
  - Cambia: modelo de configuración efectiva separado del declarado; no contiene referencias mutables ni reinterpreta el TOML al consumirlo.
  - Se comprueba con: `cargo test --locked` sobre instantáneas semánticas de configuración efectiva y prueba que modificar o reordenar la entrada después de admisión no altera el objeto ya resuelto.

- ✅ B12. Auditar las 36 pruebas trazables inicialmente asignadas a F1.2 y cerrar solo las demostrables en este hito.
  - Depende de: B1–B11.
  - Procede de: los 36 REQ inicialmente asignados a F1.2 en la matriz, DTC-044 y DTC-049.
  - Cambia: casos `TEST-REQ-*`, pruebas de integración CLI y matriz de auditoría con 16 casos superados y 20 reasignados al hito que permite demostrar el REQ completo; código, contexto y ausencia de publicación final.
  - Se comprueba con: `cargo test --workspace --all-targets --locked`, seguido de ejecución de cada fixture inválido en modo humano y JSON; JSON es un único documento y ambos renderizados conservan el mismo diagnóstico estructurado.

### Pausa obligatoria de revisión de F1.2

Revisar por separado sintaxis, semántica, dominancia y diagnósticos. La revisión
debe inspeccionar casos de máximas incomparables y reglas permutadas, no solo
casos felices. Los vectores de DTC-049 deben demostrar el orden canónico, la
suma Neumaier, el borde inclusivo y la recuperación exacta del `f64` antes de
cerrar F1.2.

**Aprobación de la pausa (9 de agosto de 2026):** F1.2 queda aprobada tras la
revisión independiente de las correcciones H-F12-01 a H-F12-04. Las puertas de
formato, Clippy, pruebas, documentación, construcción release, grafo bloqueado,
JSON Schema y `cargo audit` terminaron sin fallos. La matriz audita los 36
`TEST-REQ-*` iniciales, con 16 superados en F1.2 y 20 reasignados honestamente a
su hito consumidor o de consolidación. Los diagnósticos JSON de 17 fixtures
persistentes y 12 casos de dominio cumplen el esquema independiente. Evidencias
en `../evidence/f1.2/`. Queda autorizado comenzar C1/F1.3, pero no se inicia en
este commit.

## 7. F1.3 — Aleatoriedad reproducible

**Resultado del hito:** raíz y subsemillas de 256 bits, generadores ChaCha12
locales y comparadores de reproducibilidad demostrados mediante vectores
normativos y casos resistentes al reordenamiento.

- ✅ C1. Evaluar y aprobar las dependencias exactas de CSPRNG, BLAKE3 y ChaCha12 antes de añadirlas.
  - Depende de: A2, catálogo A5 y DTC-043.
  - Procede de: DTC-037, DTC-038 y DTC-048.
  - Cambia: evaluaciones, versiones y características mínimas; después de aprobación, `Cargo.toml` y `Cargo.lock`, manteniendo contratos propios.
  - Se comprueba con: `cargo tree --locked -e features`, revisión de código nativo/build scripts/licencias y ejecución de los vectores A5; no se admite `StdRng`.

- ✅ C2. Implementar el tipo de semilla raíz y su representación canónica.
  - Depende de: C1 y A4.
  - Procede de: DTC-037; REQ-CONF-26/27.
  - Cambia: tipo de 32 bytes, parseo y formato de exactamente 64 hexadecimales minúsculos, más origen declarado o generado.
  - Se comprueba con: `cargo test --locked` sobre longitud 63/64/65, mayúsculas, caracteres no hexadecimales, ida/vuelta de bytes y valor no truncado.

- ✅ C3. Obtener una raíz ausente desde el CSPRNG del sistema antes de simular y registrar su origen.
  - Depende de: C2.
  - Procede de: DTC-037 y DTC-044.
  - Cambia: proveedor inyectable de entropía del sistema y fallo tipado; la raíz generada entra en `EffectiveConfig`. Su serialización y publicación en manifiesto quedan en D7/E8/F1.5.
  - Se comprueba con: doble determinista de 32 bytes y doble que falla en la frontera de aplicación con `GPS-RES-001` y salida 3; el CLI rechaza una semilla declarada inválida en modos humano/JSON y no crea el destino solicitado. No se afirma una prueba CLI del fallo CSPRNG.

- ✅ C4. Implementar la codificación binaria canónica de materiales de derivación.
  - Depende de: catálogo A5 y C2.
  - Procede de: DTC-037.
  - Cambia: codificador de raíz cruda, versión y número de componentes `u32` big-endian y pares nombre ASCII/valor UTF-8 con longitud `u32` big-endian.
  - Se comprueba con: comparación byte a byte contra todos los materiales literales de A5, incluidos vacío, Unicode, longitudes y cambios de orden.

- ✅ C5. Derivar subsemillas BLAKE3 y crear generadores ChaCha12 locales por ámbito.
  - Depende de: C1, C4 y A5.
  - Procede de: DTC-037; REQ-CONF-111/112/113 y REQ-REP-22.
  - Cambia: derivador con contexto exacto `genome_population_simulator/v1/subseed`, versión explícita y fábrica que entrega estado local, nunca global o compartido.
  - Se comprueba con: vectores A5 y pruebas que añaden, omiten y reordenan consumidores ajenos sin cambiar la secuencia del ámbito observado.

- ✅ C6. Preparar metadatos RNG para que los hitos consumidores puedan reconstruir sus flujos.
  - Depende de: C3 y C5.
  - Procede de: DTC-037 y DTC-038; prepara parte de REQ-CONF-114/115 y REQ-REP-22/23 sin afirmar todavía manifiestos ni resultados reales.
  - Cambia: contrato propio de raíz/origen, algoritmo, versión de derivación, catálogo y subsemillas de evento dentro de `EffectiveConfig`; los ámbitos internos se reconstruyen desde IDs versionados.
  - Se comprueba con: prueba que reconstruye una subsemilla de evento desde la configuración efectiva y sus IDs; el manifiesto, las entradas publicadas y los consumidores reales quedan en F1.4–F9.

- ✅ C7. Implementar la infraestructura preparatoria de cuatro comparadores separados de reproducibilidad.
  - Depende de: C6 y contratos de resultado DTC-038.
  - Procede de: REQ-CONF-111 a REQ-CONF-115 y REQ-VAL-44.
  - Cambia: comparadores preparatorios para decisiones biológicas, contenido técnico, equivalencia semántica y bytes; cada nivel informa `garantizado`, `no garantizado` o `no aplicable` y sus condiciones, pero ningún estado constituye aún una afirmación sobre artefactos no implementados.
  - Se comprueba con: valores sintéticos donde el contenido lógico coincide pero cambian bytes, donde cambia una subsemilla técnica o una versión relevante; ningún nivel implica automáticamente el siguiente y dependencias, catálogo y derivación forman parte de la base RNG.

- ✅ C8. Auditar los diez REQ y cerrar únicamente la infraestructura demostrable de F1.3.
  - Depende de: C2–C7 y una configuración mínima validada de B.
  - Procede de: auditoría literal de los diez REQ inicialmente asignados a F1.3 y ACE-23.
  - Cambia: matriz con 1 requisito superado, 3 trasladados a F1.5 y 6 reasignados a F4.6/F9.7; conserva nueve pruebas parciales sin presentarlas como simulación, manifiesto, lecturas o artefactos reales.
  - Se comprueba con: `cargo test --workspace --all-targets --locked`; vectores RNG repetidos, condiciones negativas y prueba completa de que una versión distinta no recibe garantía, distinguiendo esa infraestructura de los consumidores pendientes.

### Pausa obligatoria de revisión de F1.3

No integrar RNG en generación de referencia hasta que todos los vectores A5
pasen. Revisar explícitamente que ningún ID, flujo o resultado dependa de hilo,
orden de mapa, ruta, reloj o consumo de otro ámbito. Conservar el lock y los
vectores como evidencia conjunta.

**Aprobación de la pausa (9 de agosto de 2026):** F1.3 queda aprobada tras la
revisión independiente y la corrección de H-F13-01 a H-F13-03 y
H-F13-02-R1. Las puertas de formato, Clippy, 83 pruebas, doctests,
documentación, construcción release, grafo bloqueado, JSON Schema y
`cargo audit` terminaron sin fallos. La auditoría conserva los diez casos: 1
está superado y 9 quedan reasignados con evidencia parcial y destino explícitos.
Los cuatro vectores A5, el aislamiento por ámbito y el resto de evidencias
constan en `../evidence/f1.3/`. Queda autorizado comenzar D1/F1.4, pero no se
inicia en este cambio.

## 8. F1.4 — Identidad y resultados

**Resultado del hito:** identidades tipadas, configuración efectiva y
manifiestos JSON íntegros dentro de una publicación atómica que distingue
rechazo previo, fallo controlado e interrupción.

- ✅ D1. Evaluar y aprobar dependencias exactas para JSON, JSON Schema, SHA-256 y publicación no-clobber.
  - Depende de: A2, A7 y los contratos DTC-035/040/045.
  - Procede de: DTC-035, DTC-040, DTC-043, DTC-045, DTC-047 y DTC-048.
  - Cambia: evaluaciones y, solo tras aprobación, `Cargo.toml` y `Cargo.lock`; separa productor JSON de validador independiente y oculta la llamada Linux tras adaptador propio.
  - Se comprueba con: `cargo tree --locked -e features`, auditoría del grafo y pruebas de capacidad para RFC 8259, Draft 2020-12, hash incremental y no-clobber real sin patrón check-then-rename.

- ✅ D2. Implementar tipos de ID aportado y etiqueta humana sin mezclar categorías.
  - Depende de: A6 y A4.
  - Procede de: DTC-039; REQ-CONF-38, REQ-IND-25/26.
  - Cambia: nuevos tipos por entidad, validación conceptual `[A-Za-z0-9][A-Za-z0-9._-]{0,63}`, etiqueta Unicode opcional y registro central por ámbito.
  - Se comprueba con: `cargo test --locked` sobre longitudes 1/64/65, mayúsculas, espacios, rutas, Unicode en ID/etiqueta, duplicado y colisión; un bloque documental `compile_fail` demuestra el cruce de tipos y se ejecuta con `cargo test --workspace --doc --locked`, sin añadir una dependencia de pruebas solo para esta garantía.

- ✅ D3. Implementar IDs automáticos exactamente desde el catálogo y vectores A6.
  - Depende de: C4, D1–D2 y A6.
  - Procede de: DTC-039.
  - Cambia: derivación BLAKE3 con contexto `genome_population_simulator/v1/id`, primeros 128 bits, Base32 RFC 4648 mayúscula sin relleno y prefijo por clase.
  - Se comprueba con: todos los vectores A6, longitud/alfabeto/prefijo exactos, causas de ejecución raíz+digest de configuración y rechazo de colisión inyectada.

- ✅ D4. Implementar SHA-256 incremental sobre bytes exactos y catálogo de artefactos.
  - Depende de: D1.
  - Procede de: DTC-040; REQ-MAN-08 y REQ-SAL-03.
  - Cambia: adaptador que entrega algoritmo literal `sha256`, 64 hex minúsculos, tamaño y ruta relativa por fichero; nunca sustituye contenido por metadatos o muestreo.
  - Se comprueba con: fixtures vacíos/binarios/textuales comparados con `sha256sum`; modificar un byte debe cambiar el digest y truncar debe fallar tamaño/integridad.

- ✅ D5. Definir los esquemas JSON de resultados de F1 antes de sus escritores y reutilizar el esquema diagnóstico de A4.
  - Depende de: B11, C6, D1–D4 y catálogo A4.
  - Procede de: DTC-035, DTC-036, DTC-038 y DTC-044; REQ-MAN-01 a REQ-MAN-05, REQ-MAN-08/10/32/69/72/81/83/88.
  - Cambia: esquemas Draft 2020-12 separados para `effective-config.json`, manifiesto general y manifiesto de evento, con definiciones comunes, discriminadores y `schema_version: "1"`; integra sin redefinir el contrato diagnóstico ya fijado en A4.
  - Se comprueba con: revisión campo por campo contra los REQ y fixtures válidos, truncados, corruptos, con campo desconocido y clase/versión errónea validados por la implementación independiente aprobada en A7.

- ✅ D6. Serializar `effective-config.json` con orden estable donde la semántica sea de conjunto.
  - Depende de: B11 y D5; la sección de normalización depende de B10.
  - Procede de: DTC-035 y DTC-036; REQ-CONF-123 y REQ-MAN-05/10/88.
  - Cambia: escritor legible que registra valor efectivo, regla ganadora, selectores y procedencia y enlaza la copia TOML exacta sin duplicarla.
  - Se comprueba con: validación JSON Schema independiente, permutación de entradas equivalentes y comparación semántica; no se afirma identidad binaria general.

- ✅ D7. Implementar manifiesto general y manifiestos solo para eventos explícitos.
  - Depende de: D3–D6.
  - Procede de: DTC-035 y DTC-036; REQ-CONF-101, REQ-EVO-63, REQ-MAN-01/02/03/04/08/32/69/72/83 y REQ-SAL-04/17/18.
  - Cambia: escritores con clase, versión, vía, configuración/escenarios, eventos, estados, plataforma, toolchain, dependencias y enlaces de artefacto; cero eventos no crea directorio o manifiesto de evento implícito.
  - Se comprueba con: casos cero/uno/dos eventos, IDs permutados y enlaces relativos; cada JSON supera el validador independiente y todos los enlaces resuelven dentro de la raíz.

- ✅ D8. Implementar el modelo de estado y el catálogo de unidades publicables.
  - Depende de: D4–D7.
  - Procede de: DTC-045; REQ-MAN-83 y REQ-SAL-22/26.
  - Cambia: estados funcionales `incompleto`, `completo` y `fallido`, más el marcador físico `INCOMPLETE` del provisional; artefactos y conjuntos independientes se incorporan solo tras cierre, validación, tamaño y digest.
  - Se comprueba con: máquina de estados probada exhaustivamente; ninguna transición desde artefacto abierto/inválido llega a completo y un fallo controlado enumera unidades válidas y diagnóstico.

- ✅ D9. Implementar directorio provisional hermano y temporales por artefacto.
  - Depende de: D8.
  - Procede de: DTC-045; REQ-SAL-22/26.
  - Cambia: ensamblador que trabaja en el mismo sistema de ficheros, marca `INCOMPLETE`, valida antes de catalogar y no crea destino final durante rechazo de admisión.
  - Se comprueba con: pruebas de integración para rechazo previo, fallo durante escritura y temporal truncado; se inspeccionan raíz provisional, catálogo y ausencia de destino final indebido.

- ✅ D10. Implementar publicación Linux atómica no-clobber detrás de un adaptador propio.
  - Depende de: D1 y D9.
  - Procede de: DTC-045 y DTC-042; REQ-SAL-22/26.
  - Cambia: operación equivalente a `renameat2(..., RENAME_NOREPLACE)` que rechaza `EEXIST` o sistema de ficheros sin garantía, sin fallback con comprobación previa.
  - Se comprueba con: carrera de dos procesos por el mismo destino, destino preexistente y doble que simula operación no soportada; exactamente uno publica y nunca se mezclan o sobrescriben bytes.

- ✅ D11. Completar la evidencia de manifiesto de las normalizaciones fijadas por DTC-049.
  - Depende de: B10 y D5–D7.
  - Procede de: DTC-049; REQ-MAN-81, REQ-CONF-106/107 y REQ-VAL-40.
  - Cambia: registro de tolerancia efectiva, orden canónico, algoritmo y versión, suma, delta y valores antes y después o causa de rechazo, sin reutilizar tolerancias con otra finalidad.
  - Se comprueba con: los mismos vectores independientes de B10 a través de JSON Schema y comparación semántica del manifiesto; la serialización recupera el mismo `f64` y ningún predeterminado queda implícito.

- ✅ D12. Probar rechazo, fallo controlado e interrupción no controlable como resultados distintos.
  - Depende de: D7–D10.
  - Procede de: DTC-045 y DTC-044; REQ-MAN-83 y REQ-SAL-26.
  - Cambia: pruebas de proceso para salida 2 sin destino, salida de fallo controlado con manifiesto final `fallido` y señal no controlable con provisional `INCOMPLETE`.
  - Se comprueba con: ejecución en WSL2/ext4 y CI Ubuntu 24.04; se inyectan fallos antes/después de catalogar un artefacto y se verifica que un SIGKILL no publique completo.

- ✅ D13. Auditar las 25 pruebas trazables de F1.4 con validadores independientes y cerrar solo las demostrables.
  - Depende de: D1–D12.
  - Procede de: los 25 REQ asignados a F1.4, ACE-22/25/34/93/94/105, DTC-047 y DTC-049.
  - Cambia: casos `TEST-REQ-*`, verificador de enlaces, pruebas de hash/no-clobber y matriz de evidencias.
  - Se comprueba con: `cargo test --workspace --all-targets --locked`, implementación JSON Schema independiente y `sha256sum`; corruptos, truncados y campos desconocidos deben fallar por un oráculo distinto del productor.

**Estado local de D12–D13 (10 de agosto de 2026):** la implementación y las
pruebas de proceso están completas en Ubuntu sobre WSL2/ext4: rechazo con
salida 2 y sin destino, fallo controlado con salida 1 y manifiesto `failed`,
SIGKILL con provisional `INCOMPLETE`, carrera de dos procesos, JSON Schema
independiente y `sha256sum`. La auditoría literal conserva los 25 casos con 17
superados y 8 reasignados; no se cierran D12–D13 ni la fase hasta repetir la
matriz de publicación en Ubuntu 24.04 nativo. Evidencia en `../evidence/f1.4/`.
La revisión independiente posterior detectó H-F14-R1/R2: ambos quedaron
corregidos con regresiones del marcador durante publicación y contratos JSON
anidados cerrados, sin alterar que D12–D13 siguen pendientes de la repetición
nativa. La pasada de coherencia H-F14-R3 corrigió los únicos residuos de
recuento de E13 y la pausa F1.5 a 27/61/31. La segunda revisión de esquemas
detectó H-F14-R2-R1: la regresión demostró que el contrato cerrado todavía
aceptaba contradicciones de estado, niveles y dependencias. El contrato Draft
2020-12 exige ahora los cuatro niveles y las siete dependencias exactamente una
vez, además de la coherencia entre estado general, configuración efectiva,
diagnósticos y unidades y la unidad completa de evento. En ese punto D12–D13
continuaban abiertos hasta la repetición nativa registrada a continuación.

**Aprobación de la pausa (10 de agosto de 2026):** F1.4 queda aprobada tras
la revisión independiente, la corrección de los diez hallazgos y la repetición
nativa en GitHub Actions sobre `ubuntu-24.04`. La ejecución
[`31372725778`](https://github.com/ximocanizares/genome_population_simulator/actions/runs/31372725778),
para el commit `5d12e58965561cd63eec76020c56ce7cc9d6cf8c`, terminó con éxito:
formato, Clippy, 109 pruebas, doctest y construcción release. El job
[`93405082155`](https://github.com/ximocanizares/genome_population_simulator/actions/runs/31372725778/job/93405082155)
demuestra expresamente rechazo previo, fallo controlado, carrera no-clobber y
SIGKILL con `INCOMPLETE`. D12–D13 quedan cerrados y se autoriza comenzar
E1/F1.5, sin declarar todavía F1 completa.

### Pausa obligatoria de revisión de F1.4

Inspeccionar un resultado completo, uno fallido y un provisional interrumpido.
Repetir la carrera no-clobber en WSL2/ext4 y Ubuntu nativo. Ninguna revisión
puede cerrar F1.4 mientras falte un enlace, un digest, la procedencia de un valor
efectivo o la evidencia completa de REQ-MAN-81 conforme a DTC-049.

## 9. F1.5 — Referencia mínima

**Resultado del hito:** una ejecución usa y conserva una referencia inicial
simulada o suministrada con identidad verificable, sin reescritura silenciosa
del original y sin dependencia de validadores externos en runtime.

- ✅ E1. Revalidar y decidir la aceptación final de la candidata `noodles-fasta` 0.66.0 antes de incorporarla.
  - Depende de: A2, A7 y DTC-041/DTC-043/DTC-048.
  - Procede de: DTC-041, que exigía evaluarla antes de fijarla.
  - Cambia: informe actualizado de seguridad, licencias y mantenimiento del grafo final, ausencia de código nativo no aprobado, auditoría de scripts de construcción, características mínimas síncronas y coste; solo si se acepta cambia `Cargo.toml`/`Cargo.lock`.
  - Se comprueba con: `cargo tree --locked -e features` y revisión que confirme que no se activa `async`, que cada `build.rs` tiene efectos acotados y que no ejecuta procesos, usa la red, modifica fuentes ni genera o compila código nativo no aprobado; si falla la puerta, se detienen E3–E8 y se eleva una nueva decisión, sin sustituirla silenciosamente.

  **Evidencia ejecutada (10 de agosto de 2026):** la política se corrigió para
  evaluar efectos en vez de prohibir cualquier `build.rs`. La candidata compila
  sin `async` ni código nativo. Los cuatro scripts transitivos de Crossbeam y
  Rayon solo emiten configuración para Cargo: no ejecutan procesos, usan red,
  escriben fuentes ni generan o compilan código. Se acepta y fija la versión
  exacta; el coste transitivo de BGZF, compresión y paralelismo queda registrado
  y vigilado en `evaluaciones-dependencias-f1.md`. E3–E8 quedan desbloqueados.

- ✅ E2. Implementar el contrato de dominio de referencia y secuencia independiente del adaptador FASTA.
  - Depende de: D2–D4; no depende de la candidata FASTA ni de otro adaptador externo.
  - Procede de: DTC-039 a DTC-041; REQ-REF-02/03/04/05/07/22/23 y REQ-MAN-86.
  - Cambia: tipo de secuencia DNA validada y reutilizable por futuros adaptadores, `ReferenceId`, secuencias ordenadas, ID de primer token, cabecera cruda, descripción reversible, clase cromosoma/`scaffold`, longitud, origen y catálogo de índices.
  - Se comprueba con: `cargo test --locked` sobre A/C/G/T/N en ambas cajas, símbolo inválido con posición, conservación exacta de `N`, una/dos secuencias, IDs repetidos, orden, descripción UTF-8/no UTF-8 y clasificación obligatoria, sin usar el lector FASTA como único oráculo.

  **Evidencia ejecutada (10 de agosto de 2026):** ocho pruebas de dominio y
  Clippy sin avisos cubren alfabeto, `N`, posición base 1, cabecera cruda,
  descripción reversible incluido su espacio de frontera, conservación de las
  bases normalizadas, orden, clase, origen, duplicados e índice FAI bajo demanda.
  El contrato no importa ni expone tipos de `noodles-fasta`.

- ✅ E3. Implementar validación por flujo del FASTA suministrado sin modificarlo.
  - Depende de: E1–E2, D4 y A4.
  - Procede de: DTC-041; REQ-REF-12/13/15/16/23 y REQ-VAL-45.
  - Cambia: adaptador que exige texto sin comprimir, al menos una secuencia no vacía, ID no vacío/único y alfabeto A/C/G/T/N en cualquier caja; conserva `N`, caja en la copia y cabecera cruda.
  - Se comprueba con: fixtures válido, vacío, secuencia vacía, ID vacío/repetido, descripción, CRLF, descripción no UTF-8, gzip, gap, otros IUPAC y símbolo inválido con fuente/posición; cada inválido se rechaza antes de publicación.

- ✅ E4. Copiar el FASTA suministrado byte a byte al provisional mientras calcula y verifica SHA-256.
  - Depende de: D4, D8–D10 y E3.
  - Procede de: DTC-040, DTC-041 y DTC-045; REQ-REF-17/18 y REQ-SAL-03.
  - Cambia: flujo de copia que deja intacto el origen, verifica tamaño/digest, cataloga la copia cerrada y hace que la ejecución use esa copia.
  - Se comprueba con: `cmp` entre origen y copia y `sha256sum` sobre ambos; truncado o cambio de un byte impide catalogar y publicar completo.

  **Corrección verificada H-F15-R02:** tras catalogar la copia, la ejecución
  vuelve a leerla para construir la representación consumida por configuración
  efectiva y manifiestos. La ruta original permanece únicamente como
  procedencia; una regresión muta el original antes de esa carga y demuestra
  que cabecera y bases proceden de la copia.

- ✅ E5. Implementar generación de referencia mínima según composición, longitudes y semilla efectivas.
  - Depende de: B9, B10, C5, E2 y catálogo A5; la composición general consume los valores ya validados y normalizados por B10.
  - Procede de: REQ-CONF-04/05, REQ-REF-01/02/03/04/22/23 y DTC-037/DTC-041/DTC-049.
  - Cambia: generador que crea únicamente A/C/G/T en secuencias nuevas, usa un flujo nombrado, respeta número/orden/clase/longitudes/composición efectiva y no introduce `N` al azar.
  - Se comprueba con: casos literales de longitud total y composición degenerada, más repetición con misma semilla y cambio controlado de semilla; la secuencia esperada pequeña se deriva con un oráculo de prueba independiente del generador.

  **Adaptación aprobada:** el catálogo de ámbitos v2 añade
  `reference_generation(reference_id)` sin alterar la derivación v1 ni los
  vectores del catálogo v1. Los vectores v2 prueban aislamiento y repetición.

- ✅ E6. Escribir FASTA generado de forma canónica y versionar toda elección que afecte a sus bytes.
  - Depende de: E1–E2 y E5.
  - Procede de: DTC-041; REQ-REF-03/04 y REQ-SAL-01.
  - Cambia: primero registra y aprueba el ancho exacto como parte de la versión canónica y después implementa el escritor por flujo con A/C/G/T/N mayúsculas, LF, ese ancho fijo y orden estable.
  - Se comprueba con: fixture literal de FASTA pequeño, validación interna, conteo de líneas/bases y repetición; la equivalencia semántica se comprueba separada de los bytes.

  **Elección aprobada:** FASTA canónico v1 usa 80 bases por línea, LF,
  mayúsculas y orden de catálogo; la prueba literal cubre la frontera 80/81.

- ✅ E7. Crear índices únicamente cuando una operación F1 los necesite.
  - Depende de: E1, E3 o E6.
  - Procede de: DTC-041; REQ-REF-05/06.
  - Cambia: contrato de índice bajo demanda y catálogo cero/uno/varios; F1 no genera índices anticipados sin consumidor.
  - Se comprueba con: caso sin operación indexada que conserva cero índices y caso que exige FAI; `samtools faidx` valida el FASTA y consulta secuencias conocidas en pruebas/CI, nunca durante una ejecución normal.

  **Frontera aprobada:** toda ejecución normal F1 conserva cero índices. La
  consulta indexada de aceptación crea el FAI bajo demanda fuera del resultado
  publicado y demuestra el consumidor con un oráculo `samtools` independiente.

- ✅ E8. Integrar origen, cabeceras, secuencias, índice, tamaño y digest en configuración efectiva y manifiesto.
  - Depende de: D5–D8 y E2–E7.
  - Procede de: REQ-REF-07/22, REQ-MAN-06/86 y REQ-SAL-01.
  - Cambia: representación estructurada que distingue simulada/suministrada, enlaza el FASTA único y sus índices y conserva ID, descripción, cabecera reversible, longitud y clase por secuencia.
  - Se comprueba con: JSON Schema independiente, `sha256sum`, enlaces resolubles y comparación contra el FASTA leído por un oráculo distinto.

- ✅ E9. Ejecutar los recorridos extremo a extremo suministrado y simulado, válidos e inválidos.
  - Depende de: B12, C8, D13 y E2–E8.
  - Procede de: resultado y evidencia de cierre de F1; REQ-REF-23 y ACE-01/02/03/04/05/23/64/110.
  - Cambia: `TEST-REQ-REF-23`, fixtures y pruebas CLI que producen configuración original, efectiva, semilla, referencia, manifiestos, estados y diagnósticos.
  - Se comprueba con: FASTA suministrado con `N` conservada y símbolo inválido localizado, referencia generada solo con A/C/G/T, dos ejecuciones equivalentes comparadas lógicamente, `samtools faidx`, `sha256sum`, JSON independiente y configuración inválida sin destino completo.

- ☐ E10. Establecer y ejecutar la línea base A de 100 kbp en cada cambio relevante.
  - Depende de: E9 y construcción `--release` portable.
  - Procede de: DTC-046; REQ-NOF-01 a REQ-NOF-08.
  - Cambia: escenario A con configuración y semilla fijas y captura de tiempo total/etapas, CPU fiable, RSS máximo, bytes de entrada/temporales/salida, MB/s, ficheros y contexto.
  - Se comprueba con: ejecución repetida de la orden versionada del banco en WSL2/ext4; corrección, integridad y trazabilidad son obligatorias, pero una variación temporal pequeña no provoca fallo.

  **Medición ejecutada, paso reabierto:** la corrección H-F15-R04 distingue el
  preflight del recorrido completo y el residuo temporal medido, pero el pico de
  bytes temporales sigue sin instrumentación directa y se declara
  `not_measured`. E10 no se aprueba hasta medir ese pico sin inferirlo.

- ☐ E11. Ejecutar la serie B y programar la observación C sin inventar umbrales.
  - Depende de: E10.
  - Procede de: DTC-046.
  - Cambia: referencias deterministas de 1, 10 y 100 Mb para B y caso de 1 Gb manual/programado para C, con mismas métricas y contexto; conserva resultados como línea base observada.
  - Se comprueba con: informe que muestre crecimiento y número de repeticiones en construcción optimizada portable; C no forma parte de cada cambio y ningún tamaño se presenta como promesa de producto.

  **Serie ejecutada, cierre dependiente:** B se repitió tres veces por tamaño y
  C quedó programado, pero E11 conserva abierto su cierre mientras E10 no mida
  el pico temporal exigido.

- ✅ E12. Contrastar la base de distribución F1 mediante `.tar.gz` e imagen OCI sin adelantar la publicación F9.
  - Depende de: A7, E9 y DTC-042.
  - Procede de: DTC-042 y DTC-047.
  - Cambia: prueba de empaquetado mínimo versionado y receta OCI fijada por digest que incluyen ejecutable, esquemas, licencias y ejemplos; no declara una entrega v1.
  - Se comprueba con: listado y extracción estándar de `.tar.gz`, SHA-256 de sus bytes, inspección de descriptores/configuración/capas y arranque de prueba OCI mediante los comandos exactos aprobados en A7, en Ubuntu 24.04 nativo/OCI.

- ☐ E13. Cerrar la batería trazable de F1.5 y el expediente de F1.
  - Depende de: E1–E12 y todas las pausas anteriores.
  - Procede de: los 27 REQ asignados a F1.5 tras las auditorías F1.3/F1.4, ACE-01 a ACE-05/23/64/110 y reglas de avance de la hoja de ruta.
  - Cambia: casos `TEST-REQ-*`, matriz final con comandos/salidas/digests y acta de revisión; no modifica requisitos ni marca pruebas no ejecutadas como superadas.
  - Se comprueba con: batería Cargo bloqueada por warnings, validadores independientes de F1, escenarios A/B/C según frecuencia y revisión de cero hallazgos bloqueantes o graves sin destino aprobado.

  **Auditoría ejecutada, cierre no aprobado (10 de agosto de 2026):** las 27
  filas se reauditaron y las puertas completas pasan, pero solo 25 tienen
  consumidor completo. REQ-REF-05 demuestra cero y un FAI, no «varios» con un
  segundo consumidor real; REQ-EVO-63 conserva resultados reales pendientes en
  F5.6. E13 permanece abierto y F1 no se declara cerrado.

### Pausa obligatoria de revisión de F1.5 y cierre de F1

Antes de pasar a F2 se revisan conjuntamente reproducibilidad, diagnósticos,
trazabilidad, FASTA, integridad y estados de publicación. Deben existir al menos
un recorrido suministrado, uno simulado, uno inválido de admisión, un fallo
controlado y una interrupción. F1 solo se cierra si las 61 filas asignadas
actualmente a F1 tienen evidencia ejecutada, los 31 casos reasignados a F2–F9
desde F1.2/F1.3/F1.4 conservan un destino verificable y no quedan defectos
bloqueantes o graves sin resolución o destino aprobado.

## 10. Riesgos de error silencioso y pruebas independientes

| Riesgo | Protección planificada | Oráculo independiente |
| --- | --- | --- |
| Aceptar versión, campo o duplicado desconocido y completar parcialmente | B1–B4, B12; frontera de versión anterior al parseo semántico | Fixtures negativos con bytes y ubicaciones literales |
| Conceder precedencia por orden TOML o desempatar máximas incomparables | B5–B8 | Tabla manual del retículo y permutaciones |
| Normalizar con tolerancia, orden o precisión distintos de DTC-049 | B10 y D11 | Vectores de borde, suma Neumaier y bits `f64` calculados fuera del productor |
| Cambiar resultados al añadir consumidores, hilos u observaciones | A5, C4–C8 | Vectores normativos y reordenamientos adversariales |
| Mezclar clases o ámbitos de ID o ignorar una colisión | A6, D2–D3 | Vectores literales y colisión inyectada |
| Calcular el digest tras normalizar o copiar solo parte del fichero | D4, E4 | `cmp` y `sha256sum` sobre bytes completos |
| Publicar completo antes de cerrar/validar, sobrescribir o mezclar carreras | D8–D12 | Dos procesos, fallos inyectados y SIGKILL en ext4/Ubuntu/OCI |
| Reescribir un FASTA suministrado o perder cabecera/no UTF-8/`N` | E2–E4 | `cmp`, fixtures crudos y `samtools faidx` |
| Usar el propio lector/escritor como único validador | A7, D5/D13, E3/E7/E9 | JSON Schema distinto, `sha256sum`, `samtools`, tar/OCI |
| Convertir WSL2 o `target-cpu=native` en contrato de publicación | A1/A3, E10–E12 | CI Ubuntu 24.04 nativo y OCI; inspección de flags/metadatos |
| Introducir alcance DTF o consumidores F2–F9 al preparar contratos comunes | revisiones por hito y matriz §11 | Diff contra hoja de ruta y registro DTF |

## 11. Matriz final de cobertura y evidencias

La matriz no sustituye las filas `TEST-REQ-*` del artefacto de trazabilidad. Las
auditorías de F1.2, F1.3 y F1.4 detectaron respectivamente 20, 6 y 5
asignaciones que exigían consumidores de F2–F9; además, F1.3 trasladó 3
requisitos y F1.4 trasladó otros 3 al cierre ejecutable de F1.5. Las matrices
de auditoría los conservan sin marcarlos como superados. La asignación vigente
de F1 es por tanto 16/1/17/27, 61 requisitos en total.

| Hito | Requisitos asignados por la matriz | ACE principales | Decisiones | Pasos/evidencia de cierre |
| --- | --- | --- | --- | --- |
| F1.1 | Infraestructura de la hoja de ruta; sin filas REQ exclusivas | soporte de ACE-01/03/22/23/64 | DTC-042 a DTC-048 | A1–A7; toolchain exacto, lock, cinco puertas Cargo, catálogos previos y herramientas independientes aprobadas |
| F1.2 | REQ-CONF-01/04/05/07/09/10/14/29/30/31/65/106/117–119; REQ-VAL-51 | ACE-01/02/10/11/105/111 | DTC-034, DTC-036, DTC-043, DTC-044, DTC-048, DTC-049 | B1–B12; 16 `TEST-REQ-*` superados y auditoría explícita de 20 reasignados; dominancia, invalidación previa, normalización Neumaier y diagnósticos |
| F1.3 | REQ-REP-23 | ACE-01/23/64 | DTC-037, DTC-038, DTC-042/043/048 | A5, C1–C8; 1 `TEST-REQ-*` superado y auditoría de 9 reasignados con pruebas parciales; vectores, aislamiento y comparadores preparatorios |
| F1.4 | REQ-CONF-38/101/123; REQ-MAN-01–05/08/32/69/81/88; REQ-SAL-03/04/22/26 | ACE-01/10/14/22/25/34/64/93/94/105/111 | DTC-035/036/038–040/042–045/047–049 | A4/A6, D1–D13; 17 `TEST-REQ-*` superados y auditoría explícita de 8 reasignados; JSON/SHA, evidencia DTC-049, estados y carrera no-clobber |
| F1.5 | REQ-CONF-02/03/26/27; REQ-REF-01–07/12/13/15–18/22/23; REQ-REP-22; REQ-EVO-63; REQ-VAL-45; REQ-MAN-06/72/86; REQ-SAL-01/18 | ACE-01–05/23/64/110 | DTC-037, DTC-039–049 | E1–E13; 27 `TEST-REQ-*`, incluida repetición ejecutable con semilla/configuración/entrada/manifiesto; puerta `noodles-fasta`, `cmp`, `samtools faidx`, SHA/JSON independientes y escenarios A/B/C |

### 11.1. Expediente mínimo para aprobar F1

- Salidas completas de formato, Clippy sin avisos, pruebas y construcción
  `--locked` portable.
- Lock y evaluaciones de todas las dependencias y herramientas usadas.
- Catálogos y vectores normativos de diagnósticos, subsemillas e IDs.
- Resultado de los 61 `TEST-REQ-*` asignados actualmente a F1, enlazado a su
  paso y commit, más el registro auditable de los 31 reasignados a F2–F9 desde
  F1.2/F1.3/F1.4.
- Validación independiente de JSON, SHA-256, FASTA, `.tar.gz` y OCI.
- Evidencia de publicación completa, rechazo previo, fallo controlado,
  interrupción y carrera no-clobber.
- Comparación de dos ejecuciones equivalentes en los cuatro niveles, sin afirmar
  identidad binaria no demostrada.
- Línea base A, serie B y registro programado/manual de C con contexto completo.
- Actas de las cinco pausas de revisión y confirmación de que no se incorporó
  alcance DTF ni trabajo biológico de F2–F9.

## 12. Validación documental del plan

Antes de aprobar este documento se ejecutan las comprobaciones disponibles de
§4.2 y se revisa manualmente que:

- cada paso A1–E13 tiene dependencia, origen, cambio y comprobación concreta;
- todas las casillas permanecen abiertas hasta ejecutar el trabajo;
- existen cinco pausas obligatorias, una por hito;
- la distribución vigente 16 + 1 + 17 + 27 = 61 y la reasignación auditada de
  31 casos a F2–F9 desde F1.2/F1.3/F1.4 están reproducidas en pasos,
  trazabilidad y matriz;
- DTC-034 a DTC-049 están cubiertas sin convertir candidatas o DTPD en hechos;
- `noodles-fasta` 0.66.0 solo se incorpora tras superar E1;
- los catálogos A4–A6 preceden a sus consumidores;
- los comandos Cargo se declaran futuros y no ejecutados.
