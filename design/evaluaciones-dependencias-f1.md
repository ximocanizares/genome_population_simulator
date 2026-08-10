# Evaluaciones de dependencias de F1

## Estado

Registro iniciado en F1.1 el 8 de agosto de 2026. Una entrada aprobada autoriza
solo el propósito y las características indicadas; no autoriza capacidades de
fases posteriores.

## Aplicación Rust

F1.1 no incorpora crates externos. `Cargo.toml` contiene una tabla
`[dependencies]` vacía y el árbol bloqueado consta únicamente del paquete de la
aplicación. Por ello no existen licencias transitivas, código nativo, scripts de
construcción ni características de crates que evaluar en A1–A3.

Las herramientas de desarrollo (`rustfmt` y `clippy`) y el target GNU forman
parte del toolchain Rust 1.97.1 fijado en `rust-toolchain.toml`; no son
dependencias de ejecución ni se exponen en contratos públicos.

## Automatización

`actions/checkout` queda fijada al commit
`11bd71901bbe5b1630ceea73d27597364c9af683` (versión 4.2.2) exclusivamente para
obtener el repositorio en GitHub Actions. Usa licencia MIT, no entra en el
binario ni en sus resultados y no añade código nativo o scripts de construcción
al grafo Cargo. Es una acción JavaScript distribuida como artefacto cerrado, sin
características configurables relevantes para este uso. Su API de entradas es
estable en la serie 4 y el proyecto tiene mantenimiento activo. La fijación por
commit evita que una etiqueta mutable cambie la automatización sin revisión;
los permisos del flujo son de solo lectura y `persist-credentials: false` evita
conservar credenciales en la configuración Git. Su coste se limita a obtener el
repositorio una vez y no afecta al rendimiento ni a la distribución de la
aplicación. Cualquier actualización exige una revisión aislada de licencia,
procedencia, avisos de seguridad, entradas y código fijado.

Toda dependencia posterior deberá registrar: propósito y REQ atendidos, versión
exacta, licencia, mantenimiento, seguridad, estabilidad de API, rendimiento,
características activadas, árbol final, código nativo, scripts de construcción
y coste de distribución antes de modificar `Cargo.toml` o `Cargo.lock`.

## B2 — admisión TOML

**Decisión (9 de agosto de 2026): aprobada para F1.2.** Se acepta
`toml 1.1.4+spec-1.1.0` con `default-features = false` y exclusivamente las
características `parse` y `std`. La especificación Cargo fija `=1.1.4`; la
resolución completa queda fijada por `Cargo.lock` y debe comprobarse después de
incorporarla.

- **Propósito y REQ:** analizar los documentos TOML de DTC-034 y proporcionar
  rangos de bytes para B3–B4 y REQ-CONF-01/09/10/13/14. El parser detecta
  sintaxis y claves duplicadas; el adaptador propio consume su representación
  con rangos, comprueba claves desconocidas y traduce todos los errores a
  diagnósticos `GPS-*`.
- **Corrección y evidencia requerida:** implementa TOML 1.1 y conserva el
  texto numérico y los rangos de claves y valores. B2 exige fixtures que
  demuestren duplicados, desconocidos y localización; la dependencia no recibe
  autoridad sobre `DeclaredConfig`, reglas, IDs ni semántica.
- **Versión, compatibilidad y estabilidad:** `toml 1.1.4+spec-1.1.0`, publicada
  y mantenida en `toml-rs/toml`, declara MSRV 1.85, compatible con Rust 1.97.1.
  Su API de deserialización de bajo nivel es pública pero puede evolucionar; el
  riesgo queda contenido en `adapters` y en pruebas de bytes/ubicaciones.
- **Licencias:** `toml`, `serde_spanned`, `toml_datetime` y `toml_parser` usan
  MIT o Apache-2.0; `winnow` usa MIT. Son compatibles con la licencia MIT del
  proyecto y deben incluirse en el inventario de distribución posterior.
- **Seguridad y mantenimiento:** la consulta del 9 de agosto de 2026 a RustSec
  y a los avisos del repositorio no encontró avisos publicados para las
  versiones del grafo. `cargo audit` 0.22.2 (sin opción `--locked`) cargó 1198
  avisos RustSec, escaneó las 16 dependencias de `Cargo.lock` y terminó con
  código 0, sin vulnerabilidades. Las fuentes muestran mantenimiento activo y
  publicaciones coordinadas recientes de la familia `toml-rs`.
- **Rendimiento:** el parser trabaja sobre el `str` de entrada y expone valores
  prestados; F1.2 realiza una sola conversión a tipos propios. No se activa
  `preserve_order`, porque el orden físico no concede precedencia, ni `unbounded`;
  se conserva una copia exacta de bytes de forma separada.
- **Grafo y distribución previstos:** con esas características, el grafo
  esperado es `toml 1.1.4`, `serde_spanned 1.1.1`,
  `toml_datetime 1.1.1+spec-1.1.0`, `toml_parser 1.1.3+spec-1.1.0` y
  `winnow 1.0.4`. Esos cinco paquetes forman el grafo de ejecución activo;
  Cargo.lock conserva además dependencias opcionales no activadas de los
  paquetes publicados. Todos son Rust puro, declaran `build = false`, no enlazan
  bibliotecas nativas y no ejecutan scripts de construcción. No se activan
  `serde`, `display`, `debug`, `simd`, `unsafe`, `preserve_order` ni
  dependencias opcionales asociadas.
- **Coste:** cinco paquetes activos en el árbol de ejecución, sin procesos
  externos ni requisitos de runtime. La funcionalidad queda tras contratos
  propios y puede sustituirse sin cambiar el dominio.

Se rechazan para este hito `toml_edit` (edición y preservación de formato no
necesarias), `toml-span` (otra API y otro grafo para una capacidad ya cubierta)
y un parser propio (riesgo innecesario de incumplir TOML y sus duplicados). La
copia exacta se conserva como bytes originales y no requiere un árbol editable.

## C1 — CSPRNG, BLAKE3 y ChaCha12

**Decisión (9 de agosto de 2026): aprobada para F1.3.** Se aceptan
`getrandom 0.4.3`, `blake3 1.8.6` y `rand_chacha 0.10.0`, fijados en Cargo como
versiones exactas. `blake3` se compila sin características predeterminadas y
solo con `pure`; `getrandom` y `rand_chacha`, sin características
predeterminadas. La resolución final y sus características se comprueban tras
actualizar `Cargo.lock`.

- **Propósito y REQ:** `getrandom` obtiene exactamente los 32 bytes de la raíz
  ausente desde el CSPRNG del sistema operativo (REQ-CONF-26); `blake3`
  materializa exclusivamente la derivación con dominio de DTC-037; y
  `rand_chacha` proporciona `ChaCha12Rng` local por ámbito para
  REQ-CONF-111–115 y REQ-REP-22/23. Los tres quedan detrás de contratos propios
  de raíz, entropía, derivación, flujo y comparación; ningún tipo de error de
  crate forma parte del diagnóstico público.
- **Corrección:** las APIs seleccionadas ofrecen llenado completo o error,
  `derive_key` con contexto explícito y la variante concreta de 12 rondas. Los
  materiales y digests se contrastan con los cuatro vectores literales A5,
  calculados por un oráculo independiente. Se descartan `StdRng`, un RNG global,
  `seed_from_u64` y el consumo compartido.
- **Versiones, mantenimiento y estabilidad:** las tres versiones proceden de
  los proyectos oficiales BLAKE3 y rust-random, con mantenimiento activo.
  `getrandom 0.4.3` y `rand_chacha 0.10.0` declaran Rust 1.85; BLAKE3 no declara
  MSRV en sus metadatos publicados, por lo que la compatibilidad se demuestra
  construyendo con el toolchain fijado 1.97.1. Las APIs se encapsulan para que
  una actualización no cambie contratos propios; cualquier cambio repite
  vectores y comparadores.
- **Licencias:** `getrandom` y `rand_chacha` usan MIT o Apache-2.0. `blake3` usa
  CC0-1.0 o Apache-2.0, incluida la alternativa con excepción LLVM. Las
  transitivas se inventariarían y revisarían desde el árbol bloqueado antes del
  cierre; son compatibles con la distribución MIT del proyecto.
- **Seguridad:** no se usa BLAKE3 como integridad publicada ni ChaCha12 como
  fuente de la raíz. La raíz procede del CSPRNG del sistema y los flujos locales
  son deterministas desde subsemillas de 256 bits. `cargo audit` se ejecuta
  sobre el lock final; cualquier aviso bloqueante invalida esta aceptación.
- **Características, código nativo y scripts:** `blake3/pure` selecciona las
  implementaciones Rust y evita compilar ensamblador o C; su paquete sí contiene
  un `build.rs` y la dependencia de construcción `cc`, cuyo coste se acepta
  únicamente para la selección portable. `rand_chacha` y `rand_core` declaran
  `build = false`; `ppv-lite86` usa SIMD portable de Rust. `getrandom` contiene
  un `build.rs` de selección de backend y, en Linux GNU, usa `libc` para la
  frontera del sistema operativo; no enlaza una biblioteca C externa ni ejecuta
  una herramienta en runtime. No se activan `serde`, `rayon`, `mmap`, `digest`,
  `wasm_js`, `sys_rng` ni características de RNG del sistema en `rand_chacha`.
- **Rendimiento y distribución:** la entropía se consulta una sola vez por raíz
  ausente. Las derivaciones y estados locales son pequeños; la opción BLAKE3
  pura prioriza portabilidad y ausencia de código nativo sobre SIMD ensamblado,
  apropiado para los materiales cortos de F1.3. El coste exacto del grafo queda
  registrado por `cargo tree --locked -e features`.

Se rechazan el crate agregado `rand` (API y características innecesarias),
`StdRng` (algoritmo no contractual), sembrar ChaCha directamente desde el SO
(perdería la raíz auditable) y una implementación propia de ChaCha12 o BLAKE3
(riesgo criptográfico y de incompatibilidad con los vectores normativos).

**Resultado del cierre C1.** `cargo tree --locked -e features` confirma solo
`blake3/pure`, `getrandom` sin rasgos propios y `rand_chacha` sin rasgos
predeterminados, además de las transitivas descritas. El paquete BLAKE3 conserva
su `build.rs` y `cc` como dependencia de construcción, pero la compilación
observada con `pure` no compiló C ni ensamblador. `cargo audit` 0.22.2 cargó
1198 avisos RustSec, examinó las 34 dependencias del lock ampliado y terminó con
código 0, sin vulnerabilidades. Los cuatro vectores A5 pasaron tanto por el
oráculo BLAKE3 de prueba independiente como por la implementación productiva.

## D1 — JSON, JSON Schema, SHA-256 y publicación no-clobber

**Decisión (10 de agosto de 2026): aprobada para F1.4.** Se aceptan
`serde_json 1.0.149`, `sha2 0.10.9` y `rustix 1.1.4`, fijados como versiones
exactas. `serde_json` se compila sin características predeterminadas y solo con
`std` y `float_roundtrip`; `sha2`, sin características predeterminadas; y `rustix`, sin
características predeterminadas y solo con `fs` y `std` (esta última aporta la
frontera segura de rutas `Path`). JSON Schema Draft 2020-12 se
publica como esquemas estáticos propios y se contrasta con el validador externo
ya aprobado en A7; no necesita una dependencia Cargo de validación.

- **Propósito y REQ:** `serde_json` produce los documentos RFC 8259 de
  DTC-035/DTC-036; `sha2` calcula SHA-256 incremental sobre bytes exactos para
  DTC-040; `rustix` encapsula la operación Linux
  `renameat2(..., RENAME_NOREPLACE)` de DTC-045. Ningún tipo de estos crates
  cruza los contratos propios de configuración, identidad, artefacto o
  publicación.
- **Corrección y evidencia:** `serde_json` rechaza números no finitos y escribe
  `f64` con recuperación exacta; `float_roundtrip` es necesario para que su
  lector recupere también los bits exactos de los vectores DTC-049. Sin
  `preserve_order`, los mapas mantienen orden
  determinista por clave. `sha2` implementa SHA-256 incremental y se contrasta
  con GNU `sha256sum` 9.4. `rustix::fs::renameat_with` expone
  `RenameFlags::NOREPLACE` sin patrón *check-then-rename*; se prueban destino
  existente, carrera de procesos y operación no soportada mediante doble.
- **Versiones y estabilidad:** las APIs usadas son pequeñas y quedan aisladas
  en adaptadores propios. Las versiones declaran respectivamente Rust 1.68,
  1.85 y 1.63 como mínimo, compatibles con el toolchain fijado 1.97.1. Toda
  actualización será aislada y repetirá esquemas, SHA, carrera y puertas Cargo.
- **Licencias y mantenimiento:** `serde_json` y sus dependencias usan MIT o
  Apache-2.0; `sha2` y el ecosistema RustCrypto, MIT o Apache-2.0; `rustix`,
  Apache-2.0 con LLVM-exception, Apache-2.0 o MIT. Los tres proyectos tienen
  mantenimiento activo y licencias compatibles con la distribución MIT.
- **Seguridad:** el cierre ejecuta `cargo audit` sobre el lock final. El JSON se
  valida además mediante CPython/jsonschema ajeno al productor; SHA-256 no se
  sustituye por BLAKE3; y la publicación no ofrece un *fallback* con carrera.
- **Características, código nativo y scripts:** las características elegidas no
  activan `preserve_order`, `arbitrary_precision`, `raw_value`, async ni E/S de
  red. Los crates son Rust puro. `rustix` usa la ABI Linux mediante
  `linux-raw-sys`/`libc` según plataforma, sin biblioteca C externa ni proceso
  de ejecución. Los paquetes seleccionados no requieren compilar código C/C++.
- **Rendimiento y coste:** los documentos F1 son acotados y se materializan una
  vez; los artefactos se hashean por flujo sin una segunda copia en memoria; la
  publicación realiza una única llamada atómica. El árbol y las características
  finales se conservan con `cargo tree --locked -e features`.

Se rechazan un serializador JSON propio, una implementación propia de SHA-256,
usar BLAKE3 como integridad publicada, `std::fs::rename` con comprobación previa
y llamar a `mv` como proceso externo. También se rechaza integrar un validador
JSON Schema en el runtime: DTC-047 lo reserva a pruebas y aceptación.

**Corrección durante D11.** La primera prueba de ida y vuelta del manifiesto
detectó que compilar `serde_json` solo con `std` recuperaba un vector a dos ULP
del valor original. Se añadió `float_roundtrip`, se conservó el mismo contrato
público y la repetición exige igualdad de `to_bits()` para tolerancia, sumas,
deltas y todos los valores antes/después. La selección sin esta característica
queda rechazada por evidencia ejecutada.
