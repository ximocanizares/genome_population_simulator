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

**Reevaluación F1.4.5 (20 de agosto de 2026).** Se retiró `jsonschema` del
grafo Cargo después de eliminar su último consumidor. Los contratos Draft
2020-12 continúan publicados como ficheros estáticos y las pruebas los ejecutan
con CPython/jsonschema 4.10.3, ajeno al productor Rust. `Cargo.lock` quedó en 65
paquetes y las búsquedas en `Cargo.toml`, `Cargo.lock` y `cargo tree -e
features` no contienen el crate retirado. `cargo metadata --locked`, `cargo
tree --locked -e features` y `cargo audit --no-fetch` terminaron con código 0;
la auditoría examinó los 65 paquetes y solo advirtió que no pudo tomar el lock
del índice local, sin vulnerabilidades informadas.

**Corrección durante D11.** La primera prueba de ida y vuelta del manifiesto
detectó que compilar `serde_json` solo con `std` recuperaba un vector a dos ULP
del valor original. Se añadió `float_roundtrip`, se conservó el mismo contrato
público y la repetición exige igualdad de `to_bits()` para tolerancia, sumas,
deltas y todos los valores antes/después. La selección sin esta característica
queda rechazada por evidencia ejecutada.

## E1 — reevaluación de `noodles-fasta` 0.66.0

**Decisión (10 de agosto de 2026): aceptada por la puerta de F1.5 tras corregir
la política de construcción.** La versión exacta se incorpora sin características
predeterminadas en la dependencia directa y sin `async`. La compilación con
Rust 1.97.1 termina correctamente. La nueva regla no prohíbe todo `build.rs`,
sino que exige auditar sus efectos y rechaza los que introduzcan acciones no
necesarias, no acotadas o no aprobadas.

- **Propósito y contrato:** la candidata cubriría el adaptador síncrono de
  lectura/escritura FASTA y FAI de DTC-041 tras tipos propios. No tendría
  autoridad sobre alfabeto, cabeceras, identidad, diagnósticos ni conservación
  byte a byte.
- **Versión, mantenimiento y API:** `noodles-fasta 0.66.0`, publicada el 3 de
  agosto de 2026, declara edición 2024, Rust 1.89.0, licencia MIT y
  `build = false`. La serie publicó 0.63–0.66 entre julio y agosto de 2026 y su
  historial confirma cambios frecuentes de API, mitigables mediante adaptador
  propio, pero no relevantes para el bloqueo observado.
- **Características:** la raíz se declaró como
  `default-features = false`; `cargo tree --locked -e features` confirmó que
  `tokio` y la característica `async` no se activaron. La dependencia directa
  de `noodles-fasta` sobre `noodles-bgzf 0.51.0`, sin desactivar sus
  características predeterminadas, incorporó sin embargo `bytes`,
  `crossbeam-channel`, `rayon` y `zlib-rs`.
- **Licencias y seguridad:** las licencias observadas fueron MIT,
  MIT/Apache-2.0, Unlicense/MIT y Zlib, compatibles con la distribución MIT.
  `cargo audit --no-fetch` cargó 1199 avisos y examinó el lock temporal de 65
  paquetes con código 0; avisó de que no podía bloquear el índice de crates.io
  de solo lectura.
- **Código nativo y scripts:** no se encontraron fuentes C/C++, enlace nativo
  ni ejecución de procesos en `noodles-fasta`, `noodles-bgzf` o
  `noodles-core`; `zlib-rs` es Rust. El árbol sí contiene y ejecuta cuatro
  scripts: `crossbeam-utils 0.8.22/build.rs`,
  `crossbeam-epoch 0.9.20/build.rs`, `crossbeam-deque 0.8.7/build.rs` y
  `rayon-core 1.13.0/build.rs`. Los tres scripts de Crossbeam consultan el
  objetivo y la configuración de sanitizadores y emiten directivas `cfg`/
  `check-cfg`; el de Rayon solo emite `rerun-if-changed`. Ninguno ejecuta
  procesos, accede a la red, modifica fuentes, escribe artefactos generados ni
  genera o compila código nativo. Superan la regla revisada de DTC-043.
- **Rendimiento y distribución:** incluso para FASTA sin comprimir y síncrono,
  el grafo añade BGZF, compresión y paralelismo no requeridos por el recorrido
  mínimo. La API publicada por `noodles-fasta` no permite al consumidor
  desactivar las características predeterminadas de su dependencia directa
  `noodles-bgzf`, por lo que no existe una selección de características de la
  aplicación que elimine esos scripts.

El coste de BGZF, compresión, Rayon y `zlib-rs` se acepta de forma explícita
porque la API publicada no permite retirarlo desde la aplicación. `zlib-rs` es
Rust puro, contiene `unsafe` interno y selecciona optimizaciones mediante
detección de CPU en tiempo de ejecución; no enlaza una biblioteca nativa. Este
coste y superficie se vigilan en futuras actualizaciones. La aceptación no
autoriza parchear o sustituir silenciosamente la dependencia ni relajar los
contratos propios de DTC-041/DTC-043. E3–E8 quedan desbloqueados.

## F2.5.1 — evaluación DTC-043 de VCF/BGZF/índices

**Decisión (24 de agosto de 2026): aceptada para el contrato físico F2.5.1.**
Se incorporan `noodles-vcf = 0.88.0`, `noodles-csi = 0.56.0`,
`noodles-tabix = 0.62.0` y `noodles-bgzf = 0.47.0`, todos con pin exacto y
`default-features = false` en la dependencia directa. `Cargo.lock` fija la
resolución. Su único propósito es escribir VCF 4.5, BGZF y exactamente un
índice TBI o CSI desde el adaptador `vcf_contract`; no proyectan alelos, GT,
verdad ni validadores externos. La cabecera y la selección TBI/CSI se
comprueban con pruebas propias; `bcftools` y HTSlib siguen fuera del runtime.

**Ampliación F2.5.5 (25 de agosto de 2026).** `noodles-core = 0.20.0`, ya
resuelto transitivamente por la misma familia Noodles, se declara también de
forma directa con pin exacto y `default-features = false`. Solo aporta el tipo
de posición para la relectura y consulta regional local del índice ya escrito;
no cambia el grafo de `Cargo.lock`, no abre una frontera pública ni añade una
biblioteca nativa o un validador externo.

**Ampliación F3.1.4 (28 de agosto de 2026).** La salida FASTQ comprimida
reutiliza `noodles-bgzf = 0.47.0`, ya dependencia directa aprobada, mediante un
adaptador privado de escritura y lectura por flujo. No cambia `Cargo.toml`,
`Cargo.lock`, sus *features*, el grafo, licencias, scripts ni superficie nativa
ya evaluados. BGZF es un perfil de gzip interoperable: se fija DEFLATE nivel 6,
el tamaño de bloque de la versión bloqueada y `MTIME = 0`; los SHA-256 y tamaños
se calculan sobre los bytes `.fastq.gz` finales, mientras el contrato FASTQ se
valida tras descompresión en una frontera distinta. La prueba usa `gzip -cd`
como decodificador independiente de interoperabilidad. Este uso no incorpora
índices BGZF, paralelismo ni una opción de compresión configurable; la
repetibilidad binaria queda condicionada por el simulador, lock y plataforma
admitida, como exige DTC-038.

- **Estándar, API y mantenimiento:** `noodles-vcf` implementa el escritor VCF
  y declara `FileFormat::new(4, 5)`; `noodles-tabix` y `noodles-csi` escriben
  los índices TBI y CSI, y `noodles-bgzf` finaliza el flujo BGZF. La fuente
  declarada de los cuatro es el repositorio `zaeleus/noodles`; el estado
  evaluado es exactamente el de los paquetes de crates.io fijados en el lock:
  VCF 0.88.0, CSI 0.56.0, Tabix 0.62.0 y BGZF 0.47.0. Sus cuatro
  `.cargo_vcs_info.json` identifican el mismo snapshot
  `f6a9dd58400a30084226c9c1d0fc300a9468c15a`, con los directorios respectivos
  `noodles-vcf`, `noodles-csi`, `noodles-tabix` y `noodles-bgzf`. No se pudo
  verificar actividad o fecha posterior en línea desde este entorno; la
  aceptación queda condicionada a esos paquetes exactos, al lock y a repetir
  esta evaluación antes de cualquier actualización. Todos declaran licencia
  MIT, `build = false` y MSRV 1.89.0, por debajo del toolchain Rust 1.97.1.
  Sus series `0.x` no prometen estabilidad semver de API: los tipos Noodles no
  salen del contrato público y las regresiones cubren cabecera, índice y
  cierre. No se activa `async`; `cargo tree --locked --offline -e features`
  confirma el grafo efectivo.
- **Grafo, licencias y seguridad:** las seis entradas nuevas de `Cargo.lock`
  son los cuatro crates Noodles anteriores, `bit-vec 0.9.1` (Apache-2.0 OR
  MIT) y `percent-encoding 2.3.2` (MIT OR Apache-2.0). `indexmap 2.14.0`
  (Apache-2.0 OR MIT) ya estaba resuelto, pero pertenece al subgrafo nuevo de
  CSI/VCF/Tabix, no a `noodles-fasta`. `bstr 1.13.0` (MIT OR Apache-2.0),
  `bytes 1.12.1` (MIT), `crossbeam-*` (MIT OR Apache-2.0), `rayon*`
  (Apache-2.0 OR MIT) y `zlib-rs 0.6.7` (Zlib) sí estaban resueltos por
  `noodles-fasta 0.66.0`; este cambio los reutiliza, no los añade. Todas son compatibles con la
  licencia MIT del proyecto. `cargo audit --no-fetch` cargó 1.216 avisos,
  examinó las 72 dependencias del lock y terminó con código 0. `cargo-license`
  no está instalado en este entorno; las licencias anteriores se comprobaron
  en los `Cargo.toml` exactos de la caché y deben repetirse con esa herramienta
  cuando esté disponible, sin reinterpretar su ausencia como una aprobación.
- **Construcción y código nativo:** los cuatro crates Noodles nuevos y
  `bit-vec` declaran `build = false`; no introducen fuente C/C++, enlace
  nativo ni procesos externos. Sus transitivas reutilizan los cuatro
  `build.rs` ya auditados de Crossbeam/Rayon en E1: solo emiten configuración
  Cargo y no ejecutan procesos, red, generación de fuentes ni compilación
  nativa. `zlib-rs` permanece Rust puro con `unsafe` interno y detección de
  CPU; no enlaza zlib del sistema.
- **Rendimiento, distribución y duplicación:** la escritura es por flujo y no
  materializa VCF completo. El coste nuevo es real: por incompatibilidad entre
  `noodles-fasta 0.66.0` y la familia VCF, el lock conserva las dos versiones
  `noodles-bgzf 0.51.0` y `0.47.0`; ambas comparten las transitivas de
  compresión, pero duplican código compilado y superficie de actualización.
  No existe una actualización compatible dentro de F2.5.1 que lo elimine, por
  lo que se acepta acotadamente y se reevaluará al actualizar Noodles. No se
  añadieron dependencias nativas ni validadores como dependencia de ejecución.
- **Evidencia reproducible y límite pendiente:** se ejecutaron `cargo metadata
  --locked --offline --format-version 1 --no-deps`, `cargo tree --locked
  --offline -e features`, `cargo audit --no-fetch`, las pruebas focales y las
  puertas Cargo. La puerta externa DTC-066 es
  `scripts/validate_vcf_external.sh`: exige exactamente `bcftools 1.24` y
  `tabix (htslib) 1.24`, exige el sidecar exclusivo correcto, conserva en un
  temporal la salida de `view` y `norm --check-ref e`, las compara sin reparar,
  consulta por región y reconstruye el mismo tipo de índice solo sobre copias
  temporales. Sus mutaciones REF, índice opuesto y doble índice deben fallar.
  El entorno actual aporta 1.21 y 1.19, por tanto esa puerta queda **pendiente
  y no superada**; la sonda de prueba solo informa versiones y nunca convierte
  ese estado en éxito.
