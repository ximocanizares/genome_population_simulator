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
