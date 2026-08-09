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
