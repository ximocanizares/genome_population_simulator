# Validación de F1.2

> **EVIDENCIA HISTÓRICA SUPERADA.** Este documento describe la implementación
> anterior conservada en el checkpoint `66ad8b2`; sus cifras de 63 pruebas y 16
> requisitos superados no son el estado vigente y no aprueban PAUSA-F1.2. La
> reimplementación mínima usa la matriz actual de cuatro `pass` y 32
> reasignaciones. Su aplicación de correcciones y puertas actuales se registra
> en `evidence/f1-reimplementation/f1.2-corrections.md` hasta completar una nueva
> revisión independiente.

Entorno observado: Rust/Cargo 1.97.1, objetivo
`x86_64-unknown-linux-gnu`, Linux WSL2 x86_64; CPython 3.12.3 con `jsonschema`
4.10.3.

Puertas repetidas después de corregir H-F12-01 a H-F12-03:

```text
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo doc --workspace --no-deps --locked
cargo tree --locked -e features
cargo audit
git diff --check
```

Resultado: todas terminaron con código 0. Pasaron 63 pruebas de integración,
sin fallos ni pruebas ignoradas; los objetivos unitarios y documentales no
contienen todavía casos propios. La documentación se construyó sin avisos.

El árbol activo no cambió: `toml` 1.1.4 con `parse,std`, `serde_spanned` 1.1.1,
`toml_datetime` 1.1.1, `toml_parser` 1.1.3 y `winnow` 1.0.4.

`cargo audit` 0.22.2 no admite `--locked`: cargó 1198 avisos RustSec, escaneó
las 16 dependencias de `Cargo.lock` y terminó con código 0, sin
vulnerabilidades.

La integración CLI ejecutó los 17 fixtures persistentes `invalid-*.toml` y los
12 casos de `invalid-rule-values.tsv` en modo humano y JSON. Los 29 casos
devolvieron 2, conservaron el mismo código en ambos renderizados y no crearon el
destino solicitado. Cada JSON fue además validado fuera del productor mediante:

```text
python3 -m jsonschema -i <diagnóstico.json> schemas/diagnostics-v1.schema.json
```

Los 29 documentos fueron aceptados por el esquema Draft 2020-12. Una instancia
estable y su metadato reproducible se conservan en este directorio.

La matriz `test-req-matrix.tsv` audita exactamente los 36 `TEST-REQ-*`
inicialmente asignados a F1.2: 16 están superados con evidencia semántica viva y
20 figuran como `reassigned`, sin evidencia ficticia, en el hito que permite
demostrar el REQ completo. La prueba automática comprueba que cada destino
coincide con `design/trazabilidad-requisitos-0.1.md` y existe realmente en
`design/hoja-de-ruta-0.1.md`.
