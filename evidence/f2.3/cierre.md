# Cierre F2.3 — individuos directos

## Alcance comprobado

F2.3 queda cerrada como integración de cohortes iniciales sobre una matriz SNV
común: configuración e identidades estables (IDs aportados o derivados de
`founder_keys` explícitas, únicas dentro de cada población), muestreo HWE o `F` diploide,
fase independiente, eventos explícitos dirigidos, materialización temporal bajo
demanda y estadísticas de fundadores. Tanto `validate` como `run` realizan el
preflight y descartan las copias temporales; solo `run` publica estadísticas,
nunca copias ni secuencias. La UI Marimo exporta ese mismo TOML y transmite solo
las autorizaciones específicas ya implementadas.

La configuración efectiva y el manifiesto enlazan la matriz canónica con su
metadata, auditoría de normalización y estadísticas iniciales, sin copias ni
secuencias materializadas. Estas conservan
procedencia, modelo, algoritmos, objetivos y realizaciones; una deleción que
cubre un locus se declara como `missing_by_deletion` y no se autocorrige. La
revalidación vuelve a comprobar TOML, configuración efectiva, manifiesto,
esquemas, enlaces físicos y estadísticas antes de aceptar una publicación.

La evidencia está distribuida por frontera para evitar un segundo ejecutor:

- `tests/initial_populations_config.rs`, `tests/initial_genotypes.rs` y
  `tests/initial_individual_materialization.rs` cubren contratos, RNG,
  oráculos HWE/`F`, fase, dosis, aislamiento y estadísticas.
- `tests/initial_individual_integration.rs` y `tests/cli_smoke.rs` cubren los
  eventos directos, el preflight compartido y las autorizaciones limitadas.
- `tests/{effective_config_schema,manifest_minimal,published_resource_validation}.rs`
  cubren cierre de esquemas, coherencia entre documentos y revalidación física.
- `tests_ui/test_f2_3_individuals_ui.py` y los tests UI transversales cubren
  la paridad UI→TOML→Rust y los diagnósticos localizados.

## Validación ejecutada

Se superaron `cargo fmt --all -- --check`,
`cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`,
`cargo test --workspace --all-targets --all-features --locked`,
`cargo test --workspace --doc --all-features --locked` y
`cargo doc --workspace --no-deps --all-features --locked`, además de
`cargo build --locked`.

La suite UI focal también superó
`uv run --frozen python -m unittest tests/test_population_snvs_editor.py tests_ui/test_ui_contract.py tests_ui/test_f2_3_individuals_ui.py tests_ui/test_process_adapter.py`
(34 pruebas, una omitida por su condición declarada), al igual que
`MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py`
y `git diff --check`.

El contrato de claves se prueba además por CLI y revalidación publicada:
`tests/cli_smoke.rs` comprueba que `run` publica `founder_keys`, y
`tests/published_resource_validation.rs` altera una clave, actualiza el hash
del artefacto y confirma que la revalidación contra el TOML la rechaza.

## Límites que se conservan

- F2.3 no incorpora generación individual aleatoria por cantidad o densidad,
  mutaciones *de novo*, VCF, paneles haplotípicos, escenarios, bundles ni OCI.
- REQ-CONF-20/22, REQ-IND-15 y REQ-DIR-02 permanecen planificados en F6.5.a:
  cantidad/densidad, distribuciones/regiones y sustituciones por individuo de
  la creación inicial aleatoria directa. REQ-IND-64 permanece parcial: F6.5.a
  cierra esa creación y su exclusión frente a matriz, mientras F7.2.a cierra
  las mutaciones *de novo*. REQ-CONF-128 conserva evidencia parcial F2; F5.6.a
  es el destino explícito de matrices canónicas separadas por referencia.
- Marimo no ofrece controles de fases posteriores ni sustituye la validación
  semántica Rust; abrir la aplicación no ejecuta una simulación.
