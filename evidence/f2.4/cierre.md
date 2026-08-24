# Cierre F2.4 — verdad estructurada

## Alcance comprobado

F2.4 publica `outputs/truth.jsonl` v1 como proyección del estado canónico F2.3
y de la matriz canónica; no consulta ni necesita VCF. El documento conserva
cabecera versionada, identidad de referencia, copias faseadas, procedencia,
segmentos, fragmentos vigentes, historial y casos genotípicos con dosis y
portadores. El escritor usa la publicación F1 y el manifiesto cataloga el
artefacto con la ruta canónica, tamaño y SHA-256.

El reconstructor semántico es independiente de la proyección: consume solo
JSONL y la referencia autoritativa, valida la forma, el orden, la unicidad,
recuentos, digests y relaciones, y vuelve a materializar las copias y los
genotipos. La UI Marimo muestra únicamente el resumen del catálogo del
manifiesto; no carga el fichero completo.

## Trazabilidad de cierre

- **REQ-VER-31 / ACE-62:** `tests/truth_validation.rs` reconstruye por flujo
  copias y genotipos conocidos, incluidas inserciones y deleciones, sin VCF;
  cubre esquema, orden, relaciones, recuentos y digest de referencia.
- **REQ-SAL-06 / ACE-34:** las mismas pruebas demuestran que omitir VCF no
  impide reconstruir la verdad y rechazan corrupción de esquema, dosis,
  portadores, historial, fragmentos y digest.
- **Evidencia compartida de publicación:** `tests/f2_4_integration.rs` verifica
  ruta, tamaño y SHA-256 de `truth-jsonl`, reconstruye la publicación sin VCF
  y rechaza un manifiesto sin su enlace. `tests_ui/test_truth_summary.py`
  verifica el resumen pequeño y enlaces rotos en UI.

## Validación ejecutada

Se superaron las pruebas focales:

```text
cargo test --locked --test truth_validation --test f2_4_integration
```

También terminaron con código 0:

```text
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --all-features --locked -- -D warnings
cargo test --workspace --all-targets --all-features --locked
cargo test --workspace --doc --all-features --locked
cargo doc --workspace --no-deps --all-features --locked
uv run --frozen python -m unittest tests_ui/test_truth_summary.py
git diff --check
```

El gate Marimo se verificó de forma aislada con límite explícito de 60 s:

```text
$ uv run --frozen marimo --version
0.23.16
$ timeout 60s env MARIMO_SKIP_UPDATE_CHECK=1 uv run --frozen marimo check --strict ui/gps_app.py
$ printf 'marimo_check_exit=%s\n' "$?"
marimo_check_exit=0
```

La reproducción del revisor observó una duración de **0,73 s** para ese gate.

La ejecución necesita acceso a la caché local de `uv` fuera del workspace; sin
ese acceso falla antes de iniciar Marimo con error de sistema de ficheros. No
se generó VCF, complemento VCF, escenarios, servicios ni empaquetado.

## Límite y pausa

F2.4 no implementa VCF ni permite que una vista VCF reconstruya, corrija o
sustituya `truth.jsonl`. La siguiente acción requiere aprobar la pausa F2.4:
no iniciar F2.5 si la verdad no es independiente, reconstruible y derivada
únicamente del estado canónico.
