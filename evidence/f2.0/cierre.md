# Cierre local de F2.0

Fecha: 21 de agosto de 2026.

## Alcance comprobado

- El lanzador enlaza exclusivamente con `127.0.0.1`, admite un puerto local y
  no ejecuta el simulador al abrir la aplicación.
- Marimo representa todos los campos de la configuración F1 vigente y no
  anticipa escenarios, ploidía ni poblaciones.
- La importación y exportación conservan la semántica, el orden, los IDs, las
  rutas con espacios y los literales decimales normativos.
- `validate` y `run` se invocan como procesos con argumentos separados. Rust
  conserva la autoridad semántica, una edición invalida la validación y el
  destino existente se rechaza sin ejecutar.
- La aplicación muestra estados textuales, diagnósticos, logs acotados, ruta del
  resultado y manifiesto sin cargar artefactos genómicos grandes.
- Una configuración manual y su forma exportada por la UI producen la misma
  `EffectiveConfig`, descontando solo `source_path`, que identifica los dos
  temporales de entrada distintos.

## Evidencia ejecutada

- 21 pruebas Python de autoría, lanzador, adaptador y recorrido transversal:
  superadas con `uv run --frozen --offline python -m unittest discover -s
  tests_ui -v`.
- `marimo check --strict ui/gps_app.py`: superado.
- Exportación HTML de la aplicación: superada sin validar ni ejecutar.
- Arranque temporal en `127.0.0.1:2731`: respuesta HTTP 200 y cierre limpio.
- `cargo fmt`, `cargo clippy -D warnings`, pruebas de todos los targets y
  características, doctests y `cargo doc`: superados.
- `git diff --check`: superado.

## Límites

Este cierre demuestra instalación desde el repositorio y uso local. No afirma
distribución oficial, compatibilidad acreditada de Ubuntu 24.04, paquete
`.tar.gz`, OCI ni auditoría completa de navegador/accesibilidad; esas pruebas
pertenecen a F9.
