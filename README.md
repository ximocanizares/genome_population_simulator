# genome_population_simulator
simula un genoma, una estructura poblacional, bams y vcfs de la poblaion

## Interfaz local Marimo

El aprovisionamiento inicial puede descargar las dependencias exactas fijadas
en `uv.lock`:

```bash
uv sync --frozen
```

Después, la aplicación se inicia en `127.0.0.1:2718` y no requiere servicios
remotos:

```bash
uv run --frozen python scripts/run_marimo.py
```

Se puede elegir otro puerto local con `--port` y evitar que se abra el navegador
con `--no-open`. Abrir la interfaz no valida ni ejecuta una simulación.
