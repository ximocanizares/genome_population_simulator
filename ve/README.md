# Banco de validación externa (VE)

Herramienta independiente del runtime Rust del simulador. Invoca procesos
locales, lee sus salidas y publica su propio expediente. No se incorpora al
simulador ni instala o descarga herramientas.

La autoridad normativa es
`design/plan-validacion-externa-fastq-vcf-0.1.md`.

## Entorno

`VE.1.1` congela la identidad de cada ejecutable antes de ejecutar el ensayo:
ruta resuelta, tamaño, SHA-256 y versión declarada. El fichero
`environment.lock.json` es esa instantánea.

```bash
uv run --project ve python -m ve.environment freeze    # crea o actualiza el lock
uv run --project ve python -m ve.environment verify    # comprueba contra el lock
```

`verify` termina con código distinto de cero si una herramienta falta, cambia de
ruta, de tamaño, de hash o de versión.
