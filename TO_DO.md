# Trabajo pendiente

## Mejoras de rendimiento resueltas

- [x] Eliminar las validaciones físicas redundantes durante la publicación:
  una única validación completa se realiza tras escribir y releer
  `manifest.json`; antes del renombrado atómico se conservan las comprobaciones
  de catálogo, tamaño, SHA-256 y manifiesto.

## Configuración de cohortes grandes

- [ ] Añadir una forma compacta de declarar fundadores consecutivos mediante
  una raíz y numeración determinista (por ejemplo, `founder_id_prefix = "tom"`
  para producir `tom-0001`, `tom-0002`, … hasta `individual_count`).
  - Debe conservar reproducibilidad, unicidad dentro de la población y la
    derivación estable de IDs técnicos.
  - `individual_ids` y `founder_keys` explícitos deben seguir disponibles para
    nombres no consecutivos o asignados externamente.
  - La sintaxis exacta, los conflictos entre modalidades y la interfaz Marimo
    requieren especificación antes de implementarse.

## Salidas FASTQ

- [ ] Permitir publicar los FASTQ finales comprimidos con gzip (`.fastq.gz`).
  - Mantener la escritura, validación, tamaño y SHA-256 sobre los bytes finales
    publicados; no sustituir la validación de registros FASTQ por la del flujo
    comprimido.
  - Definir el nivel y versión de compresión, la reproducibilidad de los bytes,
    la compatibilidad con el validador por flujo y las pruebas de
    interoperabilidad antes de implementarlo.
