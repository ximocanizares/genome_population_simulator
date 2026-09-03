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

- [x] Publicar los FASTQ finales comprimidos como `.fastq.gz`.
  - El perfil v1 usa gzip/BGZF de `noodles-bgzf 0.47.0`, DEFLATE nivel 6,
    bloques BGZF y `MTIME = 0`; no añade dependencias ni una opción pública de
    compresión. La repetibilidad binaria se limita al mismo simulador, lock y
    plataforma admitida.
  - La escritura es por flujo; se relee y valida el FASTQ descomprimido antes
    de publicar atómicamente. Tamaño y SHA-256 se calculan y verifican sobre
    los bytes gzip finales. Las pruebas incluyen decodificación con `gzip -cd`.

## Controles de error de secuenciación PE

- [ ] Exponer en el TOML del evento técnico F4 escalas de error independientes
  para R1 y R2.
  - Deben admitir el valor `0`, que produce exactamente cero sustituciones de
    secuenciación en ese mate, y conservar `1` como perfil actual.
  - La configuración efectiva y el manifiesto del evento deben registrar los
    valores solicitados, efectivos, oportunidades y errores realizados por
    mate; `error_scale = 0` exige `errors = 0`.
  - Este control es una dependencia de VE-D02: permite un ensayo PE neutro de
    validación externa con contaminación, quimeras e *index hopping* también
    fijados a cero.
