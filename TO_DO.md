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
  - Conservarán `1` como perfil actual. El valor `0` no se usará para
    representar ausencia de errores: una probabilidad cero no tiene un Phred
    finito compatible con Q0–Q93.
- [x] Añadir `sequencing_errors = "disabled"` como modo PE neutro. Produce
  exactamente cero sustituciones en ambos mates sin usar `error_scale = 0` y
  registra de forma visible que las calidades FASTQ conservadas no representan
  una probabilidad de error realizada.
  - La configuración efectiva registra los modos solicitado/efectivo y la
    marca de calidad. El manifiesto del evento registra además, por mate,
    `quality_shift = 0`, `error_scale = 1`, `opportunities` y `errors`; estos
    contadores agregados no forman parte de la configuración efectiva y el modo
    `disabled` exige `errors = 0`.
  - Este control deja VE-D02 listo para un ensayo PE neutro que fije también a
    cero contaminación, quimeras e *index hopping*.
