# Validación ejecutada de F1.4

Fecha: 10 de agosto de 2026.

## Entorno observado

- Ubuntu sobre WSL2, kernel `6.18.33.2-microsoft-standard-WSL2`, arquitectura
  `x86_64`.
- Árbol de trabajo situado en `ext4` (`/dev/sdc`), fuera de `/mnt/c`.
- Rust y Cargo fijados por `rust-toolchain.toml` en 1.97.1; destino release
  `x86_64-unknown-linux-gnu`.
- Validador independiente `python-jsonschema` 4.10.3 y GNU `sha256sum` 9.4.

La repetición exigida por D12 y por la pausa en Ubuntu 24.04 nativo no se ha
ejecutado en este entorno y permanece como condición explícita de aprobación.

## Puertas ejecutadas

Todas las órdenes siguientes terminaron con código 0:

```text
cargo fmt --all -- --check
cargo clippy --workspace --all-targets --locked -- -D warnings
cargo test --workspace --all-targets --locked
cargo test --workspace --doc --locked
cargo doc --no-deps --locked
cargo build --release --locked --target x86_64-unknown-linux-gnu
cargo metadata --locked --format-version 1
cargo tree --locked -e features
cargo audit
git diff --check
```

La suite completa ejecutó 109 pruebas de objetivo y los doctests ejecutaron el
caso `compile_fail` que impide cruzar `EventId` e `IndividualId`. `cargo audit`
cargó 1199 avisos RustSec y examinó las 52 dependencias del lock sin informar
vulnerabilidades.

## Oráculos independientes y recorridos D1–D13

- Los documentos completos, fallidos y de evento se contrastaron con
  `python3 -m jsonschema --base-uri file://... -i ...`, separado de
  `serde_json`. El oráculo rechazó truncado, campo desconocido, clase errónea y
  versión errónea. Tras H-F14-R2 rechazó además propiedades inesperadas,
  requeridos ausentes y tipos o rangos inválidos dentro de programa,
  plataforma, toolchain, dependencias, comprobaciones, entradas, eventos,
  salidas, unidades, diagnósticos, frecuencias controladas, valores efectivos,
  selectores y procedencia. Tras H-F14-R2-R1 rechazó también manifiestos
  sintácticamente válidos pero semánticamente contradictorios: estado general
  y unidad incompatibles, diagnóstico obligatorio ausente o espurio,
  configuración efectiva ausente en `complete`, unidad de evento no completa
  o no exportable, nivel de reproducibilidad omitido/duplicado y catálogo de
  dependencias incompleto mediante duplicación con otra versión.
- Los tamaños y SHA-256 publicados se compararon fichero a fichero con
  `sha256sum`; cambiar un byte o truncar hace fallar la identidad esperada.
- Los cuatro vectores literales A6 pasan con material, prefijo, primeros 128
  bits y Base32 exactos. El ID de ejecución usa raíz y digest canónico de los
  bytes admitidos, no rutas, reloj ni orden físico.
- Se ejecutaron cero, uno y dos eventos directos y un escenario real; las
  copias general y de escenario coinciden byte a byte con lo admitido y todos
  los enlaces quedan bajo la raíz publicada.
- La carrera de dos procesos sobre un destino produce exactamente un éxito y
  un `GPS-IO-003`/salida 5. Un destino previo conserva sus bytes y el doble sin
  `RENAME_NOREPLACE` devuelve `GPS-IO-004`, sin fallback.
- Tras H-F14-R1, un doble sincronizado observa `INCOMPLETE` todavía presente al
  entrar en la publicación atómica exitosa y fallida. El éxito mueve primero
  el directorio con marcador y lo retira después; los fallos de publisher,
  retirada o restauración dejan una unidad marcada e informan el error.
- Un rechazo de configuración termina con salida 2 y sin destino; el proceso
  de fallo controlado termina con salida 1 y manifiesto final `failed`; SIGKILL
  deja únicamente un provisional hermano con `INCOMPLETE`.
- Un temporal truncado `*.gps-tmp` no entra en el catálogo y bloquea la
  publicación. Un directorio sin manifiesto general validado tampoco puede
  publicarse.
- Todos los campos de DTC-049 recuperan exactamente los bits `f64` de
  tolerancia, sumas, deltas y valores antes/después; esto requiere la
  característica `serde_json/float_roundtrip` aprobada en D1.

## Auditoría literal de requisitos

`test-req-matrix.tsv` contiene exactamente las 25 filas inicialmente asignadas
a F1.4: 17 están superadas por consumidores reales del hito y 8 están
reasignadas. No se presentan tipos preparatorios como individuos simulados,
resultados biológicos, referencias comunes ni conjuntos genómicos reales.

Las reasignaciones son REQ-IND-25/26, REQ-EVO-63, REQ-MAN-10/72/83 y
REQ-SAL-17/18; cada fila conserva evidencia parcial y destino consumidor.

La pasada H-F14-R3 confirma también la coherencia documental: F1.5 conserva 27
REQ, F1 conserva 61 filas vigentes y 31 casos están reasignados a F2–F9. D12 y
D13 continúan sin aprobación a falta de la evidencia Ubuntu 24.04 nativa.

## Regresión H-F14-R2-R1

Antes de modificar el esquema, la ejecución aislada de
`schemas_reject_semantically_contradictory_f1_4_documents` terminó con código
101: `python-jsonschema` aceptó el primer corrupto, obtenido al cambiar solo el
estado raíz de `complete` a `failed` y conservar `diagnostics=[]`. Tras cerrar
los invariantes, la misma prueba terminó con código 0; la suite específica
`results_publication` pasó 18/18 y la completa pasó 109/109. Se verificó además
la metavalidez Draft 2020-12 de los cinco esquemas con
`Draft202012Validator.check_schema`. El oráculo rechazó todos los
corruptos independientes. La revisión de colecciones equivalentes confirmó que
`checks` ya representa exactamente sus tres entradas: longitud tres,
`uniqueItems`, nombres limitados a esas tres opciones y estado único `passed`
impiden omitir una mediante duplicación. Las colecciones variables no se
convirtieron artificialmente en catálogos cerrados.
