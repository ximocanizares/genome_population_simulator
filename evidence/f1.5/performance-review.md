# Revisión de rendimiento F1.5

- Fecha: 10 de agosto de 2026.
- Contexto: Ubuntu 24.04 bajo WSL2, x86-64, ext4; construcción release portable
  `x86_64-unknown-linux-gnu`, sin `target-cpu=native`.
- Orden reproducible: `scripts/run-f1-reference-benchmarks.sh A 3` y después
  `scripts/run-f1-reference-benchmarks.sh B 3`.
- Etapas observables: `validation_preflight`, medida mediante el subcomando
  público `validate`, y `complete_run`, medida extremo a extremo. La primera no
  se resta de la segunda porque `run` repite la validación. El tiempo de pared
  usa nanosegundos del reloj y GNU time aporta CPU, RSS y cambios de contexto.
- Temporales: se mide en disco el residuo provisional tras cada recorrido
  completo (0 bytes en todas las repeticiones). El pico de bytes temporales no
  está instrumentado y se registra explícitamente como `not_measured`; no se
  sustituye por cero.

## Resultado observado

La aceptación A de 100 kbp completó tres repeticiones en 8,50–8,61 ms, con RSS
máximo de 3,3 MiB; su preflight aislado tardó 3,74–4,33 ms. La serie B creció de 1 a 10 y 100 Mb sin cambio del número de
artefactos (cuatro) ni temporales residuales. En 100 Mb se observaron
3,34–3,45 s, 29,0–29,9 MB/s y aproximadamente 193,7 MiB de RSS.

No se fija un umbral de producto. La tendencia temporal es aproximadamente
lineal en el intervalo observado. El RSS cercano a dos veces el tamaño de la
secuencia se explica por la representación de dominio completa y la copia
propietaria requerida por la API síncrona del escritor de `noodles-fasta`.

## Hallazgo

- **H-F15-PERF-01 — medio, aceptado como línea base:** pico de memoria cercano a
  2x al escribir una secuencia generada. No afecta corrección ni incumple un
  umbral aprobado, porque no existe ninguno. Evitar la copia exigiría cambiar la
  frontera propietaria del escritor o el diseño de generación/manifestación;
  hacerlo en F1.5 sin perfil adicional sería especulativo. Se conserva como
  riesgo medido para una optimización posterior, con 100 Mb como regresión.

## Escenario C

La orden manual/programada es
`scripts/run-f1-reference-benchmarks.sh C 1`. No se ejecutó en esta pasada, de
acuerdo con DTC-046; 1 Gb no forma parte de cada cambio ni constituye promesa de
escala del producto.
