# Revisión de rendimiento F1.5

- Fecha de actualización: 14 de agosto de 2026.
- Contexto: Ubuntu 24.04 bajo WSL2, x86-64, ext4; construcción release portable
  `x86_64-unknown-linux-gnu`, sin `target-cpu=native`.
- Orden reproducible: `scripts/run-f1-reference-benchmarks.sh A 3` y después
  `scripts/run-f1-reference-benchmarks.sh B 3`.
- Etapas observables: `validation_preflight`, medida mediante el subcomando
  público `validate`, y `complete_run`, medida extremo a extremo. La primera no
  se resta de la segunda porque `run` repite la validación. El tiempo de pared
  usa nanosegundos del reloj y GNU time aporta CPU, RSS y cambios de contexto.
- Temporales: cada recorrido ordinario tiene una ejecución gemela con la misma
  configuración, semilla y binario. `strace` 6.8 con `--seccomp-bpf` retrasa
  100 ms la entrada en `renameat2`; durante esa ventana el script mide
  directamente el provisional completo mediante `du -sb`. La sonda falla si no
  observa `renameat2`, si el pico es cero o si `diff -qr` detecta una diferencia
  entre los árboles publicados por ambas ejecuciones. El tiempo, CPU y RSS
  proceden solo del recorrido ordinario, sin la sobrecarga de la sonda. El
  residuo provisional se mide después y fue 0 en todas las repeticiones. La
  sonda requiere permiso `ptrace`; un sandbox que lo prohíba debe informar el
  fallo y no sustituir la medida por cero.

## Resultado observado

La aceptación A de 100 kbp completó tres repeticiones en 7,33–10,20 ms, con RSS
máximo de 3,34 MiB; su preflight aislado tardó 3,46–4,96 ms. La serie B creció
de 1 a 10 y 100 Mb sin cambio del número de artefactos (cuatro) ni temporales
residuales. En 100 Mb se observaron 2,92–5,45 s, 18,3–34,2 MB/s y
aproximadamente 193,7 MiB de RSS. La repetición de 5,45 s tuvo 58 % de CPU y se
conserva como ruido observado, no se elimina de la línea base.

En cada escenario completo, el pico temporal fue exactamente 11 bytes mayor
que el árbol final: 110 912 bytes en A, 1 022 172 en B1, 10 134 680 en B2 y
101 259 688 en B3. Los 11 bytes corresponden al contenido `incomplete\n` del
marcador obligatorio `INCOMPLETE`, presente antes de la publicación atómica y
retirado después. Esta correspondencia se observó directamente y además se
contrastó con árboles gemelos idénticos; no se obtuvo restando o suponiendo el
tamaño del marcador.

No se fija un umbral de producto. El crecimiento temporal es compatible con una
tendencia aproximadamente lineal en el intervalo, con ruido de planificación
visible. El RSS cercano a dos veces el tamaño de la secuencia se explica por la
representación de dominio completa y la copia propietaria requerida por la API
síncrona del escritor de `noodles-fasta`.

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
