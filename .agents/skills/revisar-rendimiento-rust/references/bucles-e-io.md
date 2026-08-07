# Bucles críticos y entrada/salida

## Bucles

- Sacar del bucle validaciones, conversiones y búsquedas que no cambien por
  elemento.
- Evitar asignaciones, clonaciones, formateo y construcción de errores en la
  ruta normal.
- Reducir búsquedas repetidas en mapas cuando un identificador compacto o una
  tabla indexada sea suficiente.
- Revisar ramas impredecibles en recorridos por base o lectura; no reemplazarlas
  sin perfil.
- Comparar iteradores y bucles explícitos por claridad primero y por medida
  después; Rust puede optimizar ambos.
- Evitar despacho dinámico por registro cuando el conjunto de variantes sea
  cerrado y el perfil muestre su coste.
- Probar vectorización o `unsafe` solo cuando el perfil localice el bucle y
  existan pruebas de equivalencia y límites.

## Generación aleatoria

- No reconstruir generadores ni distribuciones por individuo o lectura si sus
  parámetros permanecen constantes.
- Derivar subsemillas fuera de bucles internos y conservar independencia entre
  tareas.
- Comparar el coste de distribuciones solo con la misma calidad estadística y
  semántica.
- No sacrificar reproducibilidad entre hilos por una mejora de velocidad no
  acordada.

## Archivos

- Usar lectura y escritura con búfer y lotes suficientemente grandes.
- Evitar una llamada al sistema por base, variante o lectura.
- Transmitir registros cuando no necesiten ordenación global; materializar solo
  lo que exija ordenar, indexar o reconstruir.
- Separar tiempo de generación, compresión, ordenación, indexación y disco.
- Medir compresión y nivel elegido con datos representativos de FASTQ, BAM,
  CRAM y VCF.
- Comprobar que el cuello de botella no sea el disco o el proceso externo antes
  de optimizar Rust.

## Procesos externos

- Evitar bloqueos mutuos al escribir entrada y leer salida o `stderr`.
- Drenar flujos, comprobar códigos de salida y limitar búferes.
- Elegir archivo temporal o tubería según memoria, posibilidad de reanudar y
  necesidad de conservar FASTQ.
- No multiplicar hilos internos de minimap2, compresión y Rayon por el número
  total de núcleos simultáneamente.

## Experimento

Medir registros o bases por segundo, bytes escritos, tiempo de CPU y tiempo de
espera. Verificar tamaños, sumas y contenido de salida antes y después.
