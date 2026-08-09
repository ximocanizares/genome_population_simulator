# Catálogos normativos de derivación v1

Los ficheros `subseed-scopes-v1.tsv` e `id-classes-v1.tsv` fijan la versión 1
de la codificación definida por DTC-037 y DTC-039. Sus filas son vectores
literales, no ejemplos orientativos. Cambiar un ámbito, nombre, orden, prefijo,
causa, material o resultado requiere otra versión de derivación.

Cada entero es `u32` big-endian. Los nombres de campo son ASCII y sus valores
son UTF-8; ambos llevan delante su longitud en bytes. Las rutas, el reloj, los
hilos, el orden de iteración, la planificación y cualquier posición accidental
están prohibidos como entradas.

Los resultados se calcularon el 8 de agosto de 2026 con `blake3` para Python
1.0.9 como oráculo independiente y los contextos literales publicados. Esa
herramienta temporal no forma parte del producto, del grafo Cargo ni del
runtime. Las pruebas reconstruyen byte a byte el material publicado, comprueban
longitudes y enteros y verifican que el Base32 de cada ID representa exactamente
los primeros 128 bits del digest literal.
