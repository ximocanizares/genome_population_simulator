# Catálogos normativos de derivación v1

Los ficheros `subseed-scopes-v1.tsv`, `subseed-scopes-v2.tsv`,
`subseed-scopes-v3.tsv`, `subseed-scopes-v4.tsv`, `subseed-scopes-v5.tsv`,
`subseed-scopes-v6.tsv` e
`id-classes-v1.tsv` fijan los contratos versionados
de la codificación definida por DTC-037 y DTC-039. La versión 3 añade el ámbito
`reference_reduction` sin modificar los vectores anteriores. Sus filas son vectores
literales, no ejemplos orientativos. Cambiar un ámbito, nombre, orden, prefijo,
causa, material o resultado requiere otra versión de derivación.

La v5 es superconjunto estricto de v4: conserva todas sus filas literalmente y
añade `initial_phase_locus` para permutar una composición multialélica completa
con un flujo canónico por locus.

La v6 conserva v5 sin modificar y sustituye únicamente su fila
`reference_reduction` por el selector canónico
`direct_proportional_selection`. Es el catálogo activo del runtime; no promete
compatibilidad con el selector transitorio anterior.

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

`diagnostics-v2.tsv` amplía el catálogo inmutable v1 con los rechazos de F1.6.
Los códigos anteriores conservan significado y contexto; los nuevos códigos se
usan para decimal exacto, objetivo o disposición inviables, mapa incoherente y
coordenadas incompatibles.
