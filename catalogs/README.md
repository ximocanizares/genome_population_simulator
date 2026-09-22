# Catálogos normativos de derivación v1

Los ficheros `subseed-scopes-v1.tsv`, `subseed-scopes-v2.tsv`,
`subseed-scopes-v3.tsv`, `subseed-scopes-v4.tsv`, `subseed-scopes-v5.tsv`,
`subseed-scopes-v6.tsv`, `subseed-scopes-v7.tsv`, `subseed-scopes-v8.tsv`,
`subseed-scopes-v9.tsv`,
`id-classes-v1.tsv` e
`id-classes-v2.tsv` e `id-classes-v3.tsv` fijan los contratos versionados
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

La v7 conserva literalmente las filas de v6 y añade `technical_template`, con
`event_id`, `library_id` y el ordinal lógico de plantilla. `id-classes-v2.tsv`
conserva las filas v1 y añade `batch`, `library`, `read_group` y `template`.
Batch, biblioteca y grupo de lectura usan BLAKE3 como las demás identidades
automáticas. `TemplateId` es la excepción acotada de DTC-039: se compone como
`TPL_<library_id>_<template_ordinal_decimal>`; como la biblioteca es única en
el evento y el ordinal lo es dentro de la biblioteca, no hay colisión posible
ni se conserva un registro creciente de plantillas. Ninguna causa incorpora
ruta, orden de iteración, hilo o planificación.

`id-classes-v3.tsv` conserva literalmente las filas v2 y añade
`PreHopInstanceId` (`PH_`) y `TechnicalPairId` (`TP_`) con derivación 2. La
instancia propia depende de evento, biblioteca y plantilla; una contaminante,
de evento, biblioteca receptora e instancia fuente; y la pareja final, de
evento, biblioteca final e instancia pre-salto. Ambas clases BLAKE3 detectan
colisiones y no son ámbitos RNG.

La v8 conserva literalmente la v7 y añade `technical_library_count`, derivado
de `event_id` y `library_id`. Es el único flujo que resuelve el número Poisson
de plantillas de una biblioteca WGS; queda separado de los sorteos de origen y
orientación de cada plantilla, y por ello añadir o reordenar bibliotecas no
altera el conteo ya existente.

La v9 conserva literalmente la v8 y añade `initial_repeat_allele`, que sortea
qué alelo de un locus de repetición en tándem le toca a una copia haplotípica.
Sus componentes son identidades estables —escenario, referencia, población,
individuo, **índice de copia** e **identificador del locus**—, y deliberadamente
**no** incluyen la posición del locus en la referencia ni el orden en que
individuos o loci aparecen declarados. Esa elección es la que hace que insertar
un individuo al principio de la lista, reordenarla o añadir un locus en un
cromosoma que ordene antes que los existentes no cambie ningún genotipo ya
sorteado.

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

`diagnostics-v3.tsv` conserva v2 y añade el rechazo F4 de categorías
*paired-end* observables y alineamiento directo, con el destino aplazado en el
diagnóstico. No representa esos parámetros ni publica sus estadísticas o
salidas.

`f5-rng-catalog-v1.tsv` publica los ámbitos aleatorios propios de F5 sin
ampliar `subseed-scopes-v8`. Cada decisión deriva una subsemilla local
independiente con el mismo material longitud-prefijado v1 y ChaCha12, e incluye
las primeras cuatro palabras `u64` literales para congelar el contrato de
inicialización del flujo.

`f5-rng-catalog-v2.tsv` conserva literalmente las doce filas de v1 —material,
subsemilla y las cuatro palabras de cada decisión— y añade una sola:
`block_haplotype_count`, el ámbito del recuento haplotípico de un
bloque×subpoblación. Una ejecución publica el catálogo del que salen **sus**
decisiones: un perfil que no declara la distribución del recuento no puede
abrir ese ámbito y publica v1, con su ruta y su huella intactas.

`id-catalog-v4.tsv` conserva literalmente las once filas de
`id-classes-v3.tsv` y añade las cinco clases automáticas de F5: `BLK_`,
`SITE_`, `HAP_`, `FCP_` y `FHA_`. `FCP_` identifica una copia fundadora por
`sequence_id` y `phase_id`; `FHA_` identifica una asignación de bloque separada
que enlaza `founder_copy_id` con `haplotype_id`. Los catálogos
`id-classes-v1.tsv` a `id-classes-v3.tsv` permanecen inmutables.
