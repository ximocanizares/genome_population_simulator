# Contrato `truth.jsonl` v1

`truth.jsonl` es UTF-8 sin compresión y cada línea contiene un objeto JSON que
cumple `truth-v1.schema.json`. La primera línea es `header`; las posteriores son
`haplotype_copy` y después `genotype`. La versión es textual (`"1"`) y una
versión no reconocida se rechaza, sin intentar reinterpretarla.

La cabecera fija `gps_truth`, ejecución, referencia, convenio
`zero_based_half_open_v1`, los recuentos declarados, `reference_content_sha256`
y SHA-256 de fuentes y grupos. `reference_content_sha256` es el SHA-256
canónico del objeto `Reference` (identidad, orden, identificadores, longitudes
y bases normalizadas); se verifica al reconstruir y no sustituye el digest de
bytes del FASTA u otro artefacto declarado para la fuente `reference` dentro de
`source_digests`. Cada digest de fuente nombra un par cerrado `kind`/`source_id`; la referencia de
cabecera y toda fuente de copia, segmento o modificación deben resolver de modo
exacto a uno de esos pares. La lista de fuentes no puede estar vacía.
La identidad del artefacto de matriz es su huella SHA-256 canónica preparada,
no el `source_id` de una fila: se exige en `source_digests` aun si la matriz
vacía fue autorizada. Cada población preparada exige su digest de grupo; los
`source_id` de filas siguen identificando de forma independiente su
procedencia suministrada o reutilizada y pueden ser heterogéneos.

Una copia declara identidad `(reference_id, individual_id, sequence_id,
copy_index)`, fase, procedencia y una cadena no vacía de segmentos. Cada
segmento conserva fuente, intervalo interno `[start,end)`, longitud de fuente y
modificaciones ordenadas y no solapadas de SNP, inserción o deleción. Los
intervalos cumplen `0 <= start < end <= sequence_length`; fuente y copia deben
compartir referencia y secuencia. Un SNP cambia realmente de base y una
deleción no es vacía. La inserción contiene sus únicas bases propias. Además,
`current_pieces` es la vista autoritativa vigente: cada elemento es `source`
con sus SNP vigentes o `inserted` con sus bases y procedencia. `history`
conserva todos los eventos pequeños en el orden aplicado, incluso un SNP que
una deleción posterior eliminó de los fragmentos vigentes.

Cada modificación declara procedencia tipada y versionada: `explicit`,
`random` o `matrix`, con su fuente. La variante `matrix` además exige
`matrix_source_id`, `population_id` y `variant_id`; conserva por tanto la
distinción de ediciones aunque recaigan sobre el mismo intervalo.

Un genotipo nombra `case_id`, individuo, referencia y locus; enumera alelos,
dosis y cada copia portadora con su índice de alelo. Los casos `snp` usan
alelos `reference`/`snp`; los casos `indel` usan `reference` y una inserción o
deleción tipada, nunca el centinela textual `<deleted>`. El alelo estructural
incluye el evento fuente, los `piece_links` supervivientes y
`final_haplotype_bases`, los bytes finales materializados de la copia. Una
inserción parcialmente eliminada es válida si conserva enlaces no vacíos y
acotados a sus fragmentos vigentes. Las relaciones se verifican
en el conjunto: una fase coincide con `copy_index`, cada copia portadora existe,
es del individuo/secuencia del genotipo y no se repite; sus recuentos igualan
las dosis. F2.4.1 solo fija y valida la estructura, sus dominios y enlaces
locales; F2.4.3 materializará y comprobará independientemente la equivalencia
semántica completa entre esos bytes, los enlaces y las modificaciones.

El orden canónico es: copias por referencia, individuo, secuencia e índice;
después genotipos por referencia, individuo, secuencia, inicio, fin y caso. La
cabecera y los esquemas por línea no sustituyen estas relaciones ni los
recuentos del conjunto. El lector externo del esquema rechaza claves JSON
duplicadas a cualquier profundidad antes de aplicar JSON Schema.
