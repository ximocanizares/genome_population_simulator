# Contrato `truth.jsonl` v1

`truth.jsonl` es UTF-8 sin compresión y cada línea contiene un objeto JSON que
cumple `truth-v1.schema.json`. La primera línea es `header`; las posteriores son
`haplotype_copy` y después `genotype`. La versión es textual (`"1"`) y una
versión no reconocida se rechaza, sin intentar reinterpretarla.

La cabecera fija `gps_truth`, ejecución, referencia, convenio
`zero_based_half_open_v1`, los recuentos declarados y SHA-256 de fuentes y
grupos. Cada digest nombra un par cerrado `kind`/`source_id`; la referencia de
cabecera y toda fuente de copia, segmento o modificación deben resolver de modo
exacto a uno de esos pares. La lista de fuentes no puede estar vacía.
`group_digests` puede estar vacía cuando no hay procedencia de matriz; cada
modificación `matrix` exige en ella el grupo de su `population_id`.

Una copia declara identidad `(reference_id, individual_id, sequence_id,
copy_index)`, fase, procedencia y una cadena no vacía de segmentos. Cada
segmento conserva fuente, intervalo interno `[start,end)`, longitud de fuente y
modificaciones ordenadas y no solapadas de SNP, inserción o deleción. Los
intervalos cumplen `0 <= start < end <= sequence_length`; fuente y copia deben
compartir referencia y secuencia. Un SNP cambia realmente de base y una
deleción no es vacía. La inserción contiene sus únicas bases propias.

Cada modificación declara procedencia tipada y versionada: `explicit`,
`random` o `matrix`, con su fuente. La variante `matrix` además exige
`matrix_source_id`, `population_id` y `variant_id`; conserva por tanto la
distinción de ediciones aunque recaigan sobre el mismo intervalo.

Un genotipo nombra `case_id`, individuo, referencia y locus; enumera alelos,
dosis y cada copia portadora con su índice de alelo. Las relaciones se verifican
en el conjunto: una fase coincide con `copy_index`, cada copia portadora existe,
es del individuo/secuencia del genotipo y no se repite; sus recuentos igualan
las dosis. F2.4.1 no materializa secuencias ni comprueba los alelos contra las
modificaciones: ese reconstructor semántico independiente es F2.4.3.

El orden canónico es: copias por referencia, individuo, secuencia e índice;
después genotipos por referencia, individuo, secuencia, inicio, fin y caso. La
cabecera y los esquemas por línea no sustituyen estas relaciones ni los
recuentos del conjunto. El lector externo del esquema rechaza claves JSON
duplicadas a cualquier profundidad antes de aplicar JSON Schema.
