# Estudio de mutación de DT-001

Cada fila estropea **una sola** guarda y comprueba que la prueba que dice
cubrirla **muere**. `sha-antes` y `sha-después` son los primeros dieciséis
dígitos hexadecimales del SHA-256 del fichero mutado: demuestran que la
mutación mutó. La copia de seguridad vive fuera de git —nunca
`git checkout --`— y la restauración se verifica comparando la huella con la
de antes.

Las 132 mutaciones están en `mutaciones.json`, con el texto exacto que cada una
sustituye y la línea del árbol final a la que apunta, y el arnés que las aplica
en `arnes-de-mutacion.py`. **Ninguna guarda sobrevivió a la suya.**

El estudio se hizo en dos tandas porque la revisión adversarial obligó a
enderezar trece rechazos del verificador del VCF que se construían a mano: la
primera tanda cubre las 117 guardas que ya pasaban por una constructora y la
segunda las trece que entraron después.

## Primera tanda: 117 mutaciones

```
id                                                   veredicto          sha-antes          sha-después       
vcf_semantic_validation.rs:174 validate_vcf_against_truth :: FORMAT!=GT fijo -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   607d309e1a7d3177   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:187 validate_vcf_against_truth :: GPS_LOSS en cabecera -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   1c70e6d6d2aabdaa   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:209 validate_vcf_against_truth :: # tras #CHROM -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   c8702cfa9b984f35   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:224 validate_vcf_against_truth :: cabecera de columnas -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   41e3dadab56e0665   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:247 validate_vcf_against_truth :: cabecera no reconocida -> rechazo suprimido MUERE (buena)      aeb09b1e1ac19062   5d57d84631634128   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:256 validate_vcf_against_truth :: número de columnas -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   6f65f089bef1bb97   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:280 validate_vcf_against_truth :: CHROM+POS ausente -> búsqueda relajada a cualquier registro MUERE (buena)      aeb09b1e1ac19062   34e66a57b980c05e   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:380 validate_vcf_against_truth :: faltan registros -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   2a9c62391ed74280   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:1427 parse_gps_case :: GPS_CASE no es el primero -> guarda retirada con valor por defecto MUERE (buena)      aeb09b1e1ac19062   11e2cc7833a23ded   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:1433 parse_gps_case :: GPS_CASE vacío -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   6eef9b4063f17b89   test the_header_and_column_contract_rejects_each_corruption_by_its_own_reason ... FAILED
vcf_semantic_validation.rs:1467 validate_reference_allele :: REF vs FASTA -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   c6838f539888db49   test vcf_semantic_validation::tests::the_reference_column_is_tied_to_the_real_genome ... FAILED
vcf_semantic_validation.rs:1463 validate_reference_allele :: REF vs FASTA -> POS deja de ser 1-based MUERE (buena)      aeb09b1e1ac19062   fe3e3b1ad54af5f8   test vcf_semantic_validation::tests::the_reference_column_is_tied_to_the_real_genome ... FAILED
vcf_semantic_validation.rs:1252 validate_oracle_phases :: fases recompuestas -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   9c27fdc7bd4f52ea   test vcf_semantic_validation::tests::the_recomputed_copy_phases_must_be_contiguous ... FAILED
vcf_semantic_validation.rs:669 oracle_path :: tramos source no colineales -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   80961757b606fe17   test vcf_semantic_validation::tests::the_independent_path_rejects_a_badly_threaded_truth ... FAILED
vcf_semantic_validation.rs:716 oracle_path :: tracto no colineal -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   d2fab4c99e1acab6   test vcf_semantic_validation::tests::the_independent_path_rejects_a_badly_threaded_truth ... FAILED
vcf_semantic_validation.rs:692 oracle_path :: SNP antes del tramo -> resta saturada en vez de rechazo MUERE (buena)      aeb09b1e1ac19062   d58b11ffd664c084   test vcf_semantic_validation::tests::the_independent_path_rejects_a_badly_threaded_truth ... FAILED
vcf_semantic_validation.rs:696 oracle_path :: SNP tras el tramo -> el tramo se agranda en vez de rechazar MUERE (buena)      aeb09b1e1ac19062   b39b125e1dc3cb76   test vcf_semantic_validation::tests::the_independent_path_rejects_a_badly_threaded_truth ... FAILED
vcf_semantic_validation.rs:732 oracle_path :: inserción sin tramo previo -> ancla por defecto MUERE (buena)      aeb09b1e1ac19062   bb7021b1d44c76b7   test vcf_semantic_validation::tests::the_independent_path_rejects_a_badly_threaded_truth ... FAILED
vcf_semantic_validation.rs:817 window :: ventana que parte un tracto -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   981265fe98a0aa8a   test vcf_semantic_validation::tests::a_window_may_not_split_a_realized_repeat_tract ... FAILED
vcf_semantic_validation.rs:901 oracle_normalize_window :: indel sin ancla izquierda -> resta saturada MUERE (buena)      aeb09b1e1ac19062   4c8ebf7bd5f1bb6c   test vcf_semantic_validation::tests::an_insertion_at_the_very_start_has_no_left_anchor ... FAILED
vcf_semantic_validation.rs:414 expected_records :: secuencia FASTA ausente -> se toma la primera secuencia MUERE (buena)      aeb09b1e1ac19062   918824ca7d4f582d   test vcf_semantic_validation::tests::the_expected_records_do_not_trust_the_truth_they_read ... FAILED
vcf_semantic_validation.rs:566 expected_records :: ALT sin caso fuente -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   1a293ffd62aba80a   test vcf_semantic_validation::tests::the_expected_records_do_not_trust_the_truth_they_read ... FAILED
vcf_semantic_validation.rs:502 expected_records :: mismo ALT con tractos distintos -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   cac7368f4297be92   test vcf_semantic_validation::tests::the_expected_records_do_not_trust_the_truth_they_read ... FAILED
vcf_semantic_validation.rs:516 expected_records :: una copia con ALT incompatibles -> condición neutralizada MUERE (buena)      aeb09b1e1ac19062   96058204c0638b9d   test vcf_semantic_validation::tests::the_expected_records_do_not_trust_the_truth_they_read ... FAILED
f5_closure.rs:874 walk_founders :: copy_index no entero -> se toma cero MUERE (buena)      143e4486a63b7de1   5cc863cce956eaf8   test a_copy_index_that_is_not_an_integer_is_rejected ... FAILED
f5_closure.rs:886 walk_founders :: FCP_ del TSV sin sintaxis tipada -> análisis relajado MUERE (buena)      143e4486a63b7de1   9e39988d6cdbc5d0   test a_founder_copy_identifier_with_forbidden_syntax_is_rejected ... FAILED
f5_closure.rs:1324 read_line :: lectura no UTF-8 -> el error de lectura se traga MUERE (buena)      143e4486a63b7de1   2cb5512b51a7ecbb   test a_founders_row_that_is_not_utf8_is_rejected ... FAILED
f5_closure.rs:1163 walk_truth_membership :: cabecera de verdad ilegible -> cabecera por defecto MUERE (buena)      143e4486a63b7de1   0aa62bd593c6c91e   test an_unreadable_truth_header_is_rejected ... FAILED
f5_closure.rs:1191 walk_truth_membership :: registro de verdad ilegible -> registro por defecto MUERE (buena)      143e4486a63b7de1   1007840772300cd3   test an_unreadable_truth_record_is_rejected ... FAILED
f5_closure.rs:1232 walk_truth_membership :: FCP_ de la verdad sin sintaxis tipada -> análisis relajado MUERE (buena)      143e4486a63b7de1   a324b103a23d4189   test a_truth_membership_identifier_with_forbidden_syntax_is_rejected ... FAILED
f5_oracle.rs:1723 read_truth_header :: cabecera ilegible -> JSON nulo en vez de rechazo MUERE (buena)      89d1ac587e809ad8   73b0c61e17e00b29   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1649 exact_object :: no es un objeto -> mensaje alterado MUERE (buena)      89d1ac587e809ad8   d48f6433461b3d05   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1650 exact_object :: claves distintas -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   7a5f21bec2a4200a   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1743 read_truth_header :: record_type -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   adc6596fcb229b79   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1746 read_truth_header :: document_kind -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   9ba2c066a0ee0157   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1756 read_truth_header :: convenio de coordenadas -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   f907d1fe44875bdf   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1793 read_truth_header :: populations no es lista -> mensaje alterado MUERE (buena)      89d1ac587e809ad8   99c244bdf2419c1b   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1794 read_truth_header :: populations vacía -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   f4d9e3cd41284a03   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1810 read_truth_header :: populations desordenada -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   f9794cbb7df63d3d   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1822 read_truth_header :: source_digests no es lista -> mensaje alterado MUERE (buena)      89d1ac587e809ad8   b1128ed5ef982afe   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1823 read_truth_header :: source_digests vacía -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   29cbcff7141d4dfd   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1839 read_truth_header :: sha256 de fuente -> condición neutralizada MUERE (buena)      89d1ac587e809ad8   fe9ab7afd7481b36   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1716 matrix_source_id :: matriz de la cabecera -> se acepta cualquier fuente MUERE (buena)      89d1ac587e809ad8   65d0e59073b732a5   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1660 text :: text no es una cadena -> cadena vacía por defecto MUERE (buena)      89d1ac587e809ad8   a6fc7834afc1af56   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1668 number :: number no es entero -> cero por defecto MUERE (buena)      89d1ac587e809ad8   cfb8dbe58cb56e17   test the_truth_header_rejects_each_corruption_by_its_own_reason ... FAILED
f5_oracle.rs:1680 single_base :: "{key} no es una base canónica" MUERE (buena)      89d1ac587e809ad8   200b594c1a717295   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:1926 read_snp_provenance :: "un SNP fundador no procede de una matriz" MUERE (buena)      89d1ac587e809ad8   f40a85a782efbdb3   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:1932 read_snp_provenance :: "versión de procedencia ajena" MUERE (buena)      89d1ac587e809ad8   dc25a2cce2ba26eb   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:1949 read_snp_provenance :: "la fuente del SNP no es la matriz que la cabecera declara MUERE (buena)      89d1ac587e809ad8   1a2c442cc4d73e91   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:1964 read_snp_provenance :: "la procedencia no nombra el sitio de su evento" MUERE (buena)      89d1ac587e809ad8   fc1ae88855d96eb3   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2006 read_segment_source :: "el segmento no pertenece a la copia" MUERE (buena)      89d1ac587e809ad8   cc5ea74837f56d00   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2015 read_segment_source :: "el segmento no cubre la secuencia completa" MUERE (buena)      89d1ac587e809ad8   66fdc8d9d003bea1   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2062 read_truth_copy :: "la pertenencia nombra {population_id}, ausente de la cabe MUERE (buena)      89d1ac587e809ad8   68edf754d159befc   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2077 read_truth_copy :: "origin ajeno a F5" MUERE (buena)      89d1ac587e809ad8   d1b3f0e9162682fa   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2082 read_truth_copy :: "sources no es una lista" MUERE (buena)      89d1ac587e809ad8   13dcf8e63d29c570   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2100 read_truth_copy :: "las fuentes {declared_sources:?} no son las {expected_sou MUERE (buena)      89d1ac587e809ad8   5cbf748f61b18d7b   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2109 read_truth_copy :: "segments no es una lista" MUERE (buena)      89d1ac587e809ad8   ed88e06e90a28aff   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2131 read_truth_copy :: "modifications no es una lista" MUERE (buena)      89d1ac587e809ad8   f593a7be81ada5d3   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2148 read_truth_copy :: "F5 solo publica modificaciones SNP" MUERE (buena)      89d1ac587e809ad8   208389b0c4a508e4   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2164 read_truth_copy :: "un SNP publica el mismo alelo dos veces" MUERE (buena)      89d1ac587e809ad8   787fa0178618aaea   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2178 read_truth_copy :: "las modificaciones no llegan por posición ascendente" MUERE (buena)      89d1ac587e809ad8   9373e652801afeb8   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2198 read_truth_copy :: "current_pieces no es una lista" MUERE (buena)      89d1ac587e809ad8   22eec719a09b35e6   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2209 read_truth_copy :: "el fragmento vigente no procede de la referencia" MUERE (buena)      89d1ac587e809ad8   551056e9bc852055   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2223 read_truth_copy :: "el fragmento vigente no cubre la misma secuencia" MUERE (buena)      89d1ac587e809ad8   f1113333020a7fd6   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2229 read_truth_copy :: "snps no es una lista" MUERE (buena)      89d1ac587e809ad8   91b427ea06436526   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2233 read_truth_copy :: "el fragmento vigente y las modificaciones no traen los mi MUERE (buena)      89d1ac587e809ad8   5cbed4b51ad39467   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2273 read_truth_copy :: "history no es una lista" MUERE (buena)      89d1ac587e809ad8   86af596ae1d60810   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2287 read_truth_copy :: "F5 solo publica hitos SNP" MUERE (buena)      89d1ac587e809ad8   33257b4c4381040b   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2304 read_truth_copy :: "el historial no repite la procedencia de su modificación" MUERE (buena)      89d1ac587e809ad8   20723b5f61059faa   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2519 recount_f5_man70 :: "las subpoblaciones no declaran la misma ploidía" MUERE (buena)      89d1ac587e809ad8   64f7b6e6595ebfeb   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2523 recount_f5_man70 :: "la ploidía declarada no cabe en u8" MUERE (buena)      89d1ac587e809ad8   4a9c76a6d1c08c5d   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2525 recount_f5_man70 :: "la ploidía declarada es cero" MUERE (buena)      89d1ac587e809ad8   6005e432188683d2   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2571 recount_f5_man70 :: "registro sin record_type" MUERE (buena)      89d1ac587e809ad8   5833d64ccdb6d3f8   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2577 recount_f5_man70 :: "una copia aparece después de un genotipo" MUERE (buena)      89d1ac587e809ad8   90a36f72c0d49ced   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2591 recount_f5_man70 :: "las copias no llegan por referencia, población, individuo MUERE (buena)      89d1ac587e809ad8   3b6a42e437888aa8   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2602 recount_f5_man70 :: "un founder_copy_id se repite" MUERE (buena)      89d1ac587e809ad8   0fa40398c1db1428   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2640 recount_f5_man70 :: "los genotipos no llegan por subpoblación y bloque ascende MUERE (buena)      89d1ac587e809ad8   346427b76cd1e92c   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2651 recount_f5_man70 :: "dentro de un bloque los genotipos no llegan por site_id e MUERE (buena)      89d1ac587e809ad8   34b6c6b124f9deb8   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2717 recount_f5_man70 :: "record_type desconocido {other:?}" MUERE (buena)      89d1ac587e809ad8   f68e30d201dc5015   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2825 recount_genotype :: "F5 solo publica casos SNP" MUERE (buena)      89d1ac587e809ad8   724bfbf337c6161a   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2831 recount_genotype :: "un genotipo nombra otra referencia que la cabecera" MUERE (buena)      89d1ac587e809ad8   dfe12762ea54bc75   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2849 recount_genotype :: "el locus SNP no mide una base" MUERE (buena)      89d1ac587e809ad8   27d8be77762cf7cb   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2870 recount_genotype :: "f5_locus repite el alelo" MUERE (buena)      89d1ac587e809ad8   f8e4756a29d96975   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2876 recount_genotype :: "subpoblación {population_id} no declarada" MUERE (buena)      89d1ac587e809ad8   03a2fb4317de5280   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2884 recount_genotype :: "alleles no es una lista" MUERE (buena)      89d1ac587e809ad8   05fe17131d299095   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2890 recount_genotype :: "alelo sin kind" MUERE (buena)      89d1ac587e809ad8   9dc2f5660c17d5e4   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2900 recount_genotype :: "el alelo tipado no es el ALT del locus" MUERE (buena)      89d1ac587e809ad8   6486facdd7a8591e   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2906 recount_genotype :: "dos alelos alternativos en un sitio bialélico" MUERE (buena)      89d1ac587e809ad8   c30e30c821b5e5b9   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2913 recount_genotype :: "alelo {other:?} ajeno a F5" MUERE (buena)      89d1ac587e809ad8   929a4ada37eaf35b   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2921 recount_genotype :: "el sitio no publica su alelo alternativo" MUERE (buena)      89d1ac587e809ad8   fe10e1d18504765c   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2928 recount_genotype :: "allele_dosages no es una lista" MUERE (buena)      89d1ac587e809ad8   1031818e052a2ea8   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2933 recount_genotype :: "dosis no entera" MUERE (buena)      89d1ac587e809ad8   7e6add4c861209fc   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2939 recount_genotype :: "las dosis no acompañan a los alelos" MUERE (buena)      89d1ac587e809ad8   2f49d1bc315e4569   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2950 recount_genotype :: "carriers no es una lista" MUERE (buena)      89d1ac587e809ad8   afbcce4a6f2f612d   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2960 recount_genotype :: "allele_index fuera de la lista de alelos" MUERE (buena)      89d1ac587e809ad8   ae1b2fe6f163ad8e   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2977 recount_genotype :: "un portador no es una copia de su individuo" MUERE (buena)      89d1ac587e809ad8   879a3c58cf944f53   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2981 recount_genotype :: "un portador no tiene copia publicada" MUERE (buena)      89d1ac587e809ad8   e7adc7dd986df9d2   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:3011 recount_genotype :: "los portadores no son las copias de la ploidía" MUERE (buena)      89d1ac587e809ad8   7f41932c9b9caa5c   test each_truth_genotype_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2274 read_truth_copy :: historial de otra longitud MUERE (buena)      89d1ac587e809ad8   24cde74fb6f19c0a   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2292 read_truth_copy :: hito que no casa con su modificación MUERE (buena)      89d1ac587e809ad8   5aed52998f405ba2   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2493 recount_f5_man70 :: documento de verdad ausente MUERE (buena)      89d1ac587e809ad8   7d78f85f45714b63   test the_recount_rejects_the_document_as_bytes ... FAILED
f5_oracle.rs:2505 recount_f5_man70 :: documento vacío MUERE (buena)      89d1ac587e809ad8   75054cd3ded986a0   test the_recount_rejects_the_document_as_bytes ... FAILED
f5_oracle.rs:2507 recount_f5_man70 :: cabecera sin LF MUERE (buena)      89d1ac587e809ad8   0d8885b1bfc3e509   test the_recount_rejects_the_document_as_bytes ... FAILED
f5_oracle.rs:2564 recount_f5_man70 :: línea vacía en el cuerpo MUERE (buena)      89d1ac587e809ad8   3f44569a537c4f3d   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2567 recount_f5_man70 :: registro ilegible en el cuerpo MUERE (buena)      89d1ac587e809ad8   273c20077f345b88   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2595 recount_f5_man70 :: copy_index fuera de la ploidía MUERE (buena)      89d1ac587e809ad8   28d43c92f7e6e17b   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2722 recount_f5_man70 :: documento sin LF final MUERE (buena)      89d1ac587e809ad8   2e18f8b8f7a85ffe   test the_recount_rejects_the_document_as_bytes ... FAILED
f5_oracle.rs:2113 read_truth_copy :: dos segmentos en una copia MUERE (buena)      89d1ac587e809ad8   d536f8a27d9fa2ee   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2202 read_truth_copy :: dos fragmentos vigentes MUERE (buena)      89d1ac587e809ad8   a241c363cfb87c1c   test each_truth_copy_rejects_its_corruptions_by_their_own_reason ... FAILED
f5_oracle.rs:2634 recount_f5_man70 :: grupo repetido MUERE (buena)      89d1ac587e809ad8   90243167572e1f04   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2702 recount_f5_man70 :: sitio que cambia de forma MUERE (buena)      89d1ac587e809ad8   a2baf71221b81abe   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2742 recount_f5_man70 :: recuento del locus contra los individuos declarados MUERE (buena)      89d1ac587e809ad8   23d8b35c1c22832f   test the_recount_rejects_each_document_level_corruption ... FAILED
f5_oracle.rs:2501 recount_f5_man70 :: cabecera con codificación inválida (mensaje sustituido entero) MUERE (buena)      89d1ac587e809ad8   e9dac990fb9fdf32   test the_recount_rejects_the_document_as_bytes ... FAILED
f5_oracle.rs:2551 recount_f5_man70 :: cuerpo con codificación inválida (mensaje sustituido entero) MUERE (buena)      89d1ac587e809ad8   fc57f8b8198093fe   test the_recount_rejects_the_document_as_bytes ... FAILED
vcf_semantic_validation.rs:254 validate_vcf_against_truth :: un rechazo construido a mano vuelve a colarse MUERE (buena)      aeb09b1e1ac19062   a0717c294887eea0   test the_vcf_oracle_builds_every_rejection_through_its_constructors ... FAILED
vcf_semantic_validation.rs:1334 oracle_error :: una constructora pierde la anotación que la ata al barrido MUERE (buena)      aeb09b1e1ac19062   0398a0b41e7d6545   test the_vcf_oracle_builds_every_rejection_through_its_constructors ... FAILED
vcf_semantic_validation.rs:1344 io_reject :: una constructora deja de anotar su disparo MUERE (buena)      aeb09b1e1ac19062   b36f9f8ea0ce5602   test the_vcf_oracle_builds_every_rejection_through_its_constructors ... FAILED

total=117 fallidas=0
```

## Segunda tanda: 13 mutaciones

```
id                                                   veredicto          sha-antes          sha-después       
la verdad no se puede abrir -> se reporta como otra clase MUERE (buena)      93f75fa30b833e0e   50c09b7c658c3060   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
el VCF comprimido no se puede abrir -> otra clase    MUERE (buena)      93f75fa30b833e0e   cb08502cd6dbe20c   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
el VCF de texto no se puede abrir -> otra clase      MUERE (buena)      93f75fa30b833e0e   13d6a186a823f1b6   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
el VCF no es texto -> se culpa a otro fichero        MUERE (buena)      93f75fa30b833e0e   6651d958b697fa29   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
registro antes de #CHROM -> otro motivo              MUERE (buena)      93f75fa30b833e0e   ffa1b4f2ac7669fb   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
POS que no es entero positivo -> otro motivo         MUERE (buena)      93f75fa30b833e0e   a30b35d2badaedcd   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
CHROM ausente del FASTA en la atadura de la columna REF -> otro motivo MUERE (buena)      93f75fa30b833e0e   992a1bdaa638f027   test vcf_semantic_validation::tests::the_reference_column_is_tied_to_the_real_genome ... FAILED
el FASTA no se puede abrir -> otra clase             MUERE (buena)      93f75fa30b833e0e   96008cee90530b9c   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
el FASTA no es texto -> otra clase                   MUERE (buena)      93f75fa30b833e0e   31c3cac39436a48b   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
bases antes de la primera cabecera -> guarda retirada MUERE (buena)      93f75fa30b833e0e   638e4f19fe0019c7   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
conjunto de secuencias inválido -> otra clase        MUERE (buena)      93f75fa30b833e0e   8ef96c32808e14a9   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
cabecera FASTA ilegible -> otra clase                MUERE (buena)      93f75fa30b833e0e   ffc5841725a66c34   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED
bases que no son ADN -> otra clase                   MUERE (buena)      93f75fa30b833e0e   9ac34e86f852dc10   test the_path_oracle_rejects_missing_files_broken_encodings_and_broken_fasta ... FAILED

total=13 fallidas=0
```

## Tercera tanda: 2 mutaciones

```
id                                                   veredicto          sha-antes          sha-después       
el barrido ofrece una segunda función pública        MUERE (buena)      996113dcdbb18b34   bfd8143eaab6b84f   test the_closure_audit_cannot_derive_an_identity_by_its_signature ... FAILED
la clausura nombra el barrido por una segunda vía    MUERE (buena)      143e4486a63b7de1   2c02dfb33ed18ad2   test the_closure_audit_cannot_derive_an_identity_by_its_signature ... FAILED

total=2 fallidas=0
```
