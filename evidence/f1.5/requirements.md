# Batería trazable F1.5

Fecha: 10 de agosto de 2026. Cada identificador `TEST-REQ-*` designa un caso
observable y enlaza su prueba ejecutable. «Parcial» no se contabiliza como
superado.

Esta batería conserva la evidencia del checkpoint `66ad8b2`, pero DTC-078
reclasifica los REQ cuyo contrato completo necesita consumidores posteriores.
Solo las filas `Superado` forman la asignación vigente; todas se repetirán antes
de aprobar la pausa F1.5 del plan 0.2.

| Caso | Estado | Evidencia ejecutable |
|---|---|---|
| TEST-REQ-CONF-02 | Superado | `f1_5_acceptance::{simulated_cli_runs_repeat_logically_and_pass_independent_oracles,supplied_cli_preserves_n_and_exact_input_while_invalid_input_has_no_destination}` |
| TEST-REQ-CONF-03 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte` |
| TEST-REQ-CONF-26 | Superado | `f1_5_acceptance::simulated_cli_runs_repeat_logically_and_pass_independent_oracles` conserva raíz declarada en la configuración efectiva real. |
| TEST-REQ-CONF-27 | Parcial | F1 repite referencia y manifiesto mínimo; la simulación biológica completa cierra en F9.7. |
| TEST-REQ-REF-01 | Superado | Los dos recorridos CLI de `f1_5_acceptance` cubren referencia simulada y suministrada. |
| TEST-REQ-REF-02 | Superado | `reference_domain::reference_preserves_sequence_order_origin_and_zero_or_one_f1_index` y `reference_generation::generated_reference_respects_order_lengths_classes_and_degenerate_composition`. |
| TEST-REQ-REF-03 | Superado | `results_publication::complete_zero_one_and_two_event_results_have_resolvable_hashed_links_and_valid_json` comprueba un único FASTA enlazado por referencia. |
| TEST-REQ-REF-04 | Superado | `reference_generation::generated_fasta_v1_has_literal_headers_lf_order_and_eighty_base_lines` cubre dos secuencias en el FASTA único. |
| TEST-REQ-REF-05 | Superado | `reference_index_acceptance::normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand` demuestra el contrato v1 de cero índices normales o un único FAI bajo demanda, consumido por una consulta y validado independientemente; DTF-048 reserva los índices múltiples. |
| TEST-REQ-REF-06 | Superado | La consulta indexada de aceptación crea exactamente el FAI que necesita y ejecuta la consulta conocida. |
| TEST-REQ-REF-07 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte` y los esquemas comprueban origen y ruta suministrada; el recorrido simulado declara `generated`. |
| TEST-REQ-REF-12 | Superado | `fasta_validation::supplied_fasta_preserves_n_case_crlf_and_non_utf8_header_bytes` y rechazo prepublicación. |
| TEST-REQ-REF-13 | Superado | `fasta_validation::supplied_fasta_rejects_empty_file_and_empty_sequence` y `fasta_validation::supplied_fasta_rejects_empty_and_duplicate_ids`. |
| TEST-REQ-REF-14 | Superado | `fasta_validation::supplied_fasta_preserves_n_case_crlf_and_non_utf8_header_bytes`, `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte` y `f1_5_acceptance::supplied_cli_preserves_n_and_exact_input_while_invalid_input_has_no_destination` admiten ambas cajas, conservan N en posición/longitud y copia exacta. |
| TEST-REQ-REF-15 | Superado | `fasta_validation::supplied_fasta_rejects_gzip_gap_iupac_and_other_symbol`. |
| TEST-REQ-REF-16 | Superado | `reference_copy::invalid_supplied_fasta_is_rejected_before_creating_a_destination`. |
| TEST-REQ-REF-17 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte`, además de `cmp` en aceptación. |
| TEST-REQ-REF-18 | Superado | El mismo caso enlaza la copia catalogada desde configuración efectiva, manifiesto raíz y eventos; el cambio posterior a validación no puede catalogarse. |
| TEST-REQ-REF-22 | Superado | Las pruebas de dominio, generación y manifiesto conservan clase por secuencia. |
| TEST-REQ-REF-23 | Parcial | F1 cubre referencias; haplotipos, copias, inserciones, motivos y catálogos cierran en F9.7. |
| TEST-REQ-REP-22 | Parcial | F1 conserva entradas y referencia; perfiles, herramientas y decisiones de todas las fases cierran en F9.7. |
| TEST-REQ-EVO-63 | Parcial | El checkpoint construyó eventos ficticios; un evento real con directorio y manifiesto propios cierra en F3.5. |
| TEST-REQ-VAL-45 | Superado | `fasta_validation` prueba primer token, duplicados y alfabeto antes de publicación con diagnóstico identificable. |
| TEST-REQ-MAN-06 | Superado | Configuración efectiva y manifiesto distinguen `generated`/`supplied` y conservan ruta de origen cuando corresponde. |
| TEST-REQ-MAN-72 | Parcial | F1 conserva el manifiesto general; subsemillas y manifiestos de eventos reales cierran en F3.5. |
| TEST-REQ-MAN-86 | Superado | `reference_copy` y JSON Schema verifican ID, cabecera reversible, descripción, clase, longitud, ploidía, validación y recuento de N. |
| TEST-REQ-SAL-01 | Superado | Recorridos válidos publican exactamente un FASTA por referencia y un FAI solo en aceptación, fuera de la salida normal. |
| TEST-REQ-SAL-18 | Parcial | F1 crea el área común; los enlaces verificables desde eventos reales cierran en F3.5. |

Resultado vigente: 22 REQ completos asignados a F1.5 y 6 evidencias parciales
reasignadas. REQ-REF-05 conserva en DTF-048 los índices múltiples; ninguna fila
parcial se usa para cerrar F1.
