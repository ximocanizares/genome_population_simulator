# Batería trazable F1.5

Fecha: 10 de agosto de 2026. Cada identificador `TEST-REQ-*` designa un caso
observable y enlaza su prueba ejecutable. «Parcial» no se contabiliza como
superado.

| Caso | Estado | Evidencia ejecutable |
|---|---|---|
| TEST-REQ-CONF-02 | Superado | `f1_5_acceptance::{simulated_cli_runs_repeat_logically_and_pass_independent_oracles,supplied_cli_preserves_n_and_exact_input_while_invalid_input_has_no_destination}` |
| TEST-REQ-CONF-03 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte` |
| TEST-REQ-CONF-26 | Superado | `f1_5_acceptance::simulated_cli_runs_repeat_logically_and_pass_independent_oracles` conserva raíz declarada en la configuración efectiva real. |
| TEST-REQ-CONF-27 | Superado | El mismo caso repite configuración efectiva, manifiesto y FASTA con la misma raíz, versión y entrada. |
| TEST-REQ-REF-01 | Superado | Los dos recorridos CLI de `f1_5_acceptance` cubren referencia simulada y suministrada. |
| TEST-REQ-REF-02 | Superado | `reference_domain::reference_preserves_sequence_order_classes_lengths_and_origin` y `reference_generation::generated_reference_respects_order_lengths_classes_and_degenerate_composition`. |
| TEST-REQ-REF-03 | Superado | `results_publication::complete_zero_one_and_two_event_results_have_resolvable_hashed_links_and_valid_json` comprueba un único FASTA enlazado por referencia. |
| TEST-REQ-REF-04 | Superado | `reference_generation::generated_fasta_v1_has_literal_headers_lf_order_and_eighty_base_lines` cubre dos secuencias en el FASTA único. |
| TEST-REQ-REF-05 | Superado | `reference_index_acceptance::normal_f1_run_has_zero_indexes_and_acceptance_query_creates_fai_on_demand` demuestra el contrato v1 de cero índices normales o un único FAI bajo demanda, consumido por una consulta y validado independientemente; DTF-048 reserva los índices múltiples. |
| TEST-REQ-REF-06 | Superado | La consulta indexada de aceptación crea exactamente el FAI que necesita y ejecuta la consulta conocida. |
| TEST-REQ-REF-07 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte` y los esquemas comprueban origen y ruta suministrada; el recorrido simulado declara `generated`. |
| TEST-REQ-REF-12 | Superado | `fasta_validation::valid_fasta_preserves_headers_descriptions_n_and_order` y rechazo prepublicación. |
| TEST-REQ-REF-13 | Superado | `fasta_validation::{empty_file_and_empty_sequence_are_rejected,empty_and_duplicate_identifiers_are_rejected}`. |
| TEST-REQ-REF-15 | Superado | `fasta_validation::gzip_gap_iupac_and_other_symbols_are_rejected_with_context`. |
| TEST-REQ-REF-16 | Superado | `reference_copy::invalid_supplied_fasta_is_rejected_before_creating_a_destination`. |
| TEST-REQ-REF-17 | Superado | `reference_copy::supplied_run_copies_the_validated_fasta_byte_for_byte`, además de `cmp` en aceptación. |
| TEST-REQ-REF-18 | Superado | El mismo caso enlaza la copia catalogada desde configuración efectiva, manifiesto raíz y eventos; el cambio posterior a validación no puede catalogarse. |
| TEST-REQ-REF-22 | Superado | Las pruebas de dominio, generación y manifiesto conservan clase por secuencia. |
| TEST-REQ-REF-23 | Superado | Pruebas de dominio/FASTA conservan solo `N` y diagnostican todo símbolo no admitido con posición base 1. |
| TEST-REQ-REP-22 | Superado | Repetición CLI conserva versión del programa, plataforma, toolchain, dependencias, raíz, configuración efectiva, FASTA y manifiesto. |
| TEST-REQ-EVO-63 | Superado | `results_publication::complete_zero_one_and_two_event_results_have_resolvable_hashed_links_and_valid_json` demuestra directorio y manifiesto propios por evento y enlaces desde la raíz a configuración, referencia y subsemilla, sin inventar resultados. REQ-EVO-78 reserva los enlaces a resultados reales para F5.6. |
| TEST-REQ-VAL-45 | Superado | `fasta_validation` prueba primer token, duplicados y alfabeto antes de publicación con diagnóstico identificable. |
| TEST-REQ-MAN-06 | Superado | Configuración efectiva y manifiesto distinguen `generated`/`supplied` y conservan ruta de origen cuando corresponde. |
| TEST-REQ-MAN-72 | Superado | `results_publication` valida enlaces resolubles a configuración efectiva, referencia común, subsemilla y manifiesto de cada evento. |
| TEST-REQ-MAN-86 | Superado | `reference_copy` y JSON Schema verifican ID, cabecera reversible, descripción, clase, longitud, ploidía, validación y recuento de N. |
| TEST-REQ-SAL-01 | Superado | Recorridos válidos publican exactamente un FASTA por referencia y un FAI solo en aceptación, fuera de la salida normal. |
| TEST-REQ-SAL-18 | Superado | El resultado contiene manifiesto raíz, área `resources/references` e índices vacíos enlazados verificablemente desde eventos. |

Resultado: 27 de 27 casos de F1.5 superados. REQ-REF-05 y REQ-EVO-63 cumplen sus
contratos v1; DTF-048 y REQ-EVO-78 conservan, respectivamente, los índices
múltiples y los enlaces a resultados reales fuera de esta batería.
