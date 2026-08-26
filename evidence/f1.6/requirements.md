# Evidencia vigente de F1.6

La matriz autoritativa es
`evidence/f1-reimplementation/f1-complete-test-req-matrix.tsv`. Esta síntesis
solo nombra símbolos ejecutables actuales; `REQ-CONF-126` permanece fuera del
cierre F1.6 simplificado.

| Caso | Evidencia ejecutable actual |
|---|---|
| TEST-REQ-CONF-124 | `reference_reduction_config::declared_reduction_preserves_target_variant_literal_and_location`. |
| TEST-REQ-CONF-125 | `reference_reduction_config::{exactly_one_target_and_known_fields_are_required,zero_endpoints_source_mismatch_and_non_reduction_are_rejected}`. |
| TEST-REQ-CONF-127 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest` y `effective_config_runtime_validation::every_reduction_field_and_exact_reference_set_are_cross_checked`. |
| TEST-REQ-REF-24 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest`. |
| TEST-REQ-REF-25 | `execution::tests::reduction_selection_uses_reloaded_copy_not_the_stale_origin_object`. |
| TEST-REQ-REF-26 | `reference_reduction::tests::{direct_selector_is_deterministic_non_overlapping_and_has_one_residual,canonical_selector_reproduces_the_frozen_direct_vector}`. |
| TEST-REQ-REF-27 | `exact_decimal::tests::{exact_boundary_does_not_round_through_f64,supports_precision_beyond_fixed_width_integers}`. |
| TEST-REQ-REF-28 | `reference_reduction::tests::{direct_selector_is_deterministic_non_overlapping_and_has_one_residual,exact_multiple_has_no_residual}`. |
| TEST-REQ-REF-29 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest`. |
| TEST-REQ-REF-30 | `f1_6_acceptance::one_block_reduction_omits_unselected_sequences_without_empty_records`. |
| TEST-REQ-REF-31 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest` y `published_resource_validation::every_derived_sequence_catalog_field_is_checked_independently`. |
| TEST-REQ-REF-32 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest` y `published_resource_validation::every_derived_sequence_catalog_field_is_checked_independently`. |
| TEST-REQ-REF-33 | `f1_6_acceptance::{percentage_mode_is_exact_and_infeasible_layout_never_creates_a_derived_artifact,one_invalid_reduction_prevents_any_partial_publication}`. |
| TEST-REQ-REF-34 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest`. |
| TEST-REQ-REF-35 | `coordinates::tests::converts_first_base_and_complete_sequence_exactly` y `f1_6_acceptance::serialized_map_validator_rejects_truncation_overlap_and_out_of_range`. |
| TEST-REQ-REF-36 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest`. |
| TEST-REQ-VAL-52 | `reference_reduction_config::zero_endpoints_source_mismatch_and_non_reduction_are_rejected` y `f1_6_acceptance::one_invalid_reduction_prevents_any_partial_publication`. |
| TEST-REQ-VAL-53 | `f1_6_acceptance::serialized_map_validator_rejects_truncation_overlap_and_out_of_range`. |
| TEST-REQ-MAN-89 | `f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest` y `effective_config_runtime_validation::every_reduction_field_and_exact_reference_set_are_cross_checked`. |
| TEST-REQ-MAN-90 | `f1_6_acceptance::serialized_map_validator_rejects_truncation_overlap_and_out_of_range` y `published_resource_validation::correlated_documents_cannot_replace_physical_authority`. |
| TEST-REQ-SAL-27 | `f1_6_acceptance::{one_invalid_reduction_prevents_any_partial_publication,target_reduction_publishes_reconstructible_fasta_map_and_manifest}` y `execution::tests::physical_resource_failure_marks_reduction_failed_without_partial_artifacts`. |

`tests/f1_6_traceability.rs` exige que la lista ejecutada contenga exactamente
estos testigos vivos, más la propia prueba de trazabilidad.
