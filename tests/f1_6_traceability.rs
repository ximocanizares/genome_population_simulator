use std::collections::BTreeSet;

type Witness = (&'static str, &'static str);

fn executable_catalog() -> [(&'static str, &'static [Witness]); 22] {
    const CONFIG: &str = include_str!("reference_reduction_config.rs");
    const ACCEPTANCE: &str = include_str!("f1_6_acceptance.rs");
    const COORDINATES: &str = include_str!("../src/coordinates.rs");
    const DECIMAL: &str = include_str!("../src/exact_decimal.rs");
    const REDUCTION: &str = include_str!("../src/reference_reduction.rs");
    const EXECUTION: &str = include_str!("../src/execution.rs");
    const EFFECTIVE: &str = include_str!("effective_config.rs");
    const PUBLICATION: &str = include_str!("../src/publication.rs");
    [
        (
            "CONF-124",
            &[(
                "reference_reduction_config::declared_reduction_preserves_target_variant_literal_and_location",
                CONFIG,
            )],
        ),
        (
            "CONF-125",
            &[(
                "reference_reduction_config::exactly_one_target_and_known_fields_are_required",
                CONFIG,
            )],
        ),
        (
            "CONF-126",
            &[
                (
                    "reference_reduction_config::effective_reduction_uses_exact_decimal_and_own_subseed",
                    CONFIG,
                ),
                (
                    "f1_6_acceptance::omitted_sequence_and_unrelated_consumers_do_not_change_the_reduction_stream",
                    ACCEPTANCE,
                ),
            ],
        ),
        (
            "CONF-127",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-24",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-25",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-26",
            &[(
                "reference_reduction::tests::safe_candidate_groups_match_exhaustive_oracle_for_small_gaps",
                REDUCTION,
            )],
        ),
        (
            "REF-27",
            &[
                (
                    "exact_decimal::tests::exact_boundary_does_not_round_through_f64",
                    DECIMAL,
                ),
                (
                    "exact_decimal::tests::enforces_toml_decimal_grammar_before_removing_separators",
                    DECIMAL,
                ),
            ],
        ),
        (
            "REF-28",
            &[(
                "reference_reduction::tests::creates_exact_non_overlapping_plan_with_one_residual",
                REDUCTION,
            )],
        ),
        (
            "REF-29",
            &[(
                "reference_reduction::tests::materialization_uses_fasta_order_and_accumulates_per_sequence",
                REDUCTION,
            )],
        ),
        (
            "REF-30",
            &[(
                "f1_6_acceptance::omitted_sequence_and_unrelated_consumers_do_not_change_the_reduction_stream",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-31",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-32",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-33",
            &[(
                "f1_6_acceptance::percentage_mode_is_exact_and_infeasible_layout_never_creates_a_derived_artifact",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-34",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "REF-35",
            &[
                (
                    "coordinates::tests::converts_first_base_and_complete_sequence_exactly",
                    COORDINATES,
                ),
                (
                    "f1_6_acceptance::serialized_map_validator_rejects_truncation_overlap_and_out_of_range",
                    ACCEPTANCE,
                ),
            ],
        ),
        (
            "REF-36",
            &[(
                "f1_6_acceptance::target_reduction_publishes_reconstructible_fasta_map_and_manifest",
                ACCEPTANCE,
            )],
        ),
        (
            "VAL-52",
            &[(
                "reference_reduction_config::zero_endpoints_source_mismatch_and_non_reduction_are_rejected",
                CONFIG,
            )],
        ),
        (
            "VAL-53",
            &[(
                "f1_6_acceptance::serialized_map_validator_rejects_truncation_overlap_and_out_of_range",
                ACCEPTANCE,
            )],
        ),
        (
            "MAN-89",
            &[(
                "execution::tests::reduction_fault_matrix_preserves_prior_complete_and_fails_target_and_later_units",
                EXECUTION,
            )],
        ),
        (
            "MAN-90",
            &[
                (
                    "execution::tests::semantic_manifest_gate_rejects_cross_collection_contradictions_from_reloaded_bytes",
                    EXECUTION,
                ),
                (
                    "execution::tests::global_manifest_links_must_match_catalog_exactly_without_reductions",
                    EXECUTION,
                ),
                (
                    "execution::tests::global_manifest_links_match_catalog_for_complete_failed_zero_one_and_two_reductions",
                    EXECUTION,
                ),
                (
                    "execution::tests::manifest_graph_rejects_every_dynamic_edge_contradiction",
                    EXECUTION,
                ),
                (
                    "execution::tests::manifest_dynamic_graph_rejects_diagnostic_and_reference_inconsistencies",
                    EXECUTION,
                ),
                (
                    "execution::tests::manifest_bundle_closes_exact_percent_seed_header_normalization_and_event_edges",
                    EXECUTION,
                ),
                (
                    "execution::tests::manifest_bundle_closes_toml_fasta_and_effective_value_edges",
                    EXECUTION,
                ),
                (
                    "execution::tests::resolution_evidence_toctou_is_rejected_before_publication",
                    EXECUTION,
                ),
                (
                    "effective_config::request_context_drives_combined_dominance_independently_of_rule_order",
                    EFFECTIVE,
                ),
                (
                    "effective_config::foreign_object_and_incomparable_maxima_are_rejected_before_snapshot",
                    EFFECTIVE,
                ),
                (
                    "effective_config::parameter_axes_reject_biological_rule_in_a_technical_event_context",
                    EFFECTIVE,
                ),
                (
                    "effective_config::complete_population_and_temporal_hierarchies_reach_ruleset_dominance",
                    EFFECTIVE,
                ),
                (
                    "effective_config::complete_locus_position_hierarchy_reaches_ruleset_dominance",
                    EFFECTIVE,
                ),
                (
                    "effective_config::orphan_and_incoherent_spatial_coordinates_are_rejected",
                    EFFECTIVE,
                ),
                (
                    "effective_config::spatial_level_and_primary_object_follow_parameter_metadata",
                    EFFECTIVE,
                ),
                (
                    "effective_config::every_parameter_is_checked_at_sequence_region_and_locus_resolution_levels",
                    EFFECTIVE,
                ),
                (
                    "effective_config::hierarchical_contexts_are_not_inferred_from_flat_identifier_lists",
                    EFFECTIVE,
                ),
                (
                    "effective_config::hierarchical_contexts_must_match_parent_interval_and_coordinate_evidence",
                    EFFECTIVE,
                ),
                (
                    "effective_config::phase_ids_are_globally_unique_across_scenarios_independently_of_order",
                    EFFECTIVE,
                ),
                (
                    "effective_config::parameter_order_v1_is_literal_and_independent_of_rule_and_request_permutations",
                    EFFECTIVE,
                ),
            ],
        ),
        (
            "SAL-27",
            &[
                (
                    "publication::tests::intermediate_directory_symlink_is_rejected_before_publisher_and_preserves_incomplete",
                    PUBLICATION,
                ),
                (
                    "publication::tests::canonical_relative_path_rejects_every_textual_alias_directly",
                    PUBLICATION,
                ),
                (
                    "execution::tests::manifest_graph_rejects_every_dynamic_edge_contradiction",
                    EXECUTION,
                ),
                (
                    "execution::tests::manifest_bundle_closes_toml_fasta_and_effective_value_edges",
                    EXECUTION,
                ),
            ],
        ),
    ]
}

fn status_count(tsv: &str, status: &str) -> usize {
    tsv.lines()
        .skip(1)
        .filter(|line| line.split('\t').nth(2) == Some(status))
        .count()
}

#[test]
fn f1_6_catalog_maps_22_requirements_to_existing_executable_tests() {
    let matrix = include_str!("../evidence/f1.6/requirements.md");
    let executed = include_str!("../evidence/f1.6/executed-test-list.txt");
    assert!(executed.contains("total_tests=196"));
    let rows = matrix
        .lines()
        .filter(|line| line.starts_with("| TEST-REQ-"))
        .collect::<Vec<_>>();
    assert_eq!(rows.len(), 22);
    assert_eq!(rows.iter().copied().collect::<BTreeSet<_>>().len(), 22);
    assert!(rows.iter().all(|row| row.contains("| Superado |")));

    let expected = BTreeSet::from([
        "CONF-124", "CONF-125", "CONF-126", "CONF-127", "REF-24", "REF-25", "REF-26", "REF-27",
        "REF-28", "REF-29", "REF-30", "REF-31", "REF-32", "REF-33", "REF-34", "REF-35", "REF-36",
        "VAL-52", "VAL-53", "MAN-89", "MAN-90", "SAL-27",
    ]);
    let observed = rows
        .iter()
        .map(|row| {
            row.strip_prefix("| TEST-REQ-")
                .and_then(|rest| rest.split(" | ").next())
                .expect("TEST-REQ row")
        })
        .collect::<BTreeSet<_>>();
    assert_eq!(observed, expected);

    let catalog = executable_catalog();
    let catalog_ids = catalog
        .iter()
        .map(|(requirement, _)| *requirement)
        .collect::<BTreeSet<_>>();
    assert_eq!(catalog_ids, expected);
    assert_eq!(catalog_ids.len(), catalog.len());
    for (requirement, witnesses) in catalog {
        assert!(
            !witnesses.is_empty(),
            "{requirement} has no executable witness"
        );
        let row = rows
            .iter()
            .find(|row| row.starts_with(&format!("| TEST-REQ-{requirement} |")))
            .expect("requirement evidence row");
        for (qualified_name, source) in witnesses {
            let function_name = qualified_name
                .rsplit("::")
                .next()
                .expect("qualified test name");
            assert!(
                source.contains(&format!("fn {function_name}(")),
                "{requirement} names absent test {qualified_name}"
            );
            assert!(
                row.contains(qualified_name),
                "{requirement} evidence does not cite exact executable {qualified_name}"
            );
            assert!(
                executed.lines().any(|line| line == *qualified_name),
                "{requirement} witness absent from executed listing: {qualified_name}"
            );
        }
    }
}

#[test]
fn complete_f1_evidence_balances_exactly_83_passed_requirements() {
    let f12 = include_str!("../evidence/f1.2/test-req-matrix.tsv");
    let f13 = include_str!("../evidence/f1.3/test-req-matrix.tsv");
    let f14 = include_str!("../evidence/f1.4/test-req-matrix.tsv");
    let f15 = include_str!("../evidence/f1.5/requirements.md");
    let f16 = include_str!("../evidence/f1.6/requirements.md");
    let passed = status_count(f12, "pass")
        + status_count(f13, "pass")
        + status_count(f14, "passed")
        + f15
            .lines()
            .filter(|line| line.contains("| Superado |"))
            .count()
        + f16
            .lines()
            .filter(|line| line.contains("| Superado |"))
            .count();
    assert_eq!(passed, 83);
}

#[test]
fn approved_coordinate_and_selection_decisions_remain_linked() {
    let decisions = include_str!("../design/decisiones-tecnicas-0.1.md");
    let plan = include_str!("../design/plan-implementacion-f1-0.1.md");
    for id in ["DTPD-101", "DTPD-116", "DTC-056", "DTC-057"] {
        assert!(decisions.contains(id), "missing decision {id}");
        assert!(plan.contains(id), "missing plan link {id}");
    }
    assert!(decisions.contains("DTPD-101. Representación interna común de coordenadas (resuelta)"));
    assert!(decisions.contains(
        "DTPD-116. Selección reproducible de intervalos para referencia reducida (resuelta)"
    ));
}

#[test]
fn f1_6_expedient_claims_match_the_executed_independent_semantic_corpus() {
    let requirements = include_str!("../evidence/f1.6/requirements.md");
    let review = include_str!("../evidence/f1.6/review.md");
    let validation = include_str!("../evidence/f1.6/validation.md");
    let plan = include_str!("../design/plan-implementacion-f1-0.1.md");
    let traceability = include_str!("../design/trazabilidad-requisitos-0.1.md");
    let resume = include_str!("../design/punto-reanudacion.md");
    let executed = include_str!("../evidence/f1.6/executed-test-list.txt");
    let execution_source = include_str!("../src/execution.rs");
    let result_source = include_str!("../src/results.rs");
    let effective_source = include_str!("effective_config.rs");
    let decimal_source = include_str!("../src/exact_decimal.rs");
    let publication_source = include_str!("../src/publication.rs");
    let publication_tests = include_str!("results_publication.rs");
    let expediente = [
        requirements,
        review,
        validation,
        plan,
        traceability,
        resume,
        executed,
    ]
    .join("\n");

    for retired in [
        "final_artifact_revalidation_closes_the_controlled_toctou_window",
        "canonical_paths_and_exact_nofollow_tree_close_publication_aliases",
    ] {
        assert!(
            !expediente.contains(retired),
            "retired test name remains in the F1.6 expediente: {retired}"
        );
    }
    for (name, source) in [
        (
            "execution::tests::semantic_manifest_gate_rejects_cross_collection_contradictions_from_reloaded_bytes",
            execution_source,
        ),
        (
            "execution::tests::global_manifest_links_must_match_catalog_exactly_without_reductions",
            execution_source,
        ),
        (
            "execution::tests::global_manifest_links_match_catalog_for_complete_failed_zero_one_and_two_reductions",
            execution_source,
        ),
        (
            "execution::tests::manifest_graph_rejects_every_dynamic_edge_contradiction",
            execution_source,
        ),
        (
            "execution::tests::manifest_dynamic_graph_rejects_diagnostic_and_reference_inconsistencies",
            execution_source,
        ),
        (
            "execution::tests::manifest_bundle_closes_exact_percent_seed_header_normalization_and_event_edges",
            execution_source,
        ),
        (
            "execution::tests::manifest_bundle_closes_toml_fasta_and_effective_value_edges",
            execution_source,
        ),
        (
            "execution::tests::resolution_evidence_toctou_is_rejected_before_publication",
            execution_source,
        ),
        (
            "effective_config::phase_ids_are_globally_unique_across_scenarios_independently_of_order",
            effective_source,
        ),
        (
            "effective_config::parameter_order_v1_is_literal_and_independent_of_rule_and_request_permutations",
            effective_source,
        ),
        (
            "publication::tests::internal_hooks_reject_extra_file_symlink_and_content_before_rename",
            publication_source,
        ),
        (
            "publication::tests::canonical_relative_path_rejects_every_textual_alias_directly",
            publication_source,
        ),
        (
            "publication::tests::intermediate_directory_symlink_is_rejected_before_publisher_and_preserves_incomplete",
            publication_source,
        ),
        (
            "results_publication::manifest_capability_is_bound_to_the_exact_workspace_catalog_and_files",
            publication_tests,
        ),
    ] {
        let function = name.rsplit("::").next().expect("test function name");
        assert!(
            source.contains(&format!("fn {function}(")),
            "expediente witness is absent from source: {name}"
        );
        assert!(
            executed.lines().any(|line| line == name),
            "expediente witness is absent from executed listing: {name}"
        );
        assert!(
            review.contains(name),
            "review does not cite exact witness {name}"
        );
    }

    let semantic_mutations = [
        "wrong-fasta-digest",
        "wrong-map-path",
        "wrong-source",
        "wrong-derived-origin",
        "orphan-reduced-reference",
        "reused-reference-fasta",
        "incompatible-origin-kind",
        "crossed-derived-fasta-id",
        "crossed-map-id",
        "crossed-derived-reference-id",
    ];
    for mutation in semantic_mutations {
        assert!(
            execution_source.contains(&format!("(\"{mutation}\",")),
            "semantic mutation absent from executable corpus: {mutation}"
        );
        assert!(
            validation.contains(&format!("`{mutation}`")),
            "validation does not enumerate executed mutation: {mutation}"
        );
    }
    for positive in ["positive-complete", "positive-failed"] {
        assert!(execution_source.contains(positive));
        assert!(validation.contains(&format!("`{positive}`")));
    }
    for global_link_case in [
        "positive-zero-complete",
        "positive-zero-failed",
        "positive-one-complete",
        "positive-one-failed",
        "positive-two-complete",
        "positive-two-failed",
        "zero-complete-reference-digest",
        "zero-complete-authoritative-resource",
        "zero-complete-effective-config-digest",
        "one-complete-reference-digest",
        "one-complete-authoritative-resource",
        "one-complete-effective-config-digest",
        "one-failed-reference-digest",
        "one-failed-authoritative-resource",
        "one-failed-effective-config-digest",
        "two-complete-reference-digest",
        "two-complete-authoritative-resource",
        "two-complete-effective-config-digest",
        "two-failed-reference-digest",
        "two-failed-authoritative-resource",
        "two-failed-effective-config-digest",
        "zero-failed-effective-config-digest",
        "zero-complete-event-digest",
        "zero-reference-id",
        "zero-reference-kind",
        "zero-reference-path",
        "zero-reference-size",
        "zero-effective-config-id",
        "zero-effective-config-kind",
        "zero-effective-config-path",
        "zero-effective-config-size",
        "zero-event-manifest-id",
        "zero-event-manifest-kind",
        "zero-event-manifest-path",
        "zero-event-manifest-size",
    ] {
        assert!(
            execution_source.contains(global_link_case),
            "global-link case absent from executable corpus: {global_link_case}"
        );
        assert!(
            validation.contains(&format!("`{global_link_case}`")),
            "validation does not enumerate global-link case: {global_link_case}"
        );
    }
    for graph_case in [
        "graph-positive-zero-complete",
        "graph-positive-zero-failed",
        "graph-positive-one-complete",
        "graph-positive-one-failed",
        "graph-positive-two-complete",
        "graph-positive-two-failed",
        "graph-positive-multiple-complete",
        "graph-positive-multiple-failed",
        "graph-configuration-copy-digest",
        "graph-configuration-copy-id",
        "graph-configuration-copy-class",
        "graph-configuration-copy-path",
        "graph-configuration-copy-size",
        "graph-configuration-copy-reused",
        "graph-configuration-input-duplicate",
        "graph-configuration-kind-incompatible",
        "graph-general-configuration-missing",
        "graph-resolution-evidence-digest",
        "graph-resolution-evidence-id",
        "graph-resolution-evidence-kind",
        "graph-resolution-evidence-path",
        "graph-resolution-evidence-size",
        "graph-event-unit-crossed",
        "graph-event-unit-orphan",
        "graph-event-unit-missing",
        "graph-execution-unit-crossed",
        "graph-output-artifact-missing",
        "graph-output-wrong-category",
        "graph-artifact-path-duplicate",
        "graph-orphan-event-artifact",
        "graph-failed-configuration-copy-digest",
        "graph-failed-execution-unit-crossed",
        "graph-one-complete-output-missing",
        "graph-one-failed-output-missing",
        "graph-two-complete-output-missing",
        "graph-two-failed-output-missing",
        "graph-two-complete-artifact-path-duplicate",
    ] {
        assert!(
            execution_source.contains(graph_case),
            "manifest graph case absent from executable corpus: {graph_case}"
        );
        assert!(
            validation.contains(&format!("`{graph_case}`")),
            "validation does not enumerate manifest graph case: {graph_case}"
        );
    }
    for dynamic_case in [
        "dynamic-positive-zero-complete",
        "dynamic-positive-zero-failed",
        "dynamic-positive-one-complete",
        "dynamic-positive-one-failed",
        "dynamic-positive-two-complete",
        "dynamic-positive-two-failed",
        "dynamic-positive-percent-complete",
        "dynamic-diagnostic-execution-absent",
        "dynamic-diagnostic-reduction-absent",
        "dynamic-diagnostic-invented",
        "dynamic-diagnostic-category",
        "dynamic-diagnostic-meaning",
        "dynamic-diagnostic-required-context",
        "dynamic-diagnostic-duplicate",
        "dynamic-diagnostic-orphan",
        "dynamic-reference-sequence-duplicate",
        "dynamic-reference-total",
        "dynamic-reference-unknown-count",
        "dynamic-reduction-partition-overlap",
        "dynamic-reduction-partition-missing",
        "dynamic-reduction-selection-order",
        "dynamic-reduction-sequence-missing",
        "dynamic-reduction-out-of-range",
        "dynamic-reduction-overlap",
        "dynamic-reduction-global-sum",
        "dynamic-reduction-per-sequence-sum",
        "dynamic-reduction-eligible",
        "dynamic-reduction-declared-target",
        "dynamic-reduction-derived-target",
        "dynamic-reduction-region",
        "dynamic-reduction-residual",
        "dynamic-reduction-derived-class",
        "dynamic-reduction-derived-ploidy",
        "dynamic-reduction-percent-target",
    ] {
        assert!(
            execution_source.contains(dynamic_case),
            "dynamic invariant case absent from executable corpus: {dynamic_case}"
        );
        assert!(
            validation.contains(&format!("`{dynamic_case}`")),
            "validation does not enumerate dynamic invariant case: {dynamic_case}"
        );
    }
    for exhaustive_case in [
        "exhaustive-positive-percent-above-integer-boundary",
        "exhaustive-positive-percent-exponent-notation",
        "exhaustive-positive-percent-below-integer-boundary",
        "exhaustive-positive-non-utf8-header",
        "exhaustive-positive-multiple-events-normalizations",
        "exhaustive-positive-two-reductions",
        "exhaustive-reduction-root-derived",
        "exhaustive-reduction-subseed-derived",
        "exhaustive-derived-reference-order",
        "exhaustive-header-token",
        "exhaustive-header-description",
        "exhaustive-header-description-non-utf8",
        "exhaustive-normalization-sum-before",
        "exhaustive-normalization-delta-before",
        "exhaustive-normalization-before",
        "exhaustive-normalization-after",
        "exhaustive-normalization-sum-after",
        "exhaustive-normalization-delta-after",
        "exhaustive-normalization-decision",
        "exhaustive-normalization-order",
        "exhaustive-normalization-duplicate-value",
        "exhaustive-normalization-duplicate-parameter",
        "exhaustive-event-general-subseed",
        "exhaustive-event-effective-subseed",
        "exhaustive-event-manifest-subseed",
        "exhaustive-event-correlated-derived-subseed",
    ] {
        assert!(
            execution_source.contains(exhaustive_case),
            "exhaustive invariant case absent from executable corpus: {exhaustive_case}"
        );
        assert!(
            validation.contains(&format!("`{exhaustive_case}`")),
            "validation does not enumerate exhaustive case: {exhaustive_case}"
        );
    }
    for r17_case in [
        "r17-positive-percent-separators",
        "r17-positive-percent-plus-exponent",
        "r17-positive-percent-high-precision",
        "r17-positive-fasta-non-utf8",
        "r17-positive-effective-multiple-scopes",
        "r17-positive-two-reductions",
        "r17-percent-leading-separator",
        "r17-percent-trailing-separator",
        "r17-percent-double-separator",
        "r17-percent-separator-before-dot",
        "r17-percent-separator-after-dot",
        "r17-percent-missing-integer",
        "r17-percent-missing-fraction",
        "r17-percent-leading-zero-integer",
        "r17-percent-leading-zero-float",
        "r17-percent-signed-leading-zero",
        "r17-percent-empty-exponent",
        "r17-percent-empty-signed-exponent",
        "r17-percent-leading-exponent-separator",
        "r17-percent-trailing-exponent-separator",
        "r17-fasta-correlated-id-header",
        "r17-fasta-correlated-description",
        "r17-fasta-correlated-length",
        "r17-fasta-correlated-unknown-count",
        "r17-fasta-correlated-order",
        "r17-fasta-intermediate-symlink",
        "r17-effective-duplicate-identical",
        "r17-effective-duplicate-contradictory",
        "r17-effective-reordered",
        "r17-effective-parameter-cross",
        "r17-effective-value-cross",
        "r17-effective-winner-cross",
        "r17-effective-selectors-cross",
        "r17-effective-provenance-path",
        "r17-effective-provenance-location",
        "r17-effective-automatic-missing",
        "r17-effective-correlated-ploidy",
        "r17-effective-source-path",
        "r17-effective-controlled-duplicate",
        "r17-effective-controlled-count",
    ] {
        assert!(
            execution_source.contains(r17_case),
            "R17 case absent from executable corpus: {r17_case}"
        );
        assert!(
            validation.contains(&format!("`{r17_case}`")),
            "validation does not enumerate R17 case: {r17_case}"
        );
    }
    for r18_case in [
        "r18-effective-dominated-base-copy",
        "r18-effective-foreign-object",
    ] {
        assert!(
            execution_source.contains(r18_case),
            "R18 case absent from executable corpus: {r18_case}"
        );
        assert!(
            validation.contains(&format!("`{r18_case}`")),
            "validation does not enumerate R18 case: {r18_case}"
        );
    }
    for r19_r21_case in [
        "r19-biological-parameter-technical-context",
        "r21-orphan-position-with-individual",
    ] {
        assert!(
            execution_source.contains(r19_r21_case),
            "R19/R21 case absent from executable corpus: {r19_r21_case}"
        );
        assert!(
            validation.contains(&format!("`{r19_r21_case}`")),
            "validation does not enumerate R19/R21 case: {r19_r21_case}"
        );
    }
    for r22_r24_case in [
        "r22-coverage-locus",
        "r22-recombination-region",
        "r22-ploidy-locus",
        "r22-reference-object",
        "r23-complete-missing-phase-domain-proof",
        "r23-failed-missing-phase-domain-proof",
        "r24-complete-scenario-copy-id",
        "r24-failed-scenario-copy-id",
    ] {
        assert!(
            execution_source.contains(r22_r24_case),
            "R22-R24 case absent from executable corpus: {r22_r24_case}"
        );
        assert!(
            validation.contains(&format!("`{r22_r24_case}`")),
            "validation does not enumerate R22-R24 case: {r22_r24_case}"
        );
    }
    for r25_r27_case in [
        "r25-complete-correlated-domain-provider-interval",
        "r25-failed-correlated-domain-provider-interval",
        "r25-correlated-global-output-to-event",
        "r26-complete-duplicate-phase-across-scenarios",
        "r26-failed-duplicate-phase-across-scenarios",
        "r27-complete-parameter-order",
        "r27-failed-parameter-order",
    ] {
        assert!(
            execution_source.contains(r25_r27_case),
            "R25-R27 case absent from executable corpus: {r25_r27_case}"
        );
        assert!(
            validation.contains(&format!("`{r25_r27_case}`")),
            "validation does not enumerate R25-R27 case: {r25_r27_case}"
        );
    }
    let oracle_start = execution_source
        .find("type IndependentDynamicInvariant")
        .expect("independent oracle definition");
    let oracle_end = execution_source[oracle_start..]
        .find("fn complete_and_failed_manifest_corpus")
        .map(|offset| oracle_start + offset)
        .expect("independent oracle boundary");
    let oracle = &execution_source[oracle_start..oracle_end];
    for forbidden_dependency in [
        "validate_general_manifest_bytes",
        "python_schema_accepts",
        "validate_general_manifest_reduction_links",
        "validate_general_manifest_graph",
    ] {
        assert!(
            !oracle.contains(forbidden_dependency),
            "independent semantic oracle calls production path {forbidden_dependency}"
        );
    }
    assert!(
        execution_source.contains("independent_semantic_manifest_oracle(&path, &bytes).is_err()")
    );
    let exact_start = execution_source
        .find("fn independent_retained_after_percent")
        .expect("independent exact percentage oracle");
    let exact_end = execution_source[exact_start..]
        .find("fn independent_decimal_at_least_hundred")
        .map(|offset| exact_start + offset)
        .expect("independent exact percentage boundary");
    let exact_oracle = &execution_source[exact_start..exact_end];
    assert!(!exact_oracle.contains("parse::<f64>"));
    assert!(!exact_oracle.contains("eligible as f64"));
    assert!(execution_source.contains("30.0000000000000000001"));
    assert!(decimal_source.contains("fn parse_toml_decimal("));
    assert!(execution_source.contains("fn independent_toml_decimal_parts("));
    assert!(result_source.contains("fn read_package_file_nofollow("));
    assert!(result_source.contains("OFlags::NOFOLLOW"));
    assert!(result_source.contains("fn validate_authoritative_fastas("));
    assert!(result_source.contains("fn validate_effective_values("));
    assert!(result_source.contains("pub fn resolution_evidence_document("));
    assert!(result_source.contains("fn parse_resolution_evidence("));
    assert!(result_source.contains("fn validate_reloaded_artifact_identity("));
    assert!(result_source.contains("validate_resolution_request("));
    assert!(execution_source.contains("fn independent_authoritative_fasta_graph("));
    assert!(execution_source.contains("fn independent_effective_exclusive_graph("));
    assert!(execution_source.contains("fn independent_rule_matches("));
    assert!(execution_source.contains("fn independent_selector_dominates("));
    assert!(execution_source.contains("fn independent_allowed_resolution_axes("));
    assert!(execution_source.contains("fn independent_allowed_spatial_levels("));
    assert!(execution_source.contains("fn independent_add_scenario_resolution_ids("));
    assert!(execution_source.contains("fn independent_resolution_domain_shape("));
    assert!(execution_source.contains("fn independent_parameter_rank("));
    assert!(execution_source.contains("independent forbidden resolution axis"));
    assert!(execution_source.contains("independent subpopulation without population"));
    assert!(execution_source.contains("independent generation without phase"));
    assert!(execution_source.contains("independent orphan spatial coordinate"));
    assert!(execution_source.contains("independent position outside sequence"));
    assert!(execution_source.contains("independent forbidden resolution spatial level"));
    assert!(execution_source.contains("independent unproven population hierarchy"));
    assert!(execution_source.contains("independent unproven temporal hierarchy"));
    assert!(execution_source.contains("independent unproven locus hierarchy"));
    assert!(execution_source.contains("independent scenario copy class/id mismatch"));
    assert!(validation.contains("diez mutaciones semánticas"));
    for independent_check in [
        "authoritative resource contradicts its FASTA",
        "effective configuration contradicts catalog",
        "event manifest absent from catalog",
        "configuration copy contradicts catalog",
        "event units do not equal event links",
        "output categories do not partition artifacts",
        "artifact is unreachable from graph roots",
        "independent diagnostic contract mismatch",
        "independent orphan diagnostic or missing unit edge",
        "independent reference total mismatch",
        "independent source partition mismatch",
        "independent overlapping selection",
        "independent derived sequence mismatch",
        "independent normalization decision",
        "independent reduction subseed mismatch",
        "independent linked event identity/subseed/unit",
    ] {
        assert!(
            oracle.contains(independent_check),
            "independent oracle misses global check: {independent_check}"
        );
    }
}
