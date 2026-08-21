use genome_population_simulator::population_snvs::{
    ALT_FREQUENCY_CANONICAL_ORDER, ALT_FREQUENCY_SUM_ALGORITHM, AltFrequencyDecision,
    AltFrequencyNormalizationError, PopulationSnvsFormatError, TSV_HEADER,
    normalize_alt_frequencies, parse_canonical_tsv, parse_metadata_json,
};

fn matrix(rows: &[&str]) -> genome_population_simulator::population_snvs::PopulationSnvsTsv {
    let mut input = TSV_HEADER.join("\t");
    input.push('\n');
    input.push_str(&rows.join("\n"));
    input.push('\n');
    parse_canonical_tsv(input.as_bytes()).expect("canonical matrix shape")
}

#[test]
fn canonical_tsv_preserves_the_frozen_eight_column_shape() {
    let input = format!(
        "{}\nvar-1\tchr1\t1\tA\tC\tPOP\t1e-3\t\n",
        TSV_HEADER.join("\t")
    );
    let table = parse_canonical_tsv(input.as_bytes()).expect("canonical tsv");
    assert_eq!(table.rows().len(), 1);
    assert_eq!(table.rows()[0].fields()[6], "1e-3");
    assert_eq!(table.rows()[0].fields()[7], "");
}

#[test]
fn tsv_rejects_vcf_and_noncanonical_columns() {
    assert_eq!(
        parse_canonical_tsv(b"##fileformat=VCFv4.5\n").expect_err("VCF must fail"),
        PopulationSnvsFormatError::VcfInput
    );
    let error = parse_canonical_tsv(b"contig\tposition\tref\talt\tpopulation_id\talt_af\n")
        .expect_err("six-column input is not the canonical output");
    assert!(matches!(
        error,
        PopulationSnvsFormatError::IncorrectHeader { .. }
    ));
}

#[test]
fn metadata_requires_the_closed_v1_envelope() {
    let metadata = br#"{
        "document_kind":"population_snvs_metadata",
        "schema_version":"1",
        "coordinate_system":"base1_closed",
        "reference":{"reference_id":"ref-main","sha256":"abc"},
        "matrix":{"path":"population_snvs.tsv","sha256":"def"},
        "populations":[{"population_id":"POP","sample_size":"unknown"}],
        "provenance":"supplied"
    }"#;
    parse_metadata_json(metadata).expect("valid metadata envelope");

    let unexpected = br#"{
        "document_kind":"population_snvs_metadata", "schema_version":"1",
        "coordinate_system":"base1_closed",
        "reference":{"reference_id":"ref-main","sha256":"abc"},
        "matrix":{"path":"population_snvs.tsv","sha256":"def"},
        "populations":[], "provenance":"generated", "rows":[]
    }"#;
    assert!(matches!(
        parse_metadata_json(unexpected),
        Err(PopulationSnvsFormatError::UnexpectedMetadataField { .. })
    ));
}

#[test]
fn multiallelic_alt_sum_uses_canonical_neumaier_audit_and_normalizes_at_the_inclusive_boundary() {
    let rows = [
        "var-t\tchr1\t7\tA\tT\tPOP\t0.5000005\tz",
        "var-c\tchr1\t7\tA\tC\tPOP\t0.5\ta",
        "var-g\tchr1\t7\tA\tG\tPOP\t0\tb",
    ];
    let table = matrix(&rows);
    let original_sum = 1.0000005_f64;
    let tolerance = original_sum - 1.0;
    let normalized = normalize_alt_frequencies(&table, tolerance).expect("inclusive boundary");
    let group = &normalized.groups()[0];

    // Oráculo independiente: identidad algebraica de los valores racionales de entrada.
    assert_eq!(group.decision(), AltFrequencyDecision::Normalized);
    assert_eq!(group.original_delta().to_bits(), tolerance.to_bits());
    assert_eq!(group.original_sum().to_bits(), original_sum.to_bits());
    assert_eq!(group.values()[0].alt(), "C");
    assert_eq!(group.values()[1].alt(), "G");
    assert_eq!(group.values()[2].alt(), "T");
    assert_eq!(
        group.values()[0].effective().to_bits(),
        (0.5 / original_sum).to_bits()
    );
    assert_eq!(
        group.values()[2].effective().to_bits(),
        (0.5000005 / original_sum).to_bits()
    );
    assert!(group.effective_delta() <= tolerance);
    assert_eq!(group.canonical_order(), ALT_FREQUENCY_CANONICAL_ORDER);
    assert_eq!(group.sum_algorithm(), ALT_FREQUENCY_SUM_ALGORITHM);

    let reversed = matrix(&[rows[2], rows[1], rows[0]]);
    let reordered =
        normalize_alt_frequencies(&reversed, tolerance).expect("physical order irrelevant");
    assert_eq!(
        reordered.groups()[0].original_sum().to_bits(),
        group.original_sum().to_bits()
    );
    assert_eq!(
        reordered.groups()[0].effective_sum().to_bits(),
        group.effective_sum().to_bits()
    );
}

#[test]
fn normalized_alt_sum_has_an_exact_zero_reference_frequency_after_rounding() {
    let normalized = normalize_alt_frequencies(
        &matrix(&[
            "var-c\tchr1\t8\tA\tC\tPOP\t0.55282124170015812\t",
            "var-g\tchr1\t8\tA\tG\tPOP\t0.44717953754220685\t",
        ]),
        1e-6,
    )
    .expect("the small excess is normalizable");
    let group = &normalized.groups()[0];
    assert_eq!(group.decision(), AltFrequencyDecision::Normalized);
    assert_eq!(group.reference_frequency().to_bits(), 0.0_f64.to_bits());
    assert!(group.effective_sum().is_finite());
    assert!(group.effective_delta() <= group.tolerance());
}

#[test]
fn alt_sum_keeps_reference_remainder_and_rejects_excess_or_non_finite_values() {
    let zero = normalize_alt_frequencies(
        &matrix(&["v2\tchr1\t9\tA\tG\tP\t0\t", "v1\tchr1\t9\tA\tC\tP\t0\t"]),
        1e-9,
    )
    .expect("zero ALT is valid");
    let group = &zero.groups()[0];
    assert_eq!(group.decision(), AltFrequencyDecision::ReferenceRemainder);
    assert_eq!(group.reference_frequency().to_bits(), 1.0_f64.to_bits());
    assert_eq!(group.original_sum().to_bits(), 0.0_f64.to_bits());
    assert_eq!(group.effective_sum().to_bits(), 0.0_f64.to_bits());

    let excess = normalize_alt_frequencies(
        &matrix(&[
            "v1\tchr1\t9\tA\tC\tP\t0.6\t",
            "v2\tchr1\t9\tA\tG\tP\t0.400001\t",
        ]),
        1e-9,
    )
    .expect_err("excess beyond tolerance");
    assert!(matches!(
        excess,
        AltFrequencyNormalizationError::ExcessAltSum { .. }
    ));

    for literal in ["NaN", "inf", "-0.1", "1.1"] {
        let error = normalize_alt_frequencies(
            &matrix(&[&format!("v1\tchr1\t9\tA\tC\tP\t{literal}\t")]),
            1e-9,
        )
        .expect_err("invalid f64 must fail");
        assert!(matches!(
            error,
            AltFrequencyNormalizationError::InvalidFrequency { .. }
        ));
    }
}
