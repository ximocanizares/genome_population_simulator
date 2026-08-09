use std::collections::BTreeSet;

use toml::de::DeTable;

#[test]
fn selected_parser_rejects_duplicate_keys_with_a_source_span() {
    let source = include_str!("fixtures/config/invalid-duplicate-key.toml");
    let error = DeTable::parse(source).expect_err("the duplicate must be rejected");

    assert!(error.span().is_some());
    assert!(error.to_string().contains("duplicate key"));
}

#[test]
fn selected_parser_preserves_key_spans_for_localizable_unknown_fields() {
    let source = include_str!("fixtures/config/invalid-unknown-field.toml");
    let parsed = DeTable::parse(source).expect("the document is valid TOML");
    let accepted = BTreeSet::from(["config_kind", "schema_version", "route"]);
    let unknown = parsed
        .get_ref()
        .keys()
        .find(|key| !accepted.contains(key.get_ref().as_ref()))
        .expect("fixture must contain an unknown key");

    assert_eq!(unknown.get_ref().as_ref(), "mystery");
    assert_eq!(&source[unknown.span()], "mystery");
}
