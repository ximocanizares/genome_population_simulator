//! Regresión de retirada: los documentos F1 archivados no se admiten como
//! contratos publicados tras la transición atómica DTC-082.

use std::fs;
use std::path::{Path, PathBuf};

use genome_population_simulator::results::{
    validate_effective_config_against_manifest_bytes, validate_general_manifest_bytes,
    validate_published_resources,
};

const F1_EFFECTIVE: &[u8] = include_bytes!("fixtures/f1-e14705a/effective-config.json");
const F1_MANIFEST: &[u8] = include_bytes!("fixtures/f1-e14705a/manifest.json");
const F1_REDUCTION_EFFECTIVE: &[u8] =
    include_bytes!("fixtures/f1-e14705a-reduction/effective-config.json");
const F1_REDUCTION_MANIFEST: &[u8] = include_bytes!("fixtures/f1-e14705a-reduction/manifest.json");

#[test]
fn archived_e14705a_f1_documents_are_rejected_without_a_legacy_reader() {
    let effective_path = Path::new("effective-config.json");
    let manifest_path = Path::new("manifest.json");

    let manifest_error = validate_general_manifest_bytes(manifest_path, F1_MANIFEST)
        .expect_err("el manifiesto archivado v1 debe rechazarse");
    assert_eq!(
        manifest_error.context().get("schema").map(String::as_str),
        Some("manifest-general-v4.schema.json")
    );
    let effective_error = validate_effective_config_against_manifest_bytes(
        effective_path,
        F1_EFFECTIVE,
        Path::new("archived/effective-config.json"),
        F1_MANIFEST,
    )
    .expect_err("la configuración efectiva archivada v1 debe rechazarse");
    assert_eq!(
        effective_error.context().get("schema").map(String::as_str),
        Some("effective-config-v4.schema.json")
    );

    let effective: serde_json::Value = serde_json::from_slice(F1_EFFECTIVE).expect("golden JSON");
    let manifest: serde_json::Value = serde_json::from_slice(F1_MANIFEST).expect("golden JSON");
    assert_eq!(effective["schema_version"], "1");
    assert_eq!(manifest["schema_version"], "1");
}

#[test]
fn version_dispatch_is_closed_and_reports_the_schema_selected_by_the_document() {
    let mut broken_v1_manifest: serde_json::Value =
        serde_json::from_slice(F1_MANIFEST).expect("golden manifest");
    broken_v1_manifest["unexpected"] = serde_json::json!(true);
    let error = validate_general_manifest_bytes(
        Path::new("manifest.json"),
        &serde_json::to_vec(&broken_v1_manifest).expect("JSON"),
    )
    .expect_err("v1 unknown field");
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some("manifest-general-v4.schema.json")
    );

    let mut unsupported_effective: serde_json::Value =
        serde_json::from_slice(F1_EFFECTIVE).expect("golden effective configuration");
    unsupported_effective["schema_version"] = serde_json::json!("2");
    let error = validate_effective_config_against_manifest_bytes(
        Path::new("effective-config.json"),
        &serde_json::to_vec(&unsupported_effective).expect("JSON"),
        Path::new("manifest.json"),
        F1_MANIFEST,
    )
    .expect_err("unsupported version");
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some("effective-config-v4.schema.json")
    );
    assert!(
        error
            .context()
            .get("detail")
            .is_some_and(|detail| detail.contains("se admiten \"4\" y \"6\""))
    );

    // El despacho sigue cerrado después de abrir la versión 6 para los
    // documentos generales v4: solo esas dos versiones se admiten, y cualquier
    // otra se rechaza nombrando el esquema.
    for rejected in ["1", "3", "5", "7"] {
        let mut document: serde_json::Value =
            serde_json::from_slice(F1_EFFECTIVE).expect("golden effective configuration");
        document["schema_version"] = serde_json::json!(rejected);
        let error = validate_effective_config_against_manifest_bytes(
            Path::new("effective-config.json"),
            &serde_json::to_vec(&document).expect("JSON"),
            Path::new("manifest.json"),
            F1_MANIFEST,
        )
        .expect_err("versión no admitida");
        assert_eq!(
            error.context().get("schema").map(String::as_str),
            Some("effective-config-v4.schema.json"),
            "versión rechazada: {rejected}"
        );
    }
}

#[test]
fn old_pairs_and_physical_validation_report_the_current_schema_on_rejection() {
    let mut broken: serde_json::Value = serde_json::from_slice(F1_MANIFEST).expect("manifest");
    broken["unexpected"] = serde_json::json!(true);
    let manifest_path = Path::new("archived/manifest.json");
    assert_rejected_manifest_error(
        validate_effective_config_against_manifest_bytes(
            Path::new("archived/effective-config.json"),
            F1_EFFECTIVE,
            manifest_path,
            &serde_json::to_vec(&broken).expect("JSON"),
        )
        .expect_err("the v1 pair must reject a corrupt manifest"),
        Path::new("archived/effective-config.json"),
    );

    let root = temporary("corrupt-v1-physical");
    fs::create_dir_all(&root).expect("temporary tree");
    write_fixture_tree(&root);
    let effective_path = root.join("effective-config.json");
    let physical_manifest_path = root.join("manifest.json");
    let mut broken: serde_json::Value =
        serde_json::from_slice(F1_REDUCTION_MANIFEST).expect("reduction manifest");
    broken["unexpected"] = serde_json::json!(true);
    assert_rejected_manifest_error(
        validate_published_resources(
            &root,
            &effective_path,
            F1_REDUCTION_EFFECTIVE,
            &physical_manifest_path,
            &serde_json::to_vec(&broken).expect("JSON"),
            None,
        )
        .expect_err("the v1 physical validator must reject a corrupt manifest"),
        &effective_path,
    );
    fs::remove_dir_all(root).expect("cleanup");
}

fn assert_rejected_manifest_error(
    error: genome_population_simulator::diagnostics::Diagnostic,
    path: &Path,
) {
    let expected_path = path.display().to_string();
    assert_eq!(
        error.context().get("schema").map(String::as_str),
        Some("effective-config-v4.schema.json")
    );
    assert_eq!(
        error.context().get("path").map(String::as_str),
        Some(expected_path.as_str())
    );
}

fn write_fixture_tree(root: &Path) {
    for (relative, bytes) in [
        ("effective-config.json", F1_REDUCTION_EFFECTIVE),
        ("manifest.json", F1_REDUCTION_MANIFEST),
        (
            "resources/references/reference-small.map.json",
            include_bytes!(
                "fixtures/f1-e14705a-reduction/resources/references/reference-small.map.json"
            ) as &[u8],
        ),
        (
            "resources/references/reference-small.fa",
            include_bytes!("fixtures/f1-e14705a-reduction/resources/references/reference-small.fa")
                as &[u8],
        ),
        (
            "resources/references/reference-main.fa",
            include_bytes!("fixtures/f1-e14705a-reduction/resources/references/reference-main.fa")
                as &[u8],
        ),
        (
            "inputs/config/general.toml",
            include_bytes!("fixtures/f1-e14705a-reduction/inputs/config/general.toml") as &[u8],
        ),
    ] {
        let destination = root.join(relative);
        fs::create_dir_all(destination.parent().expect("fixture parent")).expect("fixture parent");
        fs::write(destination, bytes).expect("fixture file");
    }
}

fn temporary(name: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f1-compatibility-{}-{name}",
        std::process::id()
    ))
}
