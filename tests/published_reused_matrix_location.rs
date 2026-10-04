//! DT-022: el juez de los recursos publicados no depende del directorio de
//! trabajo cuando el paquete reutiliza la matriz de SNV.
//!
//! La matriz reutilizada no viaja en el paquete. Su ruta, escrita en el TOML
//! copiado, es relativa al directorio de trabajo de la ejecución original; el
//! juez recibe ahora la ubicación de la matriz de quien lo llama.
//!
//! **Lo que discrimina la prueba principal es mover la matriz a un sitio que
//! ninguna ruta del TOML alcanza, desde ningún directorio**; no cambiar de
//! directorio de trabajo. Una versión anterior lanzaba el juez en un subproceso
//! con otro directorio, y esa parte no aportaba nada: la mutación que quitaba
//! el cambio de directorio moría igual, y un hijo que aceptase siempre dejaba
//! pasar la prueba en vacío.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::results::validate_published_resources;

const BINARY: &str = env!("CARGO_BIN_EXE_genome-population-simulator");
static NONCE: AtomicU64 = AtomicU64::new(0);

fn temporary(label: &str) -> PathBuf {
    let root = std::env::temp_dir().join(format!(
        "gps-dt022-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("temporary root");
    root
}

fn reference_block() -> &'static str {
    "[reference]\nkind = \"simulated\"\nid = \"ref-main\"\n\
     [reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n\
     [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 60\n"
}

fn generate_config(locus_count: u32) -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"2\"\n\
         seed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n\
         [population_snvs]\nmode = \"generate\"\n[population_snvs.generate]\n\
         count_mode = \"exact\"\nlocus_count = {locus_count}\nspatial_mode = \"uniform\"\n\
         substitution_profile = \"uniform_snv12_v1\"\n\
         [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
         frequency_model = \"uniform\"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n{}",
        reference_block()
    )
}

fn reuse_config(metadata_path: &str) -> String {
    format!(
        "config_kind = \"general\"\nschema_version = \"2\"\n\
         seed = \"000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f\"\n\
         [population_snvs]\nmode = \"reuse\"\nmetadata_path = \"{metadata_path}\"\n{}",
        reference_block()
    )
}

/// Ejecuta la CLI desde `cwd` con rutas tal como se las damos.
fn run_cli(cwd: &Path, config: &str, output: &str) {
    let result = Command::new(BINARY)
        .current_dir(cwd)
        .args(["run", config, "--output", output])
        .output()
        .expect("CLI");
    assert!(
        result.status.success(),
        "la ejecución debe publicar: {}{}",
        String::from_utf8_lossy(&result.stdout),
        String::from_utf8_lossy(&result.stderr)
    );
}

fn judge(package: &Path, metadata: Option<&Path>) -> Result<(), String> {
    let effective_path = package.join("effective-config.json");
    let manifest_path = package.join("manifest.json");
    validate_published_resources(
        package,
        &effective_path,
        &fs::read(&effective_path).expect("effective bytes"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifest bytes"),
        metadata,
    )
    .map_err(|diagnostic| format!("{diagnostic:?}").replace('\n', " "))
}

/// Un paquete con la matriz reutilizada por una ruta **relativa** en el TOML.
struct RelativeReuse {
    root: PathBuf,
    package: PathBuf,
}

fn relative_reuse_package(label: &str) -> RelativeReuse {
    let root = temporary(label);
    fs::write(root.join("generate.toml"), generate_config(2)).expect("generate TOML");
    fs::write(
        root.join("reuse.toml"),
        reuse_config("matrix/population_snvs.metadata.json"),
    )
    .expect("reuse TOML");
    run_cli(&root, "generate.toml", "matrix");
    run_cli(&root, "reuse.toml", "reused");
    let package = root.join("reused");
    RelativeReuse { root, package }
}

#[test]
fn a_reused_matrix_is_judged_from_the_given_location_not_from_the_toml_path() {
    let case = relative_reuse_package("cwd");
    // Prueba de que el paquete se publicó con ruta relativa: es lo que el juez
    // antiguo resolvía contra el directorio de trabajo.
    let manifest = fs::read_to_string(case.package.join("manifest.json")).expect("manifest");
    assert!(manifest.contains("\"original_path\": \"reuse.toml\""));
    // La matriz se mueve a un sitio que ninguna ruta relativa del TOML alcanza
    // desde ningún directorio de trabajo.
    let moved = case.root.join("elsewhere/deep/moved-matrix");
    fs::create_dir_all(moved.parent().expect("parent")).expect("moved parent");
    fs::rename(case.root.join("matrix"), &moved).expect("move matrix");
    let metadata = moved.join("population_snvs.metadata.json");

    judge(&case.package, Some(&metadata))
        .expect("con la ubicación indicada se acepta, aunque el TOML apunte a otro sitio");
    // Control negativo: sin la ubicación no se acepta, de modo que lo que
    // acepta arriba es el dato y no una ruta que el juez encontrara solo.
    let detail = judge(&case.package, None).expect_err("sin ubicación se rechaza");
    assert!(detail.contains("falta indicar la ruta"), "{detail}");
}

#[test]
fn a_package_that_reuses_the_matrix_requires_its_location() {
    let root = temporary("missing");
    fs::write(root.join("generate.toml"), generate_config(2)).expect("generate TOML");
    // Ruta **absoluta** en el TOML: el juez antiguo habría encontrado la matriz
    // sin que nadie se la dijera, y esta prueba exige que aun así falte el dato.
    let metadata = root.join("matrix/population_snvs.metadata.json");
    fs::write(
        root.join("reuse.toml"),
        reuse_config(&metadata.display().to_string()),
    )
    .expect("reuse TOML");
    run_cli(&root, "generate.toml", "matrix");
    run_cli(&root, "reuse.toml", "reused");
    let package = root.join("reused");

    judge(&package, Some(&metadata)).expect("con la ubicación, se acepta");
    let detail = judge(&package, None).expect_err("sin la ubicación, se rechaza");
    assert!(
        detail.contains("falta indicar la ruta de population_snvs.metadata.json"),
        "el error debe decir qué falta: {detail}"
    );
    // Falta un argumento del llamador: no es un defecto de esquema del
    // documento publicado.
    assert!(detail.contains("GPS-INPUT-001"), "{detail}");
    assert!(!detail.contains("GPS-FMT-006"), "{detail}");
    // Un paquete que no reutiliza no necesita el dato, y si se le da, se ignora.
    let generated = root.join("matrix");
    judge(&generated, None).expect("sin reutilización no hace falta");
    judge(&generated, Some(&metadata)).expect("y si se da, no se usa");
}

#[test]
fn another_matrix_with_the_same_file_names_is_rejected_by_fingerprint() {
    let root = temporary("fingerprint");
    fs::write(root.join("generate.toml"), generate_config(2)).expect("generate TOML");
    fs::write(root.join("generate3.toml"), generate_config(3)).expect("generate TOML");
    let metadata = root.join("matrix/population_snvs.metadata.json");
    fs::write(
        root.join("reuse.toml"),
        reuse_config(&metadata.display().to_string()),
    )
    .expect("reuse TOML");
    run_cli(&root, "generate.toml", "matrix");
    run_cli(&root, "generate3.toml", "other");
    run_cli(&root, "reuse.toml", "reused");
    let package = root.join("reused");
    judge(&package, Some(&metadata)).expect("la matriz correcta se acepta");

    // Otra matriz válida para la misma referencia, con los mismos nombres de
    // fichero: la huella de los metadatos que el manifiesto registra no coincide.
    let other = root.join("other/population_snvs.metadata.json");
    let detail = judge(&package, Some(&other)).expect_err("otros metadatos se rechazan");
    assert!(
        detail.contains("enlace SNV reutilizado no coincide con el recurso físico"),
        "debe rechazarse por la huella de los metadatos: {detail}"
    );

    // La misma matriz con los metadatos reescritos con otra sangría: mismo
    // contenido y misma auditoría, otros bytes. Solo la huella de los
    // metadatos que registra el manifiesto puede distinguirlos; la auditoría
    // enlazada es idéntica, de modo que esta es la aserción que cubre
    // `metadata_sha256` por sí sola.
    let reformatted = root.join("reformatted");
    fs::create_dir_all(&reformatted).expect("reformatted");
    let original: serde_json::Value =
        serde_json::from_slice(&fs::read(&metadata).expect("metadata")).expect("metadata JSON");
    let reformatted_bytes = serde_json::to_vec(&original).expect("compact metadata");
    assert_ne!(
        reformatted_bytes,
        fs::read(&metadata).expect("metadata"),
        "el fixture debe cambiar los bytes"
    );
    fs::write(
        reformatted.join("population_snvs.metadata.json"),
        &reformatted_bytes,
    )
    .expect("reformatted metadata");
    fs::copy(
        root.join("matrix/population_snvs.tsv"),
        reformatted.join("population_snvs.tsv"),
    )
    .expect("tsv");
    let detail = judge(
        &package,
        Some(&reformatted.join("population_snvs.metadata.json")),
    )
    .expect_err("metadatos con otros bytes se rechazan");
    assert!(
        detail.contains("enlace SNV reutilizado no coincide con el recurso físico"),
        "debe rechazarse por la huella de los metadatos: {detail}"
    );

    // Los metadatos correctos con el TSV de la otra matriz: se rechaza. **Esta
    // aserción no discrimina la comprobación de `matrix.sha256` por sí sola**:
    // las comprobaciones posteriores de bytes canónicos y de auditoría lo
    // rechazarían igualmente (mutar solo la comparación de la huella del TSV
    // no hace fallar nada). Lo que prueba es que el juez no acepta un par
    // híbrido, sea cual sea la guarda que lo rechace.
    let hybrid = root.join("hybrid");
    fs::create_dir_all(&hybrid).expect("hybrid");
    fs::copy(&metadata, hybrid.join("population_snvs.metadata.json")).expect("metadata");
    fs::copy(
        root.join("other/population_snvs.tsv"),
        hybrid.join("population_snvs.tsv"),
    )
    .expect("tsv");
    let detail = judge(
        &package,
        Some(&hybrid.join("population_snvs.metadata.json")),
    )
    .expect_err("el TSV de otra matriz se rechaza");
    assert!(
        detail.contains("reuso SNV físico inválido"),
        "debe rechazarse como reutilización física inválida: {detail}"
    );
}

#[test]
fn a_missing_matrix_next_to_the_metadata_names_the_file() {
    let root = temporary("no-tsv");
    fs::write(root.join("generate.toml"), generate_config(2)).expect("generate TOML");
    let metadata = root.join("matrix/population_snvs.metadata.json");
    fs::write(
        root.join("reuse.toml"),
        reuse_config(&metadata.display().to_string()),
    )
    .expect("reuse TOML");
    run_cli(&root, "generate.toml", "matrix");
    run_cli(&root, "reuse.toml", "reused");
    let package = root.join("reused");
    let lonely = root.join("lonely");
    fs::create_dir_all(&lonely).expect("lonely");
    fs::copy(&metadata, lonely.join("population_snvs.metadata.json")).expect("metadata");

    let detail = judge(
        &package,
        Some(&lonely.join("population_snvs.metadata.json")),
    )
    .expect_err("sin el TSV junto a los metadatos se rechaza");
    let expected = lonely.join("population_snvs.tsv").display().to_string();
    assert!(
        detail.contains(&expected),
        "el error debe nombrar el fichero que falta, {expected}: {detail}"
    );
}

#[test]
fn a_forged_normalization_audit_is_rejected_even_with_the_right_matrix() {
    let root = temporary("audit");
    fs::write(root.join("generate.toml"), generate_config(2)).expect("generate TOML");
    let metadata = root.join("matrix/population_snvs.metadata.json");
    fs::write(
        root.join("reuse.toml"),
        reuse_config(&metadata.display().to_string()),
    )
    .expect("reuse TOML");
    run_cli(&root, "generate.toml", "matrix");
    run_cli(&root, "reuse.toml", "reused");
    let package = root.join("reused");

    // La auditoría se falsifica a la vez en la configuración efectiva y en el
    // manifiesto, de modo que entre ellos concuerdan y solo la comparación con
    // la auditoría recalculada desde la matriz física puede rechazarla.
    let forge = |name: &str| -> Vec<u8> {
        let mut document: serde_json::Value =
            serde_json::from_slice(&fs::read(package.join(name)).expect("document")).expect("JSON");
        assert_ne!(
            document["population_snvs"]["normalization_audit"][0]["original_sum"],
            serde_json::json!(0.0),
            "el fixture debe cambiar el valor"
        );
        document["population_snvs"]["normalization_audit"][0]["original_sum"] =
            serde_json::json!(0.0);
        serde_json::to_vec(&document).expect("forged bytes")
    };
    let error = validate_published_resources(
        &package,
        &package.join("effective-config.json"),
        &forge("effective-config.json"),
        &package.join("manifest.json"),
        &forge("manifest.json"),
        Some(&metadata),
    )
    .expect_err("una auditoría falsificada se rechaza");
    let detail = format!("{error:?}");
    assert!(
        detail.contains("enlace SNV reutilizado no coincide con el recurso físico"),
        "debe rechazarse por la auditoría enlazada: {detail}"
    );
}
