//! DT-040 y DT-041: el catálogo de azar y el perfil técnico que publica cada
//! tipo de ejecución, y que sus documentos pasan sus esquemas.
//!
//! # Qué se corrigió (12/10/2026)
//!
//! - **DT-040.** Con tartamudeo clonal, la configuración efectiva publicaba
//!   `technical_event.rng_catalog = "f4-technical-v5"` mientras el bloque de la
//!   amplificación y el manifiesto de evento decían v6, y su esquema exige el
//!   v6. Ahora la raíz del evento publica el v6 (DTC-091: el catálogo más alto
//!   cuya etapa abre; DTC-092: solo la amplificada **sin** loci sigue en v5), y
//!   el enlace del evento en el manifiesto general también.
//! - **DT-041.** Los esquemas del manifiesto general v4 y v6 fijaban el perfil
//!   sin PCR, de modo que **toda** ejecución amplificada publicaba un
//!   manifiesto general que no pasaba su esquema. Ahora admiten el perfil con
//!   PCR atado a su catálogo: en el v4, sin PCR con v3 y con PCR con v5; en el
//!   v6, sin PCR con v3 o v4 y con PCR con v5 o v6. El v4 se amplió como
//!   excepción expresa a su inmutabilidad (DTC-089, decisión del responsable
//!   del 12/10/2026).
//!
//! # Los escenarios
//!
//! Son **los mismos ficheros** que la evidencia de DT-009
//! (`evidence/dt-009-multiplicador-congelacion/inputs/`), leídos por su ruta
//! relativa a la raíz del paquete, que es el directorio de trabajo con el que
//! `cargo test` lanza las pruebas de integración. Así la configuración efectiva
//! registra las mismas rutas que la ejecución por línea de órdenes desde la
//! raíz del repositorio, y sus bytes se pueden comparar con los del binario de
//! partida.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_published_resources;
use genome_population_simulator::technical_profile::{
    GPS_HISEQX_PCR_PE150_SHA256, GPS_HISEQX_PCR_PE150_V1, GPS_HISEQX_PCRFREE_PE150_SHA256,
    GPS_HISEQX_PCRFREE_PE150_V1,
};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const EVIDENCE: &str = "evidence/dt-009-multiplicador-congelacion/inputs";

// ---------------------------------------------------------------------------
// Las huellas congeladas.
//
// **Medidas el 12/10/2026 con el binario del commit e8135ae, ANTES de tocar una
// línea del código**: una copia sacada con `git archive`, compilada con su
// propio directorio de compilación y lanzada desde la raíz del repositorio con
// cada `general.toml` de la evidencia. Si se midieran después del cambio, la
// prueba congelaría lo que el código nuevo produce y no podría fallar por el
// defecto que dice cubrir.
//
// Las ejecuciones sin PCR, sin evento técnico y amplificadas sin loci **no se
// mueven ni un byte**. En la amplificada con loci cambian a propósito dos
// ficheros, y solo dos (DT-040): la configuración efectiva, por su
// `technical_event.rng_catalog`, y el manifiesto general, por el catálogo del
// enlace del evento y por la huella de la efectiva, que cita dos veces. Las
// huellas nuevas se midieron con el binario corregido; las de antes quedan en
// el comentario de cada una.
//
// **13/10/2026 (DT-034):** el texto de la fuente de la calibración del
// tartamudeo atribuye ahora las cifras a Gymrek (2016) además de Raz (2019), de
// modo que cambian la efectiva y el manifiesto de evento de toda ejecución con
// loci de repetición (`sin-pcr-loci` y `amplificado-loci`), y por arrastre su
// manifiesto general (huella y `size_bytes` de esos dos). Con el texto antiguo
// restituido, la efectiva y el manifiesto de evento dan las huellas de partida.
//
// **30/10/2026 (entrega 3 del plan de detección):** desde entonces un locus
// declarado solo existe con la sección de microsatélites, y los TOML de la
// evidencia con loci (`sin-pcr-loci`, `amplificado-loci`,
// `amplificado-sustitucion` y `amplificado-multiplicador`) la llevan, escrita
// por `tests/support/detection_section.rs`: una región que cubre los dos loci y
// suelos que dejan fuera la `G×4` que los separa. En las dos ejecuciones con
// loci que se congelan aquí cambian a propósito la copia del TOML, la verdad
// (por el identificador de ejecución y la huella del TOML, que lleva dentro), la
// efectiva (gana `repeat_detection` y el `origin` de cada locus) y el manifiesto
// general (gana el `origin`, y cambian el identificador de ejecución y la
// huella y el tamaño de los ficheros anteriores). Comprobado contra el binario
// de `dcae406`, compilado aparte, sobre los TOML de antes: los FASTQ, el
// manifiesto de evento, la matriz y la referencia salen byte a byte iguales, y
// la verdad, la efectiva y el manifiesto salen iguales una vez quitados esos
// campos. Las huellas viejas quedan en el comentario de cada una.
//
// **01/11/2026 (entrega 4 del plan de detección):** las dos ejecuciones con
// loci publican además el catálogo de microsatélites aparte
// (`repeat_catalog.tsv` y `repeat_catalog.metadata.json`), enlazado con su
// huella desde la efectiva. Cambian a propósito la efectiva (gana
// `repeat_detection.catalog` y pierde `repeat_detection.retired_loci`, que
// pasa al metadato) y el manifiesto general (gana los dos artefactos en
// `artifacts` y en `output_categories.minimum`, y cambian la huella y el
// tamaño de la efectiva, que cita dos veces). Comprobado campo a campo contra
// el binario de `45af26b`, compilado aparte: ningún otro fichero cambia. Las
// huellas viejas quedan en el comentario de cada una.
//
// **01/11/2026 (entrega 4 bis del plan de detección):** las dos ejecuciones con
// loci publican además la tabla de alelos de los detectados
// (`repeat_alleles.tsv`), solo con su cabecera porque los dos loci son
// declarados. Cambian a propósito la efectiva (en `repeat_detection` gana el
// enlace `alleles_table`, las constantes de `frequency_normalization` y, en
// `slippage`, las versiones, el ámbito, la versión de su catálogo, las tres
// limitaciones y los linajes) y el manifiesto general (gana el artefacto en
// `artifacts` y en `output_categories.minimum`, y cambian la huella y el tamaño
// de la efectiva, que cita dos veces). Comprobado campo a campo contra el
// binario de `afda48b`, compilado aparte: ningún otro fichero cambia, y las
// ejecuciones amplificada sin loci y sin evento salen byte a byte iguales.
// ---------------------------------------------------------------------------

const UNAMPLIFIED_WITH_LOCI: [(&str, &str); 15] = [
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // 61809c4359781e56ac0d419c77c1014a43e83c42b77a34260fbc9e3f4659b975.
    // Y otra vez el 01/11/2026 (entrega 4 del plan de detección): era
    // f9d7b7dd5c2bf0a0461ae8e4322d692e637df97af4275fd1125450465ea564f2.
    // Y otra vez el 01/11/2026 (entrega 4 bis del plan de detección): era
    // b94d849a8c8fb66f668e67cff9140e6039d0c82294121dfca033029f361e3b6c.
    (
        "effective-config.json",
        "ae7e04b032ed1e28e676cb04321b5826c7b6f17133b5b60723ffc4fe55b0a975",
    ),
    (
        "events/reads-main/manifest.json",
        "75aedb1c798cc994e91684e9b0623a691e0c8d06b7cec9814fc20517ea4ecaed",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // 04b18080405ac7659985af85343f6d12bd89374626385f484dce9eef3d0b5112.
    (
        "inputs/config/general.toml",
        "1071fd4092a2e5266e701c5e711d92347ad9983c568b0a0decd4ac1f6e73e8d3",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // eaa1a6c1e9e2d2dc5db7cb9dea370951ae4a5ddad63984b0ed8b86c0635ca123.
    // Y otra vez el 01/11/2026 (entrega 4 del plan de detección): era
    // bb1498a051ae5363926b81fa9a055ee6a9f7273c32ded0dd070b9eb6fe4287dd.
    // Y otra vez el 01/11/2026 (entrega 4 bis del plan de detección): era
    // 7bfe4e2607848cf83b933f22ad19187cec4adbd823447c8eb18845f477fc281d.
    (
        "manifest.json",
        "445f2db1d767ae89fc19da8d430b1f137d2b4501e56aeec0329b86b27a3b3118",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "6109665fe5b539b71778c8256e379e64c88d4ac34c8e4f31e96890cbe814f649",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "0df5d5c66a2d958dd623d176c1a77290bf883613cc8b564e11e0bac6cc201264",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "51161276395b41e1bbb22268f668c3c0cf0f03ae46e6da6caf3a88e66096a970",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "a29364c401393696809b92c23ab1ca3089f478b0122baf0ebf7d63c3ccd67c5f",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // bca71bde0d2b215f56a11642814c9efd2eef45e1f2f172b15f51e261133ebab0.
    (
        "outputs/truth.jsonl",
        "23ed97f2ffa7c34d2a0a7d2f418b3ed87fd77ea1601fa95a8e89bac0650b165d",
    ),
    (
        "population_snvs.metadata.json",
        "c5becf380395dba899121805e2cd6f893ed2dbcf24da28fa3e94c84bf6356325",
    ),
    (
        "population_snvs.tsv",
        "ea6647904791e5af257036dcd4ae8ebe7b869e8f1c170789d907b7fd187be9ad",
    ),
    // Nueva el 01/11/2026: la tabla de alelos de los detectados (entrega 4
    // bis del plan de detección); solo la cabecera, porque los dos loci son
    // declarados.
    (
        "repeat_alleles.tsv",
        "c1dc0fe32075225c7d2d2cf9352ea8543c26033feeb54c2456ecc0b01b66ef6b",
    ),
    // Nuevos el 01/11/2026: el catálogo de microsatélites aparte (entrega 4
    // del plan de detección, REQ-MAN-98).
    (
        "repeat_catalog.metadata.json",
        "0f91d3d080249c7aa6b9ac02e06e6d84068969dce89c89faae1d45cbdbd89d26",
    ),
    (
        "repeat_catalog.tsv",
        "3b282ff10a6e39933e5c2c70456a81b08187fa852f7e0c1c12a2e47909552c29",
    ),
    (
        "resources/references/ref-main.fa",
        "227c733e42b142387655fcac85ff711d4a9ed698f5bca9a02cbed4a705980c8a",
    ),
];

const AMPLIFIED_WITHOUT_LOCI: [(&str, &str); 12] = [
    (
        "effective-config.json",
        "19c42001bbcb7500232179316c0b066f6564f39bb1cf3ebbb52721ccede183fa",
    ),
    (
        "events/reads-main/manifest.json",
        "63139a096b7007fd1c34123dbc3137038375bbd4e3a71097f0b0cc0da26e7629",
    ),
    (
        "inputs/config/general.toml",
        "0371a7b57a0d9eb4083f354c40c65ccb33de30204215068c305064c6c8780b2f",
    ),
    (
        "manifest.json",
        "ff0528a3587a62f10f32e3f61f1913430bb0c4ae956c23dd005900cbd620c30b",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "739f9c5aaa7c5275eeca6fb139101d1dd3ac294cd21ee0466bec64f8b9ca5339",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "0df5d5c66a2d958dd623d176c1a77290bf883613cc8b564e11e0bac6cc201264",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "54a0d8b90e16d747aef47b0d414c555cfe38cbcd969f282c76121f4cf89f9936",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "a29364c401393696809b92c23ab1ca3089f478b0122baf0ebf7d63c3ccd67c5f",
    ),
    (
        "outputs/truth.jsonl",
        "36246f62e4a311944235fae90637ad66dc44bbfa2941d44c0121c7cd591dcfc0",
    ),
    (
        "population_snvs.metadata.json",
        "c5becf380395dba899121805e2cd6f893ed2dbcf24da28fa3e94c84bf6356325",
    ),
    (
        "population_snvs.tsv",
        "ea6647904791e5af257036dcd4ae8ebe7b869e8f1c170789d907b7fd187be9ad",
    ),
    (
        "resources/references/ref-main.fa",
        "227c733e42b142387655fcac85ff711d4a9ed698f5bca9a02cbed4a705980c8a",
    ),
];

const WITHOUT_TECHNICAL_EVENT: [(&str, &str); 7] = [
    (
        "effective-config.json",
        "1b3bed2e56813c77554188f4f4cb9386eccc51ec86c95e583045710567949ec3",
    ),
    (
        "inputs/config/general.toml",
        "f911ab39f92c498831fa5baf0e3d7e40344a2ced69b370d36c2f790559d23cee",
    ),
    (
        "manifest.json",
        "3afc602e1a834df11f61b37c97d9f8bbfada1f6965dcdf504a62e7aedc59351e",
    ),
    (
        "outputs/truth.jsonl",
        "ae2a96f68f651840ea2a1dd46c81d85d8e1680f340e7f8cfc61fd83f917e7d18",
    ),
    (
        "population_snvs.metadata.json",
        "c5becf380395dba899121805e2cd6f893ed2dbcf24da28fa3e94c84bf6356325",
    ),
    (
        "population_snvs.tsv",
        "ea6647904791e5af257036dcd4ae8ebe7b869e8f1c170789d907b7fd187be9ad",
    ),
    (
        "resources/references/ref-main.fa",
        "227c733e42b142387655fcac85ff711d4a9ed698f5bca9a02cbed4a705980c8a",
    ),
];

const AMPLIFIED_WITH_LOCI: [(&str, &str); 15] = [
    // Cambia a propósito el 12/10/2026 (DT-040): era
    // 834848fea5745b4d8e7b7e0f41a329bcc4d50a7a91fc28dfb82fcf82392b6b64.
    // Y otra vez el 30/10/2026 (entrega 3 del plan de detección): era
    // b9c5280d983d94b50aa08c4cab453ba67e5a260521cd814361d6770272a0b37e.
    // Y otra vez el 01/11/2026 (entrega 4 del plan de detección): era
    // 44019117080453d9c9e86f31e5122c77a82a083248d797aebbbbedce9e9710eb.
    // Y otra vez el 01/11/2026 (entrega 4 bis del plan de detección): era
    // 5b4e563c9d5181d2986e7c250c03544d6893d93bf6a1c2430c2e7597f1522d17.
    (
        "effective-config.json",
        "3054b0c71cd52c9addc67e74e1fa1a37cf91e46e9e606acfe30cf4fb4cc51df8",
    ),
    (
        "events/reads-main/manifest.json",
        "bc59706c711f5bdcabefddc8f39b00ac4ddcda3cf998c4f26aabfd198ec30576",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // 4312136af5ee12661904fc6f2aa0720222925d11b7fdc2a5cc5a1597a39f36d7.
    (
        "inputs/config/general.toml",
        "b0b4d8db655594c7ab3dd0567711c4283f848c3fac43c11976eb83d2fdcb636e",
    ),
    // Cambia a propósito el 12/10/2026 (DT-040): era
    // a40a8cbb245653eacbbf04449d1db2b66a8877ac8e4ec971f80f95b91caf0de5.
    // Y otra vez el 30/10/2026 (entrega 3 del plan de detección): era
    // 679a61e0c7924f7256b049a41adeb80a5dda0aa4c1ff274c482d370acb1badf2.
    // Y otra vez el 01/11/2026 (entrega 4 del plan de detección): era
    // 8beaf3fcc5be556c1c78d6305dfb788fe57ccb833796f3044046816624bc962e.
    // Y otra vez el 01/11/2026 (entrega 4 bis del plan de detección): era
    // 9f8f6fddb5abbf2e4c8cd0dbcfc1a25d5b2fa9c0ba2d19785f603611c95b4b46.
    (
        "manifest.json",
        "13b6fae02fb3eaa7ec038b96e4010769e76e2fd88c21afcc0d63fd3666012178",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R1.fastq.gz",
        "fe3024e5da040cef13f9372c1db3eb23e592f97e5283fbd66f7c831fb5960fe8",
    ),
    (
        "outputs/fastq/reads-main/sample-A_R2.fastq.gz",
        "0df5d5c66a2d958dd623d176c1a77290bf883613cc8b564e11e0bac6cc201264",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R1.fastq.gz",
        "d5176776afce409db804d4ef6925261ac1d02a3169bfd6d912546ba41dd3a1be",
    ),
    (
        "outputs/fastq/reads-main/sample-B_R2.fastq.gz",
        "a29364c401393696809b92c23ab1ca3089f478b0122baf0ebf7d63c3ccd67c5f",
    ),
    // Cambia a propósito el 30/10/2026 (entrega 3 del plan de detección): era
    // 0ee74f98b072b70e88ed344656e250b3a978def3ac97aa04e20f8a3de2e2d7ae.
    (
        "outputs/truth.jsonl",
        "fa398a3d952478bc4318fb37e3ac0c0246475cb91ed90c2b2c84e47d2ec11154",
    ),
    (
        "population_snvs.metadata.json",
        "c5becf380395dba899121805e2cd6f893ed2dbcf24da28fa3e94c84bf6356325",
    ),
    (
        "population_snvs.tsv",
        "ea6647904791e5af257036dcd4ae8ebe7b869e8f1c170789d907b7fd187be9ad",
    ),
    // Nueva el 01/11/2026: la tabla de alelos de los detectados (entrega 4
    // bis del plan de detección); solo la cabecera, porque los dos loci son
    // declarados.
    (
        "repeat_alleles.tsv",
        "c1dc0fe32075225c7d2d2cf9352ea8543c26033feeb54c2456ecc0b01b66ef6b",
    ),
    // Nuevos el 01/11/2026: el catálogo de microsatélites aparte (entrega 4
    // del plan de detección, REQ-MAN-98).
    (
        "repeat_catalog.metadata.json",
        "0f91d3d080249c7aa6b9ac02e06e6d84068969dce89c89faae1d45cbdbd89d26",
    ),
    (
        "repeat_catalog.tsv",
        "3b282ff10a6e39933e5c2c70456a81b08187fa852f7e0c1c12a2e47909552c29",
    ),
    (
        "resources/references/ref-main.fa",
        "227c733e42b142387655fcac85ff711d4a9ed698f5bca9a02cbed4a705980c8a",
    ),
];

struct Published {
    root: PathBuf,
    destination: PathBuf,
}

impl Drop for Published {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Published {
    fn path(&self, relative: &str) -> PathBuf {
        self.destination.join(relative)
    }

    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&fs::read(self.path(relative)).expect("artefacto publicado"))
            .expect("JSON")
    }
}

/// Publica el escenario de la evidencia con ese nombre, leído por su ruta
/// relativa a la raíz del paquete.
fn publish(scenario: &str) -> Published {
    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
    assert_eq!(
        std::env::current_dir()
            .expect("directorio de trabajo")
            .canonicalize()
            .expect("ruta canónica"),
        manifest_dir.canonicalize().expect("ruta canónica"),
        "las rutas relativas de la evidencia exigen lanzar desde la raíz del paquete"
    );
    let root = std::env::temp_dir().join(format!(
        "gps-dt-040-{}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ));
    fs::create_dir_all(&root).expect("directorio");
    // El guardián que borra la carpeta existe antes de ejecutar, de modo que
    // un fallo de la publicación tampoco la deja atrás.
    let published = Published {
        destination: root.join("published"),
        root,
    };
    let toml = PathBuf::from(format!("{EVIDENCE}/{scenario}/general.toml"));
    execute_validated_config(
        &load_and_validate_config(&toml).expect("configuración válida"),
        &published.destination,
    )
    .expect("publicación");
    published
}

fn digest(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

fn external_schema(document: &Path, schema: &str) -> std::process::Output {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(document)
        .output()
        .expect("el validador externo se ejecuta")
}

fn assert_external_schema(document: &Path, schema: &str, valid: bool, what: &str) {
    let output = external_schema(document, schema);
    assert_eq!(
        output.status.success(),
        valid,
        "{what}: {schema} sobre {document:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn profile(amplified: bool) -> Value {
    if amplified {
        json!({"id": GPS_HISEQX_PCR_PE150_V1, "sha256": GPS_HISEQX_PCR_PE150_SHA256})
    } else {
        json!({"id": GPS_HISEQX_PCRFREE_PE150_V1, "sha256": GPS_HISEQX_PCRFREE_PE150_SHA256})
    }
}

/// **Cada tipo de ejecución publica el catálogo y el perfil que le
/// corresponden, y la efectiva y el manifiesto general pasan su esquema por
/// el validador externo** (DT-040, DT-041).
///
/// | Ejecución | Efectiva | Manifiesto general | Manifiesto de evento |
/// |---|---|---|---|
/// | sin PCR, con loci | v4, sin PCR | v6: v4, sin PCR | v3: v4 |
/// | amplificada sin loci | v5, con PCR | **v4**: v5, con PCR | v3: v5 |
/// | amplificada con loci | **v6**, con PCR | v6: **v6**, con PCR | v3: v6 |
///
/// **Debe fallar contra** el código anterior por tres sitios: la efectiva de
/// la amplificada con loci decía v5 y su esquema la rechazaba; su manifiesto
/// general decía v5; y los manifiestos generales de las dos amplificadas los
/// rechazaban los esquemas v4 y v6, que fijaban el perfil sin PCR. Hasta este
/// cambio **ninguna** prueba validaba por fuera el manifiesto general de una
/// ejecución amplificada.
#[test]
fn every_kind_of_run_publishes_documents_that_pass_their_schemas() {
    for (scenario, amplified, effective_catalog, manifest_version) in [
        ("sin-pcr-loci", false, "f4-technical-v4", "6"),
        ("amplificado-sin-loci", true, "f4-technical-v5", "4"),
        ("amplificado-loci", true, "f4-technical-v6", "6"),
    ] {
        let published = publish(scenario);
        let effective = published.json("effective-config.json");
        let manifest = published.json("manifest.json");
        let event = published.json("events/reads-main/manifest.json");
        assert_eq!(effective["schema_version"], "6", "{scenario}");
        assert_eq!(
            effective["technical_event"]["rng_catalog"], effective_catalog,
            "{scenario}: catálogo de la raíz del evento en la efectiva"
        );
        assert_eq!(
            effective["technical_event"]["profile"]["id"],
            profile(amplified)["id"],
            "{scenario}"
        );
        if amplified {
            assert_eq!(
                effective["technical_event"]["artifacts"]["pcr"]["rng_catalog"], effective_catalog,
                "{scenario}: la raíz y el bloque de la amplificación dicen lo mismo"
            );
        }
        assert_eq!(manifest["schema_version"], manifest_version, "{scenario}");
        assert_eq!(
            manifest["technical_events"][0]["rng_catalog"], effective_catalog,
            "{scenario}: catálogo del enlace en el manifiesto general"
        );
        assert_eq!(
            manifest["technical_events"][0]["profile"],
            profile(amplified),
            "{scenario}"
        );
        assert_eq!(event["schema_version"], "3", "{scenario}");
        assert_eq!(event["rng_catalog"], effective_catalog, "{scenario}");
        assert_external_schema(
            &published.path("effective-config.json"),
            "effective-config-v6.schema.json",
            true,
            scenario,
        );
        assert_external_schema(
            &published.path("manifest.json"),
            &format!("manifest-general-v{manifest_version}.schema.json"),
            true,
            scenario,
        );
        assert_external_schema(
            &published.path("events/reads-main/manifest.json"),
            "manifest-event-v3.schema.json",
            true,
            scenario,
        );
    }
}

/// **Lo que no tenía que moverse no se mueve, y lo que se mueve es exactamente
/// lo decidido.** Cada fichero de los cuatro paquetes contra las huellas del
/// binario de partida (ver la nota de arriba). Las ejecuciones sin PCR, sin
/// evento técnico y amplificadas sin loci son idénticas byte a byte; en la
/// amplificada con loci cambian solo la efectiva y el manifiesto general.
///
/// **Debe fallar contra** un arreglo que publicara el v6 en toda biblioteca
/// amplificada (cambiarían la efectiva y el manifiesto general de la
/// amplificada sin loci) y contra cualquier otro cambio de bytes. Contra el
/// código anterior falla por las dos huellas nuevas de la amplificada con
/// loci. El conjunto de ficheros también se compara: un fichero de más o de
/// menos rompe la prueba.
#[test]
fn published_packages_keep_their_frozen_digests() {
    for (scenario, frozen) in [
        ("sin-pcr-loci", &UNAMPLIFIED_WITH_LOCI[..]),
        ("amplificado-sin-loci", &AMPLIFIED_WITHOUT_LOCI[..]),
        ("sin-evento", &WITHOUT_TECHNICAL_EVENT[..]),
        ("amplificado-loci", &AMPLIFIED_WITH_LOCI[..]),
    ] {
        let published = publish(scenario);
        let mut files = Vec::new();
        let mut pending = vec![published.destination.clone()];
        while let Some(directory) = pending.pop() {
            for entry in fs::read_dir(&directory).expect("directorio publicado") {
                let path = entry.expect("entrada").path();
                if path.is_dir() {
                    pending.push(path);
                } else {
                    files.push(
                        path.strip_prefix(&published.destination)
                            .expect("dentro del paquete")
                            .to_string_lossy()
                            .into_owned(),
                    );
                }
            }
        }
        files.sort();
        assert_eq!(
            files,
            frozen
                .iter()
                .map(|(relative, _)| (*relative).to_owned())
                .collect::<Vec<_>>(),
            "{scenario}: ficheros publicados"
        );
        for (relative, expected) in frozen {
            assert_eq!(
                digest(&fs::read(published.path(relative)).expect("artefacto")),
                *expected,
                "{scenario}/{relative} ha cambiado de bytes; NO se actualiza el literal sin \
                 decidirlo por escrito"
            );
        }
    }
}

type Edit = Box<dyn Fn(&mut Value)>;

/// Un caso: qué se altera, en qué paquete y documento, y qué deben decir el
/// juez completo (motivo exacto) y el esquema (si lo acepta).
type Case<'a> = (&'a str, &'a Published, &'a str, Vec<Edit>, &'a str, bool);

fn set_catalog(pointer: &'static str, catalog: &'static str) -> Edit {
    Box::new(move |document: &mut Value| {
        *document.pointer_mut(pointer).expect("campo") = json!(catalog);
    })
}

fn set_profile(amplified: bool) -> Edit {
    Box::new(move |document: &mut Value| {
        document["technical_events"][0]["profile"] = profile(amplified);
    })
}

/// Lo que dicen el juez completo y el esquema externo de un documento
/// alterado. Se escribe la copia alterada fuera del paquete para el
/// validador externo; el juez la recibe en bytes.
fn judge_and_schema(published: &Published, target: &str, edits: &[&Edit]) -> (String, bool) {
    let mut effective = published.json("effective-config.json");
    let mut manifest = published.json("manifest.json");
    let document = if target == "effective-config.json" {
        &mut effective
    } else {
        &mut manifest
    };
    for edit in edits {
        edit(document);
    }
    let altered = published.root.join(format!("altered-{target}"));
    let bytes = |value: &Value| serde_json::to_vec_pretty(value).expect("JSON");
    fs::write(&altered, bytes(document)).expect("copia alterada");
    let schema = if target == "effective-config.json" {
        "effective-config-v6.schema.json".to_owned()
    } else {
        format!(
            "manifest-general-v{}.schema.json",
            manifest["schema_version"].as_str().expect("versión")
        )
    };
    let schema_accepts = external_schema(&altered, &schema).status.success();
    let verdict = validate_published_resources(
        &published.destination,
        &published.path("effective-config.json"),
        &bytes(&effective),
        &published.path("manifest.json"),
        &bytes(&manifest),
        None,
    );
    (describe(verdict), schema_accepts)
}

fn describe(result: Result<(), Diagnostic>) -> String {
    match result {
        Ok(()) => "OK".to_owned(),
        Err(error) => error
            .context()
            .get("detail")
            .cloned()
            .unwrap_or_else(|| format!("{error:?}")),
    }
}

const EFFECTIVE_CATALOG: &str = "catálogo técnico F4 inválido";
const MANIFEST_CATALOG: &str = "el catálogo de azar del inventario técnico no corresponde a su perfil o a la versión del manifiesto";
const MANIFEST_VERSUS_EFFECTIVE: &str =
    "el catálogo de azar del manifiesto general no es el de la configuración efectiva";

/// **Los jueces y los esquemas rechazan un catálogo de azar que contradice al
/// perfil o al tartamudeo clonal**, con el motivo exacto.
///
/// En la configuración efectiva el catálogo de la raíz del evento se ata a lo
/// que el evento abre (`validate_effective_technical_event`): con PCR y
/// tartamudeo clonal, v6; con PCR sin él, v5; sin PCR, v3 o v4. En el
/// manifiesto general, el lector ata el catálogo al perfil y a la versión del
/// documento (`validate_manifest_technical_events`), y el juez completo lo
/// compara además con el de la efectiva (DT-056).
///
/// **Debe fallar contra** el código anterior: el juez de la efectiva exigía
/// el v5 con tartamudeo clonal (rechazaba el v6 y aceptaba el v5), y el lector
/// del manifiesto general no ataba el catálogo al perfil.
///
/// **Con PCR, la versión 6 solo admite el v6.** Esa versión solo la publica
/// una ejecución con loci de repetición, y con PCR y loci siempre hay
/// tartamudeo clonal (la configuración exige el de base), así que el v5 de la
/// amplificación sin loci es exclusivo de la versión 4: el esquema y el lector
/// rechazan un v6 con PCR que diga v5 (corregido en la revisión del
/// 12/10/2026; la primera versión de esta prueba lo daba por aceptable).
///
/// **Lo que el esquema sí no puede ver.** Sin PCR, la versión 6 admite el v3 y
/// el v4, y el manifiesto no dice si el tartamudeo de base está encendido: un
/// enlace en v3 frente a una efectiva en v4 pasa el esquema con razón. Lo caza
/// el juez completo al compararlo con la efectiva (DT-056); el caso lo dice tal
/// cual.
#[test]
fn the_judges_and_the_schemas_reject_a_catalog_that_contradicts_the_run() {
    let unamplified = publish("sin-pcr-loci");
    let amplified_without_loci = publish("amplificado-sin-loci");
    let amplified_with_loci = publish("amplificado-loci");
    for (label, published) in [
        ("sin PCR", &unamplified),
        ("amplificada sin loci", &amplified_without_loci),
        ("amplificada con loci", &amplified_with_loci),
    ] {
        assert_eq!(
            judge_and_schema(published, "manifest.json", &[]),
            ("OK".to_owned(), true),
            "control positivo: {label}"
        );
        assert_eq!(
            judge_and_schema(published, "effective-config.json", &[]),
            ("OK".to_owned(), true),
            "control positivo: {label}"
        );
    }
    const ROOT: &str = "/technical_event/rng_catalog";
    const LINK: &str = "/technical_events/0/rng_catalog";
    let cases: Vec<Case<'_>> = vec![
        // La configuración efectiva.
        (
            "efectiva con tartamudeo clonal y v5 (DT-040)",
            &amplified_with_loci,
            "effective-config.json",
            vec![set_catalog(ROOT, "f4-technical-v5")],
            EFFECTIVE_CATALOG,
            false,
        ),
        (
            "efectiva amplificada sin tartamudeo clonal y v6",
            &amplified_without_loci,
            "effective-config.json",
            vec![set_catalog(ROOT, "f4-technical-v6")],
            EFFECTIVE_CATALOG,
            false,
        ),
        (
            "efectiva con PCR y v3",
            &amplified_without_loci,
            "effective-config.json",
            vec![set_catalog(ROOT, "f4-technical-v3")],
            EFFECTIVE_CATALOG,
            false,
        ),
        (
            "efectiva sin PCR y v5",
            &unamplified,
            "effective-config.json",
            vec![set_catalog(ROOT, "f4-technical-v5")],
            EFFECTIVE_CATALOG,
            false,
        ),
        (
            "efectiva sin PCR y v6",
            &unamplified,
            "effective-config.json",
            vec![set_catalog(ROOT, "f4-technical-v6")],
            EFFECTIVE_CATALOG,
            false,
        ),
        // El manifiesto general v4 (amplificada sin loci).
        (
            "manifiesto v4 con PCR y v3",
            &amplified_without_loci,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v3")],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v4 sin PCR y v5",
            &amplified_without_loci,
            "manifest.json",
            vec![set_profile(false)],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v4 con PCR y v6, que solo existe con loci",
            &amplified_without_loci,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v6")],
            MANIFEST_CATALOG,
            false,
        ),
        // El manifiesto general v6.
        (
            "manifiesto v6 sin PCR y v5",
            &unamplified,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v5")],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v6 sin PCR y v6",
            &unamplified,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v6")],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v6 con PCR y v3",
            &amplified_with_loci,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v3")],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v6 con PCR y v4",
            &amplified_with_loci,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v4")],
            MANIFEST_CATALOG,
            false,
        ),
        (
            "manifiesto v6 con PCR y v5, que es exclusivo de la versión 4",
            &amplified_with_loci,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v5")],
            MANIFEST_CATALOG,
            false,
        ),
        // Coherente con el perfil y con la versión, pero no con la efectiva
        // (v4): solo el juez completo lo ve, y el esquema lo acepta con razón
        // (ver el docstring).
        (
            "manifiesto v6 sin PCR en v3 frente a la efectiva en v4",
            &unamplified,
            "manifest.json",
            vec![set_catalog(LINK, "f4-technical-v3")],
            MANIFEST_VERSUS_EFFECTIVE,
            true,
        ),
    ];
    for (label, published, target, edits, reason, schema_accepts) in &cases {
        let edits = edits.iter().collect::<Vec<_>>();
        assert_eq!(
            judge_and_schema(published, target, &edits),
            ((*reason).to_owned(), *schema_accepts),
            "{label}"
        );
    }
}
