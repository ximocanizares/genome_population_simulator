//! F5.6.1. La cadena productiva de `run` y su obligación de proceso.
//!
//! Lo que esta batería fija no es que F5 sepa publicar —eso lo cierran F5.4 y
//! F5.5— sino que **el orquestador** que F5.6.1 añade publique por donde la
//! Revisión F5.5 exige y componga los dos auditores. La guarda de la
//! obligación es **léxica** y va con control positivo: sin él, una búsqueda mal
//! escrita pasaría por no encontrar nada.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{
    CliValidatedConfig, EffectiveF5Config, execute_validated_f5_config,
    execute_validated_f5_config_with_recount, f5_execution_id, load_and_validate_cli_config,
};
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::identity::derive_execution_id;
use genome_population_simulator::integrity::configuration_digest;
use genome_population_simulator::randomness::RootSeed;
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-cli-run-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Escribe una autoría F5 factible y devuelve la ruta del general.
fn write_pair(root: &Path) -> PathBuf {
    let scenario_directory = root.join("scenarios");
    let profile_directory = scenario_directory.join("profiles");
    fs::create_dir_all(&profile_directory).expect("directorio de perfil");
    fs::create_dir_all(root.join("references")).expect("directorio de referencia");

    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n",
        ),
        (
            "snp_gaps.tsv",
            "f5_snp_gap_profile",
            "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n",
        ),
        (
            "snp_spectrum.tsv",
            "f5_snp_spectrum",
            "alt_af_min\talt_af_max\tweight\n0.25\t0.75\t1\n",
        ),
        (
            "haplotype_weights.tsv",
            "f5_haplotype_weights",
            "haplotype_rank\tweight\n0\t0.5\n1\t0.5\n",
        ),
        (
            "ld_profile.tsv",
            "f5_ld_profile",
            "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t1\t0.05\nP1\t2\t2\t1\t0.05\n",
        ),
    ];
    let links = files
        .iter()
        .map(|(name, kind, contents)| {
            fs::write(profile_directory.join(name), contents).expect("TSV de perfil");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(contents.as_bytes())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(profile_directory.join("profile.metadata.json"), &metadata).expect("metadata");

    let fasta = ">chr1\nACGTAC\n";
    fs::write(root.join("references/ref.fa"), fasta).expect("referencia");
    let general_path = root.join("general.toml");
    fs::write(
        &general_path,
        format!(
            "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{SEED}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
            sha256(fasta.as_bytes())
        ),
    )
    .expect("general");
    fs::write(
        scenario_directory.join("scenario.toml"),
        format!(
            "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = 1\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = 32\nmax_swap_proposals = 256\nmax_haplotype_cells = 65536\nblock_count = 2\n\n[[subpopulations]]\npopulation_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n",
            sha256(metadata.as_bytes())
        ),
    )
    .expect("escenario");
    general_path
}

fn admitted(general: &Path) -> Box<EffectiveF5Config> {
    match load_and_validate_cli_config(general).expect("par F5 admitido") {
        CliValidatedConfig::F5(effective) => effective,
        CliValidatedConfig::Legacy(_) => panic!("el dispatcher tenía que clasificar un general v3"),
    }
}

/// La cadena productiva publica y **los dos auditores** cuentan el paquete.
///
/// El informe de clausura solo lo produce `audit_f5_identity_closure`, y los
/// resúmenes MAN-70 del manifiesto solo los produce `audit_f5_package`: si el
/// orquestador no los compusiera, ninguna de las dos cifras existiría. Las
/// cifras se contrastan además contra el propio árbol publicado, que es el
/// término independiente: el recuento de filas de `f5-founders.tsv` no sale
/// del informe.
#[test]
fn the_productive_chain_publishes_and_both_auditors_count_the_package() {
    let root = temporary("published");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    let destination = root.join("package");
    let outcome =
        execute_validated_f5_config(&admitted(&general), &destination).expect("paquete publicado");

    assert_eq!(outcome.destination(), destination);
    let closure = outcome.closure();
    // Verdad independiente, no «mayor que cero»: el fixture declara una
    // subpoblación con un individuo, ploidía 2 y una sola secuencia poblada,
    // de modo que las copias son individuos × secuencias × ploidía = 2 y el
    // individuo lo aporta `founder_ids`, no una `founder_key`.
    assert_eq!(closure.individuals(), 1);
    assert_eq!(closure.copies(), 2);
    assert_eq!(closure.supplied_individuals(), 1);
    assert_eq!(closure.derived_individuals(), 0);

    // Término independiente: las filas del TSV publicado, leídas aquí.
    let founders = fs::read_to_string(destination.join("outputs/f5-founders.tsv"))
        .expect("f5-founders.tsv publicado");
    let rows = founders
        .lines()
        .skip(1)
        .filter(|line| !line.is_empty())
        .count();
    assert_eq!(u64::try_from(rows).expect("filas"), closure.assignments());

    // `audit_f5_package` es el recuento MAN-70 que la publicación exige: sus
    // resúmenes viajan al manifiesto y no hay otra vía de producirlos.
    let manifest =
        fs::read_to_string(destination.join("manifest.json")).expect("manifiesto publicado");
    assert!(
        manifest.contains("\"initialization_summaries\""),
        "{manifest}"
    );
    let document: serde_json::Value = serde_json::from_str(&manifest).expect("manifiesto JSON");
    let summaries = document["f5"]["initialization_summaries"]
        .as_array()
        .expect("resúmenes MAN-70");
    assert!(!summaries.is_empty(), "{manifest}");

    // El efectivo publica el nombre del documento de autoría, no su ruta.
    let effective =
        fs::read_to_string(destination.join("effective-config.json")).expect("efectivo publicado");
    assert!(
        effective.contains("\"source_path\": \"general.toml\""),
        "{effective}"
    );
    fs::remove_dir_all(&root).expect("limpieza");
}

/// Dos ejecuciones del mismo par producen el mismo paquete byte a byte.
///
/// La identidad de ejecución de F5 se deriva de los SHA-256 de los **dos**
/// documentos: si dependiera solo del general, el escenario podría cambiar sin
/// cambiar `block_id` ni `site_id`. La segunda mitad de la prueba lo fija.
#[test]
fn the_execution_identity_covers_both_authored_documents() {
    let root = temporary("identity");
    fs::create_dir_all(&root).expect("raíz");
    let first_root = root.join("first");
    let second_root = root.join("second");
    let changed_root = root.join("changed");
    for copy in [&first_root, &second_root, &changed_root] {
        fs::create_dir_all(copy).expect("copia");
        write_pair(copy);
    }
    let fingerprint = |authoring: &Path, label: &str| -> String {
        let destination = authoring.join(label);
        let general = authoring.join("general.toml");
        execute_validated_f5_config(&admitted(&general), &destination).expect("paquete");
        sha256(
            &fs::read(destination.join("outputs/f5-block-layout.tsv")).expect("layout publicado"),
        )
    };
    let first = fingerprint(&first_root, "package");
    let second = fingerprint(&second_root, "package");
    assert_eq!(first, second, "el mismo par tiene que reproducirse");

    // Solo cambia el escenario —una subpoblación más, el mismo general—, y el
    // paquete tiene que dejar de coincidir.
    let scenario_path = changed_root.join("scenarios/scenario.toml");
    let scenario = fs::read_to_string(&scenario_path).expect("escenario");
    let mutated = scenario
        .replace("founder_count = 1", "founder_count = 2")
        .replace(
            "population_id = \"P1\"\nindividual_count = 1\nfounder_ids = [\"I1\"]\n",
            "population_id = \"P1\"\nindividual_count = 2\nfounder_ids = [\"I1\", \"I2\"]\n",
        );
    assert_ne!(mutated, scenario, "la mutación tiene que mutar");
    fs::write(&scenario_path, &mutated).expect("escenario mutado");
    let changed = fingerprint(&changed_root, "package");
    assert_ne!(
        first, changed,
        "otro escenario no puede producir el mismo layout bajo el mismo general"
    );
    fs::remove_dir_all(&root).expect("limpieza");
}

/// Vector del material del `execution_id`, y su rederivación independiente.
///
/// F1 a F4 rehacen la identidad de ejecución desde semilla y digest en
/// `src/results.rs`; el oráculo de F5 solo comprueba coherencia interna y
/// **nunca** la recalcula, de modo que intercambiar los dos términos del
/// material pasaba desapercibido. Aquí la rederivación se hace con las
/// primitivas públicas y sobre los bytes leídos del disco, y el orden queda
/// fijado por un vector con el identificador escrito a mano.
#[test]
fn the_execution_identity_is_the_seed_over_the_general_then_scenario_digests() {
    // Vector **congelado**: dos SHA-256 literales y el `execution_id` que
    // producen. Su valor no es derivable a mano —media BLAKE3 por medio—, de
    // modo que lo que demuestra no es la fórmula sino su **estabilidad**: si
    // el material, su orden o el camino de derivación cambian, este literal
    // deja de coincidir. La fórmula la fija la rederivación independiente de
    // más abajo, con las mismas primitivas públicas que `src/results.rs` usa
    // para rehacer la identidad de F1 a F4.
    let general_digest = "a".repeat(64);
    let scenario_digest = "b".repeat(64);
    let seed = RootSeed::parse("seed", SEED).expect("semilla canónica");
    let identity = |first: &str, second: &str| -> String {
        let material = format!("{first}{second}");
        derive_execution_id(
            seed,
            &configuration_digest(material.as_bytes()).expect("digest"),
        )
        .expect("identidad")
        .as_str()
        .to_owned()
    };
    assert_eq!(
        identity(&general_digest, &scenario_digest),
        "RUN_ZO2WZO4TKOGRM5ZDEOFVYASDJ4"
    );
    // El orden no es indiferente, y el vector lo demuestra en el mismo acto.
    assert_ne!(
        identity(&general_digest, &scenario_digest),
        identity(&scenario_digest, &general_digest)
    );

    // Y la configuración real usa **ese** material: general y después
    // escenario, sobre los bytes que el disco conserva.
    let root = temporary("execution-material");
    fs::create_dir_all(&root).expect("raíz");
    let general_path = write_pair(&root);
    let config = admitted(&general_path);
    let observed = f5_execution_id(&config).expect("identidad de la configuración");
    let general_bytes = fs::read(&general_path).expect("general");
    let scenario_bytes = fs::read(root.join("scenarios/scenario.toml")).expect("escenario");
    let expected = identity(&sha256(&general_bytes), &sha256(&scenario_bytes));
    assert_eq!(observed.as_str(), expected);
    // Control negativo con la mutación verificada: el orden contrario da otra.
    let swapped = identity(&sha256(&scenario_bytes), &sha256(&general_bytes));
    assert_ne!(expected, swapped);
    fs::remove_dir_all(&root).expect("limpieza");
}

/// Un fallo controlado de la etapa `publication` publica el paquete `failed`.
///
/// Recorre la rama que el éxito nunca toca: el provisional vivo que
/// `publish_f5_success` devuelve, `seal_and_publish_f5_failure`,
/// `audit_f5_failure_closure` y el diagnóstico del informe como valor de
/// retorno. El recuento MAN-70 se inyecta rechazando para alcanzarla sin
/// corromper el provisional desde fuera del proceso.
#[test]
fn a_controlled_publication_failure_publishes_a_failed_package() {
    let root = temporary("failed-publication");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    let destination = root.join("package");
    let config = admitted(&general);
    let rejection = Diagnostic::f5_configuration(
        "outputs/f5-founders.tsv",
        "man70",
        "el recuento MAN-70 no corresponde a la verdad publicada",
    );
    let returned = execute_validated_f5_config_with_recount(&config, &destination, {
        let rejection = rejection.clone();
        move |_payload, _declaration| Err(rejection)
    })
    .expect_err("el recuento rechazado tiene que tumbar la publicación");

    // El destino existe y es un paquete `failed`.
    let manifest = fs::read_to_string(destination.join("manifest.json")).expect("manifiesto");
    let document: serde_json::Value = serde_json::from_str(&manifest).expect("manifiesto JSON");
    assert_eq!(document["state"], "failed");
    assert_eq!(document["f5"]["f5_failure"]["stage"], "publication");
    assert!(destination.join("outputs/f5-failure.json").is_file());

    // Ninguno de los ocho ficheros biológicos llega al paquete fallido.
    for biological in [
        "outputs/f5-block-layout.tsv",
        "outputs/f5-sites.tsv",
        "outputs/f5-haplotypes.tsv",
        "outputs/f5-frequencies.tsv",
        "outputs/f5-ld-evaluation.tsv",
        "outputs/f5-founders.tsv",
        "outputs/f5-attempts.tsv",
        "outputs/truth.jsonl",
    ] {
        assert!(
            !destination.join(biological).exists(),
            "un paquete fallido no puede llevar {biological}"
        );
    }

    // El diagnóstico devuelto es el del informe, no el de la clausura.
    assert_eq!(returned.code(), rejection.code());
    assert_eq!(
        returned.context().get("detail"),
        rejection.context().get("detail")
    );
    fs::remove_dir_all(&root).expect("limpieza");
}

/// La búsqueda haplotípica agotada publica `failed` con su propia etapa.
///
/// DTC-045 solo exime de publicar al rechazo **anterior** al comienzo, y estas
/// etapas ya han abierto ámbitos RNG derivados del `execution_id`. Antes de la
/// Revisión F5.6.1 devolvían su diagnóstico sin destino, de modo que
/// `F5FailureStage` declaraba cinco variantes con una sola alcanzable.
#[test]
fn an_exhausted_haplotype_search_publishes_a_failed_package() {
    let root = temporary("failed-search");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    // El mismo par con un objetivo de LD que la búsqueda no puede alcanzar.
    let profile = root.join("scenarios/profiles/ld_profile.tsv");
    let infeasible = "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\nP1\t1\t1\t0.5\t0.05\nP1\t2\t2\t0.5\t0.05\n";
    let previous = fs::read_to_string(&profile).expect("perfil LD");
    assert_ne!(previous, infeasible, "la mutación tiene que mutar");
    fs::write(&profile, infeasible).expect("perfil LD inviable");
    rewrite_profile_metadata(&root);

    let destination = root.join("package");
    let config = admitted(&general);
    let returned = execute_validated_f5_config(&config, &destination)
        .expect_err("la búsqueda agotada no puede publicar un éxito");
    assert_eq!(returned.code(), "GPS-DOMAIN-104");

    let manifest = fs::read_to_string(destination.join("manifest.json")).expect("manifiesto");
    let document: serde_json::Value = serde_json::from_str(&manifest).expect("manifiesto JSON");
    assert_eq!(document["state"], "failed");
    assert_eq!(document["f5"]["f5_failure"]["stage"], "haplotypes");
    assert_eq!(document["f5"]["f5_failure"]["population_id"], "P1");
    assert!(destination.join("outputs/f5-failure.json").is_file());
    assert!(!destination.join("outputs/f5-haplotypes.tsv").exists());
    assert!(!destination.join("outputs/truth.jsonl").exists());
    fs::remove_dir_all(&root).expect("limpieza");
}

/// El nombre del general se publica: la regla vive en la admisión, una vez.
///
/// Antes vivía en el orquestador y en el adaptador Python, de modo que
/// `validate` aceptaba un nombre que `run` rechazaba después.
#[test]
fn an_unpublishable_general_name_is_rejected_by_validate_and_by_run() {
    let root = temporary("unpublishable-name");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    let renamed = root.join("general con espacios.toml");
    fs::rename(&general, &renamed).expect("renombrado");

    let rejection = load_and_validate_cli_config(&renamed)
        .expect_err("el nombre no publicable se rechaza en la admisión");
    assert_eq!(rejection.code(), "GPS-CONF-101");
    assert_eq!(
        rejection.context().get("fields").map(String::as_str),
        Some("general_path")
    );
    // Control positivo: el mismo par con su nombre publicable sí se admite.
    fs::rename(&renamed, &general).expect("renombrado de vuelta");
    assert!(load_and_validate_cli_config(&general).is_ok());
    fs::remove_dir_all(&root).expect("limpieza");
}

/// Reescribe el metadata del perfil tras cambiar uno de sus cinco TSV.
fn rewrite_profile_metadata(root: &Path) {
    let directory = root.join("scenarios/profiles");
    let kinds = [
        ("block_lengths.tsv", "f5_block_length_profile"),
        ("snp_gaps.tsv", "f5_snp_gap_profile"),
        ("snp_spectrum.tsv", "f5_snp_spectrum"),
        ("haplotype_weights.tsv", "f5_haplotype_weights"),
        ("ld_profile.tsv", "f5_ld_profile"),
    ];
    let links = kinds
        .iter()
        .map(|(name, kind)| {
            let contents = fs::read(directory.join(name)).expect("TSV de perfil");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(&contents)
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(directory.join("profile.metadata.json"), &metadata).expect("metadata");
    let scenario_path = root.join("scenarios/scenario.toml");
    let scenario = fs::read_to_string(&scenario_path).expect("escenario");
    let previous = scenario
        .lines()
        .find(|line| line.starts_with("profile_metadata_sha256"))
        .expect("hash del perfil")
        .to_owned();
    let rewritten = scenario.replace(
        &previous,
        &format!(
            "profile_metadata_sha256 = \"{}\"",
            sha256(metadata.as_bytes())
        ),
    );
    assert_ne!(rewritten, scenario, "la reescritura tiene que reescribir");
    fs::write(&scenario_path, rewritten).expect("escenario reescrito");
}

/// Las dos autorizaciones estrechas de F2 se rechazan en F5, no se descartan.
///
/// Eran inocuas mientras `run` F5 siempre fallaba: ninguna vía podía
/// consumirlas y ninguna ejecución existía. Con la cadena productiva
/// encadenada, aceptarlas en silencio dejaría creer que una ejecución quedó
/// autorizada cuando ninguna autorización se aplicó. La prueba invoca el
/// **binario**, que es donde vive la frontera de argumentos.
///
/// El control positivo de la otra mitad —que las rutas F1–F4 **siguen**
/// aceptándolas— no se repite aquí: lo cubren
/// `empty_population_matrix_requires_and_records_its_explicit_authorization` y
/// `direct_events_are_preflighted_before_publication_and_the_narrow_flag_is_traced`
/// en `tests/cli_smoke.rs`, que exigen éxito **y** que el efectivo y el
/// manifiesto publiquen `authorizations.*` en `true`.
#[test]
fn the_narrow_f2_authorizations_are_rejected_by_f5_validate_and_run() {
    let root = temporary("legacy-authorizations");
    fs::create_dir_all(&root).expect("raíz");
    let general = write_pair(&root);
    let binary = env!("CARGO_BIN_EXE_genome-population-simulator");
    let general_argument = general.to_str().expect("ruta del general");

    for flags in [
        vec!["--allow-empty-snv-matrix"],
        vec!["--allow-incompatible-genotypes"],
        vec!["--allow-empty-snv-matrix", "--allow-incompatible-genotypes"],
    ] {
        let validation = Command::new(binary)
            .args(["validate", general_argument])
            .args(&flags)
            .output()
            .expect("validate F5 con autorizaciones ajenas");
        let rendered = String::from_utf8_lossy(&validation.stderr).into_owned();
        assert!(!validation.status.success(), "{rendered}");
        assert!(rendered.contains("GPS-CONF-101"), "{rendered}");
        for flag in &flags {
            assert!(rendered.contains(flag), "{rendered}");
        }

        let destination = root.join(format!("run-{}", flags.len()));
        let execution = Command::new(binary)
            .args(["run", general_argument])
            .args(&flags)
            .arg("--output")
            .arg(&destination)
            .output()
            .expect("run F5 con autorizaciones ajenas");
        let rendered = String::from_utf8_lossy(&execution.stderr).into_owned();
        assert!(!execution.status.success(), "{rendered}");
        assert!(rendered.contains("GPS-CONF-101"), "{rendered}");
        for flag in &flags {
            assert!(rendered.contains(flag), "{rendered}");
        }
        // El rechazo precede al provisional: ni destino ni paquete `failed`.
        assert!(
            !destination.exists(),
            "una autorización ajena no puede abrir destino"
        );
    }

    // Control positivo: el mismo par, sin las dos banderas, valida y publica.
    let clean = Command::new(binary)
        .args(["validate", general_argument])
        .output()
        .expect("validate F5 sin autorizaciones");
    assert!(
        clean.status.success(),
        "{}",
        String::from_utf8_lossy(&clean.stderr)
    );
    let destination = root.join("run-clean");
    let published = Command::new(binary)
        .args(["run", general_argument, "--output"])
        .arg(&destination)
        .output()
        .expect("run F5 sin autorizaciones");
    assert!(
        published.status.success(),
        "{}",
        String::from_utf8_lossy(&published.stderr)
    );
    assert!(destination.join("outputs/truth.jsonl").is_file());
    fs::remove_dir_all(&root).expect("limpieza");
}

/// Obligación de proceso de la Revisión F5.5: el orquestador publica
/// **exclusivamente** por los envoltorios que sellan.
///
/// Llamar a `publish_f5_success` o a `publish_f5_failure` daría un destino sin
/// testigo y dejaría la clausura de F5.5.1 inaplicable **sin ninguna señal**.
/// La guarda es léxica y por eso lleva dos controles positivos: uno prueba que
/// el orquestador sí nombra los envoltorios, y otro que la misma búsqueda
/// **encuentra** el nombre desnudo donde realmente está.
#[test]
fn the_orchestrator_publishes_only_through_the_sealing_wrappers() {
    let orchestrator = include_str!("../src/application.rs");
    let cli = include_str!("../src/main.rs");
    for (module, source) in [("application.rs", orchestrator), ("main.rs", cli)] {
        for forbidden in ["publish_f5_success(", "publish_f5_failure("] {
            assert_eq!(
                bare_call_positions(source, forbidden),
                Vec::<usize>::new(),
                "{module} no puede publicar por {forbidden} fuera de los envoltorios que sellan"
            );
        }
    }
    // Control positivo 1: el orquestador sí publica, y por donde debe.
    assert!(orchestrator.contains("seal_and_publish_f5_success("));
    assert!(orchestrator.contains("seal_and_publish_f5_failure("));
    // Control positivo 2: la búsqueda encuentra el nombre desnudo donde vive.
    let wrappers = include_str!("../src/f5_closure.rs");
    assert!(
        !bare_call_positions(wrappers, "publish_f5_success(").is_empty(),
        "la búsqueda negativa no discrimina si no encuentra la llamada real de f5_closure.rs"
    );
    assert!(
        !bare_call_positions(wrappers, "publish_f5_failure(").is_empty(),
        "la búsqueda negativa no discrimina si no encuentra la llamada real de f5_closure.rs"
    );
}

/// El orquestador compone los dos auditores, y se lee en su texto.
#[test]
fn the_orchestrator_composes_both_audits() {
    let orchestrator = include_str!("../src/application.rs");
    for auditor in [
        "audit_f5_package(",
        "audit_f5_identity_closure(",
        "audit_f5_failure_closure(",
    ] {
        assert!(
            orchestrator.contains(auditor),
            "el orquestador tiene que componer {auditor}"
        );
    }
    // Control positivo **con contenido**: la guarda anterior protegía tres
    // `assert!(contains(x))`, que no pueden pasar en vacío. Lo que hace falta
    // es que la misma búsqueda distinga un nombre presente y cercano de uno
    // ausente y cercano.
    assert!(orchestrator.contains("audit_f5_identity_closure(destination"));
    assert!(!orchestrator.contains("audit_f5_identity_closure(payload"));
    assert!(!orchestrator.contains("audit_f5_recombination("));
}

/// Posiciones de `needle` **no** precedidas por `seal_and_`.
///
/// `seal_and_publish_f5_success(` contiene `publish_f5_success(`: una búsqueda
/// ingenua daría siempre positivo y la guarda no diría nada. Esta distingue.
fn bare_call_positions(source: &str, needle: &str) -> Vec<usize> {
    source
        .match_indices(needle)
        .filter(|(index, _)| !source[..*index].ends_with("seal_and_"))
        .map(|(index, _)| index)
        .collect()
}
