//! F4.9.1 y F4.9.2: el perfil hermano y el juez de las constantes de contrato.
//!
//! # Qué vigila este fichero
//!
//! 1. Que el recurso del perfil **sin PCR no cambie ni en un byte**. Su nombre
//!    es literal en tres esquemas publicados y su huella está congelada: meter
//!    las dos preparaciones en un solo recurso rompería los dos.
//! 2. Que las cinco constantes de contrato de la amplificación tengan **una
//!    sola** definición y que todas sus copias —el recurso hermano, los tres
//!    esquemas y la interfaz— coincidan con ella. Una copia decimal silenciosa
//!    en un solo sitio es lo que sobrevivió al cambio de contrato la vez
//!    anterior, y por eso el juez compara **cada copia por separado** en vez de
//!    compararlas entre sí.

use std::collections::BTreeSet;
use std::path::PathBuf;

use genome_population_simulator::pcr_amplification::{
    PCR_AMPLIFICATION_ALGORITHM, PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS,
    PCR_DUPLICATE_FRACTION_INVERSION_METHOD,
};
use genome_population_simulator::technical_profile::{
    F4TechnicalProfileError, GPS_HISEQX_PCR_PE150_RESOURCE, GPS_HISEQX_PCR_PE150_SHA256,
    GPS_HISEQX_PCR_PE150_V1, GPS_HISEQX_PCRFREE_PE150_RESOURCE, GPS_HISEQX_PCRFREE_PE150_SHA256,
    GPS_HISEQX_PCRFREE_PE150_V1, LibraryPreparation, PCR_MAXIMUM_CYCLES, PCR_MAXIMUM_EARLY_CYCLES,
    PCR_MAXIMUM_EFFICIENCY_CLASSES, PCR_MAXIMUM_FAMILY_SIZE, PCR_MAXIMUM_GC_ATTEMPTS,
    PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE, resolve_f4_technical_profile,
    resolve_f4_technical_profile_for, verify_initial_profile_resource, verify_pcr_profile_resource,
};
use serde_json::Value;

fn repository() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
}

fn read_json(relative: &str) -> Value {
    let bytes = std::fs::read(repository().join(relative))
        .unwrap_or_else(|error| panic!("no se puede leer {relative}: {error}"));
    serde_json::from_slice(&bytes).unwrap_or_else(|error| panic!("{relative} no es JSON: {error}"))
}

/// Recorre un camino de claves y devuelve el valor, o `None` si falta.
fn at<'a>(document: &'a Value, path: &[&str]) -> Option<&'a Value> {
    let mut current = document;
    for step in path {
        current = match step.parse::<usize>() {
            Ok(index) => current.get(index)?,
            Err(_) => current.get(*step)?,
        };
    }
    Some(current)
}

// ---------------------------------------------------------------------------
// F4.9.1 — el recurso vigente no se toca; el hermano existe.
// ---------------------------------------------------------------------------

/// **El recurso del perfil sin PCR conserva su huella y sus diez campos.**
///
/// Esta es la prueba que falla contra una edición del recurso vigente para
/// meter las dos preparaciones en uno: añadirle una clave `preparation` o un
/// bloque `amplification` cambiaría su SHA-256 y su lista de campos, y con
/// ellos los valores fijos de los tres esquemas publicados.
#[test]
fn the_pcr_free_resource_keeps_its_frozen_digest_and_fields() {
    let bytes = std::fs::read(repository().join(GPS_HISEQX_PCRFREE_PE150_RESOURCE))
        .expect("el recurso vigente existe");
    let digest = format!("{:x}", <sha2::Sha256 as sha2::Digest>::digest(&bytes));
    assert_eq!(
        digest, GPS_HISEQX_PCRFREE_PE150_SHA256,
        "el recurso sin PCR ha cambiado: su huella está congelada"
    );
    let document: Value = serde_json::from_slice(&bytes).expect("JSON");
    let fields = document
        .as_object()
        .expect("objeto")
        .keys()
        .cloned()
        .collect::<BTreeSet<_>>();
    let expected = [
        "schema_version",
        "profile_id",
        "classification",
        "description",
        "external_reference",
        "read",
        "insert",
        "quality",
        "error",
        "reproducibility",
    ]
    .into_iter()
    .map(str::to_owned)
    .collect::<BTreeSet<_>>();
    assert_eq!(
        fields, expected,
        "el recurso sin PCR no puede ganar ni perder campos"
    );
    assert!(
        !fields.contains("preparation") && !fields.contains("amplification"),
        "la preparación amplificada va en un recurso HERMANO, no dentro del vigente"
    );
    // Y sigue validándose por la vía de siempre.
    assert_eq!(verify_initial_profile_resource(), Ok(()));
}

/// **El perfil hermano se carga y declara su preparación.**
#[test]
fn the_sibling_profile_loads_and_declares_its_preparation() {
    assert_eq!(verify_pcr_profile_resource(), Ok(()));
    let profile = resolve_f4_technical_profile(GPS_HISEQX_PCR_PE150_V1).expect("perfil hermano");
    assert_eq!(profile.id(), GPS_HISEQX_PCR_PE150_V1);
    assert_eq!(profile.resource_path(), GPS_HISEQX_PCR_PE150_RESOURCE);
    assert_eq!(profile.resource_sha256(), GPS_HISEQX_PCR_PE150_SHA256);
    assert_eq!(profile.preparation(), LibraryPreparation::Pcr);
    assert_eq!(profile.preparation().as_str(), "pcr");

    let vigente = resolve_f4_technical_profile(GPS_HISEQX_PCRFREE_PE150_V1).expect("perfil");
    assert_eq!(vigente.preparation(), LibraryPreparation::PcrFree);
    assert_eq!(vigente.preparation().as_str(), "pcr_free");
    // **Idénticos en lectura, inserto y calidad**, que es lo que el contrato
    // exige del hermano. Si alguien tocara uno de esos valores en el recurso
    // nuevo, su validador lo rechazaría y esta comparación no llegaría a
    // ejecutarse; aquí se comprueba además que los tipos resueltos coinciden.
    assert_eq!(vigente.insert_length(), profile.insert_length());
    assert_eq!(vigente.common_read_length(), profile.common_read_length());
    assert_eq!(vigente.layout(), profile.layout());
    assert_eq!(vigente.orientation(), profile.orientation());
    assert_eq!(
        vigente.quality_curve(genome_population_simulator::technical_profile::ProfileMate::R1),
        profile.quality_curve(genome_population_simulator::technical_profile::ProfileMate::R1)
    );
    assert_eq!(
        vigente.quality_curve(genome_population_simulator::technical_profile::ProfileMate::R2),
        profile.quality_curve(genome_population_simulator::technical_profile::ProfileMate::R2)
    );
    // Y **distintos** en lo que los distingue, de modo que la igualdad de
    // arriba no pasa por ser dos nombres de la misma cosa.
    assert_ne!(vigente.id(), profile.id());
    assert_ne!(vigente.resource_sha256(), profile.resource_sha256());
    assert_ne!(vigente.preparation(), profile.preparation());
}

/// **La vía rápida de los hechos del perfil no puede separarse de la verificada.**
///
/// El documento efectivo y el manifiesto publican el identificador, la ruta, la
/// huella, la preparación y el inserto **sin volver a digerir el recurso**: es
/// trabajo gratis en cada escritura y, peor, entraba en pánico si el recurso
/// estuviera alterado en vez de emitir un rechazo. Esta prueba ata las dos
/// vías: lo que la rápida devuelve es exactamente lo que la verificada
/// comprueba.
#[test]
fn the_fast_profile_facts_agree_with_the_verified_profile() {
    for preparation in [LibraryPreparation::PcrFree, LibraryPreparation::Pcr] {
        let facts = genome_population_simulator::technical_profile::f4_profile_facts(preparation);
        let verified = resolve_f4_technical_profile_for(preparation).expect("perfil");
        assert_eq!(facts.id(), verified.id());
        assert_eq!(facts.resource_path(), verified.resource_path());
        assert_eq!(facts.sha256(), verified.resource_sha256());
        assert_eq!(facts.preparation(), verified.preparation());
        assert_eq!(facts.insert_length(), verified.insert_length());
        // Y la vía por identificador devuelve lo mismo.
        assert_eq!(
            genome_population_simulator::technical_profile::f4_profile_facts_by_id(facts.id()),
            Some(facts)
        );
    }
    // Control negativo con su control positivo: un identificador inventado no
    // devuelve hechos, y los dos que sí existen acaban de devolverlos.
    assert_eq!(
        genome_population_simulator::technical_profile::f4_profile_facts_by_id("gps-inventado-v1"),
        None
    );
}

/// La preparación elige su perfil, y un identificador inventado se rechaza.
#[test]
fn a_preparation_resolves_its_own_profile_and_an_unknown_id_is_rejected() {
    assert_eq!(
        resolve_f4_technical_profile_for(LibraryPreparation::PcrFree)
            .expect("perfil")
            .id(),
        GPS_HISEQX_PCRFREE_PE150_V1
    );
    assert_eq!(
        resolve_f4_technical_profile_for(LibraryPreparation::Pcr)
            .expect("perfil")
            .id(),
        GPS_HISEQX_PCR_PE150_V1
    );
    assert_eq!(
        LibraryPreparation::parse("pcr"),
        Some(LibraryPreparation::Pcr)
    );
    assert_eq!(
        LibraryPreparation::parse("pcr_free"),
        Some(LibraryPreparation::PcrFree)
    );
    assert_eq!(LibraryPreparation::parse("PCR"), None);
    assert!(matches!(
        resolve_f4_technical_profile("gps-inventado-v1"),
        Err(F4TechnicalProfileError::UnknownProfile { .. })
    ));
}

// ---------------------------------------------------------------------------
// F4.9.2 — el juez de las constantes de contrato.
// ---------------------------------------------------------------------------

/// Dónde vive cada copia de cada constante, una entrada por copia.
///
/// Se escriben **por separado** a propósito: comparar «todas las copias entre
/// sí» dejaría pasar un cambio coherente en dos sitios, y comparar solo una
/// dejaría pasar el resto.
fn constant_sites() -> Vec<(&'static str, &'static str, Vec<&'static str>, u64)> {
    let config = "schemas/config-general-v4.schema.json";
    let effective = "schemas/effective-config-v6.schema.json";
    let manifest = "schemas/manifest-event-v3.schema.json";
    let resource = "resources/f4/gps-hiseqx-pcr-pe150-v1.json";
    let config_pcr = [
        "$defs",
        "technical_event",
        "properties",
        "pcr",
        "properties",
    ];
    let effective_pcr = [
        "$defs",
        "technical_artifacts",
        "properties",
        "pcr",
        "oneOf",
        "1",
        "properties",
    ];
    let manifest_pcr = [
        "$defs",
        "technical_artifacts",
        "properties",
        "pcr",
        "oneOf",
        "1",
        "properties",
    ];
    let mut sites = Vec::new();
    let scalars: [(&str, &str, u64); 4] = [
        ("cycles", "maximum_cycles", u64::from(PCR_MAXIMUM_CYCLES)),
        (
            "early_cycles",
            "maximum_early_cycles",
            u64::from(PCR_MAXIMUM_EARLY_CYCLES),
        ),
        (
            "max_family_size",
            "maximum_family_size",
            u64::from(PCR_MAXIMUM_FAMILY_SIZE),
        ),
        (
            "gc_max_attempts",
            "maximum_gc_attempts",
            u64::from(PCR_MAXIMUM_GC_ATTEMPTS),
        ),
    ];
    for (field, resource_field, value) in scalars {
        for (file, prefix) in [
            (config, config_pcr.as_slice()),
            (effective, effective_pcr.as_slice()),
            (manifest, manifest_pcr.as_slice()),
        ] {
            let mut path = prefix.to_vec();
            path.push(field);
            path.push("maximum");
            sites.push((file, field, path, value));
        }
        sites.push((
            resource,
            resource_field,
            vec!["amplification", resource_field],
            value,
        ));
    }
    // El tope de clases es un `maxItems` en los esquemas, no un `maximum`.
    for (file, prefix) in [
        (config, config_pcr.as_slice()),
        (effective, effective_pcr.as_slice()),
        (manifest, manifest_pcr.as_slice()),
    ] {
        let mut path = prefix.to_vec();
        path.push("efficiency_classes");
        path.push("maxItems");
        sites.push((
            file,
            "efficiency_classes",
            path,
            u64::from(PCR_MAXIMUM_EFFICIENCY_CLASSES),
        ));
    }
    sites.push((
        resource,
        "maximum_efficiency_classes",
        vec!["amplification", "maximum_efficiency_classes"],
        u64::from(PCR_MAXIMUM_EFFICIENCY_CLASSES),
    ));
    // Y el tope de parejas por molécula aparece en los dos documentos que lo
    // publican.
    for (file, prefix) in [
        (effective, effective_pcr.as_slice()),
        (manifest, manifest_pcr.as_slice()),
    ] {
        let mut path = prefix.to_vec();
        path.push("mean_pairs_per_molecule");
        path.push("maximum");
        sites.push((
            file,
            "mean_pairs_per_molecule",
            path,
            PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE as u64,
        ));
    }
    sites.push((
        resource,
        "maximum_mean_pairs_per_molecule",
        vec!["amplification", "maximum_mean_pairs_per_molecule"],
        PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE as u64,
    ));
    sites
}

/// **Cada copia de cada constante coincide con su constante de Rust.**
///
/// # El control positivo, corregido
///
/// La primera redacción decía que sin él «un camino mal escrito devolvería
/// `None` en silencio», y es **falso**: un camino mal escrito revienta con
/// `panic!` en el `unwrap_or_else` de abajo. Y la comparación entre el
/// contador y el número de sitios estaba **implicada por el bucle** que la
/// precede, de modo que no podía fallar nunca.
///
/// El control verdadero es el que compara con **23**, escrito a mano: si
/// alguien retirase una fila de la lista de sitios, el bucle seguiría cuadrando
/// consigo mismo y solo esa cifra lo notaría.
#[test]
fn every_schema_and_resource_maximum_matches_its_rust_constant() {
    let sites = constant_sites();
    let mut compared = 0_usize;
    let mut files = BTreeSet::new();
    for (file, field, path, expected) in &sites {
        let document = read_json(file);
        let found = at(&document, path).unwrap_or_else(|| {
            panic!("{file} no tiene el camino {path:?} para la constante de {field}")
        });
        let found = found
            .as_u64()
            .unwrap_or_else(|| panic!("{file}:{path:?} no es un entero: {found}"));
        assert_eq!(
            found, *expected,
            "{file} declara {found} para {field} y la constante de Rust vale {expected}"
        );
        compared += 1;
        files.insert(*file);
    }
    // 4 escalares en los tres esquemas y en el recurso (16), el tope de
    // clases en los tres esquemas y en el recurso (4), y el tope de parejas
    // por molécula en los dos documentos que lo publican y en el recurso (3).
    // El usuario no declara `m`, de modo que el esquema de configuración no
    // la acota: son 23 y no 24.
    assert_eq!(compared, 23, "se esperaban 23 copias de constante");
    assert_eq!(files.len(), 4, "las copias viven en cuatro ficheros");
}

/// **La interfaz repite las mismas constantes.**
///
/// Se leen del texto de `ui/config_editor.py`, que es donde el cuaderno acota
/// sus controles. El control positivo comprueba que la búsqueda encontraría
/// una cifra distinta, de modo que una lectura rota no pasa por coincidencia.
#[test]
fn the_notebook_controls_repeat_the_same_contract_constants() {
    let source = std::fs::read_to_string(repository().join("ui/config_editor.py"))
        .expect("el editor de la interfaz existe");
    let expected: [(&str, String); 6] = [
        ("PCR_MAXIMUM_CYCLES", PCR_MAXIMUM_CYCLES.to_string()),
        (
            "PCR_MAXIMUM_EARLY_CYCLES",
            PCR_MAXIMUM_EARLY_CYCLES.to_string(),
        ),
        (
            "PCR_MAXIMUM_FAMILY_SIZE",
            PCR_MAXIMUM_FAMILY_SIZE.to_string(),
        ),
        (
            "PCR_MAXIMUM_GC_ATTEMPTS",
            PCR_MAXIMUM_GC_ATTEMPTS.to_string(),
        ),
        (
            "PCR_MAXIMUM_EFFICIENCY_CLASSES",
            PCR_MAXIMUM_EFFICIENCY_CLASSES.to_string(),
        ),
        (
            "PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE",
            format!("{PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE:?}"),
        ),
    ];
    let mut compared = 0_usize;
    for (name, value) in &expected {
        let assignment = format!("{name} = {value}\n");
        assert!(
            source.contains(&assignment),
            "ui/config_editor.py no declara `{name} = {value}`"
        );
        compared += 1;
        // Control positivo de la búsqueda: una cifra distinta **no** aparece.
        let wrong = format!("{name} = 99999\n");
        assert!(
            !source.contains(&wrong),
            "la búsqueda encuentra cualquier cosa: está rota"
        );
    }
    assert_eq!(compared, 6, "la interfaz repite seis constantes");
}

/// El nombre del algoritmo y el del método de inversión también son contrato.
#[test]
fn the_algorithm_and_inversion_names_match_across_resource_and_schemas() {
    let resource = read_json("resources/f4/gps-hiseqx-pcr-pe150-v1.json");
    assert_eq!(
        at(&resource, &["amplification", "algorithm"]).and_then(Value::as_str),
        Some(PCR_AMPLIFICATION_ALGORITHM)
    );
    assert_eq!(
        at(
            &resource,
            &["amplification", "duplicate_fraction_inversion"]
        )
        .and_then(Value::as_str),
        Some(PCR_DUPLICATE_FRACTION_INVERSION_METHOD)
    );
    assert_eq!(
        at(
            &resource,
            &["amplification", "duplicate_fraction_inversion_iterations"]
        )
        .and_then(Value::as_u64),
        Some(u64::from(PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS))
    );
    for (file, prefix) in [
        (
            "schemas/effective-config-v6.schema.json",
            [
                "$defs",
                "technical_artifacts",
                "properties",
                "pcr",
                "oneOf",
                "1",
                "properties",
            ],
        ),
        (
            "schemas/manifest-event-v3.schema.json",
            [
                "$defs",
                "technical_artifacts",
                "properties",
                "pcr",
                "oneOf",
                "1",
                "properties",
            ],
        ),
    ] {
        let document = read_json(file);
        let mut algorithm = prefix.to_vec();
        algorithm.extend(["algorithm", "const"]);
        assert_eq!(
            at(&document, &algorithm).and_then(Value::as_str),
            Some(PCR_AMPLIFICATION_ALGORITHM),
            "{file} no fija el algoritmo"
        );
        let mut method = prefix.to_vec();
        method.extend([
            "duplicate_fraction_inversion",
            "properties",
            "method",
            "const",
        ]);
        assert_eq!(
            at(&document, &method).and_then(Value::as_str),
            Some(PCR_DUPLICATE_FRACTION_INVERSION_METHOD),
            "{file} no fija el método de inversión"
        );
    }
}

/// **La tabla nueva de cada esquema está cerrada y con dominios acotados.**
///
/// Falla contra `additionalProperties: true` y contra un dominio abierto en
/// cualquiera de los parámetros.
#[test]
fn the_new_schema_tables_are_closed_and_bounded() {
    let config = read_json("schemas/config-general-v4.schema.json");
    let table = at(&config, &["$defs", "technical_event", "properties", "pcr"])
        .expect("la tabla pcr existe");
    assert_eq!(
        table.get("additionalProperties"),
        Some(&Value::Bool(false)),
        "la tabla de amplificación no puede admitir claves desconocidas"
    );
    let properties = table
        .get("properties")
        .and_then(Value::as_object)
        .expect("propiedades");
    let mut bounded = 0_usize;
    for (name, schema) in properties {
        match name.as_str() {
            "gc_relative_efficiency" | "efficiency_classes" => {
                let items = schema.get("items").expect("items");
                assert_eq!(
                    items.get("additionalProperties"),
                    Some(&Value::Bool(false)),
                    "las filas de {name} no pueden admitir claves desconocidas"
                );
                for (field, row) in items
                    .get("properties")
                    .and_then(Value::as_object)
                    .expect("filas")
                {
                    let has_low =
                        row.get("minimum").is_some() || row.get("exclusiveMinimum").is_some();
                    let has_high =
                        row.get("maximum").is_some() || row.get("exclusiveMaximum").is_some();
                    assert!(
                        has_low && has_high,
                        "{name}.{field} tiene un dominio abierto"
                    );
                    bounded += 1;
                }
            }
            _ => {
                let has_low =
                    schema.get("minimum").is_some() || schema.get("exclusiveMinimum").is_some();
                let has_high =
                    schema.get("maximum").is_some() || schema.get("exclusiveMaximum").is_some();
                assert!(has_low && has_high, "{name} tiene un dominio abierto");
                bounded += 1;
            }
        }
    }
    assert_eq!(
        bounded,
        10,
        "se esperaban diez dominios acotados: {:?}",
        properties.keys().collect::<Vec<_>>()
    );
    // **Ningún valor por omisión materializado para la fracción de
    // duplicados**: escribirlo afirmaría una calibración que nadie ha hecho.
    assert!(
        properties["duplicate_fraction"].get("default").is_none(),
        "la fracción de duplicados no puede tener valor por omisión"
    );
    assert!(
        properties["cycles"].get("default").is_none(),
        "los ciclos no pueden tener valor por omisión"
    );
    assert!(
        properties["efficiency"].get("default").is_none(),
        "la eficiencia no puede tener valor por omisión"
    );
    // Control positivo del mismo hueco: un campo del mismo esquema que **sí**
    // tiene valor por omisión, de modo que la aserción no pasa por mirar un
    // sitio donde nunca hay `default`.
    let read_length = at(
        &config,
        &[
            "$defs",
            "technical_event",
            "properties",
            "read_length_bases",
            "default",
        ],
    );
    assert!(
        read_length.is_some(),
        "el control positivo espera que read_length_bases sí tenga valor por omisión"
    );
}
