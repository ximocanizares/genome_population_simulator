//! F1.7.2: inserción de tractos en una referencia real, referencia derivada con
//! identidad propia y tabla de correspondencia de coordenadas.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::reference_tract_insertion::{
    CoordinateSegmentClass, DeclaredInsertion, coordinate_map_document, insert_tracts,
};

static NONCE: AtomicU64 = AtomicU64::new(0);

/// Genoma de partida determinista: cada base depende solo de su posición, de
/// modo que la prueba puede comprobar identidad base a base sin oráculo aparte.
fn supplied_reference(id: &str, length: u64) -> Reference {
    let bases: Vec<u8> = (0..length)
        .map(|position| b"ACGT"[(position % 4) as usize])
        .collect();
    Reference::new(
        ReferenceId::parse(id, "execution").expect("referencia"),
        ReferenceOrigin::Supplied {
            original_path: PathBuf::from("/dev/null"),
        },
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("secuencia"),
                FastaHeader::parse(b"chr1").expect("cabecera"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(&bases).expect("bases"),
            )
            .expect("secuencia de referencia"),
        ],
    )
    .expect("referencia suministrada")
}

/// Genoma de partida con **dos** secuencias, cada una con su propia base
/// característica, de modo que una fuga entre ellas sea visible a simple vista.
fn two_sequence_reference(id: &str) -> Reference {
    let sequence = |name: &str, base: u8, length: usize| {
        ReferenceSequence::new(
            SequenceId::parse(name, "test").expect("secuencia"),
            FastaHeader::parse(name.as_bytes()).expect("cabecera"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(&vec![base; length]).expect("bases"),
        )
        .expect("secuencia de referencia")
    };
    Reference::new(
        ReferenceId::parse(id, "execution").expect("referencia"),
        ReferenceOrigin::Supplied {
            original_path: PathBuf::from("/dev/null"),
        },
        vec![sequence("chr1", b'A', 200), sequence("chr2", b'T', 200)],
    )
    .expect("referencia suministrada")
}

fn bases_of(reference: &Reference) -> Vec<u8> {
    reference.sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases materializadas")
        .to_vec()
}

/// Una inserción de treinta bases en la posición 1000 desplaza la 2000 a la
/// 2030, y la base que hay allí es la misma que había en la 2000 del original.
///
/// *Cuidado:* con **una sola** inserción, una implementación que leyese las
/// posiciones declaradas en coordenadas del genoma ya modificado daría
/// exactamente el mismo resultado, porque no hay ninguna anterior que la haya
/// desplazado. Esta prueba **no discrimina** ese defecto; lo hace la de dos
/// inserciones seguidas.
#[test]
fn one_insertion_shifts_the_map_and_the_bases_follow() {
    let source = supplied_reference("ref-supplied", 4000);
    let original = bases_of(&source);
    let tract = b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG";
    assert_eq!(tract.len(), 30);
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr1",
            source_offset: 1000,
            bases: tract,
        }],
    )
    .expect("inserción de un tracto");

    assert_eq!(
        derived.map().translate_source_position("chr1", 2000),
        Some(2030)
    );
    let result = bases_of(derived.reference());
    assert_eq!(result.len(), 4030);
    assert_eq!(result[2030], original[2000]);
    assert_eq!(&result[1000..1030], tract);
    assert_eq!(&result[..1000], &original[..1000]);
    assert_eq!(&result[1030..], &original[1000..]);
}

/// **La prueba que discrimina.** Con dos inserciones seguidas en la misma
/// secuencia, la segunda posición declarada sigue refiriéndose al original y no
/// al intermedio. Una implementación que la interpretase sobre el genoma ya
/// modificado colocaría el segundo tracto treinta bases antes.
#[test]
fn the_second_declared_position_still_refers_to_the_original_genome() {
    let source = supplied_reference("ref-supplied", 4000);
    let original = bases_of(&source);
    let first = b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG";
    let second = b"GATAGATAGATAGATAGATA";
    assert_eq!(first.len(), 30);
    assert_eq!(second.len(), 20);
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[
            DeclaredInsertion {
                sequence_id: "chr1",
                source_offset: 1000,
                bases: first,
            },
            DeclaredInsertion {
                sequence_id: "chr1",
                source_offset: 2000,
                bases: second,
            },
        ],
    )
    .expect("dos inserciones seguidas");

    let result = bases_of(derived.reference());
    assert_eq!(result.len(), 4050);
    assert_eq!(&result[1000..1030], first);
    assert_eq!(
        &result[2030..2050],
        second,
        "el segundo tracto entra donde estaba la base 2000 del original, ya desplazada"
    );
    assert_eq!(
        derived.map().translate_source_position("chr1", 2000),
        Some(2050)
    );
    assert_eq!(result[2050], original[2000]);
    assert_eq!(
        derived.map().translate_source_position("chr1", 1999),
        Some(2029)
    );
    assert_eq!(result[2029], original[1999]);
}

/// La referencia derivada tiene identificador propio y declara de qué genoma
/// procede. Publicarla con el identificador del FASTA suministrado atribuiría
/// los resultados a un genoma que no existe.
#[test]
fn the_derived_reference_carries_its_own_identity() {
    let source = supplied_reference("ref-supplied", 500);
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr1",
            source_offset: 100,
            bases: b"CAGCAGCAG",
        }],
    )
    .expect("inserción");
    assert_eq!(derived.reference().id().as_str(), "ref-derived");
    assert_ne!(derived.reference().id(), source.id());
    match derived.reference().origin() {
        ReferenceOrigin::TractInsertion {
            source_reference_id,
        } => assert_eq!(source_reference_id, source.id()),
        other => panic!("origen inesperado: {other:?}"),
    }
    assert_eq!(derived.map().derived_reference_id().as_str(), "ref-derived");
    assert_eq!(derived.map().source_reference_id().as_str(), "ref-supplied");

    let reused = insert_tracts(
        &source,
        "ref-supplied",
        &[DeclaredInsertion {
            sequence_id: "chr1",
            source_offset: 100,
            bases: b"CAGCAGCAG",
        }],
    )
    .expect_err("reutilizar el identificador del FASTA suministrado");
    assert!(format!("{reused:?}").contains("must not reuse"));
}

/// El documento del mapa valida contra su esquema con el validador externo, y
/// el segmento insertado no lleva intervalo de origen mientras los conservados
/// sí lo llevan.
#[test]
fn the_published_map_validates_against_its_schema() {
    let source = supplied_reference("ref-supplied", 400);
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[
            DeclaredInsertion {
                sequence_id: "chr1",
                source_offset: 100,
                bases: b"CAGCAGCAG",
            },
            DeclaredInsertion {
                sequence_id: "chr1",
                source_offset: 200,
                bases: b"GATAGATA",
            },
        ],
    )
    .expect("inserciones");

    let classes: Vec<CoordinateSegmentClass> = derived
        .map()
        .segments()
        .iter()
        .map(|segment| segment.class())
        .collect();
    assert_eq!(
        classes,
        vec![
            CoordinateSegmentClass::Retained,
            CoordinateSegmentClass::Inserted,
            CoordinateSegmentClass::Retained,
            CoordinateSegmentClass::Inserted,
            CoordinateSegmentClass::Retained,
        ]
    );
    for segment in derived.map().segments() {
        match segment.class() {
            CoordinateSegmentClass::Retained => assert!(segment.source_interval().is_some()),
            CoordinateSegmentClass::Inserted => assert!(segment.source_interval().is_none()),
        }
    }

    let document = coordinate_map_document(derived.map()).expect("documento del mapa");
    let path = temporary("json");
    fs::write(
        &path,
        serde_json::to_vec_pretty(&document).expect("serialización"),
    )
    .expect("documento");
    assert_external_schema(&path, "reference-coordinate-map-v2.schema.json");
    fs::remove_file(&path).ok();
}

/// El esquema del mapa de **reducción** no cambia ni un byte.
///
/// La regla de la sección 10.3 del plan: los esquemas publicados no se editan;
/// la capacidad nueva va en una versión nueva.
#[test]
fn the_existing_reduction_map_schema_is_untouched() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let path = root.join("schemas/reference-reduction-map-v1.schema.json");
    let identity = genome_population_simulator::integrity::identify_file(
        &path,
        "schemas/reference-reduction-map-v1.schema.json",
    )
    .expect("huella del esquema");
    assert_eq!(
        identity.digest(),
        "30e3924a78915ce75b86b929976bac7ca07a7a72d2e09ae16a68c15382f4d027",
        "el esquema publicado del mapa de reducción no puede cambiar"
    );
}

/// Una inserción sobre una secuencia que no existe, o que rebasa el final, se
/// rechaza nombrando la secuencia.
#[test]
fn an_insertion_outside_the_genome_is_rejected() {
    let source = supplied_reference("ref-supplied", 100);
    let unknown = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr9",
            source_offset: 0,
            bases: b"CAG",
        }],
    )
    .expect_err("secuencia desconocida");
    assert!(format!("{unknown:?}").contains("chr9"));

    let past = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr1",
            source_offset: 101,
            bases: b"CAG",
        }],
    )
    .expect_err("punto fuera de la secuencia");
    assert!(format!("{past:?}").contains("runs past"));
}

/// **Las inserciones se filtran por secuencia.** Sin ese filtro, un tracto
/// declarado sobre el segundo cromosoma se escribiría también en el primero, y
/// el mapa de coordenadas desplazaría posiciones que no se han movido.
///
/// Ningún otro fixture de este fichero usa más de una secuencia, así que esta es
/// la única prueba que lo ejercita.
#[test]
fn an_insertion_declared_on_one_sequence_does_not_leak_into_another() {
    let source = two_sequence_reference("ref-supplied");
    let tract = b"CAGCAGCAGCAG";
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr2",
            source_offset: 100,
            bases: tract,
        }],
    )
    .expect("inserción en la segunda secuencia");

    let first = derived.reference().sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases del primero");
    let second = derived.reference().sequences()[1]
        .sequence()
        .as_bytes()
        .expect("bases del segundo");
    assert_eq!(
        first.len(),
        200,
        "el primer cromosoma no cambia de longitud"
    );
    assert!(
        first.iter().all(|base| *base == b'A'),
        "el primer cromosoma no puede contener ni una base del tracto del segundo"
    );
    assert_eq!(second.len(), 212);
    assert_eq!(&second[100..112], tract);

    // Y el mapa tampoco desplaza lo que no se ha movido.
    assert_eq!(
        derived.map().translate_source_position("chr1", 150),
        Some(150),
        "una posición del primer cromosoma no puede desplazarse"
    );
    assert_eq!(
        derived.map().translate_source_position("chr2", 150),
        Some(162)
    );
}

/// El esquema del mapa rechaza los documentos que violan sus condicionales: un
/// segmento conservado sin intervalo de origen, y uno insertado que lo lleva.
///
/// Son las dos reglas propias de la versión 2 y no las comprobaba nadie.
#[test]
fn the_coordinate_map_schema_rejects_the_documents_that_violate_its_conditionals() {
    let source = supplied_reference("ref-supplied", 200);
    let derived = insert_tracts(
        &source,
        "ref-derived",
        &[DeclaredInsertion {
            sequence_id: "chr1",
            source_offset: 100,
            bases: b"CAGCAGCAG",
        }],
    )
    .expect("inserción");
    let valid = coordinate_map_document(derived.map()).expect("documento del mapa");

    // Control positivo: el documento correcto se admite. Sin él, un esquema que
    // rechazase todo pasaría los dos casos negativos sin discriminar nada.
    assert!(validates(&valid, "reference-coordinate-map-v2.schema.json"));

    let mut retained_without_source = valid.clone();
    retained_without_source["segments"][0]
        .as_object_mut()
        .expect("segmento")
        .remove("source_interval");
    assert!(
        !validates(
            &retained_without_source,
            "reference-coordinate-map-v2.schema.json"
        ),
        "un segmento conservado sin intervalo de origen tiene que rechazarse"
    );

    let mut inserted_with_source = valid.clone();
    inserted_with_source["segments"][1]["source_interval"] =
        serde_json::json!({"start": 1, "end": 9});
    assert!(
        !validates(
            &inserted_with_source,
            "reference-coordinate-map-v2.schema.json"
        ),
        "un segmento insertado con intervalo de origen tiene que rechazarse"
    );

    let mut unknown_class = valid.clone();
    unknown_class["segments"][0]["segment_class"] = serde_json::json!("residual");
    assert!(
        !validates(&unknown_class, "reference-coordinate-map-v2.schema.json"),
        "la clase de segmento es un conjunto cerrado"
    );
}

/// Ejecuta el validador externo y devuelve si admitió el documento.
fn validates(document: &serde_json::Value, schema: &str) -> bool {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let path = temporary("json");
    fs::write(
        &path,
        serde_json::to_vec_pretty(document).expect("serialización"),
    )
    .expect("documento");
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(&path)
        .output()
        .expect("validador externo de JSON Schema");
    fs::remove_file(&path).ok();
    output.status.success()
}

fn assert_external_schema(path: &Path, schema: &str) {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let output = Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(path)
        .output()
        .expect("validador externo de JSON Schema");
    assert!(
        output.status.success(),
        "validador externo para {path:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
}

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f1-7-2-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}
