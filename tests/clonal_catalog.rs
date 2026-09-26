//! F4.10.1: el catálogo de ámbitos aleatorios `f4-technical-v6`.
//!
//! # Qué protege este fichero
//!
//! Que **añadir el tartamudeo clonal no desplace ninguna decisión anterior**.
//! El nombre del catálogo es además el ámbito de derivación: si las cuatro
//! claves del linaje nacieran en el v5, compartirían ámbito con las cinco de
//! la amplificación y con las tres de `fragment` con componente `attempt`; y
//! si las heredadas subieran al v6, todas sus subsemillas cambiarían y con
//! ellas todos los FASTQ publicados.
//!
//! # Lo que este fichero **no** ancla, y dónde está eso
//!
//! Aquí se congela lo que **cada constructora** produce, no **cuál elige el
//! muestreador**. Una edición que hiciera que el muestreador pidiera la clave
//! equivocada dejaría estos vectores intactos. Quien ata eso es
//! `tests/clonal_lineage.rs`, que compara ejecuciones completas, y
//! `tests/pcr_frozen_identity.rs`, con la huella congelada de un escenario que
//! no amplifica.

use genome_population_simulator::identity::{EventId, LibraryId, RepeatLocusId};
use genome_population_simulator::randomness::{
    F4_TECHNICAL_RNG_CATALOG_V3, F4_TECHNICAL_RNG_CATALOG_V4, F4_TECHNICAL_RNG_CATALOG_V5,
    F4_TECHNICAL_RNG_CATALOG_V6, F4TechnicalDecision, FragmentRetryDecision, RootSeed,
    derive_f4_technical_subseed,
};

fn root() -> RootSeed {
    RootSeed::parse(
        "root",
        "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f",
    )
    .expect("semilla")
}

fn identifiers() -> (EventId, LibraryId, RepeatLocusId) {
    (
        EventId::parse("reads-main".to_owned(), "test").expect("evento"),
        LibraryId::parse("LIB_0001".to_owned(), "test").expect("biblioteca"),
        RepeatLocusId::parse("STR-1".to_owned(), "test").expect("locus"),
    )
}

fn key(decision: F4TechnicalDecision<'_>) -> String {
    derive_f4_technical_subseed(root(), decision)
        .expect("subsemilla")
        .to_hex()
}

/// **Las cuatro decisiones del linaje nacen en el v6; nada más se mueve.**
#[test]
fn only_the_lineage_decisions_belong_to_the_new_catalog() {
    let (event_id, library_id, locus_id) = identifiers();
    assert_eq!(F4_TECHNICAL_RNG_CATALOG_V6, "f4-technical-v6");
    for lineage in [
        F4TechnicalDecision::lineage_parent(&event_id, &library_id, 17, 3),
        F4TechnicalDecision::lineage_divergence(&event_id, &library_id, 17, 3),
        F4TechnicalDecision::lineage_slip(&event_id, &library_id, 17, 3, 5, &locus_id),
        F4TechnicalDecision::lineage_step(&event_id, &library_id, 17, 3, 5, &locus_id),
    ] {
        assert_eq!(lineage.catalog(), F4_TECHNICAL_RNG_CATALOG_V6);
    }
    let attempt = std::num::NonZeroU64::new(1).expect("intento");
    for amplification in [
        F4TechnicalDecision::amplification_efficiency_class(&event_id, &library_id, 17),
        F4TechnicalDecision::amplification_early_copies(&event_id, &library_id, 17),
        F4TechnicalDecision::amplification_gc_accept(&event_id, &library_id, 17, 0),
        F4TechnicalDecision::amplification_family_size(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_retry(
            &event_id,
            &library_id,
            17,
            FragmentRetryDecision::InsertLength,
            attempt,
        ),
    ] {
        assert_eq!(amplification.catalog(), F4_TECHNICAL_RNG_CATALOG_V5);
    }
    for stutter in [
        F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, 17, &locus_id),
        F4TechnicalDecision::stutter_step(&event_id, &library_id, 17, &locus_id),
    ] {
        assert_eq!(stutter.catalog(), F4_TECHNICAL_RNG_CATALOG_V4);
    }
    for inherited in [
        F4TechnicalDecision::fragment_insert_length(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_copy(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_start(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_bernoulli(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_breakpoint(&event_id, &library_id, 17),
    ] {
        assert_eq!(inherited.catalog(), F4_TECHNICAL_RNG_CATALOG_V3);
    }
}

/// **Vectores congelados de las cuatro claves nuevas.**
///
/// Están escritos a mano —se leyeron una vez y se fijaron aquí—, no calculados
/// por el mismo código que comprueban. Cualquier cambio del orden o del
/// literal de los componentes los mueve: el nombre de la etapa, el del
/// catálogo, el `copy_index`, el `cycle` o el `locus_id`.
#[test]
fn the_lineage_stage_has_literal_subseed_vectors() {
    let (event_id, library_id, locus_id) = identifiers();
    assert_eq!(
        key(F4TechnicalDecision::lineage_parent(
            &event_id,
            &library_id,
            17,
            3
        )),
        "3d54830cbf9405777a26967cc4d38c996f8c2eb62dee1b1a674bd373d24a4caf"
    );
    assert_eq!(
        key(F4TechnicalDecision::lineage_divergence(
            &event_id,
            &library_id,
            17,
            3
        )),
        "0acbb46f7463998bf95defba0ed327aa43688553aa01f6070688ed926e75e2bf"
    );
    assert_eq!(
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            17,
            3,
            5,
            &locus_id
        )),
        "1f742d17cd9fc622a02e8e8abfabc8a7788e36d67c55e9f64a8b48263a029aac"
    );
    assert_eq!(
        key(F4TechnicalDecision::lineage_step(
            &event_id,
            &library_id,
            17,
            3,
            5,
            &locus_id
        )),
        "fffdf2cbe1a67e3e4845042b1f16e014c5d65dbb28741be5609b71b0d9052949"
    );
}

/// **Las cuatro claves son distintas entre sí y de las de la amplificación.**
///
/// Sin esto, una implementación que reutilizara la clave del tamaño de familia
/// para el progenitor ligaría dos sorteos que tienen que ser independientes, y
/// ninguna otra prueba lo vería.
#[test]
fn the_four_lineage_decisions_do_not_share_a_key() {
    let (event_id, library_id, locus_id) = identifiers();
    let keys = [
        key(F4TechnicalDecision::lineage_parent(
            &event_id,
            &library_id,
            17,
            3,
        )),
        key(F4TechnicalDecision::lineage_divergence(
            &event_id,
            &library_id,
            17,
            3,
        )),
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            17,
            3,
            5,
            &locus_id,
        )),
        key(F4TechnicalDecision::lineage_step(
            &event_id,
            &library_id,
            17,
            3,
            5,
            &locus_id,
        )),
        key(F4TechnicalDecision::amplification_family_size(
            &event_id,
            &library_id,
            17,
        )),
        key(F4TechnicalDecision::amplification_early_copies(
            &event_id,
            &library_id,
            17,
        )),
    ];
    let distinct = keys.iter().collect::<std::collections::BTreeSet<_>>();
    assert_eq!(distinct.len(), keys.len(), "dos decisiones comparten clave");
}

/// **Cada componente del ámbito separa de verdad.**
///
/// Cambiar la familia, la copia, el ciclo o el locus tiene que cambiar la
/// clave. Sin esta prueba, una implementación que olvidara meter el ciclo en
/// el material daría el mismo resbalón en los diez ciclos y el tartamudeo
/// clonal se volvería un sorteo por copia, que es justo lo contrario del
/// fenómeno.
#[test]
fn every_component_of_the_lineage_scope_separates() {
    let (event_id, library_id, locus_id) = identifiers();
    let other_locus = RepeatLocusId::parse("STR-2".to_owned(), "test").expect("locus");
    let base = key(F4TechnicalDecision::lineage_slip(
        &event_id,
        &library_id,
        17,
        3,
        5,
        &locus_id,
    ));
    assert_ne!(
        base,
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            18,
            3,
            5,
            &locus_id
        )),
        "la familia separa"
    );
    assert_ne!(
        base,
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            17,
            4,
            5,
            &locus_id
        )),
        "la copia separa"
    );
    assert_ne!(
        base,
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            17,
            3,
            6,
            &locus_id
        )),
        "el ciclo separa"
    );
    assert_ne!(
        base,
        key(F4TechnicalDecision::lineage_slip(
            &event_id,
            &library_id,
            17,
            3,
            5,
            &other_locus
        )),
        "el locus separa"
    );
}
