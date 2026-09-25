//! F4.9.8: el catálogo de ámbitos aleatorios `f4-technical-v5`.
//!
//! # Qué protege este fichero
//!
//! Que **añadir la amplificación no desplace ninguna decisión anterior**. El
//! nombre del catálogo es además el ámbito de derivación: si las claves de la
//! etapa `amplification` nacieran en el v3 o en el v4, compartirían ámbito con
//! decisiones ya publicadas; y si las heredadas subieran al v5, todas sus
//! subsemillas cambiarían y con ellas todos los FASTQ publicados.
//!
//! El otro contrato que se comprueba aquí es la separación entre el **intento
//! cero** y los reintentos: el intento cero usa literalmente las claves de hoy,
//! que es lo que hace que una biblioteca amplificada con curva neutra consuma
//! exactamente lo mismo que una sin amplificar.
//!
//! # Lo que este fichero **no** ancla, y dónde está eso
//!
//! Aquí se congela lo que **cada constructora** produce, no **cuál elige el
//! muestreador**. Una edición que hiciera que el muestreador pidiera la clave
//! equivocada dejaría estos vectores intactos. Quien ata eso es
//! `tests/pcr_frozen_identity.rs`, con la huella congelada del FASTQ publicado
//! y la de la primera pareja del muestreador; se comprobó que sin él esa
//! edición pasaba la batería entera.

use genome_population_simulator::identity::{EventId, LibraryId, RepeatLocusId};
use genome_population_simulator::randomness::{
    F4_TECHNICAL_RNG_CATALOG, F4_TECHNICAL_RNG_CATALOG_V3, F4_TECHNICAL_RNG_CATALOG_V4,
    F4_TECHNICAL_RNG_CATALOG_V5, F4TechnicalDecision, FragmentRetryDecision, RootSeed,
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

/// **La etapa nueva nace en el v5; las heredadas y el tartamudeo no se mueven.**
#[test]
fn only_the_amplification_stage_belongs_to_the_new_catalog() {
    let (event_id, library_id, locus_id) = identifiers();
    assert_eq!(F4_TECHNICAL_RNG_CATALOG, F4_TECHNICAL_RNG_CATALOG_V3);
    assert_eq!(F4_TECHNICAL_RNG_CATALOG_V5, "f4-technical-v5");
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
    for inherited in [
        F4TechnicalDecision::fragment_insert_length(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_copy(&event_id, &library_id, 17),
        F4TechnicalDecision::fragment_start(&event_id, &library_id, 17),
        F4TechnicalDecision::chimera_bernoulli(&event_id, &library_id, 17),
        F4TechnicalDecision::index_hop_bernoulli(
            &event_id,
            &library_id,
            genome_population_simulator::randomness::F4IndexHopInstance::parse("PH_ABCDEFGH")
                .expect("instancia"),
        ),
    ] {
        assert_eq!(inherited.catalog(), F4_TECHNICAL_RNG_CATALOG_V3);
    }
    for stutter in [
        F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, 17, &locus_id),
        F4TechnicalDecision::stutter_step(&event_id, &library_id, 17, &locus_id),
    ] {
        assert_eq!(stutter.catalog(), F4_TECHNICAL_RNG_CATALOG_V4);
    }
}

/// **Cada decisión de la etapa abre su propia clave.**
///
/// Sin esto, dos decisiones distintas de la misma familia podrían compartir
/// subsemilla y el sorteo de una arrastraría al de la otra.
#[test]
fn every_amplification_decision_opens_a_distinct_key() {
    let (event_id, library_id, _) = identifiers();
    let other_event = EventId::parse("reads-other".to_owned(), "test").expect("evento");
    let other_library = LibraryId::parse("LIB_0002".to_owned(), "test").expect("biblioteca");
    let mut keys = std::collections::BTreeSet::new();
    let mut total = 0_usize;
    for (event, library) in [
        (&event_id, &library_id),
        (&other_event, &library_id),
        (&event_id, &other_library),
    ] {
        for ordinal in [0_u64, 17] {
            for decision in [
                F4TechnicalDecision::amplification_efficiency_class(event, library, ordinal),
                F4TechnicalDecision::amplification_early_copies(event, library, ordinal),
                F4TechnicalDecision::amplification_gc_accept(event, library, ordinal, 0),
                F4TechnicalDecision::amplification_gc_accept(event, library, ordinal, 1),
                F4TechnicalDecision::amplification_family_size(event, library, ordinal),
            ] {
                keys.insert(key(decision));
                total += 1;
            }
        }
    }
    assert_eq!(keys.len(), total, "cada clave tiene que ser distinta");
}

/// **Vectores literales de subsemilla de la etapa nueva.**
///
/// # Por qué hacían falta
///
/// La etapa anterior —el tartamudeo— tiene sus vectores congelados y esta no
/// los tenía: solo se comprobaba que las claves fueran **distintas entre sí**,
/// que es una propiedad mucho más débil. Con los literales, cualquier cambio
/// en el orden de los componentes, en el nombre de la etapa o en el ámbito de
/// derivación muere aquí.
#[test]
fn the_amplification_keys_match_their_frozen_vectors() {
    let (event_id, library_id, _) = identifiers();
    let attempt = std::num::NonZeroU64::new(1).expect("intento");
    let vectors = [
        (
            F4TechnicalDecision::amplification_efficiency_class(&event_id, &library_id, 17),
            "2ba424972c7eefa98a78e78870974ee9196cf50e5744216a02965a0193f14178",
        ),
        (
            F4TechnicalDecision::amplification_early_copies(&event_id, &library_id, 17),
            "c92a3e587fea072c236237ec2e44d5cbc2e8c216a957e858e7e2cf94e1b7cfa6",
        ),
        (
            F4TechnicalDecision::amplification_gc_accept(&event_id, &library_id, 17, 0),
            "1f29adb2454cc3029dcc8b38bcfb43dc1ab6b190b7ffdef9f10d93f063a255e3",
        ),
        (
            F4TechnicalDecision::amplification_family_size(&event_id, &library_id, 17),
            "3a54f63d295d4da4adc504a250a53da3562dbbbff5f430a13df6f28eb6c88ea5",
        ),
        (
            F4TechnicalDecision::fragment_retry(
                &event_id,
                &library_id,
                17,
                FragmentRetryDecision::InsertLength,
                attempt,
            ),
            "52164f33459b55e5f3ebe1c1ead086b94fba7b43ad784a50aeef45dc10b48c8b",
        ),
    ];
    for (decision, expected) in vectors {
        assert_eq!(key(decision), expected);
    }
}

/// **El intento cero conserva literalmente la clave de hoy.**
///
/// # Qué afirma de verdad, y por qué la desigualdad no bastaba
///
/// La primera versión de esta prueba solo comprobaba que la clave de reintento
/// **difiriera** de la heredada, y su título prometía mucho más. Una edición
/// que hiciera que el intento cero usara ya las claves de reintento la dejaba
/// pasar tan campante —se midió—. Lo que de verdad ata el contrato es el
/// **literal**: las tres claves heredadas valen exactamente lo que valían antes
/// de que la amplificación existiera.
///
/// Y esa es la propiedad que hace posible la reducción al caso sin amplificar:
/// con curva neutra no se abre ningún reintento y el sorteo del fragmento es el
/// mismo de siempre.
#[test]
fn the_first_attempt_keeps_the_inherited_fragment_keys() {
    let (event_id, library_id, _) = identifiers();
    let attempt = std::num::NonZeroU64::new(1).expect("intento");
    let inherited_vectors = [
        (
            F4TechnicalDecision::fragment_insert_length(&event_id, &library_id, 17),
            FragmentRetryDecision::InsertLength,
            "54d91f60d873bc3250330a8d4a3c83440b1092677cefcc5979e0100d3fd1835f",
        ),
        (
            F4TechnicalDecision::fragment_copy(&event_id, &library_id, 17),
            FragmentRetryDecision::Copy,
            "fff42edbc9fd9971b6ebacc491e25853038a59f258263a4a09af814e15d2a07b",
        ),
        (
            F4TechnicalDecision::fragment_start(&event_id, &library_id, 17),
            FragmentRetryDecision::Start,
            "13e760451815fc226d4c12d0a88cfcd041ca433a0ed4d2e0691de18c57fd4a24",
        ),
    ];
    for (inherited, retried, frozen) in inherited_vectors {
        let inherited_key = key(inherited);
        assert_eq!(
            inherited_key, frozen,
            "la clave heredada del intento cero tiene que ser la de siempre"
        );
        let retry_key = key(F4TechnicalDecision::fragment_retry(
            &event_id,
            &library_id,
            17,
            retried,
            attempt,
        ));
        assert_ne!(
            inherited_key, retry_key,
            "un reintento no puede repetir la clave del intento cero"
        );
    }
}
