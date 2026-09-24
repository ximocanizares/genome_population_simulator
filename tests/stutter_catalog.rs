//! F4.8.1: el catálogo de ámbitos aleatorios `f4-technical-v4`.
//!
//! # Qué extiende y qué no toca
//!
//! El catálogo nuevo añade la etapa `stutter` con los componentes
//! `[event_id, library_id, stutter, template_ordinal, locus_id, decision]`, y
//! **no toca ninguna de las trece decisiones heredadas**: esas siguen
//! derivándose bajo `f4-technical-v3`, de modo que ninguna subsemilla ya
//! publicada cambia. Quien quiera comprobar esa parte tiene los vectores
//! literales de las trece en `tests/f4_rng_catalog.rs`, que siguen intactos:
//! **cualquier renumeración o reordenación de las etapas existentes los rompe**,
//! porque el orden de los componentes entra literalmente en el material.
//!
//! # Por qué el ámbito se indexa por locus y no por posición
//!
//! Para que declarar un locus nuevo en otro cromosoma no desplace las
//! decisiones de los que ya había. Es lo que comprueba
//! `the_scope_is_indexed_by_locus_and_not_by_position`.

use genome_population_simulator::identity::{EventId, LibraryId, RepeatLocusId};
use genome_population_simulator::randomness::{
    F4_TECHNICAL_RNG_CATALOG, F4_TECHNICAL_RNG_CATALOG_V3, F4_TECHNICAL_RNG_CATALOG_V4,
    F4TechnicalDecision, RootSeed, derive_f4_technical_subseed,
};

fn root() -> RootSeed {
    RootSeed::from_bytes([
        0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07, 0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e,
        0x0f, 0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17, 0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d,
        0x1e, 0x1f,
    ])
}

fn identifiers() -> (EventId, LibraryId, RepeatLocusId) {
    (
        EventId::parse("event-4", "test").expect("event"),
        LibraryId::parse("library-4", "test").expect("library"),
        RepeatLocusId::parse("STR-1", "test").expect("locus"),
    )
}

/// Vectores congelados de la terna evento/biblioteca/locus.
///
/// Cualquier cambio del orden o del literal de los componentes de la etapa
/// nueva mueve estas dos cifras. Están escritas a mano, no calculadas por el
/// mismo código que se comprueba, que es lo que las convierte en un vector y no
/// en un espejo.
#[test]
fn the_stutter_stage_has_literal_subseed_vectors() {
    let (event_id, library_id, locus_id) = identifiers();
    let bernoulli = derive_f4_technical_subseed(
        root(),
        F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, 17, &locus_id),
    )
    .expect("subsemilla del ensayo");
    let step = derive_f4_technical_subseed(
        root(),
        F4TechnicalDecision::stutter_step(&event_id, &library_id, 17, &locus_id),
    )
    .expect("subsemilla del paso");
    assert_eq!(
        bernoulli.to_hex(),
        "c0caba2a3c05904f30a69a71d14fef17a148b674130df15f762acc7dd594833e"
    );
    assert_eq!(
        step.to_hex(),
        "d9e6dcbe2a9ecab99b8f5f42c77bfaf3d5d0ec8860c9316e4d16c49e4d7b5514"
    );
}

/// Las dos decisiones de la misma molécula y el mismo locus son **distintas**.
///
/// Sin esto, una implementación que usara la misma clave para el ensayo y para
/// el paso pasaría desapercibida y ligaría dos sorteos que tienen que ser
/// independientes.
#[test]
fn the_trial_and_the_step_do_not_share_a_key() {
    let (event_id, library_id, locus_id) = identifiers();
    let bernoulli = derive_f4_technical_subseed(
        root(),
        F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, 17, &locus_id),
    )
    .expect("ensayo");
    let step = derive_f4_technical_subseed(
        root(),
        F4TechnicalDecision::stutter_step(&event_id, &library_id, 17, &locus_id),
    )
    .expect("paso");
    assert_ne!(bernoulli.to_hex(), step.to_hex());
}

/// El ámbito se indexa por **identificador de locus**, no por la posición del
/// tracto dentro de la molécula: dos loci distintos de la misma molécula tienen
/// claves distintas, y el mismo locus en dos moléculas también.
#[test]
fn the_scope_is_indexed_by_locus_and_not_by_position() {
    let (event_id, library_id, first) = identifiers();
    let second = RepeatLocusId::parse("STR-2", "test").expect("locus");
    let key = |locus: &RepeatLocusId, ordinal: u64| {
        derive_f4_technical_subseed(
            root(),
            F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, ordinal, locus),
        )
        .expect("subsemilla")
        .to_hex()
    };
    assert_ne!(key(&first, 17), key(&second, 17));
    assert_ne!(key(&first, 17), key(&first, 18));
}

/// **La etapa nueva nace en el v4 y las heredadas siguen en el v3.**
///
/// Es la propiedad que hace que añadir el tartamudeo no desplace ninguna
/// decisión anterior. Si el catálogo del tartamudeo fuera el v3, la clave de la
/// etapa `stutter` compartiría ámbito de derivación con las trece heredadas;
/// si el de las heredadas pasara al v4, todas sus subsemillas cambiarían.
#[test]
fn only_the_stutter_stage_belongs_to_the_new_catalog() {
    let (event_id, library_id, locus_id) = identifiers();
    assert_eq!(F4_TECHNICAL_RNG_CATALOG, F4_TECHNICAL_RNG_CATALOG_V3);
    assert_eq!(F4_TECHNICAL_RNG_CATALOG_V3, "f4-technical-v3");
    assert_eq!(F4_TECHNICAL_RNG_CATALOG_V4, "f4-technical-v4");
    assert_eq!(
        F4TechnicalDecision::stutter_bernoulli(&event_id, &library_id, 17, &locus_id).catalog(),
        F4_TECHNICAL_RNG_CATALOG_V4
    );
    assert_eq!(
        F4TechnicalDecision::stutter_step(&event_id, &library_id, 17, &locus_id).catalog(),
        F4_TECHNICAL_RNG_CATALOG_V4
    );
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
