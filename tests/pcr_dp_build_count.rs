//! **La tabla de los primeros ciclos se construye una vez por clase.**
//!
//! # Por qué este fichero está solo
//!
//! El contador de construcciones es global al proceso. Cada fichero de prueba
//! de integración es su propio proceso, de modo que dejar aquí **una sola**
//! prueba es lo que garantiza que nadie más lo toque mientras mide. Meterla con
//! las demás la haría depender del orden en que el corredor las lance, que es
//! justo la clase de prueba intermitente que este proyecto no quiere.
//!
//! # Qué cierra, y qué no cerraba la guarda anterior
//!
//! La guarda anterior era **léxica**: comprobaba que el muestreador no nombrara
//! las funciones que construyen la tabla. Se midió que **se evade importando la
//! función con otro nombre** y llamándola igual: las dos baterías seguían en
//! verde con la tabla reconstruyéndose una vez por familia. Un contador no se
//! evade renombrando.

use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::pcr_amplification::{
    DeclaredPcrParameters, PcrPlan, PcrSettings, early_copy_build_count,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::technical_profile::{
    LibraryPreparation, resolve_f4_technical_profile_for,
};
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn settings(classes: Vec<(f64, f64)>) -> PcrSettings {
    PcrSettings::new(
        &DeclaredPcrParameters {
            duplicate_fraction: 0.3,
            cycles: 10,
            efficiency: 0.9,
            early_cycles: Some(6),
            max_family_size: Some(64),
            gc_max_attempts: Some(8),
            gc_curve: None,
            efficiency_classes: classes,
        },
        1e-9,
    )
    .expect("ajustes")
}

/// Gasta una biblioteca entera y devuelve cuántas familias abrió.
fn spend(settings: &PcrSettings, pairs: u64) -> u64 {
    let bases: Vec<u8> = (0..5_000_u64)
        .map(|index| match index % 4 {
            0 => b'G',
            1 => b'A',
            2 => b'C',
            _ => b'T',
        })
        .collect();
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("secuencia"),
        FastaHeader::parse(b"chr1").expect("cabecera"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("ADN"),
    )
    .expect("secuencia de referencia");
    let reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("referencia"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("referencia");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individuo"),
        &reference,
        Ploidy::new(1).expect("ploidía"),
    );
    let identities = derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("semilla"),
        reference.id(),
        EventId::parse("wgs-pe", "test").expect("evento"),
        &[individual.id().clone()],
    )
    .expect("identidades");
    let profile =
        resolve_f4_technical_profile_for(LibraryPreparation::Pcr).expect("perfil amplificado");
    let mut sampler = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        100,
        pairs,
    )
    .expect("muestreador")
    .with_pcr_plan(PcrPlan::new(settings.clone()));
    for pair in &mut sampler {
        pair.expect("pareja");
    }
    sampler.pcr_report().families()
}

/// **Una tanda entera construye la tabla tantas veces como clases tiene.**
#[test]
fn the_early_copy_table_is_built_once_per_efficiency_class() {
    // Una clase.
    let before = early_copy_build_count();
    let one = settings(Vec::new());
    let built_for_one = early_copy_build_count() - before;
    assert_eq!(built_for_one, 1, "una clase, una construcción");

    // Y gastar la biblioteca entera **no construye ninguna más**.
    let before = early_copy_build_count();
    let families = spend(&one, 400);
    let built_while_sampling = early_copy_build_count() - before;
    // Control positivo: la tanda abrió muchas familias, de modo que un cero no
    // sale de no haber muestreado nada.
    assert!(
        families > 50,
        "la tanda tiene que abrir muchas familias: {families}"
    );
    assert_eq!(
        built_while_sampling, 0,
        "muestrear {families} familias ha reconstruido la tabla {built_while_sampling} veces; \
         se construye una vez por clase al resolver el plan, no una por familia"
    );

    // Dos clases, dos construcciones, y tampoco ninguna al muestrear.
    let before = early_copy_build_count();
    let two = settings(vec![(0.95, 0.7), (0.6, 0.3)]);
    assert_eq!(early_copy_build_count() - before, 2, "dos clases, dos");
    let before = early_copy_build_count();
    let families = spend(&two, 400);
    assert!(families > 50);
    assert_eq!(early_copy_build_count() - before, 0);
}
