//! F4.9.10: la cota de recursos **publicada, no impuesta**.
//!
//! # Qué se demuestra aquí y qué no
//!
//! Se demuestra que las unidades de operación **realmente gastadas** por la
//! amplificación caben dentro de la cota que el manifiesto publica, y que el
//! término de los intentos de aceptación por composición es **necesario**: sin
//! él la cota se queda corta sobre un fixture saturado.
//!
//! No se demuestra nada sobre el tiempo de reloj ni sobre la memoria real: la
//! cota es un valor **calculado sin ejecutar**, y el manifiesto lo publica
//! marcado como no medido para que nadie lo lea como consumo observado.
//!
//! # Por qué el fixture tiene que estar saturado
//!
//! **Con curva neutra el tope de intentos no se ejercita**: cada familia
//! acepta su molécula al primer intento, que es el trabajo que F4 ya hacía, y
//! el término de los intentos no contribuye nada. Un mutante que lo suprimiera
//! de la fórmula sobreviviría. La última prueba de este fichero lo comprueba en
//! los dos sentidos.
//!
//! # Lo que este fixture no puede ejercitar, dicho sin adornos
//!
//! El término `2·c·T` de la fórmula, el de los tractos de repetición contenidos
//! en un inserto, **vale cero en esta entrega y no hay fixture que lo mueva**:
//! la preparación amplificada rechaza por ahora los escenarios que declaran
//! loci de repetición, de modo que `T = 0` siempre. Un mutante que suprimiera
//! ese término sobreviviría a toda esta batería, y así seguirá hasta que el
//! tartamudeo clonal heredable por linaje levante ese rechazo.

use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::pcr_amplification::{
    DeclaredPcrParameters, PcrPlan, PcrSettings, resource_bound,
};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::{
    TechnicalIdentityPlan, derive_technical_identity_plan,
};
use genome_population_simulator::technical_profile::{
    LibraryPreparation, PCR_MAXIMUM_EARLY_CYCLES, PCR_MAXIMUM_FAMILY_SIZE, PCR_MAXIMUM_GC_ATTEMPTS,
    resolve_f4_technical_profile_for,
};
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
const PAIRS: u64 = 200;
const MAXIMUM_INSERT_BASES: u64 = 450;

/// Referencia de composición equilibrada: la mitad de sus bases son G o C, de
/// modo que **toda** molécula cae en el mismo punto de la curva y la aceptación
/// no depende de dónde se muestreó.
fn balanced_reference(length: u64) -> Reference {
    let bases: Vec<u8> = (0..length)
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
    Reference::new(
        ReferenceId::parse("reference-main", "test").expect("referencia"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("referencia")
}

fn identities(reference: &Reference, individual: &Individual) -> TechnicalIdentityPlan {
    derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("semilla"),
        reference.id(),
        EventId::parse("wgs-pe", "test").expect("evento"),
        &[individual.id().clone()],
    )
    .expect("identidades")
}

/// Ajustes con los tres topes **en sus máximos de contrato**.
fn saturated_settings(relative_efficiency: f64) -> PcrSettings {
    PcrSettings::new(
        &DeclaredPcrParameters {
            duplicate_fraction: 0.05,
            cycles: 10,
            efficiency: 0.9,
            early_cycles: Some(PCR_MAXIMUM_EARLY_CYCLES),
            max_family_size: Some(PCR_MAXIMUM_FAMILY_SIZE),
            gc_max_attempts: Some(PCR_MAXIMUM_GC_ATTEMPTS),
            gc_curve: Some(vec![(0.0, relative_efficiency), (1.0, relative_efficiency)]),
            efficiency_classes: Vec::new(),
        },
        1e-9,
    )
    .expect("ajustes de amplificación")
}

/// Gasta una biblioteca completa y devuelve sus unidades y sus familias.
fn spend(settings: &PcrSettings) -> (u64, u64, u64) {
    let reference = balanced_reference(5_000);
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individuo"),
        &reference,
        Ploidy::new(1).expect("ploidía"),
    );
    let identities = identities(&reference, &individual);
    let profile =
        resolve_f4_technical_profile_for(LibraryPreparation::Pcr).expect("perfil amplificado");
    let sampler = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        100,
        PAIRS,
    )
    .expect("muestreador")
    .with_pcr_plan(PcrPlan::new(settings.clone()));
    let mut sampler = sampler;
    let mut emitted = 0_u64;
    for pair in &mut sampler {
        pair.expect("pareja");
        emitted += 1;
    }
    assert_eq!(emitted, PAIRS, "el contrato externo sigue siendo N parejas");
    let report = sampler.pcr_report();
    (
        report.operation_units(),
        report.families(),
        report.gc_attempts_total(),
    )
}

/// La cota tal como la fórmula la calcula, **sin** el término de los intentos.
///
/// Es el mutante que el plan manda cazar, escrito aquí para poder compararlo.
fn bound_without_the_attempts_term(settings: &PcrSettings, pairs: u64) -> u64 {
    let classes = settings.efficiency_classes().len() as u64;
    let early_states = 1_u64 << settings.early_cycles();
    let family_slots = u64::from(settings.max_family_size());
    let cycles = u64::from(settings.cycles());
    let per_pair = classes + early_states + family_slots + cycles + 1 + family_slots;
    pairs * per_pair + classes * u64::from(settings.early_cycles()) * early_states * early_states
}

/// **Sobre un fixture saturado, las unidades reales caben en la cota.**
///
/// El control positivo es que el contador **contó**: sin él, un contador que
/// no se incrementara nunca satisfaría la desigualdad trivialmente.
#[test]
fn the_real_units_fit_inside_the_published_bound_on_a_saturated_fixture() {
    let settings = saturated_settings(0.01);
    let (units, families, attempts) = spend(&settings);
    assert!(units > 0, "control positivo: el contador tiene que contar");
    assert!(families > 0, "tiene que haber familias");
    assert!(
        attempts > families,
        "el fixture saturado tiene que gastar más de un intento por familia: \
         {attempts} intentos para {families} familias"
    );
    let bound = resource_bound(&settings, PAIRS, 1, MAXIMUM_INSERT_BASES, 0);
    assert!(
        units <= bound.operation_units(),
        "las unidades reales ({units}) superan la cota publicada ({})",
        bound.operation_units()
    );
}

/// **El término de los intentos es necesario: sin él la cota se queda corta.**
///
/// Esta es la prueba que muere si alguien suprime `A × (P + 3·L_max)` de la
/// fórmula, porque entonces la cota publicada coincide con la mutante y la
/// desigualdad de la prueba anterior se rompe.
#[test]
fn removing_the_attempts_term_makes_the_bound_too_small_on_a_saturated_fixture() {
    let settings = saturated_settings(0.01);
    let (units, _, _) = spend(&settings);
    let mutant = bound_without_the_attempts_term(&settings, PAIRS);
    assert!(
        units > mutant,
        "sobre el fixture saturado, las unidades reales ({units}) tienen que superar la cota \
         mutante ({mutant}), o el mutante sobreviviría"
    );
    let honest = resource_bound(&settings, PAIRS, 1, MAXIMUM_INSERT_BASES, 0);
    assert!(
        honest.operation_units() > mutant,
        "la cota honesta tiene que ser mayor que la mutante"
    );
}

/// **Con curva neutra el mutante sobrevive**, y decirlo es parte del contrato.
///
/// Con la curva neutra cada familia acepta su molécula al primer intento, que
/// es el trabajo que F4 ya hacía y que esta cota deja fuera a propósito. Sin
/// el fixture saturado la prueba anterior **no discriminaría nada**.
#[test]
fn with_a_neutral_curve_the_attempts_term_is_never_exercised() {
    let settings = saturated_settings(1.0);
    assert!(
        settings.gc_curve().is_neutral(),
        "la curva de este caso tiene que ser neutra"
    );
    let (units, families, attempts) = spend(&settings);
    assert_eq!(
        attempts, families,
        "con curva neutra cada familia gasta exactamente un intento"
    );
    let mutant = bound_without_the_attempts_term(&settings, PAIRS);
    assert!(
        units <= mutant,
        "con curva neutra las unidades reales ({units}) caben incluso en la cota mutante \
         ({mutant}): el mutante sobrevive, y por eso hace falta el fixture saturado"
    );
}

/// **El valor exacto de la cota, congelado.**
///
/// # Qué agujero cierra
///
/// Todas las demás pruebas de este fichero comparan la cota con **desigualdades**:
/// que las unidades reales quepan dentro, que la mutante sea menor. Ninguna
/// mira su **valor**. Se midió que borrando los términos `+ c + 1` de la fórmula
/// —los de la divergencia de linaje y el progenitor— la batería entera seguía
/// pasando con 1413 pruebas y cero fallos.
///
/// Aquí el valor está escrito a mano, término a término, de modo que quitar
/// cualquiera de ellos mata la prueba aunque la desigualdad siga cumpliéndose.
///
/// # Los términos que este fixture **no** ejercita, dicho sin adornos
///
/// `2·c·T` —los tractos de repetición— y `c + 1` —la divergencia de linaje y el
/// progenitor— corresponden a trabajo que **desde la entrega 4 sí existe**: la
/// preparación amplificada ya admite escenarios con loci declarados y el
/// tartamudeo clonal los gasta de verdad. Lo que pasa es que **este fixture no
/// declara ninguno**: su individuo se construye a mano, sin alelos de
/// repetición, de modo que su contador nunca suma esos términos y ninguna
/// desigualdad de este fichero puede ver si desaparecen de la fórmula.
///
/// Es un hueco residual y queda dicho: lo único que sujeta esos dos términos
/// aquí sigue siendo **esta igualdad exacta**, escrita a mano. Quien quiera la
/// desigualdad con tractos de verdad tiene que montar un individuo con un
/// alelo de repetición, que es trabajo de otro fixture.
#[test]
fn the_bound_has_the_exact_value_its_formula_says() {
    let settings = saturated_settings(0.01);
    // Los términos, escritos a mano y por separado.
    let classes = 1_u64; // una sola clase de eficiencia
    let early_states = 1_u64 << PCR_MAXIMUM_EARLY_CYCLES; // 2^6 = 64
    let attempts = u64::from(PCR_MAXIMUM_GC_ATTEMPTS); // 16
    let copies = 1_u64; // una copia haplotípica
    let family_slots = u64::from(PCR_MAXIMUM_FAMILY_SIZE); // 256
    let cycles = 10_u64;
    // Cero porque **este fixture** no declara ninguno: su individuo se
    // construye a mano, sin alelos de repetición. Desde la entrega 4 la
    // amplificación sí admite tractos, y la igualdad exacta de más abajo con
    // `tracts = 3` es la que sujeta su término.
    let tracts = 0_u64;
    let per_pair = classes
        + early_states
        + attempts * (copies + 3 * MAXIMUM_INSERT_BASES)
        + family_slots
        + 2 * cycles * tracts
        + cycles
        + 1
        + family_slots;
    let dynamic_programming =
        classes * u64::from(PCR_MAXIMUM_EARLY_CYCLES) * early_states * early_states;
    let expected = PAIRS * per_pair + dynamic_programming;
    let bound = resource_bound(&settings, PAIRS, copies, MAXIMUM_INSERT_BASES, tracts);
    assert_eq!(
        bound.operation_units(),
        expected,
        "la cota publicada tiene que ser exactamente la de su fórmula"
    );
    // Control positivo de los dos términos que ninguna desigualdad ve: quitarlos
    // cambia el valor, de modo que la igualdad de arriba los sujeta de verdad.
    assert_ne!(
        expected,
        PAIRS * (per_pair - cycles - 1) + dynamic_programming
    );
    assert_ne!(
        resource_bound(&settings, PAIRS, copies, MAXIMUM_INSERT_BASES, 3).operation_units(),
        expected,
        "el término de los tractos tiene que mover la cota cuando T deja de ser cero"
    );

    // **La memoria del linaje, también escrita a mano.** Lleva un factor
    // `c + 1` que la primera redacción de la cota no tenía: para que una copia
    // herede el prefijo de su progenitor hacen falta sus pasos **ciclo a
    // ciclo**, no solo su suma. Sin esta igualdad, quitar ese factor dejaría
    // la cota publicada por debajo del gasto real sin que nada muriera.
    let tracts = 3_u64;
    let memory = MAXIMUM_INSERT_BASES
        + 8 * (classes * early_states + family_slots + cycles + 1)
        + u64::from(PCR_MAXIMUM_FAMILY_SIZE) * tracts * (cycles + 1) * 8
        + 8 * (family_slots + 1)
        + 5 * tracts * 24;
    assert_eq!(
        resource_bound(&settings, PAIRS, copies, MAXIMUM_INSERT_BASES, tracts).memory_bytes(),
        memory,
        "la cota de memoria tiene que ser exactamente la de su fórmula"
    );
    // Control positivo del factor nuevo: sin él la cifra es otra.
    assert_ne!(
        memory,
        MAXIMUM_INSERT_BASES
            + 8 * (classes * early_states + family_slots + cycles + 1)
            + u64::from(PCR_MAXIMUM_FAMILY_SIZE) * tracts * 8
            + 8 * (family_slots + 1)
            + 5 * tracts * 24
    );
    // Y control positivo del término de los recuentos por locus, que es el que
    // hace cierta la frase «no depende del genoma»: sin él la cifra también es
    // otra.
    assert_ne!(
        memory,
        MAXIMUM_INSERT_BASES
            + 8 * (classes * early_states + family_slots + cycles + 1)
            + u64::from(PCR_MAXIMUM_FAMILY_SIZE) * tracts * (cycles + 1) * 8
            + 8 * (family_slots + 1)
    );
}

/// La cota de memoria es **independiente del número de parejas**, y la de
/// operaciones crece con él.
#[test]
fn the_memory_bound_does_not_grow_with_the_pair_count() {
    let settings = saturated_settings(0.01);
    let small = resource_bound(&settings, 1_000, 2, MAXIMUM_INSERT_BASES, 0);
    let large = resource_bound(&settings, 1_000_000, 2, MAXIMUM_INSERT_BASES, 0);
    assert_eq!(
        small.memory_bytes(),
        large.memory_bytes(),
        "la memoria añadida no puede depender de N"
    );
    assert!(large.operation_units() > small.operation_units());
    assert_eq!(
        small.regeneration_passes(),
        2,
        "la cota declara que el recorrido regenera cada biblioteca dos veces"
    );
}

/// **La programación dinámica no puede volver a construirse dentro de la
/// familia.**
///
/// # Qué agujero cierra
///
/// La primera versión de este módulo reconstruía la distribución de los
/// primeros ciclos **una vez por familia**, mientras la cota publicada
/// declaraba haberla calculado una sola vez por clase. Ninguna prueba podía
/// verlo, porque el contador instrumentado sumaba los términos que la fórmula
/// declara en vez de medir el trabajo: estaba definido para cuadrar con ella.
///
/// El ancla estructural es esta: el muestreador **no nombra** las funciones que
/// construyen la tabla. Si alguien volviera a llamarlas desde la familia,
/// aparecerían en su texto y esta prueba moriría.
///
/// # Esta prueba es la segunda red, no la primera, y hay que decirlo
///
/// Es una búsqueda **léxica** y por tanto **se evade importando la función con
/// otro nombre**: se midió, y con esa evasión esta prueba sobrevive tan campante
/// mientras la tabla se reconstruye una vez por familia. La red que de verdad
/// cierra el hueco es el **contador** de `tests/pcr_dp_build_count.rs`, que
/// cuenta construcciones y no nombres. Esta se conserva porque es instantánea y
/// cuesta cero, no porque baste.
#[test]
fn the_sampler_never_rebuilds_the_early_copy_tables() {
    const SOURCE: &str = include_str!("../src/wgs_sampling.rs");
    const FORBIDDEN: [&str; 3] = [
        "early_copy_distribution(",
        "size_biased(",
        "mean_early_copies(",
    ];
    // Control positivo de la búsqueda: sobre un texto que sí las contiene,
    // cada patrón aparece. Sin esto, una búsqueda rota daría el mismo
    // resultado que un muestreador limpio.
    let sample = "early_copy_distribution(3, 0.9); size_biased(&d); mean_early_copies(3, 0.9);";
    for pattern in FORBIDDEN {
        assert!(
            sample.contains(pattern),
            "el control positivo no encuentra {pattern}: la búsqueda está rota"
        );
    }
    for pattern in FORBIDDEN {
        assert!(
            !SOURCE.contains(pattern),
            "src/wgs_sampling.rs vuelve a construir la tabla de los primeros ciclos con \
             {pattern}: eso se paga una vez por clase, no una vez por familia"
        );
    }
    // Y la tabla está donde tiene que estar: el muestreador la consulta.
    assert!(
        SOURCE.contains("early_copy_table("),
        "el muestreador tiene que consultar la tabla precalculada"
    );
}

/// **El término de la programación dinámica crece con el número de clases y
/// no con el de familias.**
///
/// Es la otra mitad del contrato: el `Q` de `Q × c_early × 4^c_early` es real.
#[test]
fn the_dynamic_programming_cost_scales_with_the_classes_and_not_the_families() {
    let one = saturated_settings(1.0);
    let two = PcrSettings::new(
        &DeclaredPcrParameters {
            duplicate_fraction: 0.05,
            cycles: 10,
            efficiency: 0.9,
            early_cycles: Some(PCR_MAXIMUM_EARLY_CYCLES),
            max_family_size: Some(PCR_MAXIMUM_FAMILY_SIZE),
            gc_max_attempts: Some(PCR_MAXIMUM_GC_ATTEMPTS),
            gc_curve: None,
            efficiency_classes: vec![(0.95, 0.5), (0.6, 0.5)],
        },
        1e-9,
    )
    .expect("dos clases");
    assert_eq!(one.efficiency_classes().len(), 1);
    assert_eq!(two.efficiency_classes().len(), 2);
    assert!(one.early_copy_build_units() > 0, "control positivo");
    assert_eq!(
        two.early_copy_build_units(),
        2 * one.early_copy_build_units(),
        "dos clases cuestan exactamente el doble que una"
    );
    // Y cada clase tiene su tabla, distinta de la de la otra: una tabla
    // compartida daría la misma distribución a eficiencias distintas.
    let first = two.early_copy_table(0).expect("tabla de la primera clase");
    let second = two.early_copy_table(1).expect("tabla de la segunda clase");
    assert_ne!(first.size_biased_weights(), second.size_biased_weights());
    assert!(two.early_copy_table(2).is_none());
}
