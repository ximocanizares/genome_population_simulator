//! F4.9.4, F4.9.5 y F4.9.6: las tres piezas matemáticas de la amplificación.
//!
//! # Qué hace independiente a cada oráculo de este fichero
//!
//! - La **inversión numérica** se compara contra el cálculo **directo** de las
//!   fracciones que el diseño ya tenía calculadas a mano, escritas aquí como
//!   literales. Comparar la inversión contra la propia función directa del
//!   módulo sería un oráculo circular.
//! - La **programación dinámica** de los primeros ciclos se compara contra una
//!   **enumeración exhaustiva de caminos** escrita aquí, que no comparte ni una
//!   línea con el módulo. Reusar la propia programación dinámica como oráculo
//!   sería circular.
//! - Los **pesos del tamaño de familia** se comparan contra los valores
//!   calculados a mano, y el último contra la masa restante, que es lo que
//!   distingue una cola pinzada de una truncada.

use genome_population_simulator::identity::{EventId, LibraryId};
use genome_population_simulator::pcr_amplification::{
    PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS, PcrAmplificationError, composition_factor,
    early_copy_distribution, expected_duplicate_fraction, family_size_weights, gc_fraction,
    invert_duplicate_fraction, mean_early_copies, size_biased,
};
use genome_population_simulator::randomness::{
    F4TechnicalDecision, RootSeed, categorical_f64, derive_f4_technical_subseed,
};
use genome_population_simulator::technical_profile::PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE;

/// El texto del módulo, incorporado en tiempo de compilación para que la
/// búsqueda negativa no dependa del directorio desde el que se lance la prueba.
const MODULE_SOURCE: &str = include_str!("../src/pcr_amplification.rs");

/// Las tres parejas `(m, f(m))` que el diseño trae calculadas en el sentido
/// **directo**, con todas las cifras de un `f64`.
///
/// Son el oráculo de la inversión: se le da `f` y tiene que devolver `m`. No se
/// obtienen llamando al módulo.
const WORKED_PAIRS: [(f64, f64); 3] = [
    (0.02, 0.009_933_665_337_762_632),
    (0.1, 0.048_374_180_359_595_176),
    (0.5, 0.213_061_319_425_266_85),
];

// ---------------------------------------------------------------------------
// F4.9.4 — inversión numérica.
// ---------------------------------------------------------------------------

/// **La inversión devuelve las parejas por molécula del cálculo directo.**
#[test]
fn inverting_the_worked_fractions_recovers_their_mean_pairs_per_molecule() {
    for (expected_mean, fraction) in WORKED_PAIRS {
        let inversion = invert_duplicate_fraction(fraction).expect("fracción alcanzable");
        let observed = inversion.mean_pairs_per_molecule();
        let relative = (observed - expected_mean).abs() / expected_mean;
        assert!(
            relative < 1e-12,
            "la inversión de {fraction} dio {observed} y el cálculo directo dice {expected_mean}"
        );
    }
}

/// **La fracción cero devuelve cero exacto y sin iterar.**
///
/// No es un adorno: sin el retorno anticipado la bisección devuelve un residuo
/// del orden de `10⁻¹⁸`, y ese residuo hace que el tamaño de familia deje de
/// ser trivial, se consuma aleatoriedad y cambien **todos** los FASTQ. Es la
/// clase de defecto que no da error.
#[test]
fn a_zero_duplicate_fraction_returns_exactly_zero_without_iterating() {
    let inversion = invert_duplicate_fraction(0.0).expect("cero es alcanzable");
    assert_eq!(
        inversion.mean_pairs_per_molecule(),
        0.0,
        "la fracción cero tiene que dar exactamente cero, no un residuo"
    );
    assert_eq!(
        inversion.iterations(),
        0,
        "la fracción cero no itera: sale por retorno anticipado"
    );
    // Control positivo del mismo par de aserciones: con una fracción distinta
    // de cero las dos cifras cambian, de modo que no pasan por ser triviales.
    let other = invert_duplicate_fraction(0.05).expect("alcanzable");
    assert!(other.mean_pairs_per_molecule() > 0.0);
    assert_eq!(
        other.iterations(),
        PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS
    );
}

/// **Una fracción inalcanzable se rechaza nombrando el máximo.**
///
/// No se pinza en silencio: el usuario tiene que enterarse de que pidió algo
/// que el tope de parejas por molécula no puede producir.
#[test]
fn an_unreachable_duplicate_fraction_is_rejected_naming_the_maximum() {
    let reachable = expected_duplicate_fraction(PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE);
    let error = invert_duplicate_fraction(0.999).expect_err("0,999 no es alcanzable");
    match error {
        PcrAmplificationError::UnreachableDuplicateFraction {
            value,
            maximum_reachable,
            maximum_mean_pairs_per_molecule,
        } => {
            assert_eq!(value, 0.999);
            assert_eq!(maximum_reachable, reachable);
            assert_eq!(
                maximum_mean_pairs_per_molecule,
                PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE
            );
        }
        other => panic!("rechazo inesperado: {other:?}"),
    }
    // Control positivo: la fracción justo alcanzable **sí** se admite, de modo
    // que el rechazo no viene de que el módulo rechace cualquier cosa.
    assert!(invert_duplicate_fraction(reachable).is_ok());
}

/// **La función directa es estrictamente creciente**, que es la precondición
/// que hace única la inversa.
#[test]
fn the_direct_fraction_is_strictly_increasing_over_a_hundred_point_mesh() {
    let mut previous = expected_duplicate_fraction(0.0);
    assert_eq!(previous, 0.0);
    for step in 1..=100 {
        let mean = PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE * f64::from(step) / 100.0;
        let current = expected_duplicate_fraction(mean);
        assert!(
            current > previous,
            "f no crece en m = {mean}: {previous} -> {current}"
        );
        previous = current;
    }
    assert!(previous < 1.0, "f no puede alcanzar 1");
}

/// **El número de iteraciones es constante** para toda entrada distinta de
/// cero.
///
/// Es lo que caza una bisección con criterio de parada por tolerancia: esa
/// variaría el número de pasos con el objetivo.
#[test]
fn the_iteration_count_is_constant_for_every_non_zero_target() {
    let mut seen = std::collections::BTreeSet::new();
    for step in 1..=200_u32 {
        let target = f64::from(step) / 250.0;
        let inversion = invert_duplicate_fraction(target).expect("alcanzable");
        seen.insert(inversion.iterations());
    }
    assert_eq!(
        seen,
        std::collections::BTreeSet::from([PCR_DUPLICATE_FRACTION_INVERSION_ITERATIONS]),
        "el número de iteraciones tiene que ser el mismo para toda entrada: {seen:?}"
    );
}

/// **El módulo no llama a la exponencial de la matemática estándar.**
///
/// # Cuidado al escribir comentarios en el módulo vigilado
///
/// La búsqueda mira **todo el texto**, comentarios incluidos. Un comentario
/// que se limite a nombrar una de las funciones prohibidas —escribir
/// `.exp()` o `f64::exp` dentro de una frase explicativa— haría fallar esta
/// prueba sin que haya defecto ninguno. Es el precio de una búsqueda léxica y
/// queda dicho aquí para que quien la vea fallar sepa mirarlo primero: en el
/// módulo se habla de «la exponencial» y de `libm::exp`, nunca de las formas
/// prohibidas.
///
/// # Por qué hace falta una búsqueda en el texto y no basta una prueba numérica
///
/// En **una sola plataforma** las dos exponenciales dan lo mismo y ninguna
/// prueba numérica las distingue. Pero la bisección compara `f(mid)` con el
/// objetivo, y si el último bit de la exponencial varía entre plataformas la
/// comparación cae del otro lado y la bisección se va por otra rama **aunque
/// el número de iteraciones sea el mismo**. Una prueba de reproducibilidad
/// ejecutada dos veces en la misma máquina no discrimina nada sobre esto,
/// porque la exponencial del sistema también es determinista dentro de un
/// sistema.
#[test]
fn the_module_never_calls_the_standard_library_transcendentals() {
    const FORBIDDEN: [&str; 7] = [
        ".exp()", "f64::exp", ".ln()", ".log(", ".powf(", ".powi(", ".sqrt()",
    ];
    // Control positivo de la búsqueda misma: sobre un texto que sí las
    // contiene, cada patrón aparece. Sin esto, una búsqueda rota daría el
    // mismo resultado que un módulo limpio.
    let sample = "let a = x.exp(); let b = f64::exp(x); let c = x.ln(); let d = x.log(2.0); \
                  let e = x.powf(2.0); let f = x.powi(2); let g = x.sqrt();";
    for pattern in FORBIDDEN {
        assert!(
            sample.contains(pattern),
            "el control positivo no encuentra {pattern}: la búsqueda está rota"
        );
    }
    for pattern in FORBIDDEN {
        assert!(
            !MODULE_SOURCE.contains(pattern),
            "src/pcr_amplification.rs usa {pattern} de la matemática estándar"
        );
    }
    // Y la exponencial que sí usa es la de la librería fijada.
    assert!(
        MODULE_SOURCE.contains("libm::exp("),
        "la exponencial tiene que venir de libm"
    );
}

// ---------------------------------------------------------------------------
// F4.9.5 — copias tras los primeros ciclos.
// ---------------------------------------------------------------------------

/// Enumeración **exhaustiva** de los caminos del proceso de ramificación,
/// escrita aquí a mano.
///
/// Con `z` copias presentes, cada una se copia o no de forma independiente: son
/// `2^z` caminos, cada uno con su probabilidad. No hay ni coeficientes
/// binomiales ni programación dinámica; es la definición desnuda del proceso.
fn enumerate_paths(cycles: u32, efficiency: f64, states: usize) -> Vec<f64> {
    fn walk(copies: usize, remaining: u32, weight: f64, efficiency: f64, totals: &mut Vec<f64>) {
        if remaining == 0 {
            totals[copies - 1] += weight;
            return;
        }
        // Cada una de las `copies` moléculas decide por su cuenta.
        for mask in 0_u32..(1_u32 << copies) {
            let replicated = mask.count_ones() as usize;
            let path = efficiency.powi(replicated as i32)
                * (1.0 - efficiency).powi((copies - replicated) as i32);
            walk(
                copies + replicated,
                remaining - 1,
                weight * path,
                efficiency,
                totals,
            );
        }
    }
    let mut totals = vec![0.0_f64; states];
    walk(1, cycles, 1.0, efficiency, &mut totals);
    totals
}

/// **La programación dinámica coincide con la enumeración a mano.**
#[test]
fn the_early_copy_distribution_matches_a_hand_written_path_enumeration() {
    let distribution = early_copy_distribution(3, 0.5);
    assert_eq!(distribution.len(), 8, "tres ciclos dan estados 1..=8");
    let enumerated = enumerate_paths(3, 0.5, 8);
    for (index, (computed, hand)) in distribution.iter().zip(&enumerated).enumerate() {
        assert!(
            (computed - hand).abs() < 1e-12,
            "el estado z = {} vale {computed} y la enumeración dice {hand}",
            index + 1
        );
    }
    // Control positivo: la distribución tiene que repartirse de verdad entre
    // varios estados, o la igualdad de arriba no compararía nada.
    let occupied = distribution.iter().filter(|value| **value > 1e-12).count();
    assert!(
        occupied > 3,
        "la distribución tiene que ocupar varios estados: {distribution:?}"
    );
    let total: f64 = distribution.iter().sum();
    assert!((total - 1.0).abs() < 1e-12, "la masa total vale {total}");
}

/// **La distribución sesgada por tamaño tiene la media que le corresponde.**
///
/// Una familia observada es más probable cuantas más copias tenga su molécula,
/// de modo que la media de la sesgada es `E[Z²]/E[Z]`, calculada aquí desde la
/// enumeración a mano.
#[test]
fn the_size_biased_distribution_has_the_mean_the_hand_enumeration_predicts() {
    let enumerated = enumerate_paths(3, 0.5, 8);
    let first: f64 = enumerated
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64) * probability)
        .sum();
    let second: f64 = enumerated
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64).powi(2) * probability)
        .sum();
    let expected_mean = second / first;
    let biased = size_biased(&early_copy_distribution(3, 0.5));
    let observed_mean: f64 = biased
        .iter()
        .enumerate()
        .map(|(index, probability)| ((index + 1) as f64) * probability)
        .sum();
    assert!(
        (observed_mean - expected_mean).abs() < 1e-12,
        "la media sesgada vale {observed_mean} y la enumeración dice {expected_mean}"
    );
    // Y el sesgo es **estricto**: la media sesgada supera a la sin sesgar. Sin
    // esta comprobación, una implementación que no sesgara pasaría siempre que
    // la aritmética cuadrase consigo misma.
    assert!(
        observed_mean > first,
        "la sesgada por tamaño tiene que tener media mayor: {observed_mean} vs {first}"
    );
    // `E[Z] = (1+e)^c_early`, calculado aparte.
    assert!((first - mean_early_copies(3, 0.5)).abs() < 1e-12);
}

/// Con cero ciclos explícitos la distribución es un punto en `z = 1`.
#[test]
fn zero_early_cycles_leave_a_single_copy() {
    let distribution = early_copy_distribution(0, 0.9);
    assert_eq!(distribution, vec![1.0]);
    assert_eq!(mean_early_copies(0, 0.9), 1.0);
}

// ---------------------------------------------------------------------------
// F4.9.6 — tamaño de familia.
// ---------------------------------------------------------------------------

/// **Los pesos coinciden con los calculados a mano y el último recoge toda la
/// masa de la cola.**
///
/// El discriminador es el último peso: una Poisson truncada **sin pinzar**
/// daría `0.000_964_674_870_809_652_7` en esa celda, y la masa que falta se
/// perdería en silencio.
#[test]
fn the_family_size_weights_match_the_hand_computed_clamped_poisson() {
    const HAND_COMPUTED: [f64; 3] = [
        0.857_488_774_053_024_8,
        0.128_623_316_107_953_72,
        0.012_862_331_610_795_37,
    ];
    const UNCLAMPED_TAIL: f64 = 0.000_964_674_870_809_652_7;
    let weights = family_size_weights(0.3, 4);
    assert_eq!(weights.len(), 4);
    for (index, expected) in HAND_COMPUTED.iter().enumerate() {
        assert!(
            (weights[index] - expected).abs() < 1e-15,
            "el peso de k = {} vale {} y a mano sale {expected}",
            index + 1,
            weights[index]
        );
    }
    let tail: f64 = 1.0 - HAND_COMPUTED.iter().sum::<f64>();
    assert!(
        (weights[3] - tail).abs() < 1e-15,
        "el último peso vale {} y la masa restante es {tail}",
        weights[3]
    );
    assert!(
        (weights[3] - UNCLAMPED_TAIL).abs() > 1e-6,
        "el último peso no puede ser el término sin pinzar {UNCLAMPED_TAIL}"
    );
    let total: f64 = weights.iter().sum();
    assert!((total - 1.0).abs() < 1e-12, "los pesos suman {total}");
}

/// **Con `λ = 0` el sorteo no consume ni una palabra de generador.**
///
/// Se cuenta comparando la siguiente palabra de dos generadores gemelos: si el
/// sorteo hubiese consumido algo, los dos dejarían de estar alineados. El
/// control positivo con `λ` distinto de cero demuestra que la comparación sí
/// detecta un consumo.
#[test]
fn a_zero_lambda_family_size_consumes_no_randomness() {
    let seed = RootSeed::from_bytes([7_u8; 32]);
    let event_id = EventId::parse("reads-main".to_owned(), "test").expect("evento");
    let library_id = LibraryId::parse("LIB_0001".to_owned(), "test").expect("biblioteca");
    let subseed = derive_f4_technical_subseed(
        seed,
        F4TechnicalDecision::amplification_family_size(&event_id, &library_id, 0),
    )
    .expect("subsemilla");

    let weights = family_size_weights(0.0, 8);
    assert_eq!(weights[0], 1.0);
    assert!(weights[1..].iter().all(|weight| *weight == 0.0));
    let mut untouched = subseed.local_rng();
    let mut drawn = subseed.local_rng();
    assert_eq!(categorical_f64(&mut drawn, &weights), Ok(0));
    assert_eq!(
        untouched.next_u64(),
        drawn.next_u64(),
        "con lambda cero el sorteo no puede consumir aleatoriedad"
    );

    // Control positivo: con varios pesos positivos el sorteo **sí** consume, de
    // modo que la comparación de arriba no pasa por no detectar nada.
    let busy_weights = family_size_weights(0.5, 8);
    assert!(busy_weights.iter().filter(|weight| **weight > 0.0).count() > 1);
    let mut untouched = subseed.local_rng();
    let mut drawn = subseed.local_rng();
    let _ = categorical_f64(&mut drawn, &busy_weights);
    assert_ne!(
        untouched.next_u64(),
        drawn.next_u64(),
        "con lambda positivo el sorteo tiene que consumir"
    );
}

/// El tamaño máximo de familia es un tope de verdad.
#[test]
fn the_family_size_weights_never_exceed_the_declared_maximum() {
    for maximum in [1_u32, 2, 16, 256] {
        let weights = family_size_weights(12.0, maximum);
        assert_eq!(weights.len(), maximum as usize);
        let total: f64 = weights.iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-12,
            "con máximo {maximum} los pesos suman {total}"
        );
    }
}

// ---------------------------------------------------------------------------
// F4.9.7 — composición de bases, la parte pura.
// ---------------------------------------------------------------------------

/// La fracción de GC se cuenta **sobre A/C/G/T**, y un fragmento sin ninguna
/// recibe el valor que no distorsiona nada.
#[test]
fn the_gc_fraction_ignores_bases_outside_acgt() {
    assert_eq!(gc_fraction(b"GCGC"), 1.0);
    assert_eq!(gc_fraction(b"ATAT"), 0.0);
    assert_eq!(gc_fraction(b"GCAT"), 0.5);
    // Las `N` no cuentan ni en el numerador ni en el denominador: la misma
    // secuencia con `N` intercaladas da la misma fracción.
    assert_eq!(gc_fraction(b"GNCNANT"), 0.5);
    assert_eq!(gc_fraction(b"NNNN"), 1.0);
}

/// **Con eficiencia relativa 1 el factor vale exactamente `1.0`.**
///
/// Es lo que hace que el ensayo de aceptación no consuma ni una palabra con la
/// curva neutra, y con ello que la reducción al caso sin amplificar produzca
/// los mismos bytes. Un factor de `0.999_999_999_999_999_9` rompería las dos
/// cosas sin dar ningún error.
#[test]
fn a_neutral_relative_efficiency_gives_exactly_one() {
    for efficiency in [0.5_f64, 0.9, 1.0, 0.123_456_789] {
        for cycles in [1_u32, 7, 30] {
            assert_eq!(
                composition_factor(efficiency, 1.0, cycles),
                1.0,
                "el factor neutro tiene que ser exactamente 1 con e = {efficiency} y {cycles} \
                 ciclos"
            );
        }
    }
    // Control positivo: con eficiencia relativa menor que 1 el factor baja de
    // verdad, y baja más cuantos más ciclos.
    let few = composition_factor(0.9, 0.2, 2);
    let many = composition_factor(0.9, 0.2, 20);
    assert!(
        few < 1.0 && many < few,
        "{many} tiene que ser menor que {few}"
    );
}
