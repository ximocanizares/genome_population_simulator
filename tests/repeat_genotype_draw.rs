//! F6.4.3.2: sorteo de los genotipos de repetición desde las frecuencias.
//!
//! # Cómo están elegidos los fixtures
//!
//! Lo que estas pruebas afirman es que el ámbito de azar depende de
//! **identidades estables** y no de la posición ni del orden. Para que puedan
//! fallar hacen falta tres cosas, y las tres están elegidas a propósito:
//!
//! 1. **El corpus tiene que producir alelos distintos.** Si todas las copias
//!    recibieran el mismo alelo, «los genotipos no cambian» sería cierto para
//!    cualquier implementación, incluida una que sorteara con el ordinal. Por
//!    eso cada prueba comprueba primero que entre los sorteos observados hay
//!    **más de un alelo**, y si no lo hubiera fallaría antes de comparar nada.
//! 2. **El individuo añadido va al principio**, no al final. Añadir al final
//!    deja intactos todos los ordinales anteriores, de modo que un ámbito
//!    indexado por orden sobreviviría.
//! 3. **El locus añadido va en un cromosoma que ordena antes.** Añadirlo
//!    después no desplaza nada y tampoco cazaría un ámbito posicional.

use std::collections::BTreeMap;
use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::{EffectiveConfig, build_effective_config};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::identity::{ReferenceId, SequenceId};
use genome_population_simulator::initial_individual_integration::IncompatibleGenotypeAuthorization;
use genome_population_simulator::initial_individual_materialization::{
    InitialIndividualsMaterialization, materialize_initial_individuals,
};
use genome_population_simulator::population_snvs::{
    PreparedPopulationSnvs, TSV_HEADER, parse_canonical_tsv, prepare_population_snvs_from_reference,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
/// Cien bases: veinte `G` de flanco, sesenta de tracto `(CA)30` y veinte más.
const SEQUENCE_LENGTH: u64 = 100;

fn sequence_bases() -> Vec<u8> {
    let mut bases = vec![b'G'; 20];
    bases.extend(std::iter::repeat_n(*b"CA", 30).flatten());
    bases.extend(vec![b'G'; 20]);
    assert_eq!(bases.len() as u64, SEQUENCE_LENGTH);
    bases
}

/// Dos cromosomas, y `chr0` **ordena antes** que `chr1`.
fn reference() -> Reference {
    let bases = sequence_bases();
    let sequence = |id: &str| {
        ReferenceSequence::new(
            SequenceId::parse(id, "test").expect("sequence id"),
            FastaHeader::parse(id.as_bytes()).expect("header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(&bases).expect("DNA"),
        )
        .expect("reference sequence")
    };
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        ReferenceOrigin::Generated,
        vec![sequence("chr0"), sequence("chr1")],
    )
    .expect("reference")
}

/// Tres alelos con pesos **desiguales**, para que un reparto uniforme que
/// ignorase los pesos no produjera la misma distribución.
const ALLELES: &str = r#"
[[reference.repeat_loci.alleles]]
units = 30
frequency = 0.5

[[reference.repeat_loci.alleles]]
units = 25
frequency = 0.3

[[reference.repeat_loci.alleles]]
units = 20
frequency = 0.2
"#;

fn locus(id: &str, sequence: &str) -> String {
    format!(
        r#"
[[reference.repeat_loci]]
id = "{id}"
sequence = "{sequence}"
start = 21
end = 80
motif = "CA"
min_longest_run_units = 3
{ALLELES}"#
    )
}

fn config(individual_ids: &[&str], loci: &str) -> EffectiveConfig {
    config_with_model(individual_ids, loci, "genotype_model = \"hwe\"")
}

/// La misma configuración con el modelo genotípico que se pida, escrito tal
/// cual en la cohorte (`genotype_model` y, si hace falta, el coeficiente).
fn config_with_model(individual_ids: &[&str], loci: &str, model: &str) -> EffectiveConfig {
    let ids = individual_ids
        .iter()
        .map(|id| format!("\"{id}\""))
        .collect::<Vec<_>>()
        .join(", ");
    let source = format!(
        r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25

[[reference.sequences]]
id = "chr0"
class = "chromosome"
length = {SEQUENCE_LENGTH}

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {SEQUENCE_LENGTH}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = {}
individual_ids = [{ids}]
ploidy = 2
{model}
"#,
        individual_ids.len()
    );
    let admitted = admit_toml_header(
        Path::new("repeat-genotype-draw.toml"),
        source.into_bytes(),
        ExpectedConfigKind::General,
    )
    .expect("header");
    build_effective_config(&parse_declared_config(admitted).expect("declared")).expect("effective")
}

fn matrix(reference: &Reference) -> PreparedPopulationSnvs {
    // La fila cae en la base 1 de `chr0`, dentro del flanco y **fuera** de los
    // dos tractos: un cambio de base dentro del intervalo de un locus impide su
    // alelo a propósito, y aquí lo que se mide es el sorteo.
    let bytes = format!(
        "{}\nVAR_1\tchr0\t1\tG\tA\tP1\t0.5\tfixed\n",
        TSV_HEADER.join("\t")
    );
    prepare_population_snvs_from_reference(
        &parse_canonical_tsv(bytes.as_bytes()).expect("TSV"),
        reference,
        1e-9,
        false,
    )
    .expect("prepared matrix")
}

fn materialize(
    config: &EffectiveConfig,
    reference: &Reference,
) -> InitialIndividualsMaterialization {
    materialize_initial_individuals(
        config,
        reference,
        &matrix(reference),
        IncompatibleGenotypeAuthorization::denied(),
        |_| Ok(Vec::new()),
    )
    .expect("materialized founders")
}

/// Los sorteos de un locus, indexados por **identificador de individuo y copia**
/// y nunca por posición en la lista.
fn genotypes(
    output: &InitialIndividualsMaterialization,
    locus_id: &str,
) -> BTreeMap<(String, u8), usize> {
    output
        .repeat_draws()
        .iter()
        .filter(|draw| draw.locus_id() == locus_id)
        .map(|draw| {
            (
                (draw.individual_id().as_str().to_owned(), draw.copy_index()),
                draw.allele_index(),
            )
        })
        .collect()
}

/// Falla si todas las copias recibieron el mismo alelo: con un corpus así,
/// «los genotipos no cambian» sería cierto para cualquier implementación.
fn assert_discriminating(genotypes: &BTreeMap<(String, u8), usize>) {
    let distinct = genotypes
        .values()
        .collect::<std::collections::BTreeSet<_>>();
    assert!(
        distinct.len() > 1,
        "el corpus no sortea alelos distintos y la comparación no diría nada: {genotypes:?}"
    );
}

const FOUR: [&str; 4] = ["I2", "I3", "I4", "I5"];

/// **Insertar un individuo al principio de la lista no cambia el genotipo de
/// ninguno de los demás**, comparado por identificador y no por posición.
#[test]
fn inserting_an_individual_at_the_front_leaves_every_other_genotype_alone() {
    let reference = reference();
    let before = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let after = materialize(
        &config(&["I1", "I2", "I3", "I4", "I5"], &locus("STR-1", "chr1")),
        &reference,
    );
    let before = genotypes(&before, "STR-1");
    let mut after = genotypes(&after, "STR-1");
    assert_discriminating(&before);
    after.retain(|(individual_id, _), _| individual_id != "I1");
    assert_eq!(before, after);
}

/// **Reordenar la lista de individuos tampoco cambia ningún genotipo.**
#[test]
fn reordering_the_individual_list_leaves_every_genotype_alone() {
    let reference = reference();
    let straight = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let reversed = materialize(
        &config(&["I5", "I4", "I3", "I2"], &locus("STR-1", "chr1")),
        &reference,
    );
    let straight = genotypes(&straight, "STR-1");
    assert_discriminating(&straight);
    assert_eq!(straight, genotypes(&reversed, "STR-1"));
}

/// **Añadir un locus en un cromosoma que ordena antes que los existentes no
/// cambia los genotipos de los loci previos.**
///
/// *Qué hace que pueda fallar:* el locus nuevo va en `chr0`, que ordena antes
/// que `chr1`, de modo que **desplaza el ordinal** del que ya estaba. Un ámbito
/// que nombrara el locus por su orden en la lista lo notaría; uno que lo
/// añadiera detrás, no.
#[test]
fn adding_a_locus_that_sorts_before_the_existing_ones_leaves_them_alone() {
    let reference = reference();
    let alone = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let together = materialize(
        &config(
            &FOUR,
            &format!("{}{}", locus("STR-0", "chr0"), locus("STR-1", "chr1")),
        ),
        &reference,
    );
    let alone = genotypes(&alone, "STR-1");
    assert_discriminating(&alone);
    assert_eq!(alone, genotypes(&together, "STR-1"));
    // Y el locus nuevo sí produce sorteos propios, para que la igualdad de
    // arriba no se deba a que no se sorteó nada.
    assert_eq!(genotypes(&together, "STR-0").len(), alone.len());
}

/// **Dos loci en la misma posición de cromosomas distintos sortean por
/// separado.**
///
/// *Por qué esta prueba existe, dicho sin rodeos:* la anterior **no caza** un
/// ámbito que nombrara el locus por su posición en la referencia, y se
/// comprobó mutándolo. Con subsemilla propia por locus, añadir uno no puede
/// cambiar los sorteos de otro se llame como se llame; lo que sí delata el
/// nombre posicional es que dos loci **distintos** que comparten coordenada
/// reciban exactamente los mismos genotipos. Por eso los dos loci del fixture
/// empiezan los dos en la base 21.
#[test]
fn two_loci_at_the_same_position_on_different_chromosomes_draw_independently() {
    let reference = reference();
    let output = materialize(
        &config(
            &FOUR,
            &format!("{}{}", locus("STR-0", "chr0"), locus("STR-1", "chr1")),
        ),
        &reference,
    );
    let first = genotypes(&output, "STR-0");
    let second = genotypes(&output, "STR-1");
    assert_discriminating(&first);
    assert_eq!(
        first.keys().collect::<Vec<_>>(),
        second.keys().collect::<Vec<_>>()
    );
    assert_ne!(
        first, second,
        "dos loci distintos en la misma coordenada comparten subsemilla si el ámbito los nombra por su posición"
    );
}

/// Misma semilla, mismos genotipos.
///
/// *Sin adornos:* esta prueba **no discrimina** casi nada por sí sola —repetir
/// el mismo cálculo determinista da lo mismo—. Está porque la reproducibilidad
/// es contractual, y su valor real es cazar una implementación que metiera
/// reloj, entropía del sistema u orden de recorrido en el ámbito.
#[test]
fn the_same_seed_draws_the_same_genotypes() {
    let reference = reference();
    let first = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let second = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let first = genotypes(&first, "STR-1");
    assert_discriminating(&first);
    assert_eq!(first, genotypes(&second, "STR-1"));
}

/// Con **una frecuencia de 1** para un alelo, todas las copias lo reciben.
///
/// El «no se consume aleatoriedad» que el paso también exige se comprueba
/// contando palabras, y eso solo puede hacerse desde dentro del módulo: vive en
/// la prueba unitaria `a_single_positive_weight_consumes_no_randomness` de
/// `src/repeat_genotypes.rs`.
#[test]
fn a_frequency_of_one_gives_every_copy_the_same_allele() {
    let reference = reference();
    let single = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 21
end = 80
motif = "CA"
min_longest_run_units = 3

[[reference.repeat_loci.alleles]]
units = 25
frequency = 1.0
"#;
    let output = materialize(&config(&FOUR, single), &reference);
    let drawn = genotypes(&output, "STR-1");
    assert_eq!(drawn.len(), 8, "cuatro individuos diploides");
    assert!(drawn.values().all(|index| *index == 0));
}

/// Con dos alelos al **0,9 y 0,1**, sobre mil copias y semilla congelada, la
/// proporción observada cae dentro del intervalo binomial declarado.
///
/// *Por qué 0,9 y no 0,5:* con dos alelos al cincuenta por ciento una elección
/// **uniforme que ignorase los pesos** produce exactamente la misma
/// distribución, de modo que el mutante sobreviviría. Con 0,9 no: el intervalo
/// de abajo excluye 0,5 por más de veinte desviaciones típicas.
#[test]
fn two_alleles_at_nine_to_one_stay_inside_the_binomial_interval() {
    let reference = reference();
    let skewed = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 21
end = 80
motif = "CA"
min_longest_run_units = 3

[[reference.repeat_loci.alleles]]
units = 30
frequency = 0.9

[[reference.repeat_loci.alleles]]
units = 25
frequency = 0.1
"#;
    let ids: Vec<String> = (1..=500).map(|index| format!("I{index:03}")).collect();
    let ids: Vec<&str> = ids.iter().map(String::as_str).collect();
    let output = materialize(&config(&ids, skewed), &reference);
    let drawn = genotypes(&output, "STR-1");
    assert_eq!(drawn.len(), 1000, "quinientos individuos diploides");
    let majority = drawn.values().filter(|index| **index == 0).count();
    let observed = majority as f64 / 1000.0;
    // Intervalo declarado: 0,9 ± 4 desviaciones típicas de una binomial de
    // n = 1000 y p = 0,9, es decir 0,9 ± 4·0,00949 = [0,862, 0,938]. Un reparto
    // uniforme daría 0,5, que queda a más de cuarenta desviaciones.
    let deviation = (0.9_f64 * 0.1 / 1000.0).sqrt();
    assert!(
        (observed - 0.9).abs() <= 4.0 * deviation,
        "proporción observada {observed}, fuera de 0,9 ± {}",
        4.0 * deviation
    );
}

/// **Vectores congelados del sorteo**, terna a terna.
///
/// Las pruebas de invariancia comparan una ejecución con otra, de modo que
/// **todas seguirían pasando si el sorteo entero cambiara de valores**: lo que
/// comparan es que dos ejecuciones coincidan, no cuáles son. Esta tabla fija
/// los valores concretos y es lo que caza un cambio silencioso del algoritmo,
/// del orden de los pesos o del material de la subsemilla.
#[test]
fn the_draw_has_frozen_vectors_for_each_individual_copy_and_locus() {
    let reference = reference();
    let output = materialize(&config(&FOUR, &locus("STR-1", "chr1")), &reference);
    let drawn = genotypes(&output, "STR-1");
    let expected: BTreeMap<(String, u8), usize> = FROZEN_DRAW
        .iter()
        .map(|(individual_id, copy_index, allele_index)| {
            (((*individual_id).to_owned(), *copy_index), *allele_index)
        })
        .collect();
    assert_eq!(drawn, expected);
}

/// Congelado el 22 de septiembre de 2026 con la semilla de este fichero y los
/// tres alelos al 0,5, 0,3 y 0,2. Cambiarlo exige explicar por qué.
const FROZEN_DRAW: [(&str, u8, usize); 8] = [
    ("I2", 0, 0),
    ("I2", 1, 1),
    ("I3", 0, 2),
    ("I3", 1, 0),
    ("I4", 0, 1),
    ("I4", 1, 1),
    ("I5", 0, 0),
    ("I5", 1, 1),
];

/// Los dos alelos de cada individuo en el locus, `(copia 0, copia 1)`.
fn diploid_pairs(genotypes: &BTreeMap<(String, u8), usize>) -> BTreeMap<String, (usize, usize)> {
    let mut pairs = BTreeMap::new();
    for ((individual_id, copy_index), allele) in genotypes {
        let entry = pairs
            .entry(individual_id.clone())
            .or_insert((usize::MAX, usize::MAX));
        match copy_index {
            0 => entry.0 = *allele,
            1 => entry.1 = *allele,
            other => panic!("un diploide no tiene copia {other}"),
        }
    }
    assert!(
        pairs
            .values()
            .all(|(first, second)| *first != usize::MAX && *second != usize::MAX),
        "cada individuo tiene que tener sus dos copias sorteadas"
    );
    pairs
}

fn many_ids(count: usize) -> Vec<String> {
    (1..=count).map(|index| format!("I{index:04}")).collect()
}

/// **Con F = 1 todos los individuos son homocigotos en el microsatélite**
/// (DT-058; REQ-EVO-47: «F=1 debe producir homocigosis completa»).
///
/// Es exacta: no hay intervalo. Discrimina porque el corpus es **polimórfico**
/// —se comprueba que hay más de un alelo entre los individuos— y grande: con
/// los tres alelos al 0,5, 0,3 y 0,2 sorteados copia a copia, cada individuo
/// es homocigoto con probabilidad 0,38, y que lo sean los doscientos tiene
/// probabilidad 0,38^200, del orden de 10^-84. Contra el código anterior, que
/// ignoraba F, falla.
#[test]
fn an_inbreeding_coefficient_of_one_makes_every_repeat_genotype_homozygous() {
    let reference = reference();
    let ids = many_ids(200);
    let ids: Vec<&str> = ids.iter().map(String::as_str).collect();
    let output = materialize(
        &config_with_model(
            &ids,
            &locus("STR-1", "chr1"),
            "genotype_model = \"f_diploid\"\ninbreeding_coefficient = 1.0",
        ),
        &reference,
    );
    let drawn = genotypes(&output, "STR-1");
    assert_eq!(drawn.len(), 400, "doscientos individuos diploides");
    assert_discriminating(&drawn);
    let heterozygous = diploid_pairs(&drawn)
        .into_iter()
        .filter(|(_, (first, second))| first != second)
        .collect::<Vec<_>>();
    assert!(
        heterozygous.is_empty(),
        "con F = 1 hay heterocigotos: {heterozygous:?}"
    );
}

/// **Con F = 0,8 la proporción de heterocigotos cae alrededor de
/// `(1 − F)·He`** (DT-058; REQ-EVO-48).
///
/// # El intervalo
///
/// Con los alelos al 0,5, 0,3 y 0,2, la heterocigosidad esperada es
/// `He = 1 − (0,25 + 0,09 + 0,04) = 0,62`, y con F = 0,8 la proporción de
/// heterocigotos es `Σ_{i<j} 2 (1 − F) p_i p_j = (1 − F)·He = 0,124`. Sobre
/// n = 2000 individuos independientes el recuento es binomial, con desviación
/// típica `sqrt(0,124 · 0,876 / 2000) ≈ 0,0074`. Se acepta ±5 desviaciones,
/// `[0,087, 0,161]`: la probabilidad de que un sorteo correcto caiga fuera es
/// menor que 10^-6, y la semilla es fija, así que el resultado no varía entre
/// ejecuciones.
///
/// # Lo que caza
///
/// - **Ignorar F** (el código anterior): da `He = 0,62`, a más de sesenta
///   desviaciones.
/// - **Aplicar F a una sola de las copias**, es decir, que la copia 1 no quede
///   ligada a la 0: vuelve a dar 0,62.
/// - **Olvidar el factor `1 − F`** en los pesos de la copia 1: da
///   `He/(1 + F) ≈ 0,344`, a treinta desviaciones.
///
/// Se comprueba además que **la marginal** sigue siendo la declarada: F
/// redistribuye genotipos, no frecuencias alélicas. La frecuencia del alelo
/// mayoritario sobre 4000 copias tiene, con F = 0,8, varianza
/// `p (1 − p)(1 + F)/(2n)`, desviación ≈ 0,0106; ±5 desviaciones alrededor de
/// 0,5. Esto caza un sorteo condicionado que sesgara la copia 1 hacia un alelo
/// fijo en vez de hacia el de la copia 0.
#[test]
fn an_inbreeding_coefficient_of_point_eight_lowers_repeat_heterozygosity() {
    let reference = reference();
    let ids = many_ids(2000);
    let ids: Vec<&str> = ids.iter().map(String::as_str).collect();
    let output = materialize(
        &config_with_model(
            &ids,
            &locus("STR-1", "chr1"),
            "genotype_model = \"f_diploid\"\ninbreeding_coefficient = 0.8",
        ),
        &reference,
    );
    let drawn = genotypes(&output, "STR-1");
    assert_eq!(drawn.len(), 4000, "dos mil individuos diploides");
    let pairs = diploid_pairs(&drawn);
    let heterozygous = pairs
        .values()
        .filter(|(first, second)| first != second)
        .count();
    let observed = heterozygous as f64 / 2000.0;
    let expected = (1.0 - 0.8) * 0.62;
    let deviation = (expected * (1.0 - expected) / 2000.0_f64).sqrt();
    assert!(
        (observed - expected).abs() <= 5.0 * deviation,
        "proporción de heterocigotos {observed}, fuera de {expected} ± {}",
        5.0 * deviation
    );

    let majority = drawn.values().filter(|index| **index == 0).count() as f64 / 4000.0;
    let marginal_deviation = (0.5_f64 * 0.5 * 1.8 / 4000.0).sqrt();
    assert!(
        (majority - 0.5).abs() <= 5.0 * marginal_deviation,
        "frecuencia del alelo mayoritario {majority}, fuera de 0,5 ± {}",
        5.0 * marginal_deviation
    );

    // **La fase de un heterocigoto es simétrica**: cada orientación tiene
    // probabilidad un medio, como la de los SNV al permutar su composición.
    // Sobre unos 250 heterocigotos la desviación del recuento es ≈ 7,9; ±5
    // desviaciones caza una fase fija, que pondría siempre el alelo de menor
    // índice en la copia 0.
    let ascending = pairs
        .values()
        .filter(|(first, second)| first < second)
        .count() as f64;
    let half = heterozygous as f64 / 2.0;
    assert!(
        (ascending - half).abs() <= 5.0 * (heterozygous as f64 * 0.25).sqrt(),
        "{ascending} de {heterozygous} heterocigotos llevan el alelo menor en la copia 0"
    );
}
