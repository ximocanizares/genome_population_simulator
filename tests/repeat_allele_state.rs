//! F6.4.3: el alelo de repetición dentro del estado haplotípico.
//!
//! # Cómo están elegidos los fixtures
//!
//! La regla, escrita antes de elegir ninguno: **qué tiene que ser distinto en
//! el fixture para que la aserción pueda fallar**. Aquí lo que se afirma es que
//! el estado guarda la *secuencia* del tracto y no un contador, de modo que:
//!
//! - los flancos de la referencia son **distintos del tracto y no vacíos**, y
//!   distintos entre sí: una materialización que devolviera la referencia, o
//!   que se comiera un flanco, se nota;
//! - la pareja de alelos que demuestra que un contador no basta tiene **la
//!   misma longitud, el mismo número de copias completas y la misma pureza**, y
//!   solo cambia **dónde** está la interrupción. Con cualquier otra pareja, una
//!   representación por contador seguiría distinguiéndolos y la prueba no diría
//!   nada.

use genome_population_simulator::coordinates::{BasePosition, GenomicInterval};
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{
    EventId, IndividualId, ReferenceId, RepeatLocusId, SequenceId,
};
use genome_population_simulator::materialization::{materialize_copy, materialize_copy_interval};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    CanonicalBase, CurrentRepeatAllelePiece, ModificationProvenance, SegmentError, SnpModification,
};
use genome_population_simulator::tandem_repeat::{RepeatMotif, TandemRepeatTract};

/// Flanco izquierdo, tracto de la referencia y flanco derecho.
///
/// Los dos flancos son distintos entre sí y ninguno contiene `CA`: si la
/// materialización se desplazara una base, o devolviera la referencia en vez
/// del tracto, el resultado dejaría de coincidir byte a byte.
const LEFT: &[u8] = b"GGGTTT";
const REFERENCE_TRACT: &[u8] = b"CACACACACACA"; // (CA)6, doce bases
const RIGHT: &[u8] = b"TTTGGG";

fn reference() -> Reference {
    let mut bases = Vec::new();
    bases.extend_from_slice(LEFT);
    bases.extend_from_slice(REFERENCE_TRACT);
    bases.extend_from_slice(RIGHT);
    Reference::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        ReferenceOrigin::Generated,
        vec![
            ReferenceSequence::new(
                SequenceId::parse("chr1", "test").expect("sequence id"),
                FastaHeader::parse(b"chr1").expect("header"),
                SequenceClass::Chromosome,
                &DnaSequence::try_from_bytes(&bases).expect("DNA"),
            )
            .expect("reference sequence"),
        ],
    )
    .expect("reference")
}

fn tract_interval() -> GenomicInterval {
    let total = (LEFT.len() + REFERENCE_TRACT.len() + RIGHT.len()) as u64;
    GenomicInterval::new(
        ReferenceId::parse("ref-main", "test").expect("reference id"),
        SequenceId::parse("chr1", "test").expect("sequence id"),
        LEFT.len() as u64,
        (LEFT.len() + REFERENCE_TRACT.len()) as u64,
        total,
    )
    .expect("interval")
}

fn allele(bases: &[u8]) -> CurrentRepeatAllelePiece {
    CurrentRepeatAllelePiece::new(
        RepeatLocusId::parse("STR-1", "test").expect("locus id"),
        tract_interval(),
        TandemRepeatTract::new(RepeatMotif::new(b"CA").expect("motif"), bases).expect("tract"),
        ModificationProvenance::random(EventId::parse("STR-1", "test").expect("event")),
    )
}

fn individual() -> Individual {
    Individual::from_reference(
        IndividualId::parse("I1", "test").expect("individual id"),
        &reference(),
        Ploidy::DIPLOID,
    )
}

fn phase(index: u8) -> HaplotypePhase {
    HaplotypePhase::from_copy_index(index).expect("phase")
}

fn expected(tract: &[u8]) -> Vec<u8> {
    let mut bases = Vec::new();
    bases.extend_from_slice(LEFT);
    bases.extend_from_slice(tract);
    bases.extend_from_slice(RIGHT);
    bases
}

fn with_allele(individual: &Individual, copy_index: u8, bases: &[u8]) -> Individual {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copy = individual
        .copies_for(&sequence_id)
        .expect("copies")
        .iter()
        .find(|copy| copy.phase() == phase(copy_index))
        .expect("copy")
        .clone();
    let segments = copy
        .segments()
        .with_current_repeat_allele(allele(bases))
        .expect("the interval arrives intact");
    individual
        .with_copy_segments(&sequence_id, phase(copy_index), segments)
        .expect("copy replaced")
}

/// Las dos copias de un diploide con alelos de **5 y 7 unidades** materializan
/// secuencias distintas, y cada una es el flanco izquierdo, su tracto y el
/// flanco derecho, byte a byte.
///
/// *Qué hace que pueda fallar:* los flancos no son vacíos y no contienen el
/// motivo, de modo que una materialización que devolviera las bases de la
/// referencia para las dos copias —la variante defectuosa que el paso nombra—
/// produciría dos veces el tracto de seis unidades y las dos aserciones
/// morirían.
#[test]
fn two_copies_with_different_repeat_alleles_materialize_different_sequences() {
    let reference = reference();
    let mut individual = individual();
    individual = with_allele(&individual, 0, b"CACACACACA"); // (CA)5
    individual = with_allele(&individual, 1, b"CACACACACACACA"); // (CA)7

    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copies = individual.copies_for(&sequence_id).expect("copies");
    let first = materialize_copy(&copies[0], &reference).expect("copy 0");
    let second = materialize_copy(&copies[1], &reference).expect("copy 1");

    assert_eq!(first, expected(b"CACACACACA"));
    assert_eq!(second, expected(b"CACACACACACACA"));
    assert_ne!(first, second);
    assert_ne!(first, expected(REFERENCE_TRACT));
}

/// **Dos alelos con el mismo número de unidades y la interrupción en sitio
/// distinto materializan secuencias distintas y se distinguen en el estado.**
///
/// *Qué hace que pueda fallar:* los dos alelos miden **doce bases**, tienen
/// **cinco copias completas** de `CA` y por tanto **la misma pureza, 10/12**, y
/// la misma tirada más larga, **cuatro**. Lo único que cambia es dónde está el
/// `TT`: en uno parte el tracto tras la primera copia y en el otro tras la
/// cuarta, y en los dos la tirada mayor mide cuatro.
/// Una representación que guardara el número de unidades —la variante
/// defectuosa que este paso nombra— los colapsaría en el mismo objeto y las dos
/// aserciones morirían. Una pareja de alelos de distinta longitud **no
/// discriminaría nada**, porque un contador también los separa.
#[test]
fn two_alleles_of_equal_length_with_the_interruption_elsewhere_stay_distinct() {
    let reference = reference();
    let left_interruption = b"CATTCACACACA";
    let right_interruption = b"CACACACATTCA";
    assert_eq!(left_interruption.len(), right_interruption.len());

    let mut individual = individual();
    individual = with_allele(&individual, 0, left_interruption);
    individual = with_allele(&individual, 1, right_interruption);

    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copies = individual.copies_for(&sequence_id).expect("copies");
    assert_ne!(copies[0].segments(), copies[1].segments());
    let first = materialize_copy(&copies[0], &reference).expect("copy 0");
    let second = materialize_copy(&copies[1], &reference).expect("copy 1");
    assert_eq!(first, expected(left_interruption));
    assert_eq!(second, expected(right_interruption));
    assert_ne!(first, second);

    // Y lo derivado coincide en todo menos en la posición de la interrupción,
    // que es justo lo que hace la comparación anterior imprescindible.
    let annotation = |bases: &[u8]| {
        TandemRepeatTract::new(RepeatMotif::new(b"CA").expect("motif"), bases)
            .expect("tract")
            .annotation()
    };
    let left = annotation(left_interruption);
    let right = annotation(right_interruption);
    assert_eq!(left.units_as_f64(), right.units_as_f64());
    assert_eq!(left.purity_as_f64(), right.purity_as_f64());
    assert_eq!(left.complete_copies(), right.complete_copies());
    // La tirada más larga también coincide: es un dato que la docstring afirma,
    // de modo que se comprueba en vez de creerse.
    assert_eq!(left.longest_run_units(), 4);
    assert_eq!(right.longest_run_units(), 4);
    assert_ne!(left.interruptions(), right.interruptions());
}

/// Dos alelos sobre el **mismo tracto** de la misma copia se rechazan.
///
/// *Qué hace que pueda fallar:* el control positivo de la línea anterior —la
/// primera aplicación, que sí se admite— es lo que impide que esta prueba
/// pasara con una operación que rechazara siempre.
#[test]
fn a_second_repeat_allele_on_the_same_tract_is_rejected() {
    let individual = with_allele(&individual(), 0, b"CACACACACA");
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copy = individual
        .copies_for(&sequence_id)
        .expect("copies")
        .iter()
        .find(|copy| copy.phase() == phase(0))
        .expect("copy")
        .clone();
    assert_eq!(
        copy.segments()
            .with_current_repeat_allele(allele(b"CACACACACACACA")),
        Err(SegmentError::RepeatIntervalNotIntact)
    );
}

/// Un cambio de base **dentro** del intervalo del locus impide el alelo, y uno
/// **fuera** no.
///
/// *Qué hace que pueda fallar:* la pareja. Sin el caso de fuera, una guarda que
/// rechazara en cuanto la copia tuviera cualquier SNP pasaría igual. Las dos
/// posiciones se eligen pegadas al borde del intervalo —la última base del
/// flanco izquierdo y la primera del tracto— para que lo único que las separe
/// sea el borde.
#[test]
fn a_base_change_inside_the_interval_blocks_the_allele_and_one_outside_does_not() {
    let reference = reference();
    let sequence_length = reference.sequences()[0].length();
    let snp_at = |offset: u64, from: CanonicalBase, to: CanonicalBase| {
        SnpModification::new(
            BasePosition::new(
                ReferenceId::parse("ref-main", "test").expect("reference id"),
                SequenceId::parse("chr1", "test").expect("sequence id"),
                offset,
                sequence_length,
            )
            .expect("position"),
            from,
            to,
            ModificationProvenance::new(EventId::parse("snp-1", "test").expect("event")),
        )
        .expect("snp")
    };
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copy = individual()
        .copies_for(&sequence_id)
        .expect("copies")
        .iter()
        .find(|copy| copy.phase() == phase(0))
        .expect("copy")
        .clone();

    // Última base del flanco izquierdo: `T` en la posición 5. Está fuera.
    let outside = copy
        .segments()
        .with_snp(snp_at(5, CanonicalBase::T, CanonicalBase::A))
        .expect("snp outside the tract");
    outside
        .with_current_repeat_allele(allele(b"CACACACACA"))
        .expect("un cambio fuera del intervalo no estorba");

    // Primera base del tracto: `C` en la posición 6. Está dentro.
    let inside = copy
        .segments()
        .with_snp(snp_at(6, CanonicalBase::C, CanonicalBase::A))
        .expect("snp inside the tract");
    assert_eq!(
        inside.with_current_repeat_allele(allele(b"CACACACACA")),
        Err(SegmentError::RepeatIntervalAlreadyModified)
    );
}

/// **Un lector independiente reconstruye la secuencia y, de ella, el número de
/// unidades y la pureza.**
///
/// El lector recorta el tracto de la copia materializada por sus coordenadas y
/// lo compara con **constantes escritas a mano**, no con otra llamada al
/// módulo: doce bases, cinco copias completas de `CA` —porque `TT` rompe la
/// cuarta—, pureza 10/12 y tirada más larga de tres.
#[test]
fn an_independent_reader_recovers_the_units_and_the_purity_from_the_bases() {
    let reference = reference();
    let tract = b"CACACATTCACA";
    let individual = with_allele(&individual(), 0, tract);
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copy = individual
        .copies_for(&sequence_id)
        .expect("copies")
        .iter()
        .find(|copy| copy.phase() == phase(0))
        .expect("copy")
        .clone();
    let bases = materialize_copy(&copy, &reference).expect("materialized copy");
    let start = LEFT.len();
    let end = start + tract.len();
    assert_eq!(&bases[start..end], tract);

    let annotation =
        TandemRepeatTract::new(RepeatMotif::new(b"CA").expect("motif"), &bases[start..end])
            .expect("tract")
            .annotation();
    assert_eq!(annotation.units_as_f64(), 6.0);
    assert_eq!(annotation.complete_copies(), 5);
    assert_eq!(annotation.purity_numerator(), 10);
    assert!((annotation.purity_as_f64() - 10.0 / 12.0).abs() < 1e-12);
    assert_eq!(annotation.longest_run_units(), 3);
}

/// La reconstrucción por ventanas ve el mismo tracto que la reconstrucción
/// completa.
///
/// *Por qué hace falta:* las lecturas no materializan la copia entera, sino un
/// intervalo. Son dos recorridos distintos del mismo estado, y sin esta prueba
/// uno podría emitir el tracto y el otro las bases de la referencia.
#[test]
fn the_windowed_reconstruction_agrees_with_the_complete_one() {
    let reference = reference();
    let tract = b"CACACATTCACACA"; // siete unidades, una interrupción
    let individual = with_allele(&individual(), 0, tract);
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence id");
    let copy = individual
        .copies_for(&sequence_id)
        .expect("copies")
        .iter()
        .find(|copy| copy.phase() == phase(0))
        .expect("copy")
        .clone();
    let complete = materialize_copy(&copy, &reference).expect("complete");
    assert_eq!(complete, expected(tract));
    let window =
        materialize_copy_interval(&copy, &reference, 0, complete.len() as u64).expect("window");
    assert_eq!(window, complete);
    // Y una ventana que empieza dentro del flanco y acaba dentro del flanco
    // derecho tiene que llevar el tracto entero, no el de la referencia.
    let partial = materialize_copy_interval(
        &copy,
        &reference,
        (LEFT.len() - 2) as u64,
        (LEFT.len() + tract.len() + 2) as u64,
    )
    .expect("partial window");
    let mut expected_partial = LEFT[LEFT.len() - 2..].to_vec();
    expected_partial.extend_from_slice(tract);
    expected_partial.extend_from_slice(&RIGHT[..2]);
    assert_eq!(partial, expected_partial);
}
