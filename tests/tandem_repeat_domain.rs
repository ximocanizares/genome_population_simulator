//! F6.4.1 y F6.4.1.1: el modelo puro del locus de repetición y el paso de
//! unidades sobre la tirada ininterrumpida más larga.
//!
//! Todas las anotaciones se comprueban contra un recuento hecho a mano en el
//! propio fixture, no contra otra llamada al modelo: un oráculo que se calcule
//! con la función que valida no discrimina nada.

use genome_population_simulator::tandem_repeat::{
    RepeatMotif, TandemRepeatError, TandemRepeatTract, TractBlock,
};

fn motif(text: &str) -> RepeatMotif {
    RepeatMotif::new(text.as_bytes()).expect("motivo canónico")
}

fn tract(motif_text: &str, tract_text: &str) -> TandemRepeatTract {
    TandemRepeatTract::new(motif(motif_text), tract_text.as_bytes()).expect("tracto canónico")
}

/// Tracto puro: quince bases, cinco unidades, pureza 1, sin interrupciones y
/// tirada más larga de cinco.
///
/// *Cuidado:* este fixture es puro y **no discrimina** casi nada. Con un tracto
/// puro coinciden la pureza contada sobre bases y sobre unidades, el número de
/// unidades entero y el redondeado, y el modelo de múltiplos exactos con el
/// fraccionario. Está aquí como caso de referencia, no como guarda.
#[test]
fn pure_acg_tract_annotates_five_units_unit_purity_and_a_run_of_five() {
    let annotation = tract("ACG", "ACGACGACGACGACG").annotation();
    assert_eq!(annotation.total_bases(), 15);
    assert_eq!(annotation.exact_units(), Some((15, 3)));
    assert!((annotation.units_as_f64() - 5.0).abs() < f64::EPSILON);
    assert_eq!(annotation.purity_numerator(), 15);
    assert_eq!(annotation.purity_denominator(), 15);
    assert!(annotation.interruptions().is_empty());
    assert_eq!(annotation.longest_run_units(), 5);
}

/// `AC AC AC GG AC AC`: doce bases, **seis** unidades, pureza **10/12** sobre
/// bases, una interrupción y tirada más larga **3**, no 6.
///
/// Es el fixture que separa la tirada ininterrumpida más larga del número
/// **total** de unidades: las dos cifras son 3 y 6, y un modelo que devuelva la
/// segunda muere aquí.
#[test]
fn impure_ac_tract_separates_the_longest_run_from_the_total_unit_count() {
    let annotation = tract("AC", "ACACACGGACAC").annotation();
    assert_eq!(annotation.total_bases(), 12);
    assert_eq!(annotation.exact_units(), Some((12, 2)));
    assert!((annotation.units_as_f64() - 6.0).abs() < f64::EPSILON);
    assert_eq!(annotation.complete_copies(), 5);
    assert_eq!(annotation.purity_numerator(), 10);
    assert_eq!(annotation.purity_denominator(), 12);
    assert_eq!(annotation.interruptions().len(), 1);
    assert_eq!(annotation.interruptions()[0].start(), 6);
    assert_eq!(annotation.interruptions()[0].length(), 2);
    assert_eq!(
        annotation.longest_run_units(),
        3,
        "la tirada más larga es 3; el total de unidades, 6"
    );
}

/// El ejemplo del propio estándar VCF 4.5, `(CAG)5(CA)1(CAG)4`: veintinueve
/// bases, 29/3 unidades, pureza 27/29 y tirada más larga 5.
///
/// Este fixture es el que mata el modelo de múltiplos exactos: 29 no es
/// múltiplo de 3 y el tracto **tiene que admitirse**.
#[test]
fn the_vcf_standard_cag_example_is_admitted_with_a_fractional_unit_count() {
    let text = "CAGCAGCAGCAGCAG".to_owned() + "CA" + "CAGCAGCAGCAG";
    assert_eq!(text.len(), 29);
    let annotation = tract("CAG", &text).annotation();
    assert_eq!(annotation.total_bases(), 29);
    assert_eq!(annotation.exact_units(), Some((29, 3)));
    assert!((annotation.units_as_f64() - 29.0 / 3.0).abs() < 1e-12);
    assert!(
        annotation.units_as_f64().fract() > 0.0,
        "un tracto con unidad parcial no puede anotarse como entero"
    );
    assert_eq!(annotation.complete_copies(), 9);
    assert_eq!(annotation.purity_numerator(), 27);
    assert_eq!(annotation.purity_denominator(), 29);
    assert_eq!(annotation.longest_run_units(), 5);
}

/// `CACACATCACA`: interrupción de **una sola base**, 5,5 unidades.
///
/// Es el fixture que separa la pureza contada sobre bases de la contada sobre
/// unidades: sobre unidades daría 5/6, sobre bases da 10/11. También mata el
/// número de unidades redondeado a entero, porque 5,5 no es 5 ni 6.
#[test]
fn a_single_base_interruption_gives_five_and_a_half_units_and_purity_over_bases() {
    let annotation = tract("CA", "CACACATCACA").annotation();
    assert_eq!(annotation.total_bases(), 11);
    assert!((annotation.units_as_f64() - 5.5).abs() < 1e-12);
    assert_eq!(annotation.complete_copies(), 5);
    assert_eq!(
        (
            annotation.purity_numerator(),
            annotation.purity_denominator()
        ),
        (10, 11),
        "sobre bases 10/11; sobre unidades habría dado 5/6"
    );
    assert!((annotation.purity_as_f64() - 10.0 / 11.0).abs() < 1e-12);
    assert_eq!(annotation.interruptions().len(), 1);
    assert_eq!(annotation.interruptions()[0].start(), 6);
    assert_eq!(annotation.interruptions()[0].length(), 1);
}

/// **La dirección del barrido voraz queda fijada aquí, y solo aquí.**
///
/// Todos los demás fixtures usan motivos **no solapantes consigo mismos**, y con
/// ellos el barrido de derecha a izquierda da los mismos números y las mismas
/// posiciones de interrupción: **no discriminan la dirección**. Hace falta un
/// motivo con período interno. Con el homopolímero `AA` sobre `AAA` hay una
/// única copia completa y la base sobrante puede quedar a un lado o al otro: por
/// la izquierda la interrupción cae en el desplazamiento **2**; por la derecha
/// caería en el **0**.
///
/// Importa porque los homopolímeros son la clase de repetición más común, y
/// porque la posición publicada de la interrupción viaja a la verdad y al VCF.
#[test]
fn the_greedy_scan_runs_left_to_right_and_a_periodic_motif_proves_it() {
    let annotation = tract("AA", "AAA").annotation();
    assert_eq!(annotation.total_bases(), 3);
    assert_eq!(annotation.complete_copies(), 1);
    assert_eq!(annotation.interruptions().len(), 1);
    assert_eq!(
        annotation.interruptions()[0].start(),
        2,
        "de izquierda a derecha la base sobrante queda al final; de derecha a izquierda quedaría al principio"
    );
    assert_eq!(annotation.runs()[0].start(), 0);

    // Un caso algo mayor, para que no dependa de un tracto de tres bases: con
    // cinco `A` y motivo `AA` hay dos copias desde la izquierda y la sobrante
    // queda en el desplazamiento 4.
    let longer = tract("AA", "AAAAA").annotation();
    assert_eq!(longer.complete_copies(), 2);
    assert_eq!(longer.interruptions()[0].start(), 4);
    assert_eq!(longer.longest_run_units(), 2);
}

/// Motivo vacío y motivo con `N` se rechazan en el constructor.
#[test]
fn an_empty_motif_and_a_motif_with_n_are_rejected() {
    assert_eq!(RepeatMotif::new(b""), Err(TandemRepeatError::EmptyMotif));
    assert_eq!(
        RepeatMotif::new(b"CAN"),
        Err(TandemRepeatError::NonCanonicalMotifBase { base: b'N' })
    );
    assert!(RepeatMotif::new(b"CAG").is_ok());
}

/// Dos alelos con el **mismo número de unidades** y la interrupción en sitio
/// distinto producen secuencias distintas y anotaciones distintas.
///
/// Es el argumento biológico de la decisión S-06: en genética de poblaciones y
/// en forense esos dos alelos son alelos distintos, y un contador de unidades
/// los colapsaría en el mismo objeto.
#[test]
fn equal_length_alleles_with_the_interruption_elsewhere_are_distinguishable() {
    let left = tract("CA", "CACACATCACA");
    let right = tract("CA", "CACATCACACA");
    assert_eq!(left.length(), right.length());
    let left_annotation = left.annotation();
    let right_annotation = right.annotation();
    assert_eq!(
        left_annotation.units_as_f64(),
        right_annotation.units_as_f64(),
        "por longitud son indistinguibles: esa es la premisa del caso"
    );
    assert_eq!(
        left_annotation.purity_numerator(),
        right_annotation.purity_numerator()
    );
    assert_ne!(left.bases(), right.bases());
    assert_ne!(
        left_annotation.interruptions()[0].start(),
        right_annotation.interruptions()[0].start()
    );
    assert_ne!(left_annotation, right_annotation);
    assert_ne!(left, right);
}

/// La anotación se recalcula desde las bases después de un paso: no se hereda.
///
/// El recuento esperado está escrito a mano en el propio fixture. La variante
/// defectuosa que caza es la que guarda la anotación junto a la secuencia y la
/// propaga al aplicar un paso; con ella, el tracto de veintinueve bases seguiría
/// anunciando veintiséis y una tirada de cinco.
#[test]
fn the_annotation_is_recomputed_from_the_bytes_after_a_step() {
    let before = tract("CAG", "CAGCAGCAGCAGCAGTTCAGCAGCAG");
    let before_annotation = before.annotation();
    assert_eq!(before_annotation.total_bases(), 26);
    assert_eq!(before_annotation.longest_run_units(), 5);
    assert_eq!(before_annotation.complete_copies(), 8);

    let after = before.apply_unit_step(1).expect("expansión de una unidad");
    let after_annotation = after.annotation();
    assert_eq!(after.length(), 29);
    assert_eq!(after_annotation.total_bases(), 29);
    assert_eq!(after_annotation.longest_run_units(), 6);
    assert_eq!(after_annotation.complete_copies(), 9);
    assert_eq!(after_annotation.purity_numerator(), 27);
    assert_eq!(after_annotation.purity_denominator(), 29);
}

/// Quitar una unidad de un tracto puro de diez deja las nueve restantes
/// idénticas y la diferencia es una deleción contigua de una unidad.
///
/// *Cuidado:* con un tracto **puro** todas las convenciones de dónde resbalar
/// coinciden y esta prueba **no discrimina** entre ellas. Lo único que
/// comprueba es que el resultado es una deleción limpia.
#[test]
fn removing_one_unit_from_a_pure_tract_is_a_clean_contiguous_deletion() {
    let ten = TandemRepeatTract::pure(&motif("GATA"), 10).expect("tracto puro");
    let nine = ten.apply_unit_step(-1).expect("contracción de una unidad");
    assert_eq!(nine.length(), ten.length() - 4);
    assert_eq!(nine.bases(), &ten.bases()[..36]);
    assert_eq!(
        nine,
        TandemRepeatTract::pure(&motif("GATA"), 9).expect("tracto puro")
    );
}

/// Con tiradas de **5 y 3**, y la más larga **a la izquierda**, el paso cae en
/// la de 5 y la de 3 queda byte a byte intacta.
///
/// Es la única prueba que caza la convención del extremo derecho que la
/// decisión S-07 retiró: bajo esa convención el paso caería en la tirada de la
/// derecha y el sufijo cambiaría.
#[test]
fn the_step_lands_on_the_longest_run_even_when_it_is_not_the_rightmost() {
    let before = tract("CAG", "CAGCAGCAGCAGCAGTTCAGCAGCAG");
    let runs = before.annotation();
    assert_eq!(runs.runs().len(), 2);
    assert_eq!(runs.runs()[0].units(), 5);
    assert_eq!(runs.runs()[1].units(), 3);

    let after = before.apply_unit_step(-1).expect("contracción");
    assert_eq!(
        after.as_text(),
        "CAGCAGCAGCAGTTCAGCAGCAG",
        "la unidad se pierde en la tirada de 5"
    );
    // *Retirada:* aquí había una comprobación de que la cola queda byte a byte
    // intacta. Estaba **implicada** por el `assert_eq!` del texto completo que
    // la precede, de modo que no podía fallar por su cuenta.
    assert_eq!(after.annotation().runs()[1].units(), 3);
}

/// Con tiradas **empatadas** de 4 y 4, el paso cae en la de la izquierda.
///
/// Es la única prueba de la batería que ejercita el desempate. Con cualquier
/// otro fixture la tirada ganadora es única y el criterio no se toca.
#[test]
fn tied_runs_break_towards_the_left_and_the_result_is_reproducible() {
    let before = tract("CA", "CACACACATCACACACA");
    let annotation = before.annotation();
    assert_eq!(annotation.runs().len(), 2);
    assert_eq!(annotation.runs()[0].units(), 4);
    assert_eq!(annotation.runs()[1].units(), 4);

    let after = before.apply_unit_step(1).expect("expansión");
    assert_eq!(
        after.as_text(),
        "CACACACACATCACACACA",
        "la unidad se gana en la tirada de la izquierda"
    );
    // *Retirada:* aquí había un `assert_eq!(after, before.apply_unit_step(1))`
    // presentado como «el desempate es determinista». **No podía fallar**: es
    // llamar dos veces a una función pura sin estado. Lo que sí discrimina el
    // desempate es la comparación anterior contra el texto literal, porque un
    // desempate hacia la derecha produce otro texto.
    assert_eq!(after.annotation().runs()[0].units(), 5);
    assert_eq!(after.annotation().runs()[1].units(), 4);
}

/// Aplicar `+1` y luego `−1` devuelve la secuencia original byte a byte, tanto
/// con tiradas desiguales como con tiradas empatadas.
#[test]
fn plus_one_then_minus_one_restores_the_original_bytes() {
    for text in ["CAGCAGCAGCAGCAGTTCAGCAGCAG", "CACACACATCACACACA"] {
        let motif_text = if text.starts_with("CAG") { "CAG" } else { "CA" };
        let original = tract(motif_text, text);
        let roundtrip = original
            .apply_unit_step(1)
            .expect("expansión")
            .apply_unit_step(-1)
            .expect("contracción");
        assert_eq!(
            roundtrip.bases(),
            original.bases(),
            "ida y vuelta en {text}"
        );
    }
}

/// Un paso que vaciaría la tirada más larga se rechaza, y uno de cero también.
#[test]
fn a_step_that_would_empty_the_longest_run_is_rejected() {
    let single = tract("CA", "CATT");
    assert_eq!(single.annotation().longest_run_units(), 1);
    assert_eq!(
        single.apply_unit_step(-1),
        Err(TandemRepeatError::StepEmptiesLongestRun {
            run_units: 1,
            requested: 1,
        })
    );
    assert!(single.apply_unit_step(1).is_ok());
    assert_eq!(single.apply_unit_step(0), Err(TandemRepeatError::ZeroStep));
}

// ---------------------------------------------------------------------------
// Decisión del responsable del 20/09/2026: un locus puede tener varios motivos.
// ---------------------------------------------------------------------------

fn composite() -> TandemRepeatTract {
    // `(TCTA)5(TCTG)5(TCTA)5`, la forma de muchos marcadores forenses.
    TandemRepeatTract::pure_blocks(&[(motif("TCTA"), 5), (motif("TCTG"), 5), (motif("TCTA"), 5)])
        .expect("microsatélite compuesto")
}

/// **El fixture que un modelo de un solo motivo no puede pasar.**
///
/// Las mismas sesenta bases leídas con un solo motivo `TCTA` dan pureza 40/60,
/// porque el bloque central no es ese motivo; leídas con sus tres bloques dan
/// pureza 1, que es la verdad biológica. Con el mínimo por omisión de 0,8 la
/// primera lectura se rechaza y la segunda se admite: no es un matiz numérico,
/// es la diferencia entre poder declarar el locus y no poder.
#[test]
fn a_composite_microsatellite_is_pure_by_blocks_and_impure_read_as_one_motif() {
    let by_blocks = composite().annotation();
    assert_eq!(by_blocks.total_bases(), 60);
    assert_eq!(by_blocks.blocks().len(), 3);
    assert_eq!(by_blocks.purity_numerator(), 60);
    assert_eq!(by_blocks.purity_denominator(), 60);
    assert!((by_blocks.purity_as_f64() - 1.0).abs() < f64::EPSILON);

    let as_one = tract("TCTA", &composite().as_text()).annotation();
    assert_eq!(as_one.total_bases(), 60);
    assert_eq!(as_one.blocks().len(), 1);
    assert_eq!(
        (as_one.purity_numerator(), as_one.purity_denominator()),
        (40, 60),
        "leído con un solo motivo, el bloque central cuenta entero como interrupción"
    );
    assert!(
        as_one.purity_as_f64() < 0.8 && by_blocks.purity_as_f64() >= 0.8,
        "el mínimo por omisión separa las dos lecturas"
    );
}

/// El número total de unidades es la suma de las de cada bloque, y las
/// magnitudes por bloque son las que el VCF publica en listas paralelas.
#[test]
fn the_block_annotation_carries_the_parallel_lists_the_vcf_defines() {
    let annotation = composite().annotation();
    assert!((annotation.units_as_f64() - 15.0).abs() < f64::EPSILON);
    let motifs: Vec<String> = annotation
        .blocks()
        .iter()
        .map(|block| block.motif().as_text())
        .collect();
    assert_eq!(motifs, vec!["TCTA", "TCTG", "TCTA"]);
    for block in annotation.blocks() {
        // La correspondencia entre las listas: copias por longitud del motivo
        // son las bases del bloque. Si una lista se desalinea respecto de otra,
        // esta igualdad deja de cumplirse.
        assert!(
            (block.units_as_f64() * block.motif().length() as f64 - block.length() as f64).abs()
                < f64::EPSILON,
            "bloque descuadrado: {block:?}"
        );
    }
    let bases: u64 = annotation.blocks().iter().map(|block| block.length()).sum();
    assert_eq!(bases, annotation.total_bases());
    assert_eq!(annotation.blocks()[0].start(), 0);
    assert_eq!(annotation.blocks()[1].start(), 20);
    assert_eq!(annotation.blocks()[2].start(), 40);
}

/// **La tirada más larga gana sea del motivo que sea.**
///
/// Con bloques de 3, 7 y 2 unidades, la tirada más larga es la del **segundo**
/// motivo. Un modelo que mirase solo el primer bloque, o que prefiriese el
/// motivo declarado en primer lugar, diría 3.
#[test]
fn the_longest_run_wins_across_blocks_whatever_its_motif() {
    let tract = TandemRepeatTract::pure_blocks(&[
        (motif("TCTA"), 3),
        (motif("TCTG"), 7),
        (motif("TCTA"), 2),
    ])
    .expect("tracto compuesto");
    let annotation = tract.annotation();
    assert_eq!(annotation.longest_run_units(), 7);
    let run = annotation.longest_run().expect("hay tirada");
    assert_eq!(
        run.block_index(),
        1,
        "la tirada más larga es la del segundo motivo"
    );
    assert_eq!(run.start(), 12);
}

/// **El paso gana o pierde una unidad del bloque donde cae el resbalón**, y los
/// demás bloques quedan byte a byte intactos.
#[test]
fn the_step_lands_on_the_block_of_the_longest_run() {
    let before = TandemRepeatTract::pure_blocks(&[
        (motif("TCTA"), 3),
        (motif("TCTG"), 7),
        (motif("TCTA"), 2),
    ])
    .expect("tracto compuesto");
    let after = before.apply_unit_step(1).expect("expansión");

    assert_eq!(after.length(), before.length() + 4);
    assert_eq!(
        after.as_text(),
        "TCTA".repeat(3) + &"TCTG".repeat(8) + &"TCTA".repeat(2),
        "la unidad ganada es una del motivo del bloque donde cayó el resbalón"
    );
    let blocks: Vec<(String, u64)> = after
        .blocks()
        .iter()
        .map(|block| (block.motif().as_text(), block.length()))
        .collect();
    assert_eq!(
        blocks,
        vec![
            ("TCTA".to_owned(), 12),
            ("TCTG".to_owned(), 32),
            ("TCTA".to_owned(), 8)
        ],
        "solo cambia de longitud el bloque afectado"
    );
    assert_eq!(&after.bases()[..12], &before.bases()[..12]);
    assert_eq!(&after.bases()[44..], &before.bases()[40..]);

    // Y la vuelta devuelve el original byte a byte.
    assert_eq!(
        after.apply_unit_step(-1).expect("contracción").bases(),
        before.bases()
    );
}

/// La pureza se mide sobre **bases** y sobre el tracto completo, sumando las
/// copias completas de todos los bloques.
#[test]
fn purity_sums_the_complete_copies_of_every_block_over_the_whole_tract() {
    // Bloque puro de 5 unidades de `TCTA` y bloque impuro de 12 bases de `CA`
    // con solo 4 copias completas.
    let blocks = vec![
        TractBlock::new(motif("TCTA"), 20).expect("bloque"),
        TractBlock::new(motif("CA"), 12).expect("bloque"),
    ];
    let bases = ("TCTA".repeat(5) + "CACACACAGGGG").into_bytes();
    let tract = TandemRepeatTract::from_blocks(blocks, bases).expect("tracto");
    let annotation = tract.annotation();
    assert_eq!(annotation.total_bases(), 32);
    assert_eq!(annotation.blocks()[0].complete_copies(), 5);
    assert_eq!(annotation.blocks()[1].complete_copies(), 4);
    assert_eq!(
        (
            annotation.purity_numerator(),
            annotation.purity_denominator()
        ),
        (28, 32),
        "cinco copias de cuatro bases más cuatro de dos, sobre treinta y dos bases"
    );
    // Y la interrupción del segundo bloque se publica con desplazamiento
    // absoluto en el tracto, no relativo a su bloque.
    let interruptions = annotation.interruptions();
    assert_eq!(interruptions.len(), 1);
    assert_eq!(interruptions[0].start(), 28);
    assert_eq!(interruptions[0].length(), 4);
}

/// Los bloques tienen que cubrir el tracto exactamente.
#[test]
fn blocks_that_do_not_cover_the_tract_are_rejected() {
    let blocks = vec![TractBlock::new(motif("CAG"), 9).expect("bloque")];
    assert_eq!(
        TandemRepeatTract::from_blocks(blocks, b"CAGCAG".to_vec()),
        Err(TandemRepeatError::BlockLengthsDoNotCoverTract {
            declared: 9,
            observed: 6,
        })
    );
    assert_eq!(
        TandemRepeatTract::from_blocks(Vec::new(), b"CAG".to_vec()),
        Err(TandemRepeatError::NoBlocks)
    );
}

/// Un tracto de un solo bloque se comporta exactamente como antes de que
/// hubiera bloques: mismas unidades, misma pureza, mismo motivo accesible.
#[test]
fn a_single_block_tract_behaves_exactly_as_before() {
    let single = tract("AC", "ACACACGGACAC");
    let annotation = single.annotation();
    assert_eq!(annotation.blocks().len(), 1);
    assert_eq!(single.motif().expect("un solo motivo").as_text(), "AC");
    assert!((annotation.units_as_f64() - 6.0).abs() < f64::EPSILON);
    assert_eq!(annotation.purity_numerator(), 10);
    assert_eq!(annotation.longest_run_units(), 3);
    assert!(
        composite().motif().is_none(),
        "con varios bloques no existe «el motivo» del tracto"
    );
}

/// **El desempate de la tirada más larga ENTRE bloques.**
///
/// *Añadido tras la segunda revisión adversarial.* Los empates **dentro** de un
/// bloque estaban probados; los de **entre** bloques, que son los que la
/// decisión de varios motivos estrena, no los cubría nadie: recorriendo los
/// bloques al revés sobrevivían las 54 pruebas de los dos ficheros, y el alelo
/// resultante **sí cambia**.
///
/// Con dos bloques empatados a tres unidades y motivos distintos, gana el de la
/// izquierda. Un recorrido inverso ganaría el de la derecha y produciría otro
/// alelo; como el generador de la entrega 5 va a consumir esta operación, el no
/// determinismo se propagaría.
#[test]
fn a_tie_between_blocks_breaks_towards_the_left_block() {
    let before = TandemRepeatTract::pure_blocks(&[(motif("TCTA"), 3), (motif("TCTG"), 3)])
        .expect("dos bloques empatados");
    let annotation = before.annotation();
    assert_eq!(annotation.blocks()[0].runs()[0].units(), 3);
    assert_eq!(annotation.blocks()[1].runs()[0].units(), 3);
    let run = annotation.longest_run().expect("hay tirada");
    assert_eq!(
        run.block_index(),
        0,
        "el empate lo gana el bloque de la izquierda"
    );
    assert_eq!(run.start(), 0);

    let after = before.apply_unit_step(1).expect("expansión");
    assert_eq!(
        after.as_text(),
        "TCTA".repeat(4) + &"TCTG".repeat(3),
        "la unidad ganada es del motivo del bloque izquierdo"
    );
    let lengths: Vec<u64> = after.blocks().iter().map(|block| block.length()).collect();
    assert_eq!(lengths, vec![16, 12]);
    // Y con el desempate hacia la derecha el alelo sería este otro, que es la
    // prueba de que los dos resultados son distinguibles.
    assert_ne!(after.as_text(), "TCTA".repeat(3) + &"TCTG".repeat(4));
}

/// **Una interrupción que cruza la frontera entre dos bloques es una, no dos.**
///
/// *Añadido tras la segunda revisión adversarial.* El tipo está documentado como
/// «bases contiguas que ninguna copia cubre», lo que implica tramos maximales, y
/// con varios bloques la contabilidad por bloque partía el tramo en el punto
/// donde el usuario decidió cortar la estructura. El número de interrupciones de
/// un microsatélite tiene significado biológico y no puede depender de eso.
#[test]
fn an_interruption_crossing_a_block_boundary_counts_as_one() {
    let blocks = vec![
        TractBlock::new(motif("CA"), 8).expect("bloque"),
        TractBlock::new(motif("GT"), 8).expect("bloque"),
    ];
    // `CACA` + `TTTT` | `TT` + `GTGTGT`: el tramo sin cubrir va de la 4 a la 10
    // y cruza la frontera, que está en la 8.
    let tract = TandemRepeatTract::from_blocks(blocks, b"CACATTTTTTGTGTGT".to_vec())
        .expect("tracto de dos bloques");
    let annotation = tract.annotation();

    let merged = annotation.interruptions();
    assert_eq!(merged.len(), 1, "la interrupción es una sola: {merged:?}");
    assert_eq!(merged[0].start(), 4);
    assert_eq!(merged[0].length(), 6);

    // Y la contabilidad de cada bloque sigue siendo la suya, cortada en la
    // frontera, que es lo que su documentación dice.
    assert_eq!(annotation.blocks()[0].interruptions().len(), 1);
    assert_eq!(annotation.blocks()[0].interruptions()[0].length(), 4);
    assert_eq!(annotation.blocks()[1].interruptions().len(), 1);
    assert_eq!(annotation.blocks()[1].interruptions()[0].start(), 8);
    assert_eq!(annotation.blocks()[1].interruptions()[0].length(), 2);

    // Dos interrupciones que NO se tocan siguen siendo dos.
    let separated = TandemRepeatTract::from_blocks(
        vec![
            TractBlock::new(motif("CA"), 8).expect("bloque"),
            TractBlock::new(motif("GT"), 8).expect("bloque"),
        ],
        b"CACATTCAGTGTTTGT".to_vec(),
    )
    .expect("tracto de dos bloques");
    assert!(separated.annotation().interruptions().len() >= 2);
}

/// La fracción exacta de unidades solo existe con un bloque, y con varios se
/// dice que no existe en vez de devolver un número sin significado.
#[test]
fn the_exact_unit_fraction_is_absent_when_there_are_several_blocks() {
    let single = tract("CAG", "CAGCAGCAGCAGCAG");
    assert_eq!(single.annotation().exact_units(), Some((15, 3)));
    assert_eq!(composite().annotation().exact_units(), None);
}
