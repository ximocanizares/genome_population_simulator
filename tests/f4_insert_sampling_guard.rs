//! Guarda de admisión de DT-006: el sorteo del inserto de DTC-088 es un
//! rechazo que consume dos números de la cinta de azar por intento fallido y
//! repite hasta acertar. Cuando la copia haplotípica más larga de un
//! individuo hunde el intervalo admitido en la cola de la normal, el número
//! de intentos por sorteo se dispara a millones y la simulación parece
//! colgada sin ningún aviso.
//!
//! La primera versión de esta guarda predecía el coste a partir de la forma
//! de la distribución (cuántas desviaciones típicas por debajo de la media).
//! Una revisión adversarial encontró que eso mide la magnitud equivocada: lo
//! que cuesta tiempo es el número de parejas dividido por la probabilidad de
//! acertar, no la probabilidad sola, y un umbral por forma o bien rechaza
//! trabajo legítimo y barato (una sola pareja con una copia moderadamente
//! corta, resuelta en menos de cien milisegundos) o bien admite trabajo caro
//! por volumen (un panel grande con copias moderadamente cortas). Además
//! chocaba con la especificación: dos criterios de aceptación exigen que el
//! tamaño de secuencia no tenga límites codificados.
//!
//! El criterio vigente **cuenta intentos de verdad** dentro de
//! `f4_truncated_normal` (en `src/randomness.rs`) y aborta un sorteo cuando
//! agota su cota, en vez de predecir el coste desde la forma de la
//! distribución. Es exacto, no rechaza nada que no cueste de verdad, y no es
//! un límite de tamaño codificado.

use std::fs;
use std::path::PathBuf;
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::coordinates::{
    GenomicInterval, InsertionAnchor, InternalDeletion,
};
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::haplotype::{HaplotypePhase, Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::randomness::RootSeed;
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::segments::{
    DeletionModification, InsertionModification, ModificationProvenance, Segment, SegmentChain,
    SmallModification,
};
use genome_population_simulator::technical_identity::derive_technical_identity_plan;
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::{WgsPairedEndSampler, WgsSamplingError};

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";
static NONCE: AtomicU64 = AtomicU64::new(0);

fn single_sequence_reference(length: u64) -> Reference {
    let bases = vec![b'A'; length as usize];
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr1", "test").expect("sequence ID"),
        FastaHeader::parse(b"chr1").expect("FASTA header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&bases).expect("DNA"),
    )
    .expect("reference sequence");
    Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference ID"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference")
}

fn haploid_individual(reference: &Reference) -> Individual {
    Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual ID"),
        reference,
        Ploidy::new(1).expect("ploidy"),
    )
}

fn identities(
    reference: &Reference,
    individual: &Individual,
) -> genome_population_simulator::technical_identity::TechnicalIdentityPlan {
    derive_technical_identity_plan(
        RootSeed::parse("root", ROOT).expect("root seed"),
        reference.id(),
        EventId::parse("wgs-pe", "test").expect("event ID"),
        &[individual.id().clone()],
    )
    .expect("identities")
}

/// Construye un muestreador con una copia única de `copy_length_bases` y
/// pide `pair_count` parejas. Aísla exactamente la variable que decide la
/// guarda: el resto de la biblioteca es la misma en todos los casos.
fn sampler_for_copy_length(
    copy_length_bases: u64,
    pair_count: u64,
) -> Result<usize, WgsSamplingError> {
    let reference = single_sequence_reference(copy_length_bases);
    let individual = haploid_individual(&reference);
    let identities = identities(&reference, &individual);
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let sampler = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        profile,
        100,
        pair_count,
    )?;
    Ok(sampler.collect::<Result<Vec<_>, _>>()?.len())
}

#[test]
fn a_catastrophic_case_exhausts_the_attempt_budget_and_is_rejected() {
    // Oráculo independiente: medido en el arnés de la revisión adversarial,
    // una copia de 255 pb agota la cota de un millón de intentos por sorteo
    // con este mismo individuo, evento y semilla, ya en la primera o segunda
    // pareja. 255 se elige a propósito **distinto** de 250 (el mínimo del
    // perfil, que también aparece en el diagnóstico): si el productor
    // intercambiara `longest_copy_bases` por `minimum_insert_bases`, esta
    // aserción lo notaría, porque los dos números no coinciden.
    let error = sampler_for_copy_length(255, 2).expect_err("255 pb debe agotar la cota");
    assert_eq!(
        error,
        WgsSamplingError::ImpracticalInsertSampling {
            individual: IndividualId::parse("sample-A", "test").expect("individual ID"),
            longest_copy_bases: 255,
            minimum_insert_bases: 250,
        }
    );
}

#[test]
fn a_costly_but_survivable_case_is_admitted() {
    // Control positivo exigido por la revisión: 271 pb necesita del orden de
    // 23 000 intentos por sorteo --medido por la revisión adversarial--, muy
    // por debajo de la cota de un millón. Debe admitirse y completar las
    // parejas pedidas; si la guarda rechazara todo, esta aserción de éxito
    // fallaría.
    let produced = sampler_for_copy_length(271, 3).expect("271 pb debe admitirse");
    assert_eq!(produced, 3, "las tres parejas pedidas deben materializarse");
}

#[test]
fn the_known_safe_fixture_copy_length_from_the_stutter_test_remains_admitted() {
    // `tests/stutter_pipeline.rs` usa a propósito una copia de 280 pb para
    // esquivar el coste del extremo corto sin tocar el sorteo. Si la guarda
    // la rechazara, el criterio estaría mal puesto, no la prueba.
    let produced = sampler_for_copy_length(280, 2).expect("280 pb debe seguir admitiéndose");
    assert_eq!(produced, 2);
}

#[test]
fn the_guard_does_not_apply_when_no_pair_is_requested() {
    // Requisito de admisión: la guarda solo se ejercita cuando de verdad va a
    // haber lecturas. Con `pair_count = 0` el sorteo del inserto nunca se
    // ejecuta, así que incluso la copia en el mínimo del perfil debe
    // admitirse sin gastar ningún intento.
    let produced = sampler_for_copy_length(250, 0).expect("sin parejas pedidas no hay rechazo");
    assert_eq!(produced, 0);
}

#[test]
fn a_large_legitimate_panel_with_moderate_per_pair_cost_is_not_rejected_by_volume() {
    // La revisión adversarial construyó el caso que esta prueba reproduce en
    // miniatura: un panel grande (aquí, muchas parejas) sobre una copia cuyo
    // coste por sorteo es alto pero individualmente admisible (271 pb, del
    // orden de 23 000 intentos) no debe rechazarse solo por acumular trabajo.
    // El criterio cuenta intentos **por sorteo**, no un total agregado: un
    // trabajo grande y proporcionado a lo pedido no es un defecto.
    let produced = sampler_for_copy_length(271, 60).expect("el panel grande debe completarse");
    assert_eq!(produced, 60);
}

fn provenance() -> ModificationProvenance {
    ModificationProvenance::new(EventId::parse("edit", "test").expect("event ID"))
}

fn with_chain(reference: &Reference, modifications: Vec<SmallModification>) -> Individual {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let source_length = reference.sequences()[0].length();
    let chain = SegmentChain::new(vec![
        Segment::new(
            GenomicInterval::new(
                reference.id().clone(),
                sequence_id.clone(),
                0,
                source_length,
                source_length,
            )
            .expect("full interval"),
            modifications,
        )
        .expect("segment"),
    ])
    .expect("chain");
    haploid_individual(reference)
        .with_copy_segments(
            &sequence_id,
            HaplotypePhase::from_copy_index(0).expect("phase"),
            chain,
        )
        .expect("edited copy")
}

fn deletion_of(reference: &Reference, removed: u64) -> SmallModification {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let source_length = reference.sequences()[0].length();
    SmallModification::Deletion(DeletionModification::new(
        InternalDeletion::new(
            GenomicInterval::new(
                reference.id().clone(),
                sequence_id,
                1,
                1 + removed,
                source_length,
            )
            .expect("deletion interval"),
        )
        .expect("internal deletion"),
        provenance(),
    ))
}

fn insertion_of(reference: &Reference, added: usize) -> SmallModification {
    let sequence_id = SequenceId::parse("chr1", "test").expect("sequence ID");
    let source_length = reference.sequences()[0].length();
    SmallModification::Insertion(
        InsertionModification::new(
            InsertionAnchor::from_base_1(
                reference.id().clone(),
                sequence_id,
                source_length / 2,
                source_length,
            )
            .expect("anchor"),
            Arc::<[u8]>::from(vec![b'C'; added]),
            provenance(),
        )
        .expect("insertion"),
    )
}

fn build(reference: &Reference, individual: &Individual) -> Result<usize, WgsSamplingError> {
    let identities = identities(reference, individual);
    let profile = resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile");
    let sampler = WgsPairedEndSampler::new_for_pair_count(
        reference,
        individual,
        &identities,
        profile,
        100,
        2,
    )?;
    Ok(sampler.collect::<Result<Vec<_>, _>>()?.len())
}

/// Los tres casos siguientes vienen de la revisión adversarial: ¿la guarda
/// mira la copia haplotípica realizada, o la secuencia de referencia?

#[test]
fn the_guard_rejects_the_realized_copy_not_the_reference_sequence() {
    // Referencia de 300 pb (sana por sí sola) con una deleción de 45 pb: la
    // copia realizada mide 255 pb, la que de verdad importa. Si la guarda
    // mirase la referencia (300 pb), admitiría; mirando la copia, rechaza.
    let reference = single_sequence_reference(300);
    let individual = with_chain(&reference, vec![deletion_of(&reference, 45)]);
    let realized = individual.sequences()[0].copies()[0]
        .segments()
        .current_length()
        .expect("longitud realizada");
    assert_eq!(realized, 255, "la copia realizada debe medir 255 pb");
    assert!(
        matches!(
            build(&reference, &individual),
            Err(WgsSamplingError::ImpracticalInsertSampling {
                longest_copy_bases: 255,
                ..
            })
        ),
        "la guarda debe mirar la copia (255 pb), no la referencia (300 pb)"
    );
}

#[test]
fn the_guard_admits_a_copy_grown_past_the_reference_by_insertion() {
    // Referencia de 280 pb con una inserción de 200 pb: la copia mide 480 pb,
    // por encima incluso del máximo del perfil. Debe admitirse sin rechazo.
    let reference = single_sequence_reference(280);
    let individual = with_chain(&reference, vec![insertion_of(&reference, 200)]);
    let realized = individual.sequences()[0].copies()[0]
        .segments()
        .current_length()
        .expect("longitud realizada");
    assert_eq!(realized, 480, "la copia realizada debe medir 480 pb");
    assert_eq!(
        build(&reference, &individual).expect("480 pb debe admitirse"),
        2
    );
}

#[test]
fn the_guard_admits_even_when_the_reference_alone_would_have_exhausted_the_budget() {
    // El discriminador más afilado: una referencia de 255 pb es, por sí
    // sola, exactamente el caso que la primera prueba de este fichero
    // demuestra que agota la cota. Con una inserción de 225 pb la copia
    // realizada mide 480 pb, sana. Si la guarda mirase la referencia en vez
    // de la copia, rechazaría (255 pb agota la cota); mirando la copia,
    // admite. Esto es lo que de verdad distingue mirar la copia de mirar la
    // referencia: los otros dos casos de este fichero también pasarían si el
    // productor mirara la referencia por error, porque sus referencias son
    // sanas.
    let reference = single_sequence_reference(255);
    let individual = with_chain(&reference, vec![insertion_of(&reference, 225)]);
    let realized = individual.sequences()[0].copies()[0]
        .segments()
        .current_length()
        .expect("longitud realizada");
    assert_eq!(realized, 480, "la copia realizada debe medir 480 pb");
    assert_eq!(
        build(&reference, &individual).expect("copia de 480 pb debe admitirse"),
        2
    );
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f4-insert-sampling-guard-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
    ))
}

fn configuration(reference_length_bases: u64) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"
seed = "{ROOT}"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {reference_length_bases}
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[[initial_populations]]
population_id = "P1"
individual_count = 1
individual_ids = ["sample-A"]
ploidy = 1
genotype_model = "hwe"
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 5
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
    )
}

/// El diagnóstico que ve la persona usuaria debe nombrar el individuo, su
/// copia más larga, el mínimo del perfil y una salida que exista de verdad,
/// sin ningún nombre de función ni de variante de Rust.
#[test]
fn the_rejection_names_the_individual_the_lengths_and_a_real_remedy_without_rust_identifiers() {
    let root = temporary("message");
    fs::create_dir_all(&root).expect("directorio temporal");
    let config_path = root.join("general.toml");
    // 255 pb: distinto del mínimo del perfil (250), y agota la cota de
    // intentos con este mismo individuo y evento (comprobado en la primera
    // prueba de este fichero).
    fs::write(&config_path, configuration(255)).expect("config");
    let effective = load_and_validate_config(&config_path).expect("configuración efectiva");
    let destination = root.join("published");

    let error =
        execute_validated_config(&effective, &destination).expect_err("255 pb debe rechazarse");

    assert_eq!(error.code(), "GPS-DOMAIN-003");
    assert_eq!(error.category(), "DOMAIN");
    let rendered = error.render_human();
    assert!(
        rendered.contains("sample-A"),
        "nombra el individuo: {rendered}"
    );
    assert!(
        rendered.contains("255"),
        "nombra la copia más larga: {rendered}"
    );
    assert!(
        rendered.contains("250"),
        "nombra el mínimo del perfil: {rendered}"
    );
    assert!(
        rendered.contains("referencia más larga"),
        "propone una salida que existe de verdad: {rendered}"
    );
    // V1 tiene un único perfil técnico, clavado en el código: no hay un
    // perfil alternativo con inserto mínimo menor que se pueda elegir. Una
    // versión anterior de este mensaje lo ofrecía igualmente, y no era
    // cierto.
    assert!(
        !rendered.contains("inserto mínimo menor"),
        "no debe ofrecer una salida inexistente: {rendered}"
    );
    for forbidden in [
        "ImpracticalInsertSampling",
        "WgsSamplingError",
        "TruncatedNormalAttemptsExhausted",
        "RandomSamplingError",
        "wgs_paired_sampler",
        "Diagnostic",
    ] {
        assert!(
            !rendered.contains(forbidden),
            "no debe nombrar identificadores de Rust ({forbidden}): {rendered}"
        );
    }

    fs::remove_dir_all(&root).ok();
}
