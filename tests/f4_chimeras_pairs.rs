//! Verdad manual de F4.3.3: enumera las segundas fuentes sin reutilizar el
//! selector compacto del productor y reconstruye la orientación PE esperada.

use std::num::NonZeroU64;
use std::path::Path;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::chimera::{ChimeraError, ChimeraSettings, transform_paired_end};
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::haplotype::{Individual, Ploidy};
use genome_population_simulator::identity::{EventId, IndividualId, ReferenceId, SequenceId};
use genome_population_simulator::paired_end::{PairedEndFragment, emit_paired_end};
use genome_population_simulator::randomness::{
    F4TechnicalDecision, RootSeed, derive_f4_technical_subseed, f4_uniform_index,
};
use genome_population_simulator::reference::{
    DnaSequence, FastaHeader, Reference, ReferenceOrigin, ReferenceSequence, SequenceClass,
};
use genome_population_simulator::technical_identity::{
    F4TechnicalIdentityDeriver, TemplateOrdinal, derive_technical_identity_plan,
};
use genome_population_simulator::technical_profile::resolve_f4_technical_profile;
use genome_population_simulator::wgs_sampling::WgsPairedEndSampler;

const ROOT: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'N' => b'N',
            _ => panic!("fixture DNA"),
        })
        .collect()
}

fn fixture(length: usize, ploidy: u8) -> (Reference, Individual, EventId) {
    let sequences = [
        ("chr-a", b"ACGT".repeat(length / 4 + 1)),
        ("chr-b", b"TGCA".repeat(length / 4 + 1)),
    ]
    .into_iter()
    .map(|(id, mut bases)| {
        bases.truncate(length);
        ReferenceSequence::new(
            SequenceId::parse(id, "test").expect("sequence"),
            FastaHeader::parse(id.as_bytes()).expect("header"),
            SequenceClass::Chromosome,
            &DnaSequence::try_from_bytes(&bases).expect("DNA"),
        )
        .expect("reference sequence")
    })
    .collect();
    let reference = Reference::new(
        ReferenceId::parse("reference-main", "test").expect("reference"),
        ReferenceOrigin::Generated,
        sequences,
    )
    .expect("reference");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-A", "test").expect("individual"),
        &reference,
        Ploidy::new(ploidy).expect("ploidy"),
    );
    (
        reference,
        individual,
        EventId::parse("chimera-event", "test").expect("event"),
    )
}

fn pair(
    reference: &Reference,
    individual: &Individual,
    event: &EventId,
) -> genome_population_simulator::paired_end::PairedEndPair {
    let root = RootSeed::parse("root", ROOT).expect("root");
    let identities = derive_technical_identity_plan(
        root,
        reference.id(),
        event.clone(),
        &[individual.id().clone()],
    )
    .expect("identities");
    WgsPairedEndSampler::new_for_pair_count(
        reference,
        individual,
        &identities,
        resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile"),
        150,
        1,
    )
    .expect("sampler")
    .next()
    .expect("one pair")
    .expect("pair")
}

fn pair_without_second_source() -> (
    Reference,
    Individual,
    EventId,
    genome_population_simulator::paired_end::PairedEndPair,
) {
    let sequence = ReferenceSequence::new(
        SequenceId::parse("chr-only", "test").expect("sequence"),
        FastaHeader::parse(b"chr-only").expect("header"),
        SequenceClass::Chromosome,
        &DnaSequence::try_from_bytes(&[b'A'; 250]).expect("DNA"),
    )
    .expect("reference sequence");
    let reference = Reference::new(
        ReferenceId::parse("reference-only", "test").expect("reference"),
        ReferenceOrigin::Generated,
        vec![sequence],
    )
    .expect("reference");
    let individual = Individual::from_reference(
        IndividualId::parse("sample-only", "test").expect("individual"),
        &reference,
        Ploidy::new(1).expect("ploidy"),
    );
    let event = EventId::parse("chimera-only", "test").expect("event");
    let root = RootSeed::parse("root", ROOT).expect("root");
    let identities = derive_technical_identity_plan(
        root,
        reference.id(),
        event.clone(),
        &[individual.id().clone()],
    )
    .expect("identities");
    let template = identities
        .derive_template(individual.id(), TemplateOrdinal::new(0))
        .expect("template");
    let copy = &individual.sequences()[0].copies()[0];
    let fragment = PairedEndFragment::new(template.clone(), copy.id().clone(), 0, vec![b'A'; 250])
        .expect("fragment");
    let pre_hop = F4TechnicalIdentityDeriver::new(root, event.clone())
        .derive_own_pre_hop_pair(template)
        .expect("pre-hop identity");
    let pair = emit_paired_end(fragment, pre_hop, 150).expect("pair");
    (reference, individual, event, pair)
}

#[test]
fn chimera_uses_a_manual_canonical_second_source_list_and_orients_pe_after_composition() {
    let (reference, individual, event) = fixture(500, 2);
    let root = RootSeed::parse("root", ROOT).expect("root");
    let pair = pair(&reference, &individual, &event);
    let fragment = pair.fragment();
    let ordinal = fragment.template().template_ordinal().get();
    let library = fragment.template().library_id().clone();
    let insert = fragment.sequence().len() as u64;
    let mut breakpoint_rng = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::chimera_breakpoint(&event, &library, ordinal),
    )
    .expect("breakpoint seed")
    .local_rng();
    let breakpoint = f4_uniform_index(&mut breakpoint_rng, insert - 1).expect("breakpoint") + 1;
    let suffix = insert - breakpoint;
    let first_start = fragment.start();
    let first_end = first_start + breakpoint;
    let first_copy = fragment.haplotype_copy_id().clone();
    let first_copy_ref = &first_copy;
    let first_sequence = fragment.sequence().to_vec();

    // Oráculo manual: materializa la lista, en vez de los rangos compactos del
    // productor. La exclusión incluye los dos intervalos que solo se tocan.
    let mut candidates = individual
        .sequences()
        .iter()
        .flat_map(|genotype| genotype.copies().iter())
        .flat_map(|copy| {
            let length = copy.segments().current_length().expect("copy length");
            let maximum_start = length.saturating_sub(suffix);
            (length >= suffix).then_some((copy, maximum_start))
        })
        .flat_map(|(copy, maximum_start)| {
            (0..=maximum_start).filter_map(move |start| {
                let end = start + suffix;
                (copy.id() != first_copy_ref || end < first_start || start > first_end)
                    .then_some((copy, start))
            })
        })
        .collect::<Vec<_>>();
    candidates.sort_by(|(left_copy, left_start), (right_copy, right_start)| {
        left_copy
            .id()
            .cmp(right_copy.id())
            .then(left_start.cmp(right_start))
    });
    assert!(
        !candidates.is_empty(),
        "the fixture has a manual candidate list"
    );
    let mut source_rng = derive_f4_technical_subseed(
        root,
        F4TechnicalDecision::chimera_second_source(
            &event,
            &library,
            ordinal,
            NonZeroU64::new(breakpoint).expect("positive breakpoint"),
        ),
    )
    .expect("source seed")
    .local_rng();
    let expected = candidates
        [f4_uniform_index(&mut source_rng, candidates.len() as u64).expect("source") as usize];

    let transformed = transform_paired_end(
        root,
        &event,
        &reference,
        &individual,
        pair,
        ChimeraSettings::new(1.0).expect("rate"),
    )
    .expect("chimera");
    let origin = transformed.origin().expect("transient provenance");
    assert_eq!(origin.breakpoint().get(), breakpoint);
    assert_eq!(origin.first().copy_id(), &first_copy);
    assert_eq!(origin.first().start(), first_start);
    assert_eq!(origin.first().end(), first_end);
    assert_eq!(origin.second().copy_id(), expected.0.id());
    assert_eq!(origin.second().start(), expected.1);
    assert_eq!(origin.second().end(), expected.1 + suffix);

    let expected_sequence = [
        &first_sequence[..breakpoint as usize],
        &reference
            .sequences()
            .iter()
            .find(|sequence| sequence.id() == expected.0.id().sequence_id())
            .expect("source sequence")
            .sequence()
            .as_bytes()
            .expect("bases materializadas")[expected.1 as usize..(expected.1 + suffix) as usize],
    ]
    .concat();
    assert_eq!(transformed.pair().fragment().sequence(), expected_sequence);
    assert_eq!(
        transformed.pair().r1().sequence(),
        &expected_sequence[..150]
    );
    assert_eq!(
        transformed.pair().r2().sequence(),
        reverse_complement(&expected_sequence[expected_sequence.len() - 150..])
    );
}

#[test]
fn selected_chimera_without_a_second_source_fails_and_default_rate_keeps_pairs_clean() {
    let (reference, individual, event, clean) = pair_without_second_source();
    let root = RootSeed::parse("root", ROOT).expect("root");
    let neutral = transform_paired_end(
        root,
        &event,
        &reference,
        &individual,
        clean.clone(),
        ChimeraSettings::default(),
    )
    .expect("zero rate does not need candidates");
    assert!(!neutral.transformed());
    assert_eq!(neutral.into_pair(), clean);

    assert!(matches!(
        transform_paired_end(
            root,
            &event,
            &reference,
            &individual,
            clean,
            ChimeraSettings::new(1.0).expect("rate"),
        ),
        Err(ChimeraError::NoSecondSource)
    ));
}

#[test]
fn sampler_keeps_pair_cardinality_and_exposes_only_an_aggregate_chimera_count() {
    let (reference, individual, event) = fixture(500, 2);
    let root = RootSeed::parse("root", ROOT).expect("root");
    let identities =
        derive_technical_identity_plan(root, reference.id(), event, &[individual.id().clone()])
            .expect("identities");
    let mut sampler = WgsPairedEndSampler::new_for_pair_count(
        &reference,
        &individual,
        &identities,
        resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1").expect("profile"),
        150,
        4,
    )
    .expect("sampler")
    .with_effective_technical_event(&effective_event_with_chimera_rate(1.0));
    let pairs = (&mut sampler)
        .collect::<Result<Vec<_>, _>>()
        .expect("all pairs");
    assert_eq!(pairs.len(), 4);
    assert_eq!(sampler.chimera_report().realized_pairs(), 4);
}

#[test]
fn configured_chimera_rate_is_defaulted_validated_and_reaches_the_pair_sampler() {
    let defaulted = technical_event_toml("0").replace("chimera_rate = 0\n", "");
    let defaulted = parse_declared_config(
        admit_toml_header(
            Path::new("default-chimera-rate.toml"),
            defaulted.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared");
    assert_eq!(
        build_effective_config(&defaulted)
            .expect("effective")
            .technical_event()
            .expect("event")
            .chimera_rate(),
        0.0
    );
    for rate in ["-0.1", "1.1", "nan"] {
        let text = technical_event_toml(rate);
        let declared = parse_declared_config(
            admit_toml_header(
                Path::new("invalid-chimera-rate.toml"),
                text.into_bytes(),
                ExpectedConfigKind::General,
            )
            .expect("header"),
        )
        .expect("declared");
        assert_eq!(
            build_effective_config(&declared)
                .expect_err("invalid rate")
                .code(),
            "GPS-CONF-008"
        );
    }
}

fn effective_event_with_chimera_rate(
    rate: f64,
) -> genome_population_simulator::application::EffectiveTechnicalEvent {
    let text = technical_event_toml(&rate.to_string());
    let declared = parse_declared_config(
        admit_toml_header(
            Path::new("chimera-rate.toml"),
            text.into_bytes(),
            ExpectedConfigKind::General,
        )
        .expect("header"),
    )
    .expect("declared");
    build_effective_config(&declared)
        .expect("effective")
        .technical_event()
        .expect("event")
        .clone()
}

fn technical_event_toml(chimera_rate: &str) -> String {
    format!(
        r#"config_kind = "general"
schema_version = "2"

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
length = 1000

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "chimera-event"
individual_ids = ["sample-A"]
target_coverage_wgs = 1
chimera_rate = {chimera_rate}

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
    )
}

#[test]
fn initial_profile_exposes_the_zero_chimera_rate() {
    assert_eq!(
        resolve_f4_technical_profile("gps-hiseqx-pcrfree-pe150-v1")
            .expect("profile")
            .chimera_rate(),
        0.0
    );
}

/// DT-031: una referencia de 600 bases, una sola copia y una tasa de quimeras
/// de 1 pueden dejar a una pareja sin segundo origen (DTC-084). Debe salir como
/// diagnóstico de dominio (`GPS-DOMAIN-003`, salida 1, el mismo del precedente
/// DT-006) que nombra al individuo, y no como invariante interna imposible
/// (`GPS-INTERNAL-001`, salida 70), sin reintento ni bajada silenciosa de la
/// tasa.
///
/// **Depende de la semilla.** Con esta configuración falla solo para algunas
/// semillas (el revisor midió que 19 de 20 publican): la prueba usa `ROOT`
/// porque con ella el sorteo falla, y no demuestra que 600 pb con una copia y
/// tasa 1 nunca quepan. El control positivo es la misma configuración con tasa
/// 0, que publica: lo que la rechaza es la tasa y no el resto del escenario.
#[test]
fn a_chimera_rate_that_does_not_fit_a_short_reference_is_a_domain_diagnostic() {
    use genome_population_simulator::application::load_and_validate_config;
    use genome_population_simulator::execution::execute_validated_config;

    let toml = |chimera_rate: &str| {
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
length = 600

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
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "chimera-event"
individual_ids = ["sample-A"]
target_coverage_wgs = 0.5
read_length_bases = 2
chimera_rate = {chimera_rate}

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"#
        )
    };
    let root = std::env::temp_dir().join(format!("gps-dt031-{}", std::process::id()));
    std::fs::create_dir_all(&root).expect("temporary directory");

    let control_path = root.join("control.toml");
    std::fs::write(&control_path, toml("0")).expect("control config");
    let control = load_and_validate_config(&control_path).expect("control config validates");
    execute_validated_config(&control, &root.join("control-result"))
        .expect("without chimeras the scenario publishes");

    let path = root.join("rate-one.toml");
    std::fs::write(&path, toml("1")).expect("config");
    let effective = load_and_validate_config(&path).expect("the configuration is admissible");
    let error = execute_validated_config(&effective, &root.join("result"))
        .expect_err("with this seed the chimera draw finds a pair without a second source");

    assert_eq!(error.code(), "GPS-DOMAIN-003", "{error:?}");
    assert_eq!(error.category(), "DOMAIN");
    assert_eq!(error.exit_code(), 1);
    let rendered = error.render_human();
    assert!(rendered.contains("sample-A"), "{rendered}");
    assert!(rendered.contains("segundo_origen_de_quimera"), "{rendered}");
    assert!(rendered.contains("referencia más larga"), "{rendered}");
    assert!(rendered.contains("semilla"), "{rendered}");
    for forbidden in ["NoSecondSource", "ChimeraError", "INTERNAL"] {
        assert!(!rendered.contains(forbidden), "{forbidden}: {rendered}");
    }
    std::fs::remove_dir_all(&root).ok();
}
