//! F5.1.3. Marginales `p_i` por bloque×subpoblación desde el espectro.
//!
//! Estas pruebas recorren la frontera completa —configuración, admisión
//! sellada, layout y marginales— y contrastan cada valor con un oráculo entero
//! independiente que no reutiliza el generador productivo.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveF5Config, load_and_validate_f5_config};
use genome_population_simulator::f5_layout::{
    F5LayoutRequest, F5PhysicalLayout, build_f5_physical_layout,
};
use genome_population_simulator::f5_marginals::{
    F5MarginalSpectrum, F5MarginalTable, F5SpectrumBand, build_f5_marginals,
};
use genome_population_simulator::f5_preflight::{
    F5OtherResources, SealedF5Admission, recover_sealed_f5_admission,
};
use genome_population_simulator::identity::{
    ExecutionId, F5BlockId, F5SiteId, PopulationId, ScenarioId,
};
use genome_population_simulator::randomness::{F5RandomDecision, RootSeed, derive_f5_subseed};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);

const F5_RNG_CATALOG_TSV: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");
const TWO_POW_53: u64 = 1_u64 << 53;
const EXECUTION: &str = "RUN_AQWGL4DQ5HJMZYWT3U4EDNLTJ4";
const SEED: &str = "a1a2a3a4a5a6a7a8b1b2b3b4b5b6b7b8c1c2c3c4c5c6c7c8d1d2d3d4d5d6d7d8";

/// Espectro del fixture: un intervalo degenerado en cero, uno continuo y otro
/// degenerado en uno. Los pesos son diádicos exactos para que el oráculo
/// entero no dependa de ninguna tolerancia.
const SPECTRUM_TSV: &str =
    "alt_af_min\talt_af_max\tweight\n0\t0\t0.25\n0.25\t0.75\t0.5\n1\t1\t0.25\n";
const SPECTRUM_NUMERATORS: [u64; 3] = [1, 2, 1];
const SPECTRUM_DENOMINATOR: u64 = 4;

fn sha256(bytes: &[u8]) -> String {
    let mut hasher = Sha256::new();
    hasher.update(bytes);
    format!("{:x}", hasher.finalize())
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f5-marginals-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Escribe el paquete de perfil y los dos TOML de autoría con las
/// subpoblaciones indicadas, en el orden literal recibido.
fn write_fixture(
    root: &Path,
    populations: &[&str],
    block_count: u64,
    haplotype_count: usize,
) -> PathBuf {
    let scenario_directory = root.join("authoring/scenarios");
    let profile_directory = scenario_directory.join("profiles");
    fs::create_dir_all(&profile_directory).expect("profile directory");
    fs::create_dir_all(root.join("authoring/references")).expect("reference directory");

    let mut ld_rows = String::from(
        "population_profile_id\tdistance_start_bp\tdistance_end_bp\ttarget_mean_r2\tmax_abs_error\n",
    );
    let mut sorted = populations.to_vec();
    sorted.sort_unstable();
    for population in &sorted {
        ld_rows.push_str(&format!("{population}\t1\t2\t0.5\t0.05\n"));
    }

    let files = [
        (
            "block_lengths.tsv",
            "f5_block_length_profile",
            "length_min_bp\tlength_max_bp\tweight\n3\t3\t1\n".to_owned(),
        ),
        (
            "snp_gaps.tsv",
            "f5_snp_gap_profile",
            "gap_min_bp\tgap_max_bp\tweight\n1\t1\t1\n".to_owned(),
        ),
        (
            "snp_spectrum.tsv",
            "f5_snp_spectrum",
            SPECTRUM_TSV.to_owned(),
        ),
        ("haplotype_weights.tsv", "f5_haplotype_weights", {
            let weight = 1.0 / haplotype_count as f64;
            let mut rows = String::from("haplotype_rank\tweight\n");
            for rank in 0..haplotype_count {
                rows.push_str(&format!("{rank}\t{weight}\n"));
            }
            rows
        }),
        ("ld_profile.tsv", "f5_ld_profile", ld_rows),
    ];
    let links = files
        .iter()
        .map(|(name, kind, contents)| {
            fs::write(profile_directory.join(name), contents).expect("profile file");
            format!(
                r#"{{"artifact_kind":"{kind}","schema_version":"1","path":"{name}","size_bytes":{},"sha256":"{}"}}"#,
                contents.len(),
                sha256(contents.as_bytes())
            )
        })
        .collect::<Vec<_>>()
        .join(",");
    let metadata = format!(
        r#"{{"profile_kind":"f5_generation_profile","schema_version":"1","profile_id":"profile","provenance_kind":"user_supplied","source_citation":"test","source_license":null,"source_sha256":null,"files":[{links}]}}"#
    );
    fs::write(profile_directory.join("profile.metadata.json"), &metadata).expect("metadata");

    let fasta = ">chr1\nACGTAC\n";
    fs::write(root.join("authoring/references/ref.fa"), fasta).expect("reference");
    let general = root.join("authoring/general.toml");
    fs::write(
        &general,
        format!(
            "config_kind = \"general\"\nschema_version = \"3\"\nseed = \"{SEED}\"\nploidy = 2\nreference_id = \"ref\"\nreference_path = \"references/ref.fa\"\nreference_sha256 = \"{}\"\n\n[f5_scenario]\nscenario_id = \"scenario\"\nscenario_path = \"scenarios/scenario.toml\"\n",
            sha256(fasta.as_bytes())
        ),
    )
    .expect("general");

    let mut subpopulations = String::new();
    for (index, population) in populations.iter().enumerate() {
        subpopulations.push_str(&format!(
            "\n[[subpopulations]]\npopulation_id = \"{population}\"\nindividual_count = 1\nfounder_ids = [\"I{}\"]\n",
            index + 1
        ));
    }
    fs::write(
        scenario_directory.join("scenario.toml"),
        format!(
            "config_kind = \"f5_scenario\"\nschema_version = \"1\"\nscenario_id = \"scenario\"\nfounder_count = {}\n\n[f5]\nprofile_metadata_path = \"profiles/profile.metadata.json\"\nprofile_metadata_sha256 = \"{}\"\nmax_layout_attempts = 8\nmax_haplotype_attempts = 1\nmax_swap_proposals = 1\nmax_haplotype_cells = 64\nblock_count = {block_count}\n{subpopulations}",
            populations.len(),
            sha256(metadata.as_bytes())
        ),
    )
    .expect("scenario");
    general
}

fn build(
    root: &Path,
    populations: &[&str],
) -> (
    EffectiveF5Config,
    SealedF5Admission,
    ExecutionId,
    F5PhysicalLayout,
) {
    build_with_blocks(root, populations, 2)
}

fn build_with_blocks(
    root: &Path,
    populations: &[&str],
    block_count: u64,
) -> (
    EffectiveF5Config,
    SealedF5Admission,
    ExecutionId,
    F5PhysicalLayout,
) {
    let general = write_fixture(root, populations, block_count, 2);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud de layout"),
    )
    .expect("layout físico");
    (config, sealed, execution, layout)
}

fn fixture_spectrum() -> F5MarginalSpectrum {
    F5MarginalSpectrum::new(
        &ScenarioId::parse("scenario", "test").expect("scenario"),
        vec![
            F5SpectrumBand::new(0.0, 0.0, 0.25),
            F5SpectrumBand::new(0.25, 0.75, 0.5),
            F5SpectrumBand::new(1.0, 1.0, 0.25),
        ],
    )
    .expect("espectro del fixture")
}

/// Oráculo independiente del espectro del fixture.
///
/// No usa `f64` ni el generador productivo: resuelve la fila con aritmética
/// entera exacta sobre los pesos diádicos y, para la fila continua
/// `[1/4, 3/4]`, redondea el racional exacto `(2^52+k)/2^54` a binary64 con
/// ties-to-even explícito.
fn oracle_marginal(words: &mut impl Iterator<Item = u64>) -> f64 {
    let row_bits = words.next().expect("palabra de fila") >> 11;
    let mut prefix = 0_u128;
    let mut row = SPECTRUM_NUMERATORS.len() - 1;
    for (index, numerator) in SPECTRUM_NUMERATORS.iter().enumerate() {
        prefix += u128::from(*numerator);
        if u128::from(row_bits) * u128::from(SPECTRUM_DENOMINATOR) < prefix * u128::from(TWO_POW_53)
        {
            row = index;
            break;
        }
    }
    match row {
        0 => 0.0,
        2 => 1.0,
        _ => {
            let value_bits = u128::from(words.next().expect("palabra de valor") >> 11);
            let numerator = u128::from(1_u64 << 52) + value_bits;
            if numerator < u128::from(TWO_POW_53) {
                // Cabe en 53 bits: el racional ya es un binary64 exacto.
                numerator as f64 / (1_u128 << 54) as f64
            } else {
                let mut half = numerator / 2;
                if numerator % 2 == 1 && half % 2 == 1 {
                    // Empate exacto a media unidad: ties-to-even sube al par.
                    half += 1;
                }
                half as f64 / (1_u128 << 53) as f64
            }
        }
    }
}

fn oracle_site(
    root_seed: RootSeed,
    scenario_id: &ScenarioId,
    population_id: &PopulationId,
    block_id: &F5BlockId,
    site_id: &F5SiteId,
) -> f64 {
    let mut rng = derive_f5_subseed(
        root_seed,
        F5RandomDecision::marginal_spectrum(scenario_id, population_id, block_id, site_id),
    )
    .expect("subsemilla del ámbito marginal_spectrum")
    .local_rng();
    let mut words = std::iter::from_fn(|| Some(rng.next_u64()));
    oracle_marginal(&mut words)
}

fn table_values(table: &F5MarginalTable) -> Vec<(String, u64, u64, f64)> {
    table
        .blocks()
        .iter()
        .flat_map(|block| {
            block.sites().iter().map(move |site| {
                (
                    block.population_id().as_str().to_owned(),
                    block.block_index(),
                    site.site_index(),
                    site.target_alt_af(),
                )
            })
        })
        .collect()
}

#[test]
fn marginal_spectrum_reproduces_the_published_catalog_vector_and_its_first_two_words() {
    let scenario_id = ScenarioId::parse("scenario-f5", "test").expect("scenario");
    let population_id = PopulationId::parse("population-A", "test").expect("population");
    let block_id = F5BlockId::parse("BLK_blockA", "test").expect("block");
    let site_id = F5SiteId::parse("SITE_siteA", "test").expect("site");
    let root_seed = RootSeed::from_bytes(std::array::from_fn(|index| index as u8));

    let row = F5_RNG_CATALOG_TSV
        .lines()
        .find(|line| line.starts_with("marginal_spectrum|"))
        .expect("la línea publicada del ámbito marginal_spectrum");
    let fields = row.split('|').collect::<Vec<_>>();
    assert_eq!(fields.len(), 9);
    assert_eq!(
        fields[1],
        "scenario_id=scenario-f5;population_id=population-A;block_id=BLK_blockA;site_id=SITE_siteA;stage=marginal_spectrum"
    );

    let decision =
        F5RandomDecision::marginal_spectrum(&scenario_id, &population_id, &block_id, &site_id);
    let subseed = derive_f5_subseed(root_seed, decision).expect("subsemilla");
    assert_eq!(subseed.to_hex(), fields[4], "subsemilla publicada");
    let mut rng = subseed.local_rng();
    let words: [u64; 4] = std::array::from_fn(|_| rng.next_u64());
    for (index, word) in words.iter().enumerate() {
        assert_eq!(
            format!("{word:016x}"),
            fields[index + 5],
            "palabra {index} publicada"
        );
    }

    // El vector se cierra hasta el marginal: las dos primeras palabras del
    // ámbito publicado, leídas por el oráculo entero, dan el mismo `p_i` que
    // el generador.
    let expected = oracle_marginal(&mut words.into_iter());
    let observed = fixture_spectrum()
        .sample_site(root_seed, &scenario_id, &population_id, &block_id, &site_id)
        .expect("marginal del ámbito publicado");
    assert_eq!(observed, expected);
}

#[test]
fn the_table_covers_every_block_and_site_of_every_subpopulation_in_canonical_order() {
    let root = temporary("canonical-order");
    let (config, config_sealed, config_execution, layout) = build(&root, &["P2", "P1"]);
    let table = build_f5_marginals(&config, &config_sealed, &config_execution, &layout)
        .expect("marginales");

    assert_eq!(layout.blocks().len(), 2);
    assert!(layout.blocks().iter().all(|block| block.sites().len() == 3));
    assert_eq!(table.blocks().len(), 4);
    assert_eq!(
        table
            .blocks()
            .iter()
            .map(|block| (
                block.population_id().as_str().to_owned(),
                block.block_index()
            ))
            .collect::<Vec<_>>(),
        vec![
            ("P1".to_owned(), 0),
            ("P1".to_owned(), 1),
            ("P2".to_owned(), 0),
            ("P2".to_owned(), 1),
        ],
        "el orden canónico es subpoblación y después bloque, no el del escenario"
    );

    let scenario_id = ScenarioId::parse(config.scenario_id(), "test").expect("scenario");
    let root_seed = RootSeed::parse("seed", config.seed()).expect("semilla");
    assert_eq!(
        F5MarginalSpectrum::from_profile(&scenario_id, config.profile()).expect("espectro leído"),
        fixture_spectrum(),
        "el oráculo entero describe exactamente el espectro que el perfil aporta"
    );
    for block in table.blocks() {
        let layout_block = layout
            .blocks()
            .iter()
            .find(|candidate| candidate.block_index() == block.block_index())
            .expect("bloque del layout");
        assert_eq!(block.block_id(), layout_block.block_id());
        assert_eq!(block.sites().len(), layout_block.sites().len());
        for (site, layout_site) in block.sites().iter().zip(layout_block.sites()) {
            assert_eq!(site.site_index(), layout_site.site_index());
            assert_eq!(site.site_id(), layout_site.site_id());
            assert_eq!(
                site.target_alt_af(),
                oracle_site(
                    root_seed,
                    &scenario_id,
                    block.population_id(),
                    block.block_id(),
                    site.site_id(),
                ),
                "{}/{}/{}",
                block.population_id(),
                block.block_index(),
                site.site_index()
            );
            assert!((0.0..=1.0).contains(&site.target_alt_af()));
        }
    }

    // El espectro del fixture solo admite 0, 1 o el interior de [1/4, 3/4]:
    // ningún valor puede caer fuera de esas tres regiones.
    for (_, _, _, value) in table_values(&table) {
        assert!(
            value == 0.0 || value == 1.0 || (0.25..=0.75).contains(&value),
            "{value} no pertenece a ninguna fila del espectro"
        );
    }
    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn reordering_or_adding_a_subpopulation_leaves_the_other_marginals_untouched() {
    let first_root = temporary("declared-order");
    let (first_config, first_config_sealed, first_config_execution, first_layout) =
        build(&first_root, &["P1", "P2"]);
    let first = build_f5_marginals(
        &first_config,
        &first_config_sealed,
        &first_config_execution,
        &first_layout,
    )
    .expect("marginales declarados");

    let swapped_root = temporary("swapped-order");
    let (swapped_config, swapped_config_sealed, swapped_config_execution, swapped_layout) =
        build(&swapped_root, &["P2", "P1"]);
    let swapped = build_f5_marginals(
        &swapped_config,
        &swapped_config_sealed,
        &swapped_config_execution,
        &swapped_layout,
    )
    .expect("marginales al revés");
    assert_eq!(
        table_values(&first),
        table_values(&swapped),
        "reordenar las subpoblaciones declaradas no cambia la tabla"
    );

    let added_root = temporary("added-subpopulation");
    let (added_config, added_config_sealed, added_config_execution, added_layout) =
        build(&added_root, &["P1", "P2", "P3"]);
    let added = build_f5_marginals(
        &added_config,
        &added_config_sealed,
        &added_config_execution,
        &added_layout,
    )
    .expect("marginales ampliados");
    assert_eq!(
        added.blocks().len(),
        6,
        "la tercera subpoblación añade sus propios bloques"
    );
    assert_eq!(
        table_values(&first),
        table_values(&added)
            .into_iter()
            .filter(|(population, ..)| population != "P3")
            .collect::<Vec<_>>(),
        "añadir una subpoblación no altera los marginales de las anteriores"
    );

    // Las subpoblaciones no comparten flujo: sus columnas no son copias.
    let p1 = table_values(&added)
        .into_iter()
        .filter(|(population, ..)| population == "P1")
        .map(|(.., value)| value)
        .collect::<Vec<_>>();
    let p3 = table_values(&added)
        .into_iter()
        .filter(|(population, ..)| population == "P3")
        .map(|(.., value)| value)
        .collect::<Vec<_>>();
    assert_ne!(
        p1, p3,
        "dos subpoblaciones distintas no repiten el marginal"
    );

    for root in [first_root, swapped_root, added_root] {
        fs::remove_dir_all(&root).expect("cleanup");
    }
}

#[test]
fn the_same_configuration_regenerates_the_same_marginals() {
    let root = temporary("reproducible");
    let (config, config_sealed, config_execution, layout) = build(&root, &["P1", "P2"]);
    let first = build_f5_marginals(&config, &config_sealed, &config_execution, &layout)
        .expect("primeros marginales");
    let second = build_f5_marginals(&config, &config_sealed, &config_execution, &layout)
        .expect("segundos marginales");
    assert_eq!(first, second);
    fs::remove_dir_all(&root).expect("cleanup");
}

#[test]
fn a_layout_from_another_admission_cannot_produce_marginals() {
    let two_root = temporary("two-blocks");
    let (two_config, two_config_sealed, two_config_execution, two_layout) =
        build_with_blocks(&two_root, &["P1"], 2);
    let one_root = temporary("one-block");
    let (one_config, one_config_sealed, one_config_execution, one_layout) =
        build_with_blocks(&one_root, &["P1"], 1);
    assert_eq!(two_layout.blocks().len(), 2);
    assert_eq!(one_layout.blocks().len(), 1);

    // Geometría distinta: el layout ajeno no pertenece a esta ejecución.
    for (config, sealed, execution, layout) in [
        (
            &two_config,
            &two_config_sealed,
            &two_config_execution,
            &one_layout,
        ),
        (
            &one_config,
            &one_config_sealed,
            &one_config_execution,
            &two_layout,
        ),
    ] {
        let error = build_f5_marginals(config, sealed, execution, layout)
            .expect_err("un layout de otra admisión no puede sortear marginales");
        assert_eq!(error.code(), "GPS-CONF-101");
        assert_eq!(error.context()["fields"], "preflight.block_count");
    }

    // Y el caso que de verdad importa: misma geometría, otra ejecución. Antes
    // se aceptaba, y los marginales salían de identificadores ajenos sin que
    // ninguna capa posterior pudiera notarlo, porque el ámbito
    // `marginal_spectrum` se deriva de `block_id` y `site_id`.
    let other_execution =
        ExecutionId::parse("RUN_BQWGL4DQ5HJMZYWT3U4EDNLTJ4", "test").expect("otra ejecución");
    assert_ne!(other_execution, two_config_execution);
    let other_layout = build_f5_physical_layout(
        F5LayoutRequest::new(&two_config, &two_config_sealed, &other_execution)
            .expect("solicitud de layout de otra ejecución"),
    )
    .expect("layout de otra ejecución");
    assert_eq!(other_layout.blocks().len(), two_layout.blocks().len());
    assert_ne!(
        other_layout.blocks()[0].block_id(),
        two_layout.blocks()[0].block_id(),
        "dos ejecuciones distintas emiten identificadores distintos"
    );
    let error = build_f5_marginals(
        &two_config,
        &two_config_sealed,
        &two_config_execution,
        &other_layout,
    )
    .expect_err("un layout de igual geometría pero otra ejecución no puede sortear marginales");
    assert_eq!(error.code(), "GPS-CONF-101");
    assert_eq!(error.context()["fields"], "layout.identity");
    // Y con su propia ejecución sí procede.
    assert!(
        build_f5_marginals(
            &two_config,
            &two_config_sealed,
            &other_execution,
            &other_layout
        )
        .is_ok()
    );
    assert!(
        build_f5_marginals(
            &two_config,
            &two_config_sealed,
            &two_config_execution,
            &two_layout
        )
        .is_ok()
    );
    assert!(
        build_f5_marginals(
            &one_config,
            &one_config_sealed,
            &one_config_execution,
            &one_layout
        )
        .is_ok()
    );

    for root in [two_root, one_root] {
        fs::remove_dir_all(&root).expect("cleanup");
    }
}

/// DTC-089 exige `H <= 2^S` **antes** de sortear marginales, y este hito es
/// «antes». Con bloques de tres sitios caben ocho haplotipos: nueve no.
#[test]
fn more_haplotypes_than_the_block_can_represent_are_rejected_before_sampling() {
    let root = temporary("haplotype-cardinality");
    let general = write_fixture(&root, &["P1"], 1, 9);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let execution = ExecutionId::parse(EXECUTION, "test").expect("execution");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud de layout"),
    )
    .expect("layout físico");
    assert!(layout.blocks().iter().all(|block| block.sites().len() == 3));

    let error = build_f5_marginals(&config, &sealed, &execution, &layout)
        .expect_err("nueve haplotipos no caben en tres sitios");
    assert_eq!(error.code(), "GPS-DOMAIN-105");
    assert_eq!(error.context()["population_id"], "P1");
    assert!(
        error.context()["detail"].contains('9') && error.context()["detail"].contains('8'),
        "el detalle debe citar la cardinalidad y la capacidad: {:?}",
        error.context()["detail"]
    );

    // Ocho sí caben, de modo que la frontera es exacta y no un rechazo genérico.
    let root_ok = temporary("haplotype-cardinality-exact");
    let general = write_fixture(&root_ok, &["P1"], 1, 8);
    let config = load_and_validate_f5_config(&general).expect("configuración F5 admitida");
    let sealed = recover_sealed_f5_admission(&config, F5OtherResources::default())
        .expect("admisión sellada");
    let layout = build_f5_physical_layout(
        F5LayoutRequest::new(&config, &sealed, &execution).expect("solicitud de layout"),
    )
    .expect("layout físico");
    build_f5_marginals(&config, &sealed, &execution, &layout)
        .expect("ocho haplotipos caben exactamente en tres sitios");

    for root in [root, root_ok] {
        fs::remove_dir_all(&root).expect("cleanup");
    }
}
