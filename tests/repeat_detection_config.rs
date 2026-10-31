//! Entrega 3 del plan de detección de microsatélites: **la sección de
//! microsatélites y su resolución** (pasos F6.4.a.9, F6.4.a.10b y F6.4.a.13).
//!
//! Cada prueba escribe un documento general pequeño y lo resuelve con
//! `load_and_validate_config`, sin ejecutar nada: lo que se juzga es la
//! frontera de configuración. Lo que la ejecución publica lo prueba
//! `tests/repeat_detection_run.rs`.

use std::fs;
use std::path::PathBuf;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::{EffectiveConfig, load_and_validate_config};
use genome_population_simulator::diagnostics::Diagnostic;

static NONCE: AtomicU64 = AtomicU64::new(0);

const HEAD: &str = "config_kind = \"general\"\nschema_version = \"4\"\n\
     seed = \"00000000000000000000000000000000000000000000000000000000000000c3\"\n\n\
     [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
     [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 1000\n\n\
     [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 500\n\n";

const SLIPPAGE: &str =
    "[reference.repeat_detection.slippage]\nmodel = \"stepwise\"\nmean_events = 1.0\n\n";

const LOCUS: &str = "[[reference.repeat_loci]]\nid = \"STR-1\"\nsequence = \"chr1\"\n\
     start = 101\nend = 120\nmotif = \"CA\"\n\n";

/// Resuelve `HEAD` más `body`. El FASTA no se lee al resolver.
fn resolve(body: &str) -> Result<EffectiveConfig, Diagnostic> {
    let root: PathBuf = std::env::temp_dir().join(format!(
        "gps-detection-config-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&root).expect("directorio temporal");
    let config = root.join("general.toml");
    fs::write(&config, format!("{HEAD}{body}")).expect("configuración");
    let result = load_and_validate_config(&config);
    let _ = fs::remove_dir_all(&root);
    result
}

fn rejected(body: &str) -> String {
    resolve(body)
        .map(|_| ())
        .expect_err("el documento debe rechazarse")
        .render_human()
}

/// **Sin la sección no hay microsatélites** (REQ-CONF-142): un locus
/// declarado se rechaza con un diagnóstico que remite a la sección.
#[test]
fn a_declared_locus_without_the_section_is_rejected_pointing_to_it() {
    let text = rejected(&format!(
        "{LOCUS}[[reference.repeat_loci.alleles]]\nunits = 10\nfrequency = 1.0\n\n"
    ));
    assert!(
        text.contains("reference.repeat_detection") && text.contains("mean_events"),
        "{text}"
    );
}

/// **El deslizamiento común es obligatorio y sin valores ocultos**: falta la
/// tabla, falta el modelo o falta la media, y el rechazo nombra el campo y
/// trae un ejemplo. Un valor por omisión de `mean_events` haría pasar el
/// tercer caso.
#[test]
fn the_common_slippage_is_required_without_hidden_defaults() {
    let cases = [
        (
            "[reference.repeat_detection]\nmin_copies = [8, 6, 6, 6, 5, 4]\n\n",
            "reference.repeat_detection.slippage",
        ),
        (
            "[reference.repeat_detection]\n\n[reference.repeat_detection.slippage]\n\
             mean_events = 1.0\n\n",
            "reference.repeat_detection.slippage.model",
        ),
        (
            "[reference.repeat_detection]\n\n[reference.repeat_detection.slippage]\n\
             model = \"stepwise\"\n\n",
            "reference.repeat_detection.slippage.mean_events",
        ),
    ];
    for (body, field) in cases {
        let text = rejected(body);
        assert!(text.contains(field), "falta el campo {field}: {text}");
        assert!(
            text.contains("por ejemplo") && text.contains("model = \\\"stepwise\\\""),
            "el mensaje trae un ejemplo: {text}"
        );
    }
    // Y con los dos, la sección se admite y el sesgo y el paso toman los de
    // REQ-CONF-139.
    let effective = resolve(&format!("[reference.repeat_detection]\n\n{SLIPPAGE}"))
        .expect("la sección mínima se admite");
    let detection = effective.repeat_detection().expect("sección");
    assert_eq!(detection.slippage().mean_events(), 1.0);
    assert_eq!(detection.slippage().expansion_bias(), 0.5);
}

/// **Los suelos y el tope por omisión** son `[8, 6, 6, 6, 5, 4]` y 200, y lo
/// declarado los sustituye. Una tabla de cinco valores, un suelo cero y un
/// tope cero se rechazan.
#[test]
fn the_floors_and_the_cap_have_their_defaults_and_their_bounds() {
    let effective = resolve(&format!("[reference.repeat_detection]\n\n{SLIPPAGE}"))
        .expect("sección sin suelos");
    let detection = effective.repeat_detection().expect("sección");
    assert_eq!(detection.min_copies(), &[8, 6, 6, 6, 5, 4]);
    assert_eq!(detection.max_units(), 200);
    assert!(detection.regions().is_none());
    let effective = resolve(&format!(
        "[reference.repeat_detection]\nmin_copies = [9, 7, 5, 5, 4, 3]\nmax_units = 50\n\n{SLIPPAGE}"
    ))
    .expect("suelos declarados");
    let detection = effective.repeat_detection().expect("sección");
    assert_eq!(detection.min_copies(), &[9, 7, 5, 5, 4, 3]);
    assert_eq!(detection.max_units(), 50);
    for (body, needle) in [
        (
            "min_copies = [8, 6, 6, 6, 5]\n",
            "min_copies lleva 5 valores",
        ),
        ("min_copies = [8, 6, 0, 6, 5, 4]\n", "motivos de 3 pb es 0"),
        ("max_units = 0\n", "reference.repeat_detection.max_units"),
    ] {
        let text = rejected(&format!("[reference.repeat_detection]\n{body}\n{SLIPPAGE}"));
        assert!(text.contains(needle), "{body}: {text}");
    }
}

/// **Un declarado sin alelos ni generador hereda el deslizamiento común**
/// (REQ-CONF-139); uno con generador propio conserva el suyo.
#[test]
fn a_declared_locus_without_alleles_inherits_the_common_slippage() {
    let effective = resolve(&format!(
        "[reference.repeat_detection]\n\n[reference.repeat_detection.slippage]\n\
         model = \"stepwise\"\nmean_events = 3.0\nexpansion_bias = 0.25\n\n\
         {LOCUS}\
         [[reference.repeat_loci]]\nid = \"STR-2\"\nsequence = \"chr2\"\nstart = 11\nend = 30\n\
         motif = \"GT\"\n\n[reference.repeat_loci.slippage]\nmodel = \"stepwise\"\n\
         mean_events = 0.5\n\n"
    ))
    .expect("los declarados se admiten");
    let loci = effective.repeat_loci();
    let inherited = loci[0].slippage().expect("hereda").settings();
    assert_eq!(
        (inherited.mean_events(), inherited.expansion_bias()),
        (3.0, 0.25)
    );
    let own = loci[1].slippage().expect("el suyo").settings();
    assert_eq!((own.mean_events(), own.expansion_bias()), (0.5, 0.5));
    // Una lista de alelos vacía sigue sin decir de dónde salen.
    let text = rejected(&format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}{LOCUS}alleles = []\n\n"
    ));
    assert!(text.contains("no declara ningún alelo"), "{text}");
}

/// **Detección y reducción no van en la misma ejecución** (DP-SSR-12), y el
/// diagnóstico da el camino en dos ejecuciones y el campo `fasta` de la
/// segunda.
#[test]
fn detection_with_a_reduction_in_the_same_run_points_to_two_runs() {
    let text = rejected(&format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}\
         [[reduced_references]]\nid = \"ref-red\"\nsource_reference_id = \"ref-main\"\n\
         region_size_bases = 100\nreduction_percent = 50\n\n"
    ));
    assert!(
        text.contains("dos ejecuciones") && text.contains("reference.fasta"),
        "{text}"
    );
}

/// **Las regiones** viven en la referencia final, en base uno y cerradas:
/// una secuencia desconocida, un tramo fuera de la secuencia, dos regiones
/// que se solapan y una lista vacía se rechazan; dos que solo se tocan se
/// admiten y salen ordenadas.
#[test]
fn detection_regions_are_checked_against_the_final_catalog() {
    let region = |sequence: &str, start: u64, end: u64| {
        format!(
            "[[reference.repeat_detection.regions]]\nsequence = \"{sequence}\"\nstart = {start}\n\
             end = {end}\n\n"
        )
    };
    for (regions, needle) in [
        (region("chr9", 1, 10), "chr9"),
        (region("chr1", 990, 1001), "no cabe"),
        (region("chr1", 0, 10), "no cabe"),
        (
            format!("{}{}", region("chr1", 100, 200), region("chr1", 200, 300)),
            "se solapan",
        ),
        ("regions = []\n\n".to_owned(), "vacía"),
    ] {
        let body = if regions.starts_with("regions") {
            format!("[reference.repeat_detection]\n{regions}{SLIPPAGE}")
        } else {
            format!("[reference.repeat_detection]\n\n{SLIPPAGE}{regions}")
        };
        let text = rejected(&body);
        assert!(text.contains(needle), "{needle}: {text}");
    }
    let effective = resolve(&format!(
        "[reference.repeat_detection]\n\n{SLIPPAGE}{}{}{}",
        region("chr2", 1, 10),
        region("chr1", 201, 300),
        region("chr1", 100, 200)
    ))
    .expect("regiones que solo se tocan");
    let regions = effective
        .repeat_detection()
        .and_then(|detection| detection.regions())
        .expect("regiones");
    let found: Vec<(&str, u64, u64)> = regions
        .iter()
        .map(|region| (region.sequence_id(), region.start(), region.end()))
        .collect();
    assert_eq!(
        found,
        vec![("chr1", 100, 200), ("chr1", 201, 300), ("chr2", 1, 10)]
    );
}
