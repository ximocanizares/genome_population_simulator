//! Contabilidad de bases de tracto: cota declarada, cifra realizada y tope.
//!
//! # Lo que hace que estas pruebas puedan fallar
//!
//! - **La cota y la cifra realizada tienen que diferir de verdad.** El fixture
//!   principal declara **tres** alelos de longitudes distintas sobre cuatro
//!   copias: es imposible que las cuatro saquen el más largo, de modo que lo
//!   realizado queda por debajo de la cota por aritmética y no por suerte. Si
//!   las dos cifras coincidieran, comparar la una con la otra no distinguiría
//!   nada.
//! - **Las dos cifras se calculan por caminos que no se tocan.** La cota sale
//!   de lo declarado —alelos, individuos y ploidía—; la realizada, de recorrer
//!   los documentos publicados. El juez independiente vuelve a contarla una
//!   tercera vez releyendo los bytes. Estropear una no puede arreglar la otra,
//!   y hay una prueba que lo exige.
//! - **El tope se comprueba con su pareja**: uno que deja pasar y otro que
//!   rechaza, sobre el mismo escenario. Sin el primero, un rechazo constante
//!   pasaría la mitad de la prueba.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::results::validate_published_tract_output_accounting;
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "000102030405060708090a0b0c0d0e0f101112131415161718191a1b1c1d1e1f";

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f6-4-8-1-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}

fn reference_bases() -> String {
    format!("{}{}{}", "G".repeat(100), "CA".repeat(50), "G".repeat(100))
}

/// Tres alelos de longitudes **distintas**: 90, 80 y 70 bases.
///
/// Con cuatro copias es imposible que todas saquen el de 90, que es lo que la
/// cota supone, de modo que lo realizado queda estrictamente por debajo.
const THREE_UNEVEN_ALLELES: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 0.4

[[reference.repeat_loci.alleles]]
units = 40
frequency = 0.35

[[reference.repeat_loci.alleles]]
units = 35
frequency = 0.25
"#;

fn scenario(loci: &str, cap: Option<u64>, vcf: bool) -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    let bases = reference_bases();
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let cap = cap.map_or_else(String::new, |value| {
        format!("max_tract_output_bases = {value}\n")
    });
    let vcf = if vcf {
        "\n[vcf_output]\nmode = \"text\"\n"
    } else {
        ""
    };
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            r#"config_kind = "general"
schema_version = "4"
seed = "{SEED}"

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

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"
{cap}
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = {}
{loci}
[[initial_populations]]
population_id = "P1"
individual_count = 2
individual_ids = ["I1", "I2"]
ploidy = 2
genotype_model = "hwe"
{vcf}"#,
            bases.len()
        ),
    )
    .expect("configuración");
    (root, config)
}

struct Published {
    root: PathBuf,
    destination: PathBuf,
    manifest: Value,
    effective: Value,
}

fn publish(loci: &str, cap: Option<u64>, vcf: bool) -> Published {
    let (root, config) = scenario(loci, cap, vcf);
    let effective_config = load_and_validate_config(&config).expect("configuración válida");
    let destination = root.join("published");
    execute_validated_config(&effective_config, &destination).expect("publicación");
    let read = |name: &str| -> Value {
        serde_json::from_slice(&fs::read(destination.join(name)).expect(name)).expect("JSON")
    };
    let manifest = read("manifest.json");
    let effective = read("effective-config.json");
    Published {
        root,
        destination,
        manifest,
        effective,
    }
}

fn judge(published: &Published, vcf: bool) -> Result<(), String> {
    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let vcf_path = published.destination.join("outputs/truth.vcf");
    validate_published_tract_output_accounting(
        &effective_path,
        &fs::read(&effective_path).expect("configuración efectiva"),
        &manifest_path,
        &fs::read(&manifest_path).expect("manifiesto"),
        &published.destination.join("outputs/truth.jsonl"),
        vcf.then_some(vcf_path.as_path()),
    )
    .map_err(|error| format!("{error:?}"))
}

/// Cuenta las bases de tracto releyendo los ficheros con Python, que es un
/// tercer camino ajeno por completo al productor y al juez de Rust.
fn counted_with_python(destination: &Path, vcf: bool) -> (u64, u64) {
    let script = r##"
import json, sys
truth = 0
for line in open(sys.argv[1]):
    r = json.loads(line)
    for piece in r.get("current_pieces", []) + r.get("history", []):
        if piece.get("kind") == "repeat_allele":
            truth += len(piece["bases"])
    for allele in r.get("alleles", []):
        if allele.get("kind") == "repeat":
            truth += len(allele["bases"])
vcf = 0
if len(sys.argv) > 2:
    for line in open(sys.argv[2]):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if any(k.startswith("RN=") for k in f[7].split(";")):
            vcf += len(f[3]) + sum(len(a) for a in f[4].split(","))
print(truth, vcf)
"##;
    let mut command = Command::new("python3");
    command
        .arg("-c")
        .arg(script)
        .arg(destination.join("outputs/truth.jsonl"));
    if vcf {
        command.arg(destination.join("outputs/truth.vcf"));
    }
    let output = command.output().expect("contador independiente en Python");
    assert!(
        output.status.success(),
        "el contador de Python falló: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let text = String::from_utf8(output.stdout).expect("UTF-8");
    let mut values = text.split_whitespace();
    (
        values
            .next()
            .expect("bases de verdad")
            .parse()
            .expect("entero"),
        values
            .next()
            .expect("bases de VCF")
            .parse()
            .expect("entero"),
    )
}

// ---------------------------------------------------------------------------
// Pruebas
// ---------------------------------------------------------------------------

/// **La cota se publica en la configuración efectiva y se calcula sobre lo
/// declarado.**
///
/// *La cuenta, hecha a mano:* dos individuos diploides son cuatro copias; el
/// alelo más largo mide noventa bases; la verdad publica el tracto de cada
/// copia en tres sitios. Tres por cuatro por noventa son **1080**. El VCF
/// publica el tracto de la referencia, cien bases, más los tres alelos
/// declarados: cien más noventa más ochenta más setenta son **340**.
#[test]
fn the_effective_configuration_publishes_the_bound_derived_from_the_declaration() {
    let published = publish(THREE_UNEVEN_ALLELES, None, true);
    let budget = &published.effective["repeat_output_budget"];
    assert_eq!(budget["truth_copies_per_haplotype_copy"], 3);
    assert_eq!(budget["truth_tract_bases_bound"], 3 * 4 * 90);
    assert_eq!(budget["vcf_tract_bases_bound"], 100 + 90 + 80 + 70);
    assert_eq!(budget["total_tract_bases_bound"], 1080 + 340);
    assert!(budget["declared_max_tract_output_bases"].is_null());
    assert_eq!(budget["per_locus"][0]["id"], "STR-1");
    assert_eq!(budget["per_locus"][0]["truth_tract_bases_bound"], 1080);
    fs::remove_dir_all(&published.root).ok();
}

/// **Lo realizado queda por debajo de la cota, y las dos cifras difieren de
/// verdad.**
///
/// *Qué hace que pueda fallar:* con tres alelos de longitudes distintas y solo
/// cuatro copias es **imposible** que lo realizado alcance la cota, que supone
/// que las cuatro sacan el más largo y que los tres alelos llegan al VCF. La
/// desigualdad estricta se exige, no solo la no superación: si alguien igualara
/// las dos cifras «arreglando» la diferencia, esta prueba moriría.
#[test]
fn the_realized_figures_stay_strictly_below_the_bound() {
    let published = publish(THREE_UNEVEN_ALLELES, None, true);
    let accounting = &published.manifest["repeat_loci"]["tract_output"];
    let bound = accounting["total_tract_bases_bound"]
        .as_u64()
        .expect("cota");
    let realized = accounting["realized_total_tract_bases"]
        .as_u64()
        .expect("realizado");
    assert_eq!(bound, 1080 + 340, "la cota es la que la declaración fija");
    assert!(
        realized < bound,
        "lo realizado ({realized}) tiene que quedar por debajo de la cota ({bound})"
    );
    assert!(realized > 0, "y el escenario sí publica tracto");
    assert_eq!(
        accounting["realized_truth_tract_bases"]
            .as_u64()
            .expect("verdad")
            + accounting["realized_vcf_tract_bases"]
                .as_u64()
                .expect("vcf"),
        realized,
        "el total es la suma de sus dos términos"
    );
    fs::remove_dir_all(&published.root).ok();
}

/// **El juez independiente acepta la contabilidad publicada y un contador
/// ajeno da la misma cifra.**
///
/// *Qué hace que pueda fallar:* el contador de Python relee los mismos ficheros
/// sin compartir una línea con el productor ni con el juez de Rust. Si el
/// productor contara mal, las tres cifras no coincidirían. La comparación con
/// la cota se exige también aquí, sobre los documentos ya escritos.
#[test]
fn the_independent_judge_and_a_foreign_counter_agree_with_the_manifest() {
    let published = publish(THREE_UNEVEN_ALLELES, None, true);
    judge(&published, true).expect("el juez acepta la contabilidad publicada");

    let (truth_bases, vcf_bases) = counted_with_python(&published.destination, true);
    let accounting = &published.manifest["repeat_loci"]["tract_output"];
    assert_eq!(
        accounting["realized_truth_tract_bases"].as_u64(),
        Some(truth_bases),
        "el contador ajeno no ve las mismas bases de verdad"
    );
    assert_eq!(
        accounting["realized_vcf_tract_bases"].as_u64(),
        Some(vcf_bases),
        "el contador ajeno no ve las mismas bases de VCF"
    );
    assert!(
        truth_bases > 0 && vcf_bases > 0,
        "control: hay algo que contar"
    );
    fs::remove_dir_all(&published.root).ok();
}

/// **El juez rechaza una contabilidad manipulada.**
///
/// *Qué hace que pueda fallar:* se altera el manifiesto de cuatro formas, cada
/// una dirigida a una comprobación distinta, y se exige el rechazo de cada una.
/// La primera aserción fija que el documento intacto se acepta.
#[test]
fn the_judge_rejects_a_tampered_accounting() {
    let published = publish(THREE_UNEVEN_ALLELES, None, true);
    judge(&published, true).expect("el par intacto se acepta");

    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("configuración efectiva");
    let truth_path = published.destination.join("outputs/truth.jsonl");
    let vcf_path = published.destination.join("outputs/truth.vcf");

    /// Una manipulación del manifiesto, si hay que rebajar la cota también en
    /// la configuración efectiva, y el motivo por el que se hace.
    type Tampering = (&'static str, bool, Box<dyn Fn(&mut Value)>);

    let cases: [Tampering; 5] = [
        (
            "bases de verdad realizadas infladas",
            false,
            Box::new(|document: &mut Value| {
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let truth = accounting["realized_truth_tract_bases"]
                    .as_u64()
                    .expect("verdad");
                accounting["realized_truth_tract_bases"] = serde_json::json!(truth + 2);
                let total = accounting["realized_total_tract_bases"]
                    .as_u64()
                    .expect("total");
                accounting["realized_total_tract_bases"] = serde_json::json!(total + 2);
            }),
        ),
        (
            "bases de VCF realizadas rebajadas",
            false,
            Box::new(|document: &mut Value| {
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let vcf = accounting["realized_vcf_tract_bases"]
                    .as_u64()
                    .expect("vcf");
                accounting["realized_vcf_tract_bases"] = serde_json::json!(vcf - 2);
                let total = accounting["realized_total_tract_bases"]
                    .as_u64()
                    .expect("total");
                accounting["realized_total_tract_bases"] = serde_json::json!(total - 2);
            }),
        ),
        (
            "cota distinta de la de la configuración efectiva",
            false,
            Box::new(|document: &mut Value| {
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let total = accounting["total_tract_bases_bound"]
                    .as_u64()
                    .expect("cota");
                accounting["total_tract_bases_bound"] = serde_json::json!(total + 1);
            }),
        ),
        (
            "cota rebajada por debajo de lo realizado en los dos documentos",
            true,
            Box::new(|document: &mut Value| {
                // Se rebaja la cota **en el manifiesto**; la de la
                // configuración efectiva se rebaja igual en el sitio donde se
                // usa esta manipulación, de modo que las dos siguen
                // coincidiendo y lo único que puede fallar es la comparación
                // con lo realizado.
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let realized = accounting["realized_total_tract_bases"]
                    .as_u64()
                    .expect("realizado");
                accounting["total_tract_bases_bound"] = serde_json::json!(realized - 1);
            }),
        ),
        (
            "total realizado que no es la suma de sus términos",
            false,
            Box::new(|document: &mut Value| {
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let total = accounting["realized_total_tract_bases"]
                    .as_u64()
                    .expect("total");
                accounting["realized_total_tract_bases"] = serde_json::json!(total + 1);
            }),
        ),
    ];
    for (reason, pair_effective, mutate) in cases {
        let mut document: Value =
            serde_json::from_slice(&fs::read(&manifest_path).expect("manifiesto")).expect("JSON");
        mutate(&mut document);
        // Solo la fila de la cota rebajada toca **también** la configuración
        // efectiva: si no la tocara, el rechazo vendría de que las dos no
        // coinciden y no de la comparación con lo realizado, que es la guarda
        // que esa fila protege. Las demás filas dejan el par intacto, porque
        // la guarda que protegen es otra.
        let paired = if pair_effective {
            let lowered = document["repeat_loci"]["tract_output"]["total_tract_bases_bound"]
                .as_u64()
                .expect("cota");
            let mut effective: Value =
                serde_json::from_slice(&effective_bytes).expect("configuración efectiva");
            effective["repeat_output_budget"]["total_tract_bases_bound"] =
                serde_json::json!(lowered);
            serde_json::to_vec(&effective).expect("configuración manipulada")
        } else {
            effective_bytes.clone()
        };
        let tampered = serde_json::to_vec(&document).expect("JSON manipulado");
        assert!(
            validate_published_tract_output_accounting(
                &effective_path,
                &paired,
                &manifest_path,
                &tampered,
                &truth_path,
                Some(&vcf_path),
            )
            .is_err(),
            "el juez tiene que rechazar «{reason}»"
        );
    }
    fs::remove_dir_all(&published.root).ok();
}

/// **El tope rechaza el escenario antes de simular, y su pareja lo deja
/// pasar.**
///
/// *Qué hace que pueda fallar:* los dos topes se aplican al **mismo**
/// escenario, uno por encima de la cota y otro por debajo. Sin el primero, una
/// implementación que rechazara siempre pasaría la mitad de la prueba. El
/// diagnóstico tiene que nombrar el tope, la cota y el locus que más aporta,
/// que es lo que lo hace accionable.
#[test]
fn a_declared_cap_rejects_the_scenario_before_simulating_and_its_pair_admits_it() {
    const BOUND: u64 = 1080 + 340;

    let (permissive_root, permissive) = scenario(THREE_UNEVEN_ALLELES, Some(BOUND), true);
    let effective = load_and_validate_config(&permissive).expect("un tope igual a la cota admite");
    assert_eq!(
        effective.repeat_output_budget().declared_cap(),
        Some(BOUND),
        "el tope declarado se conserva"
    );
    assert_eq!(effective.repeat_output_budget().total_bases(), BOUND);
    fs::remove_dir_all(&permissive_root).ok();

    let (strict_root, strict) = scenario(THREE_UNEVEN_ALLELES, Some(BOUND - 1), true);
    let error =
        load_and_validate_config(&strict).expect_err("un tope por debajo de la cota rechaza");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(
        detail.contains(&BOUND.to_string())
            && detail.contains(&(BOUND - 1).to_string())
            && detail.contains("STR-1"),
        "el rechazo tiene que nombrar la cota, el tope y el locus: {detail}"
    );
    assert!(
        !strict_root.join("published").exists(),
        "el rechazo llega antes de publicar nada"
    );
    fs::remove_dir_all(&strict_root).ok();
}

/// **Un escenario sin VCF sigue cuadrando: lo realizado en el VCF es cero.**
#[test]
fn a_run_without_a_vcf_reports_zero_realized_vcf_bases() {
    let published = publish(THREE_UNEVEN_ALLELES, None, false);
    let accounting = &published.manifest["repeat_loci"]["tract_output"];
    assert_eq!(accounting["realized_vcf_tract_bases"], 0);
    assert!(
        accounting["realized_truth_tract_bases"]
            .as_u64()
            .expect("verdad")
            > 0,
        "la verdad sí publica tracto aunque no haya VCF"
    );
    judge(&published, false).expect("el juez acepta una ejecución sin VCF");
    fs::remove_dir_all(&published.root).ok();
}
