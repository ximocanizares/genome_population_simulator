//! Contabilidad de bases de tracto: cota declarada, cifra realizada y tope.
//!
//! **Desde el 2 de octubre de 2026 la parte VCF de la cota se mide, no se
//! estima** (DTC-097): es lo que el VCF lleva, contado al publicarlo antes de
//! escribirlo. La parte de la verdad sigue saliendo de lo declarado.
//!
//! # Lo que hace que estas pruebas puedan fallar
//!
//! - **La cota de la verdad y su cifra realizada tienen que diferir de
//!   verdad.** El fixture principal declara **tres** alelos de longitudes
//!   distintas sobre cuatro copias: es imposible que las cuatro saquen el más
//!   largo, de modo que lo realizado queda por debajo por aritmética y no por
//!   suerte. La parte del VCF, en cambio, tiene que coincidir con lo que un
//!   contador ajeno ve en el fichero.
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

/// Un único alelo de 45 unidades: todas las copias lo sacan.
const SINGLE_ALLELE: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 1.0
"#;

fn scenario(loci: &str, cap: Option<u64>, vcf: bool) -> (PathBuf, PathBuf) {
    scenario_on(&reference_bases(), loci, cap, vcf.then_some("text"))
}

/// El escenario sobre una referencia dada y con el VCF en el modo pedido
/// (`text` o `bgzf`), o sin VCF.
fn scenario_on(
    bases: &str,
    loci: &str,
    cap: Option<u64>,
    vcf_mode: Option<&str>,
) -> (PathBuf, PathBuf) {
    let root = temporary("run");
    fs::create_dir_all(&root).expect("directorio");
    fs::write(root.join("reference.fa"), format!(">chr1\n{bases}\n")).expect("FASTA");
    let cap = cap.map_or_else(String::new, |value| {
        format!("max_tract_output_bases = {value}\n")
    });
    let vcf = vcf_mode.map_or_else(String::new, |mode| {
        format!("\n[vcf_output]\nmode = \"{mode}\"\n")
    });
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

fn publish(loci: &str, vcf: bool) -> Published {
    publish_on(&reference_bases(), loci, vcf.then_some("text"))
}

fn publish_on(bases: &str, loci: &str, vcf_mode: Option<&str>) -> Published {
    publish_scenario(scenario_on(bases, loci, None, vcf_mode))
}

fn publish_scenario((root, config): (PathBuf, PathBuf)) -> Published {
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

/// **La cota se publica en la configuración efectiva: la de la verdad sale de
/// lo declarado y la del VCF es lo que el VCF lleva.**
///
/// *La cuenta, hecha a mano:* dos individuos diploides son cuatro copias; el
/// alelo más largo mide noventa bases; la verdad publica el tracto de cada
/// copia en tres sitios. Tres por cuatro por noventa son **1080**. La parte del
/// VCF no se calcula aquí: se exige que sea la que un contador ajeno, en
/// Python, cuenta en el VCF publicado, y que no sea la estimación de antes
/// —cien del tracto más noventa, ochenta y setenta de los alelos, 340—, que
/// cuatro copias no pueden llenar.
#[test]
fn the_effective_configuration_publishes_the_truth_bound_and_the_measured_vcf_bases() {
    let published = publish(THREE_UNEVEN_ALLELES, true);
    let (_, vcf_bases) = counted_with_python(&published.destination, true);
    let budget = &published.effective["repeat_output_budget"];
    assert_eq!(budget["truth_copies_per_haplotype_copy"], 3);
    assert_eq!(budget["truth_tract_bases_bound"], 3 * 4 * 90);
    assert_eq!(budget["vcf_tract_bases_bound"].as_u64(), Some(vcf_bases));
    assert_ne!(
        vcf_bases, 340,
        "control: lo medido no es la estimación de antes"
    );
    assert!(vcf_bases > 0, "control: el VCF publica tracto");
    assert_eq!(
        budget["total_tract_bases_bound"].as_u64(),
        Some(1080 + vcf_bases)
    );
    assert!(budget["declared_max_tract_output_bases"].is_null());
    assert_eq!(budget["per_locus"][0]["id"], "STR-1");
    assert_eq!(budget["per_locus"][0]["truth_tract_bases_bound"], 1080);
    assert_eq!(
        budget["per_locus"][0]["vcf_tract_bases_bound"].as_u64(),
        Some(vcf_bases)
    );
    fs::remove_dir_all(&published.root).ok();
}

/// **Lo realizado en la verdad queda por debajo de su cota; lo del VCF es su
/// cota.**
///
/// *Qué hace que pueda fallar:* con tres alelos de longitudes distintas y solo
/// cuatro copias es **imposible** que la verdad alcance su cota, que supone que
/// las cuatro sacan el más largo; la desigualdad estricta se exige. La parte
/// del VCF se mide, de modo que lo realizado y la cota tienen que coincidir.
#[test]
fn the_realized_truth_stays_below_its_bound_and_the_vcf_part_is_measured() {
    let published = publish(THREE_UNEVEN_ALLELES, true);
    let accounting = &published.manifest["repeat_loci"]["tract_output"];
    let read = |key: &str| accounting[key].as_u64().expect(key);
    assert_eq!(read("truth_tract_bases_bound"), 1080);
    assert!(
        read("realized_truth_tract_bases") < read("truth_tract_bases_bound"),
        "lo realizado en la verdad tiene que quedar por debajo de su cota"
    );
    assert!(
        read("realized_truth_tract_bases") > 0,
        "y la verdad sí publica tracto"
    );
    assert_eq!(
        read("realized_vcf_tract_bases"),
        read("vcf_tract_bases_bound"),
        "la parte VCF de la cota es lo que el VCF lleva"
    );
    assert_eq!(
        read("total_tract_bases_bound"),
        read("truth_tract_bases_bound") + read("vcf_tract_bases_bound")
    );
    assert_eq!(
        read("realized_truth_tract_bases") + read("realized_vcf_tract_bases"),
        read("realized_total_tract_bases"),
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
    let published = publish(THREE_UNEVEN_ALLELES, true);
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
    let published = publish(THREE_UNEVEN_ALLELES, true);
    judge(&published, true).expect("el par intacto se acepta");

    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let effective_bytes = fs::read(&effective_path).expect("configuración efectiva");
    let truth_path = published.destination.join("outputs/truth.jsonl");
    let vcf_path = published.destination.join("outputs/truth.vcf");

    /// Una manipulación del manifiesto, si hay que rebajar la cota también en
    /// la configuración efectiva, y el motivo por el que se hace.
    type Tampering = (&'static str, bool, Box<dyn Fn(&mut Value)>);

    let cases: [Tampering; 6] = [
        (
            "parte VCF de la cota distinta de lo que el VCF lleva",
            true,
            Box::new(|document: &mut Value| {
                // La cota se sube en el manifiesto y en la configuración
                // efectiva por igual, y lo realizado no se toca: lo único que
                // puede fallar es que la parte VCF de la cota no sea lo que el
                // VCF publicado lleva.
                let accounting = &mut document["repeat_loci"]["tract_output"];
                let vcf = accounting["vcf_tract_bases_bound"].as_u64().expect("vcf");
                accounting["vcf_tract_bases_bound"] = serde_json::json!(vcf + 3);
                let total = accounting["total_tract_bases_bound"]
                    .as_u64()
                    .expect("cota");
                accounting["total_tract_bases_bound"] = serde_json::json!(total + 3);
            }),
        ),
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
                // con lo realizado. La parte VCF no se toca: sigue siendo lo
                // que el VCF lleva.
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
    // Las tres guardas de la parte VCF de la cota se exigen **por su motivo**,
    // cada una con una manipulación que solo ella puede rechazar: con un
    // simple «se rechaza», quitar una la dejaría cubierta por la siguiente.
    let reject = |effective: &Value, manifest: &Value| {
        format!(
            "{:?}",
            validate_published_tract_output_accounting(
                &effective_path,
                &serde_json::to_vec(effective).expect("JSON"),
                &manifest_path,
                &serde_json::to_vec(manifest).expect("JSON"),
                &truth_path,
                Some(&vcf_path),
            )
            .expect_err("la manipulación tiene que rechazarse")
        )
    };
    let pristine_manifest: Value =
        serde_json::from_slice(&fs::read(&manifest_path).expect("manifiesto")).expect("JSON");
    let pristine_effective: Value =
        serde_json::from_slice(&effective_bytes).expect("configuración efectiva");
    {
        // Las dos cifras VCF suben igual en los dos documentos, de modo que
        // coinciden entre sí; lo que falla es que no son lo que el VCF lleva.
        let mut manifest = pristine_manifest.clone();
        let mut effective = pristine_effective.clone();
        for (document, section) in [
            (&mut manifest, "/repeat_loci/tract_output"),
            (&mut effective, "/repeat_output_budget"),
        ] {
            let block = document.pointer_mut(section).expect("sección");
            for key in ["vcf_tract_bases_bound", "total_tract_bases_bound"] {
                let value = block[key].as_u64().expect("entero");
                block[key] = serde_json::json!(value + 3);
            }
        }
        let locus = &mut effective["repeat_output_budget"]["per_locus"][0];
        let value = locus["vcf_tract_bases_bound"].as_u64().expect("entero");
        locus["vcf_tract_bases_bound"] = serde_json::json!(value + 3);
        let detail = reject(&effective, &manifest);
        assert!(
            detail.contains("la parte VCF de la cota no es lo que el VCF publicado lleva"),
            "{detail}"
        );
    }
    {
        // El manifiesto queda intacto; la configuración efectiva dice otra
        // cifra VCF sin tocar su total, que sigue coincidiendo.
        let mut effective = pristine_effective.clone();
        let value = effective["repeat_output_budget"]["vcf_tract_bases_bound"]
            .as_u64()
            .expect("entero");
        effective["repeat_output_budget"]["vcf_tract_bases_bound"] = serde_json::json!(value + 3);
        let detail = reject(&effective, &pristine_manifest);
        assert!(
            detail.contains(
                "la configuración efectiva no publica las bases de tracto que el VCF lleva"
            ),
            "{detail}"
        );
    }
    {
        // Solo cambia el reparto por locus.
        let mut effective = pristine_effective.clone();
        let locus = &mut effective["repeat_output_budget"]["per_locus"][0];
        let value = locus["vcf_tract_bases_bound"].as_u64().expect("entero");
        locus["vcf_tract_bases_bound"] = serde_json::json!(value + 3);
        let detail = reject(&effective, &pristine_manifest);
        assert!(
            detail.contains("las bases de tracto del VCF por locus no suman las del VCF publicado"),
            "{detail}"
        );
    }

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
            let accounting = &document["repeat_loci"]["tract_output"];
            let mut effective: Value =
                serde_json::from_slice(&effective_bytes).expect("configuración efectiva");
            for key in ["total_tract_bases_bound", "vcf_tract_bases_bound"] {
                effective["repeat_output_budget"][key] = accounting[key].clone();
            }
            effective["repeat_output_budget"]["per_locus"][0]["vcf_tract_bases_bound"] =
                accounting["vcf_tract_bases_bound"].clone();
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

/// **El tope se aplica en dos momentos, cada uno con su pareja.**
///
/// - **Al admitir**, por precaución, sobre la cota de la verdad, que es una
///   cota superior: un tope por debajo de ella rechaza antes de simular. Una
///   ejecución así podría caber, porque lo realizado suele quedar por debajo.
/// - **Al publicar el VCF**, sobre lo que la ejecución va a publicar de
///   verdad: la verdad ya construida más lo que el VCF va a escribir. Un tope
///   una base por debajo rechaza sin escribir el VCF y deja un paquete
///   fallido; uno igual publica. Es el caso en que `validate` acepta y `run`
///   rechaza después.
///
/// Sin las parejas, una implementación que rechazara siempre pasaría la mitad
/// de la prueba. Los dos diagnósticos tienen que nombrar el tope, la cifra y
/// el locus que más aporta.
#[test]
fn a_declared_cap_is_applied_at_admission_and_when_publishing_the_vcf() {
    const TRUTH_BOUND: u64 = 1080;

    let (strict_root, strict) = scenario(THREE_UNEVEN_ALLELES, Some(TRUTH_BOUND - 1), true);
    let error =
        load_and_validate_config(&strict).expect_err("un tope por debajo de la verdad rechaza");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(
        detail.contains(&TRUTH_BOUND.to_string())
            && detail.contains(&(TRUTH_BOUND - 1).to_string())
            && detail.contains("STR-1"),
        "el rechazo tiene que nombrar la cota, el tope y el locus: {detail}"
    );
    assert!(
        !strict_root.join("published").exists(),
        "el rechazo llega antes de publicar nada"
    );
    fs::remove_dir_all(&strict_root).ok();

    // Para que el tope rechace **al publicar** hace falta que la verdad
    // realizada más el VCF superen un tope que la cota de la verdad no
    // supera. Con un único alelo de 45 unidades toda copia lo saca: la verdad
    // realizada son dos sitios por copia y uno por individuo homocigoto,
    // 2·4·90 + 2·90 = 900, y el VCF lleva el tracto y el alelo, 190. Lo que se
    // publica son 1090, por encima de la cota de la verdad, 1080.
    let measured = publish(SINGLE_ALLELE, true);
    let accounting = &measured.manifest["repeat_loci"]["tract_output"];
    let total = accounting["realized_total_tract_bases"]
        .as_u64()
        .expect("realizado");
    assert_eq!(accounting["realized_truth_tract_bases"], 900);
    assert_eq!(accounting["realized_vcf_tract_bases"], 190);
    assert_eq!(total, 1090);
    assert!(
        total > TRUTH_BOUND,
        "control: lo realizado supera la cota de la verdad"
    );
    fs::remove_dir_all(&measured.root).ok();

    let (tight_root, tight) = scenario(SINGLE_ALLELE, Some(total - 1), true);
    let effective = load_and_validate_config(&tight).expect("la verdad sola cabe en el tope");
    let destination = tight_root.join("published");
    let error = execute_validated_config(&effective, &destination)
        .expect_err("la verdad más lo medido en el VCF no cabe");
    assert_eq!(error.code(), "GPS-CONF-008");
    let detail = error.context().get("detail").cloned().unwrap_or_default();
    assert!(
        detail.contains(&total.to_string())
            && detail.contains(&(total - 1).to_string())
            && detail.contains("STR-1"),
        "el rechazo tiene que nombrar el tope, la cifra y el locus: {detail}"
    );
    assert!(
        !destination.join("outputs/truth.vcf").exists(),
        "el VCF no llega a escribirse"
    );
    let failed: Value = serde_json::from_slice(
        &fs::read(destination.join("manifest.json")).expect("manifiesto del paquete fallido"),
    )
    .expect("JSON");
    assert_eq!(failed["state"], "failed", "el paquete fallido es coherente");
    fs::remove_dir_all(&tight_root).ok();

    let (exact_root, exact) = scenario(SINGLE_ALLELE, Some(total), true);
    let effective = load_and_validate_config(&exact).expect("un tope igual a la cifra admite");
    execute_validated_config(&effective, &exact_root.join("published"))
        .expect("y publica, porque lo medido cabe justo");
    fs::remove_dir_all(&exact_root).ok();
}

/// **Un escenario sin VCF sigue cuadrando: lo realizado en el VCF es cero.**
#[test]
fn a_run_without_a_vcf_reports_zero_realized_vcf_bases() {
    let published = publish(THREE_UNEVEN_ALLELES, false);
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

/// Dos loci en la misma secuencia, con figuras distintas: `(CA)50` en
/// 101-200 y `(AT)30` en 301-360, separados por cien `G`.
fn two_loci_reference() -> String {
    format!(
        "{}{}{}{}{}",
        "G".repeat(100),
        "CA".repeat(50),
        "G".repeat(100),
        "AT".repeat(30),
        "G".repeat(100)
    )
}

const TWO_LOCI: &str = r#"
[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 200
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 45
frequency = 1.0

[[reference.repeat_loci]]
id = "STR-2"
sequence = "chr1"
start = 301
end = 360
motif = "AT"

[[reference.repeat_loci.alleles]]
units = 20
frequency = 1.0
"#;

/// El juez de contabilidad con una configuración efectiva dada en bytes y el
/// VCF publicado con el nombre dado.
fn judge_with(published: &Published, effective: &[u8], vcf_name: &str) -> Result<(), String> {
    let effective_path = published.destination.join("effective-config.json");
    let manifest_path = published.destination.join("manifest.json");
    let vcf_path = published.destination.join("outputs").join(vcf_name);
    validate_published_tract_output_accounting(
        &effective_path,
        effective,
        &manifest_path,
        &fs::read(&manifest_path).expect("manifiesto"),
        &published.destination.join("outputs/truth.jsonl"),
        Some(&vcf_path),
    )
    .map_err(|error| format!("{error:?}"))
}

fn per_locus_vcf_figures(effective: &Value) -> Vec<(String, u64)> {
    effective["repeat_output_budget"]["per_locus"]
        .as_array()
        .expect("reparto por locus")
        .iter()
        .map(|locus| {
            (
                locus["id"].as_str().expect("id").to_owned(),
                locus["vcf_tract_bases_bound"].as_u64().expect("cifra"),
            )
        })
        .collect()
}

/// **Cada locus lleva las bases de tracto de sus propios registros.**
///
/// *Qué hace que pueda fallar:* dos loci en la misma secuencia con figuras
/// distintas. Sus cifras se cuentan aquí releyendo el VCF y atribuyendo cada
/// registro de repetición al locus cuyo intervalo —escrito a mano, no leído
/// del manifiesto— contiene. Un productor que lo atribuyera todo al primer
/// locus de la secuencia conservaría la suma y fallaría la igualdad por
/// locus; y el juez tiene que rechazar las dos manipulaciones que conservan
/// la suma: cambiar las cifras de sitio y darle todo al primero.
#[test]
fn each_locus_carries_the_vcf_tract_bases_of_its_own_records() {
    let published = publish_on(&two_loci_reference(), TWO_LOCI, Some("text"));
    judge(&published, true).expect("el juez acepta el reparto publicado");

    let vcf = fs::read_to_string(published.destination.join("outputs/truth.vcf")).expect("VCF");
    let mut counted = vec![("STR-1".to_owned(), 0_u64), ("STR-2".to_owned(), 0_u64)];
    for line in vcf.lines().filter(|line| !line.starts_with('#')) {
        let fields = line.split('\t').collect::<Vec<_>>();
        if !fields[7].split(';').any(|field| field.starts_with("RN=")) {
            continue;
        }
        let position: u64 = fields[1].parse().expect("POS");
        let last = position + fields[3].len() as u64 - 1;
        let bases =
            fields[3].len() as u64 + fields[4].split(',').map(|a| a.len() as u64).sum::<u64>();
        let index = [(101, 200), (301, 360)]
            .iter()
            .position(|(start, end)| position <= *start && *end <= last)
            .expect("cada registro de repetición contiene un locus");
        counted[index].1 += bases;
    }
    let published_figures = per_locus_vcf_figures(&published.effective);
    assert_eq!(
        published_figures, counted,
        "el reparto publicado no es el de los registros"
    );
    assert!(
        counted[0].1 > 0 && counted[1].1 > 0 && counted[0].1 != counted[1].1,
        "control: dos cifras distintas y no nulas, para que cambiarlas se note: {counted:?}"
    );

    let effective_bytes =
        fs::read(published.destination.join("effective-config.json")).expect("efectiva");
    let total = counted[0].1 + counted[1].1;
    for (reason, first, second) in [
        ("cifras cambiadas de sitio", counted[1].1, counted[0].1),
        ("todo al primer locus", total, 0),
    ] {
        let mut document: Value = serde_json::from_slice(&effective_bytes).expect("JSON");
        let loci = document["repeat_output_budget"]["per_locus"]
            .as_array_mut()
            .expect("reparto");
        loci[0]["vcf_tract_bases_bound"] = serde_json::json!(first);
        loci[1]["vcf_tract_bases_bound"] = serde_json::json!(second);
        let tampered = serde_json::to_vec(&document).expect("JSON");
        let error = judge_with(&published, &tampered, "truth.vcf")
            .expect_err(&format!("el juez tiene que rechazar «{reason}»"));
        assert!(
            error.contains("del locus STR-1 no son las que llevan sus registros"),
            "«{reason}» tiene que caer por el reparto, no por otra comprobación: {error}"
        );
    }
    fs::remove_dir_all(&published.root).ok();
}

/// **El juez recuenta también un VCF comprimido en BGZF.**
///
/// *Qué hace que pueda fallar:* el mismo escenario publicado en texto y en
/// BGZF tiene que dar la misma cifra del VCF, y el juez tiene que aceptar el
/// comprimido y rechazar su configuración efectiva con la cifra alterada. Un
/// juez que no reconociera el `.gz` leería bytes comprimidos como texto y no
/// aceptaría el paquete intacto.
#[test]
fn the_judge_also_recounts_a_bgzf_vcf() {
    let text = publish(SINGLE_ALLELE, true);
    let compressed = publish_on(&reference_bases(), SINGLE_ALLELE, Some("bgzf"));
    assert!(
        compressed.destination.join("outputs/truth.vcf.gz").exists(),
        "control: el VCF se publicó comprimido"
    );
    let figure = |published: &Published| {
        published.manifest["repeat_loci"]["tract_output"]["realized_vcf_tract_bases"].as_u64()
    };
    assert_eq!(figure(&compressed), Some(190));
    assert_eq!(figure(&compressed), figure(&text));

    let effective_bytes =
        fs::read(compressed.destination.join("effective-config.json")).expect("efectiva");
    judge_with(&compressed, &effective_bytes, "truth.vcf.gz")
        .expect("el juez acepta el paquete comprimido intacto");
    let mut document: Value = serde_json::from_slice(&effective_bytes).expect("JSON");
    document["repeat_output_budget"]["vcf_tract_bases_bound"] = serde_json::json!(191);
    let error = judge_with(
        &compressed,
        &serde_json::to_vec(&document).expect("JSON"),
        "truth.vcf.gz",
    )
    .expect_err("una cifra alterada se rechaza también con BGZF");
    assert!(
        error.contains("no publica las bases de tracto que el VCF lleva"),
        "{error}"
    );
    fs::remove_dir_all(&text.root).ok();
    fs::remove_dir_all(&compressed.root).ok();
}
