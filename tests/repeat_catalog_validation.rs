//! Entrega 4 del plan de detección de microsatélites: **el catálogo aparte**
//! (F6.4.a.14, REQ-MAN-98) y **su juez independiente** (F6.4.a.15,
//! REQ-VAL-64, `src/repeat_catalog_validation.rs`).
//!
//! # Cómo se prueba el juez
//!
//! Se publica **de verdad** un paquete con la sección de microsatélites sobre
//! una referencia con tractos plantados en posiciones escritas aquí, y cada
//! caso **altera a mano** la tabla o el metadato publicados (o el enlace, el
//! manifiesto o la efectiva) en un solo punto. Salvo en los casos que prueban
//! las huellas, se rehacen las huellas que encadenan los ficheros, de modo que
//! lo único distinto es lo que el caso altera y el rechazo tiene que ser el de
//! esa comprobación: se exige el motivo exacto. Nada de lo esperado se calcula
//! con el detector ni con la conversión.
//!
//! El juez se llama directamente, con lo que el paquete publica (la efectiva,
//! la lista de artefactos del manifiesto, el FASTA final y la matriz),
//! leídos aquí sin el productor. Así cada comprobación es alcanzable aunque en
//! el juez completo otra la tape antes. Tres pruebas pasan además por el juez
//! completo (`validate_published_resources`): que el juez está encadenado, que
//! la reconstrucción con el productor ata lo que el juez independiente no
//! puede ver, y que sin la sección no se publica nada nuevo.
//!
//! # El tracto de pureza real 0,79
//!
//! `GATGATGATGAC × 6` con su flanco da un intervalo del detector de pureza
//! 0,792 (57/72 bases en copias completas) que supera todos los demás
//! criterios. El productor lo descarta por pureza; una fila suya añadida a mano
//! la rechaza el juez por pureza. La prueba de independencia del plan —el
//! productor con la medida de pureza estropeada lo publica y el juez lo
//! rechaza— se hizo en una copia del árbol y está registrada en el plan.
//!
//! # Lo que estas pruebas no ven, dicho sin adornos
//!
//! La exhaustividad: el juez no puede ver un tramo que falte, y ninguna prueba
//! de aquí finge lo contrario. La reconstrucción del juez completo sí lo ve, y
//! la prueba `the_reconstruction_ties_what_the_judge_cannot_see` lo demuestra.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::diagnostics::Diagnostic;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::repeat_catalog_validation::{
    PublishedArtifact, PublishedCatalogInputs, PublishedSequence, validate_published_repeat_catalog,
};
use genome_population_simulator::results::{json_document_bytes, validate_published_resources};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};

static NONCE: AtomicU64 = AtomicU64::new(0);
const SEED: &str = "00000000000000000000000000000000000000000000000000000000000000d1";

fn temporary(label: &str) -> PathBuf {
    let path = std::env::temp_dir().join(format!(
        "gps-repeat-catalog-{label}-{}-{}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ));
    fs::create_dir_all(&path).expect("directorio temporal");
    path
}

fn sha(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

// ---------------------------------------------------------------------------
// La referencia y el paquete publicado
// ---------------------------------------------------------------------------

/// Fondo sin repeticiones de periodo 1 ni 2: cada base distinta de las dos
/// anteriores (xorshift64).
fn background(seed: u64, length: usize) -> Vec<u8> {
    let mut state = seed;
    let mut bases: Vec<u8> = Vec::with_capacity(length);
    while bases.len() < length {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        let candidate = b"ACGT"[(state % 4) as usize];
        let n = bases.len();
        if (n >= 1 && bases[n - 1] == candidate) || (n >= 2 && bases[n - 2] == candidate) {
            continue;
        }
        bases.push(candidate);
    }
    bases
}

/// Escribe `text` desde la base `start` (base uno).
fn put(bases: &mut [u8], start: usize, text: &str) {
    bases[start - 1..start - 1 + text.len()].copy_from_slice(text.as_bytes());
}

/// Planta `motif × units` desde `start` con un flanco a cada lado que no lo
/// prolonga.
fn plant(bases: &mut [u8], start: usize, motif: &str, units: usize) {
    let flank = *b"ACGT"
        .iter()
        .find(|base| **base != motif.as_bytes()[0] && **base != motif.as_bytes()[motif.len() - 1])
        .expect("una base libre");
    put(bases, start, &motif.repeat(units));
    bases[start - 2] = flank;
    bases[start - 1 + motif.len() * units] = flank;
}

/// El flanco y el tracto de pureza 0,79 (ver la cabecera).
const IMPURE_LEFT: &str = "GTCAGTCTGATCGTAGCTGATGCATCGTGCTAGT";
const IMPURE_RIGHT: &str = "TGCATGACTGATCGATGCTAGTCGATCGTAGCTG";
/// Base uno de la primera base de `IMPURE_LEFT` en `chr2`.
const IMPURE_AT: usize = 1967;

/// Dos secuencias con todo lo que las pruebas necesitan, en posiciones
/// escritas aquí:
///
/// - `chr1`: `A×12` en 1..12 (pegado al borde: se descarta), `A×12` en
///   301..312 (detectado), `(CA)×10` en 801..820 (detectado y retirado por la
///   fila de la matriz en 810), `(GAT)×7` en 1301..1321 (detectado),
///   `(AAAG)×6` en 1801..1824 (lo declara el locus `STR-AAAG`), `(TG)×5` en
///   2401..2410 (por debajo del suelo de 6) y `(CA)×8` en 2601..2616 pegado a
///   once `A` (2617..2627), que el detector da como dos tramos que se tocan
///   (el de las `A`, retirado por la fila de la matriz en 2622);
/// - `chr2`: `(ACGTT)×5` en 401..425 (detectado), `A×210` en 1001..1210 (por
///   encima del tope de 200 unidades), `(CACAG)×3 CACA` en 1601..1619 (tirada
///   de 2 copias; el detector no lo ve como intervalo), el tracto de pureza
///   0,79 en 2000..2071, `(GT)×8` en 2401..2416 (fuera de las regiones),
///   `A×12` en 2601..2612 con una `N` delante, `A×12` en 2701..2712 con una
///   `N` detrás y `A×12` en 2989..3000, que acaba en la última base.
fn reference() -> Vec<(&'static str, Vec<u8>)> {
    let mut chr1 = background(0x5DEE_CE66_D1CE_4E5B, 3000);
    put(&mut chr1, 1, "AAAAAAAAAAAA");
    chr1[12] = b'C';
    plant(&mut chr1, 301, "A", 12);
    plant(&mut chr1, 801, "CA", 10);
    plant(&mut chr1, 1301, "GAT", 7);
    plant(&mut chr1, 1801, "AAAG", 6);
    plant(&mut chr1, 2401, "TG", 5);
    put(
        &mut chr1,
        2600,
        &format!("G{}{}G", "CA".repeat(8), "A".repeat(11)),
    );
    let mut chr2 = background(0x2545_F491_4F6C_DD1D, 3000);
    plant(&mut chr2, 401, "ACGTT", 5);
    plant(&mut chr2, 1001, "A", 210);
    put(&mut chr2, 1600, "TCACAGCACAGCACAGCACAT");
    put(
        &mut chr2,
        IMPURE_AT,
        &format!("{IMPURE_LEFT}{}{IMPURE_RIGHT}", "GATGATGATGAC".repeat(6)),
    );
    plant(&mut chr2, 2401, "GT", 8);
    put(&mut chr2, 2600, &format!("N{}C", "A".repeat(12)));
    put(&mut chr2, 2700, &format!("C{}N", "A".repeat(12)));
    put(&mut chr2, 2988, &format!("C{}", "A".repeat(12)));
    vec![("chr1", chr1), ("chr2", chr2)]
}

fn write_fasta(path: &Path, sequences: &[(&str, Vec<u8>)]) {
    let mut text = String::new();
    for (id, bases) in sequences {
        text.push_str(&format!(">{id}\n"));
        for chunk in bases.chunks(60) {
            text.push_str(std::str::from_utf8(chunk).expect("ASCII"));
            text.push('\n');
        }
    }
    fs::write(path, text).expect("FASTA");
}

/// Publica el escenario: con la sección de microsatélites (y sus regiones) si
/// `detection`, sin ella si no. La matriz suministrada lleva una fila dentro
/// del `(CA)×10` y otra fuera de todo tracto.
fn publish(label: &str, detection: bool) -> PathBuf {
    let root = temporary(label);
    let sequences = reference();
    write_fasta(&root.join("reference.fa"), &sequences);
    let mut matrix = String::from("contig\tposition\tref\talt\tpopulation_id\talt_af\n");
    for position in [810_usize, 2000, 2622] {
        let reference = sequences[0].1[position - 1];
        let alternate = *b"ACGT"
            .iter()
            .find(|base| **base != reference)
            .expect("alt");
        matrix.push_str(&format!(
            "chr1\t{position}\t{}\t{}\tP1\t0.5\n",
            char::from(reference),
            char::from(alternate)
        ));
    }
    fs::write(root.join("snvs.tsv"), matrix).expect("matriz");
    let section = if detection {
        "[reference.repeat_detection]\n\n\
         [reference.repeat_detection.slippage]\nmodel = \"stepwise\"\nmean_events = 1.0\n\n\
         [[reference.repeat_detection.regions]]\nsequence = \"chr1\"\nstart = 1\nend = 3000\n\n\
         [[reference.repeat_detection.regions]]\nsequence = \"chr2\"\nstart = 1\nend = 2300\n\n\
         [[reference.repeat_detection.regions]]\nsequence = \"chr2\"\nstart = 2500\nend = 3000\n\n\
         [[reference.repeat_loci]]\nid = \"STR-AAAG\"\nsequence = \"chr1\"\nstart = 1801\n\
         end = 1824\nmotif = \"AAAG\"\n\n"
    } else {
        ""
    };
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 3000\n\n\
             [[reference.sequences]]\nid = \"chr2\"\nclass = \"chromosome\"\nlength = 3000\n\n\
             {section}\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n"
        ),
    )
    .expect("configuración");
    let destination = root.join("published");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &destination,
    )
    .unwrap_or_else(|error| panic!("publicación: {}", error.render_human()));
    destination
}

/// Un paquete publicado para un solo caso, que se borra al soltarlo. Cada
/// caso publica el suyo (unas centésimas de segundo): así no queda ningún
/// paquete compartido en el directorio temporal al acabar.
struct Package {
    root: PathBuf,
}

impl Drop for Package {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.root);
    }
}

impl Package {
    fn copy() -> Self {
        let destination = publish("case", true);
        Self {
            root: destination.parent().expect("raíz").to_path_buf(),
        }
    }

    fn dir(&self) -> PathBuf {
        self.root.join("published")
    }

    fn json(&self, relative: &str) -> Value {
        serde_json::from_slice(&fs::read(self.dir().join(relative)).expect("JSON")).expect("JSON")
    }

    fn text(&self, relative: &str) -> String {
        fs::read_to_string(self.dir().join(relative)).expect("texto")
    }

    fn write_json(&self, relative: &str, value: &Value) -> Vec<u8> {
        let bytes = json_document_bytes(value).expect("JSON");
        fs::write(self.dir().join(relative), &bytes).expect("escritura");
        bytes
    }

    /// Escribe la tabla y el metadato dados y **rehace las huellas** que los
    /// encadenan: la de la tabla en el metadato, la del metadato en el enlace
    /// de la efectiva, y las de los dos en el manifiesto.
    fn republish(&self, table: &str, metadata: &Value, effective: &Value) {
        fs::write(self.dir().join("repeat_catalog.tsv"), table).expect("tabla");
        let mut metadata = metadata.clone();
        if metadata.is_object() {
            metadata["table"]["sha256"] = json!(sha(table.as_bytes()));
        }
        let metadata_bytes = self.write_json("repeat_catalog.metadata.json", &metadata);
        let mut effective = effective.clone();
        if effective.get("repeat_detection").is_some() {
            effective["repeat_detection"]["catalog"]["metadata_sha256"] =
                json!(sha(&metadata_bytes));
        }
        let effective_bytes = self.write_json("effective-config.json", &effective);
        let mut manifest = self.json("manifest.json");
        for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
            let bytes = match artifact["path"].as_str() {
                Some("repeat_catalog.tsv") => table.as_bytes().to_vec(),
                Some("repeat_catalog.metadata.json") => metadata_bytes.clone(),
                Some("effective-config.json") => effective_bytes.clone(),
                _ => continue,
            };
            artifact["checksum"]["digest"] = json!(sha(&bytes));
            artifact["size_bytes"] = json!(bytes.len());
        }
        self.write_json("manifest.json", &manifest);
    }

    /// El juez del catálogo sobre lo que hay en disco, con las entradas leídas
    /// aquí: la efectiva, la lista de artefactos del manifiesto, el FASTA final
    /// y la matriz publicados.
    fn judge(&self) -> Result<(), Diagnostic> {
        let effective = self.json("effective-config.json");
        let manifest = self.json("manifest.json");
        let artifacts: Vec<(String, String, String, u64, String)> = manifest["artifacts"]
            .as_array()
            .expect("artefactos")
            .iter()
            .map(|artifact| {
                (
                    artifact["id"].as_str().expect("id").to_owned(),
                    artifact["kind"].as_str().expect("kind").to_owned(),
                    artifact["path"].as_str().expect("path").to_owned(),
                    artifact["size_bytes"].as_u64().expect("size_bytes"),
                    artifact["checksum"]["digest"]
                        .as_str()
                        .expect("digest")
                        .to_owned(),
                )
            })
            .collect();
        let artifacts: Vec<PublishedArtifact<'_>> = artifacts
            .iter()
            .map(|(id, kind, path, size_bytes, digest)| PublishedArtifact {
                id,
                kind,
                path,
                size_bytes: *size_bytes,
                digest,
            })
            .collect();
        let fasta = self.text("resources/references/ref-main.fa");
        let mut sequences: Vec<(String, Vec<u8>)> = Vec::new();
        for line in fasta.lines() {
            if let Some(header) = line.strip_prefix('>') {
                sequences.push((
                    header.split_whitespace().next().expect("id").to_owned(),
                    Vec::new(),
                ));
            } else {
                sequences
                    .last_mut()
                    .expect("cabecera")
                    .1
                    .extend(line.bytes());
            }
        }
        let sequences: Vec<PublishedSequence<'_>> = sequences
            .iter()
            .map(|(id, bases)| PublishedSequence { id, bases })
            .collect();
        let matrix = self.text("population_snvs.tsv");
        let rows: Vec<(&str, u64, &str)> = matrix
            .lines()
            .skip(1)
            .map(|line| {
                let fields: Vec<&str> = line.split('\t').collect();
                (fields[1], fields[2].parse().expect("posición"), fields[0])
            })
            .collect();
        validate_published_repeat_catalog(&PublishedCatalogInputs {
            root: &self.dir(),
            effective: effective.as_object().expect("objeto"),
            artifacts: &artifacts,
            reference_id: "ref-main",
            sequences: &sequences,
            matrix_variants: &rows,
        })
    }

    /// El juez completo de los recursos publicados.
    fn judge_resources(&self) -> Result<(), Diagnostic> {
        let effective_path = self.dir().join("effective-config.json");
        let manifest_path = self.dir().join("manifest.json");
        validate_published_resources(
            &self.dir(),
            &effective_path,
            &fs::read(&effective_path).expect("efectiva"),
            &manifest_path,
            &fs::read(&manifest_path).expect("manifiesto"),
            None,
        )
    }
}

fn detail(result: Result<(), Diagnostic>) -> String {
    match result {
        Ok(()) => "aceptado".to_owned(),
        Err(error) => error
            .context()
            .get("detail")
            .cloned()
            .unwrap_or_else(|| error.render_human()),
    }
}

/// La tabla publicada del paquete, escrita a mano: `STR-AAAG` es el
/// declarado; el `(CA)×10` y las once `A` faltan porque se retiraron. Las
/// puntuaciones son +2 por cada base igual a la de un periodo antes: 11, 18,
/// 14 y 20 coincidencias.
const TABLE: &str = "locus_id\tsequence_id\tstart\tend\tmotif\tscore\torigin\n\
SSR-1-301-312-A\tchr1\t301\t312\tA\t22\tdetected\n\
SSR-1-1301-1321-GAT\tchr1\t1301\t1321\tGAT\t36\tdetected\n\
STR-AAAG\tchr1\t1801\t1824\tAAAG\t.\tdeclared\n\
SSR-1-2601-2616-CA\tchr1\t2601\t2616\tCA\t28\tdetected\n\
SSR-2-401-425-ACGTT\tchr2\t401\t425\tACGTT\t40\tdetected\n";

/// Un caso: qué se altera y el motivo exacto del rechazo.
struct Case {
    label: &'static str,
    expected: &'static str,
    table: fn(&str) -> String,
    metadata: fn(&mut Value),
    effective: fn(&mut Value),
}

fn same(table: &str) -> String {
    table.to_owned()
}

fn untouched(_: &mut Value) {}

/// Cambia el texto `from` por `to` en la tabla, exigiendo que esté.
fn swap(table: &str, from: &str, to: &str) -> String {
    assert!(table.contains(from), "la tabla no contiene {from:?}");
    table.replacen(from, to, 1)
}

/// Añade una fila de un detectado antes de la del declarado o al final.
fn with_row(table: &str, row: &str) -> String {
    format!("{table}{row}\n")
}

fn run_cases(cases: &[Case]) {
    for case in cases {
        let package = Package::copy();
        let table = package.text("repeat_catalog.tsv");
        assert_eq!(table, TABLE, "la tabla publicada no es la escrita a mano");
        let mut metadata = package.json("repeat_catalog.metadata.json");
        let mut effective = package.json("effective-config.json");
        let altered = (case.table)(&table);
        let before = (altered.clone(), metadata.clone(), effective.clone());
        (case.metadata)(&mut metadata);
        (case.effective)(&mut effective);
        assert!(
            before != (table.clone(), metadata.clone(), effective.clone()),
            "«{}»: la alteración no alteró nada",
            case.label
        );
        package.republish(&altered, &metadata, &effective);
        assert_eq!(
            detail(package.judge()),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

// ---------------------------------------------------------------------------
// F6.4.a.14: lo que se publica
// ---------------------------------------------------------------------------

/// **El catálogo publicado es el de los tractos plantados**, escrito a mano;
/// su metadato pasa su esquema por el validador externo de Python y la efectiva
/// el suyo; y **los recuentos suman lo que el barrido encontró**, comprobado
/// aquí con aritmética propia.
#[test]
fn the_published_catalog_is_the_planted_one_and_its_counts_add_up() {
    let package = Package::copy();
    assert_eq!(package.text("repeat_catalog.tsv"), TABLE);
    package
        .judge()
        .expect("control positivo: el paquete intacto se acepta");
    package
        .judge_resources()
        .expect("control positivo: el juez completo acepta el paquete intacto");
    let metadata = package.json("repeat_catalog.metadata.json");
    let effective = package.json("effective-config.json");
    let metadata_bytes =
        fs::read(package.dir().join("repeat_catalog.metadata.json")).expect("bytes");
    assert_eq!(
        effective["repeat_detection"]["catalog"],
        json!({"metadata_path": "repeat_catalog.metadata.json", "metadata_sha256": sha(&metadata_bytes)})
    );
    assert!(effective["repeat_detection"].get("retired_loci").is_none());
    assert_eq!(metadata["table"]["sha256"], json!(sha(TABLE.as_bytes())));
    let fasta = fs::read(package.dir().join("resources/references/ref-main.fa")).expect("FASTA");
    assert_eq!(metadata["reference"]["fasta_sha256"], json!(sha(&fasta)));
    // Lo que se espera de cada longitud de motivo, escrito a mano.
    let expected: [(u64, u64, u64, u64); 6] = [
        // (publicados, sustituidos, retirados, sin copias)
        (1, 0, 1, 0),
        (1, 0, 1, 0),
        (1, 0, 0, 0),
        (0, 1, 0, 0),
        (1, 0, 0, 0),
        (0, 0, 0, 0),
    ];
    let counts = metadata["counts_by_motif_length"]
        .as_array()
        .expect("recuentos");
    for (index, entry) in counts.iter().enumerate() {
        let read = |key: &str| entry[key].as_u64().expect("entero");
        let discarded: u64 = entry["discarded"]
            .as_object()
            .expect("descartes")
            .values()
            .map(|value| value.as_u64().expect("entero"))
            .sum();
        assert_eq!(read("motif_length"), index as u64 + 1);
        assert_eq!(
            read("detector_intervals"),
            discarded + read("accepted"),
            "{entry}"
        );
        assert_eq!(
            read("accepted"),
            read("published")
                + read("replaced_by_declared")
                + read("retired_by_variant")
                + read("without_haplotype_copies"),
            "{entry}"
        );
        assert_eq!(
            (
                read("published"),
                read("replaced_by_declared"),
                read("retired_by_variant"),
                read("without_haplotype_copies")
            ),
            expected[index],
            "{entry}"
        );
    }
    // Cada tracto plantado que no entra cae en su etapa.
    let stage = |length: usize, name: &str| counts[length - 1]["discarded"][name].as_u64();
    assert!(stage(1, "edge") >= Some(1), "el A×12 pegado al borde");
    assert!(stage(1, "max_units") >= Some(1), "el A×210");
    assert!(stage(2, "min_copies") >= Some(1), "el (TG)×5");
    // El `(CACAG)×3 CACA` no llega a ser un intervalo del detector (sus
    // discrepancias pesan más que sus coincidencias), así que no se cuenta en
    // ninguna etapa: solo sirve para la fila añadida a mano de la prueba de la
    // tirada.
    assert!(
        stage(3, "min_purity") >= Some(1),
        "el tracto de pureza 0,79"
    );
    assert!(
        stage(2, "outside_regions") >= Some(1),
        "el (GT)×8 fuera de las regiones"
    );
    assert_eq!(
        metadata["replaced_tracts"],
        json!([{"sequence_id": "chr1", "start": 1801, "end": 1824, "motif": "AAAG", "score": 40}])
    );
    let retired = &metadata["retired_tracts"];
    assert_eq!(retired.as_array().map(Vec::len), Some(2), "{retired}");
    assert_eq!(retired[1]["locus_id"], "SSR-1-2617-2627-A");
    assert_eq!(retired[1]["variant"]["position"], 2622);
    assert_eq!(retired[0]["locus_id"], "SSR-1-801-820-CA");
    assert_eq!(retired[0]["score"], 36);
    assert_eq!(retired[0]["variant"]["source"], "population_snvs");
    assert_eq!(retired[0]["variant"]["position"], 810);

    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    for (document, schema, accepted) in [
        (
            metadata.clone(),
            "repeat-catalog-metadata-v1.schema.json",
            true,
        ),
        (effective.clone(), "effective-config-v6.schema.json", true),
        (
            {
                let mut document = metadata.clone();
                document["unexpected"] = json!(1);
                document
            },
            "repeat-catalog-metadata-v1.schema.json",
            false,
        ),
        (
            {
                let mut document = effective.clone();
                document["repeat_detection"]["retired_loci"] = json!([]);
                document
            },
            "effective-config-v6.schema.json",
            false,
        ),
    ] {
        let path = package.root.join("instance.json");
        fs::write(&path, json_document_bytes(&document).expect("JSON")).expect("instancia");
        let status = Command::new("python3")
            .arg(root.join("tests/external_jsonschema_validator.py"))
            .arg(root.join("schemas").join(schema))
            .arg(&path)
            .status()
            .expect("validador de esquemas en Python");
        assert_eq!(
            status.success(),
            accepted,
            "{schema}: aceptado = {accepted}"
        );
    }
}

/// **Sin la sección no se publica nada nuevo**: ni la tabla, ni el metadato,
/// ni la sección en la efectiva (la huella de todo lo demás la vigila
/// `tests/repeat_detection_frozen_identity.rs`). Y el juez rechaza un catálogo
/// listado en el manifiesto de una ejecución sin la sección.
#[test]
fn without_the_section_nothing_new_is_published() {
    let destination = publish("without", false);
    assert!(!destination.join("repeat_catalog.tsv").exists());
    assert!(!destination.join("repeat_catalog.metadata.json").exists());
    let effective: Value = serde_json::from_slice(
        &fs::read(destination.join("effective-config.json")).expect("efectiva"),
    )
    .expect("JSON");
    assert!(effective.get("repeat_detection").is_none());
    let package = Package {
        root: destination.parent().expect("raíz").to_path_buf(),
    };
    package
        .judge()
        .expect("control positivo: sin sección y sin catálogo se acepta");
    let mut manifest = package.json("manifest.json");
    manifest["artifacts"]
        .as_array_mut()
        .expect("artefactos")
        .push(json!({
            "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)},
            "id": "repeat-catalog-tsv",
            "kind": "repeat_catalog_tsv",
            "path": "repeat_catalog.tsv",
            "size_bytes": 0
        }));
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge()),
        "hay un catálogo publicado y la configuración efectiva no declara la sección de microsatélites"
    );
}

// ---------------------------------------------------------------------------
// F6.4.a.15: cada comprobación del juez, disparada con lo publicado alterado
// ---------------------------------------------------------------------------

/// **Las huellas y el enlace.** Aquí las huellas **no** se rehacen después de
/// la alteración que las rompe: es lo que se prueba.
#[test]
fn the_identity_of_the_catalog_is_judged() {
    // Enlace de la efectiva.
    let package = Package::copy();
    let mut effective = package.json("effective-config.json");
    effective["repeat_detection"]["catalog"]["metadata_path"] = json!("otro.json");
    package.write_json("effective-config.json", &effective);
    assert_eq!(
        detail(package.judge()),
        "repeat_detection.catalog.metadata_path no es la ruta del metadato"
    );
    let package = Package::copy();
    let mut effective = package.json("effective-config.json");
    effective["repeat_detection"]["catalog"]["metadata_sha256"] = json!(5);
    package.write_json("effective-config.json", &effective);
    assert_eq!(
        detail(package.judge()),
        "catalog.metadata_sha256 no es un texto"
    );
    let package = Package::copy();
    let mut effective = package.json("effective-config.json");
    effective["repeat_detection"]["catalog"]["metadata_sha256"] = json!("0".repeat(64));
    package.write_json("effective-config.json", &effective);
    assert_eq!(
        detail(package.judge()),
        "la huella del metadato no es la del enlace de la efectiva o la del manifiesto"
    );
    let package = Package::copy();
    let mut effective = package.json("effective-config.json");
    effective["repeat_detection"]
        .as_object_mut()
        .expect("sección")
        .remove("catalog");
    package.write_json("effective-config.json", &effective);
    assert_eq!(
        detail(package.judge()),
        "repeat_detection.catalog no es un objeto"
    );
    // El manifiesto.
    let package = Package::copy();
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "repeat-catalog-tsv" {
            artifact["path"] = json!("otra.tsv");
        }
    }
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge()),
        "el manifiesto no lista repeat-catalog-tsv con clase repeat_catalog_tsv y ruta repeat_catalog.tsv"
    );
    let package = Package::copy();
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "repeat-catalog-tsv" {
            artifact["checksum"]["digest"] = json!("0".repeat(64));
        }
    }
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge()),
        "la huella de la tabla no es la del manifiesto"
    );
    let package = Package::copy();
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "repeat-catalog-metadata" {
            artifact["checksum"]["digest"] = json!("0".repeat(64));
        }
    }
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge()),
        "la huella del metadato no es la del enlace de la efectiva o la del manifiesto"
    );
    // Los ficheros.
    let package = Package::copy();
    fs::remove_file(package.dir().join("repeat_catalog.tsv")).expect("borrar");
    assert!(
        detail(package.judge()).starts_with("no se puede leer: "),
        "{}",
        detail(package.judge())
    );
}

/// **Lo que el metadato dice de sí mismo, de la referencia y de los
/// criterios**, y su lectura.
#[test]
fn the_metadata_is_judged() {
    // El metadato ilegible se escribe a mano, sin pasar por JSON.
    let package = Package::copy();
    let effective = package.json("effective-config.json");
    package.republish(TABLE, &json!({}), &effective);
    fs::write(
        package.dir().join("repeat_catalog.metadata.json"),
        b"no es JSON",
    )
    .expect("escritura");
    let bytes = b"no es JSON";
    let mut effective = effective.clone();
    effective["repeat_detection"]["catalog"]["metadata_sha256"] = json!(sha(bytes));
    package.write_json("effective-config.json", &effective);
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "repeat-catalog-metadata" {
            artifact["checksum"]["digest"] = json!(sha(bytes));
            artifact["size_bytes"] = json!(bytes.len());
        }
    }
    package.write_json("manifest.json", &manifest);
    assert!(
        detail(package.judge()).starts_with("JSON no legible: "),
        "{}",
        detail(package.judge())
    );
    run_cases(&[
        Case {
            label: "un metadato que no es un objeto",
            expected: "metadato no es un objeto",
            table: same,
            metadata: |m| *m = json!([1]),
            effective: untouched,
        },
        Case {
            label: "otra clase de documento",
            expected: "la clase, la versión, el sistema de coordenadas, la versión del algoritmo o la tabla (ruta y huella) del metadato no son los esperados",
            table: same,
            metadata: |m| m["document_kind"] = json!("otro"),
            effective: untouched,
        },
        Case {
            label: "otra ruta de la tabla",
            expected: "la clase, la versión, el sistema de coordenadas, la versión del algoritmo o la tabla (ruta y huella) del metadato no son los esperados",
            table: same,
            metadata: |m| m["table"]["path"] = json!("otra.tsv"),
            effective: untouched,
        },
        Case {
            label: "otra huella del FASTA final",
            expected: "la referencia del metadato no es la referencia final publicada (identificador y huella del FASTA)",
            table: same,
            metadata: |m| m["reference"]["fasta_sha256"] = json!("0".repeat(64)),
            effective: untouched,
        },
        Case {
            label: "otro identificador de referencia",
            expected: "la referencia del metadato no es la referencia final publicada (identificador y huella del FASTA)",
            table: same,
            metadata: |m| m["reference"]["reference_id"] = json!("ref-other"),
            effective: untouched,
        },
        Case {
            label: "otro tope en el metadato",
            expected: "los criterios del metadato no son los de la configuración efectiva",
            table: same,
            metadata: |m| m["criteria"]["max_units"] = json!(199),
            effective: untouched,
        },
        Case {
            label: "otra pureza mínima en la efectiva",
            expected: "los criterios del metadato no son los de la configuración efectiva",
            table: same,
            metadata: untouched,
            effective: |e| e["repeat_admission"]["min_purity"] = json!(0.7),
        },
        Case {
            label: "una pureza mínima que no es un número, en los dos",
            expected: "criteria.min_purity no es un número",
            table: same,
            metadata: |m| m["criteria"]["min_purity"] = json!("x"),
            effective: |e| e["repeat_admission"]["min_purity"] = json!("x"),
        },
        Case {
            label: "cinco suelos, en los dos",
            expected: "criteria.min_copies no tiene seis suelos",
            table: same,
            metadata: |m| m["criteria"]["min_copies"] = json!([8, 6, 6, 6, 5]),
            effective: |e| e["repeat_detection"]["min_copies"] = json!([8, 6, 6, 6, 5]),
        },
        Case {
            label: "un suelo que no es un entero, en los dos",
            expected: "criteria.min_copies no es un entero no negativo",
            table: same,
            metadata: |m| m["criteria"]["min_copies"][0] = json!("8"),
            effective: |e| e["repeat_detection"]["min_copies"][0] = json!("8"),
        },
        Case {
            label: "unas regiones que no son una lista, en los dos",
            expected: "criteria.regions no es una lista",
            table: same,
            metadata: |m| m["criteria"]["regions"] = json!("todas"),
            effective: |e| e["repeat_detection"]["regions"] = json!("todas"),
        },
        Case {
            label: "unos tramos sustituidos que no son una lista",
            expected: "replaced_tracts no es una lista",
            table: same,
            metadata: |m| m["replaced_tracts"] = json!({}),
            effective: untouched,
        },
        Case {
            label: "un tramo sustituido con una puntuación que no es un entero",
            expected: "replaced_tracts: la puntuación no es un entero",
            table: same,
            metadata: |m| m["replaced_tracts"][0]["score"] = json!("40"),
            effective: untouched,
        },
    ]);
}

/// **La lectura de la tabla.**
#[test]
fn the_table_is_read_strictly() {
    run_cases(&[
        Case {
            label: "sin salto de línea final",
            expected: "la tabla no termina en un salto de línea",
            table: |t| t.trim_end().to_owned(),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "otra cabecera",
            expected: "la cabecera de la tabla no es la del contrato",
            table: |t| swap(t, "locus_id\t", "id\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "una fila de seis columnas",
            expected: "una fila no tiene siete columnas: \"SSR-1-301-312-A\\tchr1\\t301\\t312\\tA\\t22\"",
            table: |t| swap(t, "\t22\tdetected\n", "\t22\n"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un inicio con un cero a la izquierda",
            expected: "SSR-1-301-312-A: el inicio o el fin no es un entero canónico",
            table: |t| swap(t, "\t301\t", "\t0301\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un origen desconocido",
            expected: "SSR-1-301-312-A: origen desconocido \"found\"",
            table: |t| swap(t, "\t22\tdetected", "\t22\tfound"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un detectado sin puntuación",
            expected: "SSR-1-301-312-A: la puntuación de un detectado es un entero y la de un declarado es «.»",
            table: |t| swap(t, "\t22\tdetected", "\t.\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un declarado con puntuación",
            expected: "STR-AAAG: la puntuación de un detectado es un entero y la de un declarado es «.»",
            table: |t| swap(t, "\t.\tdeclared", "\t40\tdeclared"),
            metadata: untouched,
            effective: untouched,
        },
    ]);
    // Una tabla que no es UTF-8, escrita a mano.
    let package = Package::copy();
    let metadata = package.json("repeat_catalog.metadata.json");
    let effective = package.json("effective-config.json");
    package.republish(TABLE, &metadata, &effective);
    let mut bytes = TABLE.as_bytes().to_vec();
    bytes.push(0xFF);
    fs::write(package.dir().join("repeat_catalog.tsv"), &bytes).expect("tabla");
    let mut metadata = metadata.clone();
    metadata["table"]["sha256"] = json!(sha(&bytes));
    let metadata_bytes = package.write_json("repeat_catalog.metadata.json", &metadata);
    let mut effective = effective.clone();
    effective["repeat_detection"]["catalog"]["metadata_sha256"] = json!(sha(&metadata_bytes));
    package.write_json("effective-config.json", &effective);
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        match artifact["id"].as_str() {
            Some("repeat-catalog-tsv") => {
                artifact["checksum"]["digest"] = json!(sha(&bytes));
                artifact["size_bytes"] = json!(bytes.len());
            }
            Some("repeat-catalog-metadata") => {
                artifact["checksum"]["digest"] = json!(sha(&metadata_bytes));
                artifact["size_bytes"] = json!(metadata_bytes.len());
            }
            _ => {}
        }
    }
    package.write_json("manifest.json", &manifest);
    assert_eq!(detail(package.judge()), "la tabla no es UTF-8");
}

/// **Cada criterio sobre las bases del FASTA final**, con una fila añadida o
/// cambiada a mano. Las filas añadidas son de tractos plantados que el
/// productor descartó: cada uno falla un solo criterio.
#[test]
fn each_tract_criterion_is_judged_against_the_published_bases() {
    run_cases(&[
        Case {
            label: "una secuencia que no está en el FASTA",
            expected: "la secuencia chr9 no está en el FASTA final publicado",
            table: |t| swap(t, "SSR-1-301-312-A\tchr1\t", "SSR-1-301-312-A\tchr9\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un fin fuera de la secuencia",
            expected: "chr2:401..3001 no cabe en su secuencia",
            table: |t| swap(t, "\t401\t425\t", "\t401\t3001\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un inicio cero",
            expected: "chr1:0..312 no cabe en su secuencia",
            table: |t| swap(t, "\t301\t312\t", "\t0\t312\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un detectado con dos motivos",
            expected: "SSR-1-301-312-A: un locus detectado tiene un solo motivo",
            table: |t| swap(t, "\t312\tA\t", "\t312\tA,C\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un declarado con un motivo que no es ACGT",
            expected: "STR-AAAG: un motivo declarado no es una cadena de A, C, G y T",
            table: |t| swap(t, "\tAAAG\t.", "\tAANG\t."),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un motivo con una N",
            expected: "SSR-1-301-312-A (chr1:301..312): el motivo no tiene de 1 a 6 bases A, C, G o T",
            table: |t| swap(t, "\t312\tA\t", "\t312\tN\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un motivo de siete bases",
            expected: "SSR-1-301-312-A (chr1:301..312): el motivo no tiene de 1 a 6 bases A, C, G o T",
            table: |t| swap(t, "\t312\tA\t", "\t312\tAAAAAAC\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un motivo que no es primitivo",
            expected: "SSR-1-301-312-A (chr1:301..312): el motivo no es primitivo",
            table: |t| swap(t, "\t312\tA\t", "\t312\tAA\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un motivo que no es el de las bases del inicio",
            expected: "SSR-1-301-312-A (chr1:301..312): el motivo no es la lectura de las bases del inicio del tramo",
            table: |t| swap(t, "\t312\tA\t", "\t312\tC\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el A×12 pegado al principio de la secuencia",
            expected: "SSR-1-1-12-A (chr1:1..12): no tiene una base A, C, G o T a cada lado",
            table: |t| {
                swap(
                    t,
                    "SSR-1-301",
                    "SSR-1-1-12-A\tchr1\t1\t12\tA\t22\tdetected\nSSR-1-301",
                )
            },
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un A×12 con una N delante",
            expected: "SSR-2-2601-2612-A (chr2:2601..2612): no tiene una base A, C, G o T a cada lado",
            table: |t| with_row(t, "SSR-2-2601-2612-A\tchr2\t2601\t2612\tA\t22\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un A×12 con una N detrás",
            expected: "SSR-2-2701-2712-A (chr2:2701..2712): no tiene una base A, C, G o T a cada lado",
            table: |t| with_row(t, "SSR-2-2701-2712-A\tchr2\t2701\t2712\tA\t22\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un A×12 que acaba en la última base de la secuencia",
            expected: "SSR-2-2989-3000-A (chr2:2989..3000): no tiene una base A, C, G o T a cada lado",
            table: |t| with_row(t, "SSR-2-2989-3000-A\tchr2\t2989\t3000\tA\t22\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el A×210, por encima del tope",
            expected: "SSR-2-1001-1210-A (chr2:1001..1210): rebasa el tope de 200 unidades",
            table: |t| with_row(t, "SSR-2-1001-1210-A\tchr2\t1001\t1210\tA\t418\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el (TG)×5, por debajo del suelo",
            expected: "SSR-1-2401-2410-TG (chr1:2401..2410): tiene 5 copias completas y el suelo de motivos de 2 pb es 6",
            table: |t| with_row(t, "SSR-1-2401-2410-TG\tchr1\t2401\t2410\tTG\t16\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el A×12 alargado cinco bases: pureza 13/17 (una A suelta del fondo cuenta)",
            expected: "SSR-1-301-312-A (chr1:301..317): su pureza 13/17 no alcanza 0.8",
            table: |t| swap(t, "\t301\t312\t", "\t301\t317\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el tracto de pureza real 0,79",
            expected: "SSR-2-2000-2071-TGA (chr2:2000..2071): su pureza 57/72 no alcanza 0.8",
            table: |t| with_row(t, "SSR-2-2000-2071-TGA\tchr2\t2000\t2071\tTGA\t0\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el (CACAG)×3 CACA, tirada de 2",
            expected: "SSR-2-1601-1619-CA (chr2:1601..1619): su tirada más larga es de 2 copias y el mínimo es 3",
            table: |t| with_row(t, "SSR-2-1601-1619-CA\tchr2\t1601\t1619\tCA\t0\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "el (GT)×8, fuera de las regiones",
            expected: "SSR-2-2401-2416-GT (chr2:2401..2416): no cabe entero en ninguna región de detección",
            table: |t| with_row(t, "SSR-2-2401-2416-GT\tchr2\t2401\t2416\tGT\t28\tdetected"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "una puntuación que no es la de las comparaciones",
            expected: "SSR-1-301-312-A (chr1:301..312): su puntuación publicada es 23 y la de sus comparaciones 22",
            table: |t| swap(t, "\tA\t22\t", "\tA\t23\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un identificador que no sale de las coordenadas",
            expected: "SSR-1-301-312-C: el identificador de un detectado es SSR-1-301-312-A",
            table: |t| swap(t, "SSR-1-301-312-A\t", "SSR-1-301-312-C\t"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un sustituido por debajo del suelo",
            expected: "un tramo sustituido (chr1:1801..1820): tiene 5 copias completas y el suelo de motivos de 4 pb es 6",
            table: same,
            metadata: |m| m["replaced_tracts"][0]["end"] = json!(1820),
            effective: untouched,
        },
        Case {
            label: "un retirado con otra puntuación",
            expected: "SSR-1-801-820-CA (chr1:801..820): su puntuación publicada es 0 y la de sus comparaciones 36",
            table: same,
            metadata: |m| m["retired_tracts"][0]["score"] = json!(0),
            effective: untouched,
        },
        Case {
            label: "un retirado con otro identificador",
            expected: "SSR-1-801-820-AC: el identificador de un detectado es SSR-1-801-820-CA",
            table: same,
            metadata: |m| m["retired_tracts"][0]["locus_id"] = json!("SSR-1-801-820-AC"),
            effective: untouched,
        },
    ]);
}

/// **El orden, los solapes, la coincidencia de los declarados y los
/// retirados.**
#[test]
fn order_overlaps_declared_matches_and_retirements_are_judged() {
    run_cases(&[
        Case {
            label: "dos filas fuera de orden",
            expected: "SSR-1-1301-1321-GAT y SSR-1-301-312-A no van por secuencia y por inicio",
            table: |t| {
                let mut lines: Vec<&str> = t.lines().collect();
                lines.swap(1, 2);
                format!("{}\n", lines.join("\n"))
            },
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un declarado que pisa a un detectado",
            expected: "SSR-1-1301-1321-GAT y STR-AAAG se solapan",
            table: |t| swap(t, "\t1801\t1824\tAAAG", "\t1321\t1824\tAAAG"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un sustituido que no toca a ningún declarado",
            expected: "el tramo sustituido chr1:1801..1824 no toca a ningún locus declarado",
            table: |t| swap(t, "\t1801\t1824\tAAAG", "\t2201\t2224\tAAAG"),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un declarado de un motivo que no es rotación del sustituido",
            expected: "el locus declarado STR-AAAG no coincide con ningún tramo sustituido de un motivo rotación del suyo",
            table: |t| swap(t, "\tAAAG\t.", "\tAAGG\t."),
            metadata: untouched,
            effective: untouched,
        },
        Case {
            label: "un retirado cuya variante cae fuera",
            expected: "SSR-1-801-820-CA: la variante VAR en 2000 no cae dentro del tramo",
            table: same,
            metadata: |m| {
                m["retired_tracts"][0]["variant"]["position"] = json!(2000);
                m["retired_tracts"][0]["variant"]["id"] = json!("VAR");
            },
            effective: untouched,
        },
        Case {
            label: "un retirado por una variante de la matriz que no existe",
            expected: "SSR-1-801-820-CA: no hay ninguna variante VAR-INVENTADA de origen population_snvs en chr1:810",
            table: same,
            metadata: |m| m["retired_tracts"][0]["variant"]["id"] = json!("VAR-INVENTADA"),
            effective: untouched,
        },
        Case {
            label: "un retirado por un cambio directo que no existe",
            expected: "SSR-1-801-820-CA: no hay ninguna variante snp-x de origen direct_event en chr1:810",
            table: same,
            metadata: |m| {
                m["retired_tracts"][0]["variant"]["id"] = json!("snp-x");
                m["retired_tracts"][0]["variant"]["source"] = json!("direct_event");
            },
            effective: untouched,
        },
        Case {
            label: "un detectado que pisa a un retirado",
            expected: "dos tramos aceptados (detectados, sustituidos o retirados) se solapan",
            table: |t| {
                swap(
                    t,
                    "SSR-1-1301",
                    "SSR-1-801-820-CA\tchr1\t801\t820\tCA\t36\tdetected\nSSR-1-1301",
                )
            },
            metadata: untouched,
            effective: untouched,
        },
    ]);
}

/// Publica un paquete, le aplica `alter` al metadato y a la efectiva, rehace
/// las huellas y devuelve el veredicto del juez del catálogo.
fn judged_after(alter: impl FnOnce(&Package, &mut Value, &mut Value)) -> String {
    let package = Package::copy();
    let mut metadata = package.json("repeat_catalog.metadata.json");
    let mut effective = package.json("effective-config.json");
    let before = (metadata.clone(), effective.clone());
    alter(&package, &mut metadata, &mut effective);
    assert!(
        before != (metadata.clone(), effective.clone()),
        "la alteración no alteró"
    );
    package.republish(TABLE, &metadata, &effective);
    detail(package.judge())
}

/// Un cambio de base escrito a mano en la efectiva.
fn direct_event(kind: &str, contig: &str, position: u64) -> Value {
    json!([{
        "event_id": "snp-x", "individual_id": "I1", "phase": 0, "kind": kind,
        "contig": contig, "position_1based": position, "ref": "C", "alt": "T"
    }])
}

/// El primer retirado (`SSR-1-801-820-CA`) pasa a citar el cambio directo
/// `snp-x` en `position`.
fn retire_by_direct_event(metadata: &mut Value, position: u64) {
    metadata["retired_tracts"][0]["variant"] =
        json!({"source": "direct_event", "id": "snp-x", "position": position});
}

/// **La variante de un retirado tiene que existir en su secuencia, en su
/// posición y con su identificador**, y caer dentro del tramo (primera y última
/// base incluidas). Cada caso cambia una sola de esas cosas.
#[test]
fn the_retiring_variant_is_pinned_to_its_sequence_position_and_identifier() {
    let missing =
        "SSR-1-801-820-CA: no hay ninguna variante snp-x de origen direct_event en chr1:810";
    // Lado aceptado: el cambio directo existe donde dice el retirado, también
    // en la primera y en la última base del tramo.
    for position in [810, 801, 820] {
        assert_eq!(
            judged_after(|_, metadata, effective| {
                retire_by_direct_event(metadata, position);
                effective["direct_events"] = direct_event("snp", "chr1", position);
            }),
            "aceptado",
            "en {position}"
        );
    }
    // Una base antes del tramo.
    assert_eq!(
        judged_after(|_, metadata, effective| {
            retire_by_direct_event(metadata, 800);
            effective["direct_events"] = direct_event("snp", "chr1", 800);
        }),
        "SSR-1-801-820-CA: la variante snp-x en 800 no cae dentro del tramo"
    );
    // El cambio existe, pero en otra posición del tramo.
    assert_eq!(
        judged_after(|_, metadata, effective| {
            retire_by_direct_event(metadata, 810);
            effective["direct_events"] = direct_event("snp", "chr1", 811);
        }),
        missing
    );
    // El cambio existe en la misma posición de otra secuencia.
    assert_eq!(
        judged_after(|_, metadata, effective| {
            retire_by_direct_event(metadata, 810);
            effective["direct_events"] = direct_event("snp", "chr2", 810);
        }),
        missing
    );
    // Existe un evento con ese identificador y esa posición, pero no es un
    // cambio de base.
    assert_eq!(
        judged_after(|_, metadata, effective| {
            retire_by_direct_event(metadata, 810);
            effective["direct_events"] = direct_event("insertion", "chr1", 810);
        }),
        missing
    );
    // La matriz: el identificador real de otra fila (la de chr1:2000) puesto
    // en la posición del retirado.
    let package = Package::copy();
    let other = package
        .text("population_snvs.tsv")
        .lines()
        .skip(1)
        .find(|line| line.split('\t').nth(2) == Some("2000"))
        .and_then(|line| line.split('\t').next().map(str::to_owned))
        .expect("la fila de chr1:2000");
    drop(package);
    let expected = format!(
        "SSR-1-801-820-CA: no hay ninguna variante {other} de origen population_snvs en chr1:810"
    );
    assert_eq!(
        judged_after(|_, metadata, _| {
            metadata["retired_tracts"][0]["variant"]["id"] = json!(other.clone());
        }),
        expected
    );
}

/// **La cabecera del metadato**: cada campo, alterado solo, se rechaza.
#[test]
fn each_header_field_of_the_metadata_is_judged() {
    let header = "la clase, la versión, el sistema de coordenadas, la versión del algoritmo o la tabla (ruta y huella) del metadato no son los esperados";
    for (field, value) in [
        ("schema_version", json!("2")),
        ("coordinate_system", json!("base0_half_open")),
        ("algorithm_version", json!("ssr-detection-v0")),
    ] {
        assert_eq!(
            judged_after(|_, metadata, _| metadata[field] = value.clone()),
            header,
            "{field}"
        );
    }
    // La huella de la tabla dentro del metadato: se escribe a mano después de
    // rehacer las demás huellas, que siguen cuadrando.
    let package = Package::copy();
    let mut metadata = package.json("repeat_catalog.metadata.json");
    metadata["table"]["sha256"] = json!("0".repeat(64));
    let bytes = package.write_json("repeat_catalog.metadata.json", &metadata);
    let mut effective = package.json("effective-config.json");
    effective["repeat_detection"]["catalog"]["metadata_sha256"] = json!(sha(&bytes));
    package.write_json("effective-config.json", &effective);
    let mut manifest = package.json("manifest.json");
    for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
        if artifact["id"] == "repeat-catalog-metadata" {
            artifact["checksum"]["digest"] = json!(sha(&bytes));
            artifact["size_bytes"] = json!(bytes.len());
        }
    }
    package.write_json("manifest.json", &manifest);
    assert_eq!(detail(package.judge()), header);
}

/// **Dos aceptados que comparten una sola base se solapan**: el retirado de
/// las once `A` alargado una base hacia atrás pisa la última del `(CA)×8`
/// detectado. El lado aceptado, los dos tocándose sin solaparse, es el paquete
/// intacto (lo acepta `the_published_catalog_is_the_planted_one_and_its_counts_add_up`).
#[test]
fn two_accepted_tracts_sharing_one_base_overlap() {
    assert_eq!(
        judged_after(|_, metadata, _| {
            let retired = &mut metadata["retired_tracts"][1];
            assert_eq!(retired["start"], 2617);
            retired["start"] = json!(2616);
            retired["score"] = json!(22);
            retired["locus_id"] = json!("SSR-1-2616-2627-A");
        }),
        "dos tramos aceptados (detectados, sustituidos o retirados) se solapan"
    );
}

/// **El tamaño de cada fichero del catálogo es el del manifiesto**, como en
/// la matriz de SNV.
#[test]
fn the_sizes_in_the_manifest_are_judged() {
    for (id, expected) in [
        (
            "repeat-catalog-tsv",
            "el tamaño de la tabla no es el del manifiesto",
        ),
        (
            "repeat-catalog-metadata",
            "el tamaño del metadato no es el del manifiesto",
        ),
    ] {
        let package = Package::copy();
        let mut manifest = package.json("manifest.json");
        for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
            if artifact["id"] == id {
                artifact["size_bytes"] =
                    json!(artifact["size_bytes"].as_u64().expect("tamaño") + 1);
            }
        }
        package.write_json("manifest.json", &manifest);
        assert_eq!(detail(package.judge()), expected, "{id}");
        assert_eq!(
            detail(package.judge_resources()),
            expected,
            "{id}, juez completo"
        );
    }
}

/// **Los recuentos.** Cada caso mantiene cuadradas las demás sumas, para que
/// solo falle la que se prueba.
#[test]
fn the_counts_are_judged() {
    fn entry(metadata: &mut Value, length: usize) -> &mut Value {
        &mut metadata["counts_by_motif_length"][length - 1]
    }
    fn bump(value: &mut Value, by: i64) {
        *value = json!(value.as_i64().expect("entero") + by);
    }
    run_cases(&[
        Case {
            label: "cinco entradas",
            expected: "counts_by_motif_length no tiene una entrada por longitud de motivo de 1 a 6",
            table: same,
            metadata: |m| {
                m["counts_by_motif_length"]
                    .as_array_mut()
                    .expect("lista")
                    .pop();
            },
            effective: untouched,
        },
        Case {
            label: "dos entradas cambiadas de sitio",
            expected: "la entrada 1 de counts_by_motif_length no es la de motivos de 1 pb",
            table: same,
            metadata: |m| {
                m["counts_by_motif_length"]
                    .as_array_mut()
                    .expect("lista")
                    .swap(0, 1)
            },
            effective: untouched,
        },
        Case {
            label: "una longitud que no es un entero",
            expected: "motif_length no es un entero no negativo",
            table: same,
            metadata: |m| entry(m, 1)["motif_length"] = json!("1"),
            effective: untouched,
        },
        Case {
            label: "una etapa de menos",
            expected: "motivos de 1 pb: los descartes no son los de las ocho etapas",
            table: same,
            metadata: |m| {
                entry(m, 1)["discarded"]
                    .as_object_mut()
                    .expect("etapas")
                    .remove("edge");
            },
            effective: untouched,
        },
        Case {
            label: "una etapa con otro nombre",
            expected: "motivos de 1 pb: los descartes no son los de las ocho etapas",
            table: same,
            metadata: |m| {
                let stages = entry(m, 1)["discarded"].as_object_mut().expect("etapas");
                let value = stages.remove("edge").expect("borde");
                stages.insert("border".to_owned(), value);
            },
            effective: untouched,
        },
        Case {
            label: "un intervalo de más en el detector",
            expected: "motivos de 2 pb: los intervalos del detector no son los descartes más los aceptados",
            table: same,
            metadata: |m| bump(&mut entry(m, 2)["detector_intervals"], 1),
            effective: untouched,
        },
        Case {
            label: "un aceptado de más, también en el detector",
            expected: "motivos de 2 pb: los aceptados no son los sustituidos, los retirados, los que quedan sin copias y los publicados",
            table: same,
            metadata: |m| {
                bump(&mut entry(m, 2)["accepted"], 1);
                bump(&mut entry(m, 2)["detector_intervals"], 1);
            },
            effective: untouched,
        },
        Case {
            label: "un sustituido de más y un publicado de menos",
            expected: "motivos de 1 pb: replaced_by_declared no es lo que lista replaced_tracts",
            table: same,
            metadata: |m| {
                bump(&mut entry(m, 1)["replaced_by_declared"], 1);
                bump(&mut entry(m, 1)["published"], -1);
            },
            effective: untouched,
        },
        Case {
            label: "un retirado de más y un publicado de menos",
            expected: "motivos de 1 pb: retired_by_variant no es lo que lista retired_tracts",
            table: same,
            metadata: |m| {
                bump(&mut entry(m, 1)["retired_by_variant"], 1);
                bump(&mut entry(m, 1)["published"], -1);
            },
            effective: untouched,
        },
        Case {
            label: "un publicado de más, con su aceptado y su intervalo",
            expected: "motivos de 1 pb: published dice 2 y la tabla tiene 1 loci detectados",
            table: same,
            metadata: |m| {
                bump(&mut entry(m, 1)["published"], 1);
                bump(&mut entry(m, 1)["accepted"], 1);
                bump(&mut entry(m, 1)["detector_intervals"], 1);
            },
            effective: untouched,
        },
        Case {
            label: "sin copias haplotípicas y con publicados",
            expected: "motivos de 1 pb: sin copias haplotípicas no entra ningún detectado, y hay publicados o retirados",
            table: same,
            metadata: |m| {
                bump(&mut entry(m, 1)["without_haplotype_copies"], 1);
                bump(&mut entry(m, 1)["accepted"], 1);
                bump(&mut entry(m, 1)["detector_intervals"], 1);
            },
            effective: untouched,
        },
    ]);
}

/// **La tabla dice lo mismo que `repeat_loci` de la efectiva.**
#[test]
fn the_table_says_what_the_effective_loci_say() {
    run_cases(&[
        Case {
            label: "un locus de menos en la efectiva",
            expected: "la tabla tiene 5 loci y repeat_loci de la efectiva 4",
            table: same,
            metadata: untouched,
            effective: |e| {
                e["repeat_loci"].as_array_mut().expect("loci").pop();
            },
        },
        Case {
            label: "otro origen en la efectiva",
            expected: "SSR-1-301-312-A: la fila no es el locus de repeat_loci de la efectiva",
            table: same,
            metadata: untouched,
            effective: |e| e["repeat_loci"][0]["origin"] = json!("declared"),
        },
        Case {
            label: "otro fin en la efectiva",
            expected: "SSR-2-401-425-ACGTT: la fila no es el locus de repeat_loci de la efectiva",
            table: same,
            metadata: untouched,
            effective: |e| e["repeat_loci"][4]["end"] = json!(424),
        },
    ]);
}

// ---------------------------------------------------------------------------
// El juez completo
// ---------------------------------------------------------------------------

/// **El juez del catálogo está encadenado en el juez completo**: una
/// puntuación alterada en la tabla, con las huellas rehechas, la rechaza
/// `validate_published_resources` con el motivo del juez del catálogo.
#[test]
fn the_full_judge_runs_the_catalog_judge() {
    let package = Package::copy();
    let metadata = package.json("repeat_catalog.metadata.json");
    let effective = package.json("effective-config.json");
    package.republish(
        &swap(TABLE, "\tA\t22\t", "\tA\t23\t"),
        &metadata,
        &effective,
    );
    assert_eq!(
        detail(package.judge_resources()),
        "SSR-1-301-312-A (chr1:301..312): su puntuación publicada es 23 y la de sus comparaciones 22"
    );
}

/// **Lo que el juez independiente no puede ver lo ata la reconstrucción.** Un
/// descarte que cambia de etapa (las sumas siguen cuadrando) y un detectado
/// que desaparece de la tabla y de la efectiva, contado como descarte por
/// solape: los dos los acepta el juez del catálogo —es su límite: no rehace el
/// barrido— y los rechaza la reconstrucción del juez completo, que llama al
/// productor.
#[test]
fn the_reconstruction_ties_what_the_judge_cannot_see() {
    let reconstruction = "el catálogo de microsatélites publicado no es el que el TOML, la referencia y la matriz producen";
    let package = Package::copy();
    let mut metadata = package.json("repeat_catalog.metadata.json");
    let effective = package.json("effective-config.json");
    let stages = &mut metadata["counts_by_motif_length"][1]["discarded"];
    let moved = stages["min_copies"].as_u64().expect("entero");
    assert!(moved > 0, "hace falta un descarte por suelo de periodo 2");
    stages["min_copies"] = json!(moved - 1);
    stages["min_purity"] = json!(stages["min_purity"].as_u64().expect("entero") + 1);
    package.republish(TABLE, &metadata, &effective);
    assert_eq!(detail(package.judge()), "aceptado", "el límite del juez");
    assert_eq!(detail(package.judge_resources()), reconstruction);

    let package = Package::copy();
    let mut metadata = package.json("repeat_catalog.metadata.json");
    let mut effective = package.json("effective-config.json");
    let counts = &mut metadata["counts_by_motif_length"][4];
    counts["published"] = json!(0);
    counts["accepted"] = json!(counts["accepted"].as_u64().expect("entero") - 1);
    counts["discarded"]["overlap"] =
        json!(counts["discarded"]["overlap"].as_u64().expect("entero") + 1);
    effective["repeat_loci"].as_array_mut().expect("loci").pop();
    let table = TABLE.replace(
        "SSR-2-401-425-ACGTT\tchr2\t401\t425\tACGTT\t40\tdetected\n",
        "",
    );
    package.republish(&table, &metadata, &effective);
    assert_eq!(
        detail(package.judge()),
        "aceptado",
        "el límite del juez: la exhaustividad"
    );
    // Lo caza la comparación de los loci de la efectiva con los que rehace el
    // TOML, que también detecta.
    assert_eq!(
        detail(package.judge_resources()),
        "los loci de repetición efectivos no son los que el TOML y la referencia producen"
    );
    // Y si solo falta en la tabla (un escritor de la tabla que perdiera una
    // fila, que la reconstrucción no vería porque usa el mismo escritor), lo
    // caza el juez del catálogo al compararla con la efectiva.
    let package = Package::copy();
    package.republish(&table, &metadata, &package.json("effective-config.json"));
    assert_eq!(
        detail(package.judge()),
        "la tabla tiene 4 loci y repeat_loci de la efectiva 5"
    );
}

/// **El juez no llama al productor.** Guarda léxica sobre el texto del juez:
/// no nombra el detector, la medida de pureza del productor, la conversión de
/// tramos ni el módulo de anotación de tractos. No ve una llamada a través de
/// otro módulo que los use; la prueba de independencia del plan (con la
/// pureza del productor estropeada en una copia del árbol) es la que lo
/// demuestra de verdad.
#[test]
fn the_catalog_judge_names_no_producer() {
    let source = include_str!("../src/repeat_catalog_validation.rs");
    let code: String = source
        .lines()
        .filter(|line| !line.trim_start().starts_with("//"))
        .collect::<Vec<_>>()
        .join("\n");
    for producer in [
        "repeat_detection::",
        "observe_tract",
        "classify_interval",
        "tandem_repeat::",
        "with_detected_repeat_loci",
        "results::",
        "application::",
    ] {
        assert!(
            !code.contains(producer),
            "el juez del catálogo nombra {producer}"
        );
    }
    // Control positivo del instrumento: ve lo que hay.
    assert!(code.contains("fn decompose") && code.contains("guard_trace::record"));
}
