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
    publish_with(label, detection, "")
}

/// Lo mismo que [`publish`], con `extra` añadido al final del TOML (el VCF,
/// por ejemplo).
fn publish_with(label: &str, detection: bool, extra: &str) -> PathBuf {
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
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n{extra}"
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
        Self::copy_with("")
    }

    /// Un paquete con `extra` añadido al TOML (ver [`publish_with`]).
    fn copy_with(extra: &str) -> Self {
        let destination = publish_with("case", true, extra);
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
        self.republish_with_alleles(table, self.text("repeat_alleles.tsv"), metadata, effective);
    }

    /// Lo mismo que [`Package::republish`] con la tabla de alelos dada, cuya
    /// huella rehace en el enlace de la efectiva y en el manifiesto.
    fn republish_with_alleles(
        &self,
        table: &str,
        alleles: impl AsRef<[u8]>,
        metadata: &Value,
        effective: &Value,
    ) {
        let alleles = alleles.as_ref();
        fs::write(self.dir().join("repeat_catalog.tsv"), table).expect("tabla");
        fs::write(self.dir().join("repeat_alleles.tsv"), alleles).expect("tabla de alelos");
        let mut metadata = metadata.clone();
        if metadata.is_object() {
            metadata["table"]["sha256"] = json!(sha(table.as_bytes()));
        }
        let metadata_bytes = self.write_json("repeat_catalog.metadata.json", &metadata);
        let mut effective = effective.clone();
        if effective.get("repeat_detection").is_some() {
            effective["repeat_detection"]["catalog"]["metadata_sha256"] =
                json!(sha(&metadata_bytes));
            if effective["repeat_detection"]["alleles_table"].is_object() {
                effective["repeat_detection"]["alleles_table"]["sha256"] = json!(sha(alleles));
            }
        }
        let effective_bytes = self.write_json("effective-config.json", &effective);
        let mut manifest = self.json("manifest.json");
        for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
            let bytes = match artifact["path"].as_str() {
                Some("repeat_catalog.tsv") => table.as_bytes().to_vec(),
                Some("repeat_alleles.tsv") => alleles.to_vec(),
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
    let effective = package.json("effective-config.json");
    let counts = &mut metadata["counts_by_motif_length"][4];
    counts["published"] = json!(0);
    counts["accepted"] = json!(counts["accepted"].as_u64().expect("entero") - 1);
    counts["discarded"]["overlap"] =
        json!(counts["discarded"]["overlap"].as_u64().expect("entero") + 1);
    let table = TABLE.replace(
        "SSR-2-401-425-ACGTT\tchr2\t401\t425\tACGTT\t40\tdetected\n",
        "",
    );
    let alleles = alter_rows(&package.text("repeat_alleles.tsv"), |rows| {
        rows.retain(|row| row[0] != "SSR-2-401-425-ACGTT");
    });
    package.republish_with_alleles(&table, &alleles, &metadata, &effective);
    assert_eq!(
        detail(package.judge()),
        "aceptado",
        "el límite del juez: la exhaustividad"
    );
    // Desde la entrega 4 bis la efectiva ya no lleva los detectados: lo caza
    // la reconstrucción, que también detecta.
    assert_eq!(detail(package.judge_resources()), reconstruction);
    // Y si solo falta en una de las dos tablas, lo caza el juez del catálogo
    // al compararlas. Si faltara en las dos a la vez (un escritor que perdiera
    // la misma fila en las dos), la reconstrucción no lo vería porque usa los
    // mismos escritores: lo cazan las listas escritas a mano de
    // `the_detected_loci_leave_the_json_documents` y de
    // `the_published_catalog_is_the_planted_one_and_its_counts_add_up`.
    let package = Package::copy();
    package.republish(&table, &metadata, &package.json("effective-config.json"));
    assert_eq!(
        detail(package.judge()),
        "los loci de la tabla de alelos no son los detectados del catálogo, en su orden"
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
        "repeat_slippage::",
        "repeat_genotypes::",
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

// ---------------------------------------------------------------------------
// Entrega 4 bis (F6.4.a.15c): los detectados salen de los documentos JSON y
// sus alelos se publican en `repeat_alleles.tsv`
// ---------------------------------------------------------------------------

/// Los detectados del paquete, escritos a mano en el orden del catálogo: son
/// las filas `detected` de [`TABLE`]. Un escritor que perdiera la misma fila
/// en el catálogo y en la tabla de alelos pasaría el juez (las dos tablas
/// dirían lo mismo) y la reconstrucción (usa los mismos escritores); lo caza
/// esta lista.
const DETECTED: [&str; 4] = [
    "SSR-1-301-312-A",
    "SSR-1-1301-1321-GAT",
    "SSR-1-2601-2616-CA",
    "SSR-2-401-425-ACGTT",
];

/// Cabecera de la tabla de alelos, escrita a mano (DTC-104).
const ALLELES_HEADER: &str = "locus_id\tallele_index\tsequence\torigin\tdeclared_frequency\t\
effective_frequency\trealized_copies\trealized_frequency\tdrawn_events\tapplied_events\t\
unapplied_events\ttruth_tract_bases_bound\tvcf_tract_bases_bound";

/// El VCF de texto, para los casos que atribuyen bases de tracto a los loci.
const VCF: &str = "\n[vcf_output]\nmode = \"text\"\n";

/// Las filas de la tabla de alelos, cada una como sus trece campos.
fn allele_rows(alleles: &str) -> Vec<Vec<String>> {
    alleles
        .lines()
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

/// Vuelve a escribir la tabla de alelos desde sus filas.
fn alleles_text(rows: &[Vec<String>]) -> String {
    let mut text = format!("{ALLELES_HEADER}\n");
    for row in rows {
        text.push_str(&row.join("\t"));
        text.push('\n');
    }
    text
}

/// Aplica `alter` a las filas de la tabla de alelos.
fn alter_rows(alleles: &str, alter: impl FnOnce(&mut Vec<Vec<String>>)) -> String {
    let mut rows = allele_rows(alleles);
    alter(&mut rows);
    alleles_text(&rows)
}

/// El índice de la primera fila de un locus con al menos dos alelos y ese
/// locus: varios casos necesitan dos alelos del mismo locus.
fn locus_with_two_alleles(rows: &[Vec<String>]) -> (usize, String) {
    rows.windows(2)
        .position(|pair| pair[0][0] == pair[1][0])
        .map(|index| (index, rows[index][0].clone()))
        .expect("algún detectado tiene dos alelos")
}

/// **Los detectados salen de los dos documentos JSON** (REQ-MAN-95,
/// REQ-MAN-98, DTC-104): `repeat_loci` y `repeat_output_budget.per_locus` de
/// la efectiva y `repeat_loci.loci` del manifiesto llevan solo el declarado
/// `STR-AAAG`, y ni la efectiva ni el manifiesto nombran un solo `SSR-`. Lo
/// común del deslizamiento va una vez en `repeat_detection`, con la tabla de
/// alelos enlazada por su huella; la tabla lleva los detectados escritos a
/// mano ([`DETECTED`]), en su orden, con su cabecera escrita a mano.
///
/// **Debe fallar contra** un abanico de un detectado que se quede en la
/// efectiva o en el manifiesto (los dos `SSR-` de abajo), y contra un
/// escritor que pierda la misma fila en el catálogo y en la tabla de alelos.
#[test]
fn the_detected_loci_leave_the_json_documents() {
    let package = Package::copy();
    let effective = package.json("effective-config.json");
    let manifest = package.json("manifest.json");
    let ids = |loci: &Value| -> Vec<String> {
        loci.as_array()
            .expect("lista")
            .iter()
            .map(|locus| locus["id"].as_str().expect("id").to_owned())
            .collect()
    };
    assert_eq!(ids(&effective["repeat_loci"]), ["STR-AAAG"]);
    assert_eq!(
        ids(&effective["repeat_output_budget"]["per_locus"]),
        ["STR-AAAG"]
    );
    assert_eq!(ids(&manifest["repeat_loci"]["loci"]), ["STR-AAAG"]);
    for file in ["effective-config.json", "manifest.json"] {
        assert!(
            !package.text(file).contains("SSR-"),
            "{file} nombra un locus detectado"
        );
    }
    let detection = &effective["repeat_detection"];
    let slippage = &detection["slippage"];
    assert_eq!(slippage["model"], "stepwise");
    assert_eq!(slippage["step_units"], 1);
    assert_eq!(slippage["mean_events"], 1.0);
    assert_eq!(slippage["algorithm"], "biological-slippage-v1");
    assert_eq!(slippage["rng_scope"], "initial_repeat_slippage");
    assert_eq!(slippage["limitations"].as_array().expect("lista").len(), 3);
    // Dos individuos diploides: cuatro copias haplotípicas, cuatro linajes.
    assert_eq!(slippage["lineages"], 4);
    for field in ["drawn_events", "applied_events", "unapplied_events"] {
        assert!(slippage.get(field).is_none(), "{field} es de cada locus");
    }
    assert_eq!(
        detection["frequency_normalization"],
        json!({"tolerance": effective["normalization"]["tolerance"],
               "sum_algorithm": "neumaier_sum_v1", "canonical_order": "locus_id,allele_index"})
    );
    let alleles = package.text("repeat_alleles.tsv");
    assert_eq!(
        detection["alleles_table"],
        json!({"path": "repeat_alleles.tsv", "sha256": sha(alleles.as_bytes())})
    );
    assert_eq!(alleles.lines().next(), Some(ALLELES_HEADER));
    let rows = allele_rows(&alleles);
    let mut loci: Vec<&str> = rows.iter().map(|row| row[0].as_str()).collect();
    loci.dedup();
    assert_eq!(loci, DETECTED);
    for locus in DETECTED {
        let rows: Vec<&Vec<String>> = rows.iter().filter(|row| row[0] == locus).collect();
        let copies: u64 = rows
            .iter()
            .map(|row| row[6].parse::<u64>().expect("copias"))
            .sum();
        assert_eq!(
            copies, 4,
            "{locus}: las copias realizadas suman los linajes"
        );
        assert!(rows.iter().all(|row| row[3] == "generated"), "{locus}");
    }
    let manifest_artifact = manifest["artifacts"]
        .as_array()
        .expect("artefactos")
        .iter()
        .find(|artifact| artifact["id"] == "repeat-alleles-tsv")
        .expect("la tabla figura en el manifiesto");
    assert_eq!(manifest_artifact["kind"], "repeat_alleles_tsv");
    assert_eq!(manifest_artifact["path"], "repeat_alleles.tsv");
    assert_eq!(manifest_artifact["size_bytes"], alleles.len());
    assert!(
        manifest["output_categories"]["minimum"]
            .as_array()
            .expect("mínimo")
            .contains(&json!("repeat-alleles-tsv"))
    );
    package
        .judge_resources()
        .expect("control positivo: el juez completo acepta el paquete intacto");
}

/// **Sin la sección no se publica la tabla de alelos**, y el juez rechaza
/// una tabla de alelos listada en el manifiesto de una ejecución sin ella.
#[test]
fn without_the_section_no_alleles_table_is_published() {
    let destination = publish("without-alleles", false);
    assert!(!destination.join("repeat_alleles.tsv").exists());
    let package = Package {
        root: destination.parent().expect("raíz").to_path_buf(),
    };
    let mut manifest = package.json("manifest.json");
    manifest["artifacts"]
        .as_array_mut()
        .expect("artefactos")
        .push(json!({
            "checksum": {"algorithm": "sha256", "digest": "0".repeat(64)},
            "id": "repeat-alleles-tsv",
            "kind": "repeat_alleles_tsv",
            "path": "repeat_alleles.tsv",
            "size_bytes": 0
        }));
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge()),
        "hay un catálogo publicado y la configuración efectiva no declara la sección de microsatélites"
    );
}

/// **El enlace y la identidad de la tabla de alelos.** Aquí las huellas no se
/// rehacen después de la alteración que las rompe, salvo donde se dice.
#[test]
fn the_identity_of_the_alleles_table_is_judged() {
    type Alter = fn(&Package, &mut Value, &mut Value);
    let cases: [(&str, Alter); 7] = [
        (
            "repeat_detection.alleles_table no es un objeto",
            |_, effective, _| {
                effective["repeat_detection"]
                    .as_object_mut()
                    .expect("sección")
                    .remove("alleles_table");
            },
        ),
        (
            "repeat_detection.alleles_table.path no es la ruta de la tabla de alelos",
            |_, effective, _| {
                effective["repeat_detection"]["alleles_table"]["path"] = json!("otra.tsv");
            },
        ),
        ("alleles_table.sha256 no es un texto", |_, effective, _| {
            effective["repeat_detection"]["alleles_table"]["sha256"] = json!(7);
        }),
        (
            "el manifiesto no lista repeat-alleles-tsv con clase repeat_alleles_tsv y ruta \
             repeat_alleles.tsv",
            |_, _, manifest| {
                for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
                    if artifact["id"] == "repeat-alleles-tsv" {
                        artifact["kind"] = json!("otra");
                    }
                }
            },
        ),
        (
            "la huella de la tabla de alelos no es la del enlace de la efectiva o la del manifiesto",
            |_, effective, _| {
                effective["repeat_detection"]["alleles_table"]["sha256"] = json!("0".repeat(64));
            },
        ),
        (
            "la huella de la tabla de alelos no es la del enlace de la efectiva o la del manifiesto",
            |_, _, manifest| {
                for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
                    if artifact["id"] == "repeat-alleles-tsv" {
                        artifact["checksum"]["digest"] = json!("0".repeat(64));
                    }
                }
            },
        ),
        (
            "el tamaño de la tabla de alelos no es el del manifiesto",
            |_, _, manifest| {
                for artifact in manifest["artifacts"].as_array_mut().expect("artefactos") {
                    if artifact["id"] == "repeat-alleles-tsv" {
                        artifact["size_bytes"] =
                            json!(artifact["size_bytes"].as_u64().expect("tamaño") + 1);
                    }
                }
            },
        ),
    ];
    for (expected, alter) in cases {
        let package = Package::copy();
        let mut effective = package.json("effective-config.json");
        let mut manifest = package.json("manifest.json");
        let before = (effective.clone(), manifest.clone());
        alter(&package, &mut effective, &mut manifest);
        assert!(
            before != (effective.clone(), manifest.clone()),
            "{expected}: no alteró"
        );
        // El enlace del metadato del catálogo sigue apuntando a sus bytes; la
        // efectiva y el manifiesto se escriben sin rehacer nada más.
        package.write_json("effective-config.json", &effective);
        package.write_json("manifest.json", &manifest);
        assert_eq!(detail(package.judge()), expected);
    }
    // Sin el fichero.
    let package = Package::copy();
    fs::remove_file(package.dir().join("repeat_alleles.tsv")).expect("borrado");
    assert!(
        detail(package.judge()).starts_with("no se puede leer:"),
        "{}",
        detail(package.judge())
    );
}

/// Un caso de la tabla de alelos: qué se altera y el motivo exacto.
struct AllelesCase {
    label: &'static str,
    expected: &'static str,
    alleles: fn(&str) -> String,
    effective: fn(&mut Value),
}

/// Publica un paquete por caso, altera la tabla de alelos (y la efectiva),
/// rehace las huellas y exige el motivo exacto del juez del catálogo.
fn run_alleles_cases(cases: &[AllelesCase]) {
    for case in cases {
        let package = Package::copy();
        let alleles = package.text("repeat_alleles.tsv");
        let metadata = package.json("repeat_catalog.metadata.json");
        let mut effective = package.json("effective-config.json");
        let altered = (case.alleles)(&alleles);
        let before = effective.clone();
        (case.effective)(&mut effective);
        assert!(
            altered != alleles || before != effective,
            "«{}»: la alteración no alteró nada",
            case.label
        );
        package.republish_with_alleles(TABLE, &altered, &metadata, &effective);
        assert_eq!(
            detail(package.judge()),
            case.expected,
            "«{}» se rechazó por otra razón",
            case.label
        );
    }
}

fn same_alleles(alleles: &str) -> String {
    alleles.to_owned()
}

/// **La forma de la tabla de alelos**, leída con su lector estricto.
#[test]
fn the_alleles_table_is_read_strictly() {
    run_alleles_cases(&[
        AllelesCase {
            label: "sin salto final",
            expected: "la tabla de alelos no termina en un salto de línea",
            alleles: |a| a.trim_end_matches('\n').to_owned(),
            effective: untouched,
        },
        AllelesCase {
            label: "otra cabecera",
            expected: "la cabecera de la tabla de alelos no es la del contrato",
            alleles: |a| a.replacen("allele_index", "index", 1),
            effective: untouched,
        },
        AllelesCase {
            label: "una columna de menos",
            expected: "SSR-1-301-312-A: una fila de la tabla de alelos no tiene trece columnas",
            alleles: |a| {
                alter_rows(a, |rows| {
                    rows[0].pop();
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "entero no canónico",
            expected: "SSR-1-301-312-A: una columna entera no es un entero canónico",
            alleles: |a| alter_rows(a, |rows| rows[0][6] = format!("0{}", rows[0][6])),
            effective: untouched,
        },
        AllelesCase {
            label: "cota de la verdad no entera",
            expected: "SSR-1-301-312-A: una columna entera no es un entero canónico",
            alleles: |a| alter_rows(a, |rows| rows[0][11] = "-1".to_owned()),
            effective: untouched,
        },
        AllelesCase {
            label: "frecuencia declarada fuera de [0, 1]",
            expected: "SSR-1-301-312-A: una frecuencia no es un decimal en [0, 1]",
            alleles: |a| alter_rows(a, |rows| rows[0][4] = "1.5".to_owned()),
            effective: untouched,
        },
        AllelesCase {
            label: "frecuencia efectiva no decimal",
            expected: "SSR-1-301-312-A: una frecuencia no es un decimal en [0, 1]",
            alleles: |a| alter_rows(a, |rows| rows[0][5] = "x".to_owned()),
            effective: untouched,
        },
        AllelesCase {
            label: "frecuencia realizada vacía",
            expected: "SSR-1-301-312-A: la frecuencia realizada no es un decimal en [0, 1] ni «.»",
            alleles: |a| alter_rows(a, |rows| rows[0][7] = String::new()),
            effective: untouched,
        },
    ]);
}

/// **Una tabla de alelos que no es UTF-8**, con las huellas rehechas.
#[test]
fn an_alleles_table_that_is_not_utf8_is_rejected() {
    let package = Package::copy();
    let mut bytes = package.text("repeat_alleles.tsv").into_bytes();
    let at = bytes.len() - 2;
    bytes[at] = 0x80;
    package.republish_with_alleles(
        TABLE,
        &bytes,
        &package.json("repeat_catalog.metadata.json"),
        &package.json("effective-config.json"),
    );
    assert_eq!(detail(package.judge()), "la tabla de alelos no es UTF-8");
}

/// **Los loci y la aritmética de la tabla de alelos** (REQ-VAL-64). Cada caso
/// cambia una sola cosa; donde hace falta, recompone lo demás para que solo
/// falle lo que se prueba.
#[test]
fn the_alleles_table_arithmetic_is_judged() {
    let loci = "los loci de la tabla de alelos no son los detectados del catálogo, en su orden";
    run_alleles_cases(&[
        AllelesCase {
            label: "falta un detectado entero",
            expected: loci,
            alleles: |a| alter_rows(a, |rows| rows.retain(|row| row[0] != "SSR-1-2601-2616-CA")),
            effective: untouched,
        },
        AllelesCase {
            label: "un locus que no es un detectado",
            expected: loci,
            alleles: |a| {
                alter_rows(a, |rows| {
                    for row in rows
                        .iter_mut()
                        .filter(|row| row[0] == "SSR-2-401-425-ACGTT")
                    {
                        row[0] = "STR-AAAG".to_owned();
                    }
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "los detectados en otro orden",
            expected: loci,
            alleles: |a| {
                alter_rows(a, |rows| {
                    let first = rows[0][0].clone();
                    let (head, tail): (Vec<_>, Vec<_>) =
                        rows.drain(..).partition(|row| row[0] == first);
                    rows.extend(tail);
                    rows.extend(head);
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "un hueco en los índices",
            expected: "SSR-1-301-312-A: los índices de alelo no van de 0 en 0 sin huecos",
            alleles: |a| alter_rows(a, |rows| rows[0][1] = "1".to_owned()),
            effective: untouched,
        },
        AllelesCase {
            label: "un alelo declarado",
            expected: "SSR-1-301-312-A: un alelo de un detectado no es generated",
            alleles: |a| alter_rows(a, |rows| rows[0][3] = "declared".to_owned()),
            effective: untouched,
        },
        AllelesCase {
            label: "una X en la secuencia",
            expected: "SSR-1-301-312-A: la secuencia de un alelo no es una cadena de A, C, G, T y N",
            alleles: |a| alter_rows(a, |rows| rows[0][2].replace_range(0..1, "X")),
            effective: untouched,
        },
        AllelesCase {
            label: "una secuencia vacía",
            expected: "SSR-1-301-312-A: la secuencia de un alelo no es una cadena de A, C, G, T y N",
            alleles: |a| alter_rows(a, |rows| rows[0][2] = String::new()),
            effective: untouched,
        },
        AllelesCase {
            label: "pasos que no suman",
            expected: "SSR-1-301-312-A: los pasos aplicados y no aplicados no suman los sorteados",
            alleles: |a| {
                alter_rows(a, |rows| {
                    let drawn: u64 = rows[0][8].parse().expect("entero");
                    rows[0][8] = (drawn + 1).to_string();
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "una copia de más",
            expected: "SSR-1-301-312-A: las copias realizadas suman 5 y las copias haplotípicas \
                       publicadas son 4",
            alleles: |a| {
                alter_rows(a, |rows| {
                    let copies: u64 = rows[0][6].parse().expect("entero");
                    rows[0][6] = (copies + 1).to_string();
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "copias haplotípicas que la efectiva no publica",
            expected: "la efectiva no publica las copias haplotípicas de sus poblaciones iniciales",
            alleles: same_alleles,
            effective: |e| e["initial_populations"][0]["ploidy"] = json!("dos"),
        },
        AllelesCase {
            label: "una frecuencia realizada que no es copias entre copias",
            expected: "SSR-1-301-312-A: la frecuencia realizada de un alelo no es sus copias entre \
                       las del locus, o «.» sin copias",
            alleles: |a| {
                alter_rows(a, |rows| {
                    let frequency: f64 = rows[0][7].parse().expect("decimal");
                    rows[0][7] = if frequency == 0.0 {
                        "0.25".to_owned()
                    } else {
                        (frequency / 2.0).to_string()
                    };
                })
            },
            effective: untouched,
        },
        AllelesCase {
            label: "un punto con copias",
            expected: "SSR-1-301-312-A: la frecuencia realizada de un alelo no es sus copias entre \
                       las del locus, o «.» sin copias",
            alleles: |a| alter_rows(a, |rows| rows[0][7] = ".".to_owned()),
            effective: untouched,
        },
    ]);
}

/// **Dos alelos iguales en un locus** se rechazan nombrando el locus, que se
/// busca aquí: el primero con dos alelos.
#[test]
fn two_equal_alleles_are_named_by_their_locus() {
    let package = Package::copy();
    let alleles = package.text("repeat_alleles.tsv");
    let (_, locus) = locus_with_two_alleles(&allele_rows(&alleles));
    let altered = alter_rows(&alleles, |rows| {
        let (index, _) = locus_with_two_alleles(rows);
        rows[index + 1][2] = rows[index][2].clone();
    });
    package.republish_with_alleles(
        TABLE,
        &altered,
        &package.json("repeat_catalog.metadata.json"),
        &package.json("effective-config.json"),
    );
    assert_eq!(
        detail(package.judge()),
        format!("{locus}: dos alelos tienen la misma secuencia")
    );
}

/// **Las filas declaradas del catálogo son los loci de `repeat_loci` de la
/// efectiva** (REQ-VAL-64 desde la entrega 4 bis). El segundo caso es el
/// abanico de un detectado que se queda en la efectiva: el juez lo rechaza
/// por el recuento.
#[test]
fn the_declared_rows_are_the_effective_loci() {
    run_cases(&[
        Case {
            label: "el declarado falta en la efectiva",
            expected: "la tabla tiene 1 loci declarados y repeat_loci de la efectiva 0",
            table: same,
            metadata: untouched,
            effective: |e| {
                e["repeat_loci"].as_array_mut().expect("loci").clear();
            },
        },
        Case {
            label: "un detectado se queda en la efectiva",
            expected: "la tabla tiene 1 loci declarados y repeat_loci de la efectiva 2",
            table: same,
            metadata: untouched,
            effective: |e| {
                let mut detected = e["repeat_loci"][0].clone();
                detected["id"] = json!("SSR-1-301-312-A");
                detected["origin"] = json!("detected");
                e["repeat_loci"]
                    .as_array_mut()
                    .expect("loci")
                    .insert(0, detected);
            },
        },
        Case {
            label: "otro origen en la efectiva",
            expected: "STR-AAAG: la fila no es el locus de repeat_loci de la efectiva",
            table: same,
            metadata: untouched,
            effective: |e| e["repeat_loci"][0]["origin"] = json!("detected"),
        },
        Case {
            label: "otro fin en la efectiva",
            expected: "STR-AAAG: la fila no es el locus de repeat_loci de la efectiva",
            table: same,
            metadata: untouched,
            effective: |e| e["repeat_loci"][0]["end"] = json!(1823),
        },
    ]);
}

/// **Un detectado que se queda en la efectiva o en el manifiesto lo rechaza
/// también el juez completo**, al atar los dos documentos a lo que resuelve
/// el TOML: un locus `SSR-1-301-312-A` hecho a mano con el bloque del
/// declarado, añadido a la efectiva y, aparte, al manifiesto.
#[test]
fn a_detected_locus_left_in_the_json_documents_is_rejected_by_the_full_judge() {
    let package = Package::copy();
    let mut effective = package.json("effective-config.json");
    let mut detected = effective["repeat_loci"][0].clone();
    detected["id"] = json!("SSR-1-301-312-A");
    detected["origin"] = json!("detected");
    // Sin su bloque del generador, para que el juez del abanico de los
    // declarados, que lee la efectiva antes, no lo rehaga: así se ve la
    // atadura al TOML.
    detected
        .as_object_mut()
        .expect("locus")
        .remove("slippage")
        .expect("el declarado genera sus alelos");
    effective["repeat_loci"]
        .as_array_mut()
        .expect("loci")
        .insert(0, detected);
    package.republish(
        TABLE,
        &package.json("repeat_catalog.metadata.json"),
        &effective,
    );
    assert_eq!(
        detail(package.judge_resources()),
        "los loci de repetición efectivos no son los que el TOML y la referencia producen"
    );

    let package = Package::copy();
    let mut manifest = package.json("manifest.json");
    let mut detected = manifest["repeat_loci"]["loci"][0].clone();
    detected["id"] = json!("SSR-1-301-312-A");
    detected["origin"] = json!("detected");
    manifest["repeat_loci"]["loci"]
        .as_array_mut()
        .expect("loci")
        .insert(0, detected);
    package.write_json("manifest.json", &manifest);
    assert_eq!(
        detail(package.judge_resources()),
        "la sección de loci del manifiesto no es la que el TOML y la referencia producen"
    );
}

/// Publica un paquete (con `extra` en el TOML), altera la tabla de alelos y la
/// efectiva, rehace las huellas y devuelve el veredicto del **juez completo**.
fn judged_in_full(
    extra: &str,
    alter: impl FnOnce(&mut Vec<Vec<String>>, &mut Value),
) -> (String, Vec<Vec<String>>) {
    let package = Package::copy_with(extra);
    let alleles = package.text("repeat_alleles.tsv");
    let original = allele_rows(&alleles);
    let mut rows = original.clone();
    let mut effective = package.json("effective-config.json");
    let before = effective.clone();
    alter(&mut rows, &mut effective);
    assert!(
        rows != original || effective != before,
        "la alteración no alteró"
    );
    package.republish_with_alleles(
        TABLE,
        alleles_text(&rows),
        &package.json("repeat_catalog.metadata.json"),
        &effective,
    );
    (detail(package.judge_resources()), original)
}

/// **El abanico de cada detectado lo rehace un juez independiente desde lo
/// publicado** (REQ-VAL-60 para los detectados, entrega 4 bis): su fila del
/// catálogo, el deslizamiento común de la efectiva, la referencia y la
/// semilla publicadas. Cada caso altera la tabla de alelos o el deslizamiento
/// común y se rechaza en el juez completo con el motivo de ese juez, que corre
/// antes de atar nada al TOML.
///
/// El espectro emparejado al revés (las frecuencias de dos alelos
/// intercambiadas) lo mata, con el productor estropeado, la prueba de
/// independencia registrada en el plan; aquí cada comprobación se dispara con
/// la tabla alterada a mano.
#[test]
fn the_fan_of_each_detected_locus_is_rebuilt_from_what_is_published() {
    let prefix = |locus: &str, detail: &str| format!("abanico del detectado {locus}: {detail}");

    // La frecuencia declarada de un alelo, por otra admisible.
    let (verdict, original) = judged_in_full("", |rows, _| {
        let frequency: f64 = rows[0][4].parse().expect("decimal");
        rows[0][4] = if frequency == 0.25 { "0.5" } else { "0.25" }.to_owned();
    });
    let published = if original[0][4] == "0.25" {
        "0.5"
    } else {
        "0.25"
    };
    let rebuilt: f64 = original[0][4].parse().expect("decimal");
    assert_eq!(
        verdict,
        prefix(
            &original[0][0],
            &format!(
                "la frecuencia de un alelo no es la proporción de linajes que terminan en él \
                 (alelo 0: publicada {published}, rehecha {rebuilt} = {}/4)",
                (rebuilt * 4.0) as u64
            )
        )
    );

    // Dos alelos de un locus en otro orden: el orden canónico no es el suyo.
    let (verdict, original) = judged_in_full("", |rows, _| {
        let (index, _) = locus_with_two_alleles(rows);
        let (left, right) = (rows[index][2].clone(), rows[index + 1][2].clone());
        rows[index][2] = right;
        rows[index + 1][2] = left;
    });
    let (_, locus) = locus_with_two_alleles(&original);
    assert_eq!(
        verdict,
        prefix(
            &locus,
            "la secuencia de un alelo publicado no es la que el modelo produce, o el orden \
             canónico no es el suyo"
        )
    );

    // Un alelo de menos.
    let (verdict, original) = judged_in_full("", |rows, _| {
        let (index, _) = locus_with_two_alleles(rows);
        rows.remove(index + 1);
    });
    let (_, locus) = locus_with_two_alleles(&original);
    assert_eq!(
        verdict,
        prefix(
            &locus,
            "el número de alelos publicados no es el que el modelo produce"
        )
    );

    // La contabilidad de la **segunda** fila de un locus (la primera queda
    // bien): cada fila se compara, no solo la primera.
    let (verdict, original) = judged_in_full("", |rows, _| {
        let (index, _) = locus_with_two_alleles(rows);
        let row = &mut rows[index + 1];
        for column in [8, 9] {
            let value: u64 = row[column].parse().expect("entero");
            row[column] = (value + 1).to_string();
        }
    });
    let (_, locus) = locus_with_two_alleles(&original);
    assert_eq!(
        verdict,
        prefix(
            &locus,
            "la contabilidad publicada del generador en el alelo 1 no es la que el modelo produce"
        )
    );

    // El deslizamiento común.
    for (expected, alter) in [
        (
            "los linajes del deslizamiento común no son las copias haplotípicas de las \
             poblaciones iniciales",
            (|e: &mut Value| e["repeat_detection"]["slippage"]["lineages"] = json!(5))
                as fn(&mut Value),
        ),
        (
            "la media de pasos del deslizamiento común rebasa el tope que el generador admite",
            |e| e["repeat_detection"]["slippage"]["mean_events"] = json!(1000.5),
        ),
        (
            "el deslizamiento común publicado no tiene media y sesgo decimales",
            |e| e["repeat_detection"]["slippage"]["expansion_bias"] = json!("media"),
        ),
        ("el tope de unidades de la sección no es un entero", |e| {
            e["repeat_detection"]["max_units"] = json!("doscientas")
        }),
    ] {
        let (verdict, _) = judged_in_full("", |_, effective| alter(effective));
        assert_eq!(verdict, expected);
    }
}

/// **La tabla de alelos es la que resuelve el TOML copiado, con lo medido**
/// (REQ-VAL-62 para los detectados). Lo que el juez independiente no mira —las
/// copias de cada alelo frente a la verdad, la frecuencia efectiva, la cota de
/// la verdad del locus— lo ata la reconstrucción del juez completo.
#[test]
fn the_alleles_table_is_tied_to_the_toml_and_to_the_truth() {
    // Una copia pasa de un alelo a otro del mismo locus, con las frecuencias
    // realizadas rehechas: el juez del catálogo y el del abanico lo aceptan, y
    // la verdad lo desmiente.
    let (verdict, original) = judged_in_full("", |rows, _| {
        let (_, locus) = locus_with_two_alleles(rows);
        let of_locus: Vec<usize> = (0..rows.len()).filter(|i| rows[*i][0] == locus).collect();
        let from = *of_locus
            .iter()
            .find(|i| rows[**i][6] != "0")
            .expect("un alelo con copias");
        let to = *of_locus.iter().find(|i| **i != from).expect("otro alelo");
        for (row, delta) in [(from, -1_i64), (to, 1)] {
            let copies = rows[row][6].parse::<i64>().expect("entero") + delta;
            rows[row][6] = copies.to_string();
            rows[row][7] = (copies as f64 / 4.0).to_string();
        }
    });
    let (_, locus) = locus_with_two_alleles(&original);
    assert_eq!(
        verdict,
        format!("las copias realizadas de los alelos de {locus} no son las que la verdad lleva")
    );
    let tied = "la tabla de alelos publicada no es la que el TOML, la referencia, la matriz y lo \
                medido producen";
    // La frecuencia efectiva de un alelo.
    let (verdict, _) = judged_in_full("", |rows, _| {
        rows[0][5] = if rows[0][5] == "0.25" { "0.5" } else { "0.25" }.to_owned();
    });
    assert_eq!(verdict, tied);
    // La cota de la verdad de un locus, en todas sus filas.
    let (verdict, _) = judged_in_full("", |rows, _| {
        let locus = rows[0][0].clone();
        for row in rows.iter_mut().filter(|row| row[0] == locus) {
            let bound: u64 = row[11].parse().expect("entero");
            row[11] = (bound + 1).to_string();
        }
    });
    assert_eq!(verdict, tied);
}

/// **Las bases de tracto del VCF de cada detectado son las que llevan sus
/// registros**, atribuidos al primer locus **del catálogo** que los contiene
/// (REQ-VAL-62, DTC-097). Con VCF de texto.
///
/// **Debe fallar contra** un juez que atribuya las bases solo a los
/// declarados: el paquete intacto, cuyos registros de repetición caen en los
/// detectados, se rechazaría por un registro que no contiene ningún locus.
#[test]
fn the_vcf_tract_bases_of_each_detected_locus_are_recounted() {
    let package = Package::copy_with(VCF);
    let rows = allele_rows(&package.text("repeat_alleles.tsv"));
    let carrying: Vec<&str> = rows
        .iter()
        .filter(|row| row[12] != "0")
        .map(|row| row[0].as_str())
        .collect();
    assert!(
        !carrying.is_empty(),
        "control: algún detectado lleva bases de tracto en el VCF"
    );
    package
        .judge_resources()
        .expect("control positivo: el paquete con VCF intacto se acepta");

    // Las bases de un detectado pasan a otro: la suma no cambia.
    let (verdict, original) = judged_in_full(VCF, |rows, _| {
        let from = rows
            .iter()
            .find(|row| row[12] != "0")
            .map(|row| row[0].clone())
            .expect("un detectado con bases");
        let to = rows
            .iter()
            .find(|row| row[0] != from)
            .map(|row| row[0].clone())
            .expect("otro detectado");
        let bases: u64 = rows
            .iter()
            .find(|row| row[0] == from)
            .map(|row| row[12].parse().expect("entero"))
            .expect("fila");
        let other: u64 = rows
            .iter()
            .find(|row| row[0] == to)
            .map(|row| row[12].parse().expect("entero"))
            .expect("fila");
        for row in rows.iter_mut() {
            if row[0] == from {
                row[12] = "0".to_owned();
            } else if row[0] == to {
                row[12] = (other + bases).to_string();
            }
        }
    });
    let from = original
        .iter()
        .find(|row| row[12] != "0")
        .map(|row| row[0].clone())
        .expect("fila");
    let to = original
        .iter()
        .find(|row| row[0] != from)
        .map(|row| row[0].clone())
        .expect("fila");
    let first = DETECTED
        .iter()
        .find(|locus| **locus == from || **locus == to)
        .expect("uno de los dos");
    assert_eq!(
        verdict,
        format!(
            "las bases de tracto del VCF del detectado {first} no son las que llevan sus registros"
        )
    );

    // Solo la segunda fila de un locus: la primera, que es la que se suma,
    // queda bien.
    let (verdict, original) = judged_in_full(VCF, |rows, _| {
        let (index, _) = locus_with_two_alleles(rows);
        let bases: u64 = rows[index + 1][12].parse().expect("entero");
        rows[index + 1][12] = (bases + 1).to_string();
    });
    let (_, locus) = locus_with_two_alleles(&original);
    assert_eq!(
        verdict,
        format!(
            "las bases de tracto del VCF del detectado {locus} no son las que llevan sus registros"
        )
    );

    // Sin la tabla de alelos, el juez público de la contabilidad de tracto no
    // la inventa.
    fs::remove_file(package.dir().join("repeat_alleles.tsv")).expect("borrado");
    let dir = package.dir();
    let effective_path = dir.join("effective-config.json");
    let manifest_path = dir.join("manifest.json");
    let verdict = detail(
        genome_population_simulator::results::validate_published_tract_output_accounting(
            &effective_path,
            &fs::read(&effective_path).expect("efectiva"),
            &manifest_path,
            &fs::read(&manifest_path).expect("manifiesto"),
            &dir.join("outputs/truth.jsonl"),
            Some(&dir.join("outputs/truth.vcf")),
        ),
    );
    assert_eq!(
        verdict,
        "no se puede leer el catálogo de microsatélites o la tabla de alelos publicados"
    );
}

/// Publica una referencia de una secuencia con `count` microsatélites
/// plantados (motivos de 1 a 6 bases, uno cada 200 bases) y la sección de
/// microsatélites, y devuelve el paquete.
fn publish_planted(count: usize) -> Package {
    const MOTIFS: [(&str, usize); 6] = [
        ("A", 12),
        ("CA", 10),
        ("GAT", 8),
        ("AAAG", 7),
        ("ACGTT", 6),
        ("AACGTT", 5),
    ];
    // El guardián existe antes de publicar: un fallo no deja la carpeta.
    let package = Package {
        root: temporary(&format!("planted-{count}")),
    };
    let root = package.root.clone();
    let length = 200 * count + 200;
    let mut chr1 = background(0x9E37_79B9_7F4A_7C15, length);
    for index in 0..count {
        let (motif, units) = MOTIFS[index % MOTIFS.len()];
        plant(&mut chr1, 101 + 200 * index, motif, units);
    }
    write_fasta(&root.join("reference.fa"), &[("chr1", chr1)]);
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
             [population_snvs]\nmode = \"generate\"\n\n\
             [population_snvs.generate]\ncount_mode = \"exact\"\nlocus_count = 1\n\
             spatial_mode = \"uniform\"\nsubstitution_profile = \"uniform_snv12_v1\"\n\n\
             [[population_snvs.generate.populations]]\npopulation_id = \"P1\"\n\
             frequency_model = \"uniform\"\nmin_alt_af = 0.2\nmax_alt_af = 0.8\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = {length}\n\n\
             [reference.repeat_detection]\n\n\
             [reference.repeat_detection.slippage]\nmodel = \"stepwise\"\nmean_events = 1.0\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"I1\", \"I2\"]\n\n\
             [vcf_output]\nmode = \"text\"\n"
        ),
    )
    .expect("configuración");
    execute_validated_config(
        &load_and_validate_config(&config).expect("configuración válida"),
        &root.join("published"),
    )
    .unwrap_or_else(|error| panic!("publicación: {}", error.render_human()));
    package
}

/// **La efectiva y el manifiesto no crecen con los detectados** (REQ-MAN-98,
/// entrega 4 bis): con 12 y con 48 microsatélites plantados, los dos
/// documentos cambian unos pocos bytes (los dígitos de las longitudes, de los
/// totales y de los tamaños de los ficheros), mientras el catálogo y la tabla
/// de alelos crecen con cada locus.
///
/// **Debe fallar contra** dejar los detectados en los documentos JSON: hasta
/// la entrega 4 cada uno añadía unos 4 kB a cada documento (apartado 9 de la
/// propuesta de la entrega 4 bis en el plan).
#[test]
fn the_json_documents_do_not_grow_with_the_detected_loci() {
    let small = publish_planted(12);
    let large = publish_planted(48);
    let size = |package: &Package, file: &str| {
        fs::metadata(package.dir().join(file))
            .expect("fichero publicado")
            .len()
    };
    let rows = |package: &Package| package.text("repeat_catalog.tsv").lines().count() - 1;
    // Control: el detector halla lo plantado (alguno puede quedar fuera por la
    // fila de la matriz generada o por el fondo; no importa cuál).
    let (few, many) = (rows(&small), rows(&large));
    assert!(
        few >= 10 && many >= 44,
        "control: {few} y {many} detectados"
    );
    for file in ["effective-config.json", "manifest.json"] {
        let (a, b) = (size(&small, file), size(&large, file));
        assert!(
            b.abs_diff(a) <= 32,
            "{file} crece con los detectados: {a} bytes con 12 y {b} con 48"
        );
    }
    for file in ["repeat_catalog.tsv", "repeat_alleles.tsv"] {
        let (a, b) = (size(&small, file), size(&large, file));
        assert!(b > a + (many - few) as u64 * 40, "{file}: {a} y {b} bytes");
    }
    // Con solo detectados, la lista de loci del manifiesto está vacía y los
    // dos documentos pasan su esquema; una efectiva sin el enlace a la tabla
    // de alelos, o con el deslizamiento común sin sus limitaciones, no.
    let manifest = large.json("manifest.json");
    let effective = large.json("effective-config.json");
    assert_eq!(manifest["repeat_loci"]["loci"], json!([]));
    assert_eq!(effective["repeat_loci"], json!([]));
    let mut without_link = effective.clone();
    without_link["repeat_detection"]
        .as_object_mut()
        .expect("sección")
        .remove("alleles_table");
    let mut without_limitations = effective.clone();
    without_limitations["repeat_detection"]["slippage"]
        .as_object_mut()
        .expect("deslizamiento")
        .remove("limitations");
    for (document, schema, accepted) in [
        (manifest, "manifest-general-v6.schema.json", true),
        (effective, "effective-config-v6.schema.json", true),
        (without_link, "effective-config-v6.schema.json", false),
        (
            without_limitations,
            "effective-config-v6.schema.json",
            false,
        ),
    ] {
        assert_eq!(
            schema_accepts(&large, &document, schema),
            accepted,
            "{schema}: aceptado = {accepted}"
        );
    }
}

/// Si el validador de esquemas externo (Python) acepta `document`.
fn schema_accepts(package: &Package, document: &Value, schema: &str) -> bool {
    let root = Path::new(env!("CARGO_MANIFEST_DIR"));
    let path = package.root.join("instance.json");
    fs::write(&path, json_document_bytes(document).expect("JSON")).expect("instancia");
    Command::new("python3")
        .arg(root.join("tests/external_jsonschema_validator.py"))
        .arg(root.join("schemas").join(schema))
        .arg(&path)
        .status()
        .expect("validador de esquemas en Python")
        .success()
}
