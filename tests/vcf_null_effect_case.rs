//! Un caso de efecto neto nulo existe en la verdad y no tiene registro VCF
//! (excepción de REQ-VER-25, 3 de octubre de 2026).
//!
//! Hasta esa fecha la publicación abortaba con un error interno, «caso sin
//! ALT», en cuanto los eventos de una copia se anulaban: el agrupador exigía
//! un registro para cada caso, y en ese sitio la copia es igual a la
//! referencia, de modo que no hay ALT que publicar. Ahora la verdad, que
//! registra **eventos**, conserva el caso con sus portadores, y el VCF, que
//! registra **secuencia**, no lo lista en ningún `GPS_CASE`.
//!
//! Las dos formas son las que encontró la revisión adversarial del 2 de
//! octubre, con sus ficheros tal cual:
//!
//! - `tests/fixtures/vcf-null-effect/direct-cancel/`: el escenario minimizado
//!   por el revisor. En la copia 1 de `ind-2` una inserción `CT` en el punto
//!   149 y la deleción vigente 148..150 dejan una `T` donde la referencia
//!   tiene una `T`. Con locus de repetición (`general.toml`) y sin él
//!   (`noloci.toml`).
//! - el tracto `(CA)15` declarado en 201..230 de
//!   `tests/fixtures/vcf-fusion/review/v-136/reference.fa`, seguido de un `CA`
//!   más en la referencia: la copia 0 de `ind-1` que recibe el alelo de 16
//!   unidades y además pierde ese flanco por una deleción queda igual a la
//!   referencia. Abortaba en seis de las doce primeras semillas, que son las
//!   que se ejercitan aquí.
//! - `tests/fixtures/vcf-null-effect/inserted-then-deleted/`, el escenario T1
//!   de la revisión del 3 de octubre: una inserción que una deleción
//!   posterior borra entera, en un sitio que ninguna otra copia cambia. La
//!   copia no conserva ningún cambio y no se compone ninguna ventana; la
//!   corrección del 2 de octubre seguía abortando aquí.
//!
//! # Qué hace que estas pruebas puedan fallar
//!
//! Que el juez semántico acepte no basta. Cada publicación pasa además por un
//! oráculo propio, escrito aquí, que reconstruye cada copia desde los `GT` del
//! VCF y la compara con la que reconstruye desde las piezas de la verdad; y
//! cada escenario exige, como control, que la copia portadora quede de verdad
//! igual a la referencia, leído de la verdad y no de nada que decida el VCF.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::vcf_semantic_validation::validate_published_vcf_semantics;
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

fn fixture(path: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(path)
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-null-effect-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

struct Published {
    root: PathBuf,
    out: PathBuf,
}

impl Published {
    fn truth_path(&self) -> PathBuf {
        self.out.join("outputs/truth.jsonl")
    }
    fn vcf_path(&self) -> PathBuf {
        self.out.join("outputs/truth.vcf")
    }
    fn fasta_path(&self) -> PathBuf {
        self.out.join("resources/references/ref-main.fa")
    }
    fn truth(&self) -> Vec<Value> {
        fs::read_to_string(self.truth_path())
            .expect("truth")
            .lines()
            .map(|line| serde_json::from_str(line).expect("JSON line"))
            .collect()
    }
    fn vcf(&self) -> String {
        fs::read_to_string(self.vcf_path()).expect("VCF de texto")
    }
}

impl Drop for Published {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.root).ok();
    }
}

/// Publica `config`, que se escribe en un directorio propio junto con su
/// FASTA y su matriz, con el VCF en texto.
fn publish(label: &str, config: &str, reference: &Path, matrix: &Path) -> Published {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("root");
    fs::copy(reference, root.join("reference.fa")).expect("FASTA");
    fs::copy(matrix, root.join("snvs.tsv")).expect("matrix");
    let mut config = config.replace("mode = \"bgzf\"", "mode = \"text\"");
    if !config.contains("[vcf_output]") {
        config.push_str("\n[vcf_output]\nmode = \"text\"\n");
    }
    let path = root.join("general.toml");
    fs::write(&path, config).expect("configuration");
    let out = root.join("published");
    let effective = load_and_validate_config(&path)
        .unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    execute_validated_config(&effective, &out)
        .unwrap_or_else(|error| panic!("{label}: la publicación abortó: {}", error.render_human()));
    Published { root, out }
}

fn fasta_bases(path: &Path) -> Vec<u8> {
    fs::read_to_string(path)
        .expect("FASTA")
        .lines()
        .filter(|line| !line.starts_with('>'))
        .flat_map(|line| line.trim().bytes())
        .collect()
}

/// Cada copia según la verdad, por (individuo, fase), desde sus piezas.
fn copies_from_truth(truth: &[Value], reference: &[u8]) -> BTreeMap<(String, usize), Vec<u8>> {
    let mut copies = BTreeMap::new();
    for record in truth
        .iter()
        .filter(|record| record["record_type"] == "haplotype_copy")
    {
        let mut sequence = Vec::new();
        for piece in record["current_pieces"].as_array().expect("pieces") {
            match piece["kind"].as_str().expect("kind") {
                "source" => {
                    let start = piece["source"]["start"].as_u64().expect("start") as usize;
                    let end = piece["source"]["end"].as_u64().expect("end") as usize;
                    let mut bases = reference[start..end].to_vec();
                    for snp in piece["snps"].as_array().expect("snps") {
                        let position = snp["position"].as_u64().expect("position") as usize;
                        bases[position - start] =
                            snp["alternate"].as_str().expect("alternate").as_bytes()[0];
                    }
                    sequence.extend(bases);
                }
                "inserted" | "repeat_allele" => {
                    sequence.extend(piece["bases"].as_str().expect("bases").bytes());
                }
                other => panic!("pieza desconocida {other}"),
            }
        }
        copies.insert(
            (
                record["identity"]["individual_id"]
                    .as_str()
                    .expect("individual")
                    .to_owned(),
                usize::try_from(record["phase"].as_u64().expect("phase")).expect("phase"),
            ),
            sequence,
        );
    }
    copies
}

struct Record {
    start: usize,
    reference: String,
    alternates: Vec<String>,
    cases: Vec<String>,
    genotypes: Vec<Vec<usize>>,
}

fn records(text: &str) -> Vec<Record> {
    text.lines()
        .filter(|line| !line.starts_with('#'))
        .map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            Record {
                start: fields[1].parse::<usize>().expect("POS") - 1,
                reference: fields[3].to_owned(),
                alternates: fields[4].split(',').map(str::to_owned).collect(),
                cases: fields[7]
                    .split(';')
                    .find_map(|field| field.strip_prefix("GPS_CASE="))
                    .expect("GPS_CASE")
                    .split(',')
                    .map(str::to_owned)
                    .collect(),
                genotypes: fields[9..]
                    .iter()
                    .map(|gt| {
                        gt.split('|')
                            .map(|allele| allele.parse().expect("GT"))
                            .collect()
                    })
                    .collect(),
            }
        })
        .collect()
}

/// Cada copia según el VCF: la referencia con el alelo que le da el `GT` de
/// cada registro.
fn copies_from_vcf(text: &str, reference: &[u8]) -> BTreeMap<(String, usize), Vec<u8>> {
    let records = records(text);
    let samples = text
        .lines()
        .find(|line| line.starts_with("#CHROM"))
        .expect("cabecera de columnas")
        .split('\t')
        .skip(9)
        .map(str::to_owned)
        .collect::<Vec<_>>();
    let mut copies = BTreeMap::new();
    for (index, sample) in samples.iter().enumerate() {
        let ploidy = records[0].genotypes[index].len();
        for phase in 0..ploidy {
            let mut sequence = Vec::new();
            let mut cursor = 0;
            for record in &records {
                assert!(record.start >= cursor, "dos registros se pisan");
                sequence.extend_from_slice(&reference[cursor..record.start]);
                let allele = record.genotypes[index][phase];
                let bases = if allele == 0 {
                    &record.reference
                } else {
                    &record.alternates[allele - 1]
                };
                sequence.extend_from_slice(bases.as_bytes());
                cursor = record.start + record.reference.len();
            }
            sequence.extend_from_slice(&reference[cursor..]);
            copies.insert((sample.clone(), phase), sequence);
        }
    }
    copies
}

/// El juez semántico acepta y el oráculo propio reconstruye desde el VCF
/// exactamente las copias de la verdad.
fn assert_published_exactly(published: &Published) {
    validate_published_vcf_semantics(
        &published.fasta_path(),
        &published.truth_path(),
        &published.vcf_path(),
    )
    .unwrap_or_else(|error| panic!("el juez semántico rechazó el VCF: {error}"));
    let reference = fasta_bases(&published.fasta_path());
    let expected = copies_from_truth(&published.truth(), &reference);
    let observed = copies_from_vcf(&published.vcf(), &reference);
    assert_eq!(
        observed.len(),
        expected.len(),
        "una copia por fase y muestra"
    );
    for (copy, sequence) in &expected {
        assert!(
            observed.get(copy) == Some(sequence),
            "la copia {copy:?} del VCF no es la de la verdad"
        );
    }
}

/// Casos de la verdad cuyo alelo alternativo lleva la copia dada, con su
/// tipo y su locus.
fn cases_carried_by(
    truth: &[Value],
    individual: &str,
    phase: u64,
) -> Vec<(String, String, u64, u64)> {
    truth
        .iter()
        .filter(|record| record["record_type"] == "genotype")
        .filter(|record| {
            record["carriers"]
                .as_array()
                .expect("carriers")
                .iter()
                .any(|carrier| {
                    carrier["allele_index"].as_u64() != Some(0)
                        && carrier["copy"]["individual_id"] == individual
                        && carrier["copy"]["copy_index"].as_u64() == Some(phase)
                })
        })
        .map(|record| {
            (
                record["case_id"].as_str().expect("case").to_owned(),
                record["case_kind"].as_str().expect("kind").to_owned(),
                record["locus"]["start"].as_u64().expect("start"),
                record["locus"]["end"].as_u64().expect("end"),
            )
        })
        .collect()
}

#[test]
fn two_direct_indels_that_cancel_keep_their_cases_in_the_truth_and_publish_no_record() {
    let directory = fixture("vcf-null-effect/direct-cancel");
    for config in ["general.toml", "noloci.toml"] {
        let text = fs::read_to_string(directory.join(config)).expect("config");
        let published = publish(
            config,
            &text,
            &directory.join("reference.fa"),
            &directory.join("snvs.tsv"),
        );
        let truth = published.truth();
        let reference = fasta_bases(&published.fasta_path());

        // Control: la copia 1 de ind-2 es igual a la referencia en el sitio.
        // La referencia tiene una T en 148, y la copia, una T insertada.
        assert_eq!(
            &reference[147..150],
            b"CTA",
            "{config}: contexto de la referencia"
        );
        let copy = &copies_from_truth(&truth, &reference)[&("ind-2".to_owned(), 1)];
        assert_eq!(
            &copy[140..160],
            &reference[140..160],
            "{config}: la copia queda igual a la referencia en el sitio"
        );

        // La verdad registra los dos eventos: los dos casos siguen ahí, con
        // la copia como portadora y en coordenadas de la referencia.
        let mut indels = cases_carried_by(&truth, "ind-2", 1)
            .into_iter()
            .filter(|(_, kind, ..)| kind == "indel")
            .collect::<Vec<_>>();
        indels.sort_by_key(|(_, _, start, end)| (*start, *end));
        let loci = indels
            .iter()
            .map(|(_, _, start, end)| (*start, *end))
            .collect::<Vec<_>>();
        assert_eq!(
            loci,
            [(148, 149), (149, 149)],
            "{config}: deleción e inserción"
        );

        // El VCF registra secuencia: ningún registro los nombra ni pisa el
        // sitio.
        let vcf = published.vcf();
        for record in records(&vcf) {
            for (case, ..) in &indels {
                assert!(
                    !record.cases.contains(case),
                    "{config}: un caso de efecto nulo no entra en GPS_CASE"
                );
            }
            assert!(
                record.start + record.reference.len() <= 147 || record.start >= 150,
                "{config}: ningún registro en 148..150 (base cero)"
            );
        }
        assert_published_exactly(&published);
    }
}

/// **Una inserción que una deleción posterior borra entera, en un sitio que
/// ninguna otra copia cambia.**
///
/// `tests/fixtures/vcf-null-effect/inserted-then-deleted/` es el escenario T1
/// de la revisión adversarial del 3 de octubre de 2026: la copia 0 de `ind-1`
/// recibe `GT` en el punto 150 y después la deleción vigente 150..152, que se
/// lleva justo esas dos bases. Las piezas finales son la referencia entera:
/// la copia no tiene ningún cambio, y como ninguna otra copia cambia ahí
/// tampoco se compone ninguna ventana. Hasta ese día la publicación abortaba
/// con «caso sin ALT», y el oráculo, sin la guarda del productor, lo
/// rechazaba en falso.
#[test]
fn an_insertion_deleted_entirely_where_no_copy_changes_publishes_its_deletion_case() {
    let directory = fixture("vcf-null-effect/inserted-then-deleted");
    let text = fs::read_to_string(directory.join("general.toml")).expect("config");
    let published = publish(
        "inserted-then-deleted",
        &text,
        &directory.join("reference.fa"),
        &directory.join("snvs.tsv"),
    );
    let truth = published.truth();
    let reference = fasta_bases(&published.fasta_path());

    // Control: la copia es exactamente la referencia, salvo el cambio de base
    // de la matriz en la posición 20 si le tocó.
    let copy = &copies_from_truth(&truth, &reference)[&("ind-1".to_owned(), 0)];
    assert_eq!(
        copy.len(),
        reference.len(),
        "la copia no cambia de longitud"
    );
    assert_eq!(
        &copy[30..],
        &reference[30..],
        "la copia es la referencia desde 30"
    );

    // La verdad conserva la deleción, con locus vacío en 150; la inserción,
    // borrada entera, ya no es un caso.
    let carried = cases_carried_by(&truth, "ind-1", 0)
        .into_iter()
        .filter(|(_, kind, ..)| kind == "indel")
        .collect::<Vec<_>>();
    assert_eq!(
        carried
            .iter()
            .map(|(_, _, start, end)| (*start, *end))
            .collect::<Vec<_>>(),
        [(150, 150)],
        "solo la deleción es un caso"
    );
    let vcf = published.vcf();
    for record in records(&vcf) {
        assert!(
            !record.cases.contains(&carried[0].0),
            "el caso de efecto nulo no entra en GPS_CASE"
        );
        assert!(
            record.start + record.reference.len() < 140 || record.start > 160,
            "ningún registro cerca de 150"
        );
    }
    assert_published_exactly(&published);
}

/// Las seis semillas en que la revisión del 2 de octubre vio abortar el
/// escenario, en hexadecimal de 64 cifras.
const TRACT_SEEDS: [u64; 6] = [3, 4, 5, 7, 11, 12];

fn tract_config(seed: u64) -> String {
    let base = fs::read_to_string(fixture("vcf-fusion/review/v-136/general.toml")).expect("config");
    let head = &base[..base
        .find("[[reference.repeat_loci.alleles]]")
        .expect("alleles")];
    let tail = &base[base.find("[[initial_populations]]").expect("populations")
        ..base.find("[[direct_events]]").expect("events")];
    let head = head
        .lines()
        .map(|line| {
            if line.starts_with("seed = ") {
                format!("seed = \"{seed:064x}\"")
            } else {
                line.to_owned()
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    format!(
        "{head}\n[[reference.repeat_loci.alleles]]\nunits = 15\nfrequency = 0.5\n\n\
         [[reference.repeat_loci.alleles]]\nunits = 16\nfrequency = 0.5\n\n{tail}\
         [[direct_events]]\nevent_id = \"e0\"\ncontig = \"chr1\"\nindividual_id = \"ind-1\"\n\
         phase = 0\nkind = \"deletion\"\nstart = 230\nend = 232\n"
    )
}

#[test]
fn a_repeat_allele_and_a_flank_deletion_that_rebuild_the_reference_publish_exactly() {
    let directory = fixture("vcf-fusion/review/v-136");
    for seed in TRACT_SEEDS {
        let published = publish(
            &format!("tract-{seed}"),
            &tract_config(seed),
            &directory.join("reference.fa"),
            &directory.join("snvs.tsv"),
        );
        let truth = published.truth();
        let reference = fasta_bases(&published.fasta_path());

        // Control: en esta semilla la copia 0 de ind-1 lleva el alelo de 16
        // unidades y la deleción del flanco, y queda igual a la referencia.
        assert_eq!(
            &reference[200..232],
            "CA".repeat(16).as_bytes(),
            "el tracto declarado (CA)15 va seguido de un CA más"
        );
        let carried = cases_carried_by(&truth, "ind-1", 0);
        let repeat = carried
            .iter()
            .find(|(_, kind, ..)| kind == "repeat")
            .unwrap_or_else(|| panic!("semilla {seed}: la copia 0 de ind-1 realiza el tracto"));
        let deletion = carried
            .iter()
            .find(|(_, kind, start, end)| kind == "indel" && (*start, *end) == (230, 232))
            .unwrap_or_else(|| panic!("semilla {seed}: la copia 0 de ind-1 pierde el flanco"));
        let copy = &copies_from_truth(&truth, &reference)[&("ind-1".to_owned(), 0)];
        assert_eq!(
            &copy[190..240],
            &reference[190..240],
            "semilla {seed}: la copia queda igual a la referencia"
        );

        // La deleción no aporta ALT en ningún registro; el caso de repetición
        // sigue en el registro del tracto por las demás copias que lo
        // realizan, pero la copia 0 de ind-1 lleva allí la referencia.
        let vcf = published.vcf();
        let samples = vcf
            .lines()
            .find(|line| line.starts_with("#CHROM"))
            .expect("columnas")
            .split('\t')
            .skip(9)
            .collect::<Vec<_>>();
        let ind_1 = samples
            .iter()
            .position(|sample| *sample == "ind-1")
            .expect("ind-1");
        for record in records(&vcf) {
            assert!(
                !record.cases.contains(&deletion.0),
                "semilla {seed}: la deleción de efecto nulo no entra en GPS_CASE"
            );
            if record.cases.contains(&repeat.0)
                && record.start <= 200
                && 232 <= record.start + record.reference.len()
            {
                assert_eq!(
                    record.genotypes[ind_1][0], 0,
                    "semilla {seed}: la copia 0 de ind-1 lleva la referencia en el registro del tracto"
                );
            }
        }
        assert_published_exactly(&published);
    }
}
