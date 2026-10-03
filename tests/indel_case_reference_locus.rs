//! El locus de un caso de inserción o deleción se publica en coordenadas de la
//! **referencia**, también cuando el evento se declaró después de otro *indel*
//! de la misma copia.
//!
//! Los cambios directos se declaran en coordenadas **vigentes** de la copia, y
//! así los conserva el historial. Hasta el 02/10/2026 el locus del caso las
//! copiaba tal cual: una inserción declarada en el punto vigente 150 detrás de
//! una inserción de dos bases quedaba, en la verdad, en la base 150 de la
//! referencia y no en la 148, y el VCF abortaba con «un ALT sin caso». El juez
//! de la verdad no lo veía porque nunca comprobaba ese locus.
//!
//! # Qué hace que estas pruebas puedan fallar
//!
//! Los loci esperados están **calculados a mano** desde el escenario, en la
//! tabla de [`EXPECTED_LOCI`], y no se leen de nada que produzca el programa.
//! Un error común al productor y al juez —que los dos tomaran las coordenadas
//! vigentes— lo caza esa tabla por el contenido publicado, aunque los jueces
//! aceptaran. La segunda prueba altera a propósito un locus publicado y exige
//! que el juez lo rechace **por ese motivo**.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use genome_population_simulator::truth_validation::TruthValidationError;
use genome_population_simulator::vcf_semantic_validation::{
    VcfSemanticValidationError, validate_published_vcf_semantics,
};
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const SEED: &str = "7777777777777777777777777777777777777777777777777777777777777777";

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-indel-locus-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}

/// Bases pseudoaleatorias reproducibles, sin depender del generador del
/// programa.
fn pseudo_random_bases(length: usize, mut state: u64) -> Vec<u8> {
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
        })
        .collect()
}

/// Referencia de 400 pb con las bases vecinas de cada evento fijadas para
/// que ningún *indel* pueda alinearse a la izquierda: así el registro VCF de
/// cada uno está donde dice la tabla, sin depender de la normalización.
fn reference() -> Vec<u8> {
    let mut bases = pseudo_random_bases(400, 0x1234_5678_9ABC_DEF1);
    // e0 y e3, inserciones en el punto 100 de la referencia.
    bases[99] = b'A';
    bases[100] = b'C';
    // e1 y e4, inserción en el punto 148.
    bases[147] = b'C';
    bases[148] = b'G';
    // e2, deleción de [195, 198).
    bases[194] = b'T';
    bases[195..198].copy_from_slice(b"ACG");
    bases[198] = b'T';
    // e6, deleción de [297, 300).
    bases[296] = b'T';
    bases[297..300].copy_from_slice(b"CAG");
    bases[300] = b'T';
    // e8 a e11, en la junta que deja e8 en 340..345.
    bases[339] = b'T';
    bases[340..347].copy_from_slice(b"ACATCAT");
    bases[347] = b'C';
    bases
}

/// Los eventos de la copia de fase 1 de `ind-1`, **en coordenadas vigentes**
/// y en el orden en que se aplican. El comentario de cada uno da el estado de
/// la copia que hace falta para traducirlo a la referencia.
const EVENTS: &str = r#"
# e0. Sin indels previos: vigente = referencia. Punto 100.
[[direct_events]]
event_id = "e0"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 100
bases = "GG"

# e1. Tras GG en 100, la vigente c >= 102 es la referencia c - 2: punto 148.
[[direct_events]]
event_id = "e1"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 150
bases = "TTT"

# e2. Tras GG y TTT, la vigente c >= 153 es la referencia c - 5: [195, 198).
[[direct_events]]
event_id = "e2"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 200
end = 203

# e3. Entre las dos G de e0: la última base de referencia que la precede es
# la 99, así que se ancla en 100, como e0.
[[direct_events]]
event_id = "e3"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 101
bases = "A"

# e4. Copia: ref[0,100) GAG ref[100,148) TTT ...; TTT ocupa las vigentes
# [151, 154). Elimina la T central: ninguna base de referencia, así que su
# locus es el ancla vacía de su inicio, 148.
[[direct_events]]
event_id = "e4"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 152
end = 153

# e5. Copia: ... ref[148,195) en [153, 200), ref[198,400) desde 200. La
# vigente 300 es la referencia 298.
[[direct_events]]
event_id = "e5"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 300
bases = "CC"

# e6. Elimina la vigente 299 (referencia 297), las dos C de e5 (300 y 301) y
# las vigentes 302 y 303 (referencia 298 y 299): [297, 300). e5 desaparece
# entera y deja de ser un caso.
[[direct_events]]
event_id = "e6"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 299
end = 304

# e7. Copia: ... ref[148,195) en [153, 200), ref[198,297) desde 200. El punto
# vigente 200 es la junta que dejó e2: la última base de referencia que lo
# precede es la 194, así que se ancla en 195 y no en 198, la siguiente.
[[direct_events]]
event_id = "e7"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 200
bases = "GC"

# e8. Copia: ... ref[198,297) en [202, 301), ref[300,400) desde 301: la
# vigente c >= 301 es la referencia c - 1. Borra [340, 345).
[[direct_events]]
event_id = "e8"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 341
end = 346

# e9. GGG en la junta que dejó e8 (vigente 341): la última base de
# referencia que la precede es la 339, así que se ancla en 340 y no en 345.
[[direct_events]]
event_id = "e9"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "insertion"
point = 341
bases = "GGG"

# e10. Borra la primera G de e9, una deleción de solo bases insertadas en una
# junta: su locus es el ancla izquierda vacía, 340, no la derecha, 345.
[[direct_events]]
event_id = "e10"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 341
end = 342

# e11. Copia: ref[300,340) en [301, 341), GG en [341, 343), ref[345,400)
# desde 343. Borra la segunda G y las referencias 345 y 346: una deleción
# mixta que empieza en bases insertadas en una junta. Su locus empieza en la
# primera base de referencia que borra, 345, no en el ancla, 340.
[[direct_events]]
event_id = "e11"
contig = "chr1"
individual_id = "ind-1"
phase = 1
kind = "deletion"
start = 342
end = 345
"#;

/// Locus esperado de cada caso, en coordenadas de referencia de base cero y
/// semiabiertas, calculado a mano en los comentarios de [`EVENTS`].
const EXPECTED_LOCI: [(&str, u64, u64); 11] = [
    ("e0", 100, 100),
    ("e1", 148, 148),
    ("e2", 195, 198),
    ("e3", 100, 100),
    ("e4", 148, 148),
    ("e6", 297, 300),
    ("e7", 195, 195),
    ("e8", 340, 345),
    ("e9", 340, 340),
    ("e10", 340, 340),
    ("e11", 345, 347),
];

/// Registros VCF esperados de `ind-1` (POS de base uno, REF, ALT) y los
/// eventos cuyos casos tiene que listar cada uno en `GPS_CASE`.
fn expected_records(bases: &[u8]) -> Vec<(usize, String, String, Vec<&'static str>)> {
    let text = |start: usize, end: usize| String::from_utf8(bases[start..end].to_vec()).unwrap();
    vec![
        // Ancla A (base 99) + G A G.
        (100, "A".into(), "AGAG".into(), vec!["e0", "e3"]),
        // Ancla C (base 147) + las dos T que quedan.
        (148, "C".into(), "CTT".into(), vec!["e1", "e4"]),
        // e2 quita ACG y e7 pone GC en la junta: sin ancla, porque las dos
        // longitudes son distintas de cero.
        (196, text(195, 198), "GC".into(), vec!["e2", "e7"]),
        // e8 y e11 quitan ACATCAT (340..347); de las tres G de e9 queda una.
        (
            341,
            "ACATCAT".into(),
            "G".into(),
            vec!["e8", "e9", "e10", "e11"],
        ),
        (297, text(296, 300), "T".into(), vec!["e6"]),
    ]
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
}

impl Drop for Published {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.root).ok();
    }
}

fn publish(label: &str, vcf: bool) -> (Published, Vec<u8>) {
    let root = temporary(label);
    fs::create_dir_all(&root).expect("fixture root");
    let bases = reference();
    let mut fasta = String::from(">chr1\n");
    for chunk in bases.chunks(60) {
        fasta.push_str(std::str::from_utf8(chunk).expect("ASCII"));
        fasta.push('\n');
    }
    fs::write(root.join("reference.fa"), fasta).expect("FASTA");
    let snv_reference = char::from(bases[19]);
    let snv_alternate = if snv_reference == 'A' { 'C' } else { 'A' };
    fs::write(
        root.join("snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t20\t{snv_reference}\t{snv_alternate}\tP1\t0.5\n"
        ),
    )
    .expect("matrix");
    let config = format!(
        "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
         [population_snvs]\nmode = \"supplied\"\ninput_path = \"snvs.tsv\"\n\n\
         [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
         [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = {}\n\n\
         [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 2\nploidy = 2\n\
         genotype_model = \"hwe\"\nindividual_ids = [\"ind-1\", \"ind-2\"]\n\n\
         {}{EVENTS}",
        bases.len(),
        if vcf {
            "[vcf_output]\nmode = \"text\"\n"
        } else {
            ""
        }
    );
    let path = root.join("general.toml");
    fs::write(&path, config).expect("configuration");
    let out = root.join("published");
    let effective = load_and_validate_config(&path)
        .unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    execute_validated_config(&effective, &out)
        .unwrap_or_else(|error| panic!("{label}: {}", error.render_human()));
    (Published { root, out }, bases)
}

fn truth_lines(path: &Path) -> Vec<Value> {
    fs::read_to_string(path)
        .expect("truth")
        .lines()
        .map(|line| serde_json::from_str(line).expect("JSON line"))
        .collect()
}

/// Los casos indel de `ind-1` por el evento de su alelo alternativo.
fn indel_cases(lines: &[Value]) -> Vec<(String, String, u64, u64)> {
    lines
        .iter()
        .filter(|line| line["record_type"] == "genotype" && line["case_kind"] == "indel")
        .map(|line| {
            assert_eq!(line["individual_id"], "ind-1", "solo ind-1 lleva indels");
            (
                line["alleles"][1]["event_id"]
                    .as_str()
                    .expect("event_id")
                    .to_owned(),
                line["case_id"].as_str().expect("case_id").to_owned(),
                line["locus"]["start"].as_u64().expect("start"),
                line["locus"]["end"].as_u64().expect("end"),
            )
        })
        .collect()
}

/// Los casos indel publicados son exactamente los de la tabla.
fn assert_expected_loci(cases: &[(String, String, u64, u64)]) {
    let observed = cases
        .iter()
        .map(|(event, _, start, end)| (event.as_str(), *start, *end))
        .collect::<Vec<_>>();
    let mut expected = EXPECTED_LOCI.to_vec();
    expected.sort();
    assert_eq!(
        observed, expected,
        "cada caso indel en coordenadas de la referencia; e5 no es un caso"
    );
}

/// Sin VCF solo queda el juez de la verdad. Si productor y juez compartieran
/// el error —los dos con las coordenadas vigentes—, la ejecución publicaría
/// y solo la tabla escrita a mano lo vería.
#[test]
fn without_vcf_the_truth_still_publishes_reference_loci() {
    let (published, _) = publish("reference-locus-no-vcf", false);
    assert!(!published.vcf_path().exists(), "control: no se pidió VCF");
    let mut cases = indel_cases(&truth_lines(&published.truth_path()));
    cases.sort();
    assert_expected_loci(&cases);
}

#[test]
fn an_indel_declared_after_another_indel_publishes_its_case_in_reference_coordinates() {
    let (published, bases) = publish("reference-locus", true);
    let lines = truth_lines(&published.truth_path());
    let mut cases = indel_cases(&lines);
    cases.sort();
    assert_expected_loci(&cases);

    // Control: la tabla discrimina. Diez de los once casos tienen
    // coordenadas vigentes distintas de las de referencia, de modo que
    // copiarlas —el defecto— no pasaría.
    let source_points = lines
        .iter()
        .filter(|line| line["record_type"] == "genotype" && line["case_kind"] == "indel")
        .filter(|line| {
            let event = &line["alleles"][1]["source_event"];
            let start = event
                .get("point")
                .or_else(|| event.get("start"))
                .and_then(Value::as_u64);
            start != line["locus"]["start"].as_u64()
        })
        .count();
    assert_eq!(source_points, 10, "solo e0 se declara sin desplazar");

    // El VCF publica cada caso en el registro de su sitio, y los dos jueces
    // aceptan lo publicado.
    validate_published_vcf_semantics(
        &published.fasta_path(),
        &published.truth_path(),
        &published.vcf_path(),
    )
    .unwrap_or_else(|error| panic!("el juez semántico rechazó el VCF: {error}"));
    let vcf = fs::read_to_string(published.vcf_path()).expect("VCF");
    let case_of = |event: &str| {
        cases
            .iter()
            .find(|(candidate, ..)| candidate == event)
            .map(|(_, case, ..)| case.clone())
            .expect("case of event")
    };
    for (position, reference, alternate, events) in expected_records(&bases) {
        let line = vcf
            .lines()
            .filter(|line| !line.starts_with('#'))
            .find(|line| {
                let fields = line.split('\t').collect::<Vec<_>>();
                fields[1] == position.to_string() && fields[3] == reference
            })
            .unwrap_or_else(|| panic!("falta el registro en {position}:\n{vcf}"));
        let fields = line.split('\t').collect::<Vec<_>>();
        assert_eq!(fields[4], alternate, "ALT del registro en {position}");
        assert_eq!(fields[9], "0|1", "ind-1 lo lleva en la fase 1");
        let mut expected_cases = events
            .iter()
            .map(|event| case_of(event))
            .collect::<Vec<_>>();
        expected_cases.sort();
        assert_eq!(
            fields[7],
            format!("GPS_CASE={}", expected_cases.join(",")),
            "casos del registro en {position}"
        );
    }
}

/// Reordena los genotipos como exige el orden canónico: referencia,
/// individuo, secuencia, inicio, fin y caso.
fn canonical_genotype_order(lines: &mut [Value]) {
    let first_genotype = lines
        .iter()
        .position(|line| line["record_type"] == "genotype")
        .expect("genotypes");
    lines[first_genotype..].sort_by_key(|line| {
        (
            line["reference_id"].as_str().unwrap_or_default().to_owned(),
            line["individual_id"]
                .as_str()
                .unwrap_or_default()
                .to_owned(),
            line["locus"]["sequence_id"]
                .as_str()
                .unwrap_or_default()
                .to_owned(),
            line["locus"]["start"].as_u64(),
            line["locus"]["end"].as_u64(),
            line["case_id"].as_str().unwrap_or_default().to_owned(),
        )
    });
}

#[test]
fn the_truth_judge_rejects_a_case_locus_in_copy_coordinates() {
    let (published, _) = publish("tampered-locus", true);
    let original = truth_lines(&published.truth_path());
    // Cada evento desplazado, con el locus que tendría si se copiaran sus
    // coordenadas vigentes —el defecto que esta prueba tiene que cazar—, y e7
    // anclada a la derecha de la junta, la otra regla posible.
    for (event, start, end) in [
        ("e1", 150, 150),
        ("e2", 200, 203),
        ("e3", 101, 101),
        ("e4", 152, 153),
        ("e6", 299, 304),
        ("e7", 198, 198),
        // Solo cambia el final: el juez tiene que comparar los dos extremos.
        ("e6", 297, 299),
        // Las dos reglas posibles que no son la del contrato.
        ("e10", 345, 345),
        ("e11", 340, 347),
    ] {
        let mut lines = original.clone();
        let target = lines
            .iter_mut()
            .find(|line| {
                line["record_type"] == "genotype" && line["alleles"][1]["event_id"] == event
            })
            .expect("case of the event");
        target["locus"]["start"] = start.into();
        target["locus"]["end"] = end.into();
        canonical_genotype_order(&mut lines);
        let tampered = published.out.join(format!("tampered-{event}.jsonl"));
        let mut text = String::new();
        for line in &lines {
            text.push_str(&serde_json::to_string(line).expect("JSON"));
            text.push('\n');
        }
        fs::write(&tampered, text).expect("tampered truth");
        let error = validate_published_vcf_semantics(
            &published.fasta_path(),
            &tampered,
            &published.vcf_path(),
        )
        .expect_err("un locus en coordenadas de la copia no puede aceptarse");
        assert!(
            matches!(
                &error,
                VcfSemanticValidationError::Truth(TruthValidationError::GenotypeMismatch {
                    detail: "el locus del caso no es el que su evento ocupa en la referencia",
                    ..
                })
            ),
            "{event}: rechazo inesperado {error:?}"
        );
    }
    // Control: el documento sin alterar, reescrito por el mismo camino, se
    // acepta. Así el rechazo no lo provoca la reescritura.
    let mut lines = original;
    canonical_genotype_order(&mut lines);
    let rewritten = published.out.join("rewritten.jsonl");
    let mut text = String::new();
    for line in &lines {
        text.push_str(&serde_json::to_string(line).expect("JSON"));
        text.push('\n');
    }
    fs::write(&rewritten, text).expect("rewritten truth");
    validate_published_vcf_semantics(&published.fasta_path(), &rewritten, &published.vcf_path())
        .expect("la verdad publicada, reescrita, sigue siendo válida");
}
