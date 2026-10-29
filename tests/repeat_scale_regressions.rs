//! Regresiones de corrección de la entrega 2 («escala») del plan de detección
//! de microsatélites: los arreglos de rendimiento tienen que dar **el mismo
//! diagnóstico** que el recorrido cuadrático que sustituyen.
//!
//! La huella de los bytes publicados está en
//! `tests/repeat_scale_frozen_identity.rs`; aquí van los diagnósticos, que no
//! se publican en ningún fichero y la huella no ve.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::diagnostics::Diagnostic;

const SEED: &str = "5e2a9c7b1d3f4068a2c4e6f8091b3d5f7a9c1e3f5072b4d6f8a0c2e4f6a8b0c2";
static NONCE: AtomicU64 = AtomicU64::new(0);

/// Directorio temporal: se borra al soltarse.
struct TemporaryRoot(PathBuf);

impl TemporaryRoot {
    fn new() -> Self {
        let root = std::env::temp_dir().join(format!(
            "gps-escala-regresiones-{}-{}",
            std::process::id(),
            NONCE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir_all(&root).expect("raíz");
        Self(root)
    }
}

impl Drop for TemporaryRoot {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

/// Lee, interpreta y resuelve un documento general como lo hace la línea de
/// órdenes, y devuelve el primer diagnóstico.
fn first_diagnostic(path: &Path, text: &str) -> Diagnostic {
    let admitted =
        match admit_toml_header(path, text.as_bytes().to_vec(), ExpectedConfigKind::General) {
            Ok(admitted) => admitted,
            Err(diagnostic) => return diagnostic,
        };
    let declared = match parse_declared_config(admitted) {
        Ok(declared) => declared,
        Err(diagnostic) => return diagnostic,
    };
    match build_effective_config(&declared) {
        Ok(_) => panic!("el documento tenía que rechazarse"),
        Err(diagnostic) => diagnostic,
    }
}

/// Una referencia suministrada de una secuencia de 2 000 bases, con `(CA)n`
/// plantados donde los loci los declaran, y el documento general con esos
/// loci en el orden dado. Cada locus es `(id, inicio en base 1, fin)`.
fn document_with_loci(root: &Path, loci: &[(&str, usize, usize)]) -> (PathBuf, String) {
    let mut bases = "ACGT".repeat(500).into_bytes();
    for (_, start, end) in loci {
        for (offset, base) in bases[start - 1..*end].iter_mut().enumerate() {
            *base = if offset % 2 == 0 { b'C' } else { b'A' };
        }
    }
    let fasta = root.join("reference.fa");
    fs::write(
        &fasta,
        format!(">chr1\n{}\n", std::str::from_utf8(&bases).expect("ASCII")),
    )
    .expect("FASTA");
    fs::write(
        root.join("population-snvs.tsv"),
        "contig\tposition\tref\talt\tpopulation_id\talt_af\nchr1\t2\tC\tA\tP1\t0.5\n",
    )
    .expect("matriz");
    let mut text = format!(
        "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
         [population_snvs]\nmode = \"supplied\"\ninput_path = \"{}\"\n\n\
         [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"{}\"\n\n\
         [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 2000\n\n",
        root.join("population-snvs.tsv").display(),
        fasta.display()
    );
    for (id, start, end) in loci {
        text.push_str(&format!(
            "[[reference.repeat_loci]]\nid = \"{id}\"\nsequence = \"chr1\"\nstart = {start}\n\
             end = {end}\nmotif = \"CA\"\n\n[[reference.repeat_loci.alleles]]\nunits = {}\n\
             frequency = 1.0\n\n",
            (end - start).div_ceil(2)
        ));
    }
    text.push_str(
        "[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\n\
         genotype_model = \"hwe\"\nindividual_ids = [\"S1\"]\n",
    );
    (root.join("general.toml"), text)
}

fn overlap_detail(diagnostic: &Diagnostic) -> String {
    let rendered = diagnostic.render_human();
    assert!(
        rendered.contains("se solapa con el del locus"),
        "tenía que rechazarse por solape: {rendered}"
    );
    rendered
}

/// **F6.4.a.6. Con tres loci donde el tercero solapa al primero y no al
/// segundo, el rechazo nombra el tercero y el primero, como cuando se
/// comparaba con todos.**
///
/// Una comprobación que mirara solo el vecino **declarado** inmediatamente
/// antes (comparar vecinos sin ordenar) compararía `L3` con `L2`, que no se
/// solapan, y lo dejaría pasar. El control sin solape fija que el documento es
/// válido salvo por ese solape: sin `L3` se resuelve.
#[test]
fn the_overlap_check_names_the_first_declared_locus_that_the_third_overlaps() {
    let root = TemporaryRoot::new();
    let (path, text) = document_with_loci(
        &root.0,
        &[("L1", 101, 120), ("L2", 301, 320), ("L3", 111, 130)],
    );
    let rendered = overlap_detail(&first_diagnostic(&path, &text));
    assert!(
        rendered.contains("L3") && rendered.contains("se solapa con el del locus L1"),
        "{rendered}"
    );

    let (path, text) = document_with_loci(&root.0, &[("L1", 101, 120), ("L2", 301, 320)]);
    let admitted =
        admit_toml_header(&path, text.into_bytes(), ExpectedConfigKind::General).expect("cabecera");
    let declared = parse_declared_config(admitted).expect("documento");
    build_effective_config(&declared)
        .unwrap_or_else(|diagnostic| panic!("control: {}", diagnostic.render_human()));
}

/// **F6.4.a.6. Si el locus nuevo pisa varios anteriores, el rechazo nombra el
/// primero en orden de declaración, no el primero por posición.**
///
/// `L3` pisa `L1` (declarado antes, más a la derecha) y `L2` (declarado
/// después, más a la izquierda). Una búsqueda que devolviera el primero por
/// posición nombraría `L2`.
#[test]
fn the_overlap_check_names_the_earliest_declared_of_several_overlapped_loci() {
    let root = TemporaryRoot::new();
    let (path, text) = document_with_loci(
        &root.0,
        &[("L1", 301, 320), ("L2", 101, 120), ("L3", 111, 310)],
    );
    let rendered = overlap_detail(&first_diagnostic(&path, &text));
    assert!(
        rendered.contains("se solapa con el del locus L1"),
        "{rendered}"
    );
    assert!(!rendered.contains("locus L2"), "{rendered}");
}

/// **F6.4.a.8c. Un error al final de un documento grande, en su primera línea
/// y tras finales de línea `CRLF` se sitúa en la misma línea y columna que
/// contando los saltos desde el principio.**
///
/// Las cifras esperadas se cuentan aquí, sobre el texto escrito por la prueba,
/// y no con el lector: el documento grande tiene 4 000 líneas de relleno
/// (comentarios) antes del valor erróneo, que va en la columna 10 de su línea
/// (contada desde 1).
#[test]
fn a_type_error_is_located_at_the_same_line_and_column_after_a_large_document() {
    let root = TemporaryRoot::new();
    let path = root.0.join("general.toml");
    let header = format!("config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n");
    let padding = "# relleno de la prueba de escala, con una ñ\n".repeat(4_000);

    // Al final de un documento grande: `ploidy = "dos"` en la última línea.
    let text = format!("{header}{padding}ploidy = \"dos\"\n");
    let line = text.lines().count();
    let rendered = first_diagnostic(&path, &text).render_human();
    assert!(
        rendered.contains(&format!("line {line}, column 10,")),
        "esperaba la línea {line}, columna 10: {rendered}"
    );

    // En la primera línea.
    let text = format!("ploidy = \"dos\"\n{header}");
    let rendered = first_diagnostic(&path, &text).render_human();
    assert!(
        rendered.contains("line 1, column 10,"),
        "esperaba la línea 1, columna 10: {rendered}"
    );

    // Con todas las líneas terminadas en CRLF.
    let text = format!("{header}{padding}ploidy = \"dos\"\n").replace('\n', "\r\n");
    let line = text.lines().count();
    let rendered = first_diagnostic(&path, &text).render_human();
    assert!(
        rendered.contains(&format!("line {line}, column 10,")),
        "esperaba la línea {line}, columna 10, con CRLF: {rendered}"
    );
}
