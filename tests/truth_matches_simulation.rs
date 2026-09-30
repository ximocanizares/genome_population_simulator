//! La verdad publicada describe las bases que se simularon.
//!
//! Hasta el 30 de septiembre de 2026, cada alelo de *indel* copiaba la
//! secuencia final de su copia, materializada desde el estado interno, y el
//! validador la contrastaba con la copia que reconstruía desde la verdad. Ese
//! era, de rebote, el único contraste entre las dos. Una revisión lo demostró
//! con un mutante en la conversión de bases de la proyección —si el evento
//! empieza por `ins`, cambia T por G—: sin el contraste, la verdad y el VCF
//! publicaban G mientras las lecturas llevaban T, y toda la batería seguía en
//! verde.

use std::fs;
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;
use serde_json::Value;

static NONCE: AtomicU64 = AtomicU64::new(0);

const INSERTED: &str = "TTGT";
const POINT: usize = 1200;

/// **Las bases insertadas que declara la configuración son las que publican
/// la verdad y el VCF.**
///
/// La verdad independiente es la configuración: las bases se escriben aquí, y
/// no se leen de nada que el programa produzca. Se exigen en los tres sitios
/// de la verdad donde aparecen —el fragmento vigente, el historial y el alelo
/// del caso— y en el `ALT` del VCF. El evento se llama `ins-t` y sus bases
/// llevan T a propósito: es el caso que el mutante de la revisión corrompía.
///
/// Con el contraste en memoria de la publicación, ese mutante ya no llega aquí:
/// la ejecución falla antes. Sin él, lo caza la comparación con lo declarado.
#[test]
fn published_truth_and_vcf_carry_the_declared_inserted_bases() {
    let root = temporary("truth-matches");
    fs::create_dir_all(&root).expect("fixture root");
    let config = write_fixture(&root);
    let destination = root.join("published");
    let effective = load_and_validate_config(&config).expect("valid configuration");
    execute_validated_config(&effective, &destination).expect("complete run");

    let truth = fs::read_to_string(destination.join("outputs/truth.jsonl")).expect("truth");
    let records = truth
        .lines()
        .map(|line| serde_json::from_str::<Value>(line).expect("JSON line"))
        .collect::<Vec<_>>();
    let mut seen = Vec::new();
    for record in &records {
        match record["record_type"].as_str() {
            Some("haplotype_copy") => {
                for piece in record["current_pieces"].as_array().expect("pieces") {
                    if piece["kind"] == "inserted" {
                        seen.push(("current piece", piece["bases"].clone()));
                    }
                }
                for event in record["history"].as_array().expect("history") {
                    if event["kind"] == "insertion" {
                        seen.push(("history", event["bases"].clone()));
                    }
                }
            }
            Some("genotype") => {
                for allele in record["alleles"].as_array().expect("alleles") {
                    if allele["kind"] == "insertion" {
                        seen.push(("allele", allele["bases"].clone()));
                    }
                }
            }
            _ => {}
        }
    }
    assert_eq!(
        seen.iter().map(|(place, _)| *place).collect::<Vec<_>>(),
        ["current piece", "history", "allele"],
        "one insertion, published once in each place"
    );
    for (place, bases) in &seen {
        assert_eq!(bases, INSERTED, "{place}");
    }

    // El VCF alinea la inserción a la izquierda, de modo que puede publicarla
    // rotada dentro de una repetición: lo que se exige es una rotación de las
    // bases declaradas, detrás de un ancla de una base.
    let vcf = fs::read_to_string(destination.join("outputs/truth.vcf")).expect("VCF");
    let insertions = vcf
        .lines()
        .filter(|line| !line.starts_with('#'))
        .flat_map(|line| {
            let fields = line.split('\t').collect::<Vec<_>>();
            let reference = fields[3].to_owned();
            fields[4]
                .split(',')
                .filter(|alternate| alternate.len() == reference.len() + INSERTED.len())
                .map(|alternate| (reference.clone(), alternate.to_owned()))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    assert_eq!(insertions.len(), 1, "one insertion record: {insertions:?}");
    let (anchor, alternate) = &insertions[0];
    assert_eq!(anchor.len(), 1, "{insertions:?}");
    assert!(alternate.starts_with(anchor.as_str()), "{insertions:?}");
    assert!(
        format!("{INSERTED}{INSERTED}").contains(&alternate[1..]),
        "the VCF must publish the declared insertion: {insertions:?}"
    );

    fs::remove_dir_all(&root).expect("fixture cleanup");
}

fn write_fixture(root: &Path) -> PathBuf {
    let bases = pseudo_random_bases(2000);
    let mut fasta = String::from(">chr1\n");
    for chunk in bases.chunks(60) {
        fasta.push_str(std::str::from_utf8(chunk).expect("ASCII bases"));
        fasta.push('\n');
    }
    fs::write(root.join("reference.fa"), fasta).expect("reference FASTA");
    let site = 100_usize;
    let reference_base = char::from(bases[site - 1]);
    let alternate_base = if reference_base == 'A' { 'C' } else { 'A' };
    fs::write(
        root.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n\
             chr1\t{site}\t{reference_base}\t{alternate_base}\tP1\t0.5\n"
        ),
    )
    .expect("population SNVs");
    let config = root.join("general.toml");
    fs::write(
        &config,
        format!(
            "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{seed}\"\n\n\
             [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
             [reference]\nkind = \"supplied\"\nid = \"ref-main\"\nfasta = \"reference.fa\"\n\n\
             [[reference.sequences]]\nid = \"chr1\"\nclass = \"chromosome\"\nlength = 2000\n\n\
             [[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = 1\nploidy = 2\n\
             genotype_model = \"hwe\"\nindividual_ids = [\"ind-1\"]\n\n\
             [[direct_events]]\nevent_id = \"ins-t\"\nindividual_id = \"ind-1\"\nphase = 0\n\
             kind = \"insertion\"\ncontig = \"chr1\"\npoint = {POINT}\nbases = \"{INSERTED}\"\n\n\
             [vcf_output]\nmode = \"text\"\n",
            seed = "5".repeat(64),
        ),
    )
    .expect("configuration");
    config
}

fn pseudo_random_bases(length: usize) -> Vec<u8> {
    let mut state = 0x2545_F491_4F6C_DD1D_u64;
    (0..length)
        .map(|_| {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            b"ACGT"[usize::try_from(state >> 62).expect("two bits")]
        })
        .collect()
}

fn temporary(label: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-truth-matches-{}-{}-{label}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed)
    ))
}
