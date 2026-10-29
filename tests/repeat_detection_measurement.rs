//! Medición de la entrega 1 del plan de detección de microsatélites (paso
//! F6.4.a.5).
//!
//! **No es una prueba de la batería.** Todas las funciones llevan `#[ignore]`
//! y leen sus entradas de variables de entorno; sirven para producir las
//! cifras de la sección 11 del plan con el binario optimizado, por ejemplo:
//!
//! ```text
//! cargo test --release --locked --test repeat_detection_measurement --no-run
//! GPS_SSR_FASTA=<fasta> /usr/bin/time -v <binario> --ignored --exact \
//!     measure_detection_stages --nocapture
//! ```
//!
//! Variables:
//!
//! - `GPS_SSR_FASTA`: el FASTA que se barre (todas las funciones salvo
//!   `measure_publication`).
//! - `GPS_SSR_DIR`, `GPS_SSR_LOCI`, `GPS_SSR_INDIVIDUALS`, `GPS_SSR_STUTTER`
//!   (`0`/`1`): `prepare_publication` escribe en `GPS_SSR_DIR` un documento
//!   general con `GPS_SSR_LOCI` candidatos repartidos por todo el genoma.
//! - `GPS_SSR_TOML`, `GPS_SSR_DEST`: `measure_publication` lee ese documento
//!   y publica en `GPS_SSR_DEST` (que no debe existir) por la biblioteca,
//!   saltando el techo de 64 KiB del lector de la línea de órdenes, igual que
//!   el sondeo de la sección 2.3 del plan. No cambia nada de `src/`.
//!
//! Preparar y publicar son dos procesos aparte a propósito: así la memoria
//! residente máxima de `measure_publication` es la del simulador y no la del
//! barrido que eligió los loci.

use std::fmt::Write as _;
use std::fs;
use std::path::{Path, PathBuf};
use std::time::Instant;

use genome_population_simulator::adapters::{admit_toml_header, parse_declared_config};
use genome_population_simulator::application::build_effective_config;
use genome_population_simulator::configuration::ExpectedConfigKind;
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config,
    preflight_validated_config_with_authorizations,
};
use genome_population_simulator::repeat_detection::{
    DEFAULT_PERIODS, DETECTION_MOTIF_LENGTHS, DetectionCounts, DetectionCriteria, DetectionStage,
    PeriodRange, RepeatCandidate, ScanParams, detect_repeat_candidates, find_tandem_repeats,
};
use genome_population_simulator::tandem_repeat::{RepeatMotif, observe_tract};

const SEED: &str = "1d0c5a7e9b3f42618a0d2c4e6f8091a2b3c4d5e6f708192a3b4c5d6e7f809102";

fn env(name: &str) -> String {
    std::env::var(name).unwrap_or_else(|_| panic!("falta la variable de entorno {name}"))
}

fn env_u64(name: &str) -> u64 {
    env(name)
        .parse()
        .unwrap_or_else(|_| panic!("{name} debe ser un entero"))
}

/// Lector mínimo de FASTA: nombre (primera palabra de la cabecera) y bases en
/// mayúsculas. La referencia del simulador ya está en mayúsculas.
fn read_fasta(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = fs::read(path).unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let mut sequences: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.split(|&byte| byte == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if let Some(header) = line.strip_prefix(b">") {
            let name = String::from_utf8_lossy(header)
                .split_whitespace()
                .next()
                .unwrap_or_default()
                .to_owned();
            sequences.push((name, Vec::new()));
        } else if let Some((_, bases)) = sequences.last_mut() {
            bases.extend(line.iter().map(u8::to_ascii_uppercase));
        }
    }
    sequences
}

fn total_bases(sequences: &[(String, Vec<u8>)]) -> u64 {
    sequences.iter().map(|(_, bases)| bases.len() as u64).sum()
}

/// Candidatos de todo el genoma con su índice de secuencia.
fn detect_all(
    sequences: &[(String, Vec<u8>)],
    criteria: &DetectionCriteria,
) -> (Vec<(usize, RepeatCandidate)>, DetectionCounts) {
    let mut counts = DetectionCounts::default();
    let mut candidates = Vec::new();
    for (index, (_, bases)) in sequences.iter().enumerate() {
        let detection = detect_repeat_candidates(bases, criteria);
        counts.add(detection.counts());
        candidates.extend(
            detection
                .candidates()
                .iter()
                .cloned()
                .map(|candidate| (index, candidate)),
        );
    }
    (candidates, counts)
}

fn stage_table(counts: &DetectionCounts) -> String {
    let mut out = String::from("etapa\tp1\tp2\tp3\tp4\tp5\tp6\ttotal\n");
    let mut row = |label: &str, value: &dyn Fn(usize) -> u64| {
        let values: Vec<u64> = (1..=DETECTION_MOTIF_LENGTHS).map(value).collect();
        let total: u64 = values.iter().sum();
        let cells: Vec<String> = values.iter().map(u64::to_string).collect();
        let _ = writeln!(out, "{label}\t{}\t{total}", cells.join("\t"));
    };
    row("detectados", &|length| counts.detected(length));
    for stage in DetectionStage::ALL {
        row(stage.name(), &|length| counts.discarded(stage, length));
    }
    row("aceptados", &|length| counts.accepted(length));
    out
}

/// Punto 2: solo cargar, para restar su memoria del barrido.
#[test]
#[ignore = "medición: necesita GPS_SSR_FASTA"]
fn measure_load_only() {
    let started = Instant::now();
    let sequences = read_fasta(Path::new(&env("GPS_SSR_FASTA")));
    println!(
        "load_seconds={:.3} sequences={} bases={}",
        started.elapsed().as_secs_f64(),
        sequences.len(),
        total_bases(&sequences)
    );
}

/// Punto 2: cargar y barrer con los criterios por omisión.
#[test]
#[ignore = "medición: necesita GPS_SSR_FASTA"]
fn measure_scan() {
    let started = Instant::now();
    let sequences = read_fasta(Path::new(&env("GPS_SSR_FASTA")));
    let loaded = started.elapsed().as_secs_f64();
    let scan = Instant::now();
    let (candidates, counts) = detect_all(&sequences, &DetectionCriteria::default());
    println!(
        "load_seconds={loaded:.3} scan_seconds={:.3} detected={} accepted={}",
        scan.elapsed().as_secs_f64(),
        counts.detected_total(),
        candidates.len()
    );
    assert!(counts.is_balanced());
}

/// Punto 1: recuentos por etapa y longitud de motivo con tres suelos, y las
/// cifras auxiliares de cada decisión de la pausa.
#[test]
#[ignore = "medición: necesita GPS_SSR_FASTA"]
fn measure_detection_stages() {
    let sequences = read_fasta(Path::new(&env("GPS_SSR_FASTA")));
    let genome = total_bases(&sequences);
    println!("sequences={} bases={genome}", sequences.len());
    let floors: [(&str, [u64; 6]); 3] = [
        ("bajo", [6, 4, 4, 4, 3, 3]),
        ("omision", [8, 6, 6, 6, 5, 4]),
        ("alto", [10, 8, 8, 8, 7, 6]),
    ];
    for (label, min_copies) in floors {
        let criteria = DetectionCriteria {
            min_copies,
            ..DetectionCriteria::default()
        };
        let started = Instant::now();
        let (candidates, counts) = detect_all(&sequences, &criteria);
        let covered: u64 = candidates.iter().map(|(_, c)| c.length()).sum();
        println!(
            "\n## suelo {label} {min_copies:?} ({:.2} s)\n{}bases_en_candidatos={covered} fraccion={:.5}",
            started.elapsed().as_secs_f64(),
            stage_table(&counts),
            covered as f64 / genome as f64
        );
        assert!(counts.is_balanced());
    }

    let criteria = DetectionCriteria::default();
    let (candidates, _) = detect_all(&sequences, &criteria);

    // DP-SSR-02: cuántos intervalos pasarían el suelo contando copias por
    // longitud y no lo pasan por copias completas; y de ellos, cuántos
    // superarían además pureza y tirada.
    let periods = PeriodRange::new(DEFAULT_PERIODS.0, DEFAULT_PERIODS.1).expect("1..=6");
    let mut by_length_only = [0_u64; 6];
    let mut by_length_only_pure = [0_u64; 6];
    for (_, bases) in &sequences {
        for interval in find_tandem_repeats(bases, periods, &ScanParams::default()) {
            let period = usize::from(interval.period);
            let tract = &bases[interval.start as usize..interval.end as usize];
            let Ok(motif) = RepeatMotif::new(&tract[..period]) else {
                continue;
            };
            let floor = criteria.min_copies[period - 1];
            let observed = observe_tract(tract, std::slice::from_ref(&motif));
            let by_length = tract.len() as u64 / period as u64;
            let complete = observed.covered_bases() / period as u64;
            if by_length >= floor && complete < floor {
                by_length_only[period - 1] += 1;
                if observed.purity_as_f64() >= criteria.min_purity
                    && observed.longest_run_units() >= criteria.min_longest_run_units
                {
                    by_length_only_pure[period - 1] += 1;
                }
            }
        }
    }
    println!(
        "\n## DP-SSR-02\npasan_por_longitud_y_no_por_copias_completas={by_length_only:?} total={}\n\
         de_ellos_con_pureza_y_tirada={by_length_only_pure:?} total={}",
        by_length_only.iter().sum::<u64>(),
        by_length_only_pure.iter().sum::<u64>()
    );

    // DP-SSR-04: copias parciales al final; DP-SSR-08: N dentro del tracto.
    let mut partial = [0_u64; 6];
    let mut internal_n = 0_u64;
    let mut touching = 0_u64;
    let mut purity_sum = 0.0_f64;
    for (position, (sequence, candidate)) in candidates.iter().enumerate() {
        let period = candidate.motif().length() as usize;
        if candidate.length() % period as u64 != 0 {
            partial[period - 1] += 1;
        }
        let bases = &sequences[*sequence].1[candidate.start() as usize..candidate.end() as usize];
        if bases.contains(&b'N') {
            internal_n += 1;
        }
        purity_sum += candidate.covered_bases() as f64 / candidate.length() as f64;
        if let Some((next_sequence, next)) = candidates.get(position + 1)
            && next_sequence == sequence
            && next.start() == candidate.end()
        {
            touching += 1;
        }
    }
    let lengths: Vec<u64> = candidates.iter().map(|(_, c)| c.length()).collect();
    let mut sorted = lengths.clone();
    sorted.sort_unstable();
    let percentile = |q: f64| sorted[((sorted.len() - 1) as f64 * q) as usize];
    println!(
        "\n## DP-SSR-04 y 08\ncon_copia_parcial_final={partial:?} total={}\ncon_N_interna={internal_n}\n\
         pares_que_se_tocan={touching}\npureza_media={:.4}\n\
         longitud_pb_p50={} p90={} p99={} max={}",
        partial.iter().sum::<u64>(),
        purity_sum / candidates.len() as f64,
        percentile(0.5),
        percentile(0.9),
        percentile(0.99),
        sorted[sorted.len() - 1]
    );

    // DP-SSR-07: unidades de los que superan todo menos el tope.
    let uncapped = DetectionCriteria {
        max_units: u64::MAX,
        ..DetectionCriteria::default()
    };
    let (all, _) = detect_all(&sequences, &uncapped);
    let mut above = [0_u64; 5];
    let thresholds = [50_u64, 100, 200, 500, 1000];
    let mut largest = 0_u64;
    for (_, candidate) in &all {
        let units = candidate.length() / candidate.motif().length();
        largest = largest.max(units);
        for (slot, threshold) in thresholds.iter().enumerate() {
            if candidate.length() > threshold * candidate.motif().length() {
                above[slot] += 1;
            }
        }
    }
    println!(
        "\n## DP-SSR-07\nsin_tope={} mas_de_50u={} mas_de_100u={} mas_de_200u={} mas_de_500u={} \
         mas_de_1000u={} maximo_unidades={largest}",
        all.len(),
        above[0],
        above[1],
        above[2],
        above[3],
        above[4]
    );

    // DP-SSR-18: densidad por secuencia.
    println!("\n## Por secuencia\nsecuencia\tbases\tcandidatos\tpor_Mb");
    for (index, (name, bases)) in sequences.iter().enumerate() {
        let count = candidates.iter().filter(|(s, _)| *s == index).count();
        println!(
            "{name}\t{}\t{count}\t{:.1}",
            bases.len(),
            count as f64 / (bases.len() as f64 / 1e6)
        );
    }
}

/// Punto 3, preparación: escribe el documento general con `GPS_SSR_LOCI`
/// candidatos repartidos por el genoma (uno de cada `k`, en orden de
/// referencia), una matriz de SNV de una fila fuera de todo tracto y, si se
/// pide, un evento técnico con tartamudeo sobre los diez primeros individuos.
#[test]
#[ignore = "medición: necesita GPS_SSR_FASTA, GPS_SSR_DIR, GPS_SSR_LOCI, GPS_SSR_INDIVIDUALS y GPS_SSR_STUTTER"]
fn prepare_publication() {
    let fasta = PathBuf::from(env("GPS_SSR_FASTA"));
    let fasta = fs::canonicalize(&fasta).expect("ruta del FASTA");
    let dir = PathBuf::from(env("GPS_SSR_DIR"));
    let loci = env_u64("GPS_SSR_LOCI") as usize;
    let individuals = env_u64("GPS_SSR_INDIVIDUALS") as usize;
    let stutter = env("GPS_SSR_STUTTER") == "1";
    fs::create_dir_all(&dir).expect("directorio");

    let sequences = read_fasta(&fasta);
    let (candidates, _) = detect_all(&sequences, &DetectionCriteria::default());
    assert!(
        loci <= candidates.len(),
        "solo hay {} candidatos",
        candidates.len()
    );
    let step = candidates.len() / loci;
    let chosen: Vec<&(usize, RepeatCandidate)> =
        candidates.iter().step_by(step).take(loci).collect();

    // Una fila de SNV en la primera secuencia, lejos de todo tracto elegido.
    let (snv_sequence, snv_bases) = &sequences[0];
    let occupied: Vec<(u64, u64)> = chosen
        .iter()
        .filter(|(sequence, _)| *sequence == 0)
        .map(|(_, c)| (c.start().saturating_sub(1), c.end() + 1))
        .collect();
    let position = (5_000..snv_bases.len() as u64)
        .find(|&p| {
            matches!(snv_bases[p as usize], b'A' | b'C' | b'G' | b'T')
                && !occupied.iter().any(|&(s, e)| s <= p && p < e)
        })
        .expect("posición libre");
    let reference = snv_bases[position as usize];
    let alternate = match reference {
        b'A' => 'C',
        b'C' => 'G',
        b'G' => 'T',
        _ => 'A',
    };
    fs::write(
        dir.join("population-snvs.tsv"),
        format!(
            "contig\tposition\tref\talt\tpopulation_id\talt_af\n{snv_sequence}\t{}\t{}\t{alternate}\tP1\t0.3\n",
            position + 1,
            char::from(reference)
        ),
    )
    .expect("matriz");

    let mut toml = format!(
        "config_kind = \"general\"\nschema_version = \"4\"\nseed = \"{SEED}\"\n\n\
         [population_snvs]\nmode = \"supplied\"\ninput_path = \"population-snvs.tsv\"\n\n\
         [reference]\nkind = \"supplied\"\nid = \"measured\"\nfasta = \"{}\"\n\n",
        fasta.display()
    );
    for (name, bases) in &sequences {
        let class = if name.ends_with("ch00") {
            "scaffold"
        } else {
            "chromosome"
        };
        let _ = write!(
            toml,
            "[[reference.sequences]]\nid = \"{name}\"\nclass = \"{class}\"\nlength = {}\n\n",
            bases.len()
        );
    }
    let mut tract_bases = 0_u64;
    for (index, (sequence, candidate)) in chosen.iter().enumerate() {
        tract_bases += candidate.length();
        let _ = write!(
            toml,
            "[[reference.repeat_loci]]\nid = \"SSR{index:06}\"\nsequence = \"{}\"\nstart = {}\n\
             end = {}\nmotif = \"{}\"\n\n[reference.repeat_loci.slippage]\nmodel = \"stepwise\"\n\
             mean_events = 1.0\n\n",
            sequences[*sequence].0,
            candidate.start() + 1,
            candidate.end(),
            candidate.motif().as_text()
        );
    }
    let ids: Vec<String> = (1..=individuals).map(|i| format!("I{i:03}")).collect();
    let quoted: Vec<String> = ids.iter().map(|id| format!("\"{id}\"")).collect();
    let _ = write!(
        toml,
        "[[initial_populations]]\npopulation_id = \"P1\"\nindividual_count = {individuals}\n\
         individual_ids = [{}]\nploidy = 2\ngenotype_model = \"hwe\"\n\n\
         [vcf_output]\nmode = \"text\"\n\n",
        quoted.join(", ")
    );
    if stutter {
        let event_ids = &quoted[..quoted.len().min(10)];
        let _ = write!(
            toml,
            "[technical_event]\nevent_id = \"reads\"\nindividual_ids = [{}]\n\
             target_coverage_wgs = 0.1\n\n[technical_event.stutter]\nenabled = true\n\n\
             [[technical_event.stutter.steps]]\nunits = -3\nweight = 0.05\n\n\
             [[technical_event.stutter.steps]]\nunits = -2\nweight = 0.15\n\n\
             [[technical_event.stutter.steps]]\nunits = -1\nweight = 0.6\n\n\
             [[technical_event.stutter.steps]]\nunits = 1\nweight = 0.2\n\n",
            event_ids.join(", ")
        );
        let barcode = |index: usize, salt: usize| -> String {
            (0..8)
                .map(|digit| char::from(b"ACGT"[((index + salt) >> (2 * digit)) & 3]))
                .collect()
        };
        for (index, id) in ids.iter().take(10).enumerate() {
            let _ = write!(
                toml,
                "[[technical_event.libraries]]\nindividual_id = \"{id}\"\ni7 = \"{}\"\ni5 = \"{}\"\n\n",
                barcode(index, 0),
                barcode(index, 1000)
            );
        }
    }
    fs::write(dir.join("general.toml"), &toml).expect("documento general");
    println!(
        "candidates={} step={step} declared_loci={} tract_bases={tract_bases} toml_bytes={}",
        candidates.len(),
        chosen.len(),
        toml.len()
    );
}

/// Punto 3 y 4: publica el documento de `GPS_SSR_TOML` por la biblioteca y
/// mide cada fase. `execute` repite el preflight por dentro, así que lo
/// propio de la publicación es `execute - preflight`.
#[test]
#[ignore = "medición: necesita GPS_SSR_TOML y GPS_SSR_DEST"]
fn measure_publication() {
    let path = PathBuf::from(env("GPS_SSR_TOML"));
    let destination = PathBuf::from(env("GPS_SSR_DEST"));
    let bytes = fs::read(&path).expect("documento general");
    let started = Instant::now();
    let admitted = admit_toml_header(&path, bytes, ExpectedConfigKind::General)
        .unwrap_or_else(|diagnostic| panic!("{diagnostic:?}"));
    let admit = started.elapsed().as_secs_f64();
    let declared =
        parse_declared_config(admitted).unwrap_or_else(|diagnostic| panic!("{diagnostic:?}"));
    let parse = started.elapsed().as_secs_f64() - admit;
    let build = Instant::now();
    let effective =
        build_effective_config(&declared).unwrap_or_else(|diagnostic| panic!("{diagnostic:?}"));
    let build = build.elapsed().as_secs_f64();
    let preflight = Instant::now();
    preflight_validated_config_with_authorizations(&effective, ExecutionAuthorizations::default())
        .unwrap_or_else(|diagnostic| panic!("{diagnostic:?}"));
    let preflight = preflight.elapsed().as_secs_f64();
    let execute = Instant::now();
    execute_validated_config(&effective, &destination)
        .unwrap_or_else(|diagnostic| panic!("{diagnostic:?}"));
    let execute = execute.elapsed().as_secs_f64();
    println!(
        "admit_s={admit:.3} parse_s={parse:.3} build_s={build:.3} preflight_s={preflight:.3} \
         execute_s={execute:.3} total_s={:.3}",
        started.elapsed().as_secs_f64()
    );
    let mut files = Vec::new();
    collect_files(&destination, &destination, &mut files);
    files.sort();
    for (relative, size) in files {
        println!("file\t{size}\t{relative}");
    }
}

fn collect_files(root: &Path, directory: &Path, out: &mut Vec<(String, u64)>) {
    for entry in fs::read_dir(directory).expect("directorio publicado") {
        let entry = entry.expect("entrada");
        let path = entry.path();
        if path.is_dir() {
            collect_files(root, &path, out);
        } else {
            let size = entry.metadata().expect("metadatos").len();
            let relative = path
                .strip_prefix(root)
                .expect("dentro")
                .display()
                .to_string();
            out.push((relative, size));
        }
    }
}
