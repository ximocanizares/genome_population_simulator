//! F1.7.1: generador de referencia simulada **con regiones declaradas**.
//!
//! El orden que la decisión S-03 fija —secuencias y tamaño, después las
//! regiones con su secuencia exacta, después el relleno al azar— produce
//! exactamente el **mismo FASTA** que un generador que rellenase primero y
//! sobrescribiese después. El único observable que los separa es el **consumo
//! de aleatoriedad**, y por eso la prueba del contador no es un adorno: es la
//! única que discrimina. Las comprobaciones de secuencia y de rechazo están
//! aquí porque hacen falta, no porque distingan el orden.

use std::fs;
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use genome_population_simulator::randomness::{RandomWordSource, RootSeed};
use genome_population_simulator::reference::SequenceClass;
use genome_population_simulator::reference_fasta::write_generated_fasta;
use genome_population_simulator::reference_generation::{
    BaseFrequency, DeclaredRegion, GeneratedSequenceSpec, generate_reference_from_source,
    generate_reference_with_regions,
};

const ROOT: RootSeed = RootSeed::from_bytes([
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
    26, 27, 28, 29, 30, 31,
]);

static NONCE: AtomicU64 = AtomicU64::new(0);

/// Fuente de palabras que **cuenta** cuántas se le piden.
struct CountingWords {
    inner: u64,
    consumed: u64,
}

impl CountingWords {
    const fn new() -> Self {
        Self {
            inner: 0x2545_F491_4F6C_DD1D,
            consumed: 0,
        }
    }
}

impl RandomWordSource for CountingWords {
    fn next_word(&mut self) -> u64 {
        self.consumed += 1;
        // Generador propio de la prueba: no imita a ninguno normativo, solo
        // tiene que ser determinista para que el recuento sea reproducible.
        self.inner ^= self.inner << 13;
        self.inner ^= self.inner >> 7;
        self.inner ^= self.inner << 17;
        self.inner
    }
}

fn composition(a: f64, c: f64, g: f64, t: f64) -> [BaseFrequency; 4] {
    [
        BaseFrequency {
            base: b'A',
            probability: a,
        },
        BaseFrequency {
            base: b'C',
            probability: c,
        },
        BaseFrequency {
            base: b'G',
            probability: g,
        },
        BaseFrequency {
            base: b'T',
            probability: t,
        },
    ]
}

fn specs(length: u64) -> [GeneratedSequenceSpec<'static>; 1] {
    [GeneratedSequenceSpec {
        id: "chr1",
        class: SequenceClass::Chromosome,
        length,
    }]
}

/// **La prueba que discrimina el orden.** El relleno consume exactamente una
/// palabra por base **no declarada**: la longitud de la secuencia menos las
/// bases de las regiones. Un generador que rellenase primero y sobrescribiese
/// después consumiría una por base de la secuencia entera y produciría el
/// **mismo FASTA**, de modo que ninguna comprobación sobre las bases lo caza.
#[test]
fn the_random_fill_consumes_one_word_per_undeclared_base_only() {
    let tract = b"CAGCAGCAGCAGCAG";
    let regions = [DeclaredRegion {
        sequence_id: "chr1",
        start: 40,
        bases: tract,
    }];
    let mut words = CountingWords::new();
    let reference = generate_reference_from_source(
        "reference-regions",
        &specs(200),
        &composition(0.25, 0.25, 0.25, 0.25),
        0.0,
        &regions,
        &mut words,
    )
    .expect("referencia con región declarada");
    assert_eq!(reference.total_length(), 200);
    assert_eq!(
        words.consumed,
        200 - tract.len() as u64,
        "el relleno no puede consumir aleatoriedad para las bases declaradas"
    );

    let mut without = CountingWords::new();
    generate_reference_from_source(
        "reference-regions",
        &specs(200),
        &composition(0.25, 0.25, 0.25, 0.25),
        0.0,
        &[],
        &mut without,
    )
    .expect("referencia sin regiones");
    assert_eq!(
        without.consumed, 200,
        "sin regiones el consumo sigue siendo una palabra por base"
    );
}

/// Una región declarada aparece **exactamente** en su coordenada, localizada
/// con un `grep -ob` externo al productor.
#[test]
fn a_declared_region_lands_at_its_coordinate_according_to_external_grep() {
    let tract = b"CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG";
    let regions = [DeclaredRegion {
        sequence_id: "chr1",
        start: 100,
        bases: tract,
    }];
    // Composición degenerada: todo lo que no sea la región es una `A`, de modo
    // que el `grep` externo no puede acertar por casualidad en otro sitio.
    let reference = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &specs(400),
        &composition(1.0, 0.0, 0.0, 0.0),
        0.0,
        &regions,
    )
    .expect("referencia con región declarada");

    let bases = reference.sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases materializadas");
    assert_eq!(&bases[100..130], tract);

    let path = temporary("fasta");
    let mut file = fs::File::create(&path).expect("fichero FASTA");
    let mut bytes = Vec::new();
    write_generated_fasta(&reference, &mut bytes).expect("FASTA canónico");
    file.write_all(&bytes).expect("escritura");
    drop(file);

    let output = Command::new("grep")
        .arg("-ob")
        .arg(String::from_utf8_lossy(tract).as_ref())
        .arg(&path)
        .output()
        .expect("grep externo");
    let found = String::from_utf8(output.stdout).expect("salida de grep");
    let offset: usize = found
        .split(':')
        .next()
        .expect("desplazamiento de grep")
        .trim()
        .parse()
        .expect("desplazamiento numérico");
    // Cabecera `>chr1\n` son seis bytes; después, una línea cada ochenta bases.
    assert_eq!(
        offset,
        ">chr1\n".len() + 100 + 100 / 80,
        "la región tiene que caer en su coordenada, ni una base más allá"
    );
    fs::remove_file(&path).ok();
}

/// El relleno al azar respeta las fracciones configuradas **fuera** de las
/// regiones declaradas, y no las aplica dentro.
#[test]
fn the_random_fill_honours_the_composition_outside_the_declared_region() {
    let tract = b"CCCCCCCCCC";
    let regions = [DeclaredRegion {
        sequence_id: "chr1",
        start: 10,
        bases: tract,
    }];
    let reference = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &specs(60),
        &composition(1.0, 0.0, 0.0, 0.0),
        0.0,
        &regions,
    )
    .expect("referencia con región declarada");
    let bases = reference.sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases materializadas");
    assert!(
        bases[..10].iter().all(|base| *base == b'A'),
        "fuera de la región manda la composición"
    );
    assert_eq!(&bases[10..20], tract, "dentro manda la secuencia declarada");
    assert!(bases[20..].iter().all(|base| *base == b'A'));
}

/// Misma semilla, mismo FASTA byte a byte, también con regiones declaradas.
#[test]
fn the_same_seed_yields_the_same_fasta_byte_for_byte_with_regions() {
    let regions = [DeclaredRegion {
        sequence_id: "chr1",
        start: 5,
        bases: b"GATAGATAGATA",
    }];
    let mut first = Vec::new();
    let mut second = Vec::new();
    for sink in [&mut first, &mut second] {
        let reference = generate_reference_with_regions(
            "reference-regions",
            ROOT,
            &specs(120),
            &composition(0.3, 0.2, 0.2, 0.3),
            0.0,
            &regions,
        )
        .expect("referencia con región declarada");
        write_generated_fasta(&reference, sink).expect("FASTA canónico");
    }
    assert_eq!(first, second);
}

/// Dos regiones declaradas que se solapan se rechazan **antes** de generar,
/// nombrando los dos intervalos.
#[test]
fn two_overlapping_declared_regions_are_rejected_naming_both_intervals() {
    let regions = [
        DeclaredRegion {
            sequence_id: "chr1",
            start: 10,
            bases: b"CAGCAGCAG",
        },
        DeclaredRegion {
            sequence_id: "chr1",
            start: 15,
            bases: b"GATAGATA",
        },
    ];
    let error = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &specs(100),
        &composition(0.25, 0.25, 0.25, 0.25),
        0.0,
        &regions,
    )
    .expect_err("dos regiones solapadas");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("overlap"), "{rendered}");
    assert!(rendered.contains("[10, 19)"), "{rendered}");
    assert!(rendered.contains("[15, 23)"), "{rendered}");
}

/// Una región que rebasa el final de la secuencia se rechaza.
#[test]
fn a_region_that_runs_past_the_end_of_the_sequence_is_rejected() {
    let regions = [DeclaredRegion {
        sequence_id: "chr1",
        start: 95,
        bases: b"CAGCAGCAG",
    }];
    let error = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &specs(100),
        &composition(0.25, 0.25, 0.25, 0.25),
        0.0,
        &regions,
    )
    .expect_err("región fuera de la secuencia");
    let rendered = format!("{error:?}");
    assert!(rendered.contains("runs past"), "{rendered}");
}

/// Una región sobre una secuencia que no existe se rechaza nombrándola.
#[test]
fn a_region_on_an_unknown_sequence_is_rejected() {
    let regions = [DeclaredRegion {
        sequence_id: "chr9",
        start: 0,
        bases: b"CAG",
    }];
    let error = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &specs(100),
        &composition(0.25, 0.25, 0.25, 0.25),
        0.0,
        &regions,
    )
    .expect_err("secuencia desconocida");
    assert!(format!("{error:?}").contains("chr9"));
}

/// Sin regiones declaradas, el FASTA es **el vector congelado**, byte a byte.
///
/// *Corregido tras la revisión adversarial.* La versión anterior comparaba
/// `generate_reference` contra `generate_reference_with_regions(..., &[])`, y eso
/// es un **oráculo circular**: la primera acaba llamando a la segunda con la
/// lista vacía, de modo que la prueba comparaba una llamada consigo misma y **no
/// podía fallar**. Un mutante que cambiase los bytes del relleno sin cambiar
/// cuántas palabras consume la dejaba pasar.
///
/// Ahora se congela la huella del FASTA producido. Es un valor externo al
/// productor: para moverlo hay que decidir moverlo.
#[test]
fn without_declared_regions_the_output_matches_the_frozen_digest() {
    use genome_population_simulator::reference_generation::generate_reference;
    // Vector congelado. **Comprobado contra el generador de HEAD sin tocar**, con
    // una copia limpia del repositorio compilada aparte: es el FASTA de antes del
    // cambio, no el de después. Si alguna vez hay que moverlo, hay que decidir
    // moverlo y decir por qué.
    const FROZEN_DIGEST: &str = "48f9c4793e52de643bfee894326d495be9ab82771fea5d731409ec9e6f0df372";

    let with_empty_regions = generate_reference_with_regions(
        "reference-main",
        ROOT,
        &specs(300),
        &composition(0.3, 0.2, 0.2, 0.3),
        0.0,
        &[],
    )
    .expect("referencia sin regiones");
    let mut bytes = Vec::new();
    write_generated_fasta(&with_empty_regions, &mut bytes).expect("FASTA");
    assert_eq!(
        digest_of(&bytes),
        FROZEN_DIGEST,
        "el relleno sin regiones no puede cambiar ni un byte"
    );

    // Y además el camino histórico sigue dando lo mismo. Esta comparación **no
    // discrimina por sí sola** —las dos llamadas comparten implementación—; solo
    // comprueba que la fachada anterior no se ha desviado.
    let legacy = generate_reference(
        "reference-main",
        ROOT,
        &specs(300),
        &composition(0.3, 0.2, 0.2, 0.3),
    )
    .expect("referencia por el camino anterior");
    let mut legacy_bytes = Vec::new();
    write_generated_fasta(&legacy, &mut legacy_bytes).expect("FASTA");
    assert_eq!(legacy_bytes, bytes);
}

/// **Las regiones se filtran por secuencia.** Sin ese filtro, una región
/// declarada sobre el segundo cromosoma se escribiría también dentro del
/// primero: corrupción silenciosa, porque el FASTA seguiría siendo válido.
///
/// Ningún otro fixture de este fichero usa más de una secuencia, de modo que
/// esta es la única prueba que ejercita el filtro.
#[test]
fn a_region_declared_on_one_sequence_does_not_leak_into_another() {
    let two = [
        GeneratedSequenceSpec {
            id: "chr1",
            class: SequenceClass::Chromosome,
            length: 60,
        },
        GeneratedSequenceSpec {
            id: "chr2",
            class: SequenceClass::Chromosome,
            length: 60,
        },
    ];
    let regions = [DeclaredRegion {
        sequence_id: "chr2",
        start: 10,
        bases: b"GGGGGGGGGG",
    }];
    let reference = generate_reference_with_regions(
        "reference-regions",
        ROOT,
        &two,
        &composition(1.0, 0.0, 0.0, 0.0),
        0.0,
        &regions,
    )
    .expect("dos secuencias con una región en la segunda");

    let first = reference.sequences()[0]
        .sequence()
        .as_bytes()
        .expect("bases del primero");
    let second = reference.sequences()[1]
        .sequence()
        .as_bytes()
        .expect("bases del segundo");
    assert!(
        first.iter().all(|base| *base == b'A'),
        "el primer cromosoma no puede contener ni una base de la región del segundo"
    );
    assert_eq!(&second[10..20], b"GGGGGGGGGG");
    assert!(second[..10].iter().all(|base| *base == b'A'));
}

fn digest_of(bytes: &[u8]) -> String {
    let path = temporary("fa");
    fs::write(&path, bytes).expect("volcado temporal");
    let identity = genome_population_simulator::integrity::identify_file(&path, "generated.fa")
        .expect("huella");
    fs::remove_file(&path).ok();
    identity.digest().to_owned()
}

fn temporary(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "gps-f1-7-1-{}-{}-{}.{suffix}",
        std::process::id(),
        NONCE.fetch_add(1, Ordering::Relaxed),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .expect("clock")
            .as_nanos()
    ))
}
