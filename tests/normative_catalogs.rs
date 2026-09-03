use std::collections::{BTreeMap, BTreeSet};
use std::io::Write;
use std::process::{Command, Stdio};

const DIAGNOSTIC_CATALOG_V1: &str = include_str!("../catalogs/diagnostics-v1.tsv");
const DIAGNOSTIC_CATALOG_V2: &str = include_str!("../catalogs/diagnostics-v2.tsv");
const DIAGNOSTIC_CATALOG: &str = include_str!("../catalogs/diagnostics-v3.tsv");
const DIAGNOSTIC_CATALOG_V4: &str = include_str!("../catalogs/diagnostics-v4.tsv");
const DIAGNOSTIC_SCHEMA: &str = include_str!("../schemas/diagnostics-v1.schema.json");
const SUBSEED_CATALOG_V1: &str = include_str!("../catalogs/subseed-scopes-v1.tsv");
const SUBSEED_CATALOG_V2: &str = include_str!("../catalogs/subseed-scopes-v2.tsv");
const SUBSEED_CATALOG: &str = include_str!("../catalogs/subseed-scopes-v3.tsv");
const ID_CATALOG: &str = include_str!("../catalogs/id-classes-v1.tsv");
const ID_CATALOG_V2: &str = include_str!("../catalogs/id-classes-v2.tsv");
const ID_CATALOG_V3: &str = include_str!("../catalogs/id-classes-v3.tsv");
const ID_CATALOG_V4: &str = include_str!("../catalogs/id-catalog-v4.tsv");
const F5_RNG_CATALOG: &str = include_str!("../catalogs/f5-rng-catalog-v1.tsv");
const VARIANT_ID_CATALOG: &str = include_str!("../catalogs/variant-id-v1.tsv");
const PARAMETER_CATALOG: &str = include_str!("../catalogs/parameter-families-v1.tsv");
const SUBSEED_CONTEXT: &str = "genome_population_simulator/v1/subseed";
const ID_CONTEXT: &str = "genome_population_simulator/v1/id";

// These fingerprints are the review lock for the published contracts. A
// deliberate content change must use a new versioned artifact and add its new
// fingerprint instead of replacing a historical value here.
const REVIEWED_CONTRACT_FINGERPRINTS: [(&str, &str, &str); 14] = [
    (
        "catalogs/diagnostics-v1.tsv",
        DIAGNOSTIC_CATALOG_V1,
        "4fb2be15e392c8fa9fb8791a5d8821fe14da901f88286711c2fb29d3d0a6997f",
    ),
    (
        "catalogs/diagnostics-v2.tsv",
        DIAGNOSTIC_CATALOG_V2,
        "24e6153b152e78236dfc2a3ae84e17f4c7ab733045b934d0fbf715ebcbdab432",
    ),
    (
        "catalogs/diagnostics-v3.tsv",
        DIAGNOSTIC_CATALOG,
        "6e8a8fa3c8c1306ed360bad6de402649020d109a10fc0c5b4072d47df66f64f9",
    ),
    (
        "catalogs/diagnostics-v4.tsv",
        DIAGNOSTIC_CATALOG_V4,
        "cfe0c0c1eae5f77603729bff365ae6ea070410dab32dce4e5bbb2c8276f6a473",
    ),
    (
        "schemas/diagnostics-v1.schema.json",
        DIAGNOSTIC_SCHEMA,
        "ab2b9f2cd41aa1fb28e7de3b0e72b5e8132031bd83b1adc72b09ee8b49d1855a",
    ),
    (
        "catalogs/subseed-scopes-v1.tsv",
        SUBSEED_CATALOG_V1,
        "b6c56848b1ae48d5646dcdc515709259ea65fc5c8b6a407bf403e4da8477675d",
    ),
    (
        "catalogs/subseed-scopes-v2.tsv",
        SUBSEED_CATALOG_V2,
        "c40f812d1645bb11a2e2641cc76a5b6db583f0505e246cb41a24bc0aa554ceff",
    ),
    (
        "catalogs/subseed-scopes-v3.tsv",
        SUBSEED_CATALOG,
        "415ead5f7083220953f3795f51da1cc464c41783174d11359d3e0e5759460066",
    ),
    (
        "catalogs/id-classes-v1.tsv",
        ID_CATALOG,
        "7eca3d96b64ed2c5944f06c95b5b547464ace88ca6f5b3d3479235b56d7f23cf",
    ),
    (
        "catalogs/id-classes-v2.tsv",
        ID_CATALOG_V2,
        "386267beeecd7e25638e0b5b677c8c10cca871847fe0d4566f964e042f074b57",
    ),
    (
        "catalogs/id-classes-v3.tsv",
        ID_CATALOG_V3,
        "34dafd1d47a0f96a9e6030270ac2f7c84a0bf050ffc6926d7b919bdef8abf390",
    ),
    (
        "catalogs/id-catalog-v4.tsv",
        ID_CATALOG_V4,
        "218f8fe4f6a68ac4c4da9907c8aad5c017ceb4be04dee1198b2d8b7fa9410d3b",
    ),
    (
        "catalogs/f5-rng-catalog-v1.tsv",
        F5_RNG_CATALOG,
        "5c2a1b68e744cc333c3b28a69263202b9b2cd5acdb4910db54cafcf6ff2fbc81",
    ),
    (
        "catalogs/parameter-families-v1.tsv",
        PARAMETER_CATALOG,
        "5507a763b91520449a022bd6e371cb4f03cafa5547a68ea8c0e1cfc3221b9a7a",
    ),
];

#[derive(Debug, PartialEq, Eq)]
struct CatalogError {
    line: usize,
    reason: &'static str,
}

fn sha256sum(input: &str) -> String {
    let mut child = Command::new("sha256sum")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .expect("the independently approved sha256sum validator must be available");
    child
        .stdin
        .take()
        .expect("sha256sum stdin must be piped")
        .write_all(input.as_bytes())
        .expect("catalog bytes must be writable to sha256sum");
    let output = child
        .wait_with_output()
        .expect("sha256sum must finish normally");
    assert!(output.status.success(), "sha256sum rejected catalog bytes");
    String::from_utf8(output.stdout)
        .expect("sha256sum output must be UTF-8")
        .split_whitespace()
        .next()
        .expect("sha256sum must emit a digest")
        .to_owned()
}

#[test]
fn published_contracts_match_their_reviewed_content() {
    assert!(DIAGNOSTIC_CATALOG_V2.starts_with("# diagnostic_catalog_version=2\n"));
    assert!(DIAGNOSTIC_CATALOG.starts_with("# diagnostic_catalog_version=3\n"));
    assert!(SUBSEED_CATALOG_V1.starts_with(&format!(
        "# derivation_version=1\n# context={SUBSEED_CONTEXT}\n"
    )));
    assert!(SUBSEED_CATALOG.starts_with(&format!(
        "# derivation_version=1\n# catalog_version=3\n# context={SUBSEED_CONTEXT}\n"
    )));
    assert!(ID_CATALOG.starts_with(&format!("# derivation_version=1\n# context={ID_CONTEXT}\n")));
    assert!(ID_CATALOG_V2.starts_with(&format!(
        "# derivation_version=1\n# catalog_version=2\n# context={ID_CONTEXT}\n"
    )));
    assert!(ID_CATALOG_V3.starts_with(&format!(
        "# catalog_version=3\n# context={ID_CONTEXT}\n# derivation_versions=v1,v2\n"
    )));
    assert!(ID_CATALOG_V4.starts_with("# catalog_version=4\n# catalog_id=id-catalog-v4\n"));
    assert!(F5_RNG_CATALOG.starts_with("# catalog_id=f5-rng-catalog-v1\n"));
    assert!(PARAMETER_CATALOG.starts_with("# parameter_catalog_version=1\n"));

    for (path, content, expected) in REVIEWED_CONTRACT_FINGERPRINTS {
        assert_eq!(
            sha256sum(content),
            expected,
            "{path} changed without introducing a new versioned contract"
        );
    }
}

#[test]
fn variant_id_catalog_has_an_independent_review_fingerprint() {
    assert_eq!(
        sha256sum(VARIANT_ID_CATALOG),
        "c1cc71ff9b19e5f00f62e260f49f57107098c3afc4e431c01938395263d5aff1"
    );
    let fields = VARIANT_ID_CATALOG
        .lines()
        .last()
        .expect("vector")
        .split('|')
        .collect::<Vec<_>>();
    assert_eq!(fields.len(), 5);
    assert_eq!(fields[1], "VAR_");
    assert!(
        fields[2].contains("reference_id=")
            && fields[2].contains("sequence_id=")
            && fields[2].contains("position_1based=")
    );
}

fn validate_diagnostic_catalog(catalog: &str) -> Result<usize, CatalogError> {
    const CATEGORIES: [&str; 8] = [
        "CONF", "INPUT", "FMT", "DOMAIN", "RES", "TOOL", "IO", "INTERNAL",
    ];
    const EXIT_CODES: [&str; 7] = ["0", "1", "2", "3", "4", "5", "70"];

    let mut codes = BTreeSet::new();
    let mut count = 0;
    for (index, line) in catalog.lines().enumerate() {
        let line_number = index + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split('|').collect();
        if fields.len() != 6 {
            return Err(CatalogError {
                line: line_number,
                reason: "field count",
            });
        }
        let [code, category, meaning, context, exit_code, planned_steps] = fields.as_slice() else {
            return Err(CatalogError {
                line: line_number,
                reason: "field count",
            });
        };
        if !CATEGORIES.contains(category) {
            return Err(CatalogError {
                line: line_number,
                reason: "unknown category",
            });
        }
        let expected_prefix = format!("GPS-{category}-");
        if !code.starts_with(&expected_prefix)
            || code.len() != expected_prefix.len() + 3
            || !code[expected_prefix.len()..]
                .bytes()
                .all(|byte| byte.is_ascii_digit())
        {
            return Err(CatalogError {
                line: line_number,
                reason: "invalid code",
            });
        }
        if !codes.insert(*code) {
            return Err(CatalogError {
                line: line_number,
                reason: "duplicate code",
            });
        }
        if meaning.trim().is_empty() {
            return Err(CatalogError {
                line: line_number,
                reason: "empty meaning",
            });
        }
        if context.trim().is_empty() {
            return Err(CatalogError {
                line: line_number,
                reason: "empty context",
            });
        }
        if !EXIT_CODES.contains(exit_code) {
            return Err(CatalogError {
                line: line_number,
                reason: "invalid exit code",
            });
        }
        if planned_steps.trim().is_empty() {
            return Err(CatalogError {
                line: line_number,
                reason: "empty planned steps",
            });
        }
        count += 1;
    }
    Ok(count)
}

#[test]
fn published_diagnostic_catalog_is_well_formed_and_covers_every_category() {
    assert_eq!(validate_diagnostic_catalog(DIAGNOSTIC_CATALOG_V2), Ok(46));
    assert_eq!(validate_diagnostic_catalog(DIAGNOSTIC_CATALOG), Ok(47));
    assert_eq!(validate_diagnostic_catalog(DIAGNOSTIC_CATALOG_V4), Ok(57));
    assert!(DIAGNOSTIC_CATALOG_V4.starts_with("# diagnostic_catalog_version=4\n"));
    let v3_rows = DIAGNOSTIC_CATALOG
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<BTreeSet<_>>();
    let v4_rows = DIAGNOSTIC_CATALOG_V4
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect::<BTreeSet<_>>();
    assert!(v3_rows.is_subset(&v4_rows));
    assert_eq!(v4_rows.difference(&v3_rows).count(), 10);
    assert!(
        v4_rows.difference(&v3_rows).all(|row| row.ends_with("|F5")),
        "diagnostics-v4 must add only the ten F5 rows"
    );
    for category in [
        "CONF", "INPUT", "FMT", "DOMAIN", "RES", "TOOL", "IO", "INTERNAL",
    ] {
        assert!(
            DIAGNOSTIC_CATALOG.contains(&format!("|{category}|")),
            "missing category {category}"
        );
    }
}

#[test]
fn diagnostic_catalog_rejects_duplicate_unknown_and_empty_entries() {
    let duplicate = "GPS-CONF-001|CONF|one|path|2|B3\nGPS-CONF-001|CONF|two|path|2|B3\n";
    assert_eq!(
        validate_diagnostic_catalog(duplicate),
        Err(CatalogError {
            line: 2,
            reason: "duplicate code"
        })
    );
    let unknown = "GPS-OTHER-001|OTHER|meaning|path|2|B3\n";
    assert_eq!(
        validate_diagnostic_catalog(unknown),
        Err(CatalogError {
            line: 1,
            reason: "unknown category"
        })
    );
    let empty = "GPS-CONF-001|CONF||path|2|B3\n";
    assert_eq!(
        validate_diagnostic_catalog(empty),
        Err(CatalogError {
            line: 1,
            reason: "empty meaning"
        })
    );
}

#[test]
fn human_and_json_fixtures_describe_the_same_diagnostic() {
    let human = include_str!("fixtures/diagnostics/invalid-schema-version.txt");
    let json = include_str!("fixtures/diagnostics/invalid-schema-version.json");
    for shared in [
        "GPS-CONF-005",
        "CONF",
        "La versión de esquema '2' no está admitida; se admite '1'.",
        "scenario.toml",
    ] {
        assert!(human.contains(shared));
        assert!(json.contains(shared));
    }
    assert!(json.contains("\"diagnostic_kind\": \"gps_diagnostics\""));
    assert!(json.contains("\"schema_version\": \"1\""));
}

#[test]
fn every_f1_consumer_hito_has_a_planned_diagnostic() {
    let mut by_hito = BTreeMap::new();
    for line in DIAGNOSTIC_CATALOG
        .lines()
        .filter(|line| !line.starts_with('#'))
    {
        if line.is_empty() {
            continue;
        }
        let fields: Vec<_> = line.split('|').collect();
        for step in fields[5].split(';') {
            by_hito.insert(&step[..1], true);
        }
    }
    assert_eq!(by_hito.get("B"), Some(&true));
    assert_eq!(by_hito.get("C"), Some(&true));
    assert_eq!(by_hito.get("D"), Some(&true));
    assert_eq!(by_hito.get("E"), Some(&true));
}

fn decode_hex(value: &str) -> Vec<u8> {
    assert_eq!(value.len() % 2, 0);
    value
        .as_bytes()
        .chunks_exact(2)
        .map(|pair| {
            let text = std::str::from_utf8(pair).expect("hex literal must be ASCII");
            u8::from_str_radix(text, 16).expect("hex literal must be valid")
        })
        .collect()
}

fn encode_material(root_hex: &str, components: &str) -> Vec<u8> {
    encode_material_with_version(root_hex, components, 1)
}

fn encode_material_with_version(root_hex: &str, components: &str, version: u32) -> Vec<u8> {
    let root = decode_hex(root_hex);
    assert_eq!(root.len(), 32);
    let components: Vec<_> = components
        .split(';')
        .map(|component| {
            component
                .split_once('=')
                .expect("component must contain '='")
        })
        .collect();
    let mut material = root;
    material.extend_from_slice(&version.to_be_bytes());
    material.extend_from_slice(
        &u32::try_from(components.len())
            .expect("component count must fit u32")
            .to_be_bytes(),
    );
    for (name, value) in components {
        assert!(name.is_ascii());
        material.extend_from_slice(
            &u32::try_from(name.len())
                .expect("name length must fit u32")
                .to_be_bytes(),
        );
        material.extend_from_slice(name.as_bytes());
        material.extend_from_slice(
            &u32::try_from(value.len())
                .expect("value length must fit u32")
                .to_be_bytes(),
        );
        material.extend_from_slice(value.as_bytes());
    }
    material
}

fn decode_base32_no_padding(value: &str) -> Vec<u8> {
    let mut output = Vec::new();
    let mut buffer = 0_u32;
    let mut bits = 0_u8;
    for byte in value.bytes() {
        let digit = match byte {
            b'A'..=b'Z' => u32::from(byte - b'A'),
            b'2'..=b'7' => u32::from(byte - b'2') + 26,
            _ => panic!("invalid RFC 4648 Base32 literal"),
        };
        buffer = (buffer << 5) | digit;
        bits += 5;
        if bits >= 8 {
            bits -= 8;
            output.push((buffer >> bits) as u8);
            buffer &= (1_u32 << bits) - 1;
        }
    }
    assert_eq!(buffer, 0, "unused Base32 bits must be zero");
    output
}

// Minimal test-only BLAKE3 oracle, implemented directly from the public
// compression function. It intentionally shares no crate or production code
// with the future runtime derivation implementation.
const BLAKE3_IV: [u32; 8] = [
    0x6a09_e667,
    0xbb67_ae85,
    0x3c6e_f372,
    0xa54f_f53a,
    0x510e_527f,
    0x9b05_688c,
    0x1f83_d9ab,
    0x5be0_cd19,
];
const BLAKE3_PERMUTATION: [usize; 16] = [2, 6, 3, 10, 7, 0, 4, 13, 1, 11, 12, 5, 9, 14, 15, 8];
const CHUNK_START: u32 = 1;
const CHUNK_END: u32 = 2;
const ROOT: u32 = 8;
const DERIVE_KEY_CONTEXT: u32 = 32;
const DERIVE_KEY_MATERIAL: u32 = 64;

fn blake3_mix(state: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize, x: u32, y: u32) {
    state[a] = state[a].wrapping_add(state[b]).wrapping_add(x);
    state[d] = (state[d] ^ state[a]).rotate_right(16);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] = (state[b] ^ state[c]).rotate_right(12);
    state[a] = state[a].wrapping_add(state[b]).wrapping_add(y);
    state[d] = (state[d] ^ state[a]).rotate_right(8);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] = (state[b] ^ state[c]).rotate_right(7);
}

fn blake3_round(state: &mut [u32; 16], message: &[u32; 16]) {
    blake3_mix(state, 0, 4, 8, 12, message[0], message[1]);
    blake3_mix(state, 1, 5, 9, 13, message[2], message[3]);
    blake3_mix(state, 2, 6, 10, 14, message[4], message[5]);
    blake3_mix(state, 3, 7, 11, 15, message[6], message[7]);
    blake3_mix(state, 0, 5, 10, 15, message[8], message[9]);
    blake3_mix(state, 1, 6, 11, 12, message[10], message[11]);
    blake3_mix(state, 2, 7, 8, 13, message[12], message[13]);
    blake3_mix(state, 3, 4, 9, 14, message[14], message[15]);
}

fn blake3_words(block: &[u8]) -> [u32; 16] {
    let mut padded = [0_u8; 64];
    padded[..block.len()].copy_from_slice(block);
    let mut words = [0_u32; 16];
    for (word, bytes) in words.iter_mut().zip(padded.chunks_exact(4)) {
        *word = u32::from_le_bytes(bytes.try_into().expect("word must contain four bytes"));
    }
    words
}

fn blake3_compress(
    chaining_value: [u32; 8],
    block: &[u8],
    block_len: u32,
    flags: u32,
) -> [u32; 16] {
    let mut state = [0_u32; 16];
    state[..8].copy_from_slice(&chaining_value);
    state[8..12].copy_from_slice(&BLAKE3_IV[..4]);
    state[12] = 0;
    state[13] = 0;
    state[14] = block_len;
    state[15] = flags;

    let mut message = blake3_words(block);
    for round in 0..7 {
        blake3_round(&mut state, &message);
        if round != 6 {
            message = BLAKE3_PERMUTATION.map(|index| message[index]);
        }
    }
    for index in 0..8 {
        state[index] ^= state[index + 8];
        state[index + 8] ^= chaining_value[index];
    }
    state
}

fn blake3_one_chunk(input: &[u8], key: [u32; 8], flags: u32) -> [u8; 32] {
    assert!(input.len() <= 1024, "test oracle only supports one chunk");
    let blocks: Vec<&[u8]> = if input.is_empty() {
        vec![&[]]
    } else {
        input.chunks(64).collect()
    };
    let mut chaining_value = key;
    for (index, block) in blocks.iter().enumerate() {
        let is_first = index == 0;
        let is_last = index + 1 == blocks.len();
        let block_flags = flags
            | if is_first { CHUNK_START } else { 0 }
            | if is_last { CHUNK_END | ROOT } else { 0 };
        let output = blake3_compress(
            chaining_value,
            block,
            u32::try_from(block.len()).expect("BLAKE3 block length must fit u32"),
            block_flags,
        );
        if is_last {
            let mut digest = [0_u8; 32];
            for (destination, word) in digest.chunks_exact_mut(4).zip(output[..8].iter()) {
                destination.copy_from_slice(&word.to_le_bytes());
            }
            return digest;
        }
        chaining_value.copy_from_slice(&output[..8]);
    }
    unreachable!("BLAKE3 input always has at least one block")
}

fn blake3_derive_key(context: &str, material: &[u8]) -> [u8; 32] {
    let context_key = blake3_one_chunk(context.as_bytes(), BLAKE3_IV, DERIVE_KEY_CONTEXT);
    let mut key_words = [0_u32; 8];
    for (word, bytes) in key_words.iter_mut().zip(context_key.chunks_exact(4)) {
        *word = u32::from_le_bytes(bytes.try_into().expect("key word must contain four bytes"));
    }
    blake3_one_chunk(material, key_words, DERIVE_KEY_MATERIAL)
}

fn chacha_quarter_round(state: &mut [u32; 16], a: usize, b: usize, c: usize, d: usize) {
    state[a] = state[a].wrapping_add(state[b]);
    state[d] ^= state[a];
    state[d] = state[d].rotate_left(16);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] ^= state[c];
    state[b] = state[b].rotate_left(12);
    state[a] = state[a].wrapping_add(state[b]);
    state[d] ^= state[a];
    state[d] = state[d].rotate_left(8);
    state[c] = state[c].wrapping_add(state[d]);
    state[b] ^= state[c];
    state[b] = state[b].rotate_left(7);
}

fn chacha12_first_four_u64(seed: &[u8]) -> [u64; 4] {
    assert_eq!(seed.len(), 32);
    let mut state = [
        0x6170_7865,
        0x3320_646e,
        0x7962_2d32,
        0x6b20_6574,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
        0,
    ];
    for (word, bytes) in state[4..12].iter_mut().zip(seed.chunks_exact(4)) {
        *word = u32::from_le_bytes(bytes.try_into().expect("key word must contain four bytes"));
    }
    let original = state;
    for _ in 0..6 {
        chacha_quarter_round(&mut state, 0, 4, 8, 12);
        chacha_quarter_round(&mut state, 1, 5, 9, 13);
        chacha_quarter_round(&mut state, 2, 6, 10, 14);
        chacha_quarter_round(&mut state, 3, 7, 11, 15);
        chacha_quarter_round(&mut state, 0, 5, 10, 15);
        chacha_quarter_round(&mut state, 1, 6, 11, 12);
        chacha_quarter_round(&mut state, 2, 7, 8, 13);
        chacha_quarter_round(&mut state, 3, 4, 9, 14);
    }
    for index in 0..16 {
        state[index] = state[index].wrapping_add(original[index]);
    }
    std::array::from_fn(|index| {
        let mut bytes = [0_u8; 8];
        bytes[..4].copy_from_slice(&state[index * 2].to_le_bytes());
        bytes[4..].copy_from_slice(&state[index * 2 + 1].to_le_bytes());
        u64::from_le_bytes(bytes)
    })
}

#[test]
fn subseed_vectors_publish_exact_canonical_material() {
    let rows: Vec<_> = SUBSEED_CATALOG
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(rows.len(), 6);
    let mut scopes = BTreeSet::new();
    for row in rows {
        let fields: Vec<_> = row.split('|').collect();
        assert_eq!(fields.len(), 5);
        scopes.insert(fields[0]);
        assert_eq!(encode_material(fields[2], fields[1]), decode_hex(fields[3]));
        let material = decode_hex(fields[3]);
        assert_eq!(
            blake3_derive_key(SUBSEED_CONTEXT, &material).as_slice(),
            decode_hex(fields[4])
        );
    }
    assert_eq!(
        scopes,
        BTreeSet::from([
            "biological_process",
            "initial_population",
            "reference_generation",
            "reference_reduction",
            "sampling",
            "technical_event"
        ])
    );
    for forbidden in ["path=", "clock=", "thread=", "position=", "scheduler="] {
        assert!(!SUBSEED_CATALOG.contains(forbidden));
    }
}

#[test]
fn id_vectors_separate_classes_and_encode_the_first_128_digest_bits() {
    let v1_rows: Vec<_> = ID_CATALOG
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    let rows: Vec<_> = ID_CATALOG_V2
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(&rows[..v1_rows.len()], v1_rows);
    assert_eq!(rows.len(), 8);
    let mut identifiers = BTreeSet::new();
    for row in rows {
        let fields: Vec<_> = row.split('|').collect();
        assert_eq!(fields.len(), 9);
        let [
            class,
            prefix,
            _scope,
            causes,
            components,
            root,
            material,
            digest,
            identifier,
        ] = fields.as_slice()
        else {
            panic!("ID catalog row must have nine fields");
        };
        assert!(identifier.starts_with(prefix));
        if *class == "template" {
            assert_eq!(*causes, "library_id,template_ordinal");
            assert_eq!(
                *components,
                "library_id=LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ;template_ordinal=7"
            );
            assert_eq!(*root, "-");
            assert_eq!(*material, "-");
            assert_eq!(*digest, "-");
            assert_eq!(*identifier, "TPL_LIB_EB6VHCIWUACHDO7B5IVJ5RXTCQ_7");
        } else {
            assert!(components.starts_with(&format!("class={class};")));
            assert!(causes.contains("root"));
            assert_eq!(encode_material(root, components), decode_hex(material));
            let material = decode_hex(material);
            assert_eq!(
                blake3_derive_key(ID_CONTEXT, &material).as_slice(),
                decode_hex(digest)
            );
            let encoded = &identifier[prefix.len()..];
            assert!(!encoded.contains('='));
            assert_eq!(decode_base32_no_padding(encoded), decode_hex(&digest[..32]));
        }
        assert!(identifiers.insert(*identifier));
    }
}

#[test]
fn id_catalog_v3_preserves_v2_rows_and_adds_f4_v2_derivations() {
    let v2_rows: Vec<_> = ID_CATALOG_V2
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    let rows: Vec<_> = ID_CATALOG_V3
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(&rows[..v2_rows.len()], v2_rows);
    assert_eq!(rows.len(), v2_rows.len() + 3);

    for row in &rows[v2_rows.len()..] {
        let fields: Vec<_> = row.split('|').collect();
        assert_eq!(fields.len(), 9);
        let [
            class,
            prefix,
            _scope,
            causes,
            components,
            root,
            material,
            digest,
            identifier,
        ] = fields.as_slice()
        else {
            panic!("ID catalog row must have nine fields");
        };
        assert!(matches!(*class, "PreHopInstanceId" | "TechnicalPairId"));
        assert!(matches!(*prefix, "PH_" | "TP_"));
        assert!(components.starts_with(&format!("class={class};")));
        assert!(causes.starts_with("root,event_id,"));
        assert_eq!(
            encode_material_with_version(root, components, 2),
            decode_hex(material)
        );
        assert_eq!(
            blake3_derive_key(ID_CONTEXT, &decode_hex(material)).as_slice(),
            decode_hex(digest)
        );
        assert!(identifier.starts_with(prefix));
        assert_eq!(
            decode_base32_no_padding(&identifier[prefix.len()..]),
            decode_hex(&digest[..32])
        );
    }
}

#[test]
fn f5_rng_catalog_v1_publishes_independent_blake3_and_chacha12_vectors() {
    assert!(F5_RNG_CATALOG.starts_with("# catalog_id=f5-rng-catalog-v1\n"));
    let rows: Vec<_> = F5_RNG_CATALOG
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(rows.len(), 12);
    let mut decisions = BTreeSet::new();
    for row in rows {
        let fields: Vec<_> = row.split('|').collect();
        assert_eq!(fields.len(), 9);
        let [
            decision,
            components,
            root,
            material,
            subseed,
            first_word,
            second_word,
            third_word,
            fourth_word,
        ] = fields.as_slice()
        else {
            panic!("F5 RNG row must have nine fields");
        };
        assert!(decisions.insert(*decision));
        assert!(components.contains(&format!("stage={decision}")));
        assert_eq!(encode_material(root, components), decode_hex(material));
        let expected_subseed = blake3_derive_key(SUBSEED_CONTEXT, &decode_hex(material));
        assert_eq!(expected_subseed.as_slice(), decode_hex(subseed));
        let expected_words =
            chacha12_first_four_u64(&expected_subseed).map(|word| format!("{word:016x}"));
        assert_eq!(
            expected_words.as_slice(),
            [*first_word, *second_word, *third_word, *fourth_word]
        );
    }
    assert_eq!(
        decisions,
        BTreeSet::from([
            "controlled_quota_tie",
            "founder_controlled_permutation",
            "founder_f",
            "founder_sample",
            "initial_column",
            "layout_gap",
            "layout_length",
            "layout_position",
            "marginal_spectrum",
            "polarity",
            "site_alt",
            "swap_proposal",
        ])
    );
    for forbidden in ["path=", "clock=", "thread=", "position=", "scheduler="] {
        assert!(!F5_RNG_CATALOG.contains(forbidden));
    }
}

#[test]
fn id_catalog_v4_is_a_literal_v3_superset_with_five_f5_classes() {
    assert!(ID_CATALOG_V4.starts_with("# catalog_version=4\n"));
    let legacy_rows: Vec<_> = ID_CATALOG_V3
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    let rows: Vec<_> = ID_CATALOG_V4
        .lines()
        .filter(|line| !line.is_empty() && !line.starts_with('#'))
        .collect();
    assert_eq!(&rows[..legacy_rows.len()], legacy_rows);
    assert_eq!(rows.len(), legacy_rows.len() + 5);

    let mut classes = BTreeMap::new();
    let mut identifiers = BTreeSet::new();
    for row in &rows[legacy_rows.len()..] {
        let fields: Vec<_> = row.split('|').collect();
        assert_eq!(fields.len(), 9);
        let [
            class,
            prefix,
            scope,
            _causes,
            components,
            root,
            material,
            digest,
            identifier,
        ] = fields.as_slice()
        else {
            panic!("F5 ID row must have nine fields");
        };
        assert!(identifier.starts_with(prefix));
        assert!(components.starts_with(&format!("class={class};")));
        assert_eq!(encode_material(root, components), decode_hex(material));
        assert_eq!(
            blake3_derive_key(ID_CONTEXT, &decode_hex(material)).as_slice(),
            decode_hex(digest)
        );
        assert_eq!(
            decode_base32_no_padding(&identifier[prefix.len()..]),
            decode_hex(&digest[..32])
        );
        assert!(identifiers.insert(*identifier));
        classes.insert(*class, (*prefix, *scope));
    }

    assert_eq!(
        classes,
        BTreeMap::from([
            ("f5_block", ("BLK_", "execution")),
            (
                "f5_block_haplotype_assignment",
                ("FHA_", "block_founder_copy")
            ),
            ("f5_founder_copy", ("FCP_", "founder_sequence_phase")),
            ("f5_haplotype", ("HAP_", "block_population")),
            ("f5_site", ("SITE_", "block")),
        ])
    );
    assert!(ID_CATALOG_V4.contains("sequence_id=sequence-1;phase_id=phase-0"));
    assert!(ID_CATALOG_V4.contains("founder_copy_id=FCP_"));
    assert!(!ID_CATALOG_V4.contains("copy_index_u32"));
}
