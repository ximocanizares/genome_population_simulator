use std::fs;
use std::process::Command;

use genome_population_simulator::integrity::{Sha256Writer, identify_file, verify_identity};

#[test]
fn sha256_matches_literal_vectors_incrementally() {
    let mut empty = Sha256Writer::default();
    empty.update(b"").expect("empty chunk");
    assert_eq!(
        empty.finish("empty").digest(),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );

    let mut abc = Sha256Writer::default();
    abc.update(b"a").expect("a");
    abc.update(b"bc").expect("bc");
    let identity = abc.finish("abc.txt");
    assert_eq!(identity.size_bytes(), 3);
    assert_eq!(identity.algorithm(), "sha256");
    assert_eq!(
        identity.digest(),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
}

#[test]
fn file_identity_matches_sha256sum_and_detects_change_or_truncation() {
    let root = std::env::temp_dir().join(format!("gps-integrity-{}", std::process::id()));
    fs::create_dir_all(&root).expect("temporary root");
    let path = root.join("binary.dat");
    fs::write(&path, [0_u8, 1, 2, 0xff, b'\n']).expect("fixture");
    let identity = identify_file(&path, "artifacts/binary.dat").expect("identity");
    let output = Command::new("sha256sum")
        .arg(&path)
        .output()
        .expect("sha256sum available from A7");
    assert!(output.status.success());
    let external = String::from_utf8(output.stdout).expect("ASCII output");
    assert_eq!(
        identity.digest(),
        external.split_whitespace().next().expect("digest")
    );

    fs::write(&path, [0_u8, 1, 3, 0xff, b'\n']).expect("one-byte change");
    assert_eq!(
        verify_identity(&path, &identity)
            .expect_err("change")
            .code(),
        "GPS-IO-002"
    );
    fs::write(&path, [0_u8, 1]).expect("truncate");
    assert_eq!(
        verify_identity(&path, &identity)
            .expect_err("truncate")
            .code(),
        "GPS-IO-002"
    );
    fs::remove_dir_all(&root).expect("cleanup");
}
