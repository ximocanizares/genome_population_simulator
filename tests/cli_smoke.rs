use std::process::Command;

#[test]
fn help_and_version_are_successful_cli_queries() {
    for argument in ["--help", "-h", "--version"] {
        let output = Command::new(env!("CARGO_BIN_EXE_genome-population-simulator"))
            .arg(argument)
            .output()
            .expect("el binario de integración debe arrancar");

        assert!(
            output.status.success(),
            "{argument} terminó con {:?}: {}",
            output.status.code(),
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!output.stdout.is_empty(), "{argument} no produjo salida");
    }
}
