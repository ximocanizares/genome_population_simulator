use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

use genome_population_simulator::application::load_and_validate_config;

#[derive(Clone, Copy)]
enum DiagnosticMode {
    Human,
    Json,
}

fn main() -> ExitCode {
    let (path, mode) = match parse_args(std::env::args().skip(1)) {
        Ok(parsed) => parsed,
        Err(message) => {
            eprintln!("{message}");
            return ExitCode::from(2);
        }
    };
    match load_and_validate_config(&path) {
        Ok(_) => ExitCode::SUCCESS,
        Err(diagnostic) => {
            let rendered = match mode {
                DiagnosticMode::Human => {
                    writeln!(std::io::stderr().lock(), "{}", diagnostic.render_human())
                }
                DiagnosticMode::Json => writeln!(
                    std::io::stdout().lock(),
                    "{}",
                    diagnostic.render_json_document()
                ),
            };
            if rendered.is_ok() {
                ExitCode::from(diagnostic.exit_code())
            } else {
                ExitCode::from(70)
            }
        }
    }
}

fn parse_args(
    mut args: impl Iterator<Item = String>,
) -> Result<(PathBuf, DiagnosticMode), &'static str> {
    if args.next().as_deref() != Some("validate") {
        return Err(
            "uso: genome-population-simulator validate <config.toml> [--diagnostics human|json] [--output <ruta>]",
        );
    }
    let path = PathBuf::from(
        args.next()
            .ok_or("falta <config.toml>: use genome-population-simulator validate <config.toml>")?,
    );
    let mut mode = DiagnosticMode::Human;
    while let Some(argument) = args.next() {
        match argument.as_str() {
            "--diagnostics" => {
                mode = match args.next().as_deref() {
                    Some("human") => DiagnosticMode::Human,
                    Some("json") => DiagnosticMode::Json,
                    _ => return Err("--diagnostics exige human o json"),
                };
            }
            "--output" => {
                let _unused_destination = args.next().ok_or("--output exige una ruta")?;
            }
            _ => return Err("argumento desconocido"),
        }
    }
    Ok((path, mode))
}
