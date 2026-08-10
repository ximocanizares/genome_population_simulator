use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

use genome_population_simulator::application::load_and_validate_config;
use genome_population_simulator::execution::execute_validated_config;

#[derive(Clone, Copy)]
enum DiagnosticMode {
    Human,
    Json,
}

#[derive(Clone, Copy)]
enum CommandMode {
    Validate,
    Run,
}

struct ParsedArgs {
    command: CommandMode,
    path: PathBuf,
    diagnostics: DiagnosticMode,
    output: Option<PathBuf>,
}

fn main() -> ExitCode {
    let parsed = match parse_args(std::env::args().skip(1)) {
        Ok(parsed) => parsed,
        Err(message) => {
            eprintln!("{message}");
            return ExitCode::from(2);
        }
    };
    match load_and_validate_config(&parsed.path).and_then(|effective| match parsed.command {
        CommandMode::Validate => Ok(()),
        CommandMode::Run => parsed.output.as_deref().map_or_else(
            || {
                Err(
                    genome_population_simulator::diagnostics::Diagnostic::internal_invariant(
                        "cli",
                        "run_without_output_after_argument_validation",
                    ),
                )
            },
            |destination| execute_validated_config(&effective, destination).map(|_| ()),
        ),
    }) {
        Ok(()) => ExitCode::SUCCESS,
        Err(diagnostic) => {
            let rendered = match parsed.diagnostics {
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

fn parse_args(mut args: impl Iterator<Item = String>) -> Result<ParsedArgs, &'static str> {
    let command = match args.next().as_deref() {
        Some("validate") => CommandMode::Validate,
        Some("run") => CommandMode::Run,
        _ => {
            return Err(
                "uso: genome-population-simulator <validate|run> <config.toml> [--diagnostics human|json] [--output <ruta>]",
            );
        }
    };
    let path = PathBuf::from(
        args.next()
            .ok_or("falta <config.toml>: use genome-population-simulator validate <config.toml>")?,
    );
    let mut mode = DiagnosticMode::Human;
    let mut output = None;
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
                output = Some(PathBuf::from(args.next().ok_or("--output exige una ruta")?));
            }
            _ => return Err("argumento desconocido"),
        }
    }
    if matches!(command, CommandMode::Run) && output.is_none() {
        return Err("run exige --output <ruta>");
    }
    Ok(ParsedArgs {
        command,
        path,
        diagnostics: mode,
        output,
    })
}
