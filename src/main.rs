//! Binario del simulador. Es un crate raíz distinto del de la biblioteca, de
//! modo que repite la prohibición: ver la nota de `src/lib.rs` sobre por qué
//! `forbid(unsafe_code)` vive en los raíces de `src/` y no en `Cargo.toml`.
#![forbid(unsafe_code)]

use std::io::Write;
use std::path::PathBuf;
use std::process::ExitCode;

use genome_population_simulator::application::{
    CliValidatedConfig, execute_validated_f5_config, load_and_validate_cli_config,
};
use genome_population_simulator::execution::{
    ExecutionAuthorizations, execute_validated_config_with_authorizations,
    preflight_validated_config_with_authorizations,
};

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
    authorizations: ExecutionAuthorizations,
}

fn main() -> ExitCode {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if matches!(arguments.as_slice(), [argument] if argument == "--help" || argument == "-h") {
        println!("{}", usage());
        return ExitCode::SUCCESS;
    }
    if arguments.as_slice() == ["--version"] {
        println!("{} {}", env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"));
        return ExitCode::SUCCESS;
    }
    let parsed = match parse_args(arguments.into_iter()) {
        Ok(parsed) => parsed,
        Err(message) => {
            eprintln!("{message}");
            return ExitCode::from(2);
        }
    };
    match load_and_validate_cli_config(&parsed.path).and_then(|effective| match (parsed.command, effective) {
        (CommandMode::Validate, CliValidatedConfig::Legacy(effective)) =>
            preflight_validated_config_with_authorizations(&effective, parsed.authorizations),
        // El dispatcher ya completó la admisión .4 (perfil, primera pasada
        // FASTA y única máscara) sin RNG, destino ni publicación.
        (CommandMode::Validate, CliValidatedConfig::F5(_)) =>
            reject_legacy_authorizations(&parsed.path, parsed.authorizations),
        // F5.6.1 encadena la ejecución productiva. El destino se publica
        // exclusivamente por los envoltorios que sellan.
        (CommandMode::Run, CliValidatedConfig::F5(effective)) => {
            reject_legacy_authorizations(&parsed.path, parsed.authorizations)?;
            parsed.output.as_deref().map_or_else(
                || Err(genome_population_simulator::diagnostics::Diagnostic::internal_invariant(
                    "cli", "run_without_output_after_argument_validation",
                )),
                |destination| execute_validated_f5_config(&effective, destination).map(|_| ()),
            )
        }
        (CommandMode::Run, CliValidatedConfig::Legacy(effective)) => parsed.output.as_deref().map_or_else(
            || Err(genome_population_simulator::diagnostics::Diagnostic::internal_invariant(
                "cli", "run_without_output_after_argument_validation",
            )),
            |destination| execute_validated_config_with_authorizations(
                &effective, destination, parsed.authorizations,
            ).map(|_| ()),
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
            return Err(usage());
        }
    };
    let path = PathBuf::from(
        args.next()
            .ok_or("falta <config.toml>: use genome-population-simulator validate <config.toml>")?,
    );
    let mut mode = DiagnosticMode::Human;
    let mut output = None;
    let mut authorizations = ExecutionAuthorizations::default();
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
            "--allow-empty-snv-matrix" => {
                authorizations.allow_empty_snv_matrix = true;
            }
            "--allow-incompatible-genotypes" => {
                authorizations.allow_incompatible_genotypes = true;
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
        authorizations,
    })
}

/// Las dos autorizaciones estrechas son de la matriz SNV de F2: F5 no tiene
/// ninguna vía que pueda consumirlas.
///
/// Aceptarlas y descartarlas en silencio era inocuo mientras `run` F5 siempre
/// fallaba; con la cadena productiva encadenada, dejaría creer que una
/// ejecución quedó autorizada cuando ninguna autorización se aplicó.
fn reject_legacy_authorizations(
    path: &std::path::Path,
    authorizations: ExecutionAuthorizations,
) -> Result<(), genome_population_simulator::diagnostics::Diagnostic> {
    let rejected = match (
        authorizations.allow_empty_snv_matrix,
        authorizations.allow_incompatible_genotypes,
    ) {
        (false, false) => return Ok(()),
        (true, false) => "--allow-empty-snv-matrix",
        (false, true) => "--allow-incompatible-genotypes",
        (true, true) => "--allow-empty-snv-matrix,--allow-incompatible-genotypes",
    };
    Err(
        genome_population_simulator::diagnostics::Diagnostic::f5_configuration(
            path.display().to_string(),
            "command",
            format!(
                "{rejected} pertenece a la matriz SNV de F2 y F5 no tiene ninguna vía que pueda consumirlo"
            ),
        ),
    )
}

fn usage() -> &'static str {
    "uso: genome-population-simulator <validate|run> <config.toml> [--diagnostics human|json] [--output <ruta>] [--allow-empty-snv-matrix] [--allow-incompatible-genotypes]"
}
