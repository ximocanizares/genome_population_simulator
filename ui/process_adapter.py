"""Invocación local del binario Rust mediante argumentos separados."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import subprocess
import tempfile
from typing import Any, Sequence

SUPPORTED_BINARY_NAME = "genome-population-simulator"
SUPPORTED_BINARY_VERSION = "0.1.0"
DEFAULT_BINARY = Path(__file__).resolve().parents[1] / "target" / "debug" / SUPPORTED_BINARY_NAME
PROJECT_ROOT = Path(__file__).resolve().parents[1]


class ProcessAdapterError(RuntimeError):
    pass


@dataclass(frozen=True)
class ProcessResult:
    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    diagnostics: dict[str, Any] | None
    output_directory: Path | None

    @property
    def succeeded(self) -> bool:
        return self.returncode == 0


@dataclass(frozen=True)
class ValidatedSession:
    toml: str
    binary: Path
    binary_version: str
    base_directory: Path
    allow_empty_snv_matrix: bool = False
    allow_incompatible_genotypes: bool = False

    def matches(
        self,
        toml: str,
        binary: Path,
        base_directory: Path,
        allow_empty_snv_matrix: bool = False,
        allow_incompatible_genotypes: bool = False,
    ) -> bool:
        return (
            self.toml == toml
            and self.binary == binary.resolve()
            and self.base_directory == base_directory.resolve()
            and self.allow_empty_snv_matrix == allow_empty_snv_matrix
            and self.allow_incompatible_genotypes == allow_incompatible_genotypes
        )


def check_binary_version(binary: Path) -> tuple[Path, str]:
    resolved = binary.expanduser().resolve()
    try:
        completed = subprocess.run(
            [str(resolved), "--version"],
            capture_output=True,
            text=True,
            check=False,
            timeout=10,
        )
    except (OSError, subprocess.TimeoutExpired) as error:
        raise ProcessAdapterError(f"no se puede consultar el binario: {error}") from error
    expected = f"{SUPPORTED_BINARY_NAME} {SUPPORTED_BINARY_VERSION}"
    observed = completed.stdout.strip()
    if completed.returncode != 0 or observed != expected:
        raise ProcessAdapterError(
            f"binario incompatible: se esperaba {expected!r} y se obtuvo {observed!r}"
        )
    return resolved, SUPPORTED_BINARY_VERSION


def validate_toml(
    binary: Path,
    toml: str,
    base_directory: Path,
    *,
    allow_empty_snv_matrix: bool = False,
    allow_incompatible_genotypes: bool = False,
) -> tuple[ProcessResult, ValidatedSession | None]:
    resolved_binary, version = check_binary_version(binary)
    resolved_base = _base_directory(base_directory)
    result = _invoke(
        resolved_binary,
        [
            "validate",
            *(["--allow-empty-snv-matrix"] if allow_empty_snv_matrix else []),
            *(["--allow-incompatible-genotypes"] if allow_incompatible_genotypes else []),
        ],
        toml,
        resolved_base,
        output_directory=None,
    )
    session = (
        ValidatedSession(
            toml,
            resolved_binary,
            version,
            resolved_base,
            allow_empty_snv_matrix,
            allow_incompatible_genotypes,
        )
        if result.succeeded
        else None
    )
    return result, session


def run_toml(
    session: ValidatedSession,
    toml: str,
    binary: Path,
    base_directory: Path,
    output_directory: Path,
    *,
    allow_empty_snv_matrix: bool = False,
    allow_incompatible_genotypes: bool = False,
) -> ProcessResult:
    if not session.matches(
        toml,
        binary,
        base_directory,
        allow_empty_snv_matrix,
        allow_incompatible_genotypes,
    ):
        raise ProcessAdapterError("la configuración cambió después de validar")
    declared_destination = output_directory.expanduser()
    destination = (
        declared_destination
        if declared_destination.is_absolute()
        else session.base_directory / declared_destination
    ).resolve()
    if destination.exists():
        raise ProcessAdapterError(f"el destino ya existe: {destination}")
    return _invoke(
        session.binary,
        [
            "run",
            *(["--allow-empty-snv-matrix"] if allow_empty_snv_matrix else []),
            *(["--allow-incompatible-genotypes"] if allow_incompatible_genotypes else []),
        ],
        toml,
        session.base_directory,
        output_directory=destination,
    )


def _invoke(
    binary: Path,
    action: Sequence[str],
    toml: str,
    base_directory: Path,
    output_directory: Path | None,
) -> ProcessResult:
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w",
            encoding="utf-8",
            prefix=".gps-marimo-",
            suffix=".toml",
            dir=base_directory,
            delete=False,
        ) as temporary:
            temporary.write(toml)
            temporary_path = Path(temporary.name)
        command = [
            str(binary),
            action[0],
            str(temporary_path),
            *action[1:],
            "--diagnostics",
            "json",
        ]
        if output_directory is not None:
            command.extend(["--output", str(output_directory)])
        completed = subprocess.run(
            command,
            cwd=base_directory,
            capture_output=True,
            text=True,
            check=False,
        )
        return ProcessResult(
            command=tuple(command),
            returncode=completed.returncode,
            stdout=completed.stdout,
            stderr=completed.stderr,
            diagnostics=_diagnostics(completed.stdout),
            output_directory=output_directory if completed.returncode == 0 else None,
        )
    except OSError as error:
        raise ProcessAdapterError(f"falló la invocación local: {error}") from error
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _base_directory(path: Path) -> Path:
    resolved = path.expanduser().resolve()
    if not resolved.is_dir():
        raise ProcessAdapterError(f"el directorio base no existe: {resolved}")
    return resolved


def _diagnostics(stdout: str) -> dict[str, Any] | None:
    if not stdout.strip():
        return None
    try:
        document = json.loads(stdout)
    except json.JSONDecodeError:
        return None
    return document if isinstance(document, dict) else None


def diagnostic_summary(document: dict[str, Any] | None) -> str:
    """Locate Rust diagnostics for display without interpreting their semantics."""
    if document is None:
        return ""
    diagnostics = document.get("diagnostics")
    if not isinstance(diagnostics, list):
        return ""
    lines = []
    for item in diagnostics:
        if not isinstance(item, dict):
            continue
        context = item.get("context")
        context = context if isinstance(context, dict) else {}
        field = context.get("field", context.get("fields", "documento"))
        field_text = str(field)
        tab = (
            "General/Ejecución"
            if field_text in {"config_kind", "schema_version", "seed"}
            else "Poblaciones"
            if field_text == "normalization_tolerance"
            or field_text.startswith("population_snvs.")
            else "Individuos"
            if field_text.startswith("initial_populations")
            or field_text.startswith("direct_events")
            or field_text.startswith("initial_individuals")
            else "Genoma"
        )
        location = context.get("location")
        suffix = f" · {location}" if location is not None else ""
        lines.append(
            f"{tab} · {field_text}{suffix} · "
            f"{item.get('code', 'diagnóstico')}: {item.get('message', '')}"
        )
    return "\n".join(lines)


def initial_population_statistics_summary(output_directory: Path) -> str:
    """Read only a bounded effective-config JSON summary, never genome artifacts."""
    effective_config = output_directory / "effective-config.json"
    maximum_bytes = 2_000_000
    try:
        if effective_config.stat().st_size > maximum_bytes:
            return (
                "El resumen de fundadores no se cargó porque effective-config.json "
                "supera 2 MB; los artefactos genómicos no se abren desde Marimo."
            )
        with effective_config.open("rb") as source:
            document = json.loads(source.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return f"No se pudo leer el resumen de fundadores publicado: {error}"
    if not isinstance(document, dict):
        return "La configuración efectiva publicada no es un objeto JSON."
    population_snvs = document.get("population_snvs")
    if not isinstance(population_snvs, dict):
        return "La ejecución no publicó estadísticas de fundadores F2.3.5."
    statistics = population_snvs.get("initial_population_statistics")
    if not isinstance(statistics, dict):
        return "La ejecución no publicó estadísticas de fundadores F2.3.5."
    populations = statistics.get("populations")
    if not isinstance(populations, list):
        return "El resumen de fundadores publicado no tiene poblaciones legibles."
    lines = []
    for population in populations[:20]:
        if not isinstance(population, dict):
            continue
        genotype_model = population.get("genotype_model", "?")
        if isinstance(genotype_model, dict):
            genotype_model = genotype_model.get("kind", "?")
        lines.append(
            "{population_id}: N={individual_count}, ploidía={ploidy}, "
            "modelo={genotype_model}, loci={loci}".format(
                population_id=population.get("population_id", "?"),
                individual_count=population.get("individual_count", "?"),
                ploidy=population.get("ploidy", "?"),
                genotype_model=genotype_model,
                loci=len(population.get("loci", []))
                if isinstance(population.get("loci"), list)
                else "?",
            )
        )
    if len(populations) > 20:
        lines.append(f"… y {len(populations) - 20} poblaciones más.")
    return "\n".join(lines) or "No hay poblaciones iniciales en el resumen F2.3.5."


def truth_artifact_summary(output_directory: Path) -> str:
    """Read a bounded manifest entry; `truth.jsonl` is never opened by Marimo."""
    manifest = output_directory / "manifest.json"
    maximum_bytes = 2_000_000
    try:
        if manifest.stat().st_size > maximum_bytes:
            return "El resumen de verdad no se cargó porque manifest.json supera 2 MB."
        with manifest.open("rb") as source:
            document = json.loads(source.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return f"No se pudo leer el manifiesto de resultados: {error}"
    if not isinstance(document, dict):
        return "El manifiesto publicado no es un objeto JSON."
    artifacts = document.get("artifacts")
    if not isinstance(artifacts, list):
        return "El manifiesto publicado no contiene un catálogo de artefactos."
    matches = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict)
        and artifact.get("id") == "truth-jsonl"
        and artifact.get("kind") == "truth_jsonl"
    ]
    if len(matches) != 1:
        return "El manifiesto no enlaza exactamente un artefacto truth.jsonl válido."
    artifact = matches[0]
    path = artifact.get("path")
    size = artifact.get("size_bytes")
    checksum = artifact.get("checksum")
    algorithm = checksum.get("algorithm") if isinstance(checksum, dict) else None
    digest = checksum.get("digest") if isinstance(checksum, dict) else None
    if (
        path != "outputs/truth.jsonl"
        or not isinstance(size, int)
        or size < 0
        or algorithm != "sha256"
        or not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        return "El enlace de truth.jsonl en el manifiesto está incompleto o es inválido."
    try:
        truth_path = output_directory / path
        if not truth_path.is_file() or truth_path.stat().st_size != size:
            return "El enlace de truth.jsonl en el manifiesto está incompleto o es inválido."
    except OSError:
        return "El enlace de truth.jsonl en el manifiesto está incompleto o es inválido."
    return f"truth.jsonl · ruta={path} · tamaño={size} bytes · SHA-256={digest}"


def vcf_artifact_summary(output_directory: Path) -> str:
    """Read bounded manifest links for the optional VCF without opening it."""
    manifest = output_directory / "manifest.json"
    maximum_bytes = 2_000_000
    try:
        if manifest.stat().st_size > maximum_bytes:
            return "El resumen VCF no se cargó porque manifest.json supera 2 MB."
        with manifest.open("rb") as source:
            document = json.loads(source.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return f"No se pudo leer el manifiesto de resultados: {error}"
    if not isinstance(document, dict):
        return "El manifiesto publicado no es un objeto JSON."
    outputs = document.get("vcf_outputs")
    artifacts = document.get("artifacts")
    if not isinstance(outputs, list) or not isinstance(artifacts, list):
        return "El manifiesto publicado no contiene el catálogo VCF esperado."
    if not outputs:
        return "VCF no solicitado: el manifiesto no publica una vista VCF."
    if len(outputs) != 1 or not isinstance(outputs[0], dict):
        return "El manifiesto no enlaza exactamente una salida VCF válida."

    output = outputs[0]
    mode = output.get("mode")
    vcf_id = output.get("vcf_artifact_id")
    index_kind = output.get("index_kind")
    index_id = output.get("index_artifact_id")
    if output.get("state") != "complete" or mode not in {"text", "bgzf"} or vcf_id != "truth-vcf":
        return "El enlace VCF en el manifiesto está incompleto o es inválido."
    if mode == "text" and (index_kind is not None or index_id is not None):
        return "El enlace VCF en el manifiesto está incompleto o es inválido."
    if mode == "bgzf" and (index_kind not in {"tbi", "csi"} or index_id != "truth-vcf-index"):
        return "El enlace VCF en el manifiesto está incompleto o es inválido."

    expected_vcf_path = "outputs/truth.vcf" if mode == "text" else "outputs/truth.vcf.gz"
    vcf = _cataloged_artifact(artifacts, "truth-vcf", "truth_vcf", expected_vcf_path)
    if vcf is None or not _physical_artifact_matches(output_directory, vcf):
        return "El enlace VCF en el manifiesto está incompleto o es inválido."
    summary = _artifact_summary("truth.vcf" if mode == "text" else "truth.vcf.gz", vcf)

    if mode == "text":
        return summary
    expected_index_path = f"outputs/truth.vcf.gz.{index_kind}"
    expected_index_kind = f"vcf_{index_kind}_index"
    index = _cataloged_artifact(artifacts, "truth-vcf-index", expected_index_kind, expected_index_path)
    if index is None or not _physical_artifact_matches(output_directory, index):
        return "El enlace VCF en el manifiesto está incompleto o es inválido."
    return f"{summary}\níndice {index_kind.upper()} · {_artifact_summary(Path(expected_index_path).name, index)}"


def _cataloged_artifact(
    artifacts: list[Any], artifact_id: str, kind: str, expected_path: str
) -> dict[str, Any] | None:
    matches = [
        artifact
        for artifact in artifacts
        if isinstance(artifact, dict)
        and artifact.get("id") == artifact_id
        and artifact.get("kind") == kind
        and artifact.get("path") == expected_path
    ]
    return matches[0] if len(matches) == 1 else None


def _physical_artifact_matches(output_directory: Path, artifact: dict[str, Any]) -> bool:
    size = artifact.get("size_bytes")
    checksum = artifact.get("checksum")
    algorithm = checksum.get("algorithm") if isinstance(checksum, dict) else None
    digest = checksum.get("digest") if isinstance(checksum, dict) else None
    if (
        not isinstance(size, int)
        or size < 0
        or algorithm != "sha256"
        or not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        return False
    try:
        path = output_directory / artifact["path"]
        return path.is_file() and path.stat().st_size == size
    except OSError:
        return False


def _artifact_summary(label: str, artifact: dict[str, Any]) -> str:
    checksum = artifact["checksum"]
    return (
        f"{label} · ruta={artifact['path']} · tamaño={artifact['size_bytes']} bytes "
        f"· SHA-256={checksum['digest']}"
    )
