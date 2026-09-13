"""Invocación local del binario Rust mediante argumentos separados."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import os
import re
import subprocess
import tempfile
from typing import Any, Sequence

SUPPORTED_BINARY_NAME = "genome-population-simulator"
SUPPORTED_BINARY_VERSION = "0.2.0"
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
            and self.binary == _clean_local_path(binary).expanduser().resolve()
            and self.base_directory
                == _clean_local_path(base_directory).expanduser().resolve()
            and self.allow_empty_snv_matrix == allow_empty_snv_matrix
            and self.allow_incompatible_genotypes == allow_incompatible_genotypes
        )


def check_binary_version(binary: Path) -> tuple[Path, str]:
    resolved = _clean_local_path(binary).expanduser().resolve()
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
    declared_destination = _clean_local_path(output_directory).expanduser()
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


# ---------------------------------------------------------------------------
# F5.6.1 — el par general v3 + escenario v1
# ---------------------------------------------------------------------------

#: `relative_path` de los dos esquemas F5, replicado para poder rechazar una
#: ruta de escritura antes de tocar el disco.
#:
#: La regla del **nombre publicable** del general no se replica aquí: vive en
#: la admisión Rust, que la aplica en `validate` y en `run`. Duplicarla en dos
#: lenguajes sin juez es lo que DTC-051 prohíbe, y de hecho divergía: el
#: adaptador la aplicaba y `validate` no.
_RELATIVE_PATH = re.compile(
    r"(?!/)(?!.*\\)(?!.*(?:^|/)\.{1,2}(?:/|$))[A-Za-z0-9._-]+(?:/[A-Za-z0-9._-]+)*"
)


@dataclass(frozen=True)
class ValidatedF5Session:
    """Validación vigente de un par F5. Cualquier edición posterior la invalida.

    Conserva los **bytes** de los dos documentos y no un resumen: la identidad
    de ejecución F5 sale de sus SHA-256, de modo que dos pares semánticamente
    iguales pero distintos byte a byte son ejecuciones distintas y no pueden
    compartir una validación.
    """

    general_toml: str
    scenario_toml: str
    binary: Path
    binary_version: str
    authoring_directory: Path
    general_name: str
    scenario_relative_path: str

    @property
    def general_path(self) -> Path:
        return self.authoring_directory / self.general_name

    def matches(
        self,
        general_toml: str,
        scenario_toml: str,
        binary: Path,
        authoring_directory: Path,
        general_name: str,
        scenario_relative_path: str,
    ) -> bool:
        return (
            self.general_toml == general_toml
            and self.scenario_toml == scenario_toml
            and self.binary == _clean_local_path(binary).expanduser().resolve()
            and self.authoring_directory
                == _clean_local_path(authoring_directory).expanduser().resolve()
            and self.general_name == general_name
            and self.scenario_relative_path == scenario_relative_path
        )


def write_f5_pair(
    authoring_directory: Path,
    general_name: str,
    general_toml: str,
    scenario_relative_path: str,
    scenario_toml: str,
    *,
    overwrite: bool = False,
) -> Path:
    """Escribe el par exportado y devuelve la ruta del general.

    Los dos documentos van al disco **tal cual se exportan**: el escenario se
    enlaza por una ruta relativa que Rust resuelve desde el directorio del
    general, de modo que no pueden vivir en un temporal anónimo. El nombre del
    general se comprueba aquí porque se publica en `effective-config.json` y
    Rust lo rechaza en `run`, no en `validate`.
    """
    directory = _base_directory(authoring_directory)
    if _RELATIVE_PATH.fullmatch(scenario_relative_path) is None:
        raise ProcessAdapterError(
            f"la ruta del escenario debe ser relativa y sin '.' ni '..': {scenario_relative_path!r}"
        )
    general_path = directory / general_name
    scenario_path = directory / Path(scenario_relative_path)
    # Contención, no publicabilidad: el adaptador solo escribe **dentro** del
    # directorio de autoría. Si el nombre además no es publicable, lo rechaza
    # la admisión Rust, con su diagnóstico, en `validate` y en `run`.
    if not _within(directory, general_path) or not _within(directory, scenario_path):
        raise ProcessAdapterError(
            "los dos documentos deben escribirse dentro del directorio de autoría: "
            f"{general_path} y {scenario_path}"
        )
    for path, contents, kind in (
        (general_path, general_toml, "general"),
        (scenario_path, scenario_toml, "f5_scenario"),
    ):
        if not contents.startswith(F5_DOCUMENT_HEADERS[kind]):
            raise ProcessAdapterError(
                f"el documento que se iba a escribir no es un {kind} F5: {path}"
            )
        if not overwrite and not _replaceable_f5_document(path, contents, kind):
            raise ProcessAdapterError(
                "el destino ya contiene un documento que no es el "
                f"{kind} F5 de esta autoría: {path}"
            )
    try:
        scenario_path.parent.mkdir(parents=True, exist_ok=True)
        general_path.write_text(general_toml, encoding="utf-8")
        scenario_path.write_text(scenario_toml, encoding="utf-8")
    except OSError as error:
        raise ProcessAdapterError(f"no se pudo escribir el par F5: {error}") from error
    return general_path


#: Cabecera exacta de cada rol F5: `config_kind` **y** `schema_version`.
#:
#: `config_kind = "general"` no basta para reconocer un documento de esta
#: autoría: es también el de los generales **v2** que las otras ocho pestañas
#: de la misma aplicación editan, y el destino que la interfaz ofrece por
#: omisión es `<repo>/general.toml`. Sin la versión, la primera validación F5
#: machacaría ese documento. La cabecera es la de las dos primeras líneas que
#: emite el renderizador, en su orden.
F5_DOCUMENT_HEADERS: dict[str, str] = {
    "general": 'config_kind = "general"\nschema_version = "3"\n',
    "f5_scenario": 'config_kind = "f5_scenario"\nschema_version = "1"\n',
}


def _within(directory: Path, path: Path) -> bool:
    """Whether `path` stays inside `directory` once every `..` is resolved."""
    try:
        resolved = Path(os.path.normpath(path))
    except (OSError, ValueError):
        return False
    return resolved != directory and directory in resolved.parents


def _replaceable_f5_document(path: Path, contents: str, kind: str) -> bool:
    """Whether writing here would only replace an F5 document of this exact role.

    Escribir la autoría en el árbol de quien la edita no puede convertirse en
    una vía de pérdida de datos: se admite un destino ausente, uno idéntico y
    uno que ya declara **el mismo `config_kind` y la misma `schema_version`**,
    que es el fichero que la propia interfaz escribió antes. Cualquier otro
    contenido —un general v2 incluido— exige que quien edita lo autorice a la
    vista del aviso.
    """
    header = F5_DOCUMENT_HEADERS[kind]
    try:
        if not path.exists():
            return True
        existing = path.read_text(encoding="utf-8")
    except (OSError, UnicodeDecodeError):
        return False
    return existing == contents or existing.startswith(header)


def validate_f5_pair(
    binary: Path,
    general_toml: str,
    scenario_toml: str,
    authoring_directory: Path,
    *,
    general_name: str = "general.toml",
    scenario_relative_path: str = "scenarios/scenario.toml",
    overwrite: bool = False,
) -> tuple[ProcessResult, ValidatedF5Session | None]:
    """Ejecuta el mismo `validate` de Rust sobre el par ya escrito."""
    resolved_binary, version = check_binary_version(binary)
    directory = _base_directory(authoring_directory)
    general_path = write_f5_pair(
        directory,
        general_name,
        general_toml,
        scenario_relative_path,
        scenario_toml,
        overwrite=overwrite,
    )
    result = _invoke_at(resolved_binary, ["validate"], general_path, directory, None)
    session = (
        ValidatedF5Session(
            general_toml,
            scenario_toml,
            resolved_binary,
            version,
            directory,
            general_name,
            scenario_relative_path,
        )
        if result.succeeded
        else None
    )
    return result, session


def run_f5_pair(
    session: ValidatedF5Session,
    general_toml: str,
    scenario_toml: str,
    binary: Path,
    authoring_directory: Path,
    output_directory: Path,
    *,
    general_name: str = "general.toml",
    scenario_relative_path: str = "scenarios/scenario.toml",
) -> ProcessResult:
    """Ejecuta el mismo `run` de Rust tras una validación vigente."""
    if not session.matches(
        general_toml,
        scenario_toml,
        binary,
        authoring_directory,
        general_name,
        scenario_relative_path,
    ):
        raise ProcessAdapterError("la configuración cambió después de validar")
    declared_destination = _clean_local_path(output_directory).expanduser()
    destination = (
        declared_destination
        if declared_destination.is_absolute()
        else session.authoring_directory / declared_destination
    ).resolve()
    if destination.exists():
        raise ProcessAdapterError(f"el destino ya existe: {destination}")
    return _invoke_at(
        session.binary,
        ["run"],
        session.general_path,
        session.authoring_directory,
        destination,
    )


#: Los siete TSV de un paquete F5 de éxito, en el orden del contrato, más la
#: proyección de verdad, que el manifiesto enlaza aparte en `f5.truth`.
F5_PACKAGE_TSV_ARTIFACTS: tuple[tuple[str, str], ...] = (
    ("f5-01-block-layout-v1", "outputs/f5-block-layout.tsv"),
    ("f5-02-sites-v1", "outputs/f5-sites.tsv"),
    ("f5-03-haplotypes-v1", "outputs/f5-haplotypes.tsv"),
    ("f5-04-frequencies-v1", "outputs/f5-frequencies.tsv"),
    ("f5-05-ld-evaluation-v1", "outputs/f5-ld-evaluation.tsv"),
    ("f5-06-founders-v1", "outputs/f5-founders.tsv"),
    ("f5-07-attempts-v1", "outputs/f5-attempts.tsv"),
)
F5_PACKAGE_TRUTH_ARTIFACT = ("truth-jsonl", "outputs/truth.jsonl")


def f5_package_summary(output_directory: Path) -> str:
    """Resume el paquete F5 publicado leyendo solo su manifiesto acotado.

    Ningún TSV ni `truth.jsonl` se abre desde Marimo: el resumen sale de las
    rutas, tamaños y SHA-256 que el manifiesto ya declara y que la publicación
    validó (DTC-051 y REQ-UI-15). El tamaño físico de cada enlace sí se
    comprueba, porque es lo único que distingue un manifiesto coherente de uno
    que describe un árbol que no está.
    """
    manifest = output_directory / "manifest.json"
    maximum_bytes = 2_000_000
    try:
        if manifest.stat().st_size > maximum_bytes:
            return "El resumen F5 no se cargó porque manifest.json supera 2 MB."
        with manifest.open("rb") as source:
            document = json.loads(source.read().decode("utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        return f"No se pudo leer el manifiesto del paquete F5: {error}"
    if not isinstance(document, dict):
        return "El manifiesto publicado no es un objeto JSON."
    f5 = document.get("f5")
    if not isinstance(f5, dict):
        return "El manifiesto publicado no contiene el objeto F5."
    state = document.get("state")
    if state == "failed":
        return _f5_failure_summary(f5)
    if state != "complete":
        return f"El manifiesto publicado declara un estado inesperado: {state!r}"
    artifacts = f5.get("artifacts")
    if not isinstance(artifacts, list):
        return "El manifiesto publicado no contiene un catálogo de artefactos."
    if len(artifacts) != len(F5_PACKAGE_TSV_ARTIFACTS):
        # El contrato de esta función es devolver un mensaje, nunca levantar:
        # el `zip(..., strict=True)` de más abajo rompería ante un manifiesto
        # con un enlace de más o de menos.
        return (
            "El manifiesto publicado no enlaza exactamente "
            f"{len(F5_PACKAGE_TSV_ARTIFACTS)} TSV, sino {len(artifacts)}."
        )
    lines = [
        "Paquete F5 · estado=complete · "
        f"ejecución={document.get('execution_id', '?')} · "
        f"escenario={f5.get('scenario_id', '?')}",
        _f5_attempt_line(f5),
    ]
    links = [*F5_PACKAGE_TSV_ARTIFACTS, F5_PACKAGE_TRUTH_ARTIFACT]
    catalog = [*artifacts, f5.get("truth")]
    for (artifact_id, expected_path), entry in zip(links, catalog, strict=True):
        if (
            not isinstance(entry, dict)
            or entry.get("artifact_id") != artifact_id
            or entry.get("path") != expected_path
        ):
            return f"El manifiesto no enlaza exactamente {expected_path}."
        if not _f5_link_matches(output_directory, entry):
            return f"El enlace de {expected_path} no corresponde al fichero publicado."
        lines.append(
            f"{Path(expected_path).name} · ruta={entry['path']} · "
            f"tamaño={entry['size_bytes']} bytes · SHA-256={entry['sha256']}"
        )
    return "\n".join(lines)


def _f5_attempt_line(f5: dict[str, Any]) -> str:
    summary = f5.get("attempt_summary")
    summary = summary if isinstance(summary, dict) else {}
    sampling = f5.get("sampling")
    subpopulations = f5.get("subpopulations")
    return (
        "layout: intentos={layout} · haplotipos: intentos={haplotypes}, "
        "índice máximo={maximum} · filas de muestreo={rows} · subpoblaciones={populations}"
    ).format(
        layout=summary.get("layout_attempts_used", "?"),
        haplotypes=summary.get("haplotype_attempts_total", "?"),
        maximum=summary.get("max_attempt_index", "?"),
        rows=len(sampling) if isinstance(sampling, list) else "?",
        populations=len(subpopulations) if isinstance(subpopulations, list) else "?",
    )


def _f5_failure_summary(f5: dict[str, Any]) -> str:
    failure = f5.get("f5_failure")
    if not isinstance(failure, dict):
        return "El manifiesto fallido no publica su diagnóstico."
    artifacts = failure.get("artifacts")
    if not isinstance(artifacts, list) or len(artifacts) != 1:
        return "El manifiesto fallido debe enlazar exactamente el diagnóstico."
    counters = failure.get("counters")
    counters = counters if isinstance(counters, dict) else {}
    return (
        f"Paquete F5 · estado=failed · código={failure.get('code', '?')} · "
        f"etapa={failure.get('stage', '?')} · "
        f"población={failure.get('population_id')} · "
        f"bloque={failure.get('block_id')} · "
        f"intentos de layout={counters.get('layout_attempts_used', '?')}"
    )


def _f5_link_matches(output_directory: Path, entry: dict[str, Any]) -> bool:
    size = entry.get("size_bytes")
    digest = entry.get("sha256")
    if (
        not isinstance(size, int)
        or size < 0
        or not isinstance(digest, str)
        or len(digest) != 64
        or any(character not in "0123456789abcdef" for character in digest)
    ):
        return False
    try:
        path = output_directory / str(entry["path"])
        return path.is_file() and path.stat().st_size == size
    except OSError:
        return False


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
        return _invoke_at(
            binary, action, temporary_path, base_directory, output_directory
        )
    except OSError as error:
        raise ProcessAdapterError(f"falló la invocación local: {error}") from error
    finally:
        if temporary_path is not None:
            temporary_path.unlink(missing_ok=True)


def _invoke_at(
    binary: Path,
    action: Sequence[str],
    config_path: Path,
    base_directory: Path,
    output_directory: Path | None,
) -> ProcessResult:
    """Invoke the local binary on a document that already exists on disk."""
    command = [
        str(binary),
        action[0],
        str(config_path),
        *action[1:],
        "--diagnostics",
        "json",
    ]
    if output_directory is not None:
        command.extend(["--output", str(output_directory)])
    try:
        completed = subprocess.run(
            command,
            cwd=base_directory,
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError as error:
        raise ProcessAdapterError(f"falló la invocación local: {error}") from error
    return ProcessResult(
        command=tuple(command),
        returncode=completed.returncode,
        stdout=completed.stdout,
        stderr=completed.stderr,
        diagnostics=_diagnostics(completed.stdout),
        output_directory=output_directory if completed.returncode == 0 else None,
    )


def _base_directory(path: Path) -> Path:
    resolved = _clean_local_path(path).expanduser().resolve()
    if not resolved.is_dir():
        raise ProcessAdapterError(f"el directorio base no existe: {resolved}")
    return resolved


def _clean_local_path(path: Path) -> Path:
    """Discard accidental edge whitespace while preserving spaces in path names."""
    return Path(str(path).strip())


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
