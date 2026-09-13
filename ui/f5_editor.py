"""Modelo de autoría F5 para Marimo: general v3 y escenario F5 v1.

No contiene semántica biológica ni normativa. Rust conserva la autoridad
(DTC-051): aquí solo vive la estructura de los dos documentos —tipos, enums,
las tres expresiones regulares del esquema y sus cotas— y el ida y vuelta que
REQ-UI-15 exige conservar byte a byte.

El renderizado se escribe a mano en lugar de delegarlo en `tomlkit`. La razón
es la paridad: `effective-config.json` publica el nombre del general y la
identidad de ejecución sale del SHA-256 de los dos documentos, de modo que los
bytes exportados **son** parte del contrato. Escribirlos explícitamente deja el
orden de claves, los espacios y los literales decimales bajo control de una
sola función, que es lo que una prueba puede fijar contra un TOML manual.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import re
from typing import Any, Iterable, Mapping, Sequence

import tomlkit
from tomlkit.items import Float, Integer

from config_editor import EditorError

U64_MAX = 18_446_744_073_709_551_615

F5_GENERAL_CONFIG_KIND = "general"
F5_GENERAL_SCHEMA_VERSION = "3"
F5_SCENARIO_CONFIG_KIND = "f5_scenario"
F5_SCENARIO_SCHEMA_VERSION = "1"

F5_SAMPLING_MODES = ("probabilistic", "controlled")
F5_FOUNDER_KINDS = ("founder_ids", "founder_keys")
F5_TARGET_KINDS = ("block_count", "covered_bases", "coverage_fraction")

#: Las tres expresiones regulares de los dos esquemas, sin sus anclajes: aquí
#: se aplican con `fullmatch`. Se exponen para que una prueba las contraste
#: contra los `pattern` declarados, igual que `F5_INTEGER_BOUNDS` se contrasta
#: contra `minimum`/`maximum`: una regla replicada sin juez es la que DTC-051
#: prohíbe.
F5_ID_PATTERN = re.compile(r"[A-Za-z0-9][A-Za-z0-9._-]{0,63}")
F5_SHA_PATTERN = re.compile(r"[0-9a-f]{64}")
F5_RELATIVE_PATH_PATTERN = re.compile(
    r"(?!/)(?!.*\\)(?!.*(?:^|/)\.{1,2}(?:/|$))[A-Za-z0-9._-]+(?:/[A-Za-z0-9._-]+)*"
)
_ID = F5_ID_PATTERN
_SHA = F5_SHA_PATTERN
_RELATIVE_PATH = F5_RELATIVE_PATH_PATTERN

#: Cotas enteras de los dos esquemas, en el orden en que se publican.
#:
#: Viven en una tabla y no repartidas por el código para que una prueba pueda
#: contrastarlas mecánicamente contra `minimum`/`maximum` de los ficheros
#: `schemas/config-general-v3.schema.json` y
#: `schemas/config-f5-scenario-v1.schema.json`. Duplicar una cota sin juez es
#: exactamente la regla duplicada que DTC-051 prohíbe.
F5_INTEGER_BOUNDS: dict[str, tuple[int, int]] = {
    "ploidy": (1, 10),
    "founder_count": (1, U64_MAX),
    "individual_count": (1, U64_MAX),
    "max_layout_attempts": (1, 100_000),
    "max_haplotype_attempts": (1, 10_000),
    "max_swap_proposals": (1, 10_000_000),
    "max_haplotype_cells": (1, 67_108_864),
    "max_operation_units": (1, 268_435_456),
    "max_working_memory_bytes": (1, 1_073_741_824),
    "max_output_bytes": (1, 2_147_483_648),
    "block_count": (1, 4_294_967_295),
    "covered_bases": (1, U64_MAX),
    "block_index": (0, 4_294_967_295),
}

#: Los siete presupuestos y límites algorítmicos opcionales del objeto `f5`.
F5_BUDGET_FIELDS = (
    "max_layout_attempts",
    "max_haplotype_attempts",
    "max_swap_proposals",
    "max_haplotype_cells",
    "max_operation_units",
    "max_working_memory_bytes",
    "max_output_bytes",
)

#: Matriz campo→control del general v3 (ACE-114). `config_kind` y
#: `schema_version` no tienen control: son constantes del documento.
F5_GENERAL_FIELD_CONTROLS: dict[str, str] = {
    "seed": "f5_seed_control",
    "ploidy": "f5_ploidy_control",
    "reference_id": "f5_reference_id_control",
    "reference_path": "f5_reference_path_control",
    "reference_sha256": "f5_reference_sha256_control",
    "f5_scenario": "f5_scenario_link_controls",
}

#: Matriz campo→control del escenario v1. Las tres claves del objetivo
#: exclusivo comparten el par selector+valor, que es como se autoría un
#: `oneOf`.
F5_SCENARIO_FIELD_CONTROLS: dict[str, str] = {
    "scenario_id": "f5_scenario_id_control",
    "founder_count": "f5_founder_count_control",
    "profile_metadata_path": "f5_profile_metadata_path_control",
    "profile_metadata_sha256": "f5_profile_metadata_sha256_control",
    "excluded_regions_path": "f5_excluded_regions_path_control",
    "excluded_regions_sha256": "f5_excluded_regions_sha256_control",
    "inbreeding_f": "f5_inbreeding_control",
    "max_layout_attempts": "f5_max_layout_attempts_control",
    "max_haplotype_attempts": "f5_max_haplotype_attempts_control",
    "max_swap_proposals": "f5_max_swap_proposals_control",
    "max_haplotype_cells": "f5_max_haplotype_cells_control",
    "max_operation_units": "f5_max_operation_units_control",
    "max_working_memory_bytes": "f5_max_working_memory_bytes_control",
    "max_output_bytes": "f5_max_output_bytes_control",
    "block_count": "f5_target_kind_control",
    "covered_bases": "f5_target_kind_control",
    "coverage_fraction": "f5_target_value_control",
    "subpopulations": "f5_subpopulations_editor",
    "sampling_overrides": "f5_sampling_overrides_editor",
}


@dataclass(eq=True)
class F5SubpopulationRow:
    population_id: str
    individual_count: int
    founder_kind: str
    founders: list[str]
    sampling: str | None = None
    inbreeding_f: str | None = None


@dataclass(eq=True)
class F5SamplingOverrideRow:
    population_id: str
    block_index: int
    mode: str
    controlled_max_abs_error: str | None = None


@dataclass(eq=True)
class F5GeneralForm:
    seed: str | None
    ploidy: int
    reference_id: str
    reference_path: str
    reference_sha256: str
    scenario_id: str
    scenario_path: str


@dataclass(eq=True)
class F5ScenarioForm:
    scenario_id: str
    founder_count: int
    profile_metadata_path: str
    profile_metadata_sha256: str
    target_kind: str
    target_value: str
    subpopulations: list[F5SubpopulationRow] = field(default_factory=list)
    sampling_overrides: list[F5SamplingOverrideRow] = field(default_factory=list)
    excluded_regions_path: str | None = None
    excluded_regions_sha256: str | None = None
    inbreeding_f: str | None = None
    max_layout_attempts: int | None = None
    max_haplotype_attempts: int | None = None
    max_swap_proposals: int | None = None
    max_haplotype_cells: int | None = None
    max_operation_units: int | None = None
    max_working_memory_bytes: int | None = None
    max_output_bytes: int | None = None


@dataclass(eq=True)
class F5Pair:
    """Los dos documentos de una autoría F5, siempre juntos.

    `general_source` y `scenario_source` conservan los **bytes importados**.
    REQ-UI-15 exige que importar y exportar conserve «exactamente el mismo
    TOML», y el renderizador impone un orden de claves único: un escenario
    legal escrito en otro orden no volvería igual, y como el `execution_id`
    sale de los SHA-256 observados, la misma autoría dejaría de ser la misma
    ejecución. Los bytes se reemiten mientras sigan describiendo **exactamente**
    el formulario vigente; en cuanto el formulario se edita, el renderizador
    toma el relevo. La condición se comprueba reparseando el origen, no se
    supone.
    """

    general: F5GeneralForm
    scenario: F5ScenarioForm
    general_source: str | None = None
    scenario_source: str | None = None


def default_f5_pair() -> F5Pair:
    """Par mínimo y editable. No se valida solo: Rust decide si es viable."""
    return F5Pair(
        general=F5GeneralForm(
            seed=None,
            ploidy=2,
            reference_id="reference-main",
            reference_path="references/reference.fa",
            reference_sha256="0" * 64,
            scenario_id="scenario-1",
            scenario_path="scenarios/scenario.toml",
        ),
        scenario=F5ScenarioForm(
            scenario_id="scenario-1",
            founder_count=1,
            profile_metadata_path="profiles/profile.metadata.json",
            profile_metadata_sha256="0" * 64,
            target_kind="block_count",
            target_value="1",
            subpopulations=[
                F5SubpopulationRow(
                    population_id="population-1",
                    individual_count=1,
                    founder_kind="founder_ids",
                    founders=["individual-1"],
                )
            ],
            sampling_overrides=[],
        ),
    )


# ---------------------------------------------------------------------------
# Importación: exactamente un general v3 y exactamente un escenario v1
# ---------------------------------------------------------------------------


def classify_f5_document(text: str, name: str) -> str:
    """Devuelve `general_v3`, `scenario_v1` o levanta el rechazo del contrato.

    Un general v2 no es una variante de autoría F5: pertenece a la pestaña que
    lo consume. Una v4 no existe. Los dos se nombran por su versión para que el
    rechazo sea accionable.
    """
    try:
        document = tomlkit.parse(text)
    except Exception as error:  # noqa: BLE001 - tomlkit no expone un tipo único
        raise EditorError(name, f"TOML inválido: {error}") from error
    kind = document.get("config_kind")
    version = document.get("schema_version")
    if not isinstance(kind, str) or not isinstance(version, str):
        raise EditorError(
            name, "faltan config_kind y schema_version como texto en la raíz"
        )
    kind, version = str(kind), str(version)
    if kind == F5_GENERAL_CONFIG_KIND and version == F5_GENERAL_SCHEMA_VERSION:
        return "general_v3"
    if kind == F5_SCENARIO_CONFIG_KIND and version == F5_SCENARIO_SCHEMA_VERSION:
        return "scenario_v1"
    if kind == F5_GENERAL_CONFIG_KIND:
        raise EditorError(
            name,
            f'F5 exige un general v3 y este declara la versión "{version}"',
        )
    if kind == F5_SCENARIO_CONFIG_KIND:
        raise EditorError(
            name,
            f'F5 exige un escenario v1 y este declara la versión "{version}"',
        )
    raise EditorError(name, f'config_kind no admitido en F5: "{kind}"')


def import_f5_pair(documents: Sequence[tuple[str, str]]) -> F5Pair:
    """Importa exactamente un general v3 con exactamente un escenario v1.

    Cero o varios escenarios son errores, no variantes de autoría (REQ-UI-15).
    La cardinalidad se juzga **antes** de parsear el contenido para que el
    diagnóstico nombre el problema real y no el primer campo que falle.
    """
    if not documents:
        raise EditorError("import", "se exigen dos documentos: general v3 y escenario v1")
    generals: list[tuple[str, str]] = []
    scenarios: list[tuple[str, str]] = []
    for name, text in documents:
        role = classify_f5_document(text, name)
        (generals if role == "general_v3" else scenarios).append((name, text))
    _require_exactly_one(generals, "general v3")
    _require_exactly_one(scenarios, "escenario F5 v1")
    general = parse_f5_general_toml(generals[0][1])
    scenario = parse_f5_scenario_toml(scenarios[0][1])
    return F5Pair(
        general=general,
        scenario=scenario,
        general_source=generals[0][1],
        scenario_source=scenarios[0][1],
    )


def _require_exactly_one(found: Sequence[tuple[str, str]], label: str) -> None:
    if len(found) == 1:
        return
    if not found:
        raise EditorError("import", f"falta el {label}: se exige exactamente uno")
    names = ", ".join(sorted(name for name, _ in found))
    raise EditorError(
        "import", f"se exige exactamente un {label} y se recibieron {len(found)}: {names}"
    )


# ---------------------------------------------------------------------------
# Parseo
# ---------------------------------------------------------------------------


def parse_f5_general_toml(text: str) -> F5GeneralForm:
    document = _document(text, "general")
    _keys(
        document,
        allowed={
            "config_kind",
            "schema_version",
            "seed",
            "ploidy",
            "reference_id",
            "reference_path",
            "reference_sha256",
            "f5_scenario",
        },
        required={
            "config_kind",
            "schema_version",
            "ploidy",
            "reference_id",
            "reference_path",
            "reference_sha256",
            "f5_scenario",
        },
        path="general",
    )
    if _text(document["config_kind"], "general.config_kind") != F5_GENERAL_CONFIG_KIND:
        raise EditorError("general.config_kind", 'debe ser "general"')
    if (
        _text(document["schema_version"], "general.schema_version")
        != F5_GENERAL_SCHEMA_VERSION
    ):
        raise EditorError("general.schema_version", 'debe ser "3"')
    link = _table(document["f5_scenario"], "general.f5_scenario")
    _keys(
        link,
        allowed={"scenario_id", "scenario_path"},
        required={"scenario_id", "scenario_path"},
        path="general.f5_scenario",
    )
    return F5GeneralForm(
        seed=_sha256(document["seed"], "general.seed") if "seed" in document else None,
        ploidy=_bounded_integer(document["ploidy"], "general.ploidy", "ploidy"),
        reference_id=_identifier(document["reference_id"], "general.reference_id"),
        reference_path=_relative_path(
            document["reference_path"], "general.reference_path"
        ),
        reference_sha256=_sha256(
            document["reference_sha256"], "general.reference_sha256"
        ),
        scenario_id=_identifier(link["scenario_id"], "general.f5_scenario.scenario_id"),
        scenario_path=_relative_path(
            link["scenario_path"], "general.f5_scenario.scenario_path"
        ),
    )


def parse_f5_scenario_toml(text: str) -> F5ScenarioForm:
    document = _document(text, "scenario")
    _keys(
        document,
        allowed={
            "config_kind",
            "schema_version",
            "scenario_id",
            "founder_count",
            "f5",
            "subpopulations",
            "sampling_overrides",
        },
        required={
            "config_kind",
            "schema_version",
            "scenario_id",
            "founder_count",
            "f5",
            "subpopulations",
        },
        path="scenario",
    )
    if _text(document["config_kind"], "scenario.config_kind") != F5_SCENARIO_CONFIG_KIND:
        raise EditorError("scenario.config_kind", 'debe ser "f5_scenario"')
    if (
        _text(document["schema_version"], "scenario.schema_version")
        != F5_SCENARIO_SCHEMA_VERSION
    ):
        raise EditorError("scenario.schema_version", 'debe ser "1"')
    f5 = _table(document["f5"], "scenario.f5")
    _keys(
        f5,
        allowed={
            "profile_metadata_path",
            "profile_metadata_sha256",
            "excluded_regions_path",
            "excluded_regions_sha256",
            "inbreeding_f",
            *F5_BUDGET_FIELDS,
            *F5_TARGET_KINDS,
        },
        required={"profile_metadata_path", "profile_metadata_sha256"},
        path="scenario.f5",
    )
    target_kind, target_value = _target(f5)
    excluded_path = (
        _relative_path(f5["excluded_regions_path"], "scenario.f5.excluded_regions_path")
        if "excluded_regions_path" in f5
        else None
    )
    excluded_sha = (
        _sha256(f5["excluded_regions_sha256"], "scenario.f5.excluded_regions_sha256")
        if "excluded_regions_sha256" in f5
        else None
    )
    if (excluded_path is None) != (excluded_sha is None):
        raise EditorError(
            "scenario.f5.excluded_regions_path,scenario.f5.excluded_regions_sha256",
            "la máscara exige ruta y hash juntos",
        )
    budgets = {
        name: (
            _bounded_integer(f5[name], f"scenario.f5.{name}", name)
            if name in f5
            else None
        )
        for name in F5_BUDGET_FIELDS
    }
    return F5ScenarioForm(
        scenario_id=_identifier(document["scenario_id"], "scenario.scenario_id"),
        founder_count=_bounded_integer(
            document["founder_count"], "scenario.founder_count", "founder_count"
        ),
        profile_metadata_path=_relative_path(
            f5["profile_metadata_path"], "scenario.f5.profile_metadata_path"
        ),
        profile_metadata_sha256=_sha256(
            f5["profile_metadata_sha256"], "scenario.f5.profile_metadata_sha256"
        ),
        target_kind=target_kind,
        target_value=target_value,
        subpopulations=_parse_subpopulations(document["subpopulations"]),
        sampling_overrides=_parse_sampling_overrides(
            document.get("sampling_overrides", [])
        ),
        excluded_regions_path=excluded_path,
        excluded_regions_sha256=excluded_sha,
        inbreeding_f=(
            _unit_interval(f5["inbreeding_f"], "scenario.f5.inbreeding_f")
            if "inbreeding_f" in f5
            else None
        ),
        **budgets,
    )


def _target(f5: Mapping[str, Any]) -> tuple[str, str]:
    declared = [kind for kind in F5_TARGET_KINDS if kind in f5]
    if len(declared) != 1:
        raise EditorError(
            "scenario.f5." + "|".join(F5_TARGET_KINDS),
            f"se exige exactamente un objetivo y se declararon {len(declared)}",
        )
    kind = declared[0]
    path = f"scenario.f5.{kind}"
    if kind == "coverage_fraction":
        literal = _number_literal(f5[kind], path)
        value = float(literal)
        if not 0 < value <= 1:
            raise EditorError(path, "debe pertenecer al intervalo (0, 1]")
        return kind, literal
    return kind, str(_bounded_integer(f5[kind], path, kind))


def _parse_subpopulations(value: Any) -> list[F5SubpopulationRow]:
    rows = _array(value, "scenario.subpopulations")
    if not rows:
        raise EditorError("scenario.subpopulations", "se exige al menos una subpoblación")
    parsed: list[F5SubpopulationRow] = []
    for index, item in enumerate(rows):
        path = f"scenario.subpopulations[{index}]"
        table = _table(item, path)
        _keys(
            table,
            allowed={
                "population_id",
                "individual_count",
                "sampling",
                "inbreeding_f",
                *F5_FOUNDER_KINDS,
            },
            required={"population_id", "individual_count"},
            path=path,
        )
        declared = [kind for kind in F5_FOUNDER_KINDS if kind in table]
        if len(declared) != 1:
            raise EditorError(
                f"{path}.founder_ids,{path}.founder_keys",
                "se exige exactamente uno",
            )
        founder_kind = declared[0]
        parsed.append(
            F5SubpopulationRow(
                population_id=_identifier(table["population_id"], f"{path}.population_id"),
                individual_count=_bounded_integer(
                    table["individual_count"],
                    f"{path}.individual_count",
                    "individual_count",
                ),
                founder_kind=founder_kind,
                founders=_identifier_list(
                    table[founder_kind], f"{path}.{founder_kind}"
                ),
                sampling=(
                    _sampling_mode(table["sampling"], f"{path}.sampling")
                    if "sampling" in table
                    else None
                ),
                inbreeding_f=(
                    _unit_interval(table["inbreeding_f"], f"{path}.inbreeding_f")
                    if "inbreeding_f" in table
                    else None
                ),
            )
        )
    return parsed


def _parse_sampling_overrides(value: Any) -> list[F5SamplingOverrideRow]:
    rows = _array(value, "scenario.sampling_overrides")
    parsed: list[F5SamplingOverrideRow] = []
    for index, item in enumerate(rows):
        path = f"scenario.sampling_overrides[{index}]"
        table = _table(item, path)
        _keys(
            table,
            allowed={
                "population_id",
                "block_index",
                "mode",
                "controlled_max_abs_error",
            },
            required={"population_id", "block_index", "mode"},
            path=path,
        )
        mode = _sampling_mode(table["mode"], f"{path}.mode")
        tolerance = (
            _unit_interval(
                table["controlled_max_abs_error"], f"{path}.controlled_max_abs_error"
            )
            if "controlled_max_abs_error" in table
            else None
        )
        if mode == "probabilistic" and tolerance is not None:
            raise EditorError(
                f"{path}.controlled_max_abs_error",
                "probabilistic no admite tolerancia",
            )
        parsed.append(
            F5SamplingOverrideRow(
                population_id=_identifier(table["population_id"], f"{path}.population_id"),
                block_index=_bounded_integer(
                    table["block_index"], f"{path}.block_index", "block_index"
                ),
                mode=mode,
                controlled_max_abs_error=tolerance,
            )
        )
    return parsed


# ---------------------------------------------------------------------------
# Renderizado
# ---------------------------------------------------------------------------


def render_f5_general_toml(form: F5GeneralForm) -> str:
    lines = [
        f'config_kind = "{F5_GENERAL_CONFIG_KIND}"',
        f'schema_version = "{F5_GENERAL_SCHEMA_VERSION}"',
    ]
    if form.seed is not None:
        lines.append(f'seed = "{_sha256(form.seed, "general.seed")}"')
    lines.append(f"ploidy = {_bounded_integer(form.ploidy, 'general.ploidy', 'ploidy')}")
    lines.append(
        f'reference_id = "{_identifier(form.reference_id, "general.reference_id")}"'
    )
    lines.append(
        f'reference_path = "{_relative_path(form.reference_path, "general.reference_path")}"'
    )
    lines.append(
        f'reference_sha256 = "{_sha256(form.reference_sha256, "general.reference_sha256")}"'
    )
    lines.append("")
    lines.append("[f5_scenario]")
    lines.append(
        f'scenario_id = "{_identifier(form.scenario_id, "general.f5_scenario.scenario_id")}"'
    )
    lines.append(
        f'scenario_path = "{_relative_path(form.scenario_path, "general.f5_scenario.scenario_path")}"'
    )
    return "\n".join(lines) + "\n"


def effective_inbreeding(form: F5ScenarioForm, population_id: str) -> str | None:
    """`F` de la subpoblación, con el del escenario como valor de respaldo."""
    for row in form.subpopulations:
        if row.population_id == population_id:
            return row.inbreeding_f if row.inbreeding_f is not None else form.inbreeding_f
    return form.inbreeding_f


def render_f5_scenario_toml(form: F5ScenarioForm) -> str:
    if not form.subpopulations:
        raise EditorError(
            "scenario.subpopulations", "se exige al menos una subpoblación"
        )
    # REQ-UI-15: `controlled` **se desactiva** cuando F es mayor que cero. No
    # basta anunciarlo: un documento con esa combinación no llega a existir, y
    # Rust lo rechaza por su cuenta —en la subpoblación y por unidad—, de modo
    # que la interfaz no inventa ninguna regla propia.
    for index, row in enumerate(form.subpopulations):
        if row.sampling == "controlled" and not controlled_is_available(
            effective_inbreeding(form, row.population_id)
        ):
            raise EditorError(
                f"scenario.subpopulations[{index}].sampling",
                "controlled se desactiva cuando F es mayor que cero",
            )
    for index, row in enumerate(form.sampling_overrides):
        if row.mode == "controlled" and not controlled_is_available(
            effective_inbreeding(form, row.population_id)
        ):
            raise EditorError(
                f"scenario.sampling_overrides[{index}].mode",
                "controlled se desactiva cuando F es mayor que cero",
            )
    lines = [
        f'config_kind = "{F5_SCENARIO_CONFIG_KIND}"',
        f'schema_version = "{F5_SCENARIO_SCHEMA_VERSION}"',
        f'scenario_id = "{_identifier(form.scenario_id, "scenario.scenario_id")}"',
        "founder_count = "
        + str(
            _bounded_integer(
                form.founder_count, "scenario.founder_count", "founder_count"
            )
        ),
        "",
        "[f5]",
        f'profile_metadata_path = "{_relative_path(form.profile_metadata_path, "scenario.f5.profile_metadata_path")}"',
        f'profile_metadata_sha256 = "{_sha256(form.profile_metadata_sha256, "scenario.f5.profile_metadata_sha256")}"',
    ]
    if (form.excluded_regions_path is None) != (form.excluded_regions_sha256 is None):
        raise EditorError(
            "scenario.f5.excluded_regions_path,scenario.f5.excluded_regions_sha256",
            "la máscara exige ruta y hash juntos",
        )
    if form.excluded_regions_path is not None:
        lines.append(
            f'excluded_regions_path = "{_relative_path(form.excluded_regions_path, "scenario.f5.excluded_regions_path")}"'
        )
        lines.append(
            f'excluded_regions_sha256 = "{_sha256(form.excluded_regions_sha256, "scenario.f5.excluded_regions_sha256")}"'
        )
    if form.inbreeding_f is not None:
        lines.append(
            "inbreeding_f = "
            + _unit_interval_literal(form.inbreeding_f, "scenario.f5.inbreeding_f")
        )
    for name in F5_BUDGET_FIELDS:
        value = getattr(form, name)
        if value is not None:
            lines.append(
                f"{name} = {_bounded_integer(value, f'scenario.f5.{name}', name)}"
            )
    lines.append(_render_target(form))
    for index, row in enumerate(form.subpopulations):
        lines.extend(_render_subpopulation(row, index))
    for index, row in enumerate(form.sampling_overrides):
        lines.extend(_render_override(row, index))
    return "\n".join(lines) + "\n"


def _render_target(form: F5ScenarioForm) -> str:
    if form.target_kind not in F5_TARGET_KINDS:
        raise EditorError(
            "scenario.f5.target", f"objetivo no admitido: {form.target_kind}"
        )
    path = f"scenario.f5.{form.target_kind}"
    if form.target_kind == "coverage_fraction":
        literal = _literal(form.target_value, path)
        value = float(literal)
        if not 0 < value <= 1:
            raise EditorError(path, "debe pertenecer al intervalo (0, 1]")
        return f"{form.target_kind} = {literal}"
    parsed = _bounded_integer(form.target_value, path, form.target_kind)
    return f"{form.target_kind} = {parsed}"


def _render_subpopulation(row: F5SubpopulationRow, index: int) -> list[str]:
    path = f"scenario.subpopulations[{index}]"
    if row.founder_kind not in F5_FOUNDER_KINDS:
        raise EditorError(
            f"{path}.founder_kind", "debe ser founder_ids o founder_keys"
        )
    lines = [
        "",
        "[[subpopulations]]",
        f'population_id = "{_identifier(row.population_id, f"{path}.population_id")}"',
        "individual_count = "
        + str(
            _bounded_integer(
                row.individual_count, f"{path}.individual_count", "individual_count"
            )
        ),
    ]
    if row.sampling is not None:
        lines.append(
            f'sampling = "{_sampling_mode(row.sampling, f"{path}.sampling")}"'
        )
    if row.inbreeding_f is not None:
        lines.append(
            "inbreeding_f = "
            + _unit_interval_literal(row.inbreeding_f, f"{path}.inbreeding_f")
        )
    founders = _identifier_list(row.founders, f"{path}.{row.founder_kind}")
    rendered = ", ".join(f'"{founder}"' for founder in founders)
    lines.append(f"{row.founder_kind} = [{rendered}]")
    return lines


def _render_override(row: F5SamplingOverrideRow, index: int) -> list[str]:
    path = f"scenario.sampling_overrides[{index}]"
    mode = _sampling_mode(row.mode, f"{path}.mode")
    if mode == "probabilistic" and row.controlled_max_abs_error is not None:
        raise EditorError(
            f"{path}.controlled_max_abs_error", "probabilistic no admite tolerancia"
        )
    lines = [
        "",
        "[[sampling_overrides]]",
        f'population_id = "{_identifier(row.population_id, f"{path}.population_id")}"',
        "block_index = "
        + str(_bounded_integer(row.block_index, f"{path}.block_index", "block_index")),
        f'mode = "{mode}"',
    ]
    if row.controlled_max_abs_error is not None:
        lines.append(
            "controlled_max_abs_error = "
            + _unit_interval_literal(
                row.controlled_max_abs_error, f"{path}.controlled_max_abs_error"
            )
        )
    return lines


# ---------------------------------------------------------------------------
# Reglas de presentación de REQ-UI-15
# ---------------------------------------------------------------------------


def inbreeding_is_available(ploidy: int) -> bool:
    """`F` solo se muestra con ploidía 2 (REQ-UI-15)."""
    return ploidy == 2


def controlled_is_available(inbreeding_f: str | None) -> bool:
    """`controlled` se desactiva cuando `F` es mayor que cero (REQ-UI-15)."""
    if inbreeding_f is None:
        return True
    try:
        return float(inbreeding_f) <= 0
    except ValueError:
        return False


def available_sampling_modes(inbreeding_f: str | None) -> tuple[str, ...]:
    if controlled_is_available(inbreeding_f):
        return F5_SAMPLING_MODES
    return ("probabilistic",)


def pair_without_inbreeding(pair: F5Pair) -> F5Pair:
    """Retira `F` cuando la ploidía deja de ser 2, en vez de publicarlo oculto."""
    if inbreeding_is_available(pair.general.ploidy):
        return pair
    return replace(
        pair,
        scenario=replace(
            pair.scenario,
            inbreeding_f=None,
            subpopulations=[
                replace(row, inbreeding_f=None) for row in pair.scenario.subpopulations
            ],
        ),
    )


# ---------------------------------------------------------------------------
# Filas de editor
# ---------------------------------------------------------------------------


def f5_subpopulation_editor_rows(form: F5ScenarioForm) -> list[dict[str, Any]]:
    return [
        {
            "order": index,
            "population_id": row.population_id,
            "individual_count": row.individual_count,
            "founder_kind": row.founder_kind,
            "founders": " ".join(row.founders),
            "sampling": row.sampling or "",
            "inbreeding_f": row.inbreeding_f or "",
        }
        for index, row in enumerate(form.subpopulations, start=1)
    ]


def f5_subpopulation_from_editor_row(
    row: Mapping[str, Any], path: str
) -> F5SubpopulationRow:
    founder_kind = str(row.get("founder_kind", "founder_ids")).strip()
    sampling = str(row.get("sampling", "")).strip()
    inbreeding = str(row.get("inbreeding_f", "")).strip()
    return F5SubpopulationRow(
        population_id=_identifier(row.get("population_id"), f"{path}.population_id"),
        individual_count=_bounded_integer(
            row.get("individual_count"), f"{path}.individual_count", "individual_count"
        ),
        founder_kind=founder_kind,
        founders=str(row.get("founders", "")).split(),
        sampling=sampling or None,
        inbreeding_f=inbreeding or None,
    )


def f5_override_editor_rows(form: F5ScenarioForm) -> list[dict[str, Any]]:
    return [
        {
            "order": index,
            "population_id": row.population_id,
            "block_index": row.block_index,
            "mode": row.mode,
            "controlled_max_abs_error": row.controlled_max_abs_error or "",
        }
        for index, row in enumerate(form.sampling_overrides, start=1)
    ]


def f5_override_from_editor_row(
    row: Mapping[str, Any], path: str
) -> F5SamplingOverrideRow:
    tolerance = str(row.get("controlled_max_abs_error", "")).strip()
    return F5SamplingOverrideRow(
        population_id=_identifier(row.get("population_id"), f"{path}.population_id"),
        block_index=_bounded_integer(
            row.get("block_index"), f"{path}.block_index", "block_index"
        ),
        mode=_sampling_mode(row.get("mode"), f"{path}.mode"),
        controlled_max_abs_error=tolerance or None,
    )


# ---------------------------------------------------------------------------
# Utilidades locales
# ---------------------------------------------------------------------------


def _document(text: str, path: str) -> Mapping[str, Any]:
    try:
        return tomlkit.parse(text)
    except Exception as error:  # noqa: BLE001 - tomlkit no expone un tipo único
        raise EditorError(path, f"TOML inválido: {error}") from error


def _keys(
    table: Mapping[str, Any], *, allowed: set[str], required: set[str], path: str
) -> None:
    unknown = set(table) - allowed
    if unknown:
        raise EditorError(path, f"campo desconocido: {sorted(unknown)[0]}")
    missing = required - set(table)
    if missing:
        raise EditorError(path, f"falta el campo: {sorted(missing)[0]}")


def _table(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise EditorError(path, "se esperaba una tabla")
    return value


def _array(value: Any, path: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise EditorError(path, "se esperaba una lista de tablas")
    return value


def _text(value: Any, path: str) -> str:
    if not isinstance(value, str):
        raise EditorError(path, "se esperaba texto")
    return str(value)


def _match(pattern: re.Pattern[str], value: Any, path: str, detail: str) -> str:
    text = _text(value, path)
    if pattern.fullmatch(text) is None:
        raise EditorError(path, detail)
    return text


def _identifier(value: Any, path: str) -> str:
    return _match(
        _ID, value, path, "debe empezar por alfanumérico y usar [A-Za-z0-9._-] (1..64)"
    )


def _sha256(value: Any, path: str) -> str:
    return _match(_SHA, value, path, "debe contener 64 hexadecimales minúsculos")


def _relative_path(value: Any, path: str) -> str:
    return _match(
        _RELATIVE_PATH,
        value,
        path,
        "debe ser relativa, sin '\\', sin '.' ni '..' y usar [A-Za-z0-9._-]",
    )


def _identifier_list(value: Any, path: str) -> list[str]:
    items = _array(value, path)
    if not items:
        raise EditorError(path, "se exige al menos un identificador")
    founders = [_identifier(item, path) for item in items]
    if len(set(founders)) != len(founders):
        raise EditorError(path, "los identificadores deben ser distintos")
    return founders


def _sampling_mode(value: Any, path: str) -> str:
    text = _text(value, path)
    if text not in F5_SAMPLING_MODES:
        raise EditorError(path, "debe ser probabilistic o controlled")
    return text


def _bounded_integer(value: Any, path: str, bound: str) -> int:
    minimum, maximum = F5_INTEGER_BOUNDS[bound]
    if isinstance(value, bool):
        raise EditorError(path, "se esperaba un entero")
    if isinstance(value, str):
        text = value.strip()
        if not text.isdecimal():
            raise EditorError(path, "se esperaba un entero")
        value = int(text)
    if isinstance(value, float):
        if not value.is_integer():
            raise EditorError(path, "se esperaba un entero sin redondeo")
        value = int(value)
    if not isinstance(value, int):
        raise EditorError(path, "se esperaba un entero")
    if not minimum <= value <= maximum:
        raise EditorError(path, f"debe pertenecer a {minimum}..={maximum}")
    return int(value)


def _literal(value: Any, path: str) -> str:
    """Literal decimal conservado tal y como se escribió."""
    if isinstance(value, (Integer, Float)):
        return value.as_string()
    text = str(value).strip()
    if not text:
        raise EditorError(path, "falta un número decimal")
    try:
        parsed = tomlkit.parse(f"value = {text}\n")["value"]
    except Exception as error:  # noqa: BLE001 - tomlkit no expone un tipo único
        raise EditorError(path, "literal decimal inválido") from error
    if not isinstance(parsed, (Integer, Float)):
        raise EditorError(path, "se esperaba un número decimal")
    return parsed.as_string()


def _number_literal(value: Any, path: str) -> str:
    if not isinstance(value, (Integer, Float)):
        raise EditorError(path, "se esperaba un número decimal")
    return value.as_string()


def _unit_interval(value: Any, path: str) -> str:
    literal = _number_literal(value, path)
    return _check_unit_interval(literal, path)


def _unit_interval_literal(value: Any, path: str) -> str:
    return _check_unit_interval(_literal(value, path), path)


def _check_unit_interval(literal: str, path: str) -> str:
    try:
        parsed = float(literal)
    except ValueError as error:
        raise EditorError(path, "debe pertenecer al intervalo [0, 1]") from error
    if not 0 <= parsed <= 1:
        raise EditorError(path, "debe pertenecer al intervalo [0, 1]")
    return literal


def f5_pair_documents(pair: F5Pair) -> tuple[str, str]:
    """Los dos documentos exportables, en orden de rol.

    Devuelve los bytes importados cuando siguen describiendo exactamente el
    formulario vigente, y el renderizado canónico en cuanto no lo hacen. La
    comprobación es un reparseo del propio origen: no se confía en que nadie
    haya editado el formulario, se verifica.
    """
    return (
        _preserved(pair.general_source, pair.general, parse_f5_general_toml)
        or render_f5_general_toml(pair.general),
        _preserved(pair.scenario_source, pair.scenario, parse_f5_scenario_toml)
        or render_f5_scenario_toml(pair.scenario),
    )


def _preserved(source: str | None, form: Any, parse: Any) -> str | None:
    if source is None:
        return None
    try:
        return source if parse(source) == form else None
    except EditorError:
        return None


def f5_control_names() -> Iterable[str]:
    yield from F5_GENERAL_FIELD_CONTROLS.values()
    yield from F5_SCENARIO_FIELD_CONTROLS.values()
