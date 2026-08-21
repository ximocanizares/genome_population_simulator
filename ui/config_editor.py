"""Modelo de autoría F1 para Marimo, sin semántica biológica propia."""

from __future__ import annotations

from dataclasses import dataclass
import math
import re
from typing import Any, Mapping, Sequence

import tomlkit
from tomlkit.items import Float, Integer

U64_MAX = 18_446_744_073_709_551_615
BASES = ("A", "C", "G", "T")
SEQUENCE_CLASSES = ("chromosome", "scaffold")


class EditorError(ValueError):
    """Error estructural localizable anterior a la validación normativa Rust."""

    def __init__(self, path: str, detail: str) -> None:
        super().__init__(f"{path}: {detail}")
        self.path = path
        self.detail = detail


@dataclass(eq=True)
class SequenceRow:
    id: str
    sequence_class: str
    length: int


@dataclass(eq=True)
class ReferenceForm:
    kind: str
    id: str
    fasta: str | None
    composition: dict[str, str] | None
    sequences: list[SequenceRow]


@dataclass(eq=True)
class ReductionRow:
    id: str
    source_reference_id: str
    region_size_bases: int
    target_kind: str
    target_literal: str


@dataclass(eq=True)
class PopulationSnvsForm:
    mode: str
    input_path: str | None
    metadata_path: str | None


@dataclass(eq=True)
class GeneralForm:
    seed: str | None
    ploidy: int
    normalization_tolerance: str
    reference: ReferenceForm
    reduced_references: list[ReductionRow]
    population_snvs: PopulationSnvsForm | None


def default_form() -> GeneralForm:
    return GeneralForm(
        seed=None,
        ploidy=2,
        normalization_tolerance="1e-9",
        reference=ReferenceForm(
            kind="simulated",
            id="reference-main",
            fasta=None,
            composition={base: "0.25" for base in BASES},
            sequences=[SequenceRow("chr1", "chromosome", 1_000)],
        ),
        reduced_references=[],
        population_snvs=None,
    )


def editor_integer(value: Any, path: str) -> int:
    """Accept only values that represent an integer exactly, without truncation."""
    if isinstance(value, bool):
        raise EditorError(path, "se esperaba un entero")
    if isinstance(value, int):
        return value
    if isinstance(value, float) and math.isfinite(value) and value.is_integer():
        return int(value)
    if isinstance(value, str):
        text = value.strip()
        if text and (text.isdecimal() or (text[0] == "-" and text[1:].isdecimal())):
            return int(text)
    raise EditorError(path, "se esperaba un entero sin redondeo ni truncamiento")


def ordered_editor_rows(rows: Sequence[Mapping[str, Any]], path: str) -> list[Mapping[str, Any]]:
    """Return rows in explicit positive order, rejecting ambiguous positions."""
    positioned: list[tuple[int, Mapping[str, Any]]] = []
    seen: set[int] = set()
    for index, row in enumerate(rows):
        order = editor_integer(row.get("order"), f"{path}[{index}].order")
        if order <= 0:
            raise EditorError(f"{path}[{index}].order", "debe ser positivo")
        if order in seen:
            raise EditorError(path, f"order repetido: {order}")
        seen.add(order)
        positioned.append((order, row))
    positioned.sort(key=lambda positioned_row: positioned_row[0])
    return [row for _, row in positioned]


def parse_general_toml(text: str) -> GeneralForm:
    try:
        document = tomlkit.parse(text)
    except Exception as error:
        raise EditorError("document", f"TOML inválido: {error}") from error

    _keys(
        document,
        allowed={
            "config_kind",
            "schema_version",
            "seed",
            "ploidy",
            "normalization_tolerance",
            "reference",
            "reduced_references",
            "population_snvs",
        },
        required={"config_kind", "schema_version", "reference"},
        path="document",
    )
    if _string(document["config_kind"], "config_kind") != "general":
        raise EditorError("config_kind", 'debe ser "general"')
    if _string(document["schema_version"], "schema_version") != "1":
        raise EditorError("schema_version", 'debe ser "1"')
    seed = _seed(document["seed"]) if "seed" in document else None
    ploidy = _ploidy(document["ploidy"]) if "ploidy" in document else 2
    normalization_tolerance = (
        _number_literal(document["normalization_tolerance"], "normalization_tolerance")
        if "normalization_tolerance" in document
        else "1e-9"
    )
    reference = _parse_reference(
        _mapping(document["reference"], "reference"), "reference"
    )
    reductions = _parse_reductions(document.get("reduced_references", []))
    population_snvs = (
        _parse_population_snvs(
            _mapping(document["population_snvs"], "population_snvs")
        )
        if "population_snvs" in document
        else None
    )
    return GeneralForm(
        seed, ploidy, normalization_tolerance, reference, reductions, population_snvs
    )


def _parse_population_snvs(table: Mapping[str, Any]) -> PopulationSnvsForm:
    _keys(
        table,
        allowed={"mode", "input_path", "metadata_path"},
        required={"mode"},
        path="population_snvs",
    )
    mode = _string(table["mode"], "population_snvs.mode")
    input_path = table.get("input_path")
    metadata_path = table.get("metadata_path")
    if mode == "supplied":
        if metadata_path is not None:
            raise EditorError("population_snvs.metadata_path", "no admite supplied")
        return PopulationSnvsForm(
            mode, _non_empty_string(input_path, "population_snvs.input_path"), None
        )
    if mode == "reuse":
        if input_path is not None:
            raise EditorError("population_snvs.input_path", "no admite reuse")
        return PopulationSnvsForm(
            mode, None, _non_empty_string(metadata_path, "population_snvs.metadata_path")
        )
    if mode == "generate":
        if input_path is not None or metadata_path is not None:
            raise EditorError("population_snvs", "generate no admite rutas")
        return PopulationSnvsForm(mode, None, None)
    raise EditorError("population_snvs.mode", "debe ser supplied, reuse o generate")


def _parse_reference(table: Mapping[str, Any], path: str) -> ReferenceForm:
    kind = _string(table.get("kind"), f"{path}.kind")
    common = {"kind", "id", "sequences"}
    if kind == "simulated":
        _keys(
            table,
            allowed=common | {"composition"},
            required=common | {"composition"},
            path=path,
        )
        composition_table = _mapping(table["composition"], f"{path}.composition")
        _keys(
            composition_table,
            allowed=set(BASES),
            required=set(BASES),
            path=f"{path}.composition",
        )
        composition = {
            base: _number_literal(
                composition_table[base], f"{path}.composition.{base}"
            )
            for base in BASES
        }
        fasta = None
    elif kind == "supplied":
        _keys(
            table,
            allowed=common | {"fasta"},
            required=common | {"fasta"},
            path=path,
        )
        composition = None
        fasta = _string(table["fasta"], f"{path}.fasta")
    else:
        raise EditorError(f"{path}.kind", "debe ser simulated o supplied")

    sequences_value = _sequence(table["sequences"], f"{path}.sequences")
    if not sequences_value:
        raise EditorError(f"{path}.sequences", "se exige al menos una fila")
    sequences: list[SequenceRow] = []
    for index, value in enumerate(sequences_value):
        row_path = f"{path}.sequences[{index}]"
        row = _mapping(value, row_path)
        _keys(
            row,
            allowed={"id", "class", "length"},
            required={"id", "class", "length"},
            path=row_path,
        )
        sequence_class = _string(row["class"], f"{row_path}.class")
        if sequence_class not in SEQUENCE_CLASSES:
            raise EditorError(
                f"{row_path}.class", "debe ser chromosome o scaffold"
            )
        sequences.append(
            SequenceRow(
                _string(row["id"], f"{row_path}.id"),
                sequence_class,
                _positive_u64(row["length"], f"{row_path}.length"),
            )
        )
    return ReferenceForm(
        kind=kind,
        id=_string(table["id"], f"{path}.id"),
        fasta=fasta,
        composition=composition,
        sequences=sequences,
    )


def _parse_reductions(value: Any) -> list[ReductionRow]:
    rows = _sequence(value, "reduced_references")
    reductions: list[ReductionRow] = []
    for index, entry in enumerate(rows):
        path = f"reduced_references[{index}]"
        row = _mapping(entry, path)
        allowed = {
            "id",
            "source_reference_id",
            "region_size_bases",
            "target_bases",
            "reduction_percent",
        }
        required = {"id", "source_reference_id", "region_size_bases"}
        _keys(row, allowed=allowed, required=required, path=path)
        targets = [name for name in ("target_bases", "reduction_percent") if name in row]
        if len(targets) != 1:
            raise EditorError(
                path, "se exige exactamente target_bases o reduction_percent"
            )
        target_kind = targets[0]
        target_literal = (
            str(_positive_u64(row[target_kind], f"{path}.{target_kind}"))
            if target_kind == "target_bases"
            else _number_literal(row[target_kind], f"{path}.{target_kind}")
        )
        reductions.append(
            ReductionRow(
                id=_string(row["id"], f"{path}.id"),
                source_reference_id=_string(
                    row["source_reference_id"], f"{path}.source_reference_id"
                ),
                region_size_bases=_positive_u64(
                    row["region_size_bases"], f"{path}.region_size_bases"
                ),
                target_kind=target_kind,
                target_literal=target_literal,
            )
        )
    return reductions


def render_general_toml(form: GeneralForm) -> str:
    document = tomlkit.document()
    document.add("config_kind", "general")
    document.add("schema_version", "1")
    if form.seed is not None:
        document.add("seed", _seed(form.seed))
    document.add("ploidy", _ploidy(form.ploidy))
    document.add(
        "normalization_tolerance",
        _literal_item(form.normalization_tolerance, "normalization_tolerance"),
    )

    reference = tomlkit.table()
    reference.add("kind", form.reference.kind)
    reference.add("id", form.reference.id)
    if form.reference.kind == "simulated":
        if form.reference.composition is None:
            raise EditorError("reference.composition", "falta para simulated")
        composition = tomlkit.table()
        for base in BASES:
            composition.add(
                base,
                _literal_item(
                    form.reference.composition.get(base),
                    f"reference.composition.{base}",
                ),
            )
        reference.add("composition", composition)
    elif form.reference.kind == "supplied":
        if not form.reference.fasta:
            raise EditorError("reference.fasta", "falta para supplied")
        reference.add("fasta", form.reference.fasta)
    else:
        raise EditorError("reference.kind", "debe ser simulated o supplied")

    if not form.reference.sequences:
        raise EditorError("reference.sequences", "se exige al menos una fila")
    sequences = tomlkit.aot()
    for index, sequence in enumerate(form.reference.sequences):
        row = tomlkit.table()
        row.add("id", sequence.id)
        if sequence.sequence_class not in SEQUENCE_CLASSES:
            raise EditorError(
                f"reference.sequences[{index}].class",
                "debe ser chromosome o scaffold",
            )
        row.add("class", sequence.sequence_class)
        row.add(
            "length",
            _checked_positive_u64(
                sequence.length, f"reference.sequences[{index}].length"
            ),
        )
        sequences.append(row)
    reference.add("sequences", sequences)
    document.add("reference", reference)

    if form.reduced_references:
        reductions = tomlkit.aot()
        for index, reduction in enumerate(form.reduced_references):
            path = f"reduced_references[{index}]"
            row = tomlkit.table()
            row.add("id", reduction.id)
            row.add("source_reference_id", reduction.source_reference_id)
            row.add(
                "region_size_bases",
                _checked_positive_u64(
                    reduction.region_size_bases, f"{path}.region_size_bases"
                ),
            )
            if reduction.target_kind == "target_bases":
                row.add(
                    "target_bases",
                    _checked_positive_u64(
                        int(reduction.target_literal), f"{path}.target_bases"
                    ),
                )
            elif reduction.target_kind == "reduction_percent":
                row.add(
                    "reduction_percent",
                    _literal_item(
                        reduction.target_literal, f"{path}.reduction_percent"
                    ),
                )
            else:
                raise EditorError(
                    f"{path}.target_kind",
                    "debe ser target_bases o reduction_percent",
                )
            reductions.append(row)
        document.add("reduced_references", reductions)
    if form.population_snvs is not None:
        population_snvs = form.population_snvs
        table = tomlkit.table()
        if population_snvs.mode == "supplied":
            table.add("mode", "supplied")
            table.add(
                "input_path",
                _non_empty_string(
                    population_snvs.input_path, "population_snvs.input_path"
                ),
            )
        elif population_snvs.mode == "reuse":
            table.add("mode", "reuse")
            table.add(
                "metadata_path",
                _non_empty_string(
                    population_snvs.metadata_path, "population_snvs.metadata_path"
                ),
            )
        elif population_snvs.mode == "generate":
            if population_snvs.input_path is not None or population_snvs.metadata_path is not None:
                raise EditorError("population_snvs", "generate no admite rutas")
            table.add("mode", "generate")
        else:
            raise EditorError(
                "population_snvs.mode", "debe ser supplied, reuse o generate"
            )
        document.add("population_snvs", table)
    return tomlkit.dumps(document)


def _keys(
    table: Mapping[str, Any],
    *,
    allowed: set[str],
    required: set[str],
    path: str,
) -> None:
    unknown = set(table) - allowed
    if unknown:
        raise EditorError(path, f"campo desconocido: {sorted(unknown)[0]}")
    missing = required - set(table)
    if missing:
        raise EditorError(path, f"falta el campo: {sorted(missing)[0]}")


def _mapping(value: Any, path: str) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise EditorError(path, "se esperaba una tabla")
    return value


def _sequence(value: Any, path: str) -> Sequence[Any]:
    if not isinstance(value, Sequence) or isinstance(value, (str, bytes)):
        raise EditorError(path, "se esperaba una lista de tablas")
    return value


def _string(value: Any, path: str) -> str:
    if not isinstance(value, str):
        raise EditorError(path, "se esperaba texto")
    return value


def _non_empty_string(value: Any, path: str) -> str:
    text = _string(value, path)
    if not text:
        raise EditorError(path, "no puede estar vacío")
    return text


def _seed(value: Any) -> str:
    seed = _string(value, "seed")
    if re.fullmatch(r"[0-9a-f]{64}", seed) is None:
        raise EditorError("seed", "debe contener 64 hexadecimales minúsculos")
    return seed


def _positive_u64(value: Any, path: str) -> int:
    if not isinstance(value, Integer):
        raise EditorError(path, "se esperaba un entero")
    return _checked_positive_u64(int(value), path)


def _ploidy(value: Any) -> int:
    if isinstance(value, bool) or not isinstance(value, (Integer, int)):
        raise EditorError("ploidy", "se esperaba un entero")
    return int(value)


def _checked_positive_u64(value: int, path: str) -> int:
    if isinstance(value, bool) or not 1 <= value <= U64_MAX:
        raise EditorError(path, "debe pertenecer a 1..=u64::MAX")
    return value


def _number_literal(value: Any, path: str) -> str:
    if not isinstance(value, (Integer, Float)):
        raise EditorError(path, "se esperaba un número decimal")
    return value.as_string()


def _literal_item(literal: str | None, path: str) -> Integer | Float:
    if not literal:
        raise EditorError(path, "falta un número decimal")
    try:
        value = tomlkit.parse(f"value = {literal}\n")["value"]
    except Exception as error:
        raise EditorError(path, "literal decimal inválido") from error
    if not isinstance(value, (Integer, Float)):
        raise EditorError(path, "se esperaba un número decimal")
    return value
