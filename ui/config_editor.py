"""Modelo de autoría F2 para Marimo, sin semántica biológica propia."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
import re
from typing import Any, Mapping, Sequence

import tomlkit
from tomlkit.items import Float, Integer

U64_MAX = 18_446_744_073_709_551_615
BASES = ("A", "C", "G", "T")
COMPOSITION_BASES = (*BASES, "N")
SEQUENCE_CLASSES = ("chromosome", "scaffold")
COMPOSITION_TOLERANCE = 1e-9
REDUCTION_OBJECTIVE_REMOVE_PERCENT = "Eliminar porcentaje"
REDUCTION_OBJECTIVE_KEEP_BASES = "Conservar bases"


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


def read_fasta_sequence_rows(path: Path) -> list[SequenceRow]:
    """Read FASTA identifiers and lengths for authoring, without modifying the input.

    Classes are intentionally initialized as ``scaffold``: a FASTA header does
    not reliably distinguish chromosomes from unplaced or auxiliary sequences.
    The Rust preflight remains the authority for the complete FASTA contract.
    """
    rows: list[SequenceRow] = []
    current_id: str | None = None
    current_length = 0

    def finish_current() -> None:
        nonlocal current_id, current_length
        if current_id is None:
            return
        if current_length == 0:
            raise EditorError(
                "reference.fasta", f"la secuencia {current_id} está vacía"
            )
        rows.append(SequenceRow(current_id, "scaffold", current_length))

    try:
        with path.open("rb") as source:
            for line_number, raw_line in enumerate(source, start=1):
                line = raw_line.strip()
                if not line:
                    continue
                if line.startswith(b">"):
                    finish_current()
                    identifier = line[1:].split(maxsplit=1)[0]
                    if not identifier:
                        raise EditorError(
                            "reference.fasta", f"cabecera vacía en línea {line_number}"
                        )
                    try:
                        current_id = identifier.decode("ascii")
                    except UnicodeDecodeError as error:
                        raise EditorError(
                            "reference.fasta",
                            f"ID FASTA no ASCII en línea {line_number}",
                        ) from error
                    if any(row.id == current_id for row in rows):
                        raise EditorError(
                            "reference.fasta", f"ID FASTA repetido: {current_id}"
                        )
                    current_length = 0
                else:
                    if current_id is None:
                        raise EditorError(
                            "reference.fasta",
                            f"secuencia sin cabecera en línea {line_number}",
                        )
                    current_length += len(line)
    except OSError as error:
        raise EditorError("reference.fasta", f"no se puede leer: {error}") from error

    finish_current()
    if not rows:
        raise EditorError("reference.fasta", "no contiene secuencias FASTA")
    if len({row.id for row in rows}) != len(rows):
        raise EditorError("reference.fasta", "contiene IDs FASTA repetidos")
    return rows


@dataclass(eq=True)
class RepeatInterruptionRow:
    """Una interrupción declarada: dónde empieza y qué bases la forman.

    Se declara y no se deriva: dadas unas bases sueltas no hay forma de saber
    dónde el usuario quería la impureza. Por eso la interfaz la transporta tal
    cual y no la reconstruye a partir del número de unidades.
    """

    offset: int
    bases: str


@dataclass(eq=True)
class RepeatAlleleRow:
    """Un alelo del locus: su tamaño, su frecuencia y sus interrupciones.

    El tamaño se declara de una de dos formas excluyentes: ``units``, con un
    valor por bloque del locus, o ``sequence``, con las bases exactas. La
    pureza **no aparece**: es siempre un valor calculado sobre la secuencia y
    declararla está prohibido.
    """

    frequency: str
    units: list[int] | None = None
    sequence: str | None = None
    interruptions: list[RepeatInterruptionRow] = field(default_factory=list)


@dataclass(eq=True)
class RepeatLocusRow:
    """Un locus de repetición: dónde está, con qué motivos y con qué alelos."""

    id: str
    sequence: str
    start: int
    end: int
    motifs: list[str]
    alleles: list[RepeatAlleleRow]
    max_units: int | None = None
    min_purity: str | None = None
    min_longest_run_units: int | None = None


@dataclass(eq=True)
class RepeatAdmissionForm:
    """Mínimos generales de la guarda de admisión."""

    min_purity: str | None = None
    min_longest_run_units: int | None = None


@dataclass(eq=True)
class ReferenceForm:
    kind: str
    id: str
    fasta: str | None
    composition: dict[str, str] | None
    sequences: list[SequenceRow]
    repeat_loci: list[RepeatLocusRow] = field(default_factory=list)
    repeat_admission: RepeatAdmissionForm | None = None
    max_tract_output_bases: int | None = None


@dataclass(eq=True)
class ReductionRow:
    id: str
    source_reference_id: str
    region_size_bases: int
    target_kind: str
    target_literal: str


def reduction_editor_row(row: ReductionRow, order: int) -> dict[str, Any]:
    """Present one normative reduction row with author-facing column names."""
    objectives = {
        "target_bases": REDUCTION_OBJECTIVE_KEEP_BASES,
        "reduction_percent": REDUCTION_OBJECTIVE_REMOVE_PERCENT,
    }
    try:
        objective = objectives[row.target_kind]
    except KeyError as error:
        raise EditorError(
            f"reduced_references[{order - 1}].target_kind",
            "debe ser target_bases o reduction_percent",
        ) from error
    return {
        "order": order,
        "id": row.id,
        "source_reference_id": row.source_reference_id,
        "region_size_bases": row.region_size_bases,
        "Objetivo": objective,
        "Valor": row.target_literal,
    }


def reduction_from_editor_row(row: Mapping[str, Any], path: str) -> ReductionRow:
    """Translate the readable reduction choice back to the normative TOML branch."""
    targets = {
        REDUCTION_OBJECTIVE_KEEP_BASES: "target_bases",
        REDUCTION_OBJECTIVE_REMOVE_PERCENT: "reduction_percent",
    }
    objective = str(row.get("Objetivo", ""))
    try:
        target_kind = targets[objective]
    except KeyError as error:
        raise EditorError(
            f"{path}.Objetivo",
            "debe ser Eliminar porcentaje o Conservar bases",
        ) from error
    return ReductionRow(
        id=str(row["id"]),
        source_reference_id=str(row["source_reference_id"]),
        region_size_bases=editor_integer(
            row["region_size_bases"], f"{path}.region_size_bases"
        ),
        target_kind=target_kind,
        target_literal=str(row.get("Valor", "")),
    )


@dataclass(eq=True)
class PopulationSnvsForm:
    mode: str
    input_path: str | None
    metadata_path: str | None
    generate: dict[str, Any] | None = None


@dataclass(eq=True)
class InitialPopulationRow:
    population_id: str
    individual_count: int
    ploidy: int
    genotype_model: str
    inbreeding_coefficient: str | None
    individual_ids: list[str] | None
    founder_keys: list[str] | None


@dataclass(eq=True)
class DirectEventRow:
    event_id: str
    individual_id: str
    phase: int
    kind: str
    contig: str
    position_1based: int | None = None
    reference: str | None = None
    alternate: str | None = None
    point: int | None = None
    bases: str | None = None
    start: int | None = None
    end: int | None = None


@dataclass(eq=True)
class TechnicalLibraryRow:
    individual_id: str
    i7: str
    i5: str


# Los cuatro pasos admitidos, en orden canónico. **La estructura no se
# configura**: la fija el algoritmo. La interfaz no ofrece ningún control para
# un paso `+2` ni para un modelo distinto, porque declararlo se rechaza.
STUTTER_ADMITTED_STEPS = (-3, -2, -1, 1)
# Tramo de longitud de motivo que la calibración cubre, en pares de base.
STUTTER_MOTIF_LENGTHS = (1, 2, 3, 4, 5, 6)
# Las cifras calibradas de la columna sin PCR, **las mismas que Rust lleva
# dentro**. Aquí solo sirven para que la tabla del cuaderno aparezca rellena con
# lo que el simulador usaría de todos modos; la autoridad sigue siendo Rust, y
# `tests_ui/test_stutter_ui.py` comprueba que las dos copias coinciden.
STUTTER_CALIBRATED_MOTIF_LENGTH_RATES = {
    1: 0.027,
    2: 0.011,
    3: 0.0029,
    4: 0.0016,
    5: 0.0015,
    6: 0.0010,
}


# Constantes de contrato de la amplificación, **las mismas que Rust lleva
# dentro**. Aquí solo acotan los controles del cuaderno; la autoridad sigue
# siendo Rust, y `tests/pcr_contract.rs` compara las dos copias junto con las
# de los tres esquemas. Una copia decimal silenciosa en un solo sitio es lo que
# sobrevivió al cambio de contrato la vez anterior.
PCR_MAXIMUM_CYCLES = 30
PCR_MAXIMUM_EARLY_CYCLES = 6
PCR_MAXIMUM_FAMILY_SIZE = 256
PCR_MAXIMUM_GC_ATTEMPTS = 16
PCR_MAXIMUM_EFFICIENCY_CLASSES = 4
PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE = 64.0
# Las dos preparaciones admitidas, en orden canónico. Ausente significa la
# primera, que es la única que existía.
LIBRARY_PREPARATIONS = ("pcr_free", "pcr")


@dataclass(eq=True)
class StutterMotifRateRow:
    motif_length_bases: int
    rate: str


@dataclass(eq=True)
class StutterStepRow:
    units: int
    weight: str


@dataclass(eq=True)
class StutterLocusRateRow:
    locus_id: str
    rate: str


@dataclass(eq=True)
class StutterForm:
    """Tartamudeo de base declarado desde la interfaz.

    Los **pesos de los cuatro pasos** son obligatorios cuando el tartamudeo
    está activo y no tienen valor por omisión: la fuente publica el método de
    ajuste, no los valores. La interfaz no los rellena sola.
    """

    enabled: bool = False
    motif_length_rates: list[StutterMotifRateRow] = field(default_factory=list)
    run_length_coefficient: str | None = None
    steps: list[StutterStepRow] = field(default_factory=list)
    locus_rates: list[StutterLocusRateRow] = field(default_factory=list)


@dataclass(eq=True)
class GcCurveKnotRow:
    gc_fraction: str
    relative_efficiency: str


@dataclass(eq=True)
class EfficiencyClassRow:
    efficiency: str
    weight: str


@dataclass(eq=True)
class PcrForm:
    """Preparación amplificada declarada desde la interfaz.

    Los tres primeros campos son **obligatorios y sin valor por omisión**: la
    fracción de duplicados depende de la masa de ADN de entrada, del número de
    ciclos y del laboratorio. La interfaz no los rellena sola, y la mejor
    fuente para el primero son los informes de calidad del propio laboratorio.

    Y lo que se declara es la fracción **esperada de una biblioteca
    homogénea**, no la realizada: la realizada es mayor, y el manifiesto la
    publica aparte.
    """

    duplicate_fraction: str
    cycles: int
    efficiency: str
    early_cycles: int | None = None
    max_family_size: int | None = None
    gc_max_attempts: int | None = None
    gc_relative_efficiency: list[GcCurveKnotRow] = field(default_factory=list)
    efficiency_classes: list[EfficiencyClassRow] = field(default_factory=list)


@dataclass(eq=True)
class TechnicalEventForm:
    event_id: str
    individual_ids: list[str]
    coverage_wgs: str
    libraries: list[TechnicalLibraryRow]
    read_length_bases: int | None = None
    stutter: StutterForm | None = None
    # `None` conserva la omisión normativa: la preparación sin amplificar.
    library_preparation: str | None = None
    pcr: PcrForm | None = None


@dataclass(eq=True)
class GeneralForm:
    seed: str | None
    ploidy: int
    normalization_tolerance: str
    reference: ReferenceForm
    reduced_references: list[ReductionRow]
    population_snvs: PopulationSnvsForm | None
    initial_populations: list[InitialPopulationRow]
    direct_events: list[DirectEventRow]
    technical_event: TechnicalEventForm | None = None
    # `None` conserva la omisión normativa: no se publica una vista VCF.
    vcf_output_mode: str | None = None


def default_form() -> GeneralForm:
    return GeneralForm(
        seed=None,
        ploidy=2,
        normalization_tolerance="1e-9",
        reference=ReferenceForm(
            kind="simulated",
            id="reference-main",
            fasta=None,
            composition={**{base: "0.25" for base in BASES}, "N": "0"},
            sequences=[SequenceRow("chr1", "chromosome", 1_000)],
        ),
        reduced_references=[],
        population_snvs=None,
        initial_populations=[],
        direct_events=[],
        technical_event=None,
        vcf_output_mode=None,
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


def population_snvs_population_from_editor_row(
    row: Mapping[str, Any], components_by_population: Mapping[str, Sequence[Mapping[str, Any]]]
) -> dict[str, Any]:
    """Shape a population editor row without assigning scientific meaning in Python."""
    population_id = str(row["population_id"]).strip()
    frequency_model = str(row.get("frequency_model", "uniform"))
    population: dict[str, Any] = {
        "population_id": population_id,
        "frequency_model": frequency_model,
    }
    if frequency_model == "uniform":
        population["min_alt_af"] = float(row["min_alt_af"])
        population["max_alt_af"] = float(row["max_alt_af"])
    else:
        population["frequency_mode"] = str(row.get("frequency_mode", "stochastic"))
        if frequency_model == "spectrum":
            population["spectrum"] = list(components_by_population.get(population_id, []))
    return population


def initial_population_from_editor_row(row: Mapping[str, Any], path: str) -> InitialPopulationRow:
    """Turn one editor row into TOML-shaped data without biological validation."""
    genotype_model = str(row.get("genotype_model", "hwe"))
    coefficient = str(row.get("inbreeding_coefficient", "")).strip()
    ids_text = str(row.get("individual_ids", "")).strip()
    individual_ids = [item.strip() for item in ids_text.split(",")] if ids_text else None
    founder_keys_text = str(row.get("founder_keys", "")).strip()
    founder_keys = [item.strip() for item in founder_keys_text.split(",")] if founder_keys_text else None
    if individual_ids is not None and any(not item for item in individual_ids):
        raise EditorError(f"{path}.individual_ids", "contiene un ID vacío")
    if founder_keys is not None and any(not item for item in founder_keys):
        raise EditorError(f"{path}.founder_keys", "contiene una clave vacía")
    if individual_ids is not None and founder_keys is not None:
        raise EditorError(path, "individual_ids y founder_keys son mutuamente excluyentes")
    if individual_ids is None and founder_keys is None:
        raise EditorError(f"{path}.founder_keys", "es obligatorio si se omite individual_ids")
    return InitialPopulationRow(
        population_id=str(row.get("population_id", "")).strip(),
        individual_count=editor_integer(row.get("individual_count"), f"{path}.individual_count"),
        ploidy=editor_integer(row.get("ploidy"), f"{path}.ploidy"),
        genotype_model=genotype_model,
        inbreeding_coefficient=coefficient or None,
        individual_ids=individual_ids,
        founder_keys=founder_keys,
    )


def direct_event_from_editor_row(row: Mapping[str, Any], path: str) -> DirectEventRow:
    """Turn one explicit-event editor row into its selected TOML branch."""
    kind = str(row.get("kind", ""))
    common = {
        "event_id": str(row.get("event_id", "")).strip(),
        "individual_id": str(row.get("individual_id", "")).strip(),
        "phase": editor_integer(row.get("phase"), f"{path}.phase"),
        "kind": kind,
        "contig": str(row.get("contig", "")).strip(),
    }
    if kind == "snp":
        return DirectEventRow(
            **common,
            position_1based=editor_integer(row.get("position_1based"), f"{path}.position_1based"),
            reference=str(row.get("ref", "")).strip(),
            alternate=str(row.get("alt", "")).strip(),
        )
    if kind == "insertion":
        return DirectEventRow(
            **common,
            point=editor_integer(row.get("point"), f"{path}.point"),
            bases=str(row.get("bases", "")).strip(),
        )
    if kind == "deletion":
        return DirectEventRow(
            **common,
            start=editor_integer(row.get("start"), f"{path}.start"),
            end=editor_integer(row.get("end"), f"{path}.end"),
        )
    raise EditorError(f"{path}.kind", "debe ser snp, insertion o deletion")


# ---------------------------------------------------------------------------
# Filas de los editores de repeticiones en tándem
# ---------------------------------------------------------------------------
#
# La interfaz de Marimo edita tablas planas, de modo que un locus, sus alelos y
# las interrupciones de cada alelo viven en tres tablas enlazadas por
# identificador y por orden, igual que ya hacían las poblaciones y sus
# espectros. Lo que **no** hay aquí es ninguna fila para la pureza objetivo ni
# para generar loci al azar: ni una ni otra pertenecen a este corte, y una
# columna que las pidiera sería una promesa que el simulador no cumple.


def stutter_motif_rate_editor_rows(
    stutter: StutterForm | None,
) -> list[dict[str, Any]]:
    """Las seis longitudes de motivo calibradas, **siempre las seis**.

    La tabla no crece ni encoge: la calibración cubre de 1 a 6 pares de base y
    ninguna otra longitud tiene tasa. Presentarlas todas, con la calibrada ya
    puesta, evita que alguien declare una fila de 7 pb que el simulador
    rechazará después.
    """

    declared = {row.motif_length_bases: row.rate for row in (stutter.motif_length_rates if stutter else [])}
    return [
        {
            "order": order,
            "motif_length_bases": length,
            "rate": declared.get(length, str(STUTTER_CALIBRATED_MOTIF_LENGTH_RATES[length])),
        }
        for order, length in enumerate(STUTTER_MOTIF_LENGTHS, start=1)
    ]


def stutter_step_editor_rows(stutter: StutterForm | None) -> list[dict[str, Any]]:
    """Los cuatro pasos admitidos, **siempre los cuatro y solo esos**.

    La estructura la fija el algoritmo, así que la tabla no admite añadir ni
    quitar filas: lo que el usuario declara son los pesos. No hay fila para un
    paso `+2` porque declararlo se rechaza.

    El peso llega **vacío** cuando no se ha declarado, y no con un reparto
    sugerido: los pesos son obligatorios y no tienen valor por omisión, y
    rellenarlos aquí sería inventar una calibración que nadie ha hecho.
    """

    declared = {row.units: row.weight for row in (stutter.steps if stutter else [])}
    return [
        {"order": order, "units": units, "weight": declared.get(units, "")}
        for order, units in enumerate(STUTTER_ADMITTED_STEPS, start=1)
    ]


def stutter_locus_rate_editor_rows(
    stutter: StutterForm | None,
) -> list[dict[str, Any]]:
    rows = [
        {"order": order, "locus_id": row.locus_id, "rate": row.rate}
        for order, row in enumerate(stutter.locus_rates if stutter else [], start=1)
    ]
    return rows or [{"order": 1, "locus_id": "", "rate": ""}]


def gc_curve_editor_rows(pcr: PcrForm | None) -> list[dict[str, Any]]:
    """Nudos de la curva de composición.

    Arranca con los **dos nudos neutros**, que es la curva que no distingue
    ninguna composición y con la que la aceptación no consume aleatoriedad.
    Presentarlos evita que alguien declare un solo nudo, que el simulador
    rechaza porque una curva necesita al menos dos.
    """

    rows = [
        {
            "order": order,
            "gc_fraction": row.gc_fraction,
            "relative_efficiency": row.relative_efficiency,
        }
        for order, row in enumerate(pcr.gc_relative_efficiency if pcr else [], start=1)
    ]
    return rows or [
        {"order": 1, "gc_fraction": "0.0", "relative_efficiency": "1.0"},
        {"order": 2, "gc_fraction": "1.0", "relative_efficiency": "1.0"},
    ]


def efficiency_class_editor_rows(pcr: PcrForm | None) -> list[dict[str, Any]]:
    """Clases de eficiencia declaradas.

    Llegan **vacías** cuando no se han declarado, no con un reparto sugerido:
    la heterogeneidad entre moléculas es opcional y no tiene valores por
    omisión.
    """

    rows = [
        {"order": order, "efficiency": row.efficiency, "weight": row.weight}
        for order, row in enumerate(pcr.efficiency_classes if pcr else [], start=1)
    ]
    return rows or [{"order": 1, "efficiency": "", "weight": ""}]


def pcr_from_editor_rows(
    library_preparation: str,
    duplicate_fraction: str,
    cycles: int,
    efficiency: str,
    early_cycles: int | None,
    max_family_size: int | None,
    gc_max_attempts: int | None,
    gc_curve_rows: Sequence[Mapping[str, Any]],
    efficiency_class_rows: Sequence[Mapping[str, Any]],
) -> PcrForm | None:
    """Reconstruye la amplificación desde las tablas del cuaderno.

    Vive aquí, y no dentro del cuaderno, por la misma razón que la del
    tartamudeo: una función con nombre puede probarse; una expresión dentro de
    una celda de Marimo, no, y así fue como el tartamudeo se quedó fuera sin
    que nada lo notara.

    Con preparación sin amplificar devuelve `None` **aunque haya filas
    escritas**: la tabla solo se admite con preparación amplificada, y
    escribirla con la otra se rechaza al leer el documento.
    """

    if library_preparation != "pcr":
        return None
    knots = [
        GcCurveKnotRow(
            gc_fraction=str(row.get("gc_fraction", "")).strip(),
            relative_efficiency=str(row.get("relative_efficiency", "")).strip(),
        )
        for row in gc_curve_rows
        if str(row.get("gc_fraction", "")).strip()
        and str(row.get("relative_efficiency", "")).strip()
    ]
    classes = [
        EfficiencyClassRow(
            efficiency=str(row.get("efficiency", "")).strip(),
            weight=str(row.get("weight", "")).strip(),
        )
        for row in efficiency_class_rows
        if str(row.get("efficiency", "")).strip() and str(row.get("weight", "")).strip()
    ]
    # La curva neutra **no se declara**: es el valor por omisión, y escribirla
    # cambiaría los bytes del documento sin cambiar nada del resultado.
    if _is_neutral_gc_curve(knots):
        knots = []
    return PcrForm(
        duplicate_fraction=str(duplicate_fraction).strip(),
        cycles=int(cycles),
        efficiency=str(efficiency).strip(),
        early_cycles=None if early_cycles is None else int(early_cycles),
        max_family_size=None if max_family_size is None else int(max_family_size),
        gc_max_attempts=None if gc_max_attempts is None else int(gc_max_attempts),
        gc_relative_efficiency=knots,
        efficiency_classes=classes,
    )


def _is_neutral_gc_curve(knots: Sequence[GcCurveKnotRow]) -> bool:
    """¿Esta curva repite, sin tocarla, la neutra que ya es el valor omitido?"""

    if not knots:
        return True
    try:
        return all(float(knot.relative_efficiency) == 1.0 for knot in knots)
    except ValueError:
        return False


def pcr_mean_pairs_per_molecule(duplicate_fraction: str) -> float | None:
    """La `m` derivada, **con el mismo método que Rust**.

    Bisección sobre `[0, 64]` con 64 iteraciones fijas y retorno anticipado
    exacto para la fracción cero. La interfaz la muestra junto a la fracción
    declarada para que el usuario vea qué parámetro interno está pidiendo; la
    autoridad sigue siendo Rust, y una prueba compara este número con el que la
    configuración efectiva publica para la misma fracción.

    Devuelve `None` cuando la fracción no es un número o está fuera de
    dominio, en vez de inventar un valor.
    """

    try:
        target = float(duplicate_fraction)
    except (TypeError, ValueError):
        return None
    if not 0.0 <= target < 1.0:
        return None
    if target == 0.0:
        return 0.0
    if target > _pcr_expected_duplicate_fraction(PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE):
        return None
    low = 0.0
    high = PCR_MAXIMUM_MEAN_PAIRS_PER_MOLECULE
    for _ in range(64):
        midpoint = 0.5 * (low + high)
        if _pcr_expected_duplicate_fraction(midpoint) < target:
            low = midpoint
        else:
            high = midpoint
    return high


def _pcr_expected_duplicate_fraction(mean_pairs_per_molecule: float) -> float:
    """`f(m) = 1 - (1 - e^-m)/m`, con `f(0) = 0`."""

    if mean_pairs_per_molecule <= 0.0:
        return 0.0
    return 1.0 - (1.0 - math.exp(-mean_pairs_per_molecule)) / mean_pairs_per_molecule


def _is_calibrated_stutter_rate(row: StutterMotifRateRow) -> bool:
    """¿Esta fila repite, sin tocarla, la cifra que la calibración ya fija?"""

    calibrated = STUTTER_CALIBRATED_MOTIF_LENGTH_RATES.get(row.motif_length_bases)
    if calibrated is None:
        return False
    try:
        return float(row.rate) == calibrated
    except ValueError:
        return False


def stutter_from_editor_rows(
    enabled: bool,
    run_length_coefficient: str,
    motif_rate_rows: Sequence[Mapping[str, Any]],
    step_rows: Sequence[Mapping[str, Any]],
    locus_rate_rows: Sequence[Mapping[str, Any]],
) -> StutterForm | None:
    """Reconstruye el tartamudeo desde las tablas del cuaderno.

    Vive aquí, y no dentro del cuaderno, **a propósito**: la reconstrucción del
    evento técnico se escribía campo a campo dentro de una celda de Marimo, y
    así fue como el tartamudeo se quedó fuera sin que nada lo notara. Una
    función con nombre puede probarse; una expresión dentro de una celda, no.
    """

    # **Una tasa que coincide con la calibrada no cuenta como declarada.** La
    # tabla del cuaderno viene con las seis cifras ya escritas, para que quien
    # la mire sepa qué usará el simulador; si eso bastara para dar el bloque
    # por declarado, **todo** escenario importado saldría con tartamudeo, y un
    # documento sin microsatélites saldría con un bloque que el lector rechaza.
    # Declarar la cifra calibrada a mano no cambia nada: el resultado efectivo
    # es el mismo con ella y sin ella.
    motif_length_rates = [
        row_value
        for row_value in (
            StutterMotifRateRow(
                motif_length_bases=int(row["motif_length_bases"]),
                rate=str(row.get("rate", "")).strip(),
            )
            for row in motif_rate_rows
            if str(row.get("rate", "")).strip()
        )
        if not _is_calibrated_stutter_rate(row_value)
    ]
    steps = [
        StutterStepRow(units=int(row["units"]), weight=str(row.get("weight", "")).strip())
        for row in step_rows
        if str(row.get("weight", "")).strip()
    ]
    locus_rates = [
        StutterLocusRateRow(
            locus_id=str(row.get("locus_id", "")).strip(),
            rate=str(row.get("rate", "")).strip(),
        )
        for row in locus_rate_rows
        if str(row.get("locus_id", "")).strip()
    ]
    coefficient = str(run_length_coefficient).strip() or None
    if (
        not enabled
        and not motif_length_rates
        and not steps
        and not locus_rates
        and coefficient is None
    ):
        return None
    return StutterForm(
        enabled=enabled,
        motif_length_rates=motif_length_rates,
        run_length_coefficient=coefficient,
        steps=steps,
        locus_rates=locus_rates,
    )


def repeat_locus_editor_rows(loci: Sequence[RepeatLocusRow]) -> list[dict[str, Any]]:
    return [
        {
            "order": order,
            "id": locus.id,
            "sequence": locus.sequence,
            "start": locus.start,
            "end": locus.end,
            "motifs": ",".join(locus.motifs),
            "max_units": "" if locus.max_units is None else str(locus.max_units),
            "min_purity": locus.min_purity or "",
            "min_longest_run_units": (
                "" if locus.min_longest_run_units is None else str(locus.min_longest_run_units)
            ),
        }
        for order, locus in enumerate(loci, start=1)
    ]


def repeat_allele_editor_rows(loci: Sequence[RepeatLocusRow]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for locus in loci:
        for order, allele in enumerate(locus.alleles, start=1):
            rows.append(
                {
                    "order": order,
                    "locus_id": locus.id,
                    "units": (
                        "" if allele.units is None else ",".join(str(value) for value in allele.units)
                    ),
                    "sequence": allele.sequence or "",
                    "frequency": allele.frequency,
                }
            )
    return rows


def repeat_interruption_editor_rows(
    loci: Sequence[RepeatLocusRow],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for locus in loci:
        for allele_order, allele in enumerate(locus.alleles, start=1):
            for order, interruption in enumerate(allele.interruptions, start=1):
                rows.append(
                    {
                        "order": order,
                        "locus_id": locus.id,
                        "allele_order": allele_order,
                        "offset": interruption.offset,
                        "bases": interruption.bases,
                    }
                )
    return rows


def repeat_loci_from_editor_rows(
    locus_rows: Sequence[Mapping[str, Any]],
    allele_rows: Sequence[Mapping[str, Any]],
    interruption_rows: Sequence[Mapping[str, Any]],
) -> list[RepeatLocusRow]:
    """Reúne las tres tablas planas en los loci que el TOML declara.

    Una fila sin identificador de locus se omite: es el borrador vacío con el
    que la tabla arranca. Un alelo o una interrupción que nombren un locus
    inexistente **no** se omiten, se denuncian: omitirlos en silencio perdería
    justo lo que el usuario acaba de escribir.
    """
    loci: list[RepeatLocusRow] = []
    ordered_loci = ordered_editor_rows(
        [row for row in locus_rows if str(row.get("id", "")).strip()],
        "reference.repeat_loci",
    )
    known = {str(row["id"]).strip() for row in ordered_loci}
    for index, row in enumerate(ordered_loci):
        path = f"reference.repeat_loci[{index}]"
        locus_id = str(row["id"]).strip()
        motifs = [
            motif.strip()
            for motif in str(row.get("motifs", "")).split(",")
            if motif.strip()
        ]
        if not motifs:
            raise EditorError(f"{path}.motifs", "se exige al menos un motivo")
        alleles = _repeat_alleles_from_editor_rows(
            locus_id, allele_rows, interruption_rows, f"{path}.alleles", len(motifs)
        )
        loci.append(
            RepeatLocusRow(
                id=locus_id,
                sequence=_non_empty_string(
                    str(row.get("sequence", "")).strip(), f"{path}.sequence"
                ),
                start=_checked_positive_u64(
                    editor_integer(row.get("start"), f"{path}.start"), f"{path}.start"
                ),
                end=_checked_positive_u64(
                    editor_integer(row.get("end"), f"{path}.end"), f"{path}.end"
                ),
                motifs=[_repeat_motif(motif, f"{path}.motifs") for motif in motifs],
                alleles=alleles,
                max_units=_optional_editor_u64(row.get("max_units"), f"{path}.max_units"),
                min_purity=_optional_editor_literal(
                    row.get("min_purity"), f"{path}.min_purity"
                ),
                min_longest_run_units=_optional_editor_u64(
                    row.get("min_longest_run_units"), f"{path}.min_longest_run_units"
                ),
            )
        )
    for row in allele_rows:
        locus_id = str(row.get("locus_id", "")).strip()
        if locus_id and locus_id not in known:
            raise EditorError(
                "reference.repeat_loci.alleles",
                f"el alelo nombra un locus que no existe: {locus_id}",
            )
    for row in interruption_rows:
        locus_id = str(row.get("locus_id", "")).strip()
        if locus_id and locus_id not in known:
            raise EditorError(
                "reference.repeat_loci.alleles.interruptions",
                f"la interrupción nombra un locus que no existe: {locus_id}",
            )
    return loci


def _repeat_alleles_from_editor_rows(
    locus_id: str,
    allele_rows: Sequence[Mapping[str, Any]],
    interruption_rows: Sequence[Mapping[str, Any]],
    path: str,
    block_count: int,
) -> list[RepeatAlleleRow]:
    selected = ordered_editor_rows(
        [
            row
            for row in allele_rows
            if str(row.get("locus_id", "")).strip() == locus_id
        ],
        path,
    )
    if not selected:
        raise EditorError(path, f"el locus {locus_id} no declara ningún alelo")
    alleles: list[RepeatAlleleRow] = []
    for index, row in enumerate(selected):
        allele_path = f"{path}[{index}]"
        allele_order = editor_integer(row["order"], f"{allele_path}.order")
        units_text = str(row.get("units", "")).strip()
        sequence_text = str(row.get("sequence", "")).strip().upper()
        if bool(units_text) == bool(sequence_text):
            raise EditorError(
                allele_path, "se exige exactamente units o sequence"
            )
        units = None
        if units_text:
            values = [item.strip() for item in units_text.split(",") if item.strip()]
            units = [
                _checked_positive_u64(
                    editor_integer(value, f"{allele_path}.units[{value_index}]"),
                    f"{allele_path}.units[{value_index}]",
                )
                for value_index, value in enumerate(values)
            ]
            if len(units) != block_count:
                raise EditorError(
                    f"{allele_path}.units",
                    f"se exige un valor por bloque: {block_count} esperados,"
                    f" {len(units)} declarados",
                )
        interruptions = [
            RepeatInterruptionRow(
                offset=_checked_non_negative_u64(
                    editor_integer(
                        item.get("offset"), f"{allele_path}.interruptions.offset"
                    ),
                    f"{allele_path}.interruptions.offset",
                ),
                bases=_repeat_motif(
                    str(item.get("bases", "")).strip().upper(),
                    f"{allele_path}.interruptions.bases",
                ),
            )
            for item in ordered_editor_rows(
                [
                    item
                    for item in interruption_rows
                    if str(item.get("locus_id", "")).strip() == locus_id
                    and editor_integer(
                        item.get("allele_order", 0),
                        f"{allele_path}.interruptions.allele_order",
                    )
                    == allele_order
                ],
                f"{allele_path}.interruptions",
            )
        ]
        if sequence_text and interruptions:
            raise EditorError(
                f"{allele_path}.interruptions",
                "un alelo declarado por su secuencia ya lleva sus interrupciones dentro",
            )
        alleles.append(
            RepeatAlleleRow(
                frequency=_editor_non_negative_literal(
                    row.get("frequency"), f"{allele_path}.frequency"
                ),
                units=units,
                sequence=_repeat_allele_sequence(
                    sequence_text, f"{allele_path}.sequence"
                )
                if sequence_text
                else None,
                interruptions=interruptions,
            )
        )
    return alleles


def _optional_editor_u64(value: Any, path: str) -> int | None:
    text = "" if value is None else str(value).strip()
    if not text:
        return None
    return _checked_positive_u64(editor_integer(text, path), path)


def _optional_editor_literal(value: Any, path: str) -> str | None:
    text = "" if value is None else str(value).strip()
    return None if not text else _editor_non_negative_literal(text, path)


def _editor_non_negative_literal(value: Any, path: str) -> str:
    """El literal decimal tal y como se escribió en la tabla.

    Se conserva el texto y no su valor en coma flotante: el documento tiene que
    publicar lo que el usuario escribió, que es lo que la configuración efectiva
    llama frecuencia **declarada**.
    """
    text = "" if value is None else str(value).strip()
    if not text:
        raise EditorError(path, "falta el número")
    try:
        parsed = float(text)
    except ValueError as error:
        raise EditorError(
            path, "debe ser un número finito mayor o igual que cero"
        ) from error
    if not math.isfinite(parsed) or parsed < 0:
        raise EditorError(path, "debe ser un número finito mayor o igual que cero")
    return text


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
            "initial_populations",
            "direct_events",
            "technical_event",
            "vcf_output",
        },
        required={"config_kind", "schema_version", "reference"},
        path="document",
    )
    if _string(document["config_kind"], "config_kind") != "general":
        raise EditorError("config_kind", 'debe ser "general"')
    declared_version = _string(document["schema_version"], "schema_version")
    if declared_version not in {"2", "4"}:
        raise EditorError("schema_version", 'debe ser "2" o "4"')
    seed = _seed(document["seed"]) if "seed" in document else None
    ploidy = _ploidy(document["ploidy"]) if "ploidy" in document else 2
    normalization_tolerance = (
        _number_literal(document["normalization_tolerance"], "normalization_tolerance")
        if "normalization_tolerance" in document
        else "1e-9"
    )
    reference = _parse_reference(
        _mapping(document["reference"], "reference"),
        "reference",
        declares_repeats=declared_version == "4",
    )
    reductions = _parse_reductions(document.get("reduced_references", []))
    population_snvs = (
        _parse_population_snvs(
            _mapping(document["population_snvs"], "population_snvs")
        )
        if "population_snvs" in document
        else None
    )
    initial_populations = _parse_initial_populations(document.get("initial_populations", []))
    direct_events = _parse_direct_events(document.get("direct_events", []))
    technical_event = (
        _parse_technical_event(
            _mapping(document["technical_event"], "technical_event")
        )
        if "technical_event" in document
        else None
    )
    vcf_output_mode = (
        _parse_vcf_output(_mapping(document["vcf_output"], "vcf_output"))
        if "vcf_output" in document
        else None
    )
    return GeneralForm(
        seed,
        ploidy,
        normalization_tolerance,
        reference,
        reductions,
        population_snvs,
        initial_populations,
        direct_events,
        technical_event,
        vcf_output_mode,
    )


def _parse_vcf_output(table: Mapping[str, Any]) -> str:
    _keys(table, allowed={"mode"}, required={"mode"}, path="vcf_output")
    mode = _string(table["mode"], "vcf_output.mode")
    if mode not in {"text", "bgzf"}:
        raise EditorError("vcf_output.mode", "debe ser text o bgzf")
    return mode


def _parse_population_snvs(table: Mapping[str, Any]) -> PopulationSnvsForm:
    _keys(
        table,
        allowed={"mode", "input_path", "metadata_path", "generate"},
        required={"mode"},
        path="population_snvs",
    )
    mode = _string(table["mode"], "population_snvs.mode")
    input_path = table.get("input_path")
    metadata_path = table.get("metadata_path")
    if mode == "supplied":
        if metadata_path is not None or "generate" in table:
            raise EditorError("population_snvs", "supplied no admite metadata_path ni generate")
        return PopulationSnvsForm(
            mode, _non_empty_string(input_path, "population_snvs.input_path"), None
        )
    if mode == "reuse":
        if input_path is not None or "generate" in table:
            raise EditorError("population_snvs", "reuse no admite input_path ni generate")
        return PopulationSnvsForm(
            mode, None, _non_empty_string(metadata_path, "population_snvs.metadata_path")
        )
    if mode == "generate":
        if input_path is not None or metadata_path is not None:
            raise EditorError("population_snvs", "generate no admite rutas")
        generate = _mapping(table.get("generate"), "population_snvs.generate")
        return PopulationSnvsForm(mode, None, None, dict(generate))
    raise EditorError("population_snvs.mode", "debe ser supplied, reuse o generate")


def _parse_initial_populations(value: Any) -> list[InitialPopulationRow]:
    populations: list[InitialPopulationRow] = []
    for index, value_row in enumerate(_sequence(value, "initial_populations")):
        path = f"initial_populations[{index}]"
        row = _mapping(value_row, path)
        _keys(
            row,
            allowed={
                "population_id",
                "individual_count",
                "ploidy",
                "genotype_model",
                "inbreeding_coefficient",
                "individual_ids",
                "founder_keys",
            },
            required={"population_id", "individual_count", "ploidy", "genotype_model"},
            path=path,
        )
        model = _string(row["genotype_model"], f"{path}.genotype_model")
        coefficient = (
            _number_literal(row["inbreeding_coefficient"], f"{path}.inbreeding_coefficient")
            if "inbreeding_coefficient" in row
            else None
        )
        if model == "hwe" and coefficient is not None:
            raise EditorError(f"{path}.inbreeding_coefficient", "hwe no admite inbreeding_coefficient")
        if model == "f_diploid" and coefficient is None:
            raise EditorError(f"{path}.inbreeding_coefficient", "f_diploid exige inbreeding_coefficient")
        if model not in {"hwe", "f_diploid"}:
            raise EditorError(f"{path}.genotype_model", "debe ser hwe o f_diploid")
        ids = (
            [_string(item, f"{path}.individual_ids") for item in _sequence(row["individual_ids"], f"{path}.individual_ids")]
            if "individual_ids" in row
            else None
        )
        if ids is not None and not ids:
            raise EditorError(f"{path}.individual_ids", "no puede estar vacío")
        founder_keys = (
            [_string(item, f"{path}.founder_keys") for item in _sequence(row["founder_keys"], f"{path}.founder_keys")]
            if "founder_keys" in row
            else None
        )
        if founder_keys is not None and not founder_keys:
            raise EditorError(f"{path}.founder_keys", "no puede estar vacío")
        if ids is not None and founder_keys is not None:
            raise EditorError(path, "individual_ids y founder_keys son mutuamente excluyentes")
        if ids is None and founder_keys is None:
            raise EditorError(f"{path}.founder_keys", "es obligatorio si se omite individual_ids")
        count = _positive_u64(row["individual_count"], f"{path}.individual_count")
        if founder_keys is not None and len(founder_keys) != count:
            raise EditorError(f"{path}.founder_keys", "debe contener individual_count claves")
        populations.append(
            InitialPopulationRow(
                _string(row["population_id"], f"{path}.population_id"),
                count,
                _ploidy(row["ploidy"], f"{path}.ploidy"),
                model,
                coefficient,
                ids,
                founder_keys,
            )
        )
    return populations


def _parse_direct_events(value: Any) -> list[DirectEventRow]:
    events: list[DirectEventRow] = []
    common = {"event_id", "individual_id", "phase", "kind", "contig"}
    fields_by_kind = {
        "snp": {"position_1based", "ref", "alt"},
        "insertion": {"point", "bases"},
        "deletion": {"start", "end"},
    }
    for index, value_row in enumerate(_sequence(value, "direct_events")):
        path = f"direct_events[{index}]"
        row = _mapping(value_row, path)
        kind = _string(row.get("kind"), f"{path}.kind")
        branch = fields_by_kind.get(kind)
        if branch is None:
            raise EditorError(f"{path}.kind", "debe ser snp, insertion o deletion")
        _keys(row, allowed=common | branch, required=common | branch, path=path)
        common_values = {
            "event_id": _string(row["event_id"], f"{path}.event_id"),
            "individual_id": _string(row["individual_id"], f"{path}.individual_id"),
            "phase": _non_negative_u64(row["phase"], f"{path}.phase"),
            "kind": kind,
            "contig": _string(row["contig"], f"{path}.contig"),
        }
        if kind == "snp":
            events.append(DirectEventRow(
                **common_values,
                position_1based=_positive_u64(row["position_1based"], f"{path}.position_1based"),
                reference=_string(row["ref"], f"{path}.ref"),
                alternate=_string(row["alt"], f"{path}.alt"),
            ))
        elif kind == "insertion":
            events.append(DirectEventRow(
                **common_values,
                point=_non_negative_u64(row["point"], f"{path}.point"),
                bases=_string(row["bases"], f"{path}.bases"),
            ))
        else:
            events.append(DirectEventRow(
                **common_values,
                start=_non_negative_u64(row["start"], f"{path}.start"),
                end=_non_negative_u64(row["end"], f"{path}.end"),
            ))
    return events


def _parse_technical_event(table: Mapping[str, Any]) -> TechnicalEventForm:
    _keys(
        table,
        allowed={
            "event_id",
            "individual_ids",
            "target_coverage_wgs",
            "read_length_bases",
            "stutter",
            "library_preparation",
            "pcr",
            "libraries",
        },
        required={"event_id", "individual_ids", "target_coverage_wgs", "libraries"},
        path="technical_event",
    )
    individual_ids = [
        _string(value, "technical_event.individual_ids")
        for value in _sequence(table["individual_ids"], "technical_event.individual_ids")
    ]
    if not individual_ids:
        raise EditorError("technical_event.individual_ids", "no puede estar vacío")
    if len(individual_ids) != len(set(individual_ids)):
        raise EditorError("technical_event.individual_ids", "contiene IDs repetidos")
    libraries: list[TechnicalLibraryRow] = []
    for index, value in enumerate(_sequence(table["libraries"], "technical_event.libraries")):
        path = f"technical_event.libraries[{index}]"
        row = _mapping(value, path)
        _keys(row, allowed={"individual_id", "i7", "i5"}, required={"individual_id", "i7", "i5"}, path=path)
        i7 = _udi(row["i7"], f"{path}.i7")
        i5 = _udi(row["i5"], f"{path}.i5")
        libraries.append(TechnicalLibraryRow(_string(row["individual_id"], f"{path}.individual_id"), i7, i5))
    if not libraries:
        raise EditorError("technical_event.libraries", "no puede estar vacío")
    return TechnicalEventForm(
        event_id=_string(table["event_id"], "technical_event.event_id"),
        individual_ids=individual_ids,
        coverage_wgs=_non_negative_number_literal(
            table["target_coverage_wgs"], "technical_event.target_coverage_wgs"
        ),
        libraries=libraries,
        read_length_bases=(
            _read_length(table["read_length_bases"], "technical_event.read_length_bases")
            if "read_length_bases" in table
            else None
        ),
        stutter=(
            _parse_stutter(_mapping(table["stutter"], "technical_event.stutter"))
            if "stutter" in table
            else None
        ),
        library_preparation=(
            _library_preparation(table["library_preparation"])
            if "library_preparation" in table
            else None
        ),
        pcr=(
            _parse_pcr(_mapping(table["pcr"], "technical_event.pcr"))
            if "pcr" in table
            else None
        ),
    )


def _library_preparation(value: Any) -> str:
    preparation = _string(value, "technical_event.library_preparation")
    if preparation not in LIBRARY_PREPARATIONS:
        raise EditorError(
            "technical_event.library_preparation",
            "las preparaciones admitidas son " + " y ".join(LIBRARY_PREPARATIONS),
        )
    return preparation


def _parse_pcr(table: Mapping[str, Any]) -> PcrForm:
    path = "technical_event.pcr"
    _keys(
        table,
        allowed={
            "duplicate_fraction",
            "cycles",
            "efficiency",
            "early_cycles",
            "max_family_size",
            "gc_max_attempts",
            "gc_relative_efficiency",
            "efficiency_classes",
        },
        required={"duplicate_fraction", "cycles", "efficiency"},
        path=path,
    )
    knots: list[GcCurveKnotRow] = []
    for index, value in enumerate(
        _sequence(table.get("gc_relative_efficiency", []), f"{path}.gc_relative_efficiency")
    ):
        row_path = f"{path}.gc_relative_efficiency[{index}]"
        row = _mapping(value, row_path)
        _keys(
            row,
            allowed={"gc_fraction", "relative_efficiency"},
            required={"gc_fraction", "relative_efficiency"},
            path=row_path,
        )
        knots.append(
            GcCurveKnotRow(
                _fraction_literal(row["gc_fraction"], f"{row_path}.gc_fraction"),
                _fraction_literal(
                    row["relative_efficiency"], f"{row_path}.relative_efficiency"
                ),
            )
        )
    classes: list[EfficiencyClassRow] = []
    for index, value in enumerate(
        _sequence(table.get("efficiency_classes", []), f"{path}.efficiency_classes")
    ):
        row_path = f"{path}.efficiency_classes[{index}]"
        row = _mapping(value, row_path)
        _keys(
            row,
            allowed={"efficiency", "weight"},
            required={"efficiency", "weight"},
            path=row_path,
        )
        classes.append(
            EfficiencyClassRow(
                _fraction_literal(row["efficiency"], f"{row_path}.efficiency"),
                _fraction_literal(row["weight"], f"{row_path}.weight"),
            )
        )
    if len(classes) > PCR_MAXIMUM_EFFICIENCY_CLASSES:
        raise EditorError(
            f"{path}.efficiency_classes",
            f"el máximo de clases de eficiencia es {PCR_MAXIMUM_EFFICIENCY_CLASSES}",
        )
    return PcrForm(
        duplicate_fraction=_fraction_literal(
            table["duplicate_fraction"], f"{path}.duplicate_fraction"
        ),
        cycles=_bounded_count(table["cycles"], f"{path}.cycles", 1, PCR_MAXIMUM_CYCLES),
        efficiency=_fraction_literal(table["efficiency"], f"{path}.efficiency"),
        early_cycles=(
            _bounded_count(
                table["early_cycles"], f"{path}.early_cycles", 0, PCR_MAXIMUM_EARLY_CYCLES
            )
            if "early_cycles" in table
            else None
        ),
        max_family_size=(
            _bounded_count(
                table["max_family_size"],
                f"{path}.max_family_size",
                1,
                PCR_MAXIMUM_FAMILY_SIZE,
            )
            if "max_family_size" in table
            else None
        ),
        gc_max_attempts=(
            _bounded_count(
                table["gc_max_attempts"],
                f"{path}.gc_max_attempts",
                1,
                PCR_MAXIMUM_GC_ATTEMPTS,
            )
            if "gc_max_attempts" in table
            else None
        ),
        gc_relative_efficiency=knots,
        efficiency_classes=classes,
    )


def _bounded_count(value: Any, path: str, minimum: int, maximum: int) -> int:
    """Un entero dentro de su dominio de contrato, nombrando los extremos."""

    if not isinstance(value, Integer):
        raise EditorError(path, "se esperaba un entero")
    number = int(value)
    if not minimum <= number <= maximum:
        raise EditorError(path, f"tiene que estar entre {minimum} y {maximum}")
    return number


def _parse_stutter(table: Mapping[str, Any]) -> StutterForm:
    path = "technical_event.stutter"
    _keys(
        table,
        allowed={
            "enabled",
            "motif_length_rates",
            "run_length_coefficient",
            "steps",
            "locus_rates",
        },
        required=set(),
        path=path,
    )
    enabled = bool(table.get("enabled", False))
    rates: list[StutterMotifRateRow] = []
    for index, value in enumerate(
        _sequence(table.get("motif_length_rates", []), f"{path}.motif_length_rates")
    ):
        row_path = f"{path}.motif_length_rates[{index}]"
        row = _mapping(value, row_path)
        _keys(
            row,
            allowed={"motif_length_bases", "rate"},
            required={"motif_length_bases", "rate"},
            path=row_path,
        )
        length = _positive_u64(row["motif_length_bases"], f"{row_path}.motif_length_bases")
        if length not in STUTTER_MOTIF_LENGTHS:
            raise EditorError(
                f"{row_path}.motif_length_bases",
                "la calibración solo cubre longitudes de motivo de 1 a 6 pb",
            )
        rates.append(
            StutterMotifRateRow(length, _fraction_literal(row["rate"], f"{row_path}.rate"))
        )
    steps: list[StutterStepRow] = []
    for index, value in enumerate(_sequence(table.get("steps", []), f"{path}.steps")):
        row_path = f"{path}.steps[{index}]"
        row = _mapping(value, row_path)
        _keys(row, allowed={"units", "weight"}, required={"units", "weight"}, path=row_path)
        units = row["units"]
        if not isinstance(units, Integer):
            raise EditorError(f"{row_path}.units", "se esperaba un entero de unidades")
        units = int(units)
        if units not in STUTTER_ADMITTED_STEPS:
            raise EditorError(
                f"{row_path}.units",
                "la estructura del paso la fija el algoritmo: solo -3, -2, -1 y +1",
            )
        steps.append(
            StutterStepRow(
                units, _non_negative_number_literal(row["weight"], f"{row_path}.weight")
            )
        )
    if enabled and {row.units for row in steps} != set(STUTTER_ADMITTED_STEPS):
        raise EditorError(
            f"{path}.steps",
            "con el tartamudeo activo hay que declarar el peso de los cuatro pasos admitidos",
        )
    locus_rates: list[StutterLocusRateRow] = []
    for index, value in enumerate(
        _sequence(table.get("locus_rates", []), f"{path}.locus_rates")
    ):
        row_path = f"{path}.locus_rates[{index}]"
        row = _mapping(value, row_path)
        _keys(
            row, allowed={"locus_id", "rate"}, required={"locus_id", "rate"}, path=row_path
        )
        locus_rates.append(
            StutterLocusRateRow(
                _non_empty_string(row["locus_id"], f"{row_path}.locus_id"),
                _fraction_literal(row["rate"], f"{row_path}.rate"),
            )
        )
    return StutterForm(
        enabled=enabled,
        motif_length_rates=rates,
        run_length_coefficient=(
            _non_negative_number_literal(
                table["run_length_coefficient"], f"{path}.run_length_coefficient"
            )
            if "run_length_coefficient" in table
            else None
        ),
        steps=steps,
        locus_rates=locus_rates,
    )


def _parse_reference(
    table: Mapping[str, Any], path: str, declares_repeats: bool = False
) -> ReferenceForm:
    """Lee la referencia declarada.

    Las claves de repetición pertenecen al documento **v4 y solo a él**: en un
    v2 son campos desconocidos, exactamente igual que en el lector de Rust. Sin
    esa condición la interfaz admitiría un documento que el simulador rechaza.
    """
    kind = _string(table.get("kind"), f"{path}.kind")
    common = {"kind", "id", "sequences"}
    repeat_keys: set[str] = (
        {"repeat_loci", "repeat_admission", "max_tract_output_bases"}
        if declares_repeats
        else set()
    )
    if kind == "simulated":
        _keys(
            table,
            allowed=common | {"composition"} | repeat_keys,
            required=common | {"composition"},
            path=path,
        )
        composition_table = _mapping(table["composition"], f"{path}.composition")
        _keys(
            composition_table,
            allowed=set(COMPOSITION_BASES),
            required=set(BASES),
            path=f"{path}.composition",
        )
        composition = {
            base: _number_literal(
                composition_table[base], f"{path}.composition.{base}"
            )
            for base in BASES
        }
        if "N" in composition_table:
            composition["N"] = _number_literal(
                composition_table["N"], f"{path}.composition.N"
            )
        fasta = None
    elif kind == "supplied":
        _keys(
            table,
            allowed=common | {"fasta"} | repeat_keys,
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
    repeat_loci = (
        _parse_repeat_loci(table.get("repeat_loci", []), f"{path}.repeat_loci")
        if declares_repeats
        else []
    )
    repeat_admission = (
        _parse_repeat_admission(
            _mapping(table["repeat_admission"], f"{path}.repeat_admission"),
            f"{path}.repeat_admission",
        )
        if declares_repeats and "repeat_admission" in table
        else None
    )
    return ReferenceForm(
        kind=kind,
        id=_string(table["id"], f"{path}.id"),
        fasta=fasta,
        composition=composition,
        sequences=sequences,
        repeat_loci=repeat_loci,
        repeat_admission=repeat_admission,
        max_tract_output_bases=(
            _positive_u64(
                table["max_tract_output_bases"], f"{path}.max_tract_output_bases"
            )
            if declares_repeats and "max_tract_output_bases" in table
            else None
        ),
    )


def _parse_repeat_admission(
    table: Mapping[str, Any], path: str
) -> RepeatAdmissionForm:
    _keys(
        table,
        allowed={"min_purity", "min_longest_run_units"},
        required=set(),
        path=path,
    )
    return RepeatAdmissionForm(
        min_purity=(
            _number_literal(table["min_purity"], f"{path}.min_purity")
            if "min_purity" in table
            else None
        ),
        min_longest_run_units=(
            _positive_u64(
                table["min_longest_run_units"], f"{path}.min_longest_run_units"
            )
            if "min_longest_run_units" in table
            else None
        ),
    )


def _parse_repeat_loci(value: Any, path: str) -> list[RepeatLocusRow]:
    rows = _sequence(value, path)
    loci: list[RepeatLocusRow] = []
    for index, entry in enumerate(rows):
        locus_path = f"{path}[{index}]"
        locus = _mapping(entry, locus_path)
        _keys(
            locus,
            allowed={
                "id",
                "sequence",
                "start",
                "end",
                "motif",
                "motifs",
                "max_units",
                "min_purity",
                "min_longest_run_units",
                "alleles",
            },
            required={"id", "sequence", "start", "end", "alleles"},
            path=locus_path,
        )
        motifs = _parse_repeat_motifs(locus, locus_path)
        alleles = _parse_repeat_alleles(
            locus["alleles"], f"{locus_path}.alleles", len(motifs)
        )
        loci.append(
            RepeatLocusRow(
                id=_non_empty_string(locus["id"], f"{locus_path}.id"),
                sequence=_non_empty_string(locus["sequence"], f"{locus_path}.sequence"),
                start=_positive_u64(locus["start"], f"{locus_path}.start"),
                end=_positive_u64(locus["end"], f"{locus_path}.end"),
                motifs=motifs,
                alleles=alleles,
                max_units=(
                    _positive_u64(locus["max_units"], f"{locus_path}.max_units")
                    if "max_units" in locus
                    else None
                ),
                min_purity=(
                    _number_literal(locus["min_purity"], f"{locus_path}.min_purity")
                    if "min_purity" in locus
                    else None
                ),
                min_longest_run_units=(
                    _positive_u64(
                        locus["min_longest_run_units"],
                        f"{locus_path}.min_longest_run_units",
                    )
                    if "min_longest_run_units" in locus
                    else None
                ),
            )
        )
    return loci


def _parse_repeat_motifs(locus: Mapping[str, Any], path: str) -> list[str]:
    """`motif` y `motifs` son excluyentes, igual que en el lector de Rust."""
    singular = "motif" in locus
    plural = "motifs" in locus
    if singular and plural:
        raise EditorError(f"{path}.motif", "motif y motifs son excluyentes")
    if singular:
        return [_repeat_motif(locus["motif"], f"{path}.motif")]
    if not plural:
        raise EditorError(f"{path}.motif", "se exige motif o motifs")
    values = _sequence(locus["motifs"], f"{path}.motifs")
    if not values:
        raise EditorError(f"{path}.motifs", "la lista de motivos no puede estar vacía")
    return [
        _repeat_motif(value, f"{path}.motifs[{index}]")
        for index, value in enumerate(values)
    ]


def _repeat_motif(value: Any, path: str) -> str:
    """Un motivo son bases A, C, G o T. `N` queda fuera a propósito: un motivo
    ambiguo no define ninguna copia exacta y ni la pureza ni la tirada más larga
    serían computables."""
    motif = _non_empty_string(value, path)
    if any(base not in BASES for base in motif):
        raise EditorError(path, "el motivo solo admite A, C, G o T")
    return motif


def _parse_repeat_alleles(
    value: Any, path: str, block_count: int
) -> list[RepeatAlleleRow]:
    rows = _sequence(value, path)
    if not rows:
        raise EditorError(path, "el locus no declara ningún alelo")
    alleles: list[RepeatAlleleRow] = []
    for index, entry in enumerate(rows):
        allele_path = f"{path}[{index}]"
        allele = _mapping(entry, allele_path)
        _keys(
            allele,
            allowed={"units", "sequence", "frequency", "interruptions"},
            required={"frequency"},
            path=allele_path,
        )
        units = (
            _parse_repeat_units(allele["units"], f"{allele_path}.units", block_count)
            if "units" in allele
            else None
        )
        sequence = (
            _repeat_allele_sequence(allele["sequence"], f"{allele_path}.sequence")
            if "sequence" in allele
            else None
        )
        if (units is None) == (sequence is None):
            raise EditorError(
                allele_path, "se exige exactamente units o sequence"
            )
        if sequence is not None and allele.get("interruptions"):
            # Un alelo declarado por su secuencia exacta ya lleva dentro sus
            # impurezas: declararlas además sería decir dos veces lo mismo y
            # poder contradecirse. El lector de Rust lo rechaza igual.
            raise EditorError(
                f"{allele_path}.interruptions",
                "un alelo declarado por su secuencia ya lleva sus interrupciones dentro",
            )
        interruptions = [
            _parse_repeat_interruption(
                item, f"{allele_path}.interruptions[{item_index}]"
            )
            for item_index, item in enumerate(
                _sequence(allele.get("interruptions", []), f"{allele_path}.interruptions")
            )
        ]
        alleles.append(
            RepeatAlleleRow(
                frequency=_non_negative_number_literal(
                    allele["frequency"], f"{allele_path}.frequency"
                ),
                units=units,
                sequence=sequence,
                interruptions=interruptions,
            )
        )
    return alleles


def _parse_repeat_units(value: Any, path: str, block_count: int) -> list[int]:
    """`units` admite un entero cuando el locus tiene un solo motivo y una lista
    con un valor por bloque cuando tiene varios."""
    if isinstance(value, (list, tuple)) or (
        hasattr(value, "__iter__") and not isinstance(value, (str, bytes, int))
    ):
        units = [
            _positive_u64(item, f"{path}[{index}]")
            for index, item in enumerate(_sequence(value, path))
        ]
    else:
        units = [_positive_u64(value, path)]
    if len(units) != block_count:
        raise EditorError(
            path,
            f"se exige un valor por bloque: {block_count} esperados, {len(units)} declarados",
        )
    return units


def _repeat_allele_sequence(value: Any, path: str) -> str:
    sequence = _non_empty_string(value, path)
    if any(base not in BASES for base in sequence):
        raise EditorError(path, "la secuencia del alelo solo admite A, C, G o T")
    return sequence


def _parse_repeat_interruption(value: Any, path: str) -> RepeatInterruptionRow:
    table = _mapping(value, path)
    _keys(table, allowed={"offset", "bases"}, required={"offset", "bases"}, path=path)
    bases = _non_empty_string(table["bases"], f"{path}.bases")
    if any(base not in BASES for base in bases):
        raise EditorError(f"{path}.bases", "la interrupción solo admite A, C, G o T")
    return RepeatInterruptionRow(
        offset=_non_negative_u64(table["offset"], f"{path}.offset"),
        bases=bases,
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


def normalize_repeat_frequencies(
    literals: Sequence[str], tolerance_literal: str, path: str
) -> list[float]:
    """Normaliza el reparto de frecuencias de un locus y **rechaza el inválido**.

    Es la misma regla que aplica el simulador (DTC-049): se suman los valores
    con compensación de Neumaier en un orden estable, se exige que la distancia
    a uno no pase de la tolerancia, y solo entonces se divide cada valor por la
    suma. Aquí se hace **en la interfaz** para que un reparto imposible se vea
    antes de lanzar nada, no para sustituir al lector: quien manda sigue siendo
    Rust, y por eso el TOML conserva los literales que el usuario escribió.

    **Se comprueba al escribir y no al leer, a propósito.** Importar un
    documento con un reparto inválido tiene que funcionar: si la lectura lo
    rechazara, el usuario no podría abrir en la interfaz precisamente el
    fichero que necesita corregir, y tendría que arreglarlo en un editor de
    texto. Como toda ida y vuelta pasa por la escritura, ningún documento
    inválido llega a publicarse; lo que se gana es poder corregirlo donde se
    ve.
    """
    try:
        tolerance = float(tolerance_literal)
    except (TypeError, ValueError) as error:
        raise EditorError(path, "la tolerancia no es un decimal válido") from error
    if not math.isfinite(tolerance) or not 0.0 <= tolerance <= 1e-6:
        raise EditorError(path, "la tolerancia debe pertenecer a [0, 1e-6]")
    if not literals:
        raise EditorError(path, "el locus no declara ningún alelo")
    values: list[float] = []
    for index, literal in enumerate(literals):
        try:
            value = float(literal)
        except (TypeError, ValueError) as error:
            raise EditorError(
                f"{path}[{index}].frequency", "no es un decimal válido"
            ) from error
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise EditorError(
                f"{path}[{index}].frequency", "debe ser finita y pertenecer a [0,1]"
            )
        values.append(value)
    total = _neumaier_sum(values)
    if not math.isfinite(total) or total <= 0.0:
        raise EditorError(path, "la suma de las frecuencias no es positiva y finita")
    if abs(total - 1.0) > tolerance:
        raise EditorError(
            path,
            f"las frecuencias suman {total!r} y su distancia a uno supera la tolerancia {tolerance!r}",
        )
    return [value / total for value in values]


def _neumaier_sum(values: Sequence[float]) -> float:
    """Suma compensada de Neumaier, la misma que el simulador usa."""
    total = 0.0
    compensation = 0.0
    for value in values:
        candidate = total + value
        if abs(total) >= abs(value):
            compensation += (total - candidate) + value
        else:
            compensation += (value - candidate) + total
        total = candidate
    return total + compensation


def render_general_toml(form: GeneralForm) -> str:
    document = tomlkit.document()
    document.add("config_kind", "general")
    # La versión del documento la decide **lo declarado**: un escenario sin
    # loci de repetición sigue escribiendo el documento de siempre, byte a
    # byte, y solo uno que los declare sube a la versión que los conoce.
    declares_repeats = (
        bool(form.reference.repeat_loci)
        or form.reference.repeat_admission is not None
        or form.reference.max_tract_output_bases is not None
        # **El bloque de tartamudeo cuenta.** La tabla del evento técnico solo
        # existe en el documento que puede declarar microsatélites, así que un
        # documento que la lleve y salga en la versión anterior es un fichero
        # que el simulador rechaza como campo desconocido. Mirarlo aquí es lo
        # que impide escribir ese fichero.
        or (
            form.technical_event is not None
            and form.technical_event.stutter is not None
        )
        # **Y la preparación amplificada también.** Sus dos claves viven en la
        # misma versión que el tartamudeo, por la misma razón: el lector de
        # Rust las admite solo en el documento que puede declarar
        # microsatélites, y escribirlas en el anterior produce un fichero que
        # el simulador rechaza como campo desconocido.
        or (
            form.technical_event is not None
            and (
                form.technical_event.library_preparation is not None
                or form.technical_event.pcr is not None
            )
        )
    )
    document.add("schema_version", "4" if declares_repeats else "2")
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
        validate_simulated_composition(form.reference.composition)
        composition = tomlkit.table()
        for base in BASES:
            composition.add(
                base,
                _literal_item(
                    form.reference.composition.get(base),
                    f"reference.composition.{base}",
                ),
            )
        if "N" in form.reference.composition:
            composition.add(
                "N",
                _literal_item(
                    form.reference.composition["N"], "reference.composition.N"
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
    if form.reference.repeat_admission is not None:
        admission = form.reference.repeat_admission
        table = tomlkit.table()
        if admission.min_purity is not None:
            table.add(
                "min_purity",
                _literal_item(admission.min_purity, "reference.repeat_admission.min_purity"),
            )
        if admission.min_longest_run_units is not None:
            table.add(
                "min_longest_run_units",
                _checked_positive_u64(
                    admission.min_longest_run_units,
                    "reference.repeat_admission.min_longest_run_units",
                ),
            )
        if not table:
            raise EditorError(
                "reference.repeat_admission",
                "no declara ningún mínimo; omítala en vez de dejarla vacía",
            )
        reference.add("repeat_admission", table)
    if form.reference.max_tract_output_bases is not None:
        reference.add(
            "max_tract_output_bases",
            _checked_positive_u64(
                form.reference.max_tract_output_bases,
                "reference.max_tract_output_bases",
            ),
        )
    if form.reference.repeat_loci:
        reference.add("repeat_loci", _render_repeat_loci(form))
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
            generate = population_snvs.generate
            if generate is None:
                raise EditorError("population_snvs.generate", "falta la configuración de generación")
            generated = tomlkit.table()
            for key, value in generate.items():
                generated.add(key, value)
            table.add("generate", generated)
        else:
            raise EditorError(
                "population_snvs.mode", "debe ser supplied, reuse o generate"
            )
        document.add("population_snvs", table)
    if form.initial_populations:
        populations = tomlkit.aot()
        for index, population in enumerate(form.initial_populations):
            path = f"initial_populations[{index}]"
            row = tomlkit.table()
            row.add("population_id", population.population_id)
            row.add(
                "individual_count",
                _checked_positive_u64(population.individual_count, f"{path}.individual_count"),
            )
            row.add("ploidy", _checked_ploidy(population.ploidy, f"{path}.ploidy"))
            if population.genotype_model not in {"hwe", "f_diploid"}:
                raise EditorError(f"{path}.genotype_model", "debe ser hwe o f_diploid")
            row.add("genotype_model", population.genotype_model)
            if population.genotype_model == "f_diploid":
                row.add(
                    "inbreeding_coefficient",
                    _literal_item(population.inbreeding_coefficient, f"{path}.inbreeding_coefficient"),
                )
            elif population.inbreeding_coefficient is not None:
                raise EditorError(f"{path}.inbreeding_coefficient", "hwe no admite inbreeding_coefficient")
            if population.individual_ids is not None:
                if not population.individual_ids:
                    raise EditorError(f"{path}.individual_ids", "no puede estar vacío")
                row.add("individual_ids", population.individual_ids)
            elif population.founder_keys is not None:
                if not population.founder_keys:
                    raise EditorError(f"{path}.founder_keys", "no puede estar vacío")
                if len(population.founder_keys) != population.individual_count:
                    raise EditorError(f"{path}.founder_keys", "debe contener individual_count claves")
                row.add("founder_keys", population.founder_keys)
            else:
                raise EditorError(f"{path}.founder_keys", "es obligatorio si se omite individual_ids")
            populations.append(row)
        document.add("initial_populations", populations)
    if form.direct_events:
        events = tomlkit.aot()
        for index, event in enumerate(form.direct_events):
            path = f"direct_events[{index}]"
            row = tomlkit.table()
            row.add("event_id", event.event_id)
            row.add("individual_id", event.individual_id)
            row.add("phase", _checked_non_negative_u64(event.phase, f"{path}.phase"))
            row.add("kind", event.kind)
            row.add("contig", event.contig)
            if event.kind == "snp":
                row.add("position_1based", _checked_positive_u64(event.position_1based, f"{path}.position_1based"))
                row.add("ref", _required_text(event.reference, f"{path}.ref"))
                row.add("alt", _required_text(event.alternate, f"{path}.alt"))
            elif event.kind == "insertion":
                row.add("point", _checked_non_negative_u64(event.point, f"{path}.point"))
                row.add("bases", _required_text(event.bases, f"{path}.bases"))
            elif event.kind == "deletion":
                row.add("start", _checked_non_negative_u64(event.start, f"{path}.start"))
                row.add("end", _checked_non_negative_u64(event.end, f"{path}.end"))
            else:
                raise EditorError(f"{path}.kind", "debe ser snp, insertion o deletion")
            events.append(row)
        document.add("direct_events", events)
    if form.technical_event is not None:
        event = form.technical_event
        if not event.individual_ids:
            raise EditorError("technical_event.individual_ids", "no puede estar vacío")
        if len(event.individual_ids) != len(set(event.individual_ids)):
            raise EditorError("technical_event.individual_ids", "contiene IDs repetidos")
        if not event.libraries:
            raise EditorError("technical_event.libraries", "no puede estar vacío")
        technical_event = tomlkit.table()
        technical_event.add("event_id", event.event_id)
        technical_event.add("individual_ids", event.individual_ids)
        technical_event.add(
            "target_coverage_wgs",
            _non_negative_literal_item(
                event.coverage_wgs, "technical_event.target_coverage_wgs"
            ),
        )
        if event.read_length_bases is not None:
            technical_event.add(
                "read_length_bases",
                _checked_read_length(
                    event.read_length_bases, "technical_event.read_length_bases"
                ),
            )
        if event.stutter is not None:
            technical_event.add("stutter", _render_stutter(event.stutter))
        if event.library_preparation is not None:
            if event.library_preparation not in LIBRARY_PREPARATIONS:
                raise EditorError(
                    "technical_event.library_preparation",
                    "las preparaciones admitidas son " + " y ".join(LIBRARY_PREPARATIONS),
                )
            technical_event.add("library_preparation", event.library_preparation)
        if event.pcr is not None:
            if event.library_preparation != "pcr":
                raise EditorError(
                    "technical_event.pcr",
                    "la tabla de amplificación solo se admite con library_preparation = \"pcr\"",
                )
            technical_event.add("pcr", _render_pcr(event.pcr))
        elif event.library_preparation == "pcr":
            raise EditorError(
                "technical_event.pcr",
                "la preparación pcr exige la tabla technical_event.pcr, cuyos tres primeros "
                "parámetros son obligatorios y no tienen valor por omisión",
            )
        libraries = tomlkit.aot()
        for index, library in enumerate(event.libraries):
            path = f"technical_event.libraries[{index}]"
            row = tomlkit.table()
            row.add("individual_id", library.individual_id)
            row.add("i7", _udi(library.i7, f"{path}.i7"))
            row.add("i5", _udi(library.i5, f"{path}.i5"))
            libraries.append(row)
        technical_event.add("libraries", libraries)
        document.add("technical_event", technical_event)
    if form.vcf_output_mode is not None:
        if form.vcf_output_mode not in {"text", "bgzf"}:
            raise EditorError("vcf_output.mode", "debe ser text o bgzf")
        vcf_output = tomlkit.table()
        vcf_output.add("mode", form.vcf_output_mode)
        document.add("vcf_output", vcf_output)
    return tomlkit.dumps(document)


def _render_stutter(stutter: StutterForm) -> Any:
    """Vuelca el tartamudeo al TOML conservando los literales declarados.

    Los pesos se escriben **tal cual se declararon**, sin reordenar ni
    completar: si faltara alguno, la ida y la vuelta lo delataría en vez de
    taparlo con un reparto implícito.
    """

    path = "technical_event.stutter"
    if stutter.enabled and {row.units for row in stutter.steps} != set(
        STUTTER_ADMITTED_STEPS
    ):
        raise EditorError(
            f"{path}.steps",
            "con el tartamudeo activo hay que declarar el peso de los cuatro pasos admitidos",
        )
    table = tomlkit.table()
    table.add("enabled", stutter.enabled)
    if stutter.motif_length_rates:
        rates = tomlkit.aot()
        for index, row in enumerate(stutter.motif_length_rates):
            row_path = f"{path}.motif_length_rates[{index}]"
            if row.motif_length_bases not in STUTTER_MOTIF_LENGTHS:
                raise EditorError(
                    f"{row_path}.motif_length_bases",
                    "la calibración solo cubre longitudes de motivo de 1 a 6 pb",
                )
            entry = tomlkit.table()
            entry.add("motif_length_bases", row.motif_length_bases)
            entry.add("rate", _fraction_literal_item(row.rate, f"{row_path}.rate"))
            rates.append(entry)
        table.add("motif_length_rates", rates)
    if stutter.run_length_coefficient is not None:
        table.add(
            "run_length_coefficient",
            _non_negative_literal_item(
                stutter.run_length_coefficient, f"{path}.run_length_coefficient"
            ),
        )
    if stutter.steps:
        steps = tomlkit.aot()
        for index, row in enumerate(stutter.steps):
            row_path = f"{path}.steps[{index}]"
            if row.units not in STUTTER_ADMITTED_STEPS:
                raise EditorError(
                    f"{row_path}.units",
                    "la estructura del paso la fija el algoritmo: solo -3, -2, -1 y +1",
                )
            entry = tomlkit.table()
            entry.add("units", row.units)
            entry.add(
                "weight", _non_negative_literal_item(row.weight, f"{row_path}.weight")
            )
            steps.append(entry)
        table.add("steps", steps)
    if stutter.locus_rates:
        locus_rates = tomlkit.aot()
        for index, row in enumerate(stutter.locus_rates):
            row_path = f"{path}.locus_rates[{index}]"
            entry = tomlkit.table()
            entry.add("locus_id", _required_text(row.locus_id, f"{row_path}.locus_id"))
            entry.add("rate", _fraction_literal_item(row.rate, f"{row_path}.rate"))
            locus_rates.append(entry)
        table.add("locus_rates", locus_rates)
    return table


def _render_pcr(pcr: PcrForm) -> Any:
    """Vuelca la amplificación al TOML conservando los literales declarados.

    Lo que no se declaró **no se escribe**: los cuatro parámetros opcionales y
    las dos tablas solo aparecen si el usuario los puso. Escribirlos con su
    valor por omisión cambiaría los bytes del documento sin cambiar nada del
    resultado, y la ida y la vuelta dejaría de ser una identidad.
    """

    path = "technical_event.pcr"
    table = tomlkit.table()
    table.add(
        "duplicate_fraction",
        _fraction_literal_item(pcr.duplicate_fraction, f"{path}.duplicate_fraction"),
    )
    table.add("cycles", _checked_count(pcr.cycles, f"{path}.cycles", 1, PCR_MAXIMUM_CYCLES))
    table.add("efficiency", _fraction_literal_item(pcr.efficiency, f"{path}.efficiency"))
    if pcr.early_cycles is not None:
        table.add(
            "early_cycles",
            _checked_count(
                pcr.early_cycles, f"{path}.early_cycles", 0, PCR_MAXIMUM_EARLY_CYCLES
            ),
        )
    if pcr.max_family_size is not None:
        table.add(
            "max_family_size",
            _checked_count(
                pcr.max_family_size,
                f"{path}.max_family_size",
                1,
                PCR_MAXIMUM_FAMILY_SIZE,
            ),
        )
    if pcr.gc_max_attempts is not None:
        table.add(
            "gc_max_attempts",
            _checked_count(
                pcr.gc_max_attempts,
                f"{path}.gc_max_attempts",
                1,
                PCR_MAXIMUM_GC_ATTEMPTS,
            ),
        )
    if pcr.gc_relative_efficiency:
        if len(pcr.gc_relative_efficiency) < 2:
            raise EditorError(
                f"{path}.gc_relative_efficiency",
                "una curva de composición necesita al menos dos nudos",
            )
        knots = tomlkit.aot()
        for index, row in enumerate(pcr.gc_relative_efficiency):
            row_path = f"{path}.gc_relative_efficiency[{index}]"
            entry = tomlkit.table()
            entry.add(
                "gc_fraction",
                _fraction_literal_item(row.gc_fraction, f"{row_path}.gc_fraction"),
            )
            entry.add(
                "relative_efficiency",
                _fraction_literal_item(
                    row.relative_efficiency, f"{row_path}.relative_efficiency"
                ),
            )
            knots.append(entry)
        table.add("gc_relative_efficiency", knots)
    if pcr.efficiency_classes:
        if len(pcr.efficiency_classes) > PCR_MAXIMUM_EFFICIENCY_CLASSES:
            raise EditorError(
                f"{path}.efficiency_classes",
                f"el máximo de clases de eficiencia es {PCR_MAXIMUM_EFFICIENCY_CLASSES}",
            )
        classes = tomlkit.aot()
        for index, row in enumerate(pcr.efficiency_classes):
            row_path = f"{path}.efficiency_classes[{index}]"
            entry = tomlkit.table()
            entry.add(
                "efficiency", _fraction_literal_item(row.efficiency, f"{row_path}.efficiency")
            )
            entry.add("weight", _fraction_literal_item(row.weight, f"{row_path}.weight"))
            classes.append(entry)
        table.add("efficiency_classes", classes)
    return table


def _checked_count(value: int, path: str, minimum: int, maximum: int) -> int:
    if not minimum <= value <= maximum:
        raise EditorError(path, f"tiene que estar entre {minimum} y {maximum}")
    return value


def _render_repeat_loci(form: GeneralForm) -> Any:
    """Escribe los loci declarados.

    El motivo se escribe en singular cuando hay uno solo, que es el atajo del
    contrato y lo que escribe todo documento que no sea compuesto. La ida y la
    vuelta conservan el significado, no la forma: `motifs = ["CA"]` vuelve como
    `motif = "CA"`, y son el mismo locus.
    """
    loci = tomlkit.aot()
    for index, locus in enumerate(form.reference.repeat_loci):
        path = f"reference.repeat_loci[{index}]"
        table = tomlkit.table()
        table.add("id", _non_empty_string(locus.id, f"{path}.id"))
        table.add("sequence", _non_empty_string(locus.sequence, f"{path}.sequence"))
        table.add("start", _checked_positive_u64(locus.start, f"{path}.start"))
        table.add("end", _checked_positive_u64(locus.end, f"{path}.end"))
        if locus.start > locus.end:
            raise EditorError(f"{path}.start", "el inicio no puede pasar del final")
        if not locus.motifs:
            raise EditorError(f"{path}.motif", "el locus no declara ningún motivo")
        motifs = [
            _repeat_motif(motif, f"{path}.motifs[{motif_index}]")
            for motif_index, motif in enumerate(locus.motifs)
        ]
        if len(motifs) == 1:
            table.add("motif", motifs[0])
        else:
            table.add("motifs", motifs)
        if locus.max_units is not None:
            table.add(
                "max_units",
                _checked_positive_u64(locus.max_units, f"{path}.max_units"),
            )
        if locus.min_purity is not None:
            table.add(
                "min_purity", _literal_item(locus.min_purity, f"{path}.min_purity")
            )
        if locus.min_longest_run_units is not None:
            table.add(
                "min_longest_run_units",
                _checked_positive_u64(
                    locus.min_longest_run_units, f"{path}.min_longest_run_units"
                ),
            )
        if not locus.alleles:
            raise EditorError(f"{path}.alleles", "el locus no declara ningún alelo")
        # El reparto se comprueba **aquí**, antes de escribir nada: un TOML con
        # frecuencias que no suman uno no llega a existir.
        normalize_repeat_frequencies(
            [allele.frequency for allele in locus.alleles],
            form.normalization_tolerance,
            f"{path}.alleles",
        )
        alleles = tomlkit.aot()
        for allele_index, allele in enumerate(locus.alleles):
            allele_path = f"{path}.alleles[{allele_index}]"
            row = tomlkit.table()
            if (allele.units is None) == (allele.sequence is None):
                raise EditorError(
                    allele_path, "se exige exactamente units o sequence"
                )
            if allele.units is not None:
                if len(allele.units) != len(motifs):
                    raise EditorError(
                        f"{allele_path}.units",
                        f"se exige un valor por bloque: {len(motifs)} esperados,"
                        f" {len(allele.units)} declarados",
                    )
                if len(allele.units) == 1:
                    row.add(
                        "units",
                        _checked_positive_u64(
                            allele.units[0], f"{allele_path}.units"
                        ),
                    )
                else:
                    row.add(
                        "units",
                        [
                            _checked_positive_u64(
                                value, f"{allele_path}.units[{value_index}]"
                            )
                            for value_index, value in enumerate(allele.units)
                        ],
                    )
            else:
                if allele.interruptions:
                    raise EditorError(
                        f"{allele_path}.interruptions",
                        "un alelo declarado por su secuencia ya lleva sus interrupciones dentro",
                    )
                row.add(
                    "sequence",
                    _repeat_allele_sequence(
                        allele.sequence, f"{allele_path}.sequence"
                    ),
                )
            row.add(
                "frequency",
                _non_negative_literal_item(
                    allele.frequency, f"{allele_path}.frequency"
                ),
            )
            if allele.interruptions:
                interruptions = tomlkit.aot()
                for item_index, interruption in enumerate(allele.interruptions):
                    item_path = f"{allele_path}.interruptions[{item_index}]"
                    item = tomlkit.table()
                    item.add(
                        "offset",
                        _checked_non_negative_u64(
                            interruption.offset, f"{item_path}.offset"
                        ),
                    )
                    item.add("bases", _repeat_motif(interruption.bases, f"{item_path}.bases"))
                    interruptions.append(item)
                row.add("interruptions", interruptions)
            alleles.append(row)
        table.add("alleles", alleles)
        loci.append(table)
    return loci


def validate_simulated_composition(composition: Mapping[str, Any]) -> None:
    """Validate the local composition contract before calling the Rust binary.

    `N` is intentionally outside the A/C/G/T mass: it represents an independent
    fraction of ambiguous positions in the final simulated sequence.
    """
    values: list[float] = []
    for base in BASES:
        path = f"reference.composition.{base}"
        try:
            value = float(str(composition.get(base, "")).strip())
        except (TypeError, ValueError) as error:
            raise EditorError(path, "debe ser un decimal finito entre 0 y 1") from error
        if not math.isfinite(value) or not 0.0 <= value <= 1.0:
            raise EditorError(path, "debe ser un decimal finito entre 0 y 1")
        values.append(value)

    total = math.fsum(values)
    if total > 1.0 + COMPOSITION_TOLERANCE:
        raise EditorError(
            "reference.composition",
            f"A/C/G/T suman {total:g}; no pueden superar 1",
        )
    if abs(total - 1.0) > COMPOSITION_TOLERANCE:
        raise EditorError(
            "reference.composition",
            f"A/C/G/T deben sumar 1 (suman {total:g})",
        )

    if "N" not in composition:
        return
    try:
        unknown_fraction = float(str(composition["N"]).strip())
    except (TypeError, ValueError) as error:
        raise EditorError("reference.composition.N", "debe ser un decimal finito entre 0 y 1") from error
    if not math.isfinite(unknown_fraction) or not 0.0 <= unknown_fraction <= 1.0:
        raise EditorError("reference.composition.N", "debe ser un decimal finito entre 0 y 1")


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
    # tomlkit usa una subclase de `str` que Marimo no puede serializar como
    # valor inicial de sus controles. El formulario solo conserva texto plano.
    return str(value)


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


def _ploidy(value: Any, path: str = "ploidy") -> int:
    if isinstance(value, bool) or not isinstance(value, (Integer, int)):
        raise EditorError(path, "se esperaba un entero")
    return int(value)


def _checked_positive_u64(value: int, path: str) -> int:
    if isinstance(value, bool) or not 1 <= value <= U64_MAX:
        raise EditorError(path, "debe pertenecer a 1..=u64::MAX")
    return value


def _non_negative_u64(value: Any, path: str) -> int:
    if not isinstance(value, Integer):
        raise EditorError(path, "se esperaba un entero")
    return _checked_non_negative_u64(int(value), path)


def _checked_non_negative_u64(value: int | None, path: str) -> int:
    if isinstance(value, bool) or value is None or not 0 <= value <= U64_MAX:
        raise EditorError(path, "debe pertenecer a 0..=u64::MAX")
    return value


def _checked_ploidy(value: int, path: str) -> int:
    if isinstance(value, bool) or not 1 <= value <= 10:
        raise EditorError(path, "debe pertenecer a 1..=10")
    return value


def _required_text(value: str | None, path: str) -> str:
    if value is None:
        raise EditorError(path, "falta el texto")
    return value


def _number_literal(value: Any, path: str) -> str:
    if not isinstance(value, (Integer, Float)):
        raise EditorError(path, "se esperaba un número decimal")
    return value.as_string()


def _non_negative_number_literal(value: Any, path: str) -> str:
    literal = _number_literal(value, path)
    try:
        parsed = float(literal)
    except ValueError as error:
        raise EditorError(path, "debe ser un número finito mayor o igual que cero") from error
    if not math.isfinite(parsed) or parsed < 0:
        raise EditorError(path, "debe ser un número finito mayor o igual que cero")
    return literal


def _fraction_literal(value: Any, path: str) -> str:
    """Una fracción declarada: finita y dentro de `[0,1]`."""

    literal = _non_negative_number_literal(value, path)
    if float(literal) > 1:
        raise EditorError(path, "debe ser una fracción entre 0 y 1")
    return literal


def _fraction_literal_item(literal: str | None, path: str) -> Integer | Float:
    item = _non_negative_literal_item(literal, path)
    _fraction_literal(item, path)
    return item


def _non_negative_literal_item(literal: str | None, path: str) -> Integer | Float:
    item = _literal_item(literal, path)
    _non_negative_number_literal(item, path)
    return item


def _read_length(value: Any, path: str) -> int:
    return _checked_read_length(_positive_u64(value, path), path)


def _checked_read_length(value: int, path: str) -> int:
    if not 1 <= value <= 150:
        raise EditorError(path, "debe pertenecer a 1..=150")
    return value


def _udi(value: Any, path: str) -> str:
    text = _string(value, path)
    if re.fullmatch(r"[ACGT]{8}", text) is None:
        raise EditorError(path, "debe contener exactamente ocho bases A/C/G/T")
    return text


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
