"""F4.9.13: la preparación amplificada en la interfaz Marimo.

Lo que estas pruebas comprueban, y qué se perdería sin cada una:

- **La ida y la vuelta completa.** Un TOML con la tabla de amplificación se lee
  a formulario y se vuelve a escribir sin perder nada. Sin ella, una interfaz
  que mostrara los parámetros y los tirara al guardar pasaría desapercibida.
- **Que la tabla no aparece con la preparación sin PCR.** Escribirla con la
  otra preparación se rechaza en la interfaz, no solo en el lector de Rust.
- **Que la `m` que la interfaz muestra es la que el simulador publica.** Es la
  comprobación que caza una interfaz que la calculara con otro método.
- **Que los límites de los controles son las constantes de contrato.** Un
  control cuyo máximo no coincidiera dejaría al usuario declarar algo que el
  simulador rechaza después.
"""

from __future__ import annotations

import ast
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import (  # noqa: E402
    LIBRARY_PREPARATIONS,
    PCR_MAXIMUM_CYCLES,
    PCR_MAXIMUM_EARLY_CYCLES,
    PCR_MAXIMUM_EFFICIENCY_CLASSES,
    PCR_MAXIMUM_FAMILY_SIZE,
    PCR_MAXIMUM_GC_ATTEMPTS,
    EditorError,
    EfficiencyClassRow,
    GcCurveKnotRow,
    PcrForm,
    efficiency_class_editor_rows,
    gc_curve_editor_rows,
    parse_general_toml,
    pcr_from_editor_rows,
    pcr_mean_pairs_per_molecule,
    render_general_toml,
)
from process_adapter import SUPPORTED_BINARY_NAME, run_toml, validate_toml  # noqa: E402

BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME
NOTEBOOK = ROOT / "ui" / "gps_app.py"
PANELS = ROOT / "ui" / "ui_panels.py"

SEED = "00" * 32


def _reference(length: int) -> str:
    """Bases equilibradas y sin periodo corto."""

    state = 0x9E3779B97F4A7C15
    out = []
    for _ in range(length):
        state = (state * 6364136223846793005 + 1442695040888963407) % (1 << 64)
        out.append("ACGT"[(state >> 33) % 4])
    return "".join(out)


REFERENCE = _reference(8000)

SCENARIO = f"""config_kind = "general"
schema_version = "4"
seed = "{SEED}"
ploidy = 1

[population_snvs]
mode = "supplied"
input_path = "population-snvs.tsv"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 8000

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
# La interfaz no edita `sequencing_errors`, y aquí no hace falta: esta
# prueba cuenta parejas y lee la `m` publicada, no compara duplicados
# byte a byte.
target_coverage_wgs = 6
library_preparation = "pcr"

[technical_event.pcr]
duplicate_fraction = 0.2
cycles = 12
efficiency = 0.92
early_cycles = 4
max_family_size = 32
gc_max_attempts = 5

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.0
relative_efficiency = 0.4

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.5
relative_efficiency = 1.0

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 1.0
relative_efficiency = 0.4

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.95
weight = 0.7

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.6
weight = 0.3

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"""

RUNNABLE = SCENARIO.replace(
    """
[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.0
relative_efficiency = 0.4

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 0.5
relative_efficiency = 1.0

[[technical_event.pcr.gc_relative_efficiency]]
gc_fraction = 1.0
relative_efficiency = 0.4

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.95
weight = 0.7

[[technical_event.pcr.efficiency_classes]]
efficiency = 0.6
weight = 0.3
""",
    "",
)


class PcrEditorTest(unittest.TestCase):
    def test_round_trip_keeps_every_amplification_parameter(self) -> None:
        form = parse_general_toml(SCENARIO)
        pcr = form.technical_event.pcr
        self.assertIsNotNone(pcr)
        self.assertEqual(form.technical_event.library_preparation, "pcr")
        self.assertEqual(pcr.duplicate_fraction, "0.2")
        self.assertEqual(pcr.cycles, 12)
        self.assertEqual(pcr.efficiency, "0.92")
        self.assertEqual(pcr.early_cycles, 4)
        self.assertEqual(pcr.max_family_size, 32)
        self.assertEqual(pcr.gc_max_attempts, 5)
        self.assertEqual(len(pcr.gc_relative_efficiency), 3)
        self.assertEqual(len(pcr.efficiency_classes), 2)
        # Y la vuelta conserva el formulario entero, no solo los campos de
        # arriba: una comparación campo a campo dejaría pasar el que se
        # olvidara el día que alguien añada uno.
        self.assertEqual(parse_general_toml(render_general_toml(form)), form)

    def test_the_amplification_table_needs_the_amplified_preparation(self) -> None:
        """La tabla con la preparación sin PCR se rechaza **en la interfaz**."""

        form = parse_general_toml(SCENARIO)
        form.technical_event.library_preparation = "pcr_free"
        with self.assertRaises(EditorError) as raised:
            render_general_toml(form)
        self.assertIn("technical_event.pcr", str(raised.exception))
        # Control positivo: con la preparación amplificada el mismo formulario
        # se escribe sin problemas.
        form.technical_event.library_preparation = "pcr"
        self.assertIn("[technical_event.pcr]", render_general_toml(form))

    def test_the_amplified_preparation_needs_its_table(self) -> None:
        form = parse_general_toml(SCENARIO)
        form.technical_event.pcr = None
        with self.assertRaises(EditorError) as raised:
            render_general_toml(form)
        self.assertIn("obligator", str(raised.exception))

    def test_the_editor_rows_do_not_invent_an_amplification(self) -> None:
        """Con preparación sin PCR no se reconstruye ninguna tabla.

        La variante defectuosa sería devolver un bloque cada vez que el usuario
        escribe algo en las tablas, aunque la preparación no lo admita: todo
        escenario exportado saldría con amplificación.
        """

        rows = gc_curve_editor_rows(None)
        classes = efficiency_class_editor_rows(None)
        self.assertIsNone(
            pcr_from_editor_rows(
                "pcr_free", "0.2", 10, "0.9", 3, 64, 8, rows, classes
            )
        )
        rebuilt = pcr_from_editor_rows(
            "pcr", "0.2", 10, "0.9", 3, 64, 8, rows, classes
        )
        self.assertIsNotNone(rebuilt)
        # Las filas neutras de la curva **no se declaran**: son el valor de
        # omisión, y escribirlas cambiaría los bytes sin cambiar el resultado.
        self.assertEqual(rebuilt.gc_relative_efficiency, [])
        self.assertEqual(rebuilt.efficiency_classes, [])

    def test_the_editor_rows_carry_a_declared_curve_and_classes(self) -> None:
        """Control positivo del anterior: lo declarado **sí** viaja."""

        form = parse_general_toml(SCENARIO)
        pcr = form.technical_event.pcr
        rebuilt = pcr_from_editor_rows(
            "pcr",
            pcr.duplicate_fraction,
            pcr.cycles,
            pcr.efficiency,
            pcr.early_cycles,
            pcr.max_family_size,
            pcr.gc_max_attempts,
            gc_curve_editor_rows(pcr),
            efficiency_class_editor_rows(pcr),
        )
        self.assertEqual(rebuilt, pcr)

    def test_a_curve_with_a_single_knot_is_rejected(self) -> None:
        form = parse_general_toml(SCENARIO)
        form.technical_event.pcr.gc_relative_efficiency = [
            GcCurveKnotRow("0.0", "0.4")
        ]
        with self.assertRaises(EditorError) as raised:
            render_general_toml(form)
        self.assertIn("al menos dos nudos", str(raised.exception))

    def test_more_efficiency_classes_than_the_contract_allows_are_rejected(
        self,
    ) -> None:
        form = parse_general_toml(SCENARIO)
        form.technical_event.pcr.efficiency_classes = [
            EfficiencyClassRow("0.9", "0.2")
        ] * (PCR_MAXIMUM_EFFICIENCY_CLASSES + 1)
        with self.assertRaises(EditorError) as raised:
            render_general_toml(form)
        self.assertIn(str(PCR_MAXIMUM_EFFICIENCY_CLASSES), str(raised.exception))
        # Control positivo: justo en el tope se admite.
        form.technical_event.pcr.efficiency_classes = [
            EfficiencyClassRow("0.9", "0.25")
        ] * PCR_MAXIMUM_EFFICIENCY_CLASSES
        self.assertIn("efficiency_classes", render_general_toml(form))

    def test_a_count_outside_its_contract_domain_is_rejected(self) -> None:
        for field, maximum in (
            ("cycles", PCR_MAXIMUM_CYCLES),
            ("early_cycles", PCR_MAXIMUM_EARLY_CYCLES),
            ("max_family_size", PCR_MAXIMUM_FAMILY_SIZE),
            ("gc_max_attempts", PCR_MAXIMUM_GC_ATTEMPTS),
        ):
            with self.subTest(field=field):
                document = SCENARIO.replace(
                    f"{field} = ", f"{field} = {maximum + 1}  # ", 1
                )
                with self.assertRaises(EditorError):
                    parse_general_toml(document)

    def test_an_invented_preparation_is_rejected(self) -> None:
        document = SCENARIO.replace(
            'library_preparation = "pcr"', 'library_preparation = "nanopore"'
        )
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(document)
        self.assertIn("library_preparation", str(raised.exception))
        self.assertIn(LIBRARY_PREPARATIONS[0], str(raised.exception))


class DerivedMeanPairsTest(unittest.TestCase):
    """La `m` que la interfaz muestra."""

    def test_a_zero_fraction_derives_exactly_zero(self) -> None:
        """**Cero exacto, no un residuo.**

        Es el discriminador más fino de este fichero: una bisección sin retorno
        anticipado devuelve un residuo del orden de `10⁻¹⁸`, que no da ningún
        error y cambiaría todos los FASTQ de una ejecución con amplificación
        neutra.
        """

        self.assertEqual(pcr_mean_pairs_per_molecule("0"), 0.0)
        self.assertEqual(pcr_mean_pairs_per_molecule("0.0"), 0.0)
        # Control positivo: una fracción distinta de cero sí da un número
        # positivo, de modo que la igualdad de arriba no pasa por devolver
        # siempre cero.
        self.assertGreater(pcr_mean_pairs_per_molecule("0.05"), 0.0)

    def test_an_out_of_domain_fraction_derives_nothing(self) -> None:
        for value in ("-0.1", "1", "1.5", "", "abc", "0.99"):
            with self.subTest(value=value):
                self.assertIsNone(pcr_mean_pairs_per_molecule(value))

    def test_the_derivation_inverts_the_worked_fractions(self) -> None:
        """Las tres parejas que el diseño trae calculadas en el sentido directo."""

        for expected, fraction in (
            (0.02, "0.009933665337762632"),
            (0.1, "0.048374180359595176"),
            (0.5, "0.21306131942526685"),
        ):
            with self.subTest(fraction=fraction):
                derived = pcr_mean_pairs_per_molecule(fraction)
                self.assertLess(abs(derived - expected) / expected, 1e-12)

    @unittest.skipUnless(
        BINARY.exists(), "compile the Rust binary before this integration test"
    )
    def test_the_shown_mean_pairs_match_the_ones_the_simulator_publishes(self) -> None:
        """**La comprobación que caza una interfaz con otro método.**

        Se compara con una tolerancia y no con igualdad exacta a propósito: la
        interfaz usa la exponencial del sistema y el simulador la de la
        librería fijada, y el contrato solo obliga a la segunda. Lo que la
        comparación caza es un **método distinto**, que se desviaría muchísimo
        más que esto.
        """

        toml = render_general_toml(parse_general_toml(RUNNABLE))
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "reference.fa").write_text(
                f">chr1\n{REFERENCE}\n", encoding="utf-8"
            )
            (base / "population-snvs.tsv").write_text(
                "contig\tposition\tref\talt\tpopulation_id\talt_af\n"
                f"chr1\t20\t{REFERENCE[19]}\t{'G' if REFERENCE[19] == 'T' else 'T'}\tP1\t0.5\n",
                encoding="utf-8",
            )
            validation, session = validate_toml(BINARY, toml, base)
            self.assertTrue(validation.succeeded, validation.stdout or validation.stderr)
            self.assertIsNotNone(session)
            execution = run_toml(session, toml, BINARY, base, base / "result")
            self.assertTrue(execution.succeeded, execution.stdout or execution.stderr)
            result = base / "result"
            effective = json.loads(
                (result / "effective-config.json").read_text(encoding="utf-8")
            )
            manifest = json.loads(
                (result / "events" / "reads-main" / "manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            reads = gzip.open(
                result / "outputs" / "fastq" / "reads-main" / "sample-A_R1.fastq.gz",
                "rt",
                encoding="ascii",
            ).read().splitlines()[1::4]

        published = effective["technical_event"]["artifacts"]["pcr"]
        self.assertEqual(published["applicable"], True)
        declared = published["declared_duplicate_fraction"]
        shown = pcr_mean_pairs_per_molecule(str(declared))
        self.assertIsNotNone(shown)
        self.assertLess(
            abs(shown - published["mean_pairs_per_molecule"])
            / published["mean_pairs_per_molecule"],
            1e-12,
            "la m de la interfaz no es la que publica la configuración efectiva",
        )
        # Control positivo de la ejecución: el escenario que la interfaz
        # escribió produjo de verdad duplicados, de modo que la comparación no
        # pasa por haber ejecutado algo vacío.
        pcr = manifest["artifacts"]["pcr"]
        self.assertGreater(pcr["duplicate_pairs"], 0)
        self.assertEqual(len(reads), pcr["pairs"])


class NotebookControlsTest(unittest.TestCase):
    """Los controles existen y llevan los límites del contrato."""

    def test_the_notebook_offers_a_control_for_every_amplification_parameter(
        self,
    ) -> None:
        notebook = NOTEBOOK.read_text(encoding="utf-8")
        panels = PANELS.read_text(encoding="utf-8")
        for control in (
            "library_preparation_control",
            "pcr_duplicate_fraction_control",
            "pcr_cycles_control",
            "pcr_efficiency_control",
            "pcr_early_cycles_control",
            "pcr_max_family_size_control",
            "pcr_gc_max_attempts_control",
            "pcr_gc_curve_editor",
            "pcr_efficiency_classes_editor",
        ):
            with self.subTest(control=control):
                self.assertIn(control, notebook)
                self.assertIn(control, panels)
        self.assertIn("pcr_panel", notebook)
        # Control negativo con su control positivo: el panel avisa de que la
        # fracción declarada es la esperada y no la realizada, y la búsqueda
        # encuentra de verdad lo que tiene que estar.
        self.assertIn("ESPERADA", panels)
        self.assertIn("MAYOR", panels)

    def test_the_controls_take_their_limits_from_the_contract_constants(self) -> None:
        """**Los límites son las constantes, no cifras escritas a mano.**

        Se comprueba sobre el árbol sintáctico: un control cuyo `stop` fuera un
        literal decimal pasaría una búsqueda de texto y moriría aquí.
        """

        tree = ast.parse(NOTEBOOK.read_text(encoding="utf-8"))
        expected = {
            "pcr_cycles_control": "PCR_MAXIMUM_CYCLES",
            "pcr_early_cycles_control": "PCR_MAXIMUM_EARLY_CYCLES",
            "pcr_max_family_size_control": "PCR_MAXIMUM_FAMILY_SIZE",
            "pcr_gc_max_attempts_control": "PCR_MAXIMUM_GC_ATTEMPTS",
        }
        found: dict[str, str] = {}
        for node in ast.walk(tree):
            if not isinstance(node, ast.Assign):
                continue
            target = node.targets[0]
            if not isinstance(target, ast.Name) or target.id not in expected:
                continue
            self.assertIsInstance(node.value, ast.Call)
            stop = [
                keyword.value
                for keyword in node.value.keywords
                if keyword.arg == "stop"
            ]
            self.assertEqual(len(stop), 1, target.id)
            self.assertIsInstance(stop[0], ast.Name, target.id)
            found[target.id] = stop[0].id
        self.assertEqual(found, expected)


if __name__ == "__main__":
    unittest.main()
