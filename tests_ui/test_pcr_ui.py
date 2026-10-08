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
    HOMOPOLYMER_PCR_FOLD_MAXIMUM,
    HOMOPOLYMER_PCR_FOLD_MINIMUM,
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
    PerCycleSlipRow,
    efficiency_class_editor_rows,
    gc_curve_editor_rows,
    parse_general_toml,
    pcr_from_editor_rows,
    pcr_mean_pairs_per_molecule,
    per_cycle_slip_editor_rows,
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
        librería fijada, y el contrato solo obliga a la segunda. **Se midió**:
        de sesenta fracciones barridas, dieciocho difieren en los últimos bits.
        Lo que la comparación caza es un **algoritmo distinto**, que se
        desviaría muchísimo más que esto. La autoridad es el número publicado.
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


class ClonalStutterEditorTest(unittest.TestCase):
    """F4.10.6: el control del tartamudeo clonal en la interfaz.

    Lo que se declara aquí es la probabilidad de **un ciclo**, no la fracción
    de extremo a extremo. La interfaz no la deriva ni la comprueba: la
    autoridad es Rust. Lo que sí tiene que hacer es **conservarla** en la ida y
    la vuelta y no inventarla cuando no se declaró.
    """

    def _scenario_with_slips(self) -> str:
        return SCENARIO.replace(
            "[[technical_event.libraries]]",
            "[[technical_event.pcr.stutter.per_cycle_slip]]\n"
            "motif_length_bases = 1\n"
            "probability = 0.0158\n"
            "\n"
            "[[technical_event.pcr.stutter.per_cycle_slip]]\n"
            "motif_length_bases = 2\n"
            "probability = 0.00276\n"
            "\n"
            "[[technical_event.libraries]]",
            1,
        )

    def test_round_trip_keeps_the_declared_per_cycle_slips(self) -> None:
        form = parse_general_toml(self._scenario_with_slips())
        pcr = form.technical_event.pcr
        self.assertEqual(
            pcr.per_cycle_slip,
            [PerCycleSlipRow(1, "0.0158"), PerCycleSlipRow(2, "0.00276")],
        )
        # La vuelta conserva el formulario entero, no solo este campo.
        self.assertEqual(parse_general_toml(render_general_toml(form)), form)
        # Y el literal declarado sobrevive tal cual: un `0.00276` que volviera
        # como `0.00276000000001` cambiaría los bytes del documento.
        self.assertIn("probability = 0.00276", render_general_toml(form))

    def test_an_undeclared_table_is_not_written(self) -> None:
        """Sin sustituciones, la tabla **no se escribe**.

        Escribir las seis derivadas cambiaría los bytes de todo documento ya
        publicado y, peor, el manifiesto las publicaría marcadas como del
        usuario en vez de como despejadas.
        """

        form = parse_general_toml(SCENARIO)
        self.assertEqual(form.technical_event.pcr.per_cycle_slip, [])
        self.assertNotIn("per_cycle_slip", render_general_toml(form))
        # Control positivo: con una sustitución declarada sí se escribe.
        form.technical_event.pcr.per_cycle_slip = [PerCycleSlipRow(2, "0.01")]
        self.assertIn(
            "[[technical_event.pcr.stutter.per_cycle_slip]]", render_general_toml(form)
        )

    def test_the_editor_rows_carry_the_declared_slips(self) -> None:
        form = parse_general_toml(self._scenario_with_slips())
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
            per_cycle_slip_editor_rows(pcr),
        )
        self.assertEqual(rebuilt, pcr)

    def test_the_editor_rows_do_not_invent_a_substitution(self) -> None:
        """Una fila en blanco no es una sustitución.

        La variante defectuosa sería tomar la fila vacía que el editor presenta
        como una probabilidad declarada de cero, que apagaría el tartamudeo
        clonal de esa longitud de motivo sin que nadie lo pidiera.
        """

        rows = per_cycle_slip_editor_rows(None)
        self.assertEqual(rows, [{"order": 1, "motif_length_bases": 1, "probability": ""}])
        rebuilt = pcr_from_editor_rows(
            "pcr",
            "0.2",
            10,
            "0.9",
            3,
            64,
            8,
            gc_curve_editor_rows(None),
            efficiency_class_editor_rows(None),
            rows,
        )
        self.assertEqual(rebuilt.per_cycle_slip, [])

    def test_a_motif_length_outside_the_calibrated_range_is_rejected(self) -> None:
        scenario = self._scenario_with_slips().replace(
            "motif_length_bases = 1", "motif_length_bases = 7", 1
        )
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(scenario)
        self.assertIn("motif_length_bases", str(raised.exception))

    def test_a_repeated_motif_length_is_rejected(self) -> None:
        scenario = self._scenario_with_slips().replace(
            "motif_length_bases = 1", "motif_length_bases = 2", 1
        )
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(scenario)
        self.assertIn("una sola probabilidad por ciclo", str(raised.exception))

    def test_the_panel_shows_the_control(self) -> None:
        """El control llega al panel, no solo al formulario.

        Sin esta comprobación, un control construido y nunca colocado pasaría
        desapercibido: es exactamente lo que le ocurrió al tartamudeo de base.
        """

        notebook = NOTEBOOK.read_text(encoding="utf-8")
        panels = PANELS.read_text(encoding="utf-8")
        self.assertIn("pcr_per_cycle_slip_editor = mo.ui.data_editor", notebook)
        self.assertIn("pcr_per_cycle_slip_editor,", panels)
        self.assertIn("pcr_per_cycle_slip_editor.value", notebook)



class HomopolymerFoldEditorTest(unittest.TestCase):
    """DTC-099: el multiplicador de la PCR en los homopolímeros en la interfaz.

    Antes de esta entrega el editor **rechazaba** un documento con la clave,
    porque la tabla del tartamudeo clonal solo admitía `per_cycle_slip`. La
    autoridad del dominio sigue siendo Rust; lo que la interfaz tiene que hacer
    es conservar el literal en la ida y la vuelta y no escribirlo si no se
    declaró.
    """

    def _scenario_with_fold(self, literal: str = "2.5") -> str:
        return SCENARIO.replace(
            "[[technical_event.libraries]]",
            "[technical_event.pcr.stutter]\n"
            f"homopolymer_pcr_fold = {literal}\n"
            "\n"
            "[[technical_event.libraries]]",
            1,
        )

    def test_round_trip_keeps_the_declared_fold(self) -> None:
        form = parse_general_toml(self._scenario_with_fold())
        self.assertEqual(form.technical_event.pcr.homopolymer_pcr_fold, "2.5")
        rendered = render_general_toml(form)
        self.assertIn("homopolymer_pcr_fold = 2.5", rendered)
        self.assertEqual(parse_general_toml(rendered), form)

    def test_an_undeclared_fold_is_not_written(self) -> None:
        """Sin declarar, la clave **no se escribe**: 6.30 no equivale a omitirla."""

        form = parse_general_toml(SCENARIO)
        self.assertIsNone(form.technical_event.pcr.homopolymer_pcr_fold)
        self.assertNotIn("homopolymer_pcr_fold", render_general_toml(form))
        self.assertNotIn("[technical_event.pcr.stutter]", render_general_toml(form))
        # Control positivo: declarado, sí se escribe, y sin sustituciones.
        form.technical_event.pcr.homopolymer_pcr_fold = "1"
        rendered = render_general_toml(form)
        self.assertIn("homopolymer_pcr_fold = 1", rendered)
        self.assertNotIn("per_cycle_slip", rendered)

    def test_the_fold_outside_its_domain_is_rejected(self) -> None:
        self.assertEqual(
            (HOMOPOLYMER_PCR_FOLD_MINIMUM, HOMOPOLYMER_PCR_FOLD_MAXIMUM), (1.0, 37.0)
        )
        for literal in ("0.99", "37.01", "nan", "inf"):
            with self.subTest(literal=literal):
                with self.assertRaises(EditorError) as raised:
                    parse_general_toml(self._scenario_with_fold(literal))
                self.assertIn("homopolymer_pcr_fold", str(raised.exception))
        # Control positivo: los dos extremos se admiten.
        for literal in ("1", "37"):
            with self.subTest(literal=literal):
                parse_general_toml(self._scenario_with_fold(literal))

    def test_the_editor_rows_carry_the_fold_and_a_blank_is_absent(self) -> None:
        form = parse_general_toml(self._scenario_with_fold())
        pcr = form.technical_event.pcr
        arguments = (
            "pcr",
            pcr.duplicate_fraction,
            pcr.cycles,
            pcr.efficiency,
            pcr.early_cycles,
            pcr.max_family_size,
            pcr.gc_max_attempts,
            gc_curve_editor_rows(pcr),
            efficiency_class_editor_rows(pcr),
            per_cycle_slip_editor_rows(pcr),
        )
        self.assertEqual(pcr_from_editor_rows(*arguments, " 2.5 "), pcr)
        self.assertIsNone(pcr_from_editor_rows(*arguments, "  ").homopolymer_pcr_fold)

    def test_the_panel_shows_the_control(self) -> None:
        panels = PANELS.read_text(encoding="utf-8")
        self.assertIn("pcr_homopolymer_fold_control,", panels)

    def test_the_control_value_reaches_the_form(self) -> None:
        """**El valor del control es lo que llega al formulario.**

        Lo honesto sería ejecutar la celda del cuaderno que construye el
        formulario y ver el campo en el TOML. **No se puede**: con un evento
        técnico activo esa celda falla con `NameError` sobre
        `TechnicalEventForm`, que el cuaderno nunca importa (DT-042), de modo
        que hoy ningún control del panel de PCR llega al TOML desde el cuaderno
        en ejecución. Mientras eso siga así, la prueba es **estructural**: sobre
        el árbol sintáctico exige que la única llamada a `pcr_from_editor_rows`
        reciba como argumento del multiplicador exactamente
        `str(pcr_homopolymer_fold_control.value)`, y que el control arranque
        con el literal del formulario importado. Caza un cableado que pase `""`
        o cualquier otra cosa en su lugar; no caza nada que solo se vea al
        ejecutar. Después comprueba, ya sin el cuaderno, que ese argumento llega
        al formulario y al TOML.
        """

        tree = ast.parse(NOTEBOOK.read_text(encoding="utf-8"))
        calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "pcr_from_editor_rows"
        ]
        self.assertEqual(len(calls), 1)
        self.assertEqual(len(calls[0].args), 11)
        self.assertEqual(
            ast.dump(calls[0].args[10]),
            ast.dump(
                ast.parse("str(pcr_homopolymer_fold_control.value)", mode="eval").body
            ),
        )
        defined = [
            node.value
            for node in ast.walk(tree)
            if isinstance(node, ast.Assign)
            and isinstance(node.targets[0], ast.Name)
            and node.targets[0].id == "pcr_homopolymer_fold_control"
        ]
        self.assertEqual(len(defined), 1)
        initial = [
            keyword.value for keyword in defined[0].keywords if keyword.arg == "value"
        ]
        self.assertEqual(len(initial), 1)
        self.assertIn("_selected_pcr.homopolymer_pcr_fold", ast.unparse(initial[0]))
        # Y lo que ese argumento lleva llega al formulario y al TOML.
        pcr = parse_general_toml(SCENARIO).technical_event.pcr
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
            per_cycle_slip_editor_rows(pcr),
            "2.82",
        )
        self.assertEqual(rebuilt.homopolymer_pcr_fold, "2.82")
        form = parse_general_toml(SCENARIO)
        form.technical_event.pcr = rebuilt
        self.assertIn("homopolymer_pcr_fold = 2.82", render_general_toml(form))


if __name__ == "__main__":
    unittest.main()
