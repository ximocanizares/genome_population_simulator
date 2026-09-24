"""F4.8.9: el tartamudeo en la interfaz Marimo.

Lo que estas pruebas comprueban, y qué se perdería sin cada una:

- **La ida y la vuelta completa.** Un TOML con el bloque de tartamudeo se lee a
  formulario y se vuelve a escribir sin perder nada. Sin ella, una interfaz que
  mostrara los parámetros y los tirara al guardar pasaría desapercibida.
- **Que no hay control para lo que el algoritmo no admite.** Un paso `+2` o una
  longitud de motivo fuera del tramo calibrado se rechazan **en la interfaz**,
  no solo en el lector de Rust.
- **Lo que la contradicción C-6 motiva:** partiendo solo de la interfaz se puede
  construir un escenario completo —referencia, locus, alelos, frecuencias y
  tasa— que produce lecturas tartamudeadas. Una interfaz con parámetros de
  tartamudeo y sin forma de declarar el locus no serviría de nada, y esa es la
  variante defectuosa que la última prueba caza.
"""

from __future__ import annotations

import ast
import dataclasses
import gzip
import importlib
import json
import re
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import (  # noqa: E402
    STUTTER_ADMITTED_STEPS,
    STUTTER_CALIBRATED_MOTIF_LENGTH_RATES,
    EditorError,
    TechnicalEventForm,
    parse_general_toml,
    render_general_toml,
    stutter_from_editor_rows,
    stutter_locus_rate_editor_rows,
    stutter_motif_rate_editor_rows,
    stutter_step_editor_rows,
)
from process_adapter import SUPPORTED_BINARY_NAME, run_toml, validate_toml  # noqa: E402

BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

SEED = "00" * 32

SCENARIO = f"""config_kind = "general"
schema_version = "4"
seed = "{SEED}"

[reference]
kind = "supplied"
id = "ref-main"
fasta = "reference.fa"

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 400

[[reference.repeat_loci]]
id = "STR-1"
sequence = "chr1"
start = 101
end = 130
motif = "CA"

[[reference.repeat_loci.alleles]]
units = 15
frequency = 0.6

[[reference.repeat_loci.alleles]]
units = 12
frequency = 0.4

[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 1
genotype_model = "hwe"
individual_ids = ["sample-A"]

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 30

[technical_event.stutter]
enabled = true
run_length_coefficient = 0.0

[[technical_event.stutter.motif_length_rates]]
motif_length_bases = 2
rate = 0.5

[[technical_event.stutter.steps]]
units = -3
weight = 0.05

[[technical_event.stutter.steps]]
units = -2
weight = 0.15

[[technical_event.stutter.steps]]
units = -1
weight = 0.6

[[technical_event.stutter.steps]]
units = 1
weight = 0.2

[[technical_event.stutter.locus_rates]]
locus_id = "STR-1"
rate = 0.25

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"""


# Cuatrocientas bases: cien de flanco, treinta de tracto `(CA)15` y el resto de
# flanco. El tracto queda holgadamente dentro de una lectura de 150.
RUNNABLE_REFERENCE = "G" * 100 + "CA" * 15 + "G" * 270

# El mismo escenario, con la matriz que exige tener población inicial, la tasa
# a uno y un solo paso con peso. La sustitución por locus también va a uno, para
# que lo que gobierna sea inequívoco.
RUNNABLE_SCENARIO = (
    SCENARIO.replace('fasta = "reference.fa"', 'fasta = "reference.fa"')
    .replace("length = 400", "length = 400")
    .replace("rate = 0.5", "rate = 1.0")
    .replace('locus_id = "STR-1"\nrate = 0.25', 'locus_id = "STR-1"\nrate = 1.0')
    .replace("units = -1\nweight = 0.6", "units = -1\nweight = 1.0")
    .replace("units = -3\nweight = 0.05", "units = -3\nweight = 0")
    .replace("units = -2\nweight = 0.15", "units = -2\nweight = 0")
    .replace("units = 1\nweight = 0.2", "units = 1\nweight = 0")
    .replace(
        'seed = "{seed}"'.format(seed=SEED),
        'seed = "{seed}"\n\n[population_snvs]\nmode = "supplied"\n'
        'input_path = "population-snvs.tsv"'.format(seed=SEED),
    )
)


# **Cada sustitución tiene que haber sustituido.** Una llamada a `replace` que
# no encuentra su patrón devuelve la cadena intacta y en silencio: es la misma
# trampa que una mutación que no muta, y aquí dejaría el escenario con la tasa
# de partida mientras la prueba creería estar ejecutando otra cosa.
assert "rate = 1.0" in RUNNABLE_SCENARIO
assert "rate = 0.25" not in RUNNABLE_SCENARIO
assert "rate = 0.5" not in RUNNABLE_SCENARIO
assert "units = -1\nweight = 1.0" in RUNNABLE_SCENARIO
assert "units = 1\nweight = 0\n" in RUNNABLE_SCENARIO
assert "[population_snvs]" in RUNNABLE_SCENARIO
assert RUNNABLE_SCENARIO != SCENARIO


def _read_sequences(path: Path) -> list[str]:
    """Las secuencias de un FASTQ publicado, sin tocar nada más."""

    with gzip.open(path, "rt", encoding="ascii") as handle:
        return handle.read().splitlines()[1::4]


def _reconstructed_through_the_notebook(form):
    """El formulario tal y como el cuaderno lo devuelve tras sus tablas.

    Recorre exactamente las funciones que `ui/gps_app.py` llama, de modo que lo
    que aquí se comprueba es el camino real y no una simulación suya.
    """

    event = form.technical_event
    return dataclasses.replace(
        form,
        technical_event=dataclasses.replace(
            event,
            stutter=stutter_from_editor_rows(
                bool(event.stutter and event.stutter.enabled),
                (event.stutter.run_length_coefficient if event.stutter else "") or "",
                stutter_motif_rate_editor_rows(event.stutter),
                stutter_step_editor_rows(event.stutter),
                stutter_locus_rate_editor_rows(event.stutter),
            ),
        ),
    )


def _fully_visible_runs(read: str) -> list[int]:
    """Unidades de cada tirada de `CA` **flanqueada por `G` a los dos lados**.

    El flanqueo es lo que garantiza que la tirada se ve entera: una tirada que
    toque el borde de la lectura estaría cortada y su cuenta no significaría
    nada.
    """

    return [(len(match.group(0)) - 2) // 2 for match in re.finditer(r"G(?:CA)+G", read)]


class StutterEditorTest(unittest.TestCase):
    def test_round_trip_keeps_every_stutter_parameter(self) -> None:
        form = parse_general_toml(SCENARIO)
        stutter = form.technical_event.stutter
        self.assertIsNotNone(stutter)
        self.assertTrue(stutter.enabled)
        self.assertEqual(stutter.run_length_coefficient, "0.0")
        self.assertEqual(
            [(row.motif_length_bases, row.rate) for row in stutter.motif_length_rates],
            [(2, "0.5")],
        )
        self.assertEqual(
            [row.units for row in stutter.steps], list(STUTTER_ADMITTED_STEPS)
        )
        self.assertEqual(
            [(row.locus_id, row.rate) for row in stutter.locus_rates],
            [("STR-1", "0.25")],
        )
        rendered = render_general_toml(form)
        self.assertEqual(parse_general_toml(rendered), form)

    def test_a_plus_two_step_has_no_control_and_is_rejected(self) -> None:
        document = SCENARIO.replace("units = 1\nweight = 0.2", "units = 2\nweight = 0.2")
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(document)
        self.assertIn("units", raised.exception.path)

    def test_a_motif_length_outside_the_calibration_is_rejected(self) -> None:
        document = SCENARIO.replace("motif_length_bases = 2", "motif_length_bases = 7")
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(document)
        self.assertIn("motif_length_bases", raised.exception.path)

    def test_enabled_stutter_without_the_four_weights_is_rejected(self) -> None:
        """Los pesos no tienen valor por omisión: faltar uno no se completa solo."""

        document = SCENARIO.replace(
            "[[technical_event.stutter.steps]]\nunits = -3\nweight = 0.05\n\n", ""
        )
        with self.assertRaises(EditorError) as raised:
            parse_general_toml(document)
        self.assertIn("steps", raised.exception.path)

    def test_a_rate_above_one_is_rejected(self) -> None:
        document = SCENARIO.replace("rate = 0.5", "rate = 1.5")
        with self.assertRaises(EditorError):
            parse_general_toml(document)

    def test_the_interface_alone_declares_every_piece_of_the_scenario(self) -> None:
        """El documento que la interfaz produce lleva las cinco piezas.

        Referencia, locus con su motivo, alelos, frecuencias y tasa. Una
        interfaz con parámetros de tartamudeo y sin forma de declarar el locus
        produciría un documento sin `reference.repeat_loci`.

        **Esto es una comprobación de texto y nada más**, y por eso no basta
        para cerrar la contradicción C-6: pasaría igual si Rust rechazara el
        documento. Quien cierra C-6 es la prueba de integración de abajo.
        """

        form = parse_general_toml(SCENARIO)
        rendered = render_general_toml(form)
        self.assertIn("[[reference.repeat_loci]]", rendered)
        self.assertIn('motif = "CA"', rendered)
        self.assertIn("[[reference.repeat_loci.alleles]]", rendered)
        self.assertIn("frequency", rendered)
        self.assertIn("[technical_event.stutter]", rendered)
        self.assertIn("enabled = true", rendered)
        locus_ids = [locus.id for locus in form.reference.repeat_loci]
        self.assertEqual(locus_ids, ["STR-1"])
        self.assertEqual(
            [row.locus_id for row in form.technical_event.stutter.locus_rates], locus_ids
        )

    @unittest.skipUnless(
        BINARY.exists(), "compile the Rust binary before this integration test"
    )
    def test_the_interface_alone_produces_stuttered_reads(self) -> None:
        """**La comprobación que cierra la contradicción C-6.**

        Partiendo solo del **modelo de autoría** se construye un escenario
        completo —referencia, locus, alelos, frecuencias y tasa— y se
        **ejecuta**: el documento que ese modelo escribe pasa por el lector de
        Rust y produce lecturas tartamudeadas de verdad.

        **Lo que esta prueba no toca, y hay que decirlo:** el cuaderno. Recorre
        `ui/config_editor.py` y nada más, de modo que **no puede fallar** porque
        `ui/gps_app.py` se deje el bloque fuera; de hecho no falló cuando se lo
        dejaba. Esa mitad la cubre `NotebookReconstructionTest`, más abajo.

        La tasa está a uno y el único paso con peso es el de menos una unidad,
        de modo que **toda** molécula pierde exactamente una unidad. Por eso la
        aserción puede ser tan dura: entre los tractos que se ven enteros en una
        lectura —flanqueados por `G` a los dos lados— **no puede aparecer
        ninguno del tamaño declarado**, ni 15 ni 12 unidades. Si el documento no
        llegara al simulador, o llegara sin el bloque de tartamudeo, aparecerían
        y la prueba moriría.
        """

        document = parse_general_toml(RUNNABLE_SCENARIO)
        toml = render_general_toml(document)
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "reference.fa").write_text(
                f">chr1\n{RUNNABLE_REFERENCE}\n", encoding="utf-8"
            )
            (base / "population-snvs.tsv").write_text(
                "contig\tposition\tref\talt\tpopulation_id\talt_af\n"
                "chr1\t20\tG\tT\tP1\t0.5\n",
                encoding="utf-8",
            )
            validation, session = validate_toml(BINARY, toml, base)
            self.assertTrue(validation.succeeded, validation.stdout or validation.stderr)
            self.assertIsNotNone(session)
            execution = run_toml(session, toml, BINARY, base, base / "result")
            self.assertTrue(execution.succeeded, execution.stdout or execution.stderr)
            result = base / "result"
            manifest = json.loads(
                (result / "events" / "reads-main" / "manifest.json").read_text(
                    encoding="utf-8"
                )
            )
            reads = _read_sequences(
                result / "outputs" / "fastq" / "reads-main" / "sample-A_R1.fastq.gz"
            )

        stutter = manifest["artifacts"]["stutter"]
        self.assertTrue(stutter["enabled"])
        self.assertGreater(stutter["stuttered_tracts"], 0)

        visible = [units for read in reads for units in _fully_visible_runs(read)]
        self.assertTrue(visible, "alguna lectura tiene que ver el tracto entero")
        self.assertNotIn(15, visible, "ningún tracto sale con su tamaño declarado")
        self.assertNotIn(12, visible, "ni con el del otro alelo declarado")
        self.assertTrue(
            {14, 11} & set(visible),
            f"y alguno tiene que salir con una unidad menos: {sorted(set(visible))}",
        )


if __name__ == "__main__":
    unittest.main()


class NotebookReconstructionTest(unittest.TestCase):
    """**El cuaderno tiene que reconstruir el evento técnico entero.**

    El tartamudeo llegó a `ui/config_editor.py` —el lector y escritor del
    TOML— y **no** al cuaderno. La consecuencia fue doble y silenciosa: desde
    la interfaz no se podía declarar, y quien importara un fichero con
    tartamudeo y lo volviera a exportar se llevaba uno **sin** tartamudeo, sin
    aviso ninguno.

    La prueba de ida y vuelta de arriba no podía cazarlo: recorre el modelo de
    `config_editor.py` y nunca toca `ui/gps_app.py`, de modo que con una
    interfaz sin ningún control pasaba igual.
    """

    NOTEBOOK = ROOT / "ui" / "gps_app.py"

    # Los siete formularios que el cuaderno reconstruye. La invariante vale
    # para todos: **el olvido de un campo es el mismo defecto lo reconstruya
    # quien lo reconstruya**, y exigirla solo al del evento técnico dejaría los
    # otros seis donde nadie lo notaría. Hoy los siete la cumplen, así que
    # generalizarla no cambia nada más que la red.
    RECONSTRUCTED_FORMS = (
        "F5GeneralForm",
        "F5ScenarioForm",
        "GeneralForm",
        "PopulationSnvsForm",
        "ReferenceForm",
        "RepeatAdmissionForm",
        "TechnicalEventForm",
    )

    def _calls_to(self, name: str) -> list[ast.Call]:
        tree = ast.parse(self.NOTEBOOK.read_text(encoding="utf-8"))
        return [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == name
        ]

    def test_the_notebook_passes_every_field_of_every_form_it_rebuilds(self) -> None:
        """No busca `stutter`: busca **todos** los campos de **todos** los formularios.

        Comprobar solo el que se olvidó esta vez dejaría el siguiente olvido
        igual de invisible. Lo que se exige es la invariante: cada campo del
        formulario aparece como argumento con nombre en la llamada que el
        cuaderno hace. Contra el cuaderno anterior a esta corrección, `stutter`
        falta en el evento técnico y la prueba muere.

        # Y el valor no puede ser una constante

        Pasar `stutter=None` satisface «el argumento está» y **reintroduce el
        defecto original palabra por palabra**: el campo se pierde igual. Se ha
        medido: sin esta segunda mitad, esa edición pasa. Por eso se exige
        además que el valor no sea una constante literal, que es la forma que
        toma «lo paso para callar a la prueba».
        """

        module = importlib.import_module("config_editor")
        f5_module = importlib.import_module("f5_editor")
        for name in self.RECONSTRUCTED_FORMS:
            form = getattr(module, name, None) or getattr(f5_module, name)
            declared = {field.name for field in dataclasses.fields(form)}
            calls = self._calls_to(name)
            with self.subTest(form=name):
                self.assertTrue(calls, f"el cuaderno tiene que reconstruir {name}")
                for call in calls:
                    passed = {
                        keyword.arg for keyword in call.keywords if keyword.arg is not None
                    }
                    self.assertEqual(
                        declared - passed,
                        set(),
                        f"{name} pierde estos campos al reconstruirse",
                    )
                    constant = {
                        keyword.arg
                        for keyword in call.keywords
                        if keyword.arg is not None
                        and isinstance(keyword.value, ast.Constant)
                    }
                    self.assertEqual(
                        constant,
                        set(),
                        f"{name} pasa estos campos con una constante, que es "
                        "perder el campo con otro nombre",
                    )

    def test_the_notebook_offers_a_control_for_every_stutter_parameter(self) -> None:
        """Y los controles existen de verdad, en la interfaz y en su panel."""

        notebook = self.NOTEBOOK.read_text(encoding="utf-8")
        panels = (ROOT / "ui" / "ui_panels.py").read_text(encoding="utf-8")
        for control in (
            "stutter_mode_control",
            "stutter_run_length_coefficient_control",
            "stutter_motif_rates_editor",
            "stutter_steps_editor",
            "stutter_locus_rates_editor",
        ):
            with self.subTest(control=control):
                self.assertIn(control, notebook)
                self.assertIn(control, panels)
        self.assertIn("stutter_panel", notebook)
        # Control negativo con su control positivo: no hay control para un paso
        # de dos unidades, y la búsqueda sí encuentra los que deben existir.
        self.assertNotIn("units = 2", panels)
        self.assertIn("Tartamudeo de base", panels)

    def test_importing_and_exporting_keeps_the_block_through_the_notebook_tables(
        self,
    ) -> None:
        """La ida y la vuelta **por las tablas del cuaderno**, no por el TOML.

        Es el camino que perdía el bloque: formulario → tablas de la interfaz →
        formulario. Se ejercita la misma función que el cuaderno llama.
        """

        original = parse_general_toml(SCENARIO).technical_event.stutter
        rebuilt = stutter_from_editor_rows(
            original.enabled,
            original.run_length_coefficient or "",
            stutter_motif_rate_editor_rows(original),
            stutter_step_editor_rows(original),
            stutter_locus_rate_editor_rows(original),
        )
        self.assertIsNotNone(rebuilt)
        self.assertEqual(rebuilt.enabled, original.enabled)
        self.assertEqual(rebuilt.run_length_coefficient, original.run_length_coefficient)
        self.assertEqual(
            [(row.units, row.weight) for row in rebuilt.steps],
            [(row.units, row.weight) for row in original.steps],
        )
        self.assertEqual(
            [(row.locus_id, row.rate) for row in rebuilt.locus_rates],
            [(row.locus_id, row.rate) for row in original.locus_rates],
        )
        # **Solo sobrevive la tasa que de verdad se declaró.** Las otras cinco
        # coinciden con la calibración, y una cifra que repite la calibrada sin
        # tocarla no es una declaración: contarla haría que todo escenario
        # importado saliera con un bloque de tartamudeo que nadie escribió.
        declared = {row.motif_length_bases: row.rate for row in rebuilt.motif_length_rates}
        self.assertEqual(declared, {2: "0.5"})

    def test_an_empty_form_does_not_invent_a_stutter_block(self) -> None:
        """Sin nada declarado no aparece el bloque.

        # Por qué el fixture sale de las funciones y no de un literal

        La primera versión de esta prueba pasaba a mano una fila con la tasa
        **vacía**, y eso es exactamente lo que el cuaderno **nunca** produce:
        la tabla viene con las seis cifras calibradas ya escritas. El fixture
        era degenerado y la prueba no podía fallar por el defecto que decía
        cubrir — de hecho no falló cuando el defecto estaba ahí.

        Ahora las filas salen de **las mismas funciones que el cuaderno llama**
        con un escenario sin tartamudeo, que es lo que de verdad llega a la
        reconstrucción.

        # Qué rompería si esto fallara

        La tabla del tartamudeo solo existe en el documento que puede declarar
        microsatélites. Un bloque inventado en un escenario sin loci produce un
        fichero que el simulador **rechaza** como campo desconocido, y es justo
        el escenario de la congelación acotada: el más corriente de todos.
        """

        self.assertIsNone(
            stutter_from_editor_rows(
                False,
                "",
                stutter_motif_rate_editor_rows(None),
                stutter_step_editor_rows(None),
                stutter_locus_rate_editor_rows(None),
            ),
            "la tabla precargada con la calibración no declara un bloque",
        )

    def test_an_imported_document_without_stutter_round_trips_unchanged(self) -> None:
        """**El escenario más corriente: ni bloque, ni versión equivocada.**

        Un documento que declara microsatélites pero **no** tartamudeo se
        importa, se reconstruye por las tablas del cuaderno y se vuelve a
        escribir sin bloque. Y uno que ni siquiera declara microsatélites
        conserva su versión: sin esto, el cuaderno escribía la versión 2 con
        una tabla que esa versión no conoce, y el simulador lo rechazaba.
        """

        without = SCENARIO
        for block in (
            "\n[technical_event.stutter]\nenabled = true\nrun_length_coefficient = 0.0\n",
        ):
            without = without.replace(block, "\n")
        start = without.index("[[technical_event.stutter.motif_length_rates]]")
        end = without.index("[[technical_event.libraries]]")
        without = without[:start] + without[end:]
        form = parse_general_toml(without)
        self.assertIsNone(form.technical_event.stutter)
        rendered = render_general_toml(_reconstructed_through_the_notebook(form))
        self.assertNotIn("technical_event.stutter", rendered)
        self.assertIn('schema_version = "4"', rendered)

    def test_a_document_without_repeats_keeps_its_version_and_stays_runnable(self) -> None:
        """Y sin microsatélites, el documento sigue siendo el de siempre."""

        plain = """config_kind = "general"
schema_version = "2"
seed = "{seed}"

[reference]
kind = "simulated"
id = "ref-main"

[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
N = 0

[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 400

[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 1

[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
""".format(seed=SEED)
        form = parse_general_toml(plain)
        self.assertIsNone(form.technical_event.stutter)
        rendered = render_general_toml(_reconstructed_through_the_notebook(form))
        self.assertNotIn("stutter", rendered)
        self.assertIn('schema_version = "2"', rendered)

    def test_a_declared_stutter_block_forces_the_version_that_knows_it(self) -> None:
        """Y al revés: con bloque, la versión **tiene** que subir.

        Es la mitad que faltaba. Sin ella, quitar de la decisión de versión la
        mirada al bloque de tartamudeo pasaría desapercibido en cuanto el
        escenario declarase además un locus, que es lo normal.
        """

        form = parse_general_toml(SCENARIO)
        stripped = dataclasses.replace(
            form,
            reference=dataclasses.replace(
                form.reference,
                repeat_loci=[],
                repeat_admission=None,
                max_tract_output_bases=None,
            ),
        )
        self.assertIsNotNone(stripped.technical_event.stutter)
        rendered = render_general_toml(stripped)
        self.assertIn("[technical_event.stutter]", rendered)
        self.assertIn(
            'schema_version = "4"',
            rendered,
            "un documento con tartamudeo no puede salir en una versión que no lo conoce",
        )

    def test_the_prefilled_rates_match_the_rust_calibration(self) -> None:
        """Las cifras que la tabla muestra son las que Rust lleva dentro.

        Una copia decimal que se separe del original es lo que sobrevive a un
        cambio de contrato sin que nadie lo note.
        """

        rust = (ROOT / "src" / "stutter.rs").read_text(encoding="utf-8")
        # Se compara el **valor**, no el texto: Rust escribe `0.0010` y Python
        # imprime `0.001`, y una comparación de cadenas moriría por eso sin que
        # ninguna cifra hubiera cambiado.
        match = re.search(
            r"STUTTER_DEFAULT_MOTIF_LENGTH_RATES:[^=]*=\s*\[([^\]]*)\]", rust
        )
        self.assertIsNotNone(match, "la tabla calibrada tiene que estar en src/stutter.rs")
        in_rust = [float(value) for value in match.group(1).split(",")]
        in_notebook = [
            STUTTER_CALIBRATED_MOTIF_LENGTH_RATES[length] for length in (1, 2, 3, 4, 5, 6)
        ]
        self.assertEqual(in_rust, in_notebook)
