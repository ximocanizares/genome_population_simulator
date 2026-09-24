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

import gzip
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
    EditorError,
    parse_general_toml,
    render_general_toml,
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

        Partiendo solo de la interfaz se construye un escenario completo
        —referencia, locus, alelos, frecuencias y tasa— y se **ejecuta**: el
        documento que la interfaz escribe pasa por el lector de Rust y produce
        lecturas tartamudeadas de verdad.

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
