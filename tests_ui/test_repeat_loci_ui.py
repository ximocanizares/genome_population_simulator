"""F6.4.9.1: declarar repeticiones en tándem desde la interfaz.

Qué hace que estas pruebas puedan fallar:

- la ida y la vuelta se comprueba **con un alelo que tiene interrupciones**, y
  además se exige que la vuelta las conserve una a una. Una interfaz que
  pensara en número de unidades las perdería y moriría aquí;
- el reparto de frecuencias inválido se rechaza **en la interfaz**, y se
  comprueba con su pareja válida al lado: sin ella, una interfaz que rechazara
  siempre pasaría la mitad de la prueba;
- la ausencia de controles se comprueba con un **control positivo**: los
  controles que sí deben existir se buscan con la misma búsqueda, de modo que
  un patrón incapaz de encontrar nada no pueda hacerse pasar por una ausencia;
- la equivalencia con la configuración efectiva se cierra ejecutando el
  simulador sobre el TOML que la interfaz escribe y leyendo lo que publica.
"""

from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import (  # noqa: E402
    EditorError,
    RepeatAlleleRow,
    RepeatInterruptionRow,
    RepeatLocusRow,
    default_form,
    normalize_repeat_frequencies,
    parse_general_toml,
    render_general_toml,
    repeat_allele_editor_rows,
    repeat_interruption_editor_rows,
    repeat_loci_from_editor_rows,
    repeat_locus_editor_rows,
)
from process_adapter import SUPPORTED_BINARY_NAME, run_toml, validate_toml  # noqa: E402

REFERENCE_BASES = "G" * 100 + "CA" * 50 + "G" * 100

LOCUS_ROWS = [
    {
        "order": 1,
        "id": "STR-1",
        "sequence": "chr1",
        "start": 101,
        "end": 200,
        "motifs": "CA",
        "max_units": "120",
        "min_purity": "0.9",
        "min_longest_run_units": "",
    }
]
ALLELE_ROWS = [
    {"order": 1, "locus_id": "STR-1", "units": "45", "sequence": "", "frequency": "0.7"},
    {"order": 2, "locus_id": "STR-1", "units": "40", "sequence": "", "frequency": "0.3"},
]
INTERRUPTION_ROWS = [
    {"order": 1, "locus_id": "STR-1", "allele_order": 2, "offset": 4, "bases": "TT"}
]


def supplied_form(loci: list[RepeatLocusRow]):
    form = default_form()
    form.reference.kind = "supplied"
    form.reference.fasta = "reference.fa"
    form.reference.composition = None
    form.reference.sequences[0].length = len(REFERENCE_BASES)
    form.reference.repeat_loci = loci
    return form


class RepeatLociUiTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_editor_rows_round_trip_through_toml_and_back(self) -> None:
        loci = repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, INTERRUPTION_ROWS)
        rendered = render_general_toml(supplied_form(loci))
        self.assertIn('schema_version = "4"', rendered)
        parsed = parse_general_toml(rendered)
        self.assertEqual(parsed.reference.repeat_loci, loci)
        self.assertEqual(repeat_locus_editor_rows(parsed.reference.repeat_loci), LOCUS_ROWS)
        self.assertEqual(
            repeat_allele_editor_rows(parsed.reference.repeat_loci), ALLELE_ROWS
        )
        self.assertEqual(
            repeat_interruption_editor_rows(parsed.reference.repeat_loci),
            INTERRUPTION_ROWS,
        )

    def test_declared_interruptions_survive_the_round_trip(self) -> None:
        loci = repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, INTERRUPTION_ROWS)
        interruptions = loci[0].alleles[1].interruptions
        self.assertEqual(
            interruptions, [RepeatInterruptionRow(offset=4, bases="TT")]
        )
        rendered = render_general_toml(supplied_form(loci))
        self.assertIn("[[reference.repeat_loci.alleles.interruptions]]", rendered)
        self.assertIn("offset = 4", rendered)
        self.assertIn('bases = "TT"', rendered)
        recovered = parse_general_toml(rendered).reference.repeat_loci[0].alleles[1]
        self.assertEqual(recovered.interruptions, interruptions)

        # Control: el mismo locus sin interrupciones escribe un documento
        # distinto, de modo que la comparación anterior no puede salir bien
        # porque las dos partes estén vacías.
        without = repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, [])
        self.assertEqual(without[0].alleles[1].interruptions, [])
        self.assertNotIn(
            "[[reference.repeat_loci.alleles.interruptions]]",
            render_general_toml(supplied_form(without)),
        )

    def test_a_sequence_allele_cannot_also_declare_interruptions(self) -> None:
        """El lector de Rust lo rechaza; la interfaz no puede ofrecerlo."""
        rows = [
            {
                "order": 1,
                "locus_id": "STR-1",
                "units": "",
                "sequence": "CA" * 45,
                "frequency": "1.0",
            }
        ]
        with self.assertRaises(EditorError):
            repeat_loci_from_editor_rows(
                LOCUS_ROWS,
                rows,
                [
                    {
                        "order": 1,
                        "locus_id": "STR-1",
                        "allele_order": 1,
                        "offset": 0,
                        "bases": "TT",
                    }
                ],
            )

    def test_the_interface_normalizes_frequencies_and_rejects_an_invalid_split(
        self,
    ) -> None:
        normalized = normalize_repeat_frequencies(
            ["0.7", "0.3"], "1e-9", "reference.repeat_loci[0].alleles"
        )
        self.assertEqual(len(normalized), 2)
        self.assertAlmostEqual(sum(normalized), 1.0, places=12)
        self.assertAlmostEqual(normalized[0], 0.7, places=12)

        invalid = [dict(row) for row in ALLELE_ROWS]
        invalid[1]["frequency"] = "0.2"
        loci = repeat_loci_from_editor_rows(LOCUS_ROWS, invalid, [])
        with self.assertRaises(EditorError) as rejected:
            render_general_toml(supplied_form(loci))
        self.assertIn("frecuencias", str(rejected.exception))

        # Y la pareja válida sí se escribe: el rechazo discrimina.
        valid = repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, [])
        self.assertIn(
            "[[reference.repeat_loci]]", render_general_toml(supplied_form(valid))
        )

    def test_an_invalid_split_can_be_imported_but_not_written_back(self) -> None:
        """La comprobación es asimétrica **a propósito**, y aquí queda fijada.

        Leer un documento con un reparto imposible tiene que funcionar, porque
        si no el usuario no podría abrir en la interfaz justo el fichero que
        necesita corregir. Escribirlo, no. Como toda ida y vuelta pasa por la
        escritura, ningún documento inválido llega a publicarse.
        """
        valid = render_general_toml(
            supplied_form(repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, []))
        )
        broken = valid.replace("frequency = 0.3", "frequency = 0.2")
        self.assertNotEqual(broken, valid)
        imported = parse_general_toml(broken)
        self.assertEqual(len(imported.reference.repeat_loci), 1)
        with self.assertRaises(EditorError):
            render_general_toml(imported)

    def test_a_negative_or_out_of_range_frequency_is_rejected(self) -> None:
        for literal in ("-0.1", "1.5", "nan"):
            with self.subTest(literal=literal):
                with self.assertRaises(EditorError):
                    normalize_repeat_frequencies(
                        [literal, "0.5"], "1e-9", "reference.repeat_loci[0].alleles"
                    )

    def test_the_declared_tract_cap_survives_the_round_trip(self) -> None:
        """El tope de bases de tracto sobrevive a abrir, tocar y guardar.

        Sin esta prueba el defecto es mudo: un escenario que declara el tope se
        abre bien —la interfaz lo lee—, y al guardarlo el tope desaparece del
        documento sin aviso.
        """
        loci = repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, INTERRUPTION_ROWS)
        form = supplied_form(loci)
        form.reference.max_tract_output_bases = 4096
        rendered = render_general_toml(form)
        self.assertIn("max_tract_output_bases = 4096", rendered)
        recovered = parse_general_toml(rendered)
        self.assertEqual(recovered.reference.max_tract_output_bases, 4096)
        # La vuelta completa: lo recuperado se vuelve a escribir igual.
        self.assertEqual(render_general_toml(recovered), rendered)

        # Control: el mismo escenario sin tope no escribe la clave, de modo que
        # la comprobación anterior no puede salir bien por estar siempre.
        without = supplied_form(loci)
        self.assertIsNone(without.reference.max_tract_output_bases)
        self.assertNotIn("max_tract_output_bases", render_general_toml(without))

    def test_the_notebook_hands_the_declared_tract_cap_to_the_form(self) -> None:
        """El cuaderno tiene el control **y se lo pasa** al formulario.

        Comparar cadenas diría solo que el nombre aparece en alguna parte. Lo
        que se comprueba aquí es la llamada: se busca en el árbol sintáctico la
        construcción de `ReferenceForm` y se exige que el tope viaje en ella.
        Un cuaderno que ofreciera el control y no lo pasara —que es justo el
        defecto que hubo— moriría aquí.
        """
        import ast

        source = (ROOT / "ui" / "gps_app.py").read_text(encoding="utf-8")
        self.assertIn("repeat_max_tract_output_control", source)
        constructions = [
            node
            for node in ast.walk(ast.parse(source))
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "ReferenceForm"
        ]
        self.assertEqual(len(constructions), 1, "el cuaderno construye un ReferenceForm")
        passed = {keyword.arg for keyword in constructions[0].keywords}
        # Control positivo: los campos que sí viajaban se encuentran con la
        # misma búsqueda, de modo que una búsqueda incapaz de hallar nada no
        # pueda hacerse pasar por un campo ausente.
        for present in ("repeat_loci", "repeat_admission", "sequences"):
            self.assertIn(present, passed)
        self.assertIn("max_tract_output_bases", passed)

    def test_the_interface_offers_no_target_purity_or_random_locus_generation(
        self,
    ) -> None:
        sources = (ROOT / "ui" / "gps_app.py").read_text(encoding="utf-8") + (
            ROOT / "ui" / "ui_panels.py"
        ).read_text(encoding="utf-8")
        # Control positivo: los controles que sí existen se encuentran con la
        # misma búsqueda con la que se niegan los otros.
        for present in (
            "repeat_loci_editor",
            "repeat_alleles_editor",
            "repeat_interruptions_editor",
            "repeat_min_purity_control",
        ):
            self.assertIn(present, sources)
        for absent in (
            "target_purity_control",
            "repeat_target_purity",
            "random_repeat_loci",
            "repeat_random_generation",
            "generate_repeat_loci",
        ):
            self.assertNotIn(absent, sources)

    def test_the_editor_rejects_a_declared_purity_or_a_random_generation_key(self) -> None:
        """La ausencia de controles se comprueba además **por comportamiento**.

        Buscar nombres en el código fuente dice que no hay un control; esto
        dice que el modelo no admite el campo aunque alguien lo escriba. La
        pureza es siempre una salida calculada y la generación aleatoria de
        loci no pertenece a este corte.
        """
        valid = render_general_toml(
            supplied_form(
                repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, INTERRUPTION_ROWS)
            )
        )
        # Control: el documento válido se lee sin protestar.
        self.assertEqual(len(parse_general_toml(valid).reference.repeat_loci), 1)
        for injected in (
            ('motif = "CA"', 'motif = "CA"\npurity = 0.95'),
            ('motif = "CA"', 'motif = "CA"\ntarget_purity = 0.95'),
            ('motif = "CA"', 'motif = "CA"\nrandom_count = 3'),
            ("units = 45", "units = 45\npurity = 0.95"),
        ):
            with self.subTest(injected=injected[1]):
                with self.assertRaises(EditorError):
                    parse_general_toml(valid.replace(*injected, 1))

    def test_a_document_without_loci_keeps_the_previous_schema_version(self) -> None:
        rendered = render_general_toml(supplied_form([]))
        self.assertIn('schema_version = "2"', rendered)
        self.assertNotIn("repeat_loci", rendered)

    def test_a_version_two_document_rejects_the_repeat_keys(self) -> None:
        rendered = render_general_toml(
            supplied_form(
                repeat_loci_from_editor_rows(LOCUS_ROWS, ALLELE_ROWS, INTERRUPTION_ROWS)
            )
        )
        downgraded = rendered.replace('schema_version = "4"', 'schema_version = "2"')
        with self.assertRaises(EditorError):
            parse_general_toml(downgraded)

    @unittest.skipUnless(
        BINARY.exists(), "compile the Rust binary before this integration test"
    )
    def test_the_authored_toml_reaches_the_effective_configuration(self) -> None:
        locus_rows = [dict(LOCUS_ROWS[0])]
        locus_rows[0]["min_purity"] = "0.8"
        loci = repeat_loci_from_editor_rows(locus_rows, ALLELE_ROWS, INTERRUPTION_ROWS)
        toml = render_general_toml(supplied_form(loci))
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            (base / "reference.fa").write_text(
                f">chr1\n{REFERENCE_BASES}\n", encoding="utf-8"
            )
            validation, session = validate_toml(self.BINARY, toml, base)
            self.assertTrue(validation.succeeded, validation.stdout or validation.stderr)
            self.assertIsNotNone(session)
            execution = run_toml(session, toml, self.BINARY, base, base / "result")
            self.assertTrue(execution.succeeded, execution.stdout or execution.stderr)
            effective = json.loads(
                (base / "result" / "effective-config.json").read_text(encoding="utf-8")
            )
        self.assertEqual(effective["schema_version"], "6")
        published = effective["repeat_loci"][0]
        self.assertEqual(published["id"], "STR-1")
        self.assertEqual(published["motifs"], ["CA"])
        self.assertEqual(published["start"], 101)
        self.assertEqual(published["end"], 200)
        self.assertEqual(published["max_units"], 120)
        self.assertEqual(published["admission"]["min_purity"], 0.8)
        # La interrupción declarada en la interfaz llega a la secuencia del
        # alelo, y de ahí a la pureza publicada, que nadie declaró.
        interrupted = next(
            allele
            for allele in published["alleles"]
            if allele["declared_frequency"] == 0.3
        )
        self.assertEqual(interrupted["sequence"][4:6], "TT")
        self.assertLess(interrupted["purity"], 1.0)
        self.assertEqual(
            interrupted["interruptions"], [{"offset": 4, "length": 2}]
        )


if __name__ == "__main__":
    unittest.main()
