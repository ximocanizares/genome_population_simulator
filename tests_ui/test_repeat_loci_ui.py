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

from ui_panels import (  # noqa: E402
    REPEAT_ALLELE_SOURCES,
    REPEAT_ALLELE_SOURCE_DECLARED,
    REPEAT_ALLELE_SOURCE_GENERATED,
    repeats_panel,
)
from config_editor import (  # noqa: E402
    EditorError,
    InitialPopulationRow,
    PopulationSnvsForm,
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
    repeat_slippage_editor_rows,
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

GENERATED_LOCUS_ROWS = [
    {
        "order": 1,
        "id": "STR-GEN",
        "sequence": "chr1",
        "start": 101,
        "end": 200,
        "motifs": "CA",
        "max_units": "120",
        "min_purity": "",
        "min_longest_run_units": "",
    }
]

SLIPPAGE_ROWS = [
    {
        "order": 1,
        "locus_id": "STR-GEN",
        "model": "stepwise",
        "mean_events": "3.0",
        "expansion_bias": "0.4",
        "step_units": "1",
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


class RecordingMo:
    """Doble de la biblioteca de la interfaz que registra qué hijos recibe."""

    def __init__(self) -> None:
        self.vstack_children: list = []

    @staticmethod
    def md(value: str):
        return ("markdown", value)

    @staticmethod
    def callout(value: str, kind: str = "info"):
        return ("callout", kind, value)

    def vstack(self, children):
        self.vstack_children = list(children)
        return ("vstack", children)


class _SourceControl:
    """Doble del selector: solo tiene el valor que el panel mira."""

    def __init__(self, value: str) -> None:
        self.value = value


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

    # -----------------------------------------------------------------------
    # F6.4.13 — el generador de alelos por deslizamiento en la interfaz
    # -----------------------------------------------------------------------

    def test_the_generator_survives_the_round_trip_through_toml_and_back(self) -> None:
        """La ida y la vuelta de un locus **generado**, con sus cuatro campos.

        El fixture declara los cuatro —modelo, media, sesgo y tamaño de paso—
        porque los dos últimos son opcionales: si solo se escribieran los
        obligatorios, una interfaz que perdiera los opcionales pasaría igual.
        """

        loci = repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], SLIPPAGE_ROWS)
        self.assertEqual(loci[0].alleles, [])
        self.assertIsNotNone(loci[0].slippage)
        rendered = render_general_toml(supplied_form(loci))
        self.assertIn("[reference.repeat_loci.slippage]", rendered)
        self.assertNotIn("[[reference.repeat_loci.alleles]]", rendered)
        parsed = parse_general_toml(rendered)
        self.assertEqual(parsed.reference.repeat_loci, loci)
        self.assertEqual(
            repeat_locus_editor_rows(parsed.reference.repeat_loci),
            GENERATED_LOCUS_ROWS,
        )

    def test_the_panel_shows_the_generator_or_the_allele_tables_but_never_both(
        self,
    ) -> None:
        """**El selector oculta de verdad la tabla que no toca.**

        Es lo que el paso F6.4.13 pide con esas palabras, y lo que antes no
        existía: había cuatro columnas más en la tabla de loci, ningún selector
        y nada que se ocultara. La prueba que llevaba este nombre comprobaba
        otra cosa —los rechazos del lector— y el nombre afirmaba lo que no
        verificaba.

        Se ejercita el panel **de verdad**, con un doble de la biblioteca de la
        interfaz que registra qué hijos recibe, y se mira en los dos sentidos.
        """

        panels = RecordingMo()
        declared = repeats_panel(
            panels,
            "LOCI",
            _SourceControl(REPEAT_ALLELE_SOURCE_DECLARED),
            "ALELOS",
            "INTERRUPCIONES",
            "GENERADOR",
            "ADMISION",
        )
        self.assertEqual(declared[0], "vstack")
        self.assertIn("ALELOS", panels.vstack_children)
        self.assertIn("INTERRUPCIONES", panels.vstack_children)
        self.assertNotIn("GENERADOR", panels.vstack_children)

        panels = RecordingMo()
        repeats_panel(
            panels,
            "LOCI",
            _SourceControl(REPEAT_ALLELE_SOURCE_GENERATED),
            "ALELOS",
            "INTERRUPCIONES",
            "GENERADOR",
            "ADMISION",
        )
        self.assertIn("GENERADOR", panels.vstack_children)
        self.assertNotIn("ALELOS", panels.vstack_children)
        self.assertNotIn("INTERRUPCIONES", panels.vstack_children)

        # En los dos sentidos se ven el locus, el propio selector y los mínimos:
        # el selector oculta las tablas de alelos, no el panel entero.
        for value in REPEAT_ALLELE_SOURCES:
            panels = RecordingMo()
            control = _SourceControl(value)
            repeats_panel(
                panels, "LOCI", control, "ALELOS", "INTERRUPCIONES", "GENERADOR", "ADMISION"
            )
            self.assertIn("LOCI", panels.vstack_children)
            self.assertIn(control, panels.vstack_children)
            self.assertIn("ADMISION", panels.vstack_children)

    def test_declaring_the_alleles_and_the_generator_at_once_is_rejected(self) -> None:
        """**Son excluyentes en los dos sentidos, y con control positivo.**

        Declarar las dos cosas se rechaza, y no declarar ninguna también: un
        locus sin alelos y sin generador no dice de dónde salen sus alelos. Las
        dos formas válidas —la declarada y la generada— se construyen al lado,
        de modo que una interfaz que rechazara siempre no pasaría.
        """

        both_alleles = [
            {
                "order": 1,
                "locus_id": "STR-GEN",
                "units": "50",
                "sequence": "",
                "frequency": "1.0",
            }
        ]
        with self.assertRaises(EditorError) as declared_and_generated:
            repeat_loci_from_editor_rows(
                GENERATED_LOCUS_ROWS, both_alleles, [], SLIPPAGE_ROWS
            )
        self.assertIn("generen", str(declared_and_generated.exception))

        with self.assertRaises(EditorError) as without_either:
            repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], [])
        self.assertIn("ningún alelo", str(without_either.exception))

        # **Controles positivos**: cada una por separado se acepta. Sin ellos,
        # una interfaz que rechazara cualquier locus pasaría las dos aserciones
        # anteriores sin distinguir nada.
        only_generated = repeat_loci_from_editor_rows(
            GENERATED_LOCUS_ROWS, [], [], SLIPPAGE_ROWS
        )
        self.assertEqual(only_generated[0].alleles, [])
        self.assertIsNotNone(only_generated[0].slippage)
        only_declared = repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, both_alleles, [], [])
        self.assertEqual(len(only_declared[0].alleles), 1)
        self.assertIsNone(only_declared[0].slippage)

    def test_a_generator_without_its_model_or_its_mean_is_rejected(self) -> None:
        """Media sin modelo y modelo sin media se denuncian, no se completan.

        Y **dos generadores para el mismo locus** también: quedarse con el
        primero perdería en silencio lo que el usuario acaba de escribir.
        """

        without_model = [dict(SLIPPAGE_ROWS[0])]
        without_model[0]["model"] = ""
        with self.assertRaises(EditorError) as missing_model:
            repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], without_model)
        self.assertIn("modelo", str(missing_model.exception))

        without_mean = [dict(SLIPPAGE_ROWS[0])]
        without_mean[0]["mean_events"] = ""
        with self.assertRaises(EditorError) as missing_mean:
            repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], without_mean)
        self.assertIn("media de pasos", str(missing_mean.exception))

        duplicated = [dict(SLIPPAGE_ROWS[0]), dict(SLIPPAGE_ROWS[0])]
        with self.assertRaises(EditorError) as twice:
            repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], duplicated)
        self.assertIn("más de un generador", str(twice.exception))

        # Una fila del generador que nombra un locus inexistente se denuncia,
        # igual que una de alelos. El locus real lleva aquí sus alelos, para que
        # lo que falle sea la fila descolgada y no la falta de alelos.
        unknown = [dict(SLIPPAGE_ROWS[0])]
        unknown[0]["locus_id"] = "STR-INEXISTENTE"
        with self.assertRaises(EditorError) as stray:
            repeat_loci_from_editor_rows(
                GENERATED_LOCUS_ROWS,
                [
                    {
                        "order": 1,
                        "locus_id": "STR-GEN",
                        "units": "50",
                        "sequence": "",
                        "frequency": "1.0",
                    }
                ],
                [],
                unknown,
            )
        self.assertIn("no existe", str(stray.exception))

    def test_the_authored_generator_reaches_the_effective_configuration(self) -> None:
        """El TOML que **escribe la interfaz** produce un abanico generado.

        Se ejecuta el binario de verdad y se leen los bytes publicados: es la
        equivalencia interfaz → TOML → configuración efectiva que la regla de
        cierre de la hoja de ruta exige.
        """

        loci = repeat_loci_from_editor_rows(GENERATED_LOCUS_ROWS, [], [], SLIPPAGE_ROWS)
        form = supplied_form(loci)
        # **El generador simula un linaje por copia haplotípica**, de modo que
        # sin población inicial no habría ni un linaje. Cuatro individuos
        # diploides dan ocho, que es lo que la prueba comprueba después.
        form.population_snvs = PopulationSnvsForm(
            mode="generate",
            input_path=None,
            metadata_path=None,
            generate={
                "count_mode": "exact",
                "locus_count": 1,
                # **La variante se confina al flanco.** Una posición sorteada
                # dentro del tracto declarado rechaza el escenario por la
                # decisión S-32, y con semilla de entropía eso ocurría una de
                # cada tres ejecuciones: la prueba fallaba por el fixture y no
                # por el producto.
                "spatial_mode": "weighted_regions",
                "outside_weight": 0.0,
                "regions": [
                    {"contig": "chr1", "start": 1, "end": 100, "weight": 1.0}
                ],
                "substitution_profile": "uniform_snv12_v1",
                "populations": [
                    {
                        "population_id": "P1",
                        "frequency_model": "uniform",
                        "min_alt_af": 0.2,
                        "max_alt_af": 0.8,
                    }
                ],
            },
        )
        form.initial_populations = [
            InitialPopulationRow(
                population_id="P1",
                individual_count=4,
                ploidy=2,
                genotype_model="hwe",
                inbreeding_coefficient=None,
                individual_ids=["I1", "I2", "I3", "I4"],
                founder_keys=None,
            )
        ]
        toml = render_general_toml(form)
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
        published = effective["repeat_loci"][0]
        self.assertEqual(published["id"], "STR-GEN")
        self.assertEqual(published["slippage"]["model"], "stepwise")
        self.assertEqual(published["slippage"]["mean_events"], 3.0)
        self.assertEqual(published["slippage"]["expansion_bias"], 0.4)
        self.assertEqual(published["slippage"]["step_units"], 1)
        self.assertEqual(published["slippage"]["lineages"], 8)
        self.assertTrue(published["alleles"])
        for allele in published["alleles"]:
            self.assertEqual(allele["origin"], "generated")

    def test_every_option_this_work_added_to_the_schema_has_a_control(self) -> None:
        """**Recorrido completo de opciones configurables, y regla de cierre.**

        Recorre las claves declarables que este trabajo —las cinco entregas de
        repeticiones y de preparación amplificada— añadió al documento general,
        una a una, y exige que **cada una** tenga un control con nombre en la
        interfaz. Es la única comprobación que caza de forma sistemática una
        opción añadida sin control, y es la que permite dar por cerrada la regla
        de cierre de la hoja de ruta.

        # Por qué la tabla de correspondencia es el corazón de la prueba

        La lista de opciones sale **del esquema**, no de una lista escrita a
        mano, y la tabla dice qué control atiende a cada una. Se exige que las
        dos coincidan **exactamente**: una opción nueva en el esquema sin
        entrada en la tabla rompe la prueba, y una entrada de la tabla que ya no
        corresponda a ninguna opción también. Así no se puede añadir una opción
        y olvidarla, que es justo el fallo que esto existe para impedir.

        **Control positivo**: la lista recorrida no está vacía y contiene las
        claves que se saben presentes. Sin él, un barrido que no encontrara
        ninguna opción pasaría sin mirar nada.

        # Lo que este barrido NO comprueba, dicho sin adornos

        Comprueba que el control **exista con ese nombre y esté definido** en el
        cuaderno o en su panel, no que **atienda** a esa opción. Una entrada de
        la tabla que apuntara a un control ajeno que ya existe pasaría: se midió
        y pasa. Es un suelo contra el olvido —añadir una opción sin darle control
        rompe la prueba, porque la tabla tiene que cubrir el esquema entero—, no
        una demostración de que el control hace lo que su fila dice. Atarlo de
        verdad exigiría ejecutar el cuaderno y observar el formulario que
        produce, que es trabajo de otra clase. Anotado como DT-012.
        """

        schema = json.loads(
            (ROOT / "schemas" / "config-general-v4.schema.json").read_text(
                encoding="utf-8"
            )
        )
        definitions = schema["$defs"]
        options: set[str] = set()
        for definition in (
            "repeat_locus",
            "repeat_allele",
            "repeat_interruption",
            "repeat_admission",
            "repeat_slippage",
        ):
            options.update(definitions[definition]["properties"])
        options.add("max_tract_output_bases")
        technical = definitions["technical_event"]["properties"]
        for block in ("stutter", "pcr"):
            node = technical[block]
            options.update(node["properties"])
            for child in node["properties"].values():
                if child.get("type") == "array" and "properties" in child.get(
                    "items", {}
                ):
                    options.update(child["items"]["properties"])

        # Qué control atiende a cada opción. Para una celda de tabla es el
        # nombre literal de su columna, entre comillas; para un bloque entero,
        # el nombre del editor o del control que lo edita.
        controls = {
            # Locus de repetición
            "id": '"id"',
            "sequence": '"sequence"',
            "start": '"start"',
            "end": '"end"',
            "motif": '"motifs"',
            "motifs": '"motifs"',
            "max_units": '"max_units"',
            "min_purity": '"min_purity"',
            "min_longest_run_units": '"min_longest_run_units"',
            "alleles": "repeat_alleles_editor",
            "interruptions": "repeat_interruptions_editor",
            "max_tract_output_bases": "max_tract_output_bases",
            # Alelos e interrupciones
            "units": '"units"',
            "frequency": '"frequency"',
            "offset": '"offset"',
            "bases": '"bases"',
            # Generador por deslizamiento. La elección entre declarar y generar
            # la atiende el **selector**, que es lo que el paso F6.4.13 pide.
            "slippage": "repeat_allele_source_control",
            "model": '"model"',
            "mean_events": '"mean_events"',
            "expansion_bias": '"expansion_bias"',
            "step_units": '"step_units"',
            # Tartamudeo
            "enabled": "stutter_mode_control",
            "run_length_coefficient": "stutter_run_length_coefficient_control",
            "motif_length_rates": "stutter_motif_rates_editor",
            "motif_length_bases": "stutter_motif_rates_editor",
            "rate": "stutter_motif_rates_editor",
            "steps": "stutter_steps_editor",
            "weight": "stutter_steps_editor",
            "locus_rates": "stutter_locus_rates_editor",
            "locus_id": "stutter_locus_rates_editor",
            # Preparación amplificada
            "duplicate_fraction": "pcr_duplicate_fraction_control",
            "cycles": "pcr_cycles_control",
            "efficiency": "pcr_efficiency_control",
            "early_cycles": "pcr_early_cycles_control",
            "max_family_size": "pcr_max_family_size_control",
            "gc_max_attempts": "pcr_gc_max_attempts_control",
            "gc_relative_efficiency": "pcr_gc_curve_editor",
            "gc_fraction": "pcr_gc_curve_editor",
            "relative_efficiency": "pcr_gc_curve_editor",
            "efficiency_classes": "pcr_efficiency_classes_editor",
            "stutter": "stutter_mode_control",
        }

        # Control positivo: el barrido encontró algo, y encontró lo que se sabe
        # que está.
        self.assertGreaterEqual(len(options), 30)
        for known in ("motif", "max_units", "min_purity", "model", "mean_events"):
            self.assertIn(known, options)
        # Y la tabla cubre **exactamente** las opciones del esquema.
        self.assertEqual(
            set(controls),
            options,
            "la tabla de controles y las opciones del esquema han dejado de"
            " coincidir: hay una opción sin control declarado o un control sin"
            " opción",
        )

        sources = "\n".join(
            (ROOT / "ui" / name).read_text(encoding="utf-8")
            for name in ("gps_app.py", "ui_panels.py")
        )
        for option in sorted(options):
            control = controls[option]
            with self.subTest(option=option):
                self.assertIn(
                    control,
                    sources,
                    f"la opción {option} no tiene control en la interfaz",
                )
                # Un control con nombre propio tiene que estar **definido**, no
                # solo nombrado. Sin esto, renombrar la asignación y dejar las
                # menciones colgando pasaba el barrido con el cuaderno roto: se
                # midió y pasaba. Una columna de tabla se busca entrecomillada y
                # no tiene asignación que comprobar.
                if not control.startswith('"'):
                    self.assertIn(
                        f"{control} = ",
                        sources,
                        f"el control {control} se nombra pero no se define",
                    )

if __name__ == "__main__":
    unittest.main()
