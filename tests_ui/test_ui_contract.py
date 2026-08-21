from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import parse_general_toml, render_general_toml  # noqa: E402
from process_adapter import (  # noqa: E402
    SUPPORTED_BINARY_NAME,
    run_toml,
    validate_toml,
)


class UiContractTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_current_schema_options_have_controls_and_future_options_do_not(self) -> None:
        schema = json.loads(
            (ROOT / "schemas" / "config-general-v1.schema.json").read_text()
        )
        self.assertEqual(
            set(schema["properties"]),
            {
                "config_kind",
                "schema_version",
                "seed",
                "reference",
                "reduced_references",
            },
        )

        source = (ROOT / "ui" / "gps_app.py").read_text()
        for control in (
            "seed_control",
            "reference_kind_control",
            "reference_id_control",
            "fasta_control",
            "composition_controls",
            "sequence_editor",
            "reduction_editor",
        ):
            with self.subTest(control=control):
                self.assertIn(control, source)
        for future_option in ("scenario_id", "ploidy_control", "population_control"):
            with self.subTest(future_option=future_option):
                self.assertNotIn(future_option, source)

    def test_nine_tabs_labels_keyboard_and_textual_states_are_present(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text()
        tabs = (
            "General/Ejecución",
            "Genoma",
            "Individuos",
            "Poblaciones",
            "Evolución",
            "Familias",
            "Eventos técnicos",
            "Salidas",
            "Revisar y ejecutar",
        )
        positions = [source.index(f'"{tab}"') for tab in tabs]
        self.assertEqual(positions, sorted(positions))
        self.assertIn('keyboard_shortcut="Ctrl-Shift-V"', source)
        self.assertIn('keyboard_shortcut="Ctrl-Shift-R"', source)
        for label in (
            "Importar configuración TOML general",
            "TOML general previsualizado",
            "Diagnósticos y logs del proceso",
            "Diagnósticos localizados por pestaña y campo",
        ):
            self.assertIn(f'label="{label}"', source)
        for textual_state in (
            "Pendiente de validación con Rust.",
            "Configuración validada por el binario Rust.",
            "Ejecución no iniciada.",
            "Ejecución completada.",
        ):
            self.assertIn(textual_state, source)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_ui_and_manual_toml_produce_the_same_effective_configuration(self) -> None:
        manual = """schema_version = "1"
config_kind = "general"
seed = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"

[reference]
id = "reference-main"
kind = "simulated"

[reference.composition]
T = 0.25
G = 0.25
C = 0.25
A = 0.25

[[reference.sequences]]
length = 1000
class = "chromosome"
id = "chr1"
"""
        ui_toml = render_general_toml(parse_general_toml(manual))

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            results = []
            for name, document in (("manual", manual), ("ui", ui_toml)):
                validation, session = validate_toml(self.BINARY, document, base)
                self.assertTrue(validation.succeeded, validation.stdout)
                self.assertIsNotNone(session)
                execution = run_toml(
                    session,
                    document,
                    self.BINARY,
                    base,
                    base / name,
                )
                self.assertTrue(execution.succeeded, execution.stdout)
                effective = json.loads(
                    (base / name / "effective-config.json").read_text()
                )
                effective.pop("source_path")
                results.append(effective)

        self.assertEqual(results[0], results[1])


if __name__ == "__main__":
    unittest.main()
