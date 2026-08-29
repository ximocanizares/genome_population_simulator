from __future__ import annotations

import json
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import default_form, parse_general_toml, render_general_toml  # noqa: E402
from process_adapter import (  # noqa: E402
    SUPPORTED_BINARY_NAME,
    run_toml,
    validate_toml,
)
from ui_panels import individuals_panel  # noqa: E402


class RecordingMo:
    def __init__(self) -> None:
        self.vstack_children = None

    @staticmethod
    def md(value: str) -> tuple[str, str]:
        return ("markdown", value)

    def vstack(self, children):
        self.vstack_children = children
        return ("vstack", children)


class UiContractTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_current_schema_options_have_controls_and_future_options_do_not(self) -> None:
        schema = json.loads(
            (ROOT / "schemas" / "config-general-v2.schema.json").read_text()
        )
        self.assertEqual(
            set(schema["properties"]),
            {
                "config_kind",
                "direct_events",
                "initial_populations",
                "normalization_tolerance",
                "ploidy",
                "population_snvs",
                "schema_version",
                "seed",
                "technical_event",
                "reference",
                "reduced_references",
                "vcf_output",
            },
        )

        source = (ROOT / "ui" / "gps_app.py").read_text()
        for control in (
            "seed_control",
            "load_default_config_button",
            "ploidy_control",
            "normalization_tolerance_control",
            "reference_kind_control",
            "reference_id_control",
            "fasta_control",
            "fasta_browser",
            "composition_controls",
            "sequence_editor",
            "reduction_editor",
            "population_snvs_mode_control",
            "population_snvs_input_path_control",
            "population_snvs_metadata_path_control",
            "initial_populations_editor",
            "direct_events_editor",
            "technical_event_mode_control",
            "technical_event_coverage_control",
            "technical_libraries_editor",
            "allow_incompatible_genotypes_control",
            "base_directory_browser",
            "output_parent_directory_browser",
            "output_name_control",
            "vcf_output_mode_control",
        ):
            with self.subTest(control=control):
                self.assertIn(control, source)
        for future_option in (
            "scenario_id",
            "population_control",
            "haplotype_construction",
            "recombination_probability_per_interbase",
        ):
            with self.subTest(future_option=future_option):
                self.assertNotIn(future_option, source)

    def test_reference_composition_notice_is_shown_once_and_path_browsers_stay_local(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text()

        self.assertEqual(source.count("A/C/G/T deben sumar 1"), 1)
        self.assertIn('selection_mode="file"', source)
        self.assertIn('selection_mode="directory"', source)
        self.assertIn("puede navegar fuera del proyecto", source)

    def test_default_project_toml_requires_an_explicit_load_action(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text()

        self.assertIn("Cargar configuración predeterminada del proyecto", source)
        self.assertIn("tests/fixtures/config/general-minimal.toml", source)
        self.assertIn("elif load_default_config_button.value", source)

    def test_supplied_reference_reduction_uses_readable_objectives_and_explains_scope(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text()
        editor_source = (ROOT / "ui" / "config_editor.py").read_text()

        for label in (
            "Cómo crear una referencia reducida",
            "Objetivo",
            "Valor",
            "Eliminar porcentaje",
            "Conservar bases",
            "global entre cromosomas y scaffolds elegibles",
            "bloques no solapados de `region_size_bases`",
            "valor 10 elimina el 10 % y retiene el 90 %",
            "Para conservar el 10 %, indique «Eliminar porcentaje» y valor 90",
            "solo está disponible para una referencia FASTA suministrada",
        ):
            with self.subTest(label=label):
                self.assertIn(label, source + editor_source)
        self.assertNotIn('row["target_kind"]', source)
        self.assertNotIn('row["target_literal"]', source)

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
            "Resumen de fundadores F2.3.5 (sin cargar secuencias)",
            "Resumen del VCF del manifiesto (sin cargar el VCF)",
        ):
            self.assertIn(f'label="{label}"', source)
        for textual_state in (
            "Pendiente de validación con Rust.",
            "Configuración validada por el binario Rust.",
            "Ejecución no iniciada.",
            "Ejecución completada.",
        ):
            self.assertIn(textual_state, source)

    def test_individual_controls_are_limited_to_f2_3_binary_fields(self) -> None:
        source = (ROOT / "ui" / "gps_app.py").read_text()
        for field in (
            "individual_count",
            "genotype_model",
            "inbreeding_coefficient",
            "individual_ids",
            "founder_keys",
            "position_1based",
            "--allow-incompatible-genotypes",
        ):
            with self.subTest(field=field):
                self.assertIn(field, source if field != "--allow-incompatible-genotypes" else (ROOT / "ui" / "process_adapter.py").read_text())

    def test_individuals_panel_contains_the_real_ploidy_control(self) -> None:
        mo = RecordingMo()
        control = object()

        panel = individuals_panel(mo, control)

        self.assertEqual(panel[0], "vstack")
        self.assertEqual(mo.vstack_children[0], ("markdown", "## Individuos"))
        self.assertIs(mo.vstack_children[1], control)

    def test_ploidy_boundary_values_are_exported_to_toml(self) -> None:
        for ploidy in (1, 10):
            with self.subTest(ploidy=ploidy):
                rendered = render_general_toml(replace(default_form(), ploidy=ploidy))
                self.assertIn(f"ploidy = {ploidy}", rendered)
                self.assertEqual(parse_general_toml(rendered).ploidy, ploidy)

    def test_normalization_tolerance_preserves_its_decimal_literal_on_import_export(self) -> None:
        literal = "0.000000001"
        rendered = render_general_toml(
            replace(default_form(), normalization_tolerance=literal)
        )
        self.assertIn(f"normalization_tolerance = {literal}", rendered)
        self.assertEqual(parse_general_toml(rendered).normalization_tolerance, literal)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_ui_technical_event_v2_round_trips_to_rust_validation(self) -> None:
        manual = """config_kind = "general"
schema_version = "2"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 1000
[[initial_populations]]
population_id = "P1"
individual_count = 1
ploidy = 2
genotype_model = "hwe"
individual_ids = ["sample-A"]
[technical_event]
event_id = "reads-main"
individual_ids = ["sample-A"]
target_coverage_wgs = 0
read_length_bases = 150
[[technical_event.libraries]]
individual_id = "sample-A"
i7 = "ACGTACGT"
i5 = "TGCATGCA"
"""
        exported = render_general_toml(parse_general_toml(manual))
        self.assertEqual(parse_general_toml(exported), parse_general_toml(manual))

        with tempfile.TemporaryDirectory() as directory:
            validation, session = validate_toml(
                self.BINARY, exported, Path(directory)
            )
        self.assertFalse(validation.succeeded)
        self.assertIsNone(session)
        self.assertIn("F4.1 rechaza la ejecución", validation.stdout)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_population_snv_import_export_and_rust_validation_cover_every_mode(self) -> None:
        base = """config_kind = "general"
schema_version = "2"
[reference]
kind = "simulated"
id = "reference-main"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
"""
        cases = {
            "supplied": 'input_path = "input.tsv"',
            "reuse": 'metadata_path = "population_snvs.metadata.json"',
            "generate": '''[population_snvs.generate]
count_mode = "exact"
locus_count = 1
spatial_mode = "uniform"
substitution_profile = "uniform_snv12_v1"
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.9''',
        }
        with tempfile.TemporaryDirectory() as directory:
            base_path = Path(directory)
            for mode, path in cases.items():
                with self.subTest(mode=mode):
                    source = f'{base}\n[population_snvs]\nmode = "{mode}"\n{path}\n'
                    exported = render_general_toml(parse_general_toml(source))
                    self.assertEqual(parse_general_toml(exported), parse_general_toml(source))
                    result, session = validate_toml(self.BINARY, exported, base_path)
                    if mode == "generate":
                        self.assertTrue(result.succeeded, result.stdout or result.stderr)
                        self.assertIsNotNone(session)
                    else:
                        self.assertFalse(result.succeeded)
                        self.assertIsNone(session)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_ui_and_manual_toml_produce_the_same_effective_configuration(self) -> None:
        manual = """schema_version = "2"
config_kind = "general"
seed = "0123456789abcdef0123456789abcdef0123456789abcdef0123456789abcdef"
normalization_tolerance = 0.000000001

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
                self.assertTrue((base / name / "outputs" / "truth.jsonl").is_file())
                effective = json.loads(
                    (base / name / "effective-config.json").read_text()
                )
                effective.pop("source_path")
                results.append(effective)

        self.assertEqual(results[0], results[1])


if __name__ == "__main__":
    unittest.main()
