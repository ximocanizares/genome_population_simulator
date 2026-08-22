from __future__ import annotations

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import (  # noqa: E402
    EditorError,
    editor_integer,
    ordered_editor_rows,
    parse_general_toml,
    render_general_toml,
)


class ConfigEditorTests(unittest.TestCase):
    def test_minimal_and_reordered_fixtures_round_trip_semantically(self) -> None:
        fixtures = ROOT / "tests" / "fixtures" / "config"
        first = parse_general_toml((fixtures / "general-minimal.toml").read_text())
        reordered = parse_general_toml(
            (fixtures / "general-minimal-reordered.toml").read_text()
        )

        self.assertEqual(first, reordered)
        self.assertEqual(parse_general_toml(render_general_toml(first)), first)

    def test_supplied_paths_with_spaces_and_exact_percent_survive_round_trip(self) -> None:
        source = """config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "reference-main"
fasta = "inputs/reference with spaces.fa"
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
[[reduced_references]]
id = "reference-small"
source_reference_id = "reference-main"
region_size_bases = 4
reduction_percent = 33.333333333333333333
"""

        parsed = parse_general_toml(source)
        rendered = render_general_toml(parsed)
        reparsed = parse_general_toml(rendered)

        self.assertEqual(reparsed, parsed)
        self.assertIn('fasta = "inputs/reference with spaces.fa"', rendered)
        self.assertIn("33.333333333333333333", rendered)

    def test_sequence_and_reduction_order_are_normative(self) -> None:
        source = """config_kind = "general"
schema_version = "1"
[reference]
kind = "supplied"
id = "reference-main"
fasta = "source.fa"
[[reference.sequences]]
id = "chr2"
class = "chromosome"
length = 20
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 10
[[reduced_references]]
id = "second"
source_reference_id = "reference-main"
region_size_bases = 4
target_bases = 8
[[reduced_references]]
id = "first"
source_reference_id = "reference-main"
region_size_bases = 2
target_bases = 5
"""

        parsed = parse_general_toml(source)
        rendered = render_general_toml(parsed)

        self.assertLess(rendered.index('id = "chr2"'), rendered.index('id = "chr1"'))
        self.assertLess(rendered.index('id = "second"'), rendered.index('id = "first"'))

    def test_population_snv_modes_round_trip_without_future_parameters(self) -> None:
        base = """config_kind = "general"
schema_version = "1"
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
            "supplied": 'input_path = "input matrix.tsv"',
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
        for mode, path in cases.items():
            with self.subTest(mode=mode):
                source = f'{base}\n[population_snvs]\nmode = "{mode}"\n{path}\n'
                parsed = parse_general_toml(source)
                rendered = render_general_toml(parsed)
                self.assertEqual(parse_general_toml(rendered), parsed)
                self.assertIn(f'mode = "{mode}"', rendered)

    def test_population_snv_paths_must_be_non_empty(self) -> None:
        base = (ROOT / "tests" / "fixtures" / "config" / "general-minimal.toml").read_text()
        for mode, field in (("supplied", "input_path"), ("reuse", "metadata_path")):
            source = f'{base}\n[population_snvs]\nmode = "{mode}"\n{field} = ""\n'
            with self.subTest(mode=mode):
                with self.assertRaises(EditorError) as raised:
                    parse_general_toml(source)
                self.assertEqual(raised.exception.path, f"population_snvs.{field}")

    def test_future_and_branch_incompatible_fields_are_visible_errors(self) -> None:
        minimal = (ROOT / "tests" / "fixtures" / "config" / "general-minimal.toml").read_text()
        cases = [
            minimal.replace("[reference]", "unexpected = true\n[reference]", 1),
            minimal.replace(
                'kind = "simulated"', 'kind = "simulated"\nfasta = "source.fa"', 1
            ),
            minimal.replace("length = 1000", 'length = "1000"', 1),
            minimal.replace(
                "[reference]", 'seed = ""\n[reference]', 1
            ),
            minimal.replace(
                "[reference]",
                f'seed = "{"A" * 64}"\n[reference]',
                1,
            ),
        ]
        for source in cases:
            with self.subTest(source=source):
                with self.assertRaises(EditorError):
                    parse_general_toml(source)

    def test_editor_integers_never_truncate_and_row_order_is_explicit(self) -> None:
        self.assertEqual(editor_integer(2.0, "length"), 2)
        self.assertEqual(editor_integer("3", "length"), 3)
        for value in [1.5, float("nan"), "2.0", True, None]:
            with self.subTest(value=value):
                with self.assertRaises(EditorError):
                    editor_integer(value, "length")

        rows = [{"order": 2, "id": "second"}, {"order": 1, "id": "first"}]
        self.assertEqual(
            [row["id"] for row in ordered_editor_rows(rows, "rows")],
            ["first", "second"],
        )
        for invalid in [
            [{"order": 1}, {"order": 1}],
            [{"order": 0}],
            [{"order": 1.5}],
        ]:
            with self.subTest(rows=invalid):
                with self.assertRaises(EditorError):
                    ordered_editor_rows(invalid, "rows")


if __name__ == "__main__":
    unittest.main()
