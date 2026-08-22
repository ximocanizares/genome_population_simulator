import sys
import unittest
import json
import importlib.util
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ui"))

from config_editor import EditorError, parse_general_toml, render_general_toml


class PopulationSnvsEditorTest(unittest.TestCase):
    def test_draft_schema_rejects_generate_with_supplied_or_reuse(self):
        if importlib.util.find_spec("jsonschema") is None:
            self.skipTest("jsonschema no está instalado en el entorno UI")
        from jsonschema import Draft202012Validator

        schema = json.loads((Path(__file__).resolve().parents[1] / "schemas/config-general-v1.schema.json").read_text())
        validator = Draft202012Validator(schema)
        generation = {
            "count_mode": "exact",
            "locus_count": 1,
            "spatial_mode": "uniform",
            "substitution_profile": "uniform_snv12_v1",
            "populations": [{
                "population_id": "P1",
                "frequency_model": "uniform",
                "min_alt_af": 0.1,
                "max_alt_af": 0.9,
            }],
        }
        reference = {
            "kind": "simulated",
            "id": "ref",
            "composition": {base: 0.25 for base in "ACGT"},
            "sequences": [{"id": "chr1", "class": "chromosome", "length": 5}],
        }
        for mode, path in (("supplied", "input_path"), ("reuse", "metadata_path")):
            document = {
                "config_kind": "general",
                "schema_version": "1",
                "reference": reference,
                "population_snvs": {"mode": mode, path: "matrix", "generate": generation},
            }
            self.assertTrue(list(validator.iter_errors(document)))

        inverted_af = {
            "config_kind": "general",
            "schema_version": "1",
            "reference": reference,
            "population_snvs": {
                "mode": "generate",
                "generate": {**generation, "populations": [{**generation["populations"][0], "min_alt_af": 0.9, "max_alt_af": 0.1}]},
            },
        }
        self.assertFalse(list(validator.iter_errors(inverted_af)))

        inverted_region = {
            "config_kind": "general",
            "schema_version": "1",
            "reference": reference,
            "population_snvs": {
                "mode": "generate",
                "generate": {
                    **generation,
                    "spatial_mode": "weighted_regions",
                    "outside_weight": 1.0,
                    "regions": [{"contig": "chr1", "start": 4, "end": 3, "weight": 1.0}],
                },
            },
        }
        self.assertFalse(list(validator.iter_errors(inverted_region)))
    def test_generate_round_trip_preserves_poisson_regions_and_populations(self):
        source = '''config_kind = "general"
schema_version = "1"
[population_snvs]
mode = "generate"
[population_snvs.generate]
count_mode = "poisson"
mean_density_per_base = 0.2
spatial_mode = "weighted_regions"
substitution_profile = "uniform_snv12_v1"
outside_weight = 0.5
[[population_snvs.generate.regions]]
contig = "chr1"
start = 1
end = 2
weight = 3.0
[[population_snvs.generate.populations]]
population_id = "P1"
frequency_model = "uniform"
min_alt_af = 0.1
max_alt_af = 0.7
[[population_snvs.generate.populations]]
population_id = "P2"
frequency_model = "uniform"
min_alt_af = 0.2
max_alt_af = 0.8
[reference]
kind = "simulated"
id = "ref"
[reference.composition]
A = 0.25
C = 0.25
G = 0.25
T = 0.25
[[reference.sequences]]
id = "chr1"
class = "chromosome"
length = 2
'''
        form = parse_general_toml(source)
        rendered = render_general_toml(form)
        reparsed = parse_general_toml(rendered)
        self.assertEqual(form.population_snvs, reparsed.population_snvs)

    def test_generate_requires_its_nested_table(self):
        with self.assertRaises(EditorError):
            parse_general_toml('''config_kind = "general"\nschema_version = "1"\n[population_snvs]\nmode = "generate"\n[reference]\nkind = "simulated"\nid = "ref"\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n[[reference.sequences]]\nid = "chr1"\nclass = "chromosome"\nlength = 1\n''')

    def test_supplied_and_reuse_reject_a_generate_table(self):
        base = '''config_kind = "general"\nschema_version = "1"\n[population_snvs]\nmode = "{mode}"\n{path} = "matrix"\n[population_snvs.generate]\ncount_mode = "exact"\nlocus_count = 1\nspatial_mode = "uniform"\nsubstitution_profile = "uniform_snv12_v1"\n[[population_snvs.generate.populations]]\npopulation_id = "P1"\nfrequency_model = "uniform"\nmin_alt_af = 0.1\nmax_alt_af = 0.9\n[reference]\nkind = "simulated"\nid = "ref"\n[reference.composition]\nA = 0.25\nC = 0.25\nG = 0.25\nT = 0.25\n[[reference.sequences]]\nid = "chr1"\nclass = "chromosome"\nlength = 1\n'''
        for mode, path in (("supplied", "input_path"), ("reuse", "metadata_path")):
            with self.assertRaises(EditorError) as raised:
                parse_general_toml(base.format(mode=mode, path=path))
            self.assertIn("generate", raised.exception.detail)


if __name__ == "__main__":
    unittest.main()
