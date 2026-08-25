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
    vcf_artifact_summary,
)


class F25VcfUiRouteTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_ui_toml_routes_configuration_to_truth_and_optional_vcf(self) -> None:
        ui_toml = render_general_toml(replace(default_form(), vcf_output_mode="bgzf"))
        self.assertEqual(parse_general_toml(ui_toml).vcf_output_mode, "bgzf")

        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            validation, session = validate_toml(self.BINARY, ui_toml, base)
            self.assertTrue(validation.succeeded, validation.stdout or validation.stderr)
            self.assertIsNotNone(session)
            execution = run_toml(session, ui_toml, self.BINARY, base, base / "result")
            self.assertTrue(execution.succeeded, execution.stdout or execution.stderr)

            result = base / "result"
            self.assertTrue((result / "outputs" / "truth.jsonl").is_file())
            self.assertTrue((result / "outputs" / "truth.vcf.gz").is_file())
            manifest = json.loads((result / "manifest.json").read_text())
            self.assertEqual(manifest["vcf_outputs"][0]["mode"], "bgzf")
            self.assertIn("truth.vcf.gz", vcf_artifact_summary(result))


if __name__ == "__main__":
    unittest.main()
