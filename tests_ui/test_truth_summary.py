from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from process_adapter import truth_artifact_summary  # noqa: E402


class TruthArtifactSummaryTests(unittest.TestCase):
    def test_reads_only_the_small_manifest_entry(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            truth = output / "outputs" / "truth.jsonl"
            truth.parent.mkdir()
            truth.write_bytes(b"x" * 123)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "artifacts": [
                            {
                                "id": "truth-jsonl",
                                "kind": "truth_jsonl",
                                "path": "outputs/truth.jsonl",
                                "size_bytes": 123,
                                "checksum": {"algorithm": "sha256", "digest": "a" * 64},
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )

            summary = truth_artifact_summary(output)

        self.assertIn("ruta=outputs/truth.jsonl", summary)
        self.assertIn("tamaño=123 bytes", summary)
        self.assertIn("SHA-256=" + "a" * 64, summary)

    def test_rejects_a_missing_canonical_truth_file(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "artifacts": [
                            {
                                "id": "truth-jsonl",
                                "kind": "truth_jsonl",
                                "path": "outputs/truth.jsonl",
                                "size_bytes": 123,
                                "checksum": {"algorithm": "sha256", "digest": "a" * 64},
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )

            summary = truth_artifact_summary(output)

        self.assertIn("incompleto o es inválido", summary)

    def test_rejects_a_truth_file_with_a_divergent_size(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            truth = output / "outputs" / "truth.jsonl"
            truth.parent.mkdir()
            truth.write_bytes(b"x" * 122)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "artifacts": [
                            {
                                "id": "truth-jsonl",
                                "kind": "truth_jsonl",
                                "path": "outputs/truth.jsonl",
                                "size_bytes": 123,
                                "checksum": {"algorithm": "sha256", "digest": "a" * 64},
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )

            summary = truth_artifact_summary(output)

        self.assertIn("incompleto o es inválido", summary)

    def test_reports_a_broken_manifest_link_without_opening_truth_jsonl(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "outputs").mkdir()
            (output / "outputs" / "truth.jsonl").write_text("not JSONL", encoding="utf-8")
            (output / "manifest.json").write_text(
                json.dumps({"artifacts": []}), encoding="utf-8"
            )

            summary = truth_artifact_summary(output)

        self.assertIn("no enlaza exactamente", summary)

    def test_rejects_a_noncanonical_truth_path(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "artifacts": [
                            {
                                "id": "truth-jsonl",
                                "kind": "truth_jsonl",
                                "path": "outputs/not-truth.jsonl",
                                "size_bytes": 1,
                                "checksum": {"algorithm": "sha256", "digest": "a" * 64},
                            }
                        ]
                    }
                ),
                encoding="utf-8",
            )

            summary = truth_artifact_summary(output)

        self.assertIn("incompleto o es inválido", summary)


if __name__ == "__main__":
    unittest.main()
