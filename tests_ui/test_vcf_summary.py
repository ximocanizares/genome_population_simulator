from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from process_adapter import vcf_artifact_summary  # noqa: E402


def artifact(artifact_id: str, kind: str, path: str, size: int) -> dict[str, object]:
    return {
        "id": artifact_id,
        "kind": kind,
        "path": path,
        "size_bytes": size,
        "checksum": {"algorithm": "sha256", "digest": "a" * 64},
    }


class VcfArtifactSummaryTests(unittest.TestCase):
    def test_reports_the_cataloged_text_vcf_without_opening_it(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            vcf = output / "outputs" / "truth.vcf"
            vcf.parent.mkdir()
            vcf.write_bytes(b"x" * 123)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "vcf_outputs": [
                            {
                                "state": "complete",
                                "mode": "text",
                                "vcf_artifact_id": "truth-vcf",
                                "index_kind": None,
                                "index_artifact_id": None,
                            }
                        ],
                        "artifacts": [
                            artifact("truth-vcf", "truth_vcf", "outputs/truth.vcf", 123)
                        ],
                    }
                ),
                encoding="utf-8",
            )

            summary = vcf_artifact_summary(output)

        self.assertIn("truth.vcf", summary)
        self.assertIn("ruta=outputs/truth.vcf", summary)
        self.assertIn("tamaño=123 bytes", summary)

    def test_reports_bgzf_and_exactly_the_cataloged_index(self) -> None:
        for index_kind in ("tbi", "csi"):
            with self.subTest(index_kind=index_kind), tempfile.TemporaryDirectory() as directory:
                output = Path(directory)
                vcf = output / "outputs" / "truth.vcf.gz"
                index = output / "outputs" / f"truth.vcf.gz.{index_kind}"
                vcf.parent.mkdir()
                vcf.write_bytes(b"v" * 123)
                index.write_bytes(b"i" * 45)
                (output / "manifest.json").write_text(
                    json.dumps(
                        {
                            "vcf_outputs": [
                                {
                                    "state": "complete",
                                    "mode": "bgzf",
                                    "vcf_artifact_id": "truth-vcf",
                                    "index_kind": index_kind,
                                    "index_artifact_id": "truth-vcf-index",
                                }
                            ],
                            "artifacts": [
                                artifact("truth-vcf", "truth_vcf", "outputs/truth.vcf.gz", 123),
                                artifact(
                                    "truth-vcf-index",
                                    f"vcf_{index_kind}_index",
                                    f"outputs/truth.vcf.gz.{index_kind}",
                                    45,
                                ),
                            ],
                        }
                    ),
                    encoding="utf-8",
                )

                summary = vcf_artifact_summary(output)

            self.assertIn("truth.vcf.gz", summary)
            self.assertIn(f"índice {index_kind.upper()}", summary)
            self.assertIn(f"truth.vcf.gz.{index_kind}", summary)

    def test_omitted_vcf_is_reported_without_inventing_an_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "manifest.json").write_text(
                json.dumps({"vcf_outputs": [], "artifacts": []}), encoding="utf-8"
            )

            summary = vcf_artifact_summary(output)

        self.assertIn("VCF no solicitado", summary)

    def test_rejects_a_bgzf_link_with_the_wrong_index_kind(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            (output / "manifest.json").write_text(
                json.dumps(
                    {
                        "vcf_outputs": [
                            {
                                "state": "complete",
                                "mode": "bgzf",
                                "vcf_artifact_id": "truth-vcf",
                                "index_kind": "tbi",
                                "index_artifact_id": "truth-vcf-index",
                            }
                        ],
                        "artifacts": [],
                    }
                ),
                encoding="utf-8",
            )

            summary = vcf_artifact_summary(output)

        self.assertIn("incompleto o es inválido", summary)


if __name__ == "__main__":
    unittest.main()
