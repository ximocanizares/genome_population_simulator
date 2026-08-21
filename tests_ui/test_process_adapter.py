from __future__ import annotations

import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "ui"))

from config_editor import default_form, render_general_toml  # noqa: E402
from process_adapter import (  # noqa: E402
    ProcessAdapterError,
    SUPPORTED_BINARY_NAME,
    SUPPORTED_BINARY_VERSION,
    ValidatedSession,
    check_binary_version,
    diagnostic_summary,
    run_toml,
    validate_toml,
)


class ProcessAdapterTests(unittest.TestCase):
    BINARY = ROOT / "target" / "debug" / SUPPORTED_BINARY_NAME

    def test_rust_diagnostics_are_located_without_reinterpreting_them(self) -> None:
        document = {
            "diagnostics": [
                {
                    "code": "GPS-CONF-005",
                    "message": "versión no admitida",
                    "context": {"field": "schema_version", "location": "2:1"},
                },
                {
                    "code": "GPS-CONF-006",
                    "message": "campo desconocido",
                    "context": {"field": "ploidy", "location": "4:1"},
                },
            ]
        }

        summary = diagnostic_summary(document)

        self.assertIn("General/Ejecución · schema_version · 2:1", summary)
        self.assertIn("Genoma · ploidy · 4:1", summary)
        self.assertIn("GPS-CONF-006: campo desconocido", summary)

    def test_version_check_accepts_only_the_supported_literal(self) -> None:
        binary = Path("target/debug/genome-population-simulator")
        expected = f"{SUPPORTED_BINARY_NAME} {SUPPORTED_BINARY_VERSION}\n"
        with mock.patch.object(
            subprocess,
            "run",
            return_value=subprocess.CompletedProcess([], 0, stdout=expected, stderr=""),
        ):
            resolved, version = check_binary_version(binary)
        self.assertEqual(resolved, binary.resolve())
        self.assertEqual(version, SUPPORTED_BINARY_VERSION)

        for stdout in ["wrong 0.1.0\n", f"{SUPPORTED_BINARY_NAME} 9.0.0\n"]:
            with self.subTest(stdout=stdout):
                with mock.patch.object(
                    subprocess,
                    "run",
                    return_value=subprocess.CompletedProcess(
                        [], 0, stdout=stdout, stderr=""
                    ),
                ):
                    with self.assertRaises(ProcessAdapterError):
                        check_binary_version(binary)

    def test_validated_session_is_invalidated_by_each_authoritative_input(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        base = ROOT.resolve()
        session = ValidatedSession("a = 1\n", binary, SUPPORTED_BINARY_VERSION, base)

        self.assertTrue(session.matches("a = 1\n", binary, base))
        self.assertFalse(session.matches("a = 2\n", binary, base))
        self.assertFalse(session.matches("a = 1\n", binary.with_name("other"), base))
        self.assertFalse(session.matches("a = 1\n", binary, base.parent))

    def test_run_rejects_existing_destination_before_invoking_a_process(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        session = ValidatedSession("a = 1\n", binary, SUPPORTED_BINARY_VERSION, ROOT)
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory)
            with mock.patch.object(subprocess, "run") as run:
                with self.assertRaises(ProcessAdapterError):
                    run_toml(session, "a = 1\n", binary, ROOT, destination)
        run.assert_not_called()

    def test_run_rejects_an_edited_document_before_invoking_a_process(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        session = ValidatedSession("a = 1\n", binary, SUPPORTED_BINARY_VERSION, ROOT)
        with mock.patch.object(subprocess, "run") as run:
            with self.assertRaises(ProcessAdapterError):
                run_toml(
                    session,
                    "a = 2\n",
                    binary,
                    ROOT,
                    ROOT / "new-output",
                )
        run.assert_not_called()

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_real_binary_validates_and_runs_without_loading_result_artifacts(self) -> None:
        toml = render_general_toml(default_form())
        with tempfile.TemporaryDirectory() as directory:
            base = Path(directory)
            validation, session = validate_toml(self.BINARY, toml, base)
            self.assertTrue(validation.succeeded)
            self.assertIsNotNone(session)

            destination = base / "result"
            execution = run_toml(
                session,
                toml,
                self.BINARY,
                base,
                destination,
            )

            self.assertTrue(execution.succeeded)
            self.assertEqual(execution.output_directory, destination)
            self.assertTrue((destination / "manifest.json").is_file())
            self.assertFalse(hasattr(execution, "artifacts"))
            self.assertEqual(execution.command[1], "run")
            self.assertIn("--output", execution.command)

            relative_destination = Path("second-result")
            second_execution = run_toml(
                session,
                toml,
                self.BINARY,
                base,
                relative_destination,
            )
            self.assertEqual(second_execution.output_directory, base / relative_destination)

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_real_binary_returns_the_structured_rust_diagnostic(self) -> None:
        invalid = render_general_toml(default_form()).replace(
            "[reference]", "ploidy = 2\n[reference]", 1
        )
        with tempfile.TemporaryDirectory() as directory:
            result, session = validate_toml(self.BINARY, invalid, Path(directory))

        self.assertEqual(result.returncode, 2)
        self.assertIsNone(session)
        self.assertIsNotNone(result.diagnostics)
        self.assertEqual(
            result.diagnostics["diagnostics"][0]["code"], "GPS-CONF-006"
        )


if __name__ == "__main__":
    unittest.main()
