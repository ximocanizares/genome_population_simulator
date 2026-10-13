from __future__ import annotations

import os
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
    default_binary,
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
                {
                    "code": "GPS-CONF-008",
                    "message": "ruta vacía",
                    "context": {
                        "fields": "population_snvs.input_path",
                        "location": "8:14",
                    },
                },
            ]
        }

        summary = diagnostic_summary(document)

        self.assertIn("General/Ejecución · schema_version · 2:1", summary)
        self.assertIn("Genoma · ploidy · 4:1", summary)
        self.assertIn("Poblaciones · population_snvs.input_path · 8:14", summary)
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
        self.assertFalse(session.matches("a = 1\n", binary, base, True))

    def test_validated_session_accepts_equivalent_paths_with_edge_whitespace(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        base = ROOT.resolve()
        session = ValidatedSession("a = 1\n", binary, SUPPORTED_BINARY_VERSION, base)

        self.assertTrue(
            session.matches("a = 1\n", Path(f"  {binary}  "), Path(f"  {base}  "))
        )

    def test_edge_whitespace_does_not_turn_an_absolute_output_path_into_a_relative_one(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        session = ValidatedSession("a = 1\n", binary, SUPPORTED_BINARY_VERSION, ROOT)
        with tempfile.TemporaryDirectory() as directory:
            destination = Path(directory) / "output with spaces"
            with mock.patch.object(
                subprocess,
                "run",
                return_value=subprocess.CompletedProcess([], 0, "", ""),
            ):
                result = run_toml(
                    session,
                    "a = 1\n",
                    binary,
                    ROOT,
                    Path(f"  {destination}  "),
                )
        self.assertEqual(result.output_directory, destination)

    def test_empty_matrix_authorization_is_an_argument_and_invalidates_the_session(self) -> None:
        binary = (ROOT / "target/debug/genome-population-simulator").resolve()
        session = ValidatedSession(
            "a = 1\n", binary, SUPPORTED_BINARY_VERSION, ROOT, True
        )
        self.assertTrue(session.matches("a = 1\n", binary, ROOT, True))
        with mock.patch.object(
            subprocess, "run", return_value=subprocess.CompletedProcess([], 0, "", "")
        ) as run:
            with tempfile.TemporaryDirectory() as directory:
                result = run_toml(
                    session,
                    "a = 1\n",
                    binary,
                    ROOT,
                    Path(directory) / "new-output",
                    allow_empty_snv_matrix=True,
                )
        self.assertIn("--allow-empty-snv-matrix", result.command)
        self.assertEqual(result.command[1], "run")
        self.assertLess(result.command.index("--allow-empty-snv-matrix"), result.command.index("--diagnostics"))
        run.assert_called_once()

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

    def test_e147_f1_version_is_rejected_before_any_form_action(self) -> None:
        """F1 e14705a reported 0.1.0; F2's UI must reject it at the version gate."""
        with tempfile.TemporaryDirectory() as directory:
            work = Path(directory)
            legacy = work / "genome-population-simulator"
            legacy.write_text(
                "#!/bin/sh\n"
                'echo "$1" >> "$(dirname "$0")/calls"\n'
                'if [ "$1" = "--version" ]; then\n'
                "  echo 'genome-population-simulator 0.1.0'\n"
                "  exit 0\n"
                "fi\n"
                "echo 'form action was reached' >&2\n"
                "exit 99\n",
                encoding="utf-8",
            )
            legacy.chmod(0o755)

            with self.assertRaisesRegex(ProcessAdapterError, "0.2.0"):
                validate_toml(legacy, 'config_kind = "general"', work)

            self.assertEqual((work / "calls").read_text(encoding="utf-8"), "--version\n")

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
    def test_real_binary_reports_an_unknown_top_level_field(self) -> None:
        invalid = render_general_toml(default_form()).replace(
            "[reference]", "unknown_top_level = true\n[reference]", 1
        )
        with tempfile.TemporaryDirectory() as directory:
            result, session = validate_toml(self.BINARY, invalid, Path(directory))

        self.assertEqual(result.returncode, 2)
        self.assertIsNone(session)
        self.assertIsNotNone(result.diagnostics)
        self.assertEqual(
            result.diagnostics["diagnostics"][0]["code"], "GPS-CONF-006"
        )

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_real_binary_reports_a_duplicate_top_level_field(self) -> None:
        invalid = render_general_toml(default_form()).replace(
            "[reference]", "ploidy = 2\n[reference]", 1
        )
        with tempfile.TemporaryDirectory() as directory:
            result, session = validate_toml(self.BINARY, invalid, Path(directory))

        self.assertEqual(result.returncode, 2)
        self.assertIsNone(session)
        self.assertIsNotNone(result.diagnostics)
        self.assertEqual(
            result.diagnostics["diagnostics"][0]["code"], "GPS-CONF-007"
        )

    @unittest.skipUnless(BINARY.exists(), "compile the Rust binary before this integration test")
    def test_empty_population_matrix_path_keeps_its_qualified_json_field_and_location(self) -> None:
        invalid = render_general_toml(default_form()).replace(
            "[reference]",
            '[population_snvs]\nmode = "supplied"\ninput_path = ""\n[reference]',
            1,
        )
        with tempfile.TemporaryDirectory() as directory:
            result, session = validate_toml(self.BINARY, invalid, Path(directory))

        self.assertEqual(result.returncode, 2)
        self.assertIsNone(session)
        self.assertIsNotNone(result.diagnostics)
        context = result.diagnostics["diagnostics"][0]["context"]
        self.assertEqual(context["fields"], "population_snvs.input_path")
        self.assertIn("location", context)
        self.assertIn(
            "Poblaciones · population_snvs.input_path", diagnostic_summary(result.diagnostics)
        )


class DefaultBinaryTests(unittest.TestCase):
    """El binario propuesto por omisión es el optimizado salvo que esté desfasado."""

    def _build(self, root: Path, kind: str, mtime: float) -> Path:
        directory = (
            root / "target" / "x86_64-unknown-linux-gnu" / "release"
            if kind == "release"
            else root / "target" / "debug"
        )
        directory.mkdir(parents=True, exist_ok=True)
        binary = directory / SUPPORTED_BINARY_NAME
        binary.write_bytes(b"")
        os.utime(binary, (mtime, mtime))
        return binary

    def test_the_release_binary_is_chosen_when_it_is_not_older(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            self._build(root, "debug", 1_000)
            release = self._build(root, "release", 2_000)
            self.assertEqual(default_binary(root), release)

    def test_a_release_binary_older_than_the_debug_one_is_not_chosen(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            self._build(root, "release", 1_000)
            debug = self._build(root, "debug", 2_000)
            self.assertEqual(default_binary(root), debug)

    def test_without_a_release_binary_the_debug_path_is_proposed(self) -> None:
        with tempfile.TemporaryDirectory() as raw:
            root = Path(raw)
            self.assertEqual(
                default_binary(root), root / "target" / "debug" / SUPPORTED_BINARY_NAME
            )
            release = self._build(root, "release", 1_000)
            self.assertEqual(default_binary(root), release)


if __name__ == "__main__":
    unittest.main()
