"""VE.1.1: the frozen environment must be verifiable, and tampering must fail."""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ve.environment import (
    INERT_FREEBAYES_DEFAULTS,
    LOCK_PATH,
    EnvironmentError_,
    freebayes_defaults,
    freeze,
    observe,
    verify,
)


class FrozenEnvironment(unittest.TestCase):
    def test_lock_is_present_and_well_formed(self) -> None:
        document = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
        self.assertEqual(document["document_kind"], "ve_environment_lock")
        self.assertEqual(document["schema_version"], "1")
        names = {entry["name"] for entry in document["tools"]}
        self.assertEqual(names, {"minimap2", "samtools", "bcftools", "freebayes"})
        for entry in document["tools"]:
            self.assertRegex(entry["sha256"], r"^[0-9a-f]{64}$")
            self.assertGreater(entry["size_bytes"], 0)
            self.assertTrue(Path(entry["path"]).is_absolute())

    def test_valid_lock_verifies_without_differences(self) -> None:
        self.assertEqual(verify(LOCK_PATH), [])

    def test_frozen_versions_match_the_approved_ones(self) -> None:
        """VE-D05 congela estas versiones; el banco no debe ejecutar con otras."""
        document = json.loads(LOCK_PATH.read_text(encoding="utf-8"))
        frozen = {entry["name"]: entry["version"] for entry in document["tools"]}
        self.assertEqual(frozen["minimap2"], "2.26-r1175")
        self.assertEqual(frozen["samtools"], "1.21")
        self.assertEqual(frozen["bcftools"], "1.21")
        self.assertEqual(frozen["freebayes"], "v1.3.10")


class TamperedEnvironment(unittest.TestCase):
    def _lock_copy(self, directory: str) -> Path:
        target = Path(directory) / "environment.lock.json"
        freeze(target)
        return target

    def test_altered_digest_is_reported(self) -> None:
        with TemporaryDirectory() as directory:
            lock = self._lock_copy(directory)
            document = json.loads(lock.read_text(encoding="utf-8"))
            document["tools"][0]["sha256"] = "0" * 64
            lock.write_text(json.dumps(document), encoding="utf-8")
            differences = verify(lock)
            self.assertTrue(any(".sha256:" in item for item in differences))

    def test_altered_version_is_reported(self) -> None:
        with TemporaryDirectory() as directory:
            lock = self._lock_copy(directory)
            document = json.loads(lock.read_text(encoding="utf-8"))
            document["tools"][0]["version"] = "0.0-impossible"
            lock.write_text(json.dumps(document), encoding="utf-8")
            differences = verify(lock)
            self.assertTrue(any(".version:" in item for item in differences))

    def test_missing_tool_in_lock_is_reported(self) -> None:
        with TemporaryDirectory() as directory:
            lock = self._lock_copy(directory)
            document = json.loads(lock.read_text(encoding="utf-8"))
            dropped = document["tools"].pop()["name"]
            lock.write_text(json.dumps(document), encoding="utf-8")
            self.assertIn(f"{dropped}: ausente en el lock", verify(lock))

    def test_absent_lock_is_an_error(self) -> None:
        with TemporaryDirectory() as directory:
            with self.assertRaises(EnvironmentError_):
                verify(Path(directory) / "no-existe.json")

    def test_wrong_document_kind_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            lock = Path(directory) / "environment.lock.json"
            lock.write_text(
                json.dumps({"document_kind": "otra_cosa"}), encoding="utf-8"
            )
            with self.assertRaises(EnvironmentError_):
                verify(lock)

    def test_absent_tool_is_reported_without_installing_it(self) -> None:
        """El banco no instala ni descarga: una ausencia es un error, no un arreglo."""
        with self.assertRaises(EnvironmentError_) as raised:
            observe("herramienta-que-no-existe")
        self.assertIn("no esta instalado", str(raised.exception))

    def test_tool_outside_the_catalogue_is_rejected(self) -> None:
        """Una herramienta instalada pero no declarada no entra por la puerta de atras."""
        with self.assertRaises(KeyError):
            observe("ls")


class FrozenCallerDefaults(unittest.TestCase):
    """§9 declara inertes varios filtros que el comando no pasa explicitamente.

    Un valor predeterminado que cambiase de version alteraria las metricas sin
    avisar, asi que VE.1.1 los vuelve a comprobar contra el binario congelado.
    """

    def test_inert_filters_remain_inert(self) -> None:
        observed = freebayes_defaults()
        for option, expected in INERT_FREEBAYES_DEFAULTS.items():
            with self.subTest(option=option):
                self.assertEqual(observed.get(option), expected)

    def test_mapping_quality_filter_is_the_one_declared(self) -> None:
        """`-m 1` se pasa explicito porque no es inerte: descarta MAPQ 0."""
        self.assertEqual(freebayes_defaults().get("--min-mapping-quality"), "1")


if __name__ == "__main__":
    unittest.main()
