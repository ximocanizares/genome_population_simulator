from __future__ import annotations

import argparse
import importlib.util
import sys
import unittest
from pathlib import Path
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
LAUNCHER_PATH = ROOT / "scripts" / "run_marimo.py"
SPEC = importlib.util.spec_from_file_location("gps_run_marimo", LAUNCHER_PATH)
assert SPEC is not None and SPEC.loader is not None
launcher = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(launcher)


class LauncherTests(unittest.TestCase):
    def test_default_command_is_loopback_only_and_does_not_start_a_simulation(self) -> None:
        command = launcher.marimo_command(launcher.DEFAULT_PORT, no_open=True)

        self.assertEqual(command[:4], [sys.executable, "-m", "marimo", "run"])
        self.assertEqual(command[command.index("--host") + 1], "127.0.0.1")
        self.assertEqual(command[command.index("--port") + 1], "2718")
        self.assertIn("--headless", command)
        self.assertNotIn("genome-population-simulator", command)

    def test_browser_mode_omits_headless_and_accepts_an_alternate_port(self) -> None:
        command = launcher.marimo_command(4721, no_open=False)

        self.assertNotIn("--headless", command)
        self.assertEqual(command[command.index("--port") + 1], "4721")

    def test_privileged_out_of_range_and_non_numeric_ports_are_rejected(self) -> None:
        for value in ["0", "1023", "65536", "not-a-port"]:
            with self.subTest(value=value):
                with self.assertRaises(argparse.ArgumentTypeError):
                    launcher.local_port(value)

    def test_occupied_port_is_diagnosed_without_starting_marimo(self) -> None:
        with mock.patch.object(launcher, "port_is_available", return_value=False):
            with mock.patch.object(launcher, "run_marimo") as run:
                result = launcher.main(["--port", "4721", "--no-open"])

        self.assertEqual(result, 2)
        run.assert_not_called()

    def test_available_port_starts_exactly_the_constructed_process(self) -> None:
        port = 4721
        with mock.patch.object(launcher, "port_is_available", return_value=True):
            with mock.patch.object(launcher, "run_marimo", return_value=7) as run:
                result = launcher.main(["--port", str(port), "--no-open"])

        self.assertEqual(result, 7)
        run.assert_called_once_with(launcher.marimo_command(port, no_open=True))

    def test_keyboard_interrupt_terminates_the_child_without_a_traceback(self) -> None:
        process = mock.Mock()
        process.wait.side_effect = [KeyboardInterrupt, 0]
        process.poll.return_value = None
        with mock.patch.object(launcher.subprocess, "Popen", return_value=process):
            result = launcher.run_marimo(["marimo", "run", "app.py"])

        self.assertEqual(result, 130)
        process.terminate.assert_called_once_with()
        process.kill.assert_not_called()


if __name__ == "__main__":
    unittest.main()
