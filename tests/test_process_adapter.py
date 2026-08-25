"""Regression coverage for the UI's exact binary-version gate."""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "ui"))

from process_adapter import ProcessAdapterError, validate_toml


def test_e147_f1_version_is_rejected_before_any_form_action(tmp_path: Path) -> None:
    """F1 e14705a reported 0.1.0; F2's UI must reject it at the version gate."""
    legacy = tmp_path / "genome-population-simulator"
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

    with pytest.raises(ProcessAdapterError, match="0.2.0"):
        validate_toml(legacy, 'config_kind = "general"', tmp_path)

    assert (tmp_path / "calls").read_text(encoding="utf-8") == "--version\n"
