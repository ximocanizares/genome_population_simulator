"""Freeze and verify the identity of the external tools the VE bench invokes.

VE.1.1 requires every executable to be pinned before the assay runs: resolved
path, byte size, SHA-256 and the version string the tool reports. Neither this
module nor any other part of the bench installs or downloads a tool.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

LOCK_PATH = Path(__file__).resolve().parents[2] / "environment.lock.json"

# Tools the assay actually invokes. `tabix` and `bgzip` are deliberately absent:
# VE-SNP 0.1 keeps the raw FreeBayes VCF as text and bcftools reads .vcf.gz with
# its own embedded HTSlib, so neither is on the path of the assay.
TOOLS: dict[str, list[str]] = {
    "minimap2": ["--version"],
    "samtools": ["--version"],
    "bcftools": ["--version"],
    "freebayes": ["--version"],
}

VERSION_PATTERNS: dict[str, str] = {
    "minimap2": r"^(\S+)",
    "samtools": r"^samtools\s+(\S+)",
    "bcftools": r"^bcftools\s+(\S+)",
    "freebayes": r"version:\s*(\S+)",
}


class EnvironmentError_(RuntimeError):
    """A tool is missing or does not match the frozen identity."""


@dataclass(frozen=True)
class ToolIdentity:
    name: str
    path: str
    size_bytes: int
    sha256: str
    version: str


def _digest(path: Path) -> tuple[int, str]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        while chunk := handle.read(1 << 20):
            size += len(chunk)
            digest.update(chunk)
    return size, digest.hexdigest()


def _version(name: str, args: list[str]) -> str:
    completed = subprocess.run(
        [name, *args], capture_output=True, text=True, check=False
    )
    if completed.returncode not in (0, 1):
        raise EnvironmentError_(
            f"{name} {' '.join(args)} termino con codigo {completed.returncode}"
        )
    output = completed.stdout or completed.stderr
    pattern = VERSION_PATTERNS[name]
    for line in output.splitlines():
        match = re.search(pattern, line.strip())
        if match:
            return match.group(1)
    raise EnvironmentError_(f"no se pudo leer la version de {name}")


def observe(name: str) -> ToolIdentity:
    """Resolve one tool and read its identity from the file on disk."""
    located = shutil.which(name)
    if located is None:
        raise EnvironmentError_(
            f"{name} no esta instalado; el banco no instala ni descarga herramientas"
        )
    path = Path(located).resolve()
    size, sha256 = _digest(path)
    return ToolIdentity(name, str(path), size, sha256, _version(name, TOOLS[name]))


def observe_all() -> list[ToolIdentity]:
    return [observe(name) for name in TOOLS]


def freeze(lock_path: Path = LOCK_PATH) -> list[ToolIdentity]:
    identities = observe_all()
    document = {
        "document_kind": "ve_environment_lock",
        "schema_version": "1",
        "tools": [asdict(identity) for identity in identities],
    }
    lock_path.write_text(
        json.dumps(document, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    return identities


def verify(lock_path: Path = LOCK_PATH) -> list[str]:
    """Compare the live environment against the lock. Returns the differences."""
    if not lock_path.is_file():
        raise EnvironmentError_(f"no existe {lock_path}; ejecute primero `freeze`")
    document = json.loads(lock_path.read_text(encoding="utf-8"))
    if document.get("document_kind") != "ve_environment_lock":
        raise EnvironmentError_("el fichero no es un lock de entorno VE")
    expected = {entry["name"]: entry for entry in document["tools"]}
    differences: list[str] = []
    for name in TOOLS:
        if name not in expected:
            differences.append(f"{name}: ausente en el lock")
            continue
        try:
            observed = asdict(observe(name))
        except EnvironmentError_ as error:
            differences.append(f"{name}: {error}")
            continue
        for field in ("path", "size_bytes", "sha256", "version"):
            if observed[field] != expected[name][field]:
                differences.append(
                    f"{name}.{field}: lock={expected[name][field]!r} "
                    f"observado={observed[field]!r}"
                )
    for name in expected:
        if name not in TOOLS:
            differences.append(f"{name}: en el lock pero ya no se usa")
    return differences


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ve.environment")
    parser.add_argument("action", choices=("freeze", "verify"))
    parser.add_argument("--lock", type=Path, default=LOCK_PATH)
    arguments = parser.parse_args(argv)
    try:
        if arguments.action == "freeze":
            for identity in freeze(arguments.lock):
                print(f"{identity.name:10s} {identity.version:14s} {identity.sha256}")
            return 0
        differences = verify(arguments.lock)
    except EnvironmentError_ as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    for difference in differences:
        print(f"difiere: {difference}", file=sys.stderr)
    return 1 if differences else 0


if __name__ == "__main__":
    raise SystemExit(main())


# Filters of the frozen FreeBayes command that §9 of the plan declares inert.
# They are not passed explicitly, so a default that changes across versions would
# silently alter the metrics. VE.1.1 re-checks them against the frozen binary.
INERT_FREEBAYES_DEFAULTS: dict[str, str] = {
    "--min-supporting-allele-qsum": "0",
    "--min-supporting-mapping-qsum": "0",
    "--min-alternate-qsum": "0",
    "--min-coverage": "0",
    "--read-max-mismatch-fraction": "1.0",
}


def freebayes_defaults() -> dict[str, str]:
    """Read the default value FreeBayes reports for each documented option."""
    completed = subprocess.run(
        ["freebayes", "--help"], capture_output=True, text=True, check=False
    )
    text = completed.stdout or completed.stderr
    defaults: dict[str, str] = {}
    option: str | None = None
    for line in text.splitlines():
        found = re.search(r"(--[a-z0-9-]+)", line)
        if found and not line.startswith(" " * 8):
            option = found.group(1)
        if option is not None:
            default = re.search(r"default:\s*(\S+)", line)
            if default:
                defaults[option] = default.group(1)
                option = None
    return defaults
