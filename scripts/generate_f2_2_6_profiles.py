#!/usr/bin/env python3
"""Regenera los perfiles compactos fijados por DTC-068.

No es parte del núcleo científico en Python: solamente materializa los
recursos versionados y rechaza una fuente eLife que no sea la fijada.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESOURCE = ROOT / "resources" / "f2_2_6"
INPUT = RESOURCE / "elife-24284-fig1-data1-v2.txt"
INPUT_SHA256 = "970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f"
INPUT_SIZE = 3644
POPULATIONS = ("AFR", "AMR", "EUR", "EAS", "SAS")
BASES = "ACGT"


def complement(base: str) -> str:
    return {"A": "T", "C": "G", "G": "C", "T": "A"}[base]


def write_sfs() -> None:
    rows = (
        ("interval", "0.00019968051118210862", "0.005", "64", "log_uniform"),
        ("interval", "0.005", "0.05", "12", "log_uniform"),
        ("interval", "0.05", "1", "8", "linear_uniform"),
    )
    lines = ["kind\tmin_alt_af\tmax_alt_af\tweight\tsampling"]
    lines.extend("\t".join(row) for row in rows)
    (RESOURCE / "human_1000g_global_sfs_approx_v1.tsv").write_text("\n".join(lines) + "\n")


def write_harris() -> None:
    metadata = json.loads((RESOURCE / "profile_sources.json").read_text())
    source = metadata["profiles"]["human_harris2017_snv12_approx_v1"]
    if (
        source["license"] != "CC-BY-4.0"
        or source["classification"] != "approximate_descriptive_not_universal_rate"
        or source["input_size"] != INPUT_SIZE
        or source["input_sha256"] != INPUT_SHA256
    ):
        raise SystemExit("metadatos de atribución o integridad eLife no fijados")
    data = INPUT.read_bytes()
    if len(data) != INPUT_SIZE or hashlib.sha256(data).hexdigest() != INPUT_SHA256:
        raise SystemExit("el recurso eLife no coincide con tamaño/SHA-256 fijados")
    lines = data.decode("ascii").splitlines()
    header = lines[0].split()
    contexts = header[1:]
    if len(contexts) != 96 or tuple(line.split()[0] for line in lines[1:]) != POPULATIONS:
        raise SystemExit("estructura eLife inesperada")
    class_totals: dict[str, float] = {f"{ref}>{alt}": 0.0 for ref in "CT" for alt in BASES if alt != ref}
    for line in lines[1:]:
        fields = line.split()
        values = [float(value) for value in fields[1:]]
        total = sum(values)
        for context, value in zip(contexts, values, strict=True):
            triplet, alternate = context.split("_")
            reference = triplet[1]
            if reference in "AG":
                reference, alternate = complement(reference), complement(alternate)
            class_totals[f"{reference}>{alternate}"] += value / total / len(POPULATIONS)
    directed: dict[str, float] = {f"{ref}>{alt}": 0.0 for ref in BASES for alt in BASES if alt != ref}
    for change, weight in class_totals.items():
        reference, alternate = change.split(">")
        directed[change] += weight / 2.0
        directed[f"{complement(reference)}>{complement(alternate)}"] += weight / 2.0
    total = sum(directed.values())
    order = tuple(f"{ref}>{alt}" for ref in BASES for alt in BASES if alt != ref)
    output = ["substitution\tweight"]
    output.extend(f"{change}\t{directed[change] / total:.17g}" for change in order)
    (RESOURCE / "human_harris2017_snv12_approx_v1.tsv").write_text("\n".join(output) + "\n")


if __name__ == "__main__":
    write_sfs()
    write_harris()
