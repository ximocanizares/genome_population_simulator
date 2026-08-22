#!/usr/bin/env python3
"""Oráculo independiente de los recursos F2.2.6, sin reutilizar el generador."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RESOURCE = ROOT / "resources" / "f2_2_6"
SOURCE = RESOURCE / "elife-24284-fig1-data1-v2.txt"
EXPECTED_SOURCE_SHA256 = "970cc71578fc008ef4e0d75c7c460ab874073da7c5b5226c526b387627c97a7f"
EXPECTED_SOURCE_SIZE = 3644
EXPECTED_SFS_SHA256 = "0ae67f89380293a136dd265099785790b32588a6234234e702c280f15276fc7e"
EXPECTED_HARRIS_SHA256 = "ea2f515124eafe7ed174930f7a71c4f524f4c7c9d22bfc4d3b0051dff55db0b2"


def reverse_complement_change(change: str) -> str:
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    return f"{complement[change[0]]}>{complement[change[2]]}"


def check_source() -> tuple[list[str], list[list[float]]]:
    metadata = json.loads((RESOURCE / "profile_sources.json").read_text())
    harris = metadata["profiles"]["human_harris2017_snv12_approx_v1"]
    sfs = metadata["profiles"]["human_1000g_global_sfs_approx_v1"]
    assert metadata["schema_version"] == "1"
    assert harris["license"] == "CC-BY-4.0"
    assert harris["classification"] == "approximate_descriptive_not_universal_rate"
    assert harris["input_size"] == EXPECTED_SOURCE_SIZE
    assert harris["input_sha256"] == EXPECTED_SOURCE_SHA256
    assert harris["output_sha256"] == EXPECTED_HARRIS_SHA256
    assert sfs["literals"] == {"minimum_alt_af": "1/5008", "weights": [64, 12, 8]}
    assert sfs["classification"] == "approximate_descriptive_not_universal_rate"
    assert sfs["output_sha256"] == EXPECTED_SFS_SHA256
    data = SOURCE.read_bytes()
    assert len(data) == EXPECTED_SOURCE_SIZE
    assert hashlib.sha256(data).hexdigest() == EXPECTED_SOURCE_SHA256
    lines = data.decode("ascii").splitlines()
    header = lines[0].split()
    assert len(header) == 97 and header[0] == "Pop"
    rows = [(line.split()[0], [float(value) for value in line.split()[1:]]) for line in lines[1:]]
    assert [name for name, _ in rows] == ["AFR", "AMR", "EUR", "EAS", "SAS"]
    assert all(len(values) == 96 for _, values in rows)
    return header[1:], [values for _, values in rows]


def check_sfs() -> None:
    path = RESOURCE / "human_1000g_global_sfs_approx_v1.tsv"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == EXPECTED_SFS_SHA256
    assert path.read_text().splitlines() == [
        "kind\tmin_alt_af\tmax_alt_af\tweight\tsampling",
        "interval\t0.00019968051118210862\t0.005\t64\tlog_uniform",
        "interval\t0.005\t0.05\t12\tlog_uniform",
        "interval\t0.05\t1\t8\tlinear_uniform",
    ]


def check_harris() -> None:
    contexts, population_rows = check_source()
    path = RESOURCE / "human_harris2017_snv12_approx_v1.tsv"
    assert hashlib.sha256(path.read_bytes()).hexdigest() == EXPECTED_HARRIS_SHA256
    rows = [line.split("\t") for line in path.read_text().splitlines()]
    assert rows[0] == ["substitution", "weight"] and len(rows) == 13
    weights = {change: float(weight) for change, weight in rows[1:]}
    assert tuple(weights) == tuple(f"{ref}>{alt}" for ref in "ACGT" for alt in "ACGT" if alt != ref)
    assert abs(sum(weights.values()) - 1.0) < 1e-15
    for change, weight in weights.items():
        assert abs(weight - weights[reverse_complement_change(change)]) < 1e-15
    # Recalcula independientemente: normalización por población, media igual,
    # colapso a C/T y reflexión complemento-inverso hacia las doce clases.
    complement = {"A": "T", "C": "G", "G": "C", "T": "A"}
    six = {f"{reference}>{alternate}": 0.0 for reference in "CT" for alternate in "ACGT" if alternate != reference}
    for row in population_rows:
        row_total = sum(row)
        for context, count in zip(contexts, row, strict=True):
            context_bases, alternate = context.split("_")
            reference = context_bases[1]
            if reference not in "CT":
                reference, alternate = complement[reference], complement[alternate]
            six[f"{reference}>{alternate}"] += count / row_total / 5.0
    expected = {f"{reference}>{alternate}": 0.0 for reference in "ACGT" for alternate in "ACGT" if alternate != reference}
    for change, weight in six.items():
        expected[change] += weight / 2.0
        expected[reverse_complement_change(change)] += weight / 2.0
    expected_total = sum(expected.values())
    for change in expected:
        assert abs(weights[change] - expected[change] / expected_total) < 1e-15


if __name__ == "__main__":
    check_sfs()
    check_harris()
