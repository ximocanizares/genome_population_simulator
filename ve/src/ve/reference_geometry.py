"""Publish the geometry the reduction actually produced.

VE-D14 makes the source genome and the reduction geometry parameters of each
assay, not constants of the contract, and requires VE.1.2 to publish the
realised geometry read from the FASTA index and the coordinate map rather than
assumed from the configuration.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path


class GeometryError(RuntimeError):
    """The published artefacts do not describe a consistent geometry."""


@dataclass(frozen=True)
class Contig:
    contig_id: str
    length_bases: int
    block_count: int
    source_intervals: list[dict[str, int | str]]


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_fai(path: Path) -> dict[str, int]:
    lengths: dict[str, int] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        fields = line.split("\t")
        lengths[fields[0]] = int(fields[1])
    return lengths


def describe(fasta: Path, coordinate_map: Path) -> dict[str, object]:
    """Build the realised-geometry document from the published artefacts."""
    index = fasta.with_suffix(fasta.suffix + ".fai")
    if not index.is_file():
        raise GeometryError(f"falta el indice {index}; ejecute `samtools faidx`")
    lengths = _read_fai(index)
    document = json.loads(coordinate_map.read_text(encoding="utf-8"))
    if document.get("document_kind") != "reference_reduction_map":
        raise GeometryError("el mapa no es un reference_reduction_map")

    blocks: dict[str, list[dict[str, int | str]]] = defaultdict(list)
    for segment in sorted(document["segments"], key=lambda s: s["materialization_order"]):
        blocks[segment["sequence_id"]].append(
            {
                "source_start": segment["source_interval"]["start"],
                "source_end": segment["source_interval"]["end"],
                "orientation": segment["orientation"],
                "sequence_class": segment["sequence_class"],
            }
        )

    if set(blocks) != set(lengths):
        raise GeometryError("el mapa y el indice no describen los mismos contigs")

    contigs = [
        Contig(name, lengths[name], len(blocks[name]), blocks[name])
        for name in sorted(lengths)
    ]
    for contig in contigs:
        covered = sum(
            int(b["source_end"]) - int(b["source_start"]) + 1 for b in contig.source_intervals
        )
        if covered != contig.length_bases:
            raise GeometryError(
                f"{contig.contig_id}: los bloques suman {covered} y el indice dice "
                f"{contig.length_bases}"
            )

    return {
        "document_kind": "ve_realised_reference_geometry",
        "schema_version": "1",
        "derived_reference_id": document["derived_reference_id"],
        "source_reference_id": document["source_reference_id"],
        "fasta_sha256": _digest(fasta),
        "contig_count": len(contigs),
        "total_bases": sum(c.length_bases for c in contigs),
        "block_count": sum(c.block_count for c in contigs),
        "contigs": [asdict(c) for c in contigs],
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ve.reference_geometry")
    parser.add_argument("fasta", type=Path)
    parser.add_argument("coordinate_map", type=Path)
    parser.add_argument("--out", type=Path)
    arguments = parser.parse_args(argv)
    document = describe(arguments.fasta, arguments.coordinate_map)
    text = json.dumps(document, indent=2, ensure_ascii=False) + "\n"
    if arguments.out:
        arguments.out.write_text(text, encoding="utf-8")
        print(
            f"{document['contig_count']} contigs, {document['total_bases']} bases, "
            f"{document['block_count']} bloques -> {arguments.out}"
        )
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
