"""VE.1.2: la geometria publicada debe leerse de los artefactos, no suponerse."""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ve.reference_geometry import GeometryError, describe

ASSAY = Path(__file__).resolve().parents[1] / "assays" / "ve-snp-0.1"
GEOMETRY = ASSAY / "reference-geometry.json"


class PublishedGeometry(unittest.TestCase):
    def setUp(self) -> None:
        if not GEOMETRY.is_file():
            self.skipTest("la geometria aun no se ha publicado")
        self.document = json.loads(GEOMETRY.read_text(encoding="utf-8"))

    def test_document_is_well_formed(self) -> None:
        self.assertEqual(self.document["document_kind"], "ve_realised_reference_geometry")
        self.assertEqual(self.document["derived_reference_id"], "ve-snp-tomato-500k")
        self.assertEqual(self.document["source_reference_id"], "tomato-sl4")
        self.assertRegex(self.document["fasta_sha256"], r"^[0-9a-f]{64}$")

    def test_total_is_the_exact_target(self) -> None:
        """VE-D14 exige `target_bases`; solo el objetivo exacto da 500.000."""
        self.assertEqual(self.document["total_bases"], 500_000)

    def test_contig_lengths_sum_to_the_total(self) -> None:
        total = sum(c["length_bases"] for c in self.document["contigs"])
        self.assertEqual(total, self.document["total_bases"])

    def test_every_contig_is_a_mosaic_of_declared_blocks(self) -> None:
        for contig in self.document["contigs"]:
            covered = sum(
                b["source_end"] - b["source_start"] + 1 for b in contig["source_intervals"]
            )
            with self.subTest(contig=contig["contig_id"]):
                self.assertEqual(covered, contig["length_bases"])
                self.assertEqual(len(contig["source_intervals"]), contig["block_count"])

    def test_blocks_do_not_overlap_within_a_source_sequence(self) -> None:
        for contig in self.document["contigs"]:
            spans = sorted(
                (b["source_start"], b["source_end"]) for b in contig["source_intervals"]
            )
            for (_, previous_end), (start, _) in zip(spans, spans[1:]):
                with self.subTest(contig=contig["contig_id"]):
                    self.assertGreater(start, previous_end)

    def test_contig_names_are_shared_with_the_source_genome(self) -> None:
        """VE-D15: el nombre no acredita identidad; solo la longitud distingue."""
        names = {c["contig_id"] for c in self.document["contigs"]}
        self.assertTrue(all(name.startswith("SL4.0ch") for name in names))
        self.assertTrue(all(c["length_bases"] < 100_000 for c in self.document["contigs"]))


class InconsistentGeometry(unittest.TestCase):
    def test_missing_index_is_an_error(self) -> None:
        with TemporaryDirectory() as directory:
            fasta = Path(directory) / "r.fa"
            fasta.write_text(">c1\nACGT\n", encoding="utf-8")
            coordinate_map = Path(directory) / "r.map.json"
            coordinate_map.write_text(json.dumps({"document_kind": "x"}), encoding="utf-8")
            with self.assertRaises(GeometryError):
                describe(fasta, coordinate_map)

    def test_wrong_map_kind_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            fasta = Path(directory) / "r.fa"
            fasta.write_text(">c1\nACGT\n", encoding="utf-8")
            (Path(directory) / "r.fa.fai").write_text("c1\t4\t4\t4\t5\n", encoding="utf-8")
            coordinate_map = Path(directory) / "r.map.json"
            coordinate_map.write_text(
                json.dumps({"document_kind": "otra_cosa"}), encoding="utf-8"
            )
            with self.assertRaises(GeometryError):
                describe(fasta, coordinate_map)

    def test_blocks_that_do_not_add_up_are_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            fasta = Path(directory) / "r.fa"
            fasta.write_text(">c1\nACGT\n", encoding="utf-8")
            (Path(directory) / "r.fa.fai").write_text("c1\t4\t4\t4\t5\n", encoding="utf-8")
            coordinate_map = Path(directory) / "r.map.json"
            coordinate_map.write_text(
                json.dumps(
                    {
                        "document_kind": "reference_reduction_map",
                        "derived_reference_id": "d",
                        "source_reference_id": "s",
                        "segments": [
                            {
                                "sequence_id": "c1",
                                "materialization_order": 1,
                                "orientation": "forward",
                                "sequence_class": "chromosome",
                                "source_interval": {"start": 1, "end": 99},
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaises(GeometryError):
                describe(fasta, coordinate_map)


if __name__ == "__main__":
    unittest.main()
