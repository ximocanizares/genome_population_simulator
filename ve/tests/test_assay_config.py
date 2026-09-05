"""VE.1.3: la configuracion de B se genera desde el mapa de A, no se transcribe."""

from __future__ import annotations

import json
import tomllib
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

from ve.assay_config import (
    INDIVIDUAL_COUNT,
    LOCUS_COUNT,
    MAX_ALT_AF,
    MIN_ALT_AF,
    PLOIDY,
    SUBSTITUTION_PROFILE,
    TARGET_COVERAGE,
    build,
)

ASSAY = Path(__file__).resolve().parents[1] / "assays" / "ve-snp-0.1"
GEOMETRY = ASSAY / "reference-geometry.json"
SEED = "2986c45543aaa394cd419b1a731dbe7dbbc3523943f7215515eb6f400bca596b"


class GeneratedConfiguration(unittest.TestCase):
    def setUp(self) -> None:
        if not GEOMETRY.is_file():
            self.skipTest("la geometria de VE.1.2 aun no existe")
        self.document = tomllib.loads(build(GEOMETRY, Path("/tmp/r.fa"), SEED))
        self.geometry = json.loads(GEOMETRY.read_text(encoding="utf-8"))

    def test_sequences_come_from_the_published_geometry(self) -> None:
        """Las trece longitudes son salida del muestreo de A, no valores a mano."""
        expected = {c["contig_id"]: c["length_bases"] for c in self.geometry["contigs"]}
        observed = {s["id"]: s["length"] for s in self.document["reference"]["sequences"]}
        self.assertEqual(observed, expected)

    def test_total_is_the_derived_reference_size(self) -> None:
        total = sum(s["length"] for s in self.document["reference"]["sequences"])
        self.assertEqual(total, 500_000)

    def test_seed_is_declared(self) -> None:
        """§8.1: sin semilla declarada la repeticion no significaria nada."""
        self.assertEqual(self.document["seed"], SEED)

    def test_locus_generation_matches_ve_d04(self) -> None:
        generate = self.document["population_snvs"]["generate"]
        self.assertEqual(self.document["population_snvs"]["mode"], "generate")
        self.assertEqual(generate["count_mode"], "exact")
        self.assertEqual(generate["locus_count"], LOCUS_COUNT)
        self.assertEqual(generate["spatial_mode"], "uniform")
        self.assertEqual(generate["substitution_profile"], SUBSTITUTION_PROFILE)
        population = generate["populations"][0]
        self.assertEqual(population["frequency_model"], "uniform")
        self.assertEqual(population["min_alt_af"], MIN_ALT_AF)
        self.assertEqual(population["max_alt_af"], MAX_ALT_AF)

    def test_cohort_matches_the_assay(self) -> None:
        population = self.document["initial_populations"][0]
        self.assertEqual(population["individual_count"], INDIVIDUAL_COUNT)
        self.assertEqual(population["ploidy"], PLOIDY)
        self.assertEqual(len(population["individual_ids"]), INDIVIDUAL_COUNT)
        self.assertEqual(len(set(population["individual_ids"])), INDIVIDUAL_COUNT)

    def test_technical_event_is_neutral(self) -> None:
        """VE-D02: PE neutro, con los demas artefactos tecnicos a cero."""
        event = self.document["technical_event"]
        self.assertEqual(event["sequencing_errors"], "disabled")
        self.assertEqual(event["target_coverage_wgs"], TARGET_COVERAGE)
        self.assertEqual(event["chimera_rate"], 0)
        self.assertEqual(event["index_hopping_rate"], 0)
        self.assertEqual(event["contamination_min_percent"], 0)
        self.assertEqual(event["contamination_max_percent"], 0)

    def test_exactly_one_library_per_individual(self) -> None:
        """§7.1: dos bibliotecas romperian la biyeccion individuo->FASTQ->RG->BAM."""
        libraries = self.document["technical_event"]["libraries"]
        individuals = self.document["initial_populations"][0]["individual_ids"]
        self.assertEqual(len(libraries), len(individuals))
        self.assertEqual({lib["individual_id"] for lib in libraries}, set(individuals))
        self.assertEqual(len({lib["i7"] for lib in libraries}), len(libraries))
        self.assertEqual(len({lib["i5"] for lib in libraries}), len(libraries))
        for library in libraries:
            self.assertRegex(library["i7"], r"^[ACGT]{8}$")
            self.assertRegex(library["i5"], r"^[ACGT]{8}$")

    def test_vcf_is_plain_text(self) -> None:
        """El ensayo no comprime, de modo que no necesita bgzip ni tabix."""
        self.assertEqual(self.document["vcf_output"]["mode"], "text")

    def test_generation_is_deterministic(self) -> None:
        again = build(GEOMETRY, Path("/tmp/r.fa"), SEED)
        self.assertEqual(again, build(GEOMETRY, Path("/tmp/r.fa"), SEED))


class RejectedGeometry(unittest.TestCase):
    def test_document_of_another_kind_is_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "g.json"
            path.write_text(json.dumps({"document_kind": "otra"}), encoding="utf-8")
            with self.assertRaises(ValueError):
                build(path, Path("/tmp/r.fa"), SEED)

    def test_lengths_that_do_not_add_up_are_rejected(self) -> None:
        with TemporaryDirectory() as directory:
            path = Path(directory) / "g.json"
            path.write_text(
                json.dumps(
                    {
                        "document_kind": "ve_realised_reference_geometry",
                        "derived_reference_id": "d",
                        "source_reference_id": "s",
                        "total_bases": 999,
                        "contigs": [
                            {
                                "contig_id": "c1",
                                "length_bases": 10,
                                "block_count": 1,
                                "source_intervals": [
                                    {
                                        "source_start": 1,
                                        "source_end": 10,
                                        "orientation": "forward",
                                        "sequence_class": "chromosome",
                                    }
                                ],
                            }
                        ],
                    }
                ),
                encoding="utf-8",
            )
            with self.assertRaises(ValueError):
                build(path, Path("/tmp/r.fa"), SEED)


if __name__ == "__main__":
    unittest.main()
