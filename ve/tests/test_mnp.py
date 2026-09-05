"""VE-D13: los nueve fixtures deterministas del descomponedor (VE.1.6).

El ensayo produce unas seis fusiones por azar, que no bastan para validar el
modulo. Estos casos fijan de antemano el resultado esperado de cada situacion.
"""

from __future__ import annotations

import unittest

from ve.mnp import Call, decompose


class Cis(unittest.TestCase):
    """E5: dos heterocigotos en la misma copia -> MNP bialelico."""

    def test_two_snps_two_bases_apart_yield_two_snps(self) -> None:
        # Es la fusion medida sobre FreeBayes v1.3.10: AGC -> CGA.
        result = decompose(Call("c1", 500, "AGC", ("CGA",), {"s1": "0/1"}))
        self.assertFalse(result.is_complex)
        self.assertEqual(len(result.snps), 2)
        first, second = result.snps
        self.assertEqual((first.pos, first.ref, first.alts), (500, "A", ("C",)))
        self.assertEqual((second.pos, second.ref, second.alts), (502, "C", ("A",)))
        self.assertEqual(first.genotypes["s1"], "0/1")
        self.assertEqual(second.genotypes["s1"], "0/1")
        self.assertTrue(all(snp.origin == "mnp" for snp in result.snps))

    def test_anchor_base_yields_nothing(self) -> None:
        """La base intermedia no cambia: es anclaje del haplotipo, no variante."""
        result = decompose(Call("c1", 500, "AGC", ("CGA",), {"s1": "0/1"}))
        self.assertNotIn(501, [snp.pos for snp in result.snps])


class Trans(unittest.TestCase):
    """E6: dos heterocigotos en copias distintas -> MNP multialelico 1/2."""

    def test_multiallelic_record_yields_two_heterozygous_snps(self) -> None:
        result = decompose(Call("c1", 600, "AGC", ("CGC", "AGA"), {"s1": "1/2"}))
        self.assertFalse(result.is_complex)
        self.assertEqual(len(result.snps), 2)
        first, second = result.snps
        self.assertEqual((first.pos, first.ref, first.alts), (600, "A", ("C",)))
        self.assertEqual((second.pos, second.ref, second.alts), (602, "C", ("A",)))
        # Cada alelo aporta su base en su posicion y la de referencia en la otra.
        self.assertEqual(first.genotypes["s1"], "1/0")
        self.assertEqual(second.genotypes["s1"], "0/1")

    def test_phase_separator_is_preserved(self) -> None:
        result = decompose(Call("c1", 600, "AGC", ("CGC", "AGA"), {"s1": "1|2"}))
        self.assertIn("|", result.snps[0].genotypes["s1"])


class Homozygous(unittest.TestCase):
    def test_homozygous_alternate_projects_as_homozygous(self) -> None:
        result = decompose(Call("c1", 700, "AG", ("CG",), {"s1": "1/1"}))
        self.assertEqual(len(result.snps), 1)
        self.assertEqual(result.snps[0].genotypes["s1"], "1/1")
        self.assertEqual(result.snps[0].pos, 700)


class LongerStretch(unittest.TestCase):
    def test_four_base_record_with_three_changes(self) -> None:
        result = decompose(Call("c1", 800, "ACGT", ("TCAA",), {"s1": "0/1"}))
        self.assertEqual([snp.pos for snp in result.snps], [800, 802, 803])
        self.assertEqual([snp.ref for snp in result.snps], ["A", "G", "T"])
        self.assertEqual([snp.alts[0] for snp in result.snps], ["T", "A", "A"])


class ComplexAllele(unittest.TestCase):
    """E7: len(REF) != len(ALT) no se descompone; es FP tecnico integro."""

    def test_deletion_is_left_whole(self) -> None:
        result = decompose(Call("c1", 900, "AT", ("A",), {"s1": "0/1"}))
        self.assertTrue(result.is_complex)
        self.assertEqual(result.undecomposed_alts, ("A",))
        self.assertEqual(result.snps, ())

    def test_mixed_record_decomposes_only_the_equal_length_allele(self) -> None:
        result = decompose(Call("c1", 950, "AG", ("CG", "A"), {"s1": "1/2"}))
        self.assertEqual(result.undecomposed_alts, ("A",))
        self.assertEqual(len(result.snps), 1)
        self.assertEqual(result.snps[0].pos, 950)
        # El alelo 2 no aporta base comparable: se marca ausente, no se inventa.
        self.assertEqual(result.snps[0].genotypes["s1"], "1/.")


class PlainSnp(unittest.TestCase):
    """E1: un SNP simple atraviesa el modulo sin alterarse."""

    def test_single_base_record_passes_through(self) -> None:
        result = decompose(Call("c1", 100, "A", ("G",), {"s1": "0|1"}))
        self.assertFalse(result.is_complex)
        self.assertEqual(len(result.snps), 1)
        snp = result.snps[0]
        self.assertEqual((snp.pos, snp.ref, snp.alts), (100, "A", ("G",)))
        self.assertEqual(snp.genotypes["s1"], "0|1")
        self.assertEqual(snp.origin, "direct")


class NoCall(unittest.TestCase):
    def test_missing_genotype_stays_missing(self) -> None:
        result = decompose(Call("c1", 500, "AGC", ("CGA",), {"s1": "./."}))
        self.assertEqual(len(result.snps), 2)
        for snp in result.snps:
            self.assertEqual(snp.genotypes["s1"], "./.")


class MultipleSamples(unittest.TestCase):
    def test_every_sample_is_projected_independently(self) -> None:
        result = decompose(
            Call("c1", 500, "AGC", ("CGA",), {"s1": "0/1", "s2": "1/1", "s3": "0/0"})
        )
        first = result.snps[0]
        self.assertEqual(first.genotypes["s1"], "0/1")
        self.assertEqual(first.genotypes["s2"], "1/1")
        self.assertEqual(first.genotypes["s3"], "0/0")


if __name__ == "__main__":
    unittest.main()
