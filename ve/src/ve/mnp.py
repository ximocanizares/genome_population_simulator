"""Decompose multinucleotide records into the SNPs they represent (VE-D13).

FreeBayes calls by local haplotype. With its default `--haplotype-length 3`, two
true SNPs three bases apart or less are emitted as a single multinucleotide
record spanning the whole stretch, not as two biallelic records. That fusion is
a property of the caller facing real evidence — both variants travel on the same
reads — so VE-SNP does not separate the true loci to avoid it. The comparator
decomposes instead.

Decomposition is mechanical: for every ALT whose length equals REF's, walk the
bases in parallel and emit one derived SNP per differing position. Positions
where they agree are haplotype anchors and yield nothing. An ALT whose length
differs from REF is a complex or indel allele, outside the biological scope of
VE-SNP 0.1: it is not decomposed and is inventoried whole as a technical false
positive.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field

GT_SEPARATOR = re.compile(r"[/|]")


@dataclass(frozen=True)
class Call:
    """One VCF record as the comparator needs it."""

    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    genotypes: dict[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class DerivedSnp:
    """A single-base variant recovered from a record, fused or not."""

    chrom: str
    pos: int
    ref: str
    alts: tuple[str, ...]
    genotypes: dict[str, str]
    origin: str  # "direct" when the record was already a SNP, else "mnp"


@dataclass(frozen=True)
class Decomposition:
    """What a record yielded: derived SNPs plus the ALTs left whole."""

    snps: tuple[DerivedSnp, ...]
    undecomposed_alts: tuple[str, ...]

    @property
    def is_complex(self) -> bool:
        return bool(self.undecomposed_alts)


def _allele_indices(genotype: str) -> list[int | None]:
    """Split a GT into allele indices, with None for a no-call."""
    if not genotype or genotype in {".", "./.", ".|."}:
        return []
    indices: list[int | None] = []
    for field_ in GT_SEPARATOR.split(genotype.split(":")[0]):
        indices.append(None if field_ == "." else int(field_))
    return indices


def _base_of(allele: int | None, offset: int, ref: str, alts: tuple[str, ...]) -> str | None:
    """The base an allele contributes at one offset of the record."""
    if allele is None:
        return None
    if allele == 0:
        return ref[offset]
    sequence = alts[allele - 1]
    if len(sequence) != len(ref):
        return None
    return sequence[offset]


def decompose(call: Call) -> Decomposition:
    """Split one record into the SNPs it represents."""
    simple = [alt for alt in call.alts if len(alt) == len(call.ref)]
    complex_ = tuple(alt for alt in call.alts if len(alt) != len(call.ref))

    # A record that is already a single-base biallelic SNP passes through intact.
    if len(call.ref) == 1 and not complex_:
        if not simple:
            return Decomposition((), complex_)
        return Decomposition(
            (
                DerivedSnp(
                    call.chrom, call.pos, call.ref, call.alts, dict(call.genotypes), "direct"
                ),
            ),
            complex_,
        )

    derived: list[DerivedSnp] = []
    for offset in range(len(call.ref)):
        reference_base = call.ref[offset]
        # Alternative bases at this offset, keeping the order the ALTs had.
        bases: list[str] = []
        for alt in simple:
            base = alt[offset]
            if base != reference_base and base not in bases:
                bases.append(base)
        if not bases:
            continue  # anchor position: nothing changes here

        genotypes: dict[str, str] = {}
        for sample, genotype in call.genotypes.items():
            separator = "|" if "|" in genotype.split(":")[0] else "/"
            projected: list[str] = []
            for allele in _allele_indices(genotype):
                base = _base_of(allele, offset, call.ref, call.alts)
                if base is None:
                    projected.append(".")
                elif base == reference_base:
                    projected.append("0")
                else:
                    projected.append(str(bases.index(base) + 1))
            genotypes[sample] = separator.join(projected) if projected else "./."

        derived.append(
            DerivedSnp(
                call.chrom,
                call.pos + offset,
                reference_base,
                tuple(bases),
                genotypes,
                "mnp",
            )
        )

    return Decomposition(tuple(derived), complex_)
