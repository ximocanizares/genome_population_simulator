"""Build the run-B configuration from the artefacts run A published.

VE-D15 requires the B configuration to be generated from A's map instead of
transcribed by hand: the thirteen derived contig lengths are an output of A's
sampling, and copying them manually is both tedious and the easiest place to
introduce a silent mismatch. The simulator would reject a wrong length, but
there is no reason to rely on that.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

# VE-D04: the loci come from the existing F2.2 generator, with these keys.
LOCUS_COUNT = 1000
MIN_ALT_AF = 0.1
MAX_ALT_AF = 0.4
# Control neutro: las doce sustituciones equiprobables. El perfil humano de
# Harris introduciria un sesgo de espectro que el ensayo no valida y que ademas
# no corresponde a tomate.
SUBSTITUTION_PROFILE = "uniform_snv12_v1"

# VE-D03 and §4: cohort, ploidy and requested coverage.
INDIVIDUAL_COUNT = 20
PLOIDY = 2
TARGET_COVERAGE = 30
READ_LENGTH = 150

# A distinct i7/i5 pair per individual, derived deterministically from its index
# so the table is reproducible and never hand-written.
_BASES = "ACGT"


def _index_tag(value: int, salt: int) -> str:
    """Encode an index as eight A/C/G/T bases."""
    digits = []
    remainder = value + salt * 97 + 1
    for _ in range(8):
        digits.append(_BASES[remainder % 4])
        remainder //= 4
    return "".join(digits)


def _derived_sequences(geometry: dict[str, object]) -> list[dict[str, object]]:
    sequences = []
    for contig in geometry["contigs"]:  # type: ignore[index]
        classes = {b["sequence_class"] for b in contig["source_intervals"]}
        if len(classes) != 1:
            raise ValueError(f"{contig['contig_id']}: clases de secuencia mezcladas")
        sequences.append(
            {
                "id": contig["contig_id"],
                "class": classes.pop(),
                "length": contig["length_bases"],
            }
        )
    return sequences


def build(
    geometry_path: Path,
    fasta_path: Path,
    seed: str,
    individual_count: int = INDIVIDUAL_COUNT,
    target_coverage: float = TARGET_COVERAGE,
    locus_count: int = LOCUS_COUNT,
) -> str:
    """Render the run-B TOML.

    The cohort, coverage and locus count are parameters so a reduced variant can
    be built for scaling measurements without hand-editing the configuration.
    """
    geometry = json.loads(geometry_path.read_text(encoding="utf-8"))
    if geometry.get("document_kind") != "ve_realised_reference_geometry":
        raise ValueError("la geometria no es la publicada por VE.1.2")
    sequences = _derived_sequences(geometry)
    total = sum(int(s["length"]) for s in sequences)
    if total != geometry["total_bases"]:
        raise ValueError("las longitudes no suman el total publicado")

    individuals = [f"ind-{index + 1:02d}" for index in range(individual_count)]
    lines: list[str] = []
    add = lines.append

    add("# Ejecucion B de VE-SNP 0.1: genera SNP, verdad, VCF y FASTQ contra la")
    add("# referencia derivada por la ejecucion A.")
    add("#")
    add("# GENERADO por ve.assay_config desde la geometria publicada en VE.1.2.")
    add("# No editar a mano: VE-D15 exige derivar las longitudes del mapa de A.")
    add("")
    add('config_kind = "general"')
    add('schema_version = "2"')
    add(f'seed = "{seed}"')
    add(f"ploidy = {PLOIDY}")
    add("")
    add("[reference]")
    add('kind = "supplied"')
    add(f'id = "{geometry["derived_reference_id"]}"')
    add(f'fasta = "{fasta_path}"')
    for sequence in sequences:
        add("")
        add("[[reference.sequences]]")
        add(f'id = "{sequence["id"]}"')
        add(f'class = "{sequence["class"]}"')
        add(f"length = {sequence['length']}")

    add("")
    add("# VE-D04: mil candidatos A/C/G/T, frecuencia ALT uniforme en [0,1; 0,4).")
    add("# Sin separacion minima entre loci: la fusion MNP se resuelve en el")
    add("# comparador por VE-D13, no alterando la verdad.")
    add("[population_snvs]")
    add('mode = "generate"')
    add("")
    add("[population_snvs.generate]")
    add('count_mode = "exact"')
    add(f"locus_count = {locus_count}")
    add('spatial_mode = "uniform"')
    add(f'substitution_profile = "{SUBSTITUTION_PROFILE}"')
    add("")
    add("[[population_snvs.generate.populations]]")
    add('population_id = "pop-ve"')
    add('frequency_model = "uniform"')
    add(f"min_alt_af = {MIN_ALT_AF}")
    add(f"max_alt_af = {MAX_ALT_AF}")

    add("")
    add("[[initial_populations]]")
    add('population_id = "pop-ve"')
    add(f"individual_count = {individual_count}")
    add(f"ploidy = {PLOIDY}")
    add('genotype_model = "hwe"')
    add("individual_ids = [")
    for name in individuals:
        add(f'  "{name}",')
    add("]")

    add("")
    add("# VE-D02: PE neutro. sequencing_errors disabled produce cero")
    add("# sustituciones; quimeras, index hopping y contaminacion quedan a cero.")
    add("[technical_event]")
    add('event_id = "ve-snp-pe-neutral"')
    add(f"target_coverage_wgs = {target_coverage}")
    add(f"read_length_bases = {READ_LENGTH}")
    add('sequencing_errors = "disabled"')
    add("chimera_rate = 0")
    add("index_hopping_rate = 0")
    add("contamination_min_percent = 0")
    add("contamination_max_percent = 0")
    add("individual_ids = [")
    for name in individuals:
        add(f'  "{name}",')
    add("]")
    add("")
    add("# Una sola biblioteca por individuo (§7.1): F4 publica un unico par de")
    add("# FASTQ por individuo, de modo que dos bibliotecas romperian la biyeccion.")
    for index, name in enumerate(individuals):
        add("")
        add("[[technical_event.libraries]]")
        add(f'individual_id = "{name}"')
        add(f'i7 = "{_index_tag(index, 0)}"')
        add(f'i5 = "{_index_tag(index, 1)}"')

    add("")
    add("# VCF de verdad en texto: VE-SNP conserva el VCF sin comprimir, de modo")
    add("# que el ensayo no necesita bgzip ni tabix (VE.1.1).")
    add("[vcf_output]")
    add('mode = "text"')
    add("")
    return "\n".join(lines)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="ve.assay_config")
    parser.add_argument("geometry", type=Path)
    parser.add_argument("fasta", type=Path)
    parser.add_argument("--seed", required=True)
    parser.add_argument("--individuals", type=int, default=INDIVIDUAL_COUNT)
    parser.add_argument("--coverage", type=float, default=TARGET_COVERAGE)
    parser.add_argument("--loci", type=int, default=LOCUS_COUNT)
    parser.add_argument("--out", type=Path)
    arguments = parser.parse_args(argv)
    text = build(
        arguments.geometry,
        arguments.fasta.resolve(),
        arguments.seed,
        arguments.individuals,
        arguments.coverage,
        arguments.loci,
    )
    if arguments.out:
        arguments.out.write_text(text, encoding="utf-8")
        print(f"escrito {arguments.out}")
    else:
        print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
