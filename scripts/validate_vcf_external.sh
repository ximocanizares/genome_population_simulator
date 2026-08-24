#!/usr/bin/env bash
# Puerta DTC-066: aceptación externa, nunca parte del runtime del simulador.
set -euo pipefail

if [[ $# -ne 5 ]]; then
  echo "uso: $0 <reference.fa> <text.vcf> <bgzf-tbi.vcf.gz> <bgzf-csi.vcf.gz> <region>" >&2
  exit 64
fi

reference=$1
text_vcf=$2
tbi_vcf=$3
csi_vcf=$4
region=$5
temporary=$(mktemp -d "${TMPDIR:-/tmp}/gps-dtc066.XXXXXX")

fail() {
  echo "DTC-066 no superada: $*" >&2
  echo "Los temporales de validación se conservan en: $temporary" >&2
  exit 1
}

cleanup() {
  rm -rf "$temporary"
}

require_version() {
  local program=$1 expected=$2 observed
  observed=$($program --version | sed -n '1p')
  [[ $observed == "$expected" ]] || fail "se requiere $expected; se observó ${observed:-ausente}."
}

require_sidecars() {
  local vcf=$1 expected=$2 opposite
  if [[ $expected == tbi ]]; then opposite=csi; else opposite=tbi; fi
  [[ -f "${vcf}.${expected}" ]] || fail "falta el sidecar exacto ${vcf}.${expected}."
  [[ ! -e "${vcf}.${opposite}" ]] || fail "aparece el índice opuesto ${vcf}.${opposite}."
}

require_no_sidecars() {
  local vcf=$1
  [[ ! -e "${vcf}.tbi" && ! -e "${vcf}.csi" ]] || fail "el VCF de texto no puede llevar índice."
}

validate_text() {
  local vcf=$1 label=$2
  bcftools view --no-version "$vcf" >"$temporary/${label}.view.vcf"
  bcftools norm --no-version --check-ref e --fasta-ref "$reference" "$vcf" >"$temporary/${label}.norm.vcf"
  cmp -- "$temporary/${label}.view.vcf" "$temporary/${label}.norm.vcf" || fail "la validación normalizó o alteró $vcf; no se permite reparar ni descartar salida."
}

validate_indexed() {
  local vcf=$1 kind=$2 label=$3 copied
  require_sidecars "$vcf" "$kind"
  validate_text "$vcf" "$label"
  bcftools index --stats "$vcf" >"$temporary/${label}.original.stats"
  bcftools view --no-version --regions "$region" "$vcf" >"$temporary/${label}.original.region.vcf"
  tabix -l "$vcf" >"$temporary/${label}.original.contigs"
  tabix "$vcf" "$region" >"$temporary/${label}.original.region.tabix"

  copied="$temporary/${label}.vcf.gz"
  cp -- "$vcf" "$copied"
  if [[ $kind == tbi ]]; then bcftools index --force --tbi "$copied"; else bcftools index --force --csi "$copied"; fi
  require_sidecars "$copied" "$kind"
  bcftools index --stats "$copied" >"$temporary/${label}.reconstructed.stats"
  bcftools view --no-version --regions "$region" "$copied" >"$temporary/${label}.reconstructed.region.vcf"
  tabix -l "$copied" >"$temporary/${label}.reconstructed.contigs"
  tabix "$copied" "$region" >"$temporary/${label}.reconstructed.region.tabix"
  cmp -- "$temporary/${label}.original.stats" "$temporary/${label}.reconstructed.stats" || fail "la reconstrucción temporal del índice $kind no reproduce sus estadísticas."
  cmp -- "$temporary/${label}.original.region.vcf" "$temporary/${label}.reconstructed.region.vcf" || fail "la consulta regional tras reconstruir el índice $kind difiere."
  cmp -- "$temporary/${label}.original.contigs" "$temporary/${label}.reconstructed.contigs" || fail "la lista HTSlib tras reconstruir el índice $kind difiere."
  cmp -- "$temporary/${label}.original.region.tabix" "$temporary/${label}.reconstructed.region.tabix" || fail "la consulta regional HTSlib tras reconstruir el índice $kind difiere."
}

expect_rejection() {
  local label=$1 diagnostic
  shift
  if diagnostic=$( ( "$@" ) 2>&1); then
    fail "la mutación $label fue aceptada."
  fi
  echo "DTC-066: rechazo esperado verificado para $label." >&2
}

exercise_mutation_helpers() {
  local fixture_ref="$temporary/mutation.fa" fixture_text="$temporary/mutation.vcf"
  local fixture_tbi="$temporary/mutation-tbi.vcf.gz" fixture_csi="$temporary/mutation-csi.vcf.gz"
  local mutated_ref="$temporary/mutated-ref.vcf" opposite="$temporary/opposite.vcf.gz" both="$temporary/both.vcf.gz"
  printf '>mutation\nA\n' >"$fixture_ref"
  printf '##fileformat=VCFv4.5\n##contig=<ID=mutation,length=1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nmutation\t1\t.\tA\tC\t.\t.\t.\n' >"$fixture_text"
  bcftools view --no-version --output-type z --output-file "$fixture_tbi" "$fixture_text"
  bcftools index --tbi "$fixture_tbi"
  bcftools view --no-version --output-type z --output-file "$fixture_csi" "$fixture_text"
  bcftools index --csi "$fixture_csi"

  sed 's/\tA\tC\t/\tG\tC\t/' "$fixture_text" >"$mutated_ref"
  if bcftools norm --no-version --check-ref e --fasta-ref "$fixture_ref" "$mutated_ref" >/dev/null 2>&1; then fail "la mutación de REF fue aceptada."; fi

  cp -- "$fixture_tbi" "$opposite"
  cp -- "${fixture_tbi}.tbi" "${opposite}.csi"
  expect_rejection índice_opuesto require_sidecars "$opposite" tbi

  cp -- "$fixture_tbi" "$both"
  cp -- "${fixture_tbi}.tbi" "${both}.tbi"
  cp -- "${fixture_csi}.csi" "${both}.csi"
  expect_rejection ambos_índices require_sidecars "$both" tbi
}

require_version bcftools 'bcftools 1.24'
require_version tabix 'tabix (htslib) 1.24'
require_no_sidecars "$text_vcf"
validate_text "$text_vcf" text
validate_indexed "$tbi_vcf" tbi tbi
validate_indexed "$csi_vcf" csi csi
exercise_mutation_helpers
cleanup
