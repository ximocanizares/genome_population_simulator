#!/usr/bin/env bash
# Puerta DTC-066: aceptación externa, nunca parte del runtime del simulador.
#
# Cada fixture es una tupla independiente:
#   <reference.fa> <input.vcf[.gz]> <none|tbi|csi> <region|-> <expected-records>
#
# `expected-records` contiene exactamente las líneas de registros (sin cabecera)
# que deben devolver bcftools/tabix. Se pasan, como mínimo, los tres positivos:
# texto/none, FASTA corto/TBI y FASTA largo/CSI. El tipo esperado no se toma de
# la cabecera VCF: se vuelve a derivar de las longitudes del FASTA autoritativo.
set -euo pipefail

readonly TBI_MAXIMUM_POSITION_EXCLUSIVE=536870912

if (( $# == 0 || $# % 5 != 0 )); then
  echo "uso: $0 (<reference.fa> <input.vcf[.gz]> <none|tbi|csi> <region|-> <expected-records>)+" >&2
  exit 64
fi

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

# La selección DTC-063 depende solamente de la longitud real de cada secuencia
# FASTA. No se acepta que un ##contig VCF ni su índice decidan este resultado.
index_kind_from_reference() {
  local reference=$1
  awk -v threshold="$TBI_MAXIMUM_POSITION_EXCLUSIVE" '
    function finish_sequence() {
      if (seen && length_bases >= threshold) exceeds = 1
    }
    /^>/ {
      finish_sequence()
      seen = 1
      length_bases = 0
      next
    }
    {
      sequence = $0
      gsub(/[[:space:]]/, "", sequence)
      length_bases += length(sequence)
    }
    END {
      finish_sequence()
      if (!seen) exit 2
      print exceeds ? "csi" : "tbi"
    }
  ' "$reference" || fail "no se pudieron derivar longitudes de $reference."
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

compare_exact_records() {
  local observed=$1 expected=$2 description=$3
  cmp -- "$expected" "$observed" || fail "$description no coincide exactamente con $expected."
}

# Los registros de repetición en tándem quedan FUERA de esta comparación, y no
# por comodidad: la decisión S-09 manda publicarlos con la secuencia completa en
# REF y ALT, mientras que DTC-064 recorta todo lo demás. `bcftools norm` recorta
# también los de repetición, de modo que exigir que no altere el fichero los
# haría fallar siempre por una diferencia de codificación buscada a propósito.
#
# No se apaga la comprobación: lo que a ellos les corresponde lo verifica una
# prueba de Rust —`the_external_normalizer_only_trims_the_repeat_record_without_
# changing_its_alleles`— que ejecuta el mismo `bcftools norm` y exige que la
# forma normalizada conserve la longitud de cada alelo respecto de su REF, es
# decir, que el juez externo solo recorte y no cambie ningún alelo.
records_without_tandem_repeats() {
  grep -v $'\tRN=\|;RN=' -- "$1" || true
}

validate_complete_view_and_normalization() {
  local reference=$1 vcf=$2 label=$3
  local parsed="$temporary/${label}.parsed.vcf" normalized="$temporary/${label}.normalized.vcf"
  bcftools view --no-version --output-type v --output "$parsed" "$vcf"
  bcftools norm --no-version --check-ref e --fasta-ref "$reference" \
    --output-type v --output "$normalized" "$vcf"
  records_without_tandem_repeats "$parsed" >"$temporary/${label}.parsed.plain"
  records_without_tandem_repeats "$normalized" >"$temporary/${label}.normalized.plain"
  # Sin esta comprobación un fichero compuesto solo de registros de repetición
  # dejaría la comparación reducida a las cabeceras, y la puerta pasaría en
  # silencio sin haber comparado ningún registro.
  local remaining
  remaining=$(grep -c -v '^#' -- "$temporary/${label}.parsed.plain" || true)
  [[ ${remaining:-0} -gt 0 ]] \
    || fail "tras excluir los registros de repetición no queda ningún registro que comparar en $vcf; la comparación con el normalizador externo sería vacua."
  cmp -- "$temporary/${label}.parsed.plain" "$temporary/${label}.normalized.plain" \
    || fail "la normalización externa alteró un registro no repetitivo de $vcf; no se permite reparar ni descartar salida."
}

validate_fixture() {
  local reference=$1 vcf=$2 expected_kind=$3 region=$4 expected_records=$5 label=$6
  local derived_kind copied
  derived_kind=$(index_kind_from_reference "$reference")
  [[ $expected_kind == none || $expected_kind == tbi || $expected_kind == csi ]] \
    || fail "$label declara un tipo de índice desconocido: $expected_kind."
  [[ -f $reference && -f $vcf && -f $expected_records ]] \
    || fail "$label no contiene todos sus artefactos de fixture."

  if [[ $expected_kind == none ]]; then
    [[ $region == - ]] || fail "$label de texto debe declarar región '-'."
    [[ $vcf != *.gz ]] || fail "$label de texto debe ser un .vcf sin compresión BGZF."
    require_no_sidecars "$vcf"
  else
    [[ $derived_kind == "$expected_kind" ]] || fail "$label declara $expected_kind, pero las longitudes de $reference requieren $derived_kind."
    [[ $region != - ]] || fail "$label indexado exige una región de consulta."
    [[ $vcf == *.vcf.gz ]] || fail "$label indexado debe ser un .vcf.gz BGZF."
    require_sidecars "$vcf" "$expected_kind"
  fi

  validate_complete_view_and_normalization "$reference" "$vcf" "$label"
  if [[ $expected_kind == none ]]; then
    bcftools view --no-version --no-header --output-type v \
      --output "$temporary/${label}.full.records" "$vcf"
    compare_exact_records "$temporary/${label}.full.records" "$expected_records" \
      "la lectura completa de bcftools para $label"
    return
  fi

  bcftools index --stats "$vcf" >"$temporary/${label}.original.stats"
  bcftools view --no-version --no-header --regions "$region" --output-type v \
    --output "$temporary/${label}.bcftools.records" "$vcf"
  tabix -l "$vcf" >"$temporary/${label}.contigs"
  tabix "$vcf" "$region" >"$temporary/${label}.tabix.records"
  compare_exact_records "$temporary/${label}.bcftools.records" "$expected_records" \
    "la consulta regional de bcftools para $label"
  compare_exact_records "$temporary/${label}.tabix.records" "$expected_records" \
    "la consulta regional de tabix para $label"

  # Se conserva una comprobación secundaria de reconstrucción, pero jamás es
  # la única evidencia: todas las salidas anteriores se compararon con la
  # verdad independiente del fixture.
  copied="$temporary/${label}.vcf.gz"
  cp -- "$vcf" "$copied"
  if [[ $expected_kind == tbi ]]; then
    bcftools index --tbi --output "${copied}.tbi" "$copied"
  else
    bcftools index --csi --output "${copied}.csi" "$copied"
  fi
  require_sidecars "$copied" "$expected_kind"
  bcftools index --stats "$copied" >"$temporary/${label}.reconstructed.stats"
  cmp -- "$temporary/${label}.original.stats" "$temporary/${label}.reconstructed.stats" \
    || fail "la reconstrucción temporal del índice $expected_kind no reproduce sus estadísticas para $label."
}

expect_rejection() {
  local label=$1 diagnostic
  shift
  if diagnostic=$( ( "$@" ) 2>&1); then
    fail "la mutación $label fue aceptada."
  fi
  echo "DTC-066: rechazo esperado verificado para $label." >&2
}

require_version bcftools 'bcftools 1.24'
require_version tabix 'tabix (htslib) 1.24'

saw_text=0
saw_tbi=0
saw_csi=0
short_reference=
short_vcf=
short_region=
short_expected_records=
fixture_count=0

while (( $# > 0 )); do
  reference=$1
  vcf=$2
  expected_kind=$3
  region=$4
  expected_records=$5
  shift 5
  label="fixture-${fixture_count}-${expected_kind}"
  validate_fixture "$reference" "$vcf" "$expected_kind" "$region" "$expected_records" "$label"
  case $expected_kind in
    none) saw_text=1 ;;
    tbi)
      saw_tbi=1
      short_reference=$reference
      short_vcf=$vcf
      short_region=$region
      short_expected_records=$expected_records
      ;;
    csi) saw_csi=1 ;;
  esac
  fixture_count=$((fixture_count + 1))
done

(( saw_text == 1 && saw_tbi == 1 && saw_csi == 1 )) \
  || fail "el corpus debe incluir positivos independientes texto/none, FASTA corto/TBI y FASTA largo/CSI."

# Mutación negativa: el mismo FASTA corto no puede cambiar a CSI solo porque
# el fixture lo declare. Esto prueba que la elección procede de longitudes, no
# del nombre del fichero, la cabecera VCF ni el sidecar existente.
expect_rejection referencia_corta_declara_csi \
  validate_fixture "$short_reference" "$short_vcf" csi "$short_region" "$short_expected_records" \
    mutation-short-csi

cleanup
