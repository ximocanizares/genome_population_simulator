#!/usr/bin/env bash
# Puerta del barrido de guardas de los jueces independientes, de los dos
# verificadores de configuración efectiva y manifiesto de evento y, desde el
# 05/10/2026, de los nueve módulos que juzgan lo publicado o lo leído (verdad,
# FASTQ, matriz de SNV, paquete F5, FASTA, mapa de reducción, identidad de
# artefactos y perfil técnico). La lista autoritativa es el array `judges`.
#
# QUÉ MIDE, Y QUÉ NO. Mide **si una guarda llega a dispararse** durante la
# batería. No mide si la prueba que la dispara **fallaría** en caso de que la
# guarda estuviera mal escrita. Son cosas distintas y la segunda es la que
# importa: una prueba que provoca un rechazo sin exigir **cuál** rechazo es
# satisface esta puerta igual que una que exige el motivo exacto. Por tanto
# esta puerta es un **suelo**, no un techo: impide que la cobertura vuelva a
# caer sin que nadie se entere, y nada más. Lo que demuestra que una guarda
# está probada sigue siendo mutarla y ver morir la prueba, a mano, una mutación
# cada vez.
#
# Uso: scripts/check-guard-coverage.sh DIRECTORIO_DE_EVIDENCIA
#
# El directorio no puede existir: la evidencia de cada ejecución es suya.
set -euo pipefail

project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)

if (( $# != 1 )); then
    echo "usage: $0 EVIDENCE_DIRECTORY" >&2
    exit 2
fi
evidence_dir=$1
if [[ -e "$evidence_dir" ]]; then
    echo "evidence directory already exists: $evidence_dir" >&2
    exit 2
fi
mkdir -p -- "$evidence_dir/trace"
evidence_dir=$(cd -- "$evidence_dir" && pwd)

baseline="$project_root/scripts/guard-coverage-baseline.tsv"

# Cada módulo vigilado con las constructoras de error por las que canaliza sus rechazos.
judges=(
    "src/f5_oracle.rs|reject|resource_reject|catalog_reject|truth_reject"
    "src/vcf_semantic_validation.rs|vcf_reject|oracle_error|vcf_error|fasta_reject|io_reject|io_error|bases_not_materialized"
    "src/f5_closure.rs|reject"
    "src/results.rs|effective_error|effective_error_for_schema|manifest_error|manifest_error_for_schema|repeat_section_error|reason|physical_error|artifact_identity_error"
    "src/event_manifest.rs|event_error"
    "src/truth_validation.rs|rejected|mismatch"
    "src/truth_schema.rs|rejected"
    "src/fastq_validation.rs|rejected|missing_canonical_eof"
    "src/population_snvs.rs|rejected|audit_error|missing|invalid"
    "src/f5_publication.rs|rejected|document_error|manifest_error|failure_mismatch"
    "src/reference_fasta.rs|rejected|empty_header|empty_fasta|invalid_symbol|configuration_mismatch|reference_error"
    "src/reference_reduction_io.rs|map_error|out_of_range"
    "src/integrity.rs|rejected"
    "src/technical_profile.rs|rejected|invalid"
)

echo "== batería con el barrido encendido ==" >&2
set +e
env GPS_GUARD_TRACE_DIR="$evidence_dir/trace" \
    cargo test --manifest-path "$project_root/Cargo.toml" \
        --workspace --all-targets --locked \
    > "$evidence_dir/cargo-test.log" 2>&1
test_status=$?
set -e
{
    grep -E '^test result:' "$evidence_dir/cargo-test.log" \
        | awk '{p+=$4; f+=$6; i+=$8} END {printf "passed=%d failed=%d ignored=%d\n", p, f, i}'
    echo "exit=$test_status"
} > "$evidence_dir/battery-summary.txt"
cat "$evidence_dir/battery-summary.txt" >&2
if (( test_status != 0 )); then
    echo "la batería no pasó: el barrido no es interpretable" >&2
    exit 1
fi

# Las guardas disparadas, sin repetir, tal como las anotó cada proceso.
cat "$evidence_dir"/trace/*.guards 2>/dev/null | sort -u > "$evidence_dir/fired.txt"

printf '| Juez | Comprobaciones | Nunca se disparan |\n|---|---|---|\n' \
    > "$evidence_dir/tabla.md"

failures=0
for judge in "${judges[@]}"; do
    module=${judge%%|*}
    names=${judge#*|}
    awk -v names="$names" -f "$project_root/scripts/guard-sites.awk" \
        "$project_root/$module" \
        | sed "s#^#$module:#" | sort -u > "$evidence_dir/sites-$(basename "$module").txt"
    comm -23 "$evidence_dir/sites-$(basename "$module").txt" "$evidence_dir/fired.txt" \
        > "$evidence_dir/never-fired-$(basename "$module").txt"
    total=$(wc -l < "$evidence_dir/sites-$(basename "$module").txt")
    dead=$(wc -l < "$evidence_dir/never-fired-$(basename "$module").txt")
    ceiling=$(awk -F'\t' -v m="$module" '$1 == m {print $3}' "$baseline")
    declared=$(awk -F'\t' -v m="$module" '$1 == m {print $2}' "$baseline")
    if [[ -z "$ceiling" || -z "$declared" ]]; then
        echo "$module no tiene fila en $baseline" >&2
        exit 2
    fi
    # Sin esto, borrar una guarda baja las dos cifras a la vez y la puerta
    # pasaría: no distinguiría **cubrir** una comprobación de **quitarla**.
    if (( total < declared )); then
        echo "REGRESIÓN: $module declaraba $declared comprobaciones y solo quedan $total: se han retirado guardas" >&2
        failures=$((failures + 1))
    elif (( total > declared )); then
        echo "AVISO: $module pasa de $declared a $total comprobaciones; las nuevas no las ejercita nadie hasta que se demuestre" >&2
    fi
    printf '| `%s` | %d | %d |\n' "$module" "$total" "$dead" >> "$evidence_dir/tabla.md"
    printf '%-34s comprobaciones=%-4d nunca=%-4d techo=%s\n' \
        "$module" "$total" "$dead" "$ceiling" >&2
    if (( dead > ceiling )); then
        echo "REGRESIÓN: $module deja $dead guardas sin disparar, y el techo es $ceiling" >&2
        failures=$((failures + 1))
    fi
done

cat "$evidence_dir"/never-fired-*.txt | sort > "$evidence_dir/never-fired.txt"
cat "$evidence_dir/tabla.md" >&2

if (( failures != 0 )); then
    echo "la cobertura de guardas ha empeorado en $failures juez(es)" >&2
    exit 1
fi
echo "OK: ningún juez supera su techo de guardas sin disparar" >&2
