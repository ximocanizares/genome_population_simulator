#!/usr/bin/env bash
set -euo pipefail

mode="${1:-A}"
repetitions="${2:-3}"
binary="target/x86_64-unknown-linux-gnu/release/genome-population-simulator"
template="benchmarks/f1/reference-template.toml"

case "$mode" in
  A) scenarios="A:100000" ;;
  B) scenarios="B1:1000000 B2:10000000 B3:100000000" ;;
  C) scenarios="C:1000000000" ;;
  *) echo "uso: $0 A|B|C [repeticiones]" >&2; exit 2 ;;
esac

test -x "$binary"
test -f "$template"
workspace="$(mktemp -d /tmp/gps-f1-benchmark.XXXXXX)"
trap 'rm -rf -- "$workspace"' EXIT

printf 'scenario\tbases\trepetition\tstage\telapsed_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kib\tinput_bytes\ttemporary_peak_bytes\ttemporary_residual_bytes\toutput_bytes\tthroughput_mb_s\tfiles\tvoluntary_context_switches\tinvoluntary_context_switches\n'
for scenario_spec in $scenarios; do
  scenario="${scenario_spec%%:*}"
  bases="${scenario_spec##*:}"
  repetition=1
  while [ "$repetition" -le "$repetitions" ]; do
    config="$workspace/${scenario}-${repetition}.toml"
    destination="$workspace/${scenario}-${repetition}-result"
    metrics="$workspace/${scenario}-${repetition}.metrics"
    validation_metrics="$workspace/${scenario}-${repetition}.validation.metrics"
    sed "s/__LENGTH__/${bases}/" "$template" > "$config"
    validation_started_ns="$(date +%s%N)"
    /usr/bin/time -f '%U\t%S\t%P\t%M\t%c\t%w' -o "$validation_metrics" \
      "$binary" validate "$config"
    validation_finished_ns="$(date +%s%N)"
    validation_elapsed="$(awk -v started="$validation_started_ns" -v finished="$validation_finished_ns" 'BEGIN { printf "%.6f", (finished - started) / 1000000000 }')"
    IFS=$'\t' read -r validation_user validation_system validation_cpu validation_rss validation_involuntary validation_voluntary < "$validation_metrics"
    input_bytes="$(stat -c %s "$config")"
    printf '%s\t%s\t%s\tvalidation_preflight\t%s\t%s\t%s\t%s\t%s\t%s\tnot_measured\t0\t0\tNA\t0\t%s\t%s\n' \
      "$scenario" "$bases" "$repetition" "$validation_elapsed" "$validation_user" \
      "$validation_system" "$validation_cpu" "$validation_rss" "$input_bytes" \
      "$validation_voluntary" "$validation_involuntary"
    started_ns="$(date +%s%N)"
    /usr/bin/time -f '%U\t%S\t%P\t%M\t%c\t%w' -o "$metrics" \
      "$binary" run "$config" --output "$destination"
    finished_ns="$(date +%s%N)"
    elapsed="$(awk -v started="$started_ns" -v finished="$finished_ns" 'BEGIN { printf "%.6f", (finished - started) / 1000000000 }')"
    IFS=$'\t' read -r user system cpu rss involuntary voluntary < "$metrics"
    output_bytes="$(du -sb "$destination" | cut -f1)"
    temporary_residual_bytes="$(find "$workspace" -maxdepth 1 -type d -name '*.gps-incomplete-*' -exec du -sb {} + | awk '{ total += $1 } END { print total + 0 }')"
    files="$(find "$destination" -type f | wc -l)"
    throughput="$(awk -v bases="$bases" -v elapsed="$elapsed" 'BEGIN { if (elapsed == 0) print "NA"; else printf "%.6f", (bases / 1000000) / elapsed }')"
    printf '%s\t%s\t%s\tcomplete_run\t%s\t%s\t%s\t%s\t%s\t%s\tnot_measured\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$scenario" "$bases" "$repetition" "$elapsed" "$user" "$system" "$cpu" "$rss" \
      "$input_bytes" "$temporary_residual_bytes" "$output_bytes" "$throughput" "$files" "$voluntary" "$involuntary"
    repetition=$((repetition + 1))
  done
done
