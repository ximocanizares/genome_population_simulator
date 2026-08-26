#!/usr/bin/env bash
set -euo pipefail

# Reproductor sintético de la ruta FASTA suministrado + reducción. No usa
# referencias reales ni escribe resultados dentro del repositorio.
repetitions="${1:-3}"
binary="${GPS_BENCHMARK_BINARY:-target/release/genome-population-simulator}"
template="benchmarks/f1/reduced-reference-template.toml"
workspace="$(mktemp -d /tmp/gps-f1-reduced-reference-benchmark.XXXXXX)"
trap 'rm -rf -- "$workspace"' EXIT

test -x "$binary"
test -f "$template"
case "$repetitions" in
  ''|*[!0-9]*|0) echo "las repeticiones deben ser un entero positivo" >&2; exit 2 ;;
esac
binary_path="$(cd "$(dirname -- "$binary")" && pwd -P)/$(basename -- "$binary")"
binary_sha256="$(sha256sum "$binary" | cut -d' ' -f1)"
binary_version="$("$binary" --version)"
template_sha256="$(sha256sum "$template" | cut -d' ' -f1)"
git_revision="$(git rev-parse --verify HEAD)"
git_tracked_diff_sha256="$(git diff --binary HEAD | sha256sum | cut -d' ' -f1)"
git_untracked_paths_sha256="$(git ls-files --others --exclude-standard -z | sha256sum | cut -d' ' -f1)"
if git diff --quiet --ignore-submodules HEAD -- && [ -z "$(git ls-files --others --exclude-standard)" ]; then
  git_worktree="clean"
else
  git_worktree="dirty"
fi

write_source() {
  local bases="$1"
  {
    printf '>chr-benchmark\n'
    head -c "$bases" /dev/zero | tr '\000' 'A'
    printf '\n'
  } > "$workspace/source.fa"
}

measure() {
  local label="$1"
  local config="$2"
  local destination="$3"
  local metrics="$4"
  local started finished elapsed
  started="$(date +%s%N)"
  if [ -n "$destination" ]; then
    /usr/bin/time -f '%U\t%S\t%P\t%M\t%c\t%w' -o "$metrics" \
      "$binary" "$label" "$config" --output "$destination" >/dev/null
  else
    /usr/bin/time -f '%U\t%S\t%P\t%M\t%c\t%w' -o "$metrics" \
      "$binary" "$label" "$config" >/dev/null
  fi
  finished="$(date +%s%N)"
  elapsed="$(awk -v started="$started" -v finished="$finished" 'BEGIN { printf "%.6f", (finished - started) / 1000000000 }')"
  printf '%s' "$elapsed"
}

# Tamaños suficientes para separar E/S y materialización; los tres casos de
# 100 MB mantienen la misma entrada y cambian solo el número de bloques. Así
# el coste del selector canónico se observa frente a 10, 100 y 1.000 bloques,
# sin un recorrido proporcional a las bases de holgura.
scenarios="${GPS_BENCHMARK_SCENARIOS:-small:1000000:500000:10000 direct-100mb-10-blocks:100000000:10000000:1000000 direct-100mb-100-blocks:100000000:10000000:100000 direct-100mb-1000-blocks:100000000:10000000:10000}"
printf '# binary_path=%s\n' "$binary_path"
printf '# binary_sha256=%s\n' "$binary_sha256"
printf '# binary_version=%s\n' "$binary_version"
printf '# git_revision=%s\n' "$git_revision"
printf '# git_worktree=%s\n' "$git_worktree"
printf '# git_tracked_diff_sha256=%s\n' "$git_tracked_diff_sha256"
printf '# git_untracked_paths_sha256=%s\n' "$git_untracked_paths_sha256"
printf '# template_sha256=%s\n' "$template_sha256"
printf 'scenario\tsource_bases\ttarget_bases\tregion_size_bases\trepetition\tstage\telapsed_seconds\tuser_seconds\tsystem_seconds\tcpu_percent\tmax_rss_kib\tinvoluntary_context_switches\tvoluntary_context_switches\tderived_fasta_sha256\tmap_sha256\n'
for scenario_spec in $scenarios; do
  IFS=: read -r scenario source_bases target_bases region_size_bases <<< "$scenario_spec"
  write_source "$source_bases"
  expected_derived_sha=''
  expected_map_sha=''
  for ((repetition = 1; repetition <= repetitions; repetition++)); do
    config="$workspace/${scenario}-${repetition}.toml"
    destination="$workspace/${scenario}-${repetition}-result"
    validation_metrics="$workspace/${scenario}-${repetition}-validate.metrics"
    run_metrics="$workspace/${scenario}-${repetition}-run.metrics"
    sed -e "s/__LENGTH__/${source_bases}/" \
        -e "s/__TARGET_BASES__/${target_bases}/" \
        -e "s/__REGION_SIZE__/${region_size_bases}/" \
        "$template" > "$config"
    validation_elapsed="$(measure validate "$config" '' "$validation_metrics")"
    IFS=$'\t' read -r validation_user validation_system validation_cpu validation_rss validation_involuntary validation_voluntary < "$validation_metrics"
    printf '%s\t%s\t%s\t%s\t%s\tvalidation_preflight\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tNA\tNA\n' \
      "$scenario" "$source_bases" "$target_bases" "$region_size_bases" "$repetition" \
      "$validation_elapsed" "$validation_user" "$validation_system" "$validation_cpu" \
      "$validation_rss" "$validation_involuntary" "$validation_voluntary"
    run_elapsed="$(measure run "$config" "$destination" "$run_metrics")"
    IFS=$'\t' read -r run_user run_system run_cpu run_rss run_involuntary run_voluntary < "$run_metrics"
    derived_sha="$(sha256sum "$destination/resources/references/reference-small.fa" | cut -d' ' -f1)"
    map_sha="$(sha256sum "$destination/resources/references/reference-small.map.json" | cut -d' ' -f1)"
    if [ -n "$expected_derived_sha" ] && { [ "$derived_sha" != "$expected_derived_sha" ] || [ "$map_sha" != "$expected_map_sha" ]; }; then
      echo "hash inestable para $scenario" >&2
      exit 1
    fi
    expected_derived_sha="$derived_sha"
    expected_map_sha="$map_sha"
    printf '%s\t%s\t%s\t%s\t%s\tcomplete_run\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$scenario" "$source_bases" "$target_bases" "$region_size_bases" "$repetition" \
      "$run_elapsed" "$run_user" "$run_system" "$run_cpu" "$run_rss" \
      "$run_involuntary" "$run_voluntary" "$derived_sha" "$map_sha"
  done
done
