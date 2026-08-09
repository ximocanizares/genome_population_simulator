#!/usr/bin/env bash
set -euo pipefail

project_root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
fixtures="$project_root/tests/fixtures/validators"

if (( $# != 1 )); then
    echo "usage: $0 EVIDENCE_DIRECTORY" >&2
    exit 2
fi
evidence_dir=$1
if [[ -e "$evidence_dir" ]]; then
    echo "evidence directory already exists: $evidence_dir" >&2
    exit 2
fi
mkdir -p -- "$evidence_dir/inputs"

samtools_version=$(samtools --version | sed -n '1s/^samtools //p')
htslib_version=$(samtools --version | sed -n '2s/^Using htslib //p')
sha256sum_version=$(sha256sum --version | sed -n '1s/^sha256sum (GNU coreutils) //p')
tar_version=$(tar --version | sed -n '1s/^tar (GNU tar) //p')
python_version=$(python3 --version | sed -n 's/^Python //p')
jsonschema_version=$(python3 -c 'import importlib.metadata; print(importlib.metadata.version("jsonschema"))')
docker_version=$(docker --version | sed -n 's/^Docker version \([^,]*\),.*/\1/p')
docker_build=$(docker --version | sed -n 's/^Docker version [^,]*, build \([^ ]*\)$/\1/p')

[[ "$samtools_version" == "1.21" ]] || { echo "unexpected samtools version: $samtools_version" >&2; exit 1; }
[[ "$htslib_version" == "1.21" ]] || { echo "unexpected htslib version: $htslib_version" >&2; exit 1; }
[[ "$sha256sum_version" == "9.4" ]] || { echo "unexpected sha256sum version: $sha256sum_version" >&2; exit 1; }
[[ "$tar_version" == "1.35" ]] || { echo "unexpected tar version: $tar_version" >&2; exit 1; }
[[ "$python_version" == "3.12.3" ]] || { echo "unexpected CPython version: $python_version" >&2; exit 1; }
[[ "$jsonschema_version" == "4.10.3" ]] || { echo "unexpected jsonschema version: $jsonschema_version" >&2; exit 1; }
[[ "$docker_version" == "29.6.1" ]] || { echo "unexpected Docker version: $docker_version" >&2; exit 1; }
[[ "$docker_build" == "8900f1d" ]] || { echo "unexpected Docker build: $docker_build" >&2; exit 1; }

cp -- \
    "$fixtures/fasta-without-header.fa" \
    "$fixtures/payload.txt" \
    "$fixtures/payload.sha256" \
    "$fixtures/not-an-archive.tar.gz" \
    "$fixtures/diagnostics-unknown-field.json" \
    "$project_root/schemas/diagnostics-v1.schema.json" \
    "$evidence_dir/inputs/"

platform_name=$(
    . /etc/os-release
    printf '%s %s; %s' "$NAME" "$VERSION_ID" "$(uname -m)"
)
if [[ $(uname -r) == *microsoft-standard-WSL2* ]]; then
    platform_name="$platform_name; WSL2"
else
    platform_name="$platform_name; native_or_OCI"
fi
filesystem=$(findmnt -n -o FSTYPE --target "$project_root")

run_and_capture() {
    local validator_class=$1
    local runner=$2
    local stdout_file="$evidence_dir/$validator_class.stdout"
    local stderr_file="$evidence_dir/$validator_class.stderr"
    if "$runner" >"$stdout_file" 2>"$stderr_file"; then
        validator_status=0
    else
        validator_status=$?
    fi
}

run_fasta() (
    cd -- "$evidence_dir/inputs"
    samtools faidx fasta-without-header.fa
)

run_sha256() (
    cd -- "$evidence_dir/inputs"
    sha256sum --check payload.sha256
)

run_tar_gz() (
    cd -- "$evidence_dir/inputs"
    tar --list --gzip --file not-an-archive.tar.gz
)

run_json_schema() (
    cd -- "$evidence_dir/inputs"
    python3 -m jsonschema -i diagnostics-unknown-field.json diagnostics-v1.schema.json
)

invalid_image="gps-invalid@sha256:0000000000000000000000000000000000000000000000000000000000000000"
run_oci() {
    docker image inspect "$invalid_image"
}

fail_unexpected() {
    local validator_class=$1
    local explanation=$2
    echo "$validator_class did not produce the expected semantic rejection: $explanation" >&2
    echo "stdout: $evidence_dir/$validator_class.stdout" >&2
    echo "stderr: $evidence_dir/$validator_class.stderr" >&2
    exit 1
}

write_metadata() {
    local validator_class=$1
    local tool=$2
    local tool_version=$3
    local command=$4
    local fixture=$5
    local artifact_sha256=$6
    local observed=$7
    local stdout_sha256 stderr_sha256
    stdout_sha256=$(sha256sum "$evidence_dir/$validator_class.stdout" | cut -d ' ' -f 1)
    stderr_sha256=$(sha256sum "$evidence_dir/$validator_class.stderr" | cut -d ' ' -f 1)
    printf '%s\n' \
        'evidence_version=1' \
        "validator_class=$validator_class" \
        "tool=$tool" \
        "tool_version=$tool_version" \
        "platform=$platform_name" \
        "filesystem=$filesystem" \
        "command=$command" \
        "fixture=$fixture" \
        'expected=rejection' \
        "observed=exit_code=$validator_status;$observed" \
        "stdout_sha256=$stdout_sha256" \
        "stderr_sha256=$stderr_sha256" \
        "artifact_sha256=$artifact_sha256" \
        >"$evidence_dir/$validator_class.metadata"
    printf 'PASS class=%s status=%s stdout_sha256=%s stderr_sha256=%s artifact_sha256=%s\n' \
        "$validator_class" "$validator_status" "$stdout_sha256" "$stderr_sha256" "$artifact_sha256"
}

run_and_capture fasta run_fasta
[[ $validator_status == 1 ]] || fail_unexpected fasta "exit code $validator_status, expected 1"
grep -Fq 'Format error' "$evidence_dir/fasta.stderr" || fail_unexpected fasta 'missing format-error diagnostic'
grep -Fq 'Could not build fai index' "$evidence_dir/fasta.stderr" || fail_unexpected fasta 'missing faidx rejection diagnostic'
[[ ! -e "$evidence_dir/inputs/fasta-without-header.fa.fai" ]] || fail_unexpected fasta 'validator left an index for the corrupt FASTA'
write_metadata \
    fasta samtools "samtools=$samtools_version;htslib=$htslib_version" \
    'samtools faidx fasta-without-header.fa' \
    'tests/fixtures/validators/fasta-without-header.fa' \
    "$(sha256sum "$fixtures/fasta-without-header.fa" | cut -d ' ' -f 1)" \
    'rejected_invalid_fasta_structure'

run_and_capture sha256 run_sha256
[[ $validator_status == 1 ]] || fail_unexpected sha256 "exit code $validator_status, expected 1"
grep -Fxq 'payload.txt: FAILED' "$evidence_dir/sha256.stdout" || fail_unexpected sha256 'missing checksum mismatch result'
grep -Fq 'computed checksum did NOT match' "$evidence_dir/sha256.stderr" || fail_unexpected sha256 'missing checksum mismatch diagnostic'
write_metadata \
    sha256 sha256sum "$sha256sum_version" \
    'sha256sum --check payload.sha256' \
    'tests/fixtures/validators/payload.sha256 + payload.txt' \
    "$(sha256sum "$fixtures/payload.sha256" | cut -d ' ' -f 1)" \
    'rejected_digest_mismatch'

run_and_capture tar_gz run_tar_gz
[[ $validator_status == 2 ]] || fail_unexpected tar_gz "exit code $validator_status, expected 2"
grep -Fq 'not in gzip format' "$evidence_dir/tar_gz.stderr" || fail_unexpected tar_gz 'missing gzip-format diagnostic'
grep -Fq 'Error is not recoverable' "$evidence_dir/tar_gz.stderr" || fail_unexpected tar_gz 'missing tar rejection diagnostic'
write_metadata \
    tar_gz tar "$tar_version" \
    'tar --list --gzip --file not-an-archive.tar.gz' \
    'tests/fixtures/validators/not-an-archive.tar.gz' \
    "$(sha256sum "$fixtures/not-an-archive.tar.gz" | cut -d ' ' -f 1)" \
    'rejected_non_gzip_archive'

run_and_capture json_schema run_json_schema
[[ $validator_status == 1 ]] || fail_unexpected json_schema "exit code $validator_status, expected 1"
grep -Fq "Additional properties are not allowed ('unknown_field' was unexpected)" \
    "$evidence_dir/json_schema.stderr" || fail_unexpected json_schema 'missing unknown-field validation diagnostic'
write_metadata \
    json_schema 'CPython/jsonschema' "CPython=$python_version;jsonschema=$jsonschema_version" \
    'python3 -m jsonschema -i diagnostics-unknown-field.json diagnostics-v1.schema.json' \
    'tests/fixtures/validators/diagnostics-unknown-field.json' \
    "$(sha256sum "$fixtures/diagnostics-unknown-field.json" | cut -d ' ' -f 1)" \
    'rejected_unknown_top_level_field'

run_and_capture oci run_oci
[[ $validator_status == 1 ]] || fail_unexpected oci "exit code $validator_status, expected 1"
if grep -Eqi 'cannot connect|permission denied|daemon.*not running' "$evidence_dir/oci.stderr"; then
    fail_unexpected oci 'Docker infrastructure is unavailable'
fi
if ! grep -Eqi 'no such (object|image)|not found' "$evidence_dir/oci.stderr"; then
    fail_unexpected oci 'Docker failed for an unexpected reason'
fi
write_metadata \
    oci docker "version=$docker_version;build=$docker_build" \
    "docker image inspect $invalid_image" \
    "$invalid_image" \
    '0000000000000000000000000000000000000000000000000000000000000000' \
    'rejected_missing_digest_pinned_image'
