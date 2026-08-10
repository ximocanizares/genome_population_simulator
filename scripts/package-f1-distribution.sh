#!/usr/bin/env bash
set -euo pipefail

output_root="${1:-/tmp/gps-f1-distribution}"
binary="target/x86_64-unknown-linux-gnu/release/genome-population-simulator"
package_name="genome-population-simulator-f1-0.1.0-x86_64-unknown-linux-gnu"
stage="$output_root/$package_name"
archive="$output_root/$package_name.tar.gz"

test -x "$binary"
rm -rf -- "$stage"
mkdir -p "$stage/bin" "$stage/schemas" "$stage/licenses" "$stage/examples"
install -m 0755 "$binary" "$stage/bin/genome-population-simulator"
cp schemas/*.json "$stage/schemas/"
cp LICENSE "$stage/licenses/LICENSE"
cp tests/fixtures/config/general-seeded.toml "$stage/examples/general-seeded.toml"
cp distribution/f1/Dockerfile "$stage/Dockerfile"

tar --create --gzip --file "$archive" --directory "$output_root" "$package_name"
sha256sum "$archive"
tar --list --gzip --file "$archive"
