#!/usr/bin/env bash
# Keep plasmotrackr's copy of the golden Stores, and of the input they were run on, equal to the engine's.
#
#   tools/plasmotrackr/golden-sync.sh          copy test/data/golden/ into R/inst/extdata/golden/ and
#                                              examples/example_input.json into R/inst/extdata/
#   tools/plasmotrackr/golden-sync.sh --check  exit 1 when the two copies differ (CI)
#
# The engine's .gitattributes is not copied: R/.gitattributes covers the package copy, and a hidden file in inst/
# would be installed with the package.
set -euo pipefail

root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
src="$root/test/data/golden"
dst="$root/R/inst/extdata/golden"
input_src="$root/examples/example_input.json"
input_dst="$root/R/inst/extdata/example_input.json"

if [[ "${1:-}" == "--check" ]]; then
  if diff -r --exclude=.gitattributes "$src" "$dst" && cmp "$input_src" "$input_dst"; then
    echo "golden copy in sync"
  else
    echo "plasmotrackr's golden copy differs from the engine's: run tools/plasmotrackr/golden-sync.sh" >&2
    exit 1
  fi
else
  rm -rf "$dst"
  mkdir -p "$dst"
  (cd "$src" && tar --exclude=.gitattributes -cf - .) | (cd "$dst" && tar -xf -)
  cp "$input_src" "$input_dst"
  echo "copied $src -> $dst and $input_src -> $input_dst"
fi
