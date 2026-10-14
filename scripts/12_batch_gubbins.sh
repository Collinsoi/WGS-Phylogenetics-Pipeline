#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd run_gubbins.py
alignment="$PROJECT_ROOT/results/phylogeny/core/core.full.aln"
[[ -s "$alignment" ]] || {
  echo "Core full alignment not found: $alignment (run 11_batch_phylogeny.sh first)." >&2
  exit 1
}
out_dir="$PROJECT_ROOT/results/phylogeny/gubbins"
mkdir -p "$out_dir" "$PROJECT_ROOT/results/final_trees"

run_gubbins.py --prefix "$out_dir/campy" "$alignment"
for tree in "$out_dir"/*.tree "$out_dir"/*.tre "$out_dir"/*.nwk; do
  [[ -f "$tree" ]] && cp "$tree" "$PROJECT_ROOT/results/final_trees/"
done
