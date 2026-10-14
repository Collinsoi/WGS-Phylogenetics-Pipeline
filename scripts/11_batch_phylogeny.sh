#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd snippy-core
[[ -s "$REFERENCE" ]] || {
  echo "Reference FASTA not found: $REFERENCE" >&2
  exit 1
}
out_root="$PROJECT_ROOT/results/phylogeny"
core_dir="$out_root/core"
snippy_root="$PROJECT_ROOT/results/snippy"
mkdir -p "$core_dir" "$out_root/gubbins" "$out_root/iqtree" "$PROJECT_ROOT/results/final_trees"
load_samples

snippy_dirs=()
for sample in "${SAMPLE_IDS[@]}"; do
  sample_dir="$snippy_root/$sample"
  [[ -d "$sample_dir" ]] || {
    echo "Snippy output missing for $sample: $sample_dir" >&2
    exit 1
  }
  snippy_dirs+=("$sample_dir")
done

snippy-core --ref "$REFERENCE" --prefix "$core_dir/core" "${snippy_dirs[@]}"

if command -v run_gubbins.py >/dev/null 2>&1; then
  run_gubbins.py --prefix "$out_root/gubbins/campy" "$core_dir/core.full.aln"
else
  echo "run_gubbins.py not found; skipping recombination-aware inference." >&2
fi

iqtree_bin="${IQTREE_BIN:-iqtree}"
if command -v "$iqtree_bin" >/dev/null 2>&1; then
  "$iqtree_bin" \
    -s "$core_dir/core.aln" \
    -st DNA \
    -m MFP+ASC \
    -B 1000 \
    -alrt 1000 \
    -T AUTO \
    --prefix "$out_root/iqtree/campy_final"
else
  echo "IQ-TREE executable '$iqtree_bin' not found; skipping IQ-TREE." >&2
fi

for tree in "$out_root"/gubbins/*.tree "$out_root"/gubbins/*.tre "$out_root"/gubbins/*.nwk \
  "$out_root"/iqtree/*.treefile "$out_root"/iqtree/*.contree; do
  [[ -f "$tree" ]] && cp "$tree" "$PROJECT_ROOT/results/final_trees/"
done
