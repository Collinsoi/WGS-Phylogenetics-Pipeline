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
mkdir -p "$core_dir"
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
