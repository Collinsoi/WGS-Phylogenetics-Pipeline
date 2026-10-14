#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd mlst
out_dir="$PROJECT_ROOT/results/mlst_batch"
mkdir -p "$out_dir"
result="$out_dir/mlst_results.tsv"
: > "$result"
load_samples

for sample in "${SAMPLE_IDS[@]}"; do
  assembly="$(assembly_for_sample "$sample")"
  [[ -s "$assembly" ]] || {
    echo "Assembly not found for $sample: $assembly" >&2
    exit 1
  }
  mlst "$assembly" >> "$result"
done
