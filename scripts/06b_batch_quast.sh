#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd quast.py
out_dir="$PROJECT_ROOT/results/quast_batch"
mkdir -p "$out_dir"
load_samples

assemblies=()
for sample in "${SAMPLE_IDS[@]}"; do
  assembly="$(assembly_for_sample "$sample")"
  [[ -s "$assembly" ]] || {
    echo "Assembly not found for $sample: $assembly" >&2
    exit 1
  }
  assemblies+=("$assembly")
done

quast.py "${assemblies[@]}" -o "$out_dir" -t "$THREADS" --min-contig 500
