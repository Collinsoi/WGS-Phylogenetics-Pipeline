#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd amrfinder
out_root="$PROJECT_ROOT/results/amrfinder_batch"
mkdir -p "$out_root"
load_samples

for sample in "${SAMPLE_IDS[@]}"; do
  assembly="$(assembly_for_sample "$sample")"
  [[ -s "$assembly" ]] || {
    echo "Assembly not found for $sample: $assembly" >&2
    exit 1
  }
  out_dir="$out_root/$sample"
  mkdir -p "$out_dir"
  amrfinder \
    -n "$assembly" \
    -O Campylobacter \
    --plus \
    --mutation_all "$out_dir/${sample}_all_mutations.tsv" \
    --output "$out_dir/${sample}_amrfinder.tsv" \
    --threads "$THREADS"
done
