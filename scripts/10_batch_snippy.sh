#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd snippy
[[ -s "$REFERENCE" ]] || {
  echo "Reference FASTA not found: $REFERENCE" >&2
  exit 1
}
out_root="$PROJECT_ROOT/results/snippy"
mkdir -p "$out_root"
load_samples

for sample in "${SAMPLE_IDS[@]}"; do
  r1="$(reads_for_sample "$sample" R1)"
  r2="$(reads_for_sample "$sample" R2)"
  [[ -s "$r1" && -s "$r2" ]] || {
    echo "Paired reads missing for $sample: $r1 / $r2" >&2
    exit 1
  }
  snippy \
    --outdir "$out_root/$sample" \
    --ref "$REFERENCE" \
    --R1 "$r1" \
    --R2 "$r2" \
    --cpus "$THREADS"
done
