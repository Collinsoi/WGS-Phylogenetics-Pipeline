#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd abricate
database="${ABRICATE_DB:-vfdb}"
out_dir="$PROJECT_ROOT/results/virulence/$database"
mkdir -p "$out_dir"
load_samples

for sample in "${SAMPLE_IDS[@]}"; do
  assembly="$(assembly_for_sample "$sample")"
  [[ -s "$assembly" ]] || {
    echo "Assembly not found for $sample: $assembly" >&2
    exit 1
  }
  abricate --db "$database" --minid 80 --mincov 80 "$assembly" \
    > "$out_dir/${sample}.tsv"
done

reports=()
for sample in "${SAMPLE_IDS[@]}"; do
  report="$out_dir/${sample}.tsv"
  [[ -s "$report" ]] && reports+=("$report")
done
if [[ ${#reports[@]} -eq 0 ]]; then
  echo "No per-sample VFDB reports were produced." >&2
  exit 1
fi
abricate --summary "${reports[@]}" > "$out_dir/vfdb_summary.tsv"
