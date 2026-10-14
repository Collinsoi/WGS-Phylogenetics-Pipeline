#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd mash
[[ -s "$REFERENCE" ]] || {
  echo "Reference FASTA not found: $REFERENCE" >&2
  exit 1
}
out_dir="$PROJECT_ROOT/results/phylogeny/reference_check"
mkdir -p "$out_dir"
load_samples

mash sketch -o "$out_dir/NCTC11168" "$REFERENCE"
: > "$out_dir/assemblies.txt"
for sample in "${SAMPLE_IDS[@]}"; do
  assembly="$(assembly_for_sample "$sample")"
  [[ -s "$assembly" ]] || {
    echo "Assembly not found for $sample: $assembly (run 06_batch_spades.sh first)." >&2
    exit 1
  }
  mash sketch -o "$out_dir/$sample" "$assembly"
  mash dist "$out_dir/NCTC11168.msh" "$out_dir/$sample.msh" \
    > "$out_dir/${sample}_to_NCTC11168.tsv"
  printf '%s\n' "$assembly" >> "$out_dir/assemblies.txt"
done

mash sketch -l "$out_dir/assemblies.txt" -o "$out_dir/all15"
mash dist "$out_dir/all15.msh" "$out_dir/all15.msh" \
  > "$out_dir/all15_pairwise_mash.tsv"
