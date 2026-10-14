#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd spades.py
mkdir -p "$ASSEMBLY_ROOT"
load_samples

for sample in "${SAMPLE_IDS[@]}"; do
  r1="$TRIMMED_DIR/${sample}_R1.fastq.gz"
  r2="$TRIMMED_DIR/${sample}_R2.fastq.gz"
  out_dir="$ASSEMBLY_ROOT/$sample"
  if [[ -s "$out_dir/contigs.fasta" ]]; then
    echo "$sample: assembly already exists; skipping."
    continue
  fi
  [[ -s "$r1" && -s "$r2" ]] || {
    echo "Missing trimmed reads for $sample; run 03_batch_fastp_qc.sh first." >&2
    exit 1
  }
  spades.py -1 "$r1" -2 "$r2" -o "$out_dir" -t "$THREADS" --isolate
done
