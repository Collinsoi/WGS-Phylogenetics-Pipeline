#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd prefetch
require_cmd fasterq-dump
mkdir -p "$RAW_DIR" "$PROJECT_ROOT/data/sra"
load_samples

for index in "${!SAMPLE_IDS[@]}"; do
  sample="${SAMPLE_IDS[$index]}"
  run="${SAMPLE_RUNS[$index]}"
  r1="$RAW_DIR/${run}_1.fastq.gz"
  r2="$RAW_DIR/${run}_2.fastq.gz"
  if [[ -s "$r1" && -s "$r2" ]]; then
    echo "$sample ($run): paired FASTQ files already exist; skipping."
    continue
  fi

  echo "Downloading $sample : $run"
  prefetch "$run" --output-directory "$PROJECT_ROOT/data/sra"
  sra="$PROJECT_ROOT/data/sra/$run/$run.sra"
  if [[ ! -s "$sra" ]]; then
    echo "SRA file not found after prefetch: $sra" >&2
    exit 1
  fi
  fasterq-dump "$sra" --split-files --threads "$THREADS" --outdir "$RAW_DIR"
  for mate in 1 2; do
    fastq="$RAW_DIR/${run}_${mate}.fastq"
    if [[ ! -s "$fastq" ]]; then
      echo "Expected paired FASTQ was not created: $fastq" >&2
      exit 1
    fi
    gzip "$fastq"
  done
done
