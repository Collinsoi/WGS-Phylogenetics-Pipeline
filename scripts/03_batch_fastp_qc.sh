#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd fastp
require_cmd fastqc
require_cmd multiqc
mkdir -p "$TRIMMED_DIR" "$PROJECT_ROOT/results/fastp" "$PROJECT_ROOT/results/fastqc"
load_samples

for index in "${!SAMPLE_IDS[@]}"; do
  sample="${SAMPLE_IDS[$index]}"
  run="${SAMPLE_RUNS[$index]}"
  r1="$RAW_DIR/${run}_1.fastq.gz"
  r2="$RAW_DIR/${run}_2.fastq.gz"
  [[ -s "$r1" && -s "$r2" ]] || {
    echo "Missing raw paired reads for $sample ($run)." >&2
    exit 1
  }
  mkdir -p "$PROJECT_ROOT/results/fastp/$sample"
  fastp \
    -i "$r1" -I "$r2" \
    -o "$TRIMMED_DIR/${sample}_R1.fastq.gz" \
    -O "$TRIMMED_DIR/${sample}_R2.fastq.gz" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html "$PROJECT_ROOT/results/fastp/$sample/fastp.html" \
    --json "$PROJECT_ROOT/results/fastp/$sample/fastp.json" \
    --thread "$THREADS"
  fastqc --threads "$THREADS" \
    "$TRIMMED_DIR/${sample}_R1.fastq.gz" \
    "$TRIMMED_DIR/${sample}_R2.fastq.gz" \
    --outdir "$PROJECT_ROOT/results/fastqc"
done

multiqc "$PROJECT_ROOT/results/fastp" "$PROJECT_ROOT/results/fastqc" \
  --outdir "$PROJECT_ROOT/results/multiqc"
