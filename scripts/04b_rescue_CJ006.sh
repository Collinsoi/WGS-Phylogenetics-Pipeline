#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd kraken2
require_cmd extract_kraken_reads.py
require_cmd spades.py
: "${KRAKEN_DB:?Set KRAKEN_DB to a Kraken2 database directory.}"
TAXID="${KRAKEN_TAXID:-194}"
sample="$RESCUE_SAMPLE"
read_dir="$PROJECT_ROOT/results/decontamination/$sample"
initial_assembly="$ASSEMBLY_ROOT/$sample/contigs.fasta"
read_report="$read_dir/${sample}_reads_report.txt"
read_output="$read_dir/${sample}_reads.kraken"
mkdir -p "$read_dir" "$PROJECT_ROOT/results/kraken2_assembly/$sample"

r1="$TRIMMED_DIR/${sample}_R1.fastq.gz"
r2="$TRIMMED_DIR/${sample}_R2.fastq.gz"
[[ -s "$r1" && -s "$r2" ]] || {
  echo "Missing trimmed paired reads for $sample." >&2
  exit 1
}
[[ -s "$initial_assembly" ]] || {
  echo "Initial assembly not found: $initial_assembly (run 06_batch_spades.sh first)." >&2
  exit 1
}

kraken2 \
  --db "$KRAKEN_DB" --threads "$THREADS" \
  --report "$read_report" --output "$read_output" \
  --paired "$r1" "$r2"
kraken2 \
  --db "$KRAKEN_DB" --threads "$THREADS" \
  --report "$PROJECT_ROOT/results/kraken2_assembly/$sample/${sample}_contigs_report.txt" \
  --output "$PROJECT_ROOT/results/kraken2_assembly/$sample/${sample}_contigs_output.txt" \
  "$initial_assembly"

extract_kraken_reads.py \
  -k "$read_output" -s1 "$r1" -s2 "$r2" \
  -o "$read_dir/${sample}_Campylobacter_R1.fastq" \
  -o2 "$read_dir/${sample}_Campylobacter_R2.fastq" \
  -t "$TAXID" --include-children \
  --report "$read_report" --fastq-output
gzip -f "$read_dir/${sample}_Campylobacter_R1.fastq" \
  "$read_dir/${sample}_Campylobacter_R2.fastq"

clean_assembly="$ASSEMBLY_ROOT/${sample}_clean"
spades.py \
  -1 "$read_dir/${sample}_Campylobacter_R1.fastq.gz" \
  -2 "$read_dir/${sample}_Campylobacter_R2.fastq.gz" \
  -o "$clean_assembly" -t "$THREADS" --isolate

mkdir -p "$PROJECT_ROOT/results/kraken2_assembly/${sample}_clean"
kraken2 \
  --db "$KRAKEN_DB" --threads "$THREADS" \
  --report "$PROJECT_ROOT/results/kraken2_assembly/${sample}_clean/${sample}_clean_report.txt" \
  --output /dev/null "$clean_assembly/contigs.fasta"
