#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/common.sh"

require_cmd kraken2
: "${KRAKEN_DB:?Set KRAKEN_DB to a Kraken2 database directory.}"
mkdir -p "$PROJECT_ROOT/results/kraken2_batch"
load_samples

summary="$PROJECT_ROOT/results/kraken2_batch/kraken2_summary.csv"
printf 'sample_id,total_reads,classified_percent,unclassified_percent\n' > "$summary"
for sample in "${SAMPLE_IDS[@]}"; do
  r1="$TRIMMED_DIR/${sample}_R1.fastq.gz"
  r2="$TRIMMED_DIR/${sample}_R2.fastq.gz"
  [[ -s "$r1" && -s "$r2" ]] || {
    echo "Missing trimmed paired reads for $sample." >&2
    exit 1
  }
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired "$r1" "$r2" \
    --report "$PROJECT_ROOT/results/kraken2_batch/${sample}_report.txt" \
    --output "$PROJECT_ROOT/results/kraken2_batch/${sample}.kraken"
  awk -v sample="$sample" '
    $4 == "R" && $5 == 1 { total = $2 }
    $4 == "U" && $5 == 0 { unclassified = $1 }
    END {
      if (total == "" || unclassified == "") {
        print "Incomplete Kraken2 report for " sample > "/dev/stderr"
        exit 1
      }
      printf "%s,%d,%.2f,%.2f\n", sample, total, 100 - unclassified, unclassified
    }
  ' "$PROJECT_ROOT/results/kraken2_batch/${sample}_report.txt" >> "$summary"
done
