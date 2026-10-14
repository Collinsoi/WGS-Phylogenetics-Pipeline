#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$PROJECT_ROOT"

THREADS="${THREADS:-4}"
METADATA="${METADATA:-$PROJECT_ROOT/metadata/sample_metadata.csv}"
RAW_DIR="${RAW_DIR:-$PROJECT_ROOT/data/raw}"
TRIMMED_DIR="${TRIMMED_DIR:-$PROJECT_ROOT/data/trimmed}"
ASSEMBLY_ROOT="${ASSEMBLY_ROOT:-$PROJECT_ROOT/results/assembly_batch}"
REFERENCE="${REFERENCE:-$PROJECT_ROOT/references/species/C_jejuni_NCTC11168.fasta}"
RESCUE_SAMPLE="${RESCUE_SAMPLE:-CJ006}"

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "Required command not found: $1" >&2
    return 1
  fi
}

load_samples() {
  local sample run
  SAMPLE_IDS=()
  SAMPLE_RUNS=()
  if [[ ! -f "$METADATA" ]]; then
    echo "Sample metadata CSV not found: $METADATA" >&2
    return 1
  fi
  while IFS=, read -r sample run; do
    [[ "$sample" == "sample_id" ]] && continue
    [[ -n "$sample" && -n "$run" ]] || continue
    SAMPLE_IDS+=("$sample")
    SAMPLE_RUNS+=("$run")
  done < "$METADATA"
  if [[ ${#SAMPLE_IDS[@]} -eq 0 ]]; then
    echo "No sample rows found in $METADATA" >&2
    return 1
  fi
}

assembly_for_sample() {
  local sample="$1"
  local rescued="$ASSEMBLY_ROOT/${sample}_clean/contigs.fasta"
  if [[ "$sample" == "$RESCUE_SAMPLE" && -f "$rescued" ]]; then
    printf '%s' "$rescued"
  else
    printf '%s' "$ASSEMBLY_ROOT/$sample/contigs.fasta"
  fi
}

reads_for_sample() {
  local sample="$1"
  local mate="$2"
  local rescued="$PROJECT_ROOT/results/decontamination/$sample/${sample}_Campylobacter_${mate}.fastq.gz"
  if [[ "$sample" == "$RESCUE_SAMPLE" && -f "$rescued" ]]; then
    printf '%s' "$rescued"
  else
    printf '%s' "$TRIMMED_DIR/${sample}_${mate}.fastq.gz"
  fi
}
