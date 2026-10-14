#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$SCRIPT_DIR"
cd "$PROJECT_ROOT"

usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh [reference_fasta]

Runs the README workflow using metadata/sample_metadata.csv and writes data
and results beneath this repository. Configure the Conda environment names
and optional stages with environment variables; see README.md.

Required:
  KRAKEN_DB             Kraken2 database directory

Optional:
  THREADS               Threads per tool (default: 4)
  RESCUE_SAMPLE         Sample to genus-filter (default: CJ006)
  KRAKEN_TAXID          Rescue taxon ID (default: 194)
  RUN_VFDB              Set to 0 to skip ABRicate/VFDB (default: 1)
  SKIP_DOWNLOAD         Set to 1 to reuse existing raw FASTQ files
  DRY_RUN               Set to 1 to print planned stages without running them
  MAIN_ENV              Main environment name (default: campy-wgs)
  SNIPPY_ENV            Snippy environment name (default: snippy-env)
  GUBBINS_ENV           Gubbins environment name (default: gubbins-env)
  IQTREE_ENV            IQ-TREE environment name (default: iqtree-env)
  ABRICATE_ENV          ABRicate environment name (default: abricate-env)
  IQTREE_BIN            IQ-TREE executable (default: iqtree)
EOF
}

if [[ $# -gt 1 ]]; then
  usage >&2
  exit 2
fi
if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
  usage
  exit 0
fi

THREADS="${THREADS:-4}"
REFERENCE="${1:-${REFERENCE:-$PROJECT_ROOT/references/species/C_jejuni_NCTC11168.fasta}}"
RESCUE_SAMPLE="${RESCUE_SAMPLE:-CJ006}"
KRAKEN_TAXID="${KRAKEN_TAXID:-194}"
RUN_VFDB="${RUN_VFDB:-1}"
SKIP_DOWNLOAD="${SKIP_DOWNLOAD:-0}"
DRY_RUN="${DRY_RUN:-0}"
MAIN_ENV="${MAIN_ENV:-campy-wgs}"
SNIPPY_ENV="${SNIPPY_ENV:-snippy-env}"
GUBBINS_ENV="${GUBBINS_ENV:-gubbins-env}"
IQTREE_ENV="${IQTREE_ENV:-iqtree-env}"
ABRICATE_ENV="${ABRICATE_ENV:-abricate-env}"
IQTREE_BIN="${IQTREE_BIN:-iqtree}"
LOG_DIR="$PROJECT_ROOT/results/logs"

if [[ ! "$THREADS" =~ ^[1-9][0-9]*$ ]]; then
  echo "THREADS must be a positive integer." >&2
  exit 2
fi
if [[ ! -s "$REFERENCE" ]]; then
  echo "Reference FASTA not found: $REFERENCE" >&2
  exit 1
fi
if [[ "$DRY_RUN" != 1 && -z "${KRAKEN_DB:-}" ]]; then
  echo "KRAKEN_DB must point to a Kraken2 database directory." >&2
  exit 1
fi
if [[ "$RUN_VFDB" != 0 && "$RUN_VFDB" != 1 ]]; then
  echo "RUN_VFDB must be 0 or 1." >&2
  exit 2
fi
if [[ "$SKIP_DOWNLOAD" != 0 && "$SKIP_DOWNLOAD" != 1 ]]; then
  echo "SKIP_DOWNLOAD must be 0 or 1." >&2
  exit 2
fi

if [[ "$DRY_RUN" != 1 ]] && ! command -v conda >/dev/null 2>&1; then
  echo "conda was not found on PATH; initialize Conda before running this script." >&2
  exit 1
fi

export THREADS REFERENCE RESCUE_SAMPLE KRAKEN_TAXID RUN_VFDB
export KRAKEN_DB="${KRAKEN_DB:-}"
export IQTREE_BIN

check_tool() {
  local env_name="$1"
  local executable="$2"
  if ! conda run --no-capture-output -n "$env_name" \
    bash -c 'command -v "$1" >/dev/null 2>&1' bash "$executable"; then
    echo "Required executable '$executable' is unavailable in Conda environment '$env_name'." >&2
    exit 1
  fi
}

run_stage() {
  local label="$1"
  local env_name="$2"
  local script="$3"

  if [[ "$DRY_RUN" == 1 ]]; then
    printf 'DRY RUN: conda run -n %q bash %q\n' "$env_name" "$script"
    return
  fi

  echo "===== $label (environment: $env_name) ====="
  if ! conda run --no-capture-output -n "$env_name" \
    bash "$PROJECT_ROOT/scripts/$script" 2>&1 \
    | tee "$LOG_DIR/run_pipeline.${script%.sh}.log"; then
    echo "Pipeline stopped: $label failed. See $LOG_DIR/run_pipeline.${script%.sh}.log" >&2
    exit 1
  fi
}

if [[ "$DRY_RUN" != 1 ]]; then
  mkdir -p "$LOG_DIR"
  check_tool "$MAIN_ENV" prefetch
  check_tool "$MAIN_ENV" fasterq-dump
  check_tool "$MAIN_ENV" fastp
  check_tool "$MAIN_ENV" fastqc
  check_tool "$MAIN_ENV" multiqc
  check_tool "$MAIN_ENV" kraken2
  check_tool "$MAIN_ENV" extract_kraken_reads.py
  check_tool "$MAIN_ENV" mash
  check_tool "$MAIN_ENV" spades.py
  check_tool "$MAIN_ENV" quast.py
  check_tool "$MAIN_ENV" mlst
  check_tool "$MAIN_ENV" amrfinder
  check_tool "$SNIPPY_ENV" snippy
  check_tool "$SNIPPY_ENV" snippy-core
  check_tool "$GUBBINS_ENV" run_gubbins.py
  check_tool "$IQTREE_ENV" "$IQTREE_BIN"
  if [[ "$RUN_VFDB" == 1 ]]; then
    check_tool "$ABRICATE_ENV" abricate
  fi
fi

if [[ "$SKIP_DOWNLOAD" != 1 ]]; then
  run_stage "Download SRA reads" "$MAIN_ENV" "02_download_remaining.sh"
else
  echo "Skipping SRA download (SKIP_DOWNLOAD=1)."
fi
run_stage "Read trimming and QC" "$MAIN_ENV" "03_batch_fastp_qc.sh"
run_stage "Read contamination screening" "$MAIN_ENV" "04_batch_kraken2.sh"
run_stage "Initial SPAdes assemblies" "$MAIN_ENV" "06_batch_spades.sh"
run_stage "CJ006 contamination rescue" "$MAIN_ENV" "04b_rescue_CJ006.sh"
run_stage "Mash reference and pairwise checks" "$MAIN_ENV" "05_batch_mash.sh"
run_stage "QUAST assembly evaluation" "$MAIN_ENV" "06b_batch_quast.sh"
run_stage "MLST typing" "$MAIN_ENV" "07_batch_mlst.sh"
run_stage "AMRFinderPlus screening" "$MAIN_ENV" "08_batch_amrfinder.sh"

if [[ "$RUN_VFDB" == 1 ]]; then
  run_stage "VFDB virulence screening" "$ABRICATE_ENV" "09_batch_vfdb.sh"
else
  echo "Skipping VFDB screening (RUN_VFDB=0)."
fi

run_stage "Snippy variant calling" "$SNIPPY_ENV" "10_batch_snippy.sh"
run_stage "Snippy core alignment" "$SNIPPY_ENV" "11_batch_phylogeny.sh"
run_stage "Gubbins recombination analysis" "$GUBBINS_ENV" "12_batch_gubbins.sh"
run_stage "IQ-TREE maximum-likelihood analysis" "$IQTREE_ENV" "13_batch_iqtree.sh"

echo "Pipeline completed. Results: $PROJECT_ROOT/results"
