#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  ./run_pipeline.sh <input_fastq_dir> <output_dir> [reference_fasta]

Example:
  ./run_pipeline.sh /path/to/raw_fastq /path/to/output /path/to/reference/C_jejuni_NCTC11168.fasta

The script auto-detects paired-end FASTQ samples from the input directory and
writes outputs to the specified output folder. Set KRAKEN_DB to enable
contamination screening and the optional genus-level rescue for RESCUE_SAMPLE
(default: CJ006). Set ABRICATE_DB (default: vfdb) to select the virulence
database.
EOF
  exit 1
}

if [[ $# -lt 2 || $# -gt 3 ]]; then
  usage
fi

INPUT_DIR="${1}"
OUTPUT_DIR="${2}"
REFERENCE="${3:-}"
THREADS="${THREADS:-8}"
RESCUE_SAMPLE="${RESCUE_SAMPLE:-CJ006}"
KRAKEN_TAXID="${KRAKEN_TAXID:-194}"
ABRICATE_DB="${ABRICATE_DB:-vfdb}"

if [[ ! -d "$INPUT_DIR" ]]; then
  echo "Error: input FASTQ directory not found: $INPUT_DIR" >&2
  exit 1
fi

mkdir -p "$OUTPUT_DIR"

# Recommended defaults for a C. jejuni pipeline.
if [[ -z "$REFERENCE" ]]; then
  REFERENCE="${REFERENCE_GENOME:-${OUTPUT_DIR}/reference/C_jejuni_NCTC11168.fasta}"
fi

if [[ ! -f "$REFERENCE" ]]; then
  echo "Warning: reference FASTA not found at $REFERENCE. If you have one, pass it as the third argument."
  echo "Snippy and phylogeny steps will be skipped or will fail unless a valid reference is provided."
fi

LOG_DIR="$OUTPUT_DIR/logs"
mkdir -p "$LOG_DIR"

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    echo "Error: required command not found: $cmd" >&2
    exit 1
  fi
}

# Keep the wrapper strict about the core tools used in the workflow.
require_cmd fastp
require_cmd fastqc
require_cmd multiqc
require_cmd kraken2
require_cmd mash
require_cmd spades.py
require_cmd quast.py
require_cmd mlst
require_cmd amrfinder
require_cmd snippy
require_cmd snippy-core

# Some environments may not have Gubbins or RAxML in PATH; keep this step optional.
HAS_GUBBINS=false
if command -v run_gubbins.py >/dev/null 2>&1; then
  HAS_GUBBINS=true
fi

HAS_IQTREE=false
IQTREE_BIN="${IQTREE_BIN:-iqtree}"
if command -v "$IQTREE_BIN" >/dev/null 2>&1; then
  HAS_IQTREE=true
fi

HAS_ABRICATE=false
if command -v abricate >/dev/null 2>&1; then
  HAS_ABRICATE=true
fi

# A robust sample detector: supports names like SRR_1.fastq.gz, *_R1_001.fastq.gz, *_1.fastq.gz.
find_fastq_pairs() {
  local r1_list=()
  while IFS= read -r f; do
    r1_list+=("$f")
  done < <(find "$INPUT_DIR" -type f \( -iname "*_R1*.fastq.gz" -o -iname "*_R1*.fastq" -o -iname "*_1.fastq.gz" -o -iname "*_1.fastq" -o -iname "*.R1.fastq.gz" -o -iname "*.R1.fastq" \) | sort)
  if [[ ${#r1_list[@]} -eq 0 ]]; then
    return 0
  fi
  printf '%s\n' "${r1_list[@]}"
}

sample_name_from_r1() {
  local file="$1"
  local base
  base="$(basename "$file")"

  if [[ "$base" =~ _R1 ]]; then
    base="${base%%_R1*}"
  elif [[ "$base" =~ _1\.fastq ]]; then
    base="${base%%_1.fastq*}"
  elif [[ "$base" =~ \.R1\.fastq ]]; then
    base="${base%%.R1.fastq*}"
  elif [[ "$base" =~ _1$ ]]; then
    base="${base%%_1}"
  else
    base="${base%.*}"
  fi

  printf '%s' "$base"
}

match_r2_for_r1() {
  local r1="$1"
  local r2=""
  if [[ "$r1" == *_R1* ]]; then
    r2="${r1/_R1/_R2}"
  elif [[ "$r1" == *_1.fastq* ]]; then
    r2="${r1/_1.fastq/_2.fastq}"
  elif [[ "$r1" == *.R1.fastq* ]]; then
    r2="${r1/.R1.fastq/.R2.fastq}"
  elif [[ "$r1" == *_1.fq* ]]; then
    r2="${r1/_1.fq/_2.fq}"
  else
    r2="${r1%.fastq.gz}.R2.fastq.gz"
  fi

  if [[ -f "$r2" ]]; then
    printf '%s' "$r2"
  else
    echo "" >&2
  fi
}

assembly_dir_for_sample() {
  local sample="$1"
  local rescued_assembly="$OUTPUT_DIR/05_assemblies/${sample}_clean"
  if [[ "$sample" == "$RESCUE_SAMPLE" && -f "$rescued_assembly/contigs.fasta" ]]; then
    printf '%s' "$rescued_assembly"
  else
    printf '%s' "$OUTPUT_DIR/05_assemblies/$sample"
  fi
}

reads_for_sample() {
  local sample="$1"
  local mate="$2"
  local rescued_reads="$OUTPUT_DIR/04_decontamination/${sample}/${sample}_Campylobacter_${mate}.fastq.gz"
  if [[ "$sample" == "$RESCUE_SAMPLE" && -f "$rescued_reads" ]]; then
    printf '%s' "$rescued_reads"
  else
    printf '%s' "$OUTPUT_DIR/01_trimmed_reads/${sample}_${mate}.fastq.gz"
  fi
}

run_step_1_fastp() {
  local sample="$1"
  local r1="$2"
  local r2="$3"
  local trim_dir="$OUTPUT_DIR/01_trimmed_reads"
  local qc_dir="$OUTPUT_DIR/01_fastp_reports"
  mkdir -p "$trim_dir" "$qc_dir/$sample"

  local out_r1="$trim_dir/${sample}_R1.fastq.gz"
  local out_r2="$trim_dir/${sample}_R2.fastq.gz"

  echo "[Step 1/10] fastp: $sample"
  fastp \
    -i "$r1" \
    -I "$r2" \
    -o "$out_r1" \
    -O "$out_r2" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --html "$qc_dir/$sample/fastp.html" \
    --json "$qc_dir/$sample/fastp.json" \
    --thread "$THREADS" >"$LOG_DIR/${sample}.fastp.log" 2>&1
}

run_step_2_fastqc_multiqc() {
  local sample="$1"
  local trim_dir="$OUTPUT_DIR/01_trimmed_reads"
  local qc_dir="$OUTPUT_DIR/02_qc_reports"
  mkdir -p "$qc_dir"

  echo "[Step 2/10] fastqc + multiqc: $sample"
  fastqc -t "$THREADS" "$trim_dir/${sample}_R1.fastq.gz" "$trim_dir/${sample}_R2.fastq.gz" -o "$qc_dir" >"$LOG_DIR/${sample}.fastqc.log" 2>&1 || true

  multiqc "$qc_dir" -o "$OUTPUT_DIR/03_multiqc" >"$LOG_DIR/${sample}.multiqc.log" 2>&1 || true
}

run_step_3_kraken2() {
  local sample="$1"
  local trim_dir="$OUTPUT_DIR/01_trimmed_reads"
  local kraken_dir="$OUTPUT_DIR/04_kraken2"
  mkdir -p "$kraken_dir"

  if [[ -z "${KRAKEN_DB:-}" ]]; then
    echo "[Step 3] kraken2: $sample (KRAKEN_DB unset; skipping contamination screen)"
    return 0
  fi

  echo "[Step 3] kraken2: $sample"
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --paired "$trim_dir/${sample}_R1.fastq.gz" "$trim_dir/${sample}_R2.fastq.gz" \
    --output "$kraken_dir/${sample}.kraken.out" \
    --report "$kraken_dir/${sample}.kraken.report" \
    >"$LOG_DIR/${sample}.kraken2.log" 2>&1
}

run_step_3b_rescue_cj006() {
  local sample="$1"
  local read_kraken_dir="$OUTPUT_DIR/04_kraken2"
  local contig_kraken_dir="$OUTPUT_DIR/04_kraken2_assembly/$sample"
  local rescue_dir="$OUTPUT_DIR/04_decontamination/$sample"
  local original_assembly="$OUTPUT_DIR/05_assemblies/$sample/contigs.fasta"
  local clean_assembly="$OUTPUT_DIR/05_assemblies/${sample}_clean"
  local r1="$OUTPUT_DIR/01_trimmed_reads/${sample}_R1.fastq.gz"
  local r2="$OUTPUT_DIR/01_trimmed_reads/${sample}_R2.fastq.gz"
  local out_r1="$rescue_dir/${sample}_Campylobacter_R1.fastq"
  local out_r2="$rescue_dir/${sample}_Campylobacter_R2.fastq"

  if [[ "$sample" != "$RESCUE_SAMPLE" ]]; then
    return 0
  fi
  if [[ -z "${KRAKEN_DB:-}" ]]; then
    echo "Skipping genus-level rescue for $sample: set KRAKEN_DB to enable it." >&2
    return 0
  fi
  if [[ ! -f "$read_kraken_dir/${sample}.kraken.out" || ! -f "$read_kraken_dir/${sample}.kraken.report" ]]; then
    echo "Cannot rescue $sample: read-level Kraken2 classification output/report is missing." >&2
    return 1
  fi
  if [[ ! -f "$original_assembly" ]]; then
    echo "Cannot rescue $sample: initial assembly is missing: $original_assembly" >&2
    return 1
  fi
  require_cmd extract_kraken_reads.py
  mkdir -p "$rescue_dir" "$contig_kraken_dir"

  echo "[Step 5b] Screening initial assembly for $sample"
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --report "$contig_kraken_dir/${sample}_contigs_report.txt" \
    --output "$contig_kraken_dir/${sample}_contigs_output.txt" \
    "$original_assembly" >"$LOG_DIR/${sample}.contigs_kraken2.log" 2>&1

  echo "[Step 5b] Extracting taxid $KRAKEN_TAXID reads for $sample"
  extract_kraken_reads.py \
    -k "$read_kraken_dir/${sample}.kraken.out" \
    -s1 "$r1" \
    -s2 "$r2" \
    -o "$out_r1" \
    -o2 "$out_r2" \
    -t "$KRAKEN_TAXID" \
    --include-children \
    --report "$read_kraken_dir/${sample}.kraken.report" \
    --fastq-output >"$LOG_DIR/${sample}.extract_kraken_reads.log" 2>&1
  gzip -f "$out_r1" "$out_r2"

  echo "[Step 5c] Reassembling genus-filtered reads for $sample"
  spades.py \
    -1 "${out_r1}.gz" \
    -2 "${out_r2}.gz" \
    -o "$clean_assembly" \
    -t "$THREADS" \
    --isolate >"$LOG_DIR/${sample}.clean_spades.log" 2>&1

  local clean_contig_report="$OUTPUT_DIR/04_kraken2_assembly/${sample}_clean"
  mkdir -p "$clean_contig_report"
  kraken2 \
    --db "$KRAKEN_DB" \
    --threads "$THREADS" \
    --report "$clean_contig_report/${sample}_clean_report.txt" \
    --output /dev/null \
    "$clean_assembly/contigs.fasta" >"$LOG_DIR/${sample}.clean_contigs_kraken2.log" 2>&1
}

run_step_4_mash_reference_check() {
  local sample="$1"
  local mash_dir="$OUTPUT_DIR/06_mash"
  local r1
  mkdir -p "$mash_dir"
  r1="$(reads_for_sample "$sample" R1)"

  if [[ ! -f "$REFERENCE" ]]; then
    echo "[Step 4/10] mash: $sample (reference missing; skipping)"
    return 0
  fi

  echo "[Step 4/10] mash: $sample"
  if [[ ! -f "$mash_dir/reference.msh" ]]; then
    mash sketch -o "$mash_dir/reference" "$REFERENCE" >"$LOG_DIR/reference.mash.log" 2>&1
  fi

  if [[ -f "$r1" ]]; then
    mash sketch -o "$mash_dir/${sample}" "$r1" >"$LOG_DIR/${sample}.mash_sketch.log" 2>&1
    mash dist "$mash_dir/reference.msh" "$mash_dir/${sample}.msh" > "$mash_dir/${sample}_to_reference.tsv" 2>"$LOG_DIR/${sample}.mash_dist.log"
  else
    echo "Trimmed reads not found for $sample; skipping mash comparison." >&2
  fi
}

run_step_5_spades_assembly() {
  local sample="$1"
  local asm_dir="$OUTPUT_DIR/05_assemblies/$sample"
  local r1="$OUTPUT_DIR/01_trimmed_reads/${sample}_R1.fastq.gz"
  local r2="$OUTPUT_DIR/01_trimmed_reads/${sample}_R2.fastq.gz"
  mkdir -p "$asm_dir"

  echo "[Step 5] SPAdes assembly: $sample"
  spades.py \
    -1 "$r1" \
    -2 "$r2" \
    -o "$asm_dir" \
    -t "$THREADS" \
    --isolate >"$LOG_DIR/${sample}.spades.log" 2>&1 || {
      echo "SPAdes failed for $sample. Check $LOG_DIR/${sample}.spades.log" >&2
      return 1
    }
}

run_step_6_quast() {
  local sample="$1"
  local asm_dir
  asm_dir="$(assembly_dir_for_sample "$sample")"
  local quast_dir="$OUTPUT_DIR/07_quast"
  mkdir -p "$quast_dir"

  echo "[Step 6/10] QUAST: $sample"
  if [[ -f "$asm_dir/contigs.fasta" ]]; then
    quast.py "$asm_dir/contigs.fasta" -o "$quast_dir/$sample" --threads "$THREADS" >"$LOG_DIR/${sample}.quast.log" 2>&1 || true
  else
    echo "No contigs.fasta for $sample; skipping QUAST." >&2
  fi
}

run_step_7_mlst() {
  local sample="$1"
  local asm_dir
  asm_dir="$(assembly_dir_for_sample "$sample")"
  local mlst_dir="$OUTPUT_DIR/08_mlst"
  mkdir -p "$mlst_dir"

  echo "[Step 7/10] MLST: $sample"
  if [[ -f "$asm_dir/contigs.fasta" ]]; then
    mlst "$asm_dir/contigs.fasta" > "$mlst_dir/${sample}.mlst.tsv" 2>"$LOG_DIR/${sample}.mlst.log" || true
  else
    echo "No contigs.fasta for $sample; skipping MLST." >&2
  fi
}

run_step_8_amrfinderplus() {
  local sample="$1"
  local asm_dir
  local amr_dir="$OUTPUT_DIR/09_amrfinderplus/$sample"
  asm_dir="$(assembly_dir_for_sample "$sample")"
  mkdir -p "$amr_dir"

  echo "[Step 8] AMRFinderPlus: $sample"
  if [[ -f "$asm_dir/contigs.fasta" ]]; then
    amrfinder \
      -n "$asm_dir/contigs.fasta" \
      -O Campylobacter \
      --plus \
      --mutation_all "$amr_dir/${sample}_all_mutations.tsv" \
      --output "$amr_dir/${sample}_amrfinder.tsv" \
      --threads "$THREADS" >"$LOG_DIR/${sample}.amrfinder.log" 2>&1
  else
    echo "No contigs.fasta for $sample; skipping AMRFinderPlus." >&2
  fi
}

run_step_8b_vfdb() {
  local sample="$1"
  local asm_dir
  local vfdb_dir="$OUTPUT_DIR/09_vfdb"
  asm_dir="$(assembly_dir_for_sample "$sample")"

  if [[ "$HAS_ABRICATE" != true ]]; then
    echo "ABRicate not found; skipping VFDB screen for $sample." >&2
    return 0
  fi
  if [[ ! -f "$asm_dir/contigs.fasta" ]]; then
    echo "No contigs.fasta for $sample; cannot run VFDB screen." >&2
    return 1
  fi
  mkdir -p "$vfdb_dir"
  echo "[Step 8b] ABRicate ($ABRICATE_DB): $sample"
  abricate --db "$ABRICATE_DB" --minid 80 --mincov 80 "$asm_dir/contigs.fasta" \
    >"$vfdb_dir/${sample}.tsv" 2>"$LOG_DIR/${sample}.abricate.log"
}

summarize_vfdb() {
  local vfdb_dir="$OUTPUT_DIR/09_vfdb"
  local report
  local -a sample_reports=()
  if [[ "$HAS_ABRICATE" != true ]]; then
    return 0
  fi
  for report in "$vfdb_dir"/*.tsv; do
    [[ -f "$report" ]] || continue
    if [[ "$(basename "$report")" != "vfdb_summary.tsv" ]]; then
      sample_reports+=("$report")
    fi
  done
  if [[ ${#sample_reports[@]} -eq 0 ]]; then
    return 0
  fi
  abricate --summary "${sample_reports[@]}" > "$vfdb_dir/vfdb_summary.tsv"
}

run_step_9_snippy() {
  local sample="$1"
  local r1
  local r2
  local snippy_dir="$OUTPUT_DIR/10_snippy"
  r1="$(reads_for_sample "$sample" R1)"
  r2="$(reads_for_sample "$sample" R2)"
  mkdir -p "$snippy_dir"

  if [[ ! -f "$REFERENCE" ]]; then
    echo "[Step 9/10] snippy: $sample (reference missing; skipping)"
    return 0
  fi

  echo "[Step 9] snippy: $sample"
  if [[ -f "$r1" && -f "$r2" ]]; then
    snippy \
      --outdir "$snippy_dir/$sample" \
      --ref "$REFERENCE" \
      --R1 "$r1" \
      --R2 "$r2" \
      --cpus "$THREADS" >"$LOG_DIR/${sample}.snippy.log" 2>&1
  else
    echo "Trimmed reads not found for $sample; skipping snippy." >&2
  fi
}

run_step_10_core_alignment_and_phylogeny() {
  local snippy_root="$OUTPUT_DIR/10_snippy"
  local phylo_dir="$OUTPUT_DIR/11_phylogeny"
  mkdir -p "$phylo_dir/core" "$phylo_dir/gubbins" "$OUTPUT_DIR/final_trees"

  if [[ ! -f "$REFERENCE" ]]; then
    echo "[Step 10] phylogeny: reference missing; skipping alignment and tree generation."
    return 0
  fi

  local -a snippy_samples=()
  for d in "$snippy_root"/*; do
    if [[ -d "$d" ]]; then
      snippy_samples+=("$d")
    fi
  done

  if [[ ${#snippy_samples[@]} -eq 0 ]]; then
    echo "[Step 10] no Snippy sample directories found; skipping final phylogeny."
    return 0
  fi

  echo "[Step 10] snippy-core + phylogeny"
  snippy-core --ref "$REFERENCE" --prefix "$phylo_dir/core/core" "${snippy_samples[@]}" >"$LOG_DIR/core_alignment.log" 2>&1

  if [[ -f "$phylo_dir/core/core.full.aln" ]]; then
    if command -v run_gubbins.py >/dev/null 2>&1; then
      mkdir -p "$phylo_dir/gubbins"
      run_gubbins.py --prefix "$phylo_dir/gubbins/campy" "$phylo_dir/core/core.full.aln" >"$LOG_DIR/gubbins.log" 2>&1
    else
      echo "run_gubbins.py not found; skipping recombination-aware tree inference." >&2
    fi

    if [[ "$HAS_IQTREE" == true && -f "$phylo_dir/core/core.aln" ]]; then
      mkdir -p "$phylo_dir/iqtree"
      "$IQTREE_BIN" \
        -s "$phylo_dir/core/core.aln" \
        -st DNA \
        -m MFP+ASC \
        -B 1000 \
        -alrt 1000 \
        -T AUTO \
        --prefix "$phylo_dir/iqtree/campy_final" >"$LOG_DIR/iqtree.log" 2>&1
    else
      echo "IQ-TREE unavailable or core SNP alignment missing; skipping additional ML tree." >&2
    fi

    find "$phylo_dir" -type f \( -iname "*.tree" -o -iname "*.tre" -o -iname "*.nwk" -o -iname "*.treefile" -o -iname "*.contree" -o -iname "RAxML_bestTree*" -o -iname "*gubbins*" \) -exec cp -f {} "$OUTPUT_DIR/final_trees/" \; 2>/dev/null || true

    local final_tree_found=false
    for tree in "$OUTPUT_DIR/final_trees"/*; do
      if [[ -f "$tree" ]]; then
        final_tree_found=true
        break
      fi
    done
    if [[ "$final_tree_found" == false ]]; then
      if [[ -f "$phylo_dir/core/core.full.aln" ]]; then
        echo "Final alignment saved at: $phylo_dir/core/core.full.aln"
      fi
    fi
  fi
}

main() {
  local -a fastq_pairs=()
  while IFS= read -r fastq; do
    if [[ -n "$fastq" ]]; then
      fastq_pairs+=("$fastq")
    fi
  done < <(find_fastq_pairs)

  if [[ ${#fastq_pairs[@]} -eq 0 ]]; then
    echo "No paired FASTQ files were found in $INPUT_DIR" >&2
    echo "Expected patterns like *_R1.fastq.gz or *_1.fastq.gz" >&2
    exit 1
  fi

  echo "Discovered ${#fastq_pairs[@]} paired FASTQ sample(s)."

  for r1 in "${fastq_pairs[@]}"; do
    r2="$(match_r2_for_r1 "$r1")"
    if [[ -z "$r2" || ! -f "$r2" ]]; then
      echo "Skipping $r1 because a matching read pair could not be found." >&2
      continue
    fi

    sample="$(sample_name_from_r1 "$r1")"

    if [[ -z "$sample" ]]; then
      echo "Could not infer sample name from $r1" >&2
      continue
    fi

    echo "========================================"
    echo "Processing sample: $sample"
    echo "R1: $r1"
    echo "R2: $r2"
    echo "========================================"

    run_step_1_fastp "$sample" "$r1" "$r2"
    run_step_2_fastqc_multiqc "$sample"
    run_step_3_kraken2 "$sample"
    run_step_5_spades_assembly "$sample"
    run_step_3b_rescue_cj006 "$sample"
    run_step_4_mash_reference_check "$sample"
    run_step_6_quast "$sample"
    run_step_7_mlst "$sample"
    run_step_8_amrfinderplus "$sample"
    run_step_8b_vfdb "$sample"
    run_step_9_snippy "$sample"
  done

  summarize_vfdb
  run_step_10_core_alignment_and_phylogeny

  echo
  echo "Pipeline complete. Results are in: $OUTPUT_DIR"
  echo "Final phylogenetic outputs are in: $OUTPUT_DIR/final_trees"
}

main "$@"
