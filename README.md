# Campylobacter WGS Bioinformatics Pipeline

This repository documents a whole-genome sequencing workflow for 15
*Campylobacter jejuni* isolates (CJ001–CJ015), from paired-end Illumina reads
through assembly, typing, AMR screening, variant calling, and phylogenetic
analysis.

## Study design and inputs

- Samples: CJ001–CJ015
- Sample-to-run mapping: `metadata/sample_metadata.csv`
- Raw paired reads: `data/raw/`
- Reference genome: `references/species/C_jejuni_NCTC11168.fasta`

## Conda environments

The main analysis environment is defined in `environment.yml`. Create and
activate it with:

```bash
conda env create -f environment.yml
conda activate campy-wgs
```

It includes SRA Toolkit, fastp, FastQC, MultiQC, SeqKit, Kraken2/KrakenTools,
Mash, SPAdes, QUAST, MLST, AMRFinderPlus, Biopython, Matplotlib, and SciPy.
QUAST is installed from PyPI in this environment because the native Apple
Silicon Conda package did not solve in the package check for this workflow.

On Apple Silicon macOS, Snippy, Gubbins, and IQ-TREE were run in dedicated
`osx-64` Conda environments under Rosetta because of package compatibility
constraints. The recorded setup was:

```bash
CONDA_SUBDIR=osx-64 conda create -n snippy-env \
  -c conda-forge -c bioconda snippy=4.6.0 -y
CONDA_SUBDIR=osx-64 conda create -n gubbins-env \
  -c conda-forge -c bioconda gubbins -y
CONDA_SUBDIR=osx-64 conda create -n iqtree-env \
  -c conda-forge -c bioconda iqtree -y
```

If installing more packages into one of these environments, scope the
platform setting to that environment:

```bash
conda activate snippy-env
conda config --env --set subdir osx-64
```

ABRicate is an optional VFDB screening tool and is not in the main
environment. Its dependencies did not solve natively in the checked Apple
Silicon Conda setup, but did solve for `osx-64`. Install it separately and
initialize its local databases:

```bash
CONDA_SUBDIR=osx-64 conda create -n abricate-env \
  -c conda-forge -c bioconda abricate -y
conda activate abricate-env
abricate --setupdb
```

## Individual scripts

The runnable stage scripts are in `scripts/`; shared paths, metadata loading,
and sample-selection helpers are in `scripts/common.sh`. They are intended to
be run from the repository root. The default thread count is 4; set
`THREADS` to override it. `METADATA`, `RAW_DIR`, `TRIMMED_DIR`,
`ASSEMBLY_ROOT`, and `REFERENCE` can also be set to override paths.

| Script | Purpose |
|---|---|
| `scripts/02_download_remaining.sh` | Download mapped SRA runs and create paired FASTQ files |
| `scripts/03_batch_fastp_qc.sh` | Trim reads with fastp; run FastQC and MultiQC |
| `scripts/04_batch_kraken2.sh` | Classify trimmed read pairs with Kraken2 and write a summary |
| `scripts/04b_rescue_CJ006.sh` | Screen, genus-filter, reassemble, and re-screen the CJ006 rescue |
| `scripts/05_batch_mash.sh` | Compare assemblies against the reference and one another |
| `scripts/06_batch_spades.sh` | Assemble trimmed paired reads |
| `scripts/06b_batch_quast.sh` | Evaluate assemblies with QUAST |
| `scripts/07_batch_mlst.sh` | Assign MLST profiles |
| `scripts/08_batch_amrfinder.sh` | Run AMRFinderPlus, including all-mutations output |
| `scripts/09_batch_vfdb.sh` | Optional ABRicate screen against VFDB |
| `scripts/10_batch_snippy.sh` | Run reference-based Snippy per sample |
| `scripts/11_batch_phylogeny.sh` | Create Snippy core alignments |
| `scripts/12_batch_gubbins.sh` | Run recombination analysis with Gubbins |
| `scripts/13_batch_iqtree.sh` | Filter to variable sites and run IQ-TREE |

### Run the complete workflow

After installing the Conda environments above, the orchestrator runs the
stages in workflow order and invokes each script in its required environment.
Provide the Kraken2 database path and, if needed, a reference FASTA path:

```bash
KRAKEN_DB=/path/to/kraken2_database \
  ./run_pipeline.sh /path/to/C_jejuni_NCTC11168.fasta
```

The reference argument can be omitted when the reference is at the default
path shown above. To reuse existing raw FASTQ files, set `SKIP_DOWNLOAD=1`.
To include the optional VFDB screen, set `RUN_VFDB=1`; this requires
`abricate-env` with its VFDB database initialized. `THREADS` defaults to 4.
Environment names can be overridden with `MAIN_ENV`, `SNIPPY_ENV`,
`GUBBINS_ENV`, `IQTREE_ENV`, and `ABRICATE_ENV`. The runner requires
`KRAKEN_DB` because it performs both the routine Kraken2 screen and CJ006
rescue. It writes per-stage logs to `results/logs/` and stops on the first
failed stage.

To print the planned stages without executing them, use
`DRY_RUN=1 ./run_pipeline.sh /path/to/reference.fasta`. This still validates
that the reference file exists.

### Run stages individually

Individual scripts can also be run manually in the Conda environment
appropriate to each stage. For the main environment, run download, read QC,
Kraken2 screening, SPAdes assembly, CJ006 rescue, Mash, QUAST, MLST, and
AMRFinderPlus in the sequence listed by the orchestrator. Run `09_batch_vfdb.sh`
from `abricate-env` when desired. Run `10_batch_snippy.sh` and
`11_batch_phylogeny.sh` in `snippy-env`, followed by `12_batch_gubbins.sh` in
`gubbins-env` and `13_batch_iqtree.sh` in `iqtree-env`.

The CJ006 rescue requires KrakenTools' `extract_kraken_reads.py` and
`KRAKEN_DB`. Downstream scripts select rescued CJ006 reads and assembly when
those outputs exist.

## Workflow

### 1. Read quality control and trimming

fastp trims adapters and filters low-quality or short reads; FastQC and
MultiQC provide per-sample and combined QC reports. The batch script writes
trimmed paired reads under `data/trimmed/`, fastp reports under
`results/fastp/`, FastQC output under `results/fastqc/`, and the combined
report under `results/multiqc/`.

Example fastp settings used:

```bash
fastp \
  -i data/raw/SRR23704524_1.fastq.gz \
  -I data/raw/SRR23704524_2.fastq.gz \
  -o data/trimmed/CJ001_R1.fastq.gz \
  -O data/trimmed/CJ001_R2.fastq.gz \
  --detect_adapter_for_pe \
  --qualified_quality_phred 20 \
  --length_required 50 \
  --html results/fastp/CJ001/fastp.html \
  --json results/fastp/CJ001/fastp.json \
  --thread 4
```

Run `./scripts/03_batch_fastp_qc.sh` for all samples.

### 2. Contamination screening and CJ006 troubleshooting

Kraken2 read-level classification is run with
`KRAKEN_DB=/path/to/database ./scripts/04_batch_kraken2.sh`. Per-sample reports
and classification files are written to `results/kraken2_batch/`, together
with `kraken2_summary.csv`.

#### CJ006 contamination discovery and rescue

The initial CJ006 assembly was 3.57 Mbp, unusually large for the expected
Campylobacter genome. Kraken2 classified 69.30% of its contigs, prompting
read-level screening and a targeted reassembly.

Kraken2 processed 976,004 paired read records (287.85 Mbp); 884,318 were
classified (90.61%). KrakenTools `extract_kraken_reads.py` was used with taxid
`194` (*Campylobacter*) and `--include-children`, retaining 867,149 read pairs.
These reads were reassembled with SPAdes into
`results/assembly_batch/CJ006_clean/`. The resulting assembly was 1.74 Mbp.

The cleaned assembly still had contigs classified to multiple
Campylobacter species, including *C. jejuni*, *C. coli*, and *C. lari*.
Because the extraction retained the whole genus, this step reduces
non-Campylobacter sequence but does not establish that CJ006 is a pure
*C. jejuni* isolate. Review the post-rescue taxonomic report before treating
it as a confirmed single-species sample.

Run the scripted rescue after initial assembly:

```bash
KRAKEN_DB=/path/to/kraken2_database ./scripts/04b_rescue_CJ006.sh
```

The script writes read and assembly classifications under
`results/kraken2_assembly/`, retained paired reads under
`results/decontamination/CJ006/`, and the rescued assembly under
`results/assembly_batch/CJ006_clean/`.

### 3. Reference similarity and strain grouping

Mash sketches the reference and assemblies, compares each assembly to
`C_jejuni_NCTC11168.fasta`, and produces pairwise distances among isolates:

```bash
./scripts/05_batch_mash.sh
```

Outputs are in `results/phylogeny/reference_check/`. The recorded
reference-similarity check found CJ015 to be the best central candidate, and
the isolates were reported as genetically close to the NCTC11168 reference.

### 4. De novo assembly

SPAdes runs in isolate mode on paired trimmed reads:

```bash
./scripts/06_batch_spades.sh
```

Typical outputs are `contigs.fasta`, `scaffolds.fasta`, and SPAdes logs under
`results/assembly_batch/CJxxx/`. When CJ006 rescued reads are present, the
rescue script instead writes its assembly to
`results/assembly_batch/CJ006_clean/`.

### 5. Assembly evaluation, typing, and AMR

QUAST evaluates assembly metrics; MLST assigns sequence types:

```bash
./scripts/06b_batch_quast.sh
./scripts/07_batch_mlst.sh
```

Run AMRFinderPlus with the options used in the successful recorded analysis:

```bash
./scripts/08_batch_amrfinder.sh
```

The per-sample command uses `-O Campylobacter`, `--plus`, `--mutation_all`,
and the long-form `--output` option. CJ006 uses the rescued assembly if it is
available. Outputs are written under `results/amrfinder_batch/CJxxx/`;
assembly metrics are in `results/quast_batch/` and MLST results in
`results/mlst_batch/mlst_results.tsv`.

### 6. Virulence screening against VFDB

The terminal history does not document completed VFDB results. This is an
additional recommended analysis, not a claim that virulence screening was
performed. AMRFinderPlus AMR profiling is not a substitute for a dedicated
VFDB search.

With ABRicate installed and its VFDB database initialized, run:

```bash
./scripts/09_batch_vfdb.sh
```

The script uses example minimum identity and coverage thresholds of 80%.
Per-sample candidate matches and the combined summary are written under
`results/virulence/vfdb/`. Record the ABRicate and database versions and the
thresholds with results. Database hits alone do not establish gene
expression, functionality, or virulence phenotype.

### 7. Reference-based variant calling

Snippy maps reads to the reference and calls/filter variants. Run:

```bash
conda activate snippy-env
./scripts/10_batch_snippy.sh
```

The script writes per-sample results under `results/snippy/CJxxx/`, including
Snippy VCF, alignment, consensus, BAM, and tabular outputs. CJ006 uses the
retained genus-filtered reads when its rescue reads are present.

### 8. Core alignment and recombination-aware phylogeny

With the Snippy environment active, create the core alignment:

```bash
./scripts/11_batch_phylogeny.sh
```

Snippy-core writes `results/phylogeny/core/core.full.aln`,
`core.aln`, `core.txt`, and `core.vcf`. The script only creates these core
alignments; Gubbins and IQ-TREE are separate stages for use in their dedicated
environments.

Run Gubbins in its dedicated environment, or use
`scripts/12_batch_gubbins.sh` to process the full core alignment and collect
its tree outputs:

```bash
conda activate gubbins-env
./scripts/12_batch_gubbins.sh
```

The additional IQ-TREE maximum-likelihood analysis was separate from
Gubbins/RAxML. An initial `+ASC` attempt on `campy_filtered_15.fasta` failed
because the alignment contained 94 invariant sites. IQ-TREE generated
`campy.varsites.phy`; the successful run used that variable-sites alignment.
The reproducible script below extracts variable nucleotide columns from
Snippy-core's `core.aln` into `campy.varsites.fasta` before running the same
model search and support tests:

```bash
conda activate iqtree-env
./scripts/13_batch_iqtree.sh
```

This used ModelFinder with ascertainment-bias correction (`+ASC`), 1,000
ultrafast bootstrap replicates (`-B 1000`), and 1,000 SH-aLRT replicates
(`-alrt 1000`). The selected model by BIC was `TVM+F+ASC+R2`. The principal
outputs were `campy_final.treefile`, `campy_final.contree`, and
`campy_final.iqtree` under `results/phylogeny/iqtree/`.

## Workflow summary

1. Download paired reads from SRA using the sample metadata.
2. Trim and QC reads with fastp, FastQC, and MultiQC.
3. Screen read pairs with Kraken2; investigate the CJ006 contamination signal.
4. Check reference similarity and pairwise distances with Mash.
5. Assemble reads using SPAdes; genus-filter and reassemble CJ006 as needed.
6. Evaluate assemblies using QUAST and assign types with MLST.
7. Profile AMR determinants and mutations with AMRFinderPlus.
8. Optionally screen candidate virulence genes against VFDB with ABRicate.
9. Call reference-based variants with Snippy and create a core alignment.
10. Run recombination-aware Gubbins/RAxML inference.
11. Run the separate IQ-TREE maximum-likelihood analysis with `+ASC` and
    branch support tests.

The Gubbins/RAxML and IQ-TREE analyses are distinct. CJ006 genus-level
filtering is a troubleshooting rescue, not proof of single-species purity;
VFDB screening is recommended but was not recorded as completed.
