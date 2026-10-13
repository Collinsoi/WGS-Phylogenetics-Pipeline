# Corrected Draft: Campylobacter WGS Bioinformatics Pipeline

This is the cleaned, working workflow extracted from the successful terminal history in the raw notes, excluding failed package installs and setup errors.

## Study design and dataset
- 15 Campylobacter jejuni isolates: CJ001–CJ015
- Paired-end Illumina WGS reads in `data/raw/`
- Reference genome: `references/species/C_jejuni_NCTC11168.fasta`

---

## Step 1: Raw read quality control and trimming
Goal: remove adapters and low-quality bases before downstream analyses.

Tools used:
- `fastp`
- `fastqc`
- `multiqc`

Representative command:
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

Outputs:
- trimmed paired reads in `data/trimmed/`
- per-sample `fastp` HTML/JSON reports
- consolidated QC report at `results/multiqc/multiqc_report.html`

Interpretation:
- reads passed quality filtering with adapter trimming
- low-quality reads and short reads were removed before further processing

---

## Step 2: Contamination screening
Goal: check whether each isolate had obvious contamination or taxonomy mismatch.

Tool used:
- `kraken2`

Representative workflow:
```bash
mkdir -p results/kraken2_batch
./scripts/04_batch_kraken2.sh
```

Outputs:
- per-sample Kraken reports in `results/kraken2_batch/`
- summary table in `results/kraken2_batch/kraken2_summary.csv`

Purpose:
- identify likely non-target or contaminated samples before assembly and phylogeny

---

## Step 3: Reference similarity / strain grouping check
Goal: assess whether samples were close to the expected C. jejuni reference and whether one isolate was a good central representative for phylogeny.

Tools used:
- `mash`

Representative workflow:
```bash
mash sketch -o results/phylogeny/reference_check/NCTC11168 references/species/C_jejuni_NCTC11168.fasta
```

Then for each assembly:
```bash
mash dist results/phylogeny/reference_check/NCTC11168.msh "$fasta"
```

And for the combined pairwise comparison:
```bash
mash sketch -l results/phylogeny/reference_check/assemblies.txt -o results/phylogeny/reference_check/all15
mash dist results/phylogeny/reference_check/all15.msh results/phylogeny/reference_check/all15.msh > results/phylogeny/reference_check/all15_pairwise_mash.tsv
```

Outputs:
- `results/phylogeny/reference_check/nctc11168_distances.tsv`
- `results/phylogeny/reference_check/reference_candidates.tsv`

Key finding from the log:
- best central candidate was `CJ015`
- all isolates were genetically close to C. jejuni NCTC11168, supporting use of this reference for SNP calling

---

## Step 4: De novo assembly
Goal: construct contig assemblies from filtered reads.

Tool used:
- `spades.py`

Representative command:
```bash
spades.py \
  -1 data/trimmed/CJ001_R1.fastq.gz \
  -2 data/trimmed/CJ001_R2.fastq.gz \
  -o results/assembly_batch/CJ001 \
  -t 4 \
  --isolate
```

Outputs:
- `results/assembly_batch/CJxxx/contigs.fasta`
- `results/assembly_batch/CJxxx/scaffolds.fasta`
- logs under each isolate directory

This was performed batch-wise across the 15 isolates.

---

## Step 5: Assembly evaluation and typing
Goal: assess assembly quality and assign sequence types / AMR profile.

Tools used:
- `quast.py`
- `mlst`
- `amrfinder`

Representative commands:
```bash
quast.py results/assembly_batch/CJ001/contigs.fasta results/assembly_batch/CJ002/contigs.fasta ... -o results/quast_batch
```

```bash
./scripts/07_batch_mlst.sh
```

Outputs:
- `results/quast_batch/report.tsv`
- `results/quast_batch/report.html`
- `results/mlst_batch/mlst_results.tsv`

This step provided assembly metrics and sequence type information for each isolate.

---

## Step 6: Reference-based variant calling with Snippy
Goal: compare each isolate to the C. jejuni reference and identify SNPs/indels relative to the reference genome.

Tool used:
- `snippy`

Representative workflow:
```bash
./scripts/10_batch_snippy.sh
```

Snippy performed for each isolate:
- read mapping with `bwa`
- sorting/indexing with `samtools`
- variant calling with `freebayes`
- filtering with `bcftools` and `vt`
- generation of per-sample SNP output, including:
  - `snps.vcf`
  - `snps.aligned.fa`
  - `snps.consensus.fa`
  - `snps.bam`
  - `snps.tab`

The successful logs show that all 15 samples were processed successfully.

---

## Step 7: Core genome alignment
Goal: combine all isolate SNP calls into a single whole-genome core alignment for phylogenetic inference.

Tool used:
- `snippy-core`

Representative command:
```bash
snippy-core \
  --ref references/species/C_jejuni_NCTC11168.fasta \
  --prefix results/phylogeny/core/core \
  results/snippy/CJ001 \
  results/snippy/CJ002 \
  ... \
  results/snippy/CJ015
```

Outputs:
- `results/phylogeny/core/core.full.aln`
- `results/phylogeny/core/core.aln`
- `results/phylogeny/core/core.txt`
- `results/phylogeny/core/core.vcf`

This produced a single alignment containing all isolates plus the reference sequence.

---

## Step 8: Recombination-aware phylogeny
Goal: generate a robust phylogenetic tree while accounting for recombination in the core SNP alignment.

Tools used:
- `Gubbins`
- `RAxML`

Command used:
```bash
run_gubbins.py --prefix results/phylogeny/gubbins/campy results/phylogeny/core/core.full.aln
```

This step:
- built a phylogenetic tree from the SNP alignment
- detected recombinant regions
- iteratively refined the alignment/tree to remove misleading phylogenetic signal

The logs show multiple RAxML iterations and completion of the Gubbins analysis.

---

## Final workflow summary
The successful end-to-end pipeline was:

1. `fastp` quality filtering and adapter trimming
2. `fastqc` and `multiqc` QC inspection
3. `kraken2` contamination screening
4. `mash` reference similarity check
5. `SPAdes` de novo assembly
6. `quast` assembly assessment
7. `mlst` typing
8. `snippy` reference-based SNP calling
9. `snippy-core` core SNP alignment
10. `Gubbins` + `RAxML` recombination-aware phylogenetic inference

This is the corrected draft version of the working pipeline and represents the actual successful workflow recorded in the raw terminal history.
