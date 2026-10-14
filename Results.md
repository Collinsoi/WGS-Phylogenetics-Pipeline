
## Key findings

This project analyzed a 15-isolate training subset of poultry-associated
*Campylobacter jejuni* genomes collected between 2008 and 2016. The objective
was to examine whether genomic lineages showed corresponding patterns in
antimicrobial-resistance determinants and virulence-gene content. Because this
is a deliberately selected training subset rather than a population-representative
sample, the results should be interpreted as comparative genomic observations,
not prevalence estimates.

### 1. Taxonomic screening and assembly quality

Kraken2 read-level screening showed that all 15 datasets were dominated by
*Campylobacter*-associated reads, although the proportion of unclassified reads
varied among isolates. Most de novo assemblies showed characteristics consistent
with *C. jejuni*, with genome sizes around the expected ~1.6–1.8 Mb range and
GC contents close to 30%.

CJ006 was a clear exception. Its original assembly was ~2.97 Mb with 43.91% GC,
1,402 contigs, and an N50 of 91,771 bp. Assembly-level Kraken2 classification
showed that 61.25% of contigs were assigned to *Bifidobacterium pullorum*,
whereas only 1.51% were assigned to *Campylobacter*.

Read-level screening of CJ006 nevertheless showed that 88.85% of reads were
assigned to the genus *Campylobacter*, while approximately 1.06% were assigned
to *B. pullorum*.

After retaining *Campylobacter*-classified reads and reassembling the sample,
the CJ006 assembly improved to approximately 1.74 Mb, 30.12% GC, 29 contigs,
N50 130,036 bp, and L50 4. Post-rescue Kraken2 classification was dominated by
*Campylobacter* and no longer showed the strong *Bifidobacterium* signal seen
in the original assembly.

This demonstrates the importance of combining read-level and assembly-level
taxonomic QC rather than relying on assembly statistics alone.

### 2. Sequence-type diversity

MLST identified 13 distinct sequence types among the 15 isolates, indicating
substantial lineage diversity within the small dataset.

ST-2116 was the only sequence type represented by more than one isolate:

- CJ008 — ST-2116
- CJ011 — ST-2116
- CJ015 — ST-2116

These three isolates shared the same seven-locus MLST allele profile and later
formed a strongly supported phylogenetic cluster, providing an internal test of
whether lineage relationships were reflected in AMR and virulence profiles.

### 3. Antimicrobial-resistance determinants

AMRFinderPlus identified both acquired resistance genes and resistance-associated
mutations across the dataset.

The quinolone-resistance-associated `gyrA_T86I` substitution was detected in
12/15 isolates, while CJ009 carried `gyrA_T86V`. The tetracycline-resistance
gene `tet(O)` was detected in 9/15 isolates.

Beta-lactamase-related determinants were also widespread and included variants
such as:

- `blaOXA-193`
- `blaOXA-184`
- `blaOXA-630`
- other `blaOXA` calls

CJ004 additionally carried the macrolide-associated `23S_A2075G` mutation.

Other detected determinants included `arsP`, `acr3`, `cmeB`, and
`rplV_A103V`; these were interpreted separately from the major antibiotic-
resistance markers because some represent metal/stress tolerance or candidate
resistance-associated variation rather than direct phenotypic resistance.

The three ST-2116 isolates showed substantial AMR-profile similarity. CJ008,
CJ011, and CJ015 all carried `gyrA_T86I` and `blaOXA-193`; however, CJ008 and
CJ011 also carried `tet(O)`, whereas CJ015 did not. This illustrates that
closely related isolates can share a common resistance background while still
differing in accessory resistance determinants.

### 4. Virulence-gene repertoire

ABRicate screening against VFDB at ≥80% identity and ≥80% coverage identified
128 distinct virulence-associated genes across the 15 isolates.

Of these:

- 83 genes were present in all 15 isolates
- 45 genes were variable
- 27 genes occurred in three or fewer isolates

This indicates a large conserved virulence-gene backbone together with a
smaller accessory component that varied among genomic lineages.

Selected variable genes included `kpsC`, `porA`, `gmhA2`, `hddA`, `hddC`,
`flaA`, `cstIII`, `neuA1`, `neuB1`, `neuC1`, and `wlaN`.

The LOS/sialylation-associated group `cstIII`, `neuA1`, `neuB1`, `neuC1`,
and `wlaN` occurred together in CJ001, CJ005, and CJ007, suggesting coordinated
presence of a lineage-associated locus rather than independent gene occurrence.

Because ABRicate calls depend on sequence identity and coverage thresholds,
absence of highly variable surface-associated loci such as `flaA` should be
interpreted cautiously and does not necessarily prove biological absence.

### 5. Reference suitability and core-genome coverage

Mash comparisons supported *C. jejuni* NCTC11168 as a reasonable common
reference for SNP analysis. Mash distances between study isolates and the
reference were small overall, while CJ015 was the most central isolate in
pairwise Mash comparisons within the 15-isolate dataset.

Snippy-core mapping showed substantial reference-aligned sequence across all
isolates. CJ004 had the lowest aligned fraction at approximately 89.9%, while
CJ007 had the highest at approximately 97.1%.

Variant counts ranged from 4,722 in CJ007 to 18,186 in CJ002, showing
considerable genomic diversity despite all isolates belonging to the same
species.

### 6. Recombination-aware phylogenetic analysis

Because homologous recombination can distort bacterial SNP phylogenies, the
Snippy core alignment was processed with Gubbins before final tree inference.

Branch-specific Gubbins statistics showed substantial heterogeneity in
recombination across the tree, with some branches carrying considerably more
recombinant sequence than others. The recombination-filtered polymorphic-site
alignment was therefore used for the final phylogenetic analysis.

After removal of the reference and invariant positions, the final IQ-TREE
alignment contained:

- 15 isolates
- 41,260 variable nucleotide sites
- 25,135 parsimony-informative sites
- 0 constant sites

The best-fit substitution model selected by ModelFinder was:

`TVM+F+ASC+R2`

Branch support was evaluated using 1,000 ultrafast bootstrap replicates and
1,000 SH-aLRT replicates.

### 7. Strong phylogenetic clusters

Several relationships were strongly supported in the final unrooted tree.

The clearest included:

- CJ010 + CJ014: approximately 100/100 SH-aLRT/UFBoot support
- CJ009 + CJ012: approximately 98.2/97
- CJ011 + CJ015: approximately 87.9/91
- CJ008 + (CJ011, CJ015): approximately 100/100

The broader relationship connecting CJ004 to the CJ010/CJ014 group had much
weaker support (~47.8/55) and was therefore not interpreted as a robust lineage
relationship.

### 8. Concordance between lineage, AMR, and virulence

The clearest integrated signal was observed for the three ST-2116 isolates,
CJ008, CJ011, and CJ015.

These isolates:

- shared the same MLST type and allele profile
- formed a strongly supported core-SNP phylogenetic cluster
- shared major AMR determinants
- displayed highly similar virulence-gene profiles

Whole-profile virulence similarity was approximately:

- CJ008 vs CJ011: 96.8%
- CJ008 vs CJ015: 97.9%
- CJ011 vs CJ015: 98.9%

Similarity remained high even when only variable virulence genes were
considered, although accessory differences were still detectable.

Together, these observations support the working hypothesis that genetically
related isolates tend to share more similar AMR and virulence profiles.
However, the differences in accessory genes such as `tet(O)` also show that
closely related isolates are not genomically identical.

### 9. Overall interpretation

The 15-isolate dataset contained substantial genomic diversity, with 13 MLST
sequence types, heterogeneous AMR determinants, and both conserved and variable
virulence loci.

Despite this diversity, closely related isolates—most clearly the ST-2116
CJ008/CJ011/CJ015 group—showed concordant phylogenetic, AMR, and virulence
patterns. This suggests that lineage background contributes to the distribution
of clinically relevant genomic traits, while accessory gene gain/loss and
recombination introduce additional within-lineage variation.

The CJ006 contamination episode additionally demonstrated why bacterial WGS
analysis requires active quality control throughout the workflow: apparently
minor read-level contamination can produce major assembly-level distortion and
substantially affect downstream genomic interpretation.

### Limitations

This analysis used a selected 15-isolate subset from a larger public dataset and
was developed primarily as a comparative-genomics and workflow-training
project. The sample is therefore not statistically representative of the full
population.

AMR findings are genotype-based and should not be interpreted as measured
phenotypic susceptibility. Likewise, VFDB matches indicate sequence similarity
to known virulence-associated genes but do not demonstrate gene expression or
virulence phenotype.

Phylogenetic relationships are based on a reference-mapped,
recombination-filtered core-SNP analysis and should be interpreted as genomic
relatedness within this dataset rather than direct evidence of transmission,
ancestry, or epidemiological linkage.
