[← Back to main page](index.md)

# Pipeline Overview

The SeqNado pipeline is built on Snakemake and handles the end-to-end processing of various sequencing assays. This page details the standard workflow and assay-specific steps.

## Usage

Run the pipeline for a given assay (e.g., ATAC-seq) via slurm using singularity containers using 16 cores:

```bash
seqnado pipeline atac --cores 16 --preset ss
```

### Presets

Presets control where and how jobs are executed. Pick the one that matches your setup:

| Preset | Profile                | Description |
|--------|------------------------|-------------|
| `le`   | `local_environment`    | Local execution, no containers (default) |
| `ls`   | `local_singularity`    | Local execution with Apptainer/Singularity |
| `lc`   | `local_conda`          | Local execution with Conda + Apptainer |
| `ld`   | `local_docker`         | Local execution with Conda + Docker |
| `ss`   | `slurm_singularity`    | SLURM cluster with Apptainer |
| `t`    | Test                   | For testing and development |

!!! tip
    **Not sure which to pick?** Use `le` if you're running on your own machine or a login node. Use `ss` if you're on an HPC cluster with SLURM. See the [HPC Clusters](cluster_config.md) guide for cluster setup.

### Common Options

| Option | Short | Description |
|--------|-------|-------------|
| `--cores INTEGER` | `-c` | Number of CPU cores for Snakemake to use (default: 1) |
| `--dry-run` | `-n` | Show what would be executed without running |
| `--unlock` | | Unlock the working directory after a failed/interrupted run |
| `--rerun-incomplete` | | Re-run jobs left incomplete from a previous run |
| `--scale-resources FLOAT` | `-s` | Scale memory/time requests (default: 1.0) |
| `--queue TEXT` | `-q` | Slurm queue/partition for the `ss` preset (default: short) |

Any additional arguments are passed directly to Snakemake (e.g., `--printshellcmds`).

For the full CLI reference, see [seqnado pipeline](cli.md#cli-seqnado-pipeline).

You can also point to a custom Snakemake profile directory with `--profile` (or `--profiles`), which overrides `--preset`.

## General Workflow

Regardless of the assay type, all SeqNado runs follow these core stages:

1.  **Quality Control**: FastQC checks per-base quality, GC content, and duplication levels. Fastq Screen detects contamination from other organisms — useful for catching sample swaps or adapter contamination early.
2.  **Adapter Trimming**: Trim Galore removes sequencing adapters and low-quality bases from read ends, which would otherwise cause misalignments or inflate duplicate rates.
3.  **Alignment**: Reads are mapped to the reference genome. DNA-based assays use Bowtie2 (fast, accurate for short reads). RNA-seq uses STAR, which is splice-aware and handles reads spanning exon-exon junctions.
4.  **Post-processing**: BAM files are filtered (remove unmapped reads, low-quality alignments), sorted, and indexed. Duplicates from PCR amplification are marked or removed so they don't inflate signal.
5.  **Signal Generation**: Normalised BigWig tracks are created for visualisation in genome browsers (IGV, UCSC). These let you visually inspect signal across the genome without loading the full BAM.
6.  **Summarization**: All QC metrics are aggregated into a single MultiQC HTML report for quick assessment of the entire run.

---

## Supported Assays

### ATAC-seq (ATAC)
Identifies regions of open chromatin by sequencing DNA accessible to the Tn5 transposase.

- **Filtering**: Mitochondrial reads are removed because Tn5 preferentially inserts into mitochondrial DNA, which would otherwise dominate the library. PCR duplicates are also removed. Reads overlapping ENCODE blacklist regions are discarded to avoid artefactual signal.
- **Tn5 Shift Correction**: Reads are shifted +4/-5 bp to account for the 9 bp staggered cut that Tn5 makes on insertion. This centres the signal on the actual insertion site rather than the fragment ends. Controlled by the `shift_for_tn5_insertion` config option.
- **Peak Calling**: MACS2 is the standard choice for most ATAC-seq experiments. LanceOtron uses a deep-learning approach and can be better at detecting broad or weak peaks. HOMER is also available.
- **Consensus Peaks**: When samples are grouped via a `consensus_group` column in the design file, BAMs are merged and peaks are re-called on the combined data to produce a robust consensus peak set.
- **QC**: TSS enrichment scores measure signal-to-noise — a high score (>7) confirms that open chromatin regions are well-captured. Fragment size distributions should show a clear nucleosomal banding pattern. FRiP (Fraction of Reads in Peaks) quantifies how much of your library falls within called peaks.

### ChIP-seq (ChIP)
Maps protein-DNA interactions by sequencing DNA fragments bound to a target protein.

- **Background Correction**: Input or IgG controls are essential to distinguish true binding from background noise. SeqNado pairs each IP sample with its control automatically from the design file.
- **Blacklist Filtering**: Reads in ENCODE blacklist regions are removed — these regions produce artefactually high signal regardless of the experiment.
- **Peak Calling**: Use MACS2 narrow mode for transcription factors (sharp, well-defined peaks) and broad mode for histone marks like H3K27me3 or H3K36me3 (wide, diffuse domains). SEACR is an alternative that works well for low-background experiments. HOMER is also available.
- **Consensus Peaks**: Merged sample groups produce consensus peak sets from the combined data, useful for defining a common set of regions across replicates.
- **Normalization**: Spike-in normalization (using a reference genome like *Drosophila*) corrects for global differences in pull-down efficiency between samples, which is critical when comparing conditions where total binding levels may change. Multiple normalization methods are available (Orlando, input-based, DESeq2, edgeR).
- **QC**: FRiP scores are calculated to measure the proportion of reads falling within peaks — a key indicator of enrichment quality.

### RNA-seq (RNA)
Quantifies gene expression across the transcriptome.

- **Alignment**: STAR is used because it handles spliced alignments — RNA-seq reads frequently span exon junctions, which a standard DNA aligner would fail to map.
- **Quantification**: featureCounts assigns aligned reads to genes using the genome annotation (GTF), producing a count matrix of genes-by-samples. Salmon is also available as an alternative that uses pseudoalignment for faster quantification without a separate alignment step.
- **Strand-Specific BigWigs**: For stranded libraries, separate plus- and minus-strand BigWig tracks are generated so you can visualise sense and antisense transcription independently.
- **Analysis**: DESeq2 performs differential expression analysis between conditions defined in your design file. It models count data with appropriate statistics and corrects for multiple testing.
- **QC**: Qualimap provides RNA-seq-specific metrics including gene body coverage, 5'/3' bias, and the proportion of reads mapping to exonic, intronic, and intergenic regions.

### CUT&Tag (CAT)
A low-input alternative to ChIP-seq that uses protein A-Tn5 fusion to tag chromatin at sites of antibody binding.

- **Tn5 Shift Correction**: Like ATAC-seq, reads are shifted +4/-5 bp to correct for the Tn5 insertion offset, centering signal on the true binding site. Controlled by the `shift_for_tn5_insertion` config option.
- **Filtering**: Filtering is tuned for the characteristically small fragment sizes generated by the targeted Tn5 insertion. CUT&Tag produces very low background. Blacklist regions are also removed.
- **Peak Calling**: SEACR is the recommended peak caller for CUT&Tag data — it was designed for the sparse, low-background signal that CUT&Tag produces, unlike MACS2 which assumes higher background levels. MACS2 and LanceOtron are also available if needed.
- **Consensus Peaks**: Merged sample groups produce consensus peak sets from the combined data.
- **Normalization**: Spike-in normalization is supported, as with ChIP-seq.

### SNP Analysis (SNP)
Identifies single nucleotide variants from sequencing data.

- **Alignment/Processing**: Standard Bowtie2 alignment followed by BAM post-processing with duplicate marking, which is especially important for variant calling to avoid false positives from PCR duplicates.
- **Variant Calling**: bcftools calls variants and produces VCF/BCF files containing identified variants, their quality scores, and genotype information. Multiallelic sites are split into separate records for downstream analysis.

### Methylation (METH)
Measures DNA methylation at single-base resolution using bisulfite sequencing.

- **Processing**: Bisulfite treatment converts unmethylated cytosines to uracil (read as thymine), so alignment requires a specialised approach that accounts for C-to-T conversions rather than treating them as mismatches. TAPS (TET-assisted pyridine borane sequencing) is also supported as an alternative chemistry — methylation values are automatically inverted to produce the correct output.
- **Spike-in BAM Splitting**: When spike-in genomes are configured (e.g., lambda DNA, pUC19), aligned reads are split by genome based on chromosome prefix. Reference and spike-in reads are then processed independently through all downstream steps.
- **Bias Correction**: MethylDackel calculates the M-bias plot for each genome (reference and spike-in separately), identifying positions with systematic bias at read ends. These biased positions are excluded during extraction to avoid skewing methylation estimates.
- **Conversion Rate QC**: Conversion rates are calculated from the M-bias data for each sample and genome, and aggregated into a summary table (`methylation_conversion.tsv`) and visualisation (`methylation_conversion.png`). For spike-in DNA (which is unmethylated), conversion should be >95% — a low spike-in conversion rate indicates incomplete bisulfite treatment and unreliable methylation calls. Read 1 and Read 2 are reported separately to detect strand-specific bias.
- **Methylation Calls**: MethylDackel extracts per-cytosine methylation levels from the aligned reads, producing CpG bedGraph files for each sample-genome combination. This gives you both the biological methylation state (from the reference genome) and the technical control (from the spike-in).

### MCC (MCC)
Micro Capture-C maps chromatin interactions at high resolution using targeted capture of specific viewpoints.

The MCC pipeline follows a multi-phase workflow:

**Phase 1: Viewpoint Preparation** 

- **Viewpoint FASTA Generation**: Viewpoint coordinates from the BED file are extracted as FASTA sequences from the reference genome using bedtools. These serve as alignment targets for identifying which capture probe pulled down each read.
- **Exclusion Regions**: Buffer zones around each viewpoint are defined (configurable via `exclusion_zone`) to exclude self-ligation and undigested fragments from the analysis.

**Phase 2: Read-to-Viewpoint Assignment** 

- **Viewpoint Alignment**: Trimmed reads (merged with FLASH for overlapping pairs) are aligned to the viewpoint FASTA using minimap2 with short-read settings (`-k 8 -w 1`). This identifies which capture probe each read originated from.
- **Read Splitting**: Reads that aligned to a viewpoint are extracted and converted back to FASTQ for re-alignment to the full genome.

**Phase 3: Genome Alignment** 

- **Primary Genome Alignment**: Viewpoint-assigned reads are aligned to the reference genome with Bowtie2 to determine the genomic location of each captured fragment.
- **Sensitive Re-alignment**: Reads that failed to map in the primary alignment are re-aligned with `--very-sensitive-local` settings to recover additional mappings. Both sets are then merged.

**Phase 4: Per-Replicate Processing** 

- **Query Name Sorting**: BAMs are sorted by read name to group paired reads for junction identification.
- **Viewpoint Annotation**: Each read pair is annotated with a viewpoint (VP) tag indicating which capture probe it originated from.
- **Deduplication**: PCR duplicates are removed using viewpoint-aware deduplication to avoid inflating contact frequencies.
- **Ligation Junction Extraction**: Read pairs representing capture-C ligation junctions are identified for each viewpoint. These are output as pairs files (chromosome, position for each end of the contact).
- **Ligation Statistics**: Per-sample cis and trans contact counts are extracted for normalisation and QC. A high cis ratio (contacts on the same chromosome as the viewpoint) indicates good capture efficiency.

**Phase 5: Group Merging and Aggregation** 

- **BAM Merging**: Replicate BAMs within each consensus group are merged for increased statistical power.
- **Grouped Junction Extraction**: Ligation junctions are re-extracted from the merged BAMs for each viewpoint, providing group-level contact data.

**Phase 6: Contact Matrix Generation** 

- **Cooler Creation**: Pairs files are loaded into Cooler format (HDF5-based) at the primary resolution defined in the config.
- **Multi-resolution Zoomification**: Each Cooler is zoomified to create `.mcool` files with multiple resolution levels for browsing at different scales.
- **Cooler Aggregation**: Per-viewpoint Cooler files are combined into a single group-level `.mcool` file.

**Phase 7: Signal Tracks** 

- **Per-Replicate BigWigs**: BigWig tracks are generated for each sample and viewpoint, normalised by the number of cis contacts (n_cis CPM) to allow comparison between samples with different capture efficiencies.
- **Group BigWigs**: Both normalised (n_cis-scaled) and raw (unscaled) BigWigs are generated for each consensus group.
- **Aggregated BigWigs**: Replicate BigWigs are averaged within each condition to produce mean signal tracks.
- **Comparison BigWigs**: Subtraction BigWigs are generated between conditions (e.g., treatment minus control) to highlight differential interactions.

**Phase 8: Peak Calling** 

- **LanceOtron-MCC**: A deep-learning peak caller identifies significant interaction peaks from the unscaled group BigWigs. This runs on GPU when available.

### CRISPR Screens (CRISPR)
Quantifies guide RNA representation across pooled CRISPR screen experiments.

- **Adapter Detection**: Automatically detects CRISPR-specific adapter sequences in the reads to ensure correct trimming.
- **Quantification**: Counts the number of reads mapping to each guide RNA in the library, producing a guide-level count matrix across samples.
- **Analysis**: MAGeCK performs statistical analysis to identify guides and genes that are enriched or depleted between conditions (e.g., treatment vs. control). Both the test module (rank-based) and MLE module (maximum likelihood) are available, accounting for the multiple guides per gene.

### Multiomics (MULTIOMICS)
Run multiple assay types together in one project for integrated outputs. This is useful when you have matched samples across assays (e.g., ATAC + RNA from the same conditions) and want to analyse them in a single coordinated run. Configure per-assay sections and enable multiomics mode in the config.

---

## Downstream Analysis

### Motif Analysis
Identifies DNA sequence motifs enriched in peak regions, helping to determine which transcription factors may be driving binding or accessibility.

- **MEME-ChIP**: Performs both de novo motif discovery and comparison against known motif databases in a single analysis.
- **HOMER**: Runs `findMotifsGenome` for de novo and known motif enrichment, with broader database support.

### Heatmaps and Metaplots
Deeptools generates signal heatmaps and average profile plots (metaplots) centred on peak regions or gene features. These provide a visual summary of signal distribution across all peaks or genes in a single figure.

### Visualisation
- **UCSC Genome Browser Hub**: When enabled, SeqNado generates a track hub that can be loaded directly into the UCSC Genome Browser for interactive exploration of BigWig tracks and peak calls.
- **Plotnado**: An interactive browser visualisation for exploring signal at specific loci.

### Merged Quantification
When consensus peaks are defined from merged sample groups, reads from individual samples are counted across the consensus peak set to produce a peaks-by-samples count matrix — analogous to the gene count matrix in RNA-seq.

---

## Technical Details

### Resource Management
SeqNado automatically calculates required cores and memory for each step based on your provided configuration and the available system resources.

### Parallelization
The pipeline leverages Snakemake's ability to run samples in parallel, scaling from local machines to large high-performance computing (HPC) clusters.

For command-line options (presets, queues, scaling), see [seqnado pipeline](cli.md#cli-seqnado-pipeline).

---

**See Also:**

- [Outputs Reference](outputs.md) - Understanding your results
- [HPC Clusters](cluster_config.md) - Configure for HPC environments
- [Troubleshooting](troubleshooting.md#pipeline-execution-seqnado-pipeline) - Pipeline execution issues
- [Output Examples](examples.md) - Example analysis results
