[← Back to main page](index.md)

# Pipeline Overview

The SeqNado pipeline is built on Snakemake and handles the end-to-end processing of various sequencing assays. This page details the standard workflow and assay-specific steps.

## Usage

Run the pipeline for a given assay (e.g., ATAC-seq):

```bash
seqnado pipeline atac --preset le
```

For all arguments and presets, see the CLI reference: [seqnado pipeline](cli.md#cli-seqnado-pipeline).

## General Workflow

Regardless of the assay type, all SeqNado runs follow these core stages:

1.  **Quality Control**: Raw FASTQ validation using FastQC and Fastq Screen.
2.  **Adapter Trimming**: Automatic detection and removal of adapters via Trim Galore.
3.  **Alignment**: Mapping reads to the reference genome (Bowtie2 for DNA, STAR for RNA).
4.  **Post-processing**: Filtering, sorting, and indexing of BAM files.
5.  **Signal Generation**: BigWig track creation for visualization.
6.  **Summarization**: Integration of all QC metrics into a final MultiQC report.

---

## Supported Assays

### ATAC-seq (ATAC)
The ATAC-seq pipeline focuses on identifying regions of open chromatin.

- **Filtering**: Removal of mitochondrial reads and duplicates.
- **Peak Calling**: Uses MACS2 or LanceOtron for accessibility peak detection.
- **QC**: Calculates TSS enrichment scores and fragment size distributions.

### ChIP-seq (ChIP)
Designed for Protein-DNA interaction mapping.

- **Background Correction**: Support for Input/IgG controls.
- **Peak Calling**: MACS2 (narrow or broad) and SEACR options.
- **Normalization**: Supports spike-in normalization if specified in the config.

### RNA-seq (RNA)
Standard transcriptomic profiling.

- **Alignment**: Splice-aware mapping using STAR.
- **Quantification**: Gene-level counts via `featureCounts`.
- **Analysis**: Optional automated DESeq2 differential expression results.

### CUT&Tag (CAT)
Low-input chromatin profiling targeting protein-DNA interactions.

- **Specialized Filtering**: Optimized for small fragment sizes.
- **Callers**: SEACR integration for sparse signal.

### SNP Analysis (SNP)
Variant detection workflows.

- **Alignment/Processing**: Standard mapping and BAM post-processing.
- **Variant Calling**: Uses established variant calling utilities; outputs VCF/BCF.

### Methylation (METH)
Bisulfite sequencing and methylation calling.

- **Processing**: Bisulfite-aware read handling.
- **Methylation Calls**: Generates cytosine methylation metrics and summaries.

### MCC (MCC)
Multiplex chromatin conformation capture.

- **Processing**: Contact map-oriented post-processing.
- **Analysis**: MCC-specific peak/enrichment outputs and QC metrics.

### CRISPR Screens (CRISPR)
Guide-level quantification and screen analysis.

- **Quantification**: Count guides/targets across samples.
- **Analysis**: Screen statistics and summary tables.

### Multiomics (MULTIOMICS)
Run multiple assay types together in one project for integrated outputs. Configure per-assay sections and enable multiomics mode in the config.

---

## Technical Details

### Resource Management
SeqNado automatically calculates required cores and memory for each step based on your provided configuration and the available system resources.

### Parallelization
The pipeline leverages Snakemake's ability to run samples in parallel, scaling from local machines to large high-performance computing (HPC) clusters.

## Next Steps
Once your pipeline is running, you can monitor the progress in the terminal. After completion, visit the [Outputs](outputs.md) page to understand the result structure.

For command-line options (presets, queues, scaling), see [seqnado pipeline](cli.md#cli-seqnado-pipeline).
