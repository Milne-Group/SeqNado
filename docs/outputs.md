[← Back to main page](index.md)

# Pipeline Outputs

All SeqNado analysis results are organized within the `seqnado_output/` directory (or your custom output directory specified during configuration). This page describes the comprehensive structure and types of files you can expect from your pipeline runs.

## 📁 General Output Structure

```
seqnado_output/
├── {assay}/                       # Assay-specific directory (atac, chip, rna, etc.)
│   ├── seqnado_report.html       # 🎯 Main interactive QC report
│   ├── aligned/                   # BAM alignment files
│   ├── bigwigs/                   # BigWig coverage tracks
│   ├── peaks/                     # Peak calling results (for applicable assays)
│   ├── qc/                        # Quality control metrics
│   ├── hub/                       # UCSC Genome Browser hub
│   ├── genome_browser_plots/      # PlotNado visualizations
│   ├── geo_submission/            # GEO submission-ready files
│   ├── trimmed/                   # Adapter-trimmed FastQ files
│   └── logs/                      # Process execution logs
```

## 🎯 Main Entry Point: SeqNado Report

The **`seqnado_report.html`** file is your primary analysis dashboard, providing:

- **Interactive QC Summary**: FastQC, alignment statistics, and quality metrics
- **Sample Overview**: All samples displayed with key metrics
- **Peak Statistics**: Number of peaks called, genomic distribution
- **Fragment Size Distributions**: For ATAC-seq and ChIP-seq analyses
- **Library Complexity**: Duplication rates and unique read counts
- **Multi-Sample Comparisons**: Side-by-side QC metrics

!!! tip "Viewing the Report"
    Open `seqnado_report.html` in any modern web browser. No server required!

## 📊 Core Output Files

### Alignment Files (`aligned/`)

Processed BAM files at various stages:

```
aligned/
├── {sample}.bam                   # Final processed BAM
├── {sample}.bam.bai              # BAM index
├── raw/                           # Initial alignments
├── sorted/                        # Coordinate-sorted BAMs
├── filtered/                      # Quality-filtered reads
├── duplicates_removed/            # PCR duplicate-removed reads
├── blacklist_regions_removed/     # Blacklist-filtered reads
└── shifted_for_tn5_insertion/     # Tn5-corrected (ATAC-seq only)
```

**File Formats:**
- **BAM**: Binary alignment format, viewable with samtools/IGV
- **BAI**: Index files for rapid random access

**Example sizes** (from test data):
- Input control: ~66 KB
- ChIP sample: ~3.9 MB

### Coverage Tracks (`bigwigs/`)

Genome-wide signal tracks for visualization:

```
bigwigs/
├── bamnado/                       # Bamnado-generated tracks
│   ├── {sample}_CPM.bw           # CPM normalized
│   ├── {sample}_RPGC.bw          # RPGC normalized
│   └── {sample}_spike_in.bw      # Spike-in normalized (if applicable)
├── deeptools/                     # DeepTools coverage
│   └── {sample}.bw
└── homer/                         # HOMER tag directories
    └── {sample}_homer.ucsc.bedGraph.gz
```

**Normalization Methods:**
- **CPM** (Counts Per Million): Standard library size normalization
- **RPGC** (Reads Per Genomic Content): Accounts for effective genome size
- **Spike-in**: Uses external control DNA for absolute quantification

### Peak Calls (`peaks/`)

Peak calling results from multiple callers:

```
peaks/
├── macs2/
│   ├── {sample}_peaks.narrowPeak  # Narrow peaks (TF/histone marks)
│   ├── {sample}_peaks.xls         # Detailed peak information
│   ├── {sample}_summits.bed       # Peak summit positions
│   └── {sample}_model.r           # Peak calling model
├── macs3/                          # MACS3 results (if enabled)
├── homer/                          # HOMER peak calls
│   └── {sample}_peaks.txt
├── lanceotron/                     # ML-based peak calling
│   └── {sample}_L-tron.bed
└── consensus/                      # Merged peaks across replicates
    └── consensus_peaks.bed
```

**Peak File Formats:**
- **narrowPeak**: ENCODE standard peak format (BED6+4)
- **broadPeak**: For histone marks with broad enrichment
- **BED**: Simple genomic coordinates

### Quality Control (`qc/`)

Comprehensive QC metrics:

```
qc/
├── fastqc_raw/                    # Pre-trimming FastQC reports
│   └── {sample}_fastqc.html
├── fastqc_trimmed/                # Post-trimming FastQC
├── fastq_screen/                  # Contamination screening
│   └── {sample}_screen.html
├── qualimap_bamqc/                # BAM quality metrics
│   ├── {sample}/
│   │   ├── qualimapReport.html
│   │   └── images_qualimapReport/
│   │       ├── genome_coverage_histogram.png
│   │       ├── genome_insert_size_histogram.png
│   │       └── genome_gc_content_per_window.png
└── preseq/                        # Library complexity estimates
    └── {sample}_complexity.txt
```

**QC Metrics Include:**
- Read quality scores per base position
- GC content distribution
- Adapter content
- Duplication rates
- Insert size distributions
- Mapping statistics
- Coverage uniformity

### UCSC Genome Browser Hub (`hub/`)

Ready-to-load UCSC track hub:

```
hub/
├── seqnado_hub.hub.txt           # Hub description file
├── seqnado_hub.genomes.txt       # Genome assemblies
├── tracknado_config.json         # TrackNado configuration
└── {assay}/
    ├── trackDb.txt                # Track definitions
    ├── seqnado_report.html        # QC report copy
    └── *.bigWig                   # Processed signal tracks
```

**Usage:**
1. Upload the `hub/` directory to a web-accessible location
2. Load in UCSC Genome Browser using the hub URL
3. Or use locally with IGV/other genome browsers

### Genome Browser Plots (`genome_browser_plots/`)

High-quality publication-ready visualizations generated with PlotNado:

```
genome_browser_plots/
├── {region}_{sample_comparison}.pdf
├── {region}_{sample_comparison}.png
└── regions.bed                    # Plotted genomic regions
```

## 🧬 Assay-Specific Outputs

### ATAC-seq

**Additional files:**
```
atac/
├── peaks/
│   ├── lanceotron/                # ML peak calling optimized for ATAC
│   └── consensus_peaks.bed        # Reproducible peaks across replicates
├── tss_enrichment/                # TSS enrichment scores
│   └── {sample}_tss_enrichment.txt
└── fragment_analysis/             # Fragment size analysis
    └── {sample}_fragments.txt
```

**Key Metrics:**
- TSS enrichment score (>7 for high-quality ATAC)
- Fragment size distribution (nucleosome periodicity)
- FRiP (Fraction of Reads in Peaks) score

### ChIP-seq / CUT&Tag

**Additional files:**
```
chip/
├── peaks/
│   ├── macs2/                     # Standard MACS2 peaks
│   ├── homer/                     # HOMER findPeaks results
│   └── lanceotron/                # ML-based peak calling
├── tag_dirs/                      # HOMER tag directories
│   └── {sample}/
└── motif_analysis/                # Motif enrichment (if enabled)
    └── {sample}_motifs/
```

**Spike-in Normalization** (if applicable):
- Normalized BigWigs account for IP efficiency
- Scaling factors stored in log files

### RNA-seq

**Additional files:**
```
rna/
├── aligned/
│   └── {sample}_Aligned.sortedByCoord.out.bam
├── counts/
│   ├── {sample}_ReadsPerGene.out.tab
│   └── combined_counts.tsv       # All samples merged
├── differential_expression/       # DESeq2 results
│   ├── deseq2_results.csv
│   ├── MA_plot.pdf
│   ├── volcano_plot.pdf
│   └── PCA_plot.pdf
└── splice_junctions/
    └── {sample}_SJ.out.tab
```

**Count Tables:**
- Gene-level quantification
- Transcript-level counts (if enabled)
- TPM/FPKM normalized values

### Methylation (METH)

```
meth/
├── methylation_calls/
│   ├── {sample}_CpG_report.txt
│   ├── {sample}.bedGraph
│   └── {sample}.bismark.cov
├── mbias/                         # M-bias plots
└── splitting_report.txt           # Methylation extraction summary
```

### CRISPR Screens

```
crispr/
├── counts/
│   ├── guide_counts.txt           # Raw guide counts
│   └── normalized_counts.txt      # Library size-normalized
├── mageck/                        # MAGeCK analysis results
│   ├── gene_summary.txt
│   ├── sgrna_summary.txt
│   └── QC_plots/
└── screen_report.html             # CRISPR-specific QC
```

## 📈 Example Output Sizes

Based on test datasets:

| Assay Type | Sample | Final BAM | Peaks | BigWig | Total Output |
|------------|--------|-----------|-------|--------|--------------|
| ATAC-seq   | Single | ~500 MB   | 50K   | 42 MB  | ~1.2 GB      |
| ChIP-seq   | + Input| 3.9 MB    | 15K   | 593 KB | ~150 MB      |
| RNA-seq    | Single | ~1.5 GB   | N/A   | 85 MB  | ~2.5 GB      |

!!! note "Actual sizes vary"
    Output sizes depend on sequencing depth, genome size, and enabled analyses.

## 🔄 Accessing Your Results

### Command Line

```bash
# Navigate to output directory
cd seqnado_output/

# View main report
firefox {assay}/seqnado_report.html &

# List peaks
ls -lh {assay}/peaks/macs2/

# Load BAM in IGV
igv {assay}/aligned/{sample}.bam
```

### Opening Reports

The HTML reports can be opened directly in your browser:

```bash
# On local machine
open seqnado_output/chip/seqnado_report.html

# Via X11 forwarding on HPC
firefox seqnado_output/chip/seqnado_report.html &

# Transfer to local machine
scp -r user@hpc:path/to/seqnado_output/ ./
```

## 📤 GEO Submission (`geo_submission/`)

Pre-formatted files for GEO/SRA submission:

```
geo_submission/
├── metadata.xlsx                  # Sample metadata spreadsheet
├── processed_files/               # Symlinks to key processed files
│   ├── bigwigs/
│   └── peaks/
└── raw_fastq/                     # Links to original FastQ files
```

## 🔍 Finding Specific Outputs

### Peak calling results
```bash
find seqnado_output/ -name "*_peaks.narrowPeak"
```

### Coverage tracks for visualization
```bash
find seqnado_output/ -name "*.bw" -o -name "*.bigWig"
```

### QC HTML reports
```bash
find seqnado_output/ -name "*qc.html" -o -name "seqnado_report.html"
```

## 💡 Next Steps

After reviewing your outputs:

1. **QC Assessment**: Check `seqnado_report.html` for sample quality
2. **Peak Analysis**: Explore peaks in `peaks/` directories
3. **Visualization**: Load BigWigs in genome browsers
4. **Downstream Analysis**: Use count tables or peaks for further analysis
5. **Publication**: Use plots from `genome_browser_plots/` and QC metrics

For pipeline rerunning or parameter adjustments, see: [seqnado pipeline](cli.md#cli-seqnado-pipeline)
