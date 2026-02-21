# Output Examples

This page showcases example outputs from SeqNado pipelines to help you understand what to expect from your analyses. All examples are based on real test data processed through the pipeline.

## SeqNado QC Report

The main `seqnado_report.html` provides a MultiQC report with comprehensive quality control metrics.

### Report Sections

#### 1. General Statistics

View high-level sample information and key metrics:

- **Total reads**: Raw sequencing depth
- **Mapped reads**: Alignment success rate
- **Duplication rate**: PCR duplicate percentage
- **GC content**: Library GC distribution
- **Insert size**: Fragment size metrics (PE data)

#### 2. FastQC Results

Quality metrics on raw reads:

- **Per base sequence quality**: Quality scores across read positions
- **Per sequence quality scores**: Overall read quality distribution
- **Per base sequence content**: Nucleotide balance
- **Sequence duplication levels**: Library complexity indicators
- **Adapter content**: Contamination levels

#### 3. Alignment Metrics

Mapping statistics from Bowtie2 (DNA assays) or STAR (RNA-seq):

- **Alignment rates**: Percentage uniquely mapped, multimapped, unmapped
- **Paired-end concordance**: Proper pair percentages
- **Insert size distribution**: Fragment size histograms

#### 4. Peak Calling Summary (ChIP/ATAC/CUT&Tag)

Peak detection results:

- **Number of peaks**: Total peaks called per caller
- **FRiP scores**: Fraction of Reads in Peaks (if enabled)
- **Peak caller comparison**: Overlap between MACS2, HOMER, LanceOtron, SEACR

#### 5. Library Complexity

Picard duplicate metrics:

- **Unique reads**: Non-duplicate read counts
- **Duplication rates**: PCR duplicate percentages
- **Library complexity estimates**: From Picard MarkDuplicates metrics

## ChIP-seq Example Output

### Directory Structure

```
seqnado_output/chip/
├── seqnado_report.html                    # Main QC dashboard
├── protocol.txt                           # Data processing protocol
├── aligned/
│   ├── chip-rx_MLL.bam                    # Final processed BAM
│   ├── chip-rx_MLL.bam.bai
│   ├── chip-rx_input.bam
│   └── chip-rx_input.bam.bai
├── bigwigs/
│   ├── bamnado/
│   │   ├── chip-rx_MLL.bigWig
│   │   └── chip-rx_input.bigWig
│   ├── deeptools/
│   │   ├── chip-rx_MLL.bigWig
│   │   └── chip-rx_input.bigWig
│   └── homer/
│       ├── chip-rx_MLL.bigWig
│       └── chip-rx_input.bigWig
├── peaks/
│   ├── macs2/
│   │   └── chip-rx_MLL.bed               # Simplified 3-column BED
│   ├── homer/
│   │   └── chip-rx_MLL.bed
│   └── lanceotron/
│       └── chip-rx_MLL.bed
├── qc/
│   ├── fastqc_raw/
│   │   ├── chip-rx_MLL_1_fastqc.html
│   │   ├── chip-rx_MLL_2_fastqc.html
│   │   ├── chip-rx_input_1_fastqc.html
│   │   └── chip-rx_input_2_fastqc.html
│   ├── fastq_screen/                      # If enabled
│   │   ├── chip-rx_MLL_1_screen.html
│   │   └── chip-rx_input_1_screen.html
│   ├── qualimap_bamqc/
│   │   ├── chip-rx_MLL/qualimapReport.html
│   │   └── chip-rx_input/qualimapReport.html
│   ├── alignment_stats.tsv
│   └── library_complexity/
│       ├── chip-rx_MLL.metrics
│       └── chip-rx_input.metrics
├── hub/
│   └── seqnado_hub.hub.txt
└── tag_dirs/
    └── chip-rx_MLL/
```

### Example Peak File Content

SeqNado outputs simplified 3-column BED files for all peak callers:

**BED format** (`chip-rx_MLL.bed`):

```
chr1    3054728    3055228
chr1    3669834    3670334
chr2    5847291    5847791
...
```

**Columns:**

1. Chromosome
2. Start position
3. End position

## ATAC-seq Example Output

### Directory Structure

```
seqnado_output/atac/
├── seqnado_report.html
├── protocol.txt
├── aligned/
│   ├── atac_sample.bam                    # Tn5-shifted, filtered BAM
│   └── atac_sample.bam.bai
├── bigwigs/
│   ├── bamnado/
│   │   └── atac_sample.bigWig
│   ├── deeptools/
│   │   └── atac_sample.bigWig
│   └── homer/
│       └── atac_sample.bigWig
├── peaks/
│   └── lanceotron/                        # Default peak caller for ATAC
│       └── atac_sample.bed
├── qc/
│   ├── fastqc_raw/
│   ├── qualimap_bamqc/
│   │   └── atac_sample/
│   │       └── qualimapReport.html
│   ├── alignment_stats.tsv
│   └── library_complexity/
│       └── atac_sample.metrics
└── hub/
    └── seqnado_hub.hub.txt
```

!!! note
    All intermediate BAM processing stages (sorting, blacklist removal, duplicate removal, Tn5 shifting, filtering) are temporary and automatically deleted. Only the final `aligned/{sample}.bam` is retained.

### ATAC-seq Quality Indicators

Key metrics to check in the MultiQC report:

| Metric | Good Quality | What to Look For |
|--------|-------------|------------------|
| Nucleosome periodicity | Clear peaks at ~200bp intervals | Visible in insert size distribution |
| Mitochondrial % | <10% | Low mitochondrial read contamination |
| FRiP score | >30% | High fraction of reads in peaks |
| Unique reads | >80% | Good library complexity |

### Fragment Size Distribution

ATAC-seq shows characteristic nucleosome-free (~150bp) and mono-nucleosome (~200bp) peaks, visible in the insert size distribution within the MultiQC report.

## Qualimap BAM QC Report

The Qualimap reports provide detailed alignment quality metrics. For RNA-seq, `qualimap_rnaseq` is used instead of `qualimap_bamqc`.

### Key Visualisations

**Coverage Histogram**

- Distribution of genome coverage depths
- Helps identify over/under-sequenced regions
- Shows sequencing uniformity

**Insert Size Distribution**

- Fragment size histogram for paired-end data
- Critical for ATAC-seq quality assessment
- Reveals nucleosome positioning

**GC Content Distribution**

- AT/GC bias detection
- Compares observed vs theoretical
- Identifies contamination or bias

**Mapping Quality**

- Distribution of MAPQ scores
- Higher scores = more confident alignments
- Helps assess multi-mapping issues

## Example FastQ Screen Results

FastQ Screen checks for contamination across reference genomes (when enabled via `run_fastq_screen`):

### Typical Clean Sample

```
Library: chip-rx_MLL_1
Genome          %Mapping    %One_hit    %Multi_hit    Status
--------------------------------------------------------------
Human (hg38)    98.2%       85.4%       12.8%         OK
Mouse (mm10)     0.8%        0.5%        0.3%         OK
E. coli          0.0%        0.0%        0.0%         OK
Adapters         0.3%        0.3%        0.0%         OK
PhiX             0.0%        0.0%        0.0%         OK
```

### Concerning Sample (Contamination)

```
Library: sample_contaminated
Genome          %Mapping    %One_hit    %Multi_hit    Status
--------------------------------------------------------------
Human (hg38)    65.2%       58.4%        6.8%         WARNING
Mouse (mm10)    32.8%       29.5%        3.3%         WARNING
E. coli          1.2%        1.1%        0.1%         WARNING
```

## HOMER Tag Directory

HOMER creates tag directories for downstream analysis:

```
tag_dirs/chip-rx_MLL/
├── tagInfo.txt                     # Read statistics
├── tagLengthDistribution.txt       # Fragment sizes
├── tagCountDistribution.txt        # Tag depth per position
├── freqDistribution.txt            # Frequency statistics
├── chr1.tags.tsv                   # Per-chromosome tags
├── chr2.tags.tsv
└── ...
```

## BigWig Coverage Tracks

BigWig files provide genome-wide signal visualisation.

### File Naming Convention

BigWig files are organised by tool, scaling method, and individual vs merged:

```
Individual samples:
- bigwigs/deeptools/chip-rx_MLL.bigWig              # Unscaled
- bigwigs/bamnado/chip-rx_MLL.bigWig
- bigwigs/homer/chip-rx_MLL.bigWig
- bigwigs/deeptools/csaw/chip-rx_MLL.bigWig         # CSAW-normalised
- bigwigs/deeptools/spikein/orlando/chip-rx_MLL.bigWig  # Spike-in normalised

Merged consensus groups:
- bigwigs/deeptools/merged/consensus_group.bigWig             # Unscaled merged
- bigwigs/deeptools/merged/csaw/consensus_group.bigWig        # CSAW-scaled merged
- bigwigs/deeptools/merged/spikein/orlando/consensus_group.bigWig  # Spike-in scaled merged
```

See [Pipeline Outputs — Normalisation factor calculation](outputs.md#normalisation-factor-calculation) for a full explanation of how per-sample and merged scale factors are derived.

For RNA-seq, stranded tracks include `_plus` and `_minus` suffixes:

```
- bigwigs/deeptools/rna_sample_plus.bigWig
- bigwigs/deeptools/rna_sample_minus.bigWig
```

### Loading in UCSC Genome Browser

1. Upload BigWig files to a web-accessible location
2. Or use the auto-generated hub at `hub/seqnado_hub.hub.txt`
3. Tracks display sample signal across genome
4. Compare multiple samples side-by-side

## GEO Submission Files

Ready-to-submit files for GEO/SRA (when enabled):

```
geo_submission/
├── samples_table.txt                       # Sample metadata (TSV format)
├── md5sums.txt                             # Combined checksums
├── raw_data_checksums.txt                  # Checksums for raw FASTQs
├── processed_data_checksums.txt            # Checksums for processed files
├── upload_instructions.txt                 # GEO upload instructions
├── chip-rx_MLL_1.fastq.gz                  # Symlinks to raw FASTQ R1
├── chip-rx_MLL_2.fastq.gz                  # Symlinks to raw FASTQ R2
├── chip-rx_MLL_deeptools_unscaled.bigWig   # Renamed processed files
├── chip-rx_MLL_macs2.bed                   
└── chip/                                   # Upload directory
```

Files are flattened from the nested directory structure into a single directory with descriptive filenames that encode the tool and scaling method.

## Genome Browser Plots (PlotNado)

Publication-ready visualisations of genomic regions (when plotting coordinates are configured):

### Output Files

```
track_plots/
├── MYC_promoter.svg             # Named region from BED file
├── chr1-1000000-1005000.svg     # Unnamed region uses coordinates
└── template.toml                # PlotNado configuration template
```

Plot filenames are derived from the Name column in the input BED file, or from `{chr}-{start}-{end}` if no name is provided. Output format can be `svg`, `png`, or `pdf` as configured.

## Tips for Exploring Outputs

### Quick Quality Check

```bash
# Check main report
firefox seqnado_output/chip/seqnado_report.html

# Count peaks called
wc -l seqnado_output/chip/peaks/macs2/*.bed

# View alignment stats
samtools flagstat seqnado_output/chip/aligned/chip-rx_MLL.bam

# Check bigwig file
bigWigInfo seqnado_output/chip/bigwigs/deeptools/chip-rx_MLL.bigWig
```

### Finding Specific Results

```bash
# All HTML reports
find seqnado_output/ -name "*.html"

# All peak files
find seqnado_output/ -name "*.bed"

# All coverage tracks
find seqnado_output/ -name "*.bigWig"
```

## Understanding File Formats

### BAM Files
- **Binary alignment format** (`.bam` extension)
- Stores aligned sequencing reads
- Includes alignment quality, CIGAR strings, and flags

### BigWig Files
- **Binary coverage track format** (`.bigWig` extension)
- Efficient genome browser visualisation
- Contains normalised signal values

### BED Files
- **Tab-delimited genomic coordinates** (`.bed` extension)
- SeqNado peak outputs use 3-column BED: chr, start, end
- Standard BED can include additional columns (name, score, strand)

### FastQ Files
- **Raw sequencing reads** (`.fastq.gz` extension)
- Four lines per read: header, sequence, +, quality scores
- Usually gzip compressed (.gz)

---

For more information on interpreting these outputs for your specific experiment, consult the [Pipeline Overview](pipeline.md).
