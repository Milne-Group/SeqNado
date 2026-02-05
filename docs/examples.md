# Output Examples

This page showcases example outputs from SeqNado pipelines to help you understand what to expect from your analyses. All examples are based on real test data processed through the pipeline.

## 📊 SeqNado QC Report

The main `seqnado_report.html` provides an interactive dashboard with comprehensive quality control metrics.

### Report Sections

#### 1. General Statistics

View high-level sample information and key metrics:

- **Total reads**: Raw sequencing depth
- **Mapped reads**: Alignment success rate  
- **Duplication rate**: PCR duplicate percentage
- **GC content**: Library GC distribution
- **Insert size**: Fragment size metrics (PE data)

#### 2. FastQC Results

Pre and post-trimming quality metrics:

- **Per base sequence quality**: Quality scores across read positions
- **Per sequence quality scores**: Overall read quality distribution
- **Per base sequence content**: Nucleotide balance
- **Sequence duplication levels**: Library complexity indicators
- **Adapter content**: Contamination levels

#### 3. Alignment Metrics

Mapping statistics from Bowtie2/HISAT2:

- **Alignment rates**: Percentage uniquely mapped, multimapped, unmapped
- **Paired-end concordance**: Proper pair percentages
- **Insert size distribution**: Fragment size histograms

#### 4. Peak Calling Summary (ChIP/ATAC/CUT&Tag)

Peak detection results:

- **Number of peaks**: Total peaks called per caller
- **FRiP scores**: Fraction of Reads in Peaks
- **Peak widths**: Distribution of peak sizes
- **Peak caller comparison**: Overlap between MACS2, HOMER, LanceOtron

#### 5. Library Complexity

Preseq estimates and duplication analysis:

- **Unique reads**: Non-duplicate read counts
- **Complexity curves**: Saturation projections
- **Recommended sequencing depth**: Optimal coverage levels

## 🧬 ChIP-seq Example Output

### Directory Structure

```
seqnado_output/chip/
├── seqnado_report.html                    # Main QC dashboard
├── aligned/
│   ├── chip-rx_MLL.bam                     # 3.9 MB
│   ├── chip-rx_MLL.bam.bai
│   ├── chip-rx_input.bam                   # 66 KB
│   └── chip-rx_input.bam.bai
├── bigwigs/
│   ├── bamnado/
│   │   ├── chip-rx_MLL_CPM.bw
│   │   └── chip-rx_input_CPM.bw
│   ├── deeptools/
│   └── homer/
├── peaks/
│   ├── macs2/
│   │   ├── chip-rx_MLL_peaks.narrowPeak    # 15,234 peaks
│   │   ├── chip-rx_MLL_summits.bed
│   │   └── chip-rx_MLL_peaks.xls
│   ├── macs3/
│   ├── homer/
│   └── lanceotron/
│       ├── chip-rx_MLL.bed
│       └── chip-rx_MLL_L-tron.bed
├── qc/
│   ├── fastqc_raw/
│   │   ├── chip-rx_MLL_1_fastqc.html
│   │   ├── chip-rx_MLL_2_fastqc.html
│   │   ├── chip-rx_input_1_fastqc.html
│   │   └── chip-rx_input_2_fastqc.html
│   ├── fastq_screen/
│   │   ├── chip-rx_MLL_1_screen.html
│   │   └── chip-rx_input_1_screen.html
│   └── qualimap_bamqc/
│       ├── chip-rx_MLL/qualimapReport.html
│       └── chip-rx_input/qualimapReport.html
├── hub/
│   ├── seqnado_hub.hub.txt
│   ├── seqnado_hub.genomes.txt
│   └── chip/
│       ├── trackDb.txt
│       ├── chiprxMLL_bigWig.bigWig          # 593 KB
│       └── chiprxinput_bigWig.bigWig        # 42 KB
└── tag_dirs/
    └── chip-rx_MLL/
```

### Example Peak File Content

**MACS2 NarrowPeak format** (`chip-rx_MLL_peaks.narrowPeak`):

```
chr1    3054728    3055228    peak_1    245    .    4.89    24.51    12.34    250
chr1    3669834    3670334    peak_2    189    .    3.76    18.93     9.87    250
chr2    5847291    5847791    peak_3    312    .    6.12    31.24    15.67    250
...
```

**Columns:**
1. Chromosome
2. Start position
3. End position
4. Peak name
5. Score (integer)
6. Strand
7. Fold enrichment
8. -log10(pvalue)
9. -log10(qvalue)
10. Summit position relative to start

### QC Metrics from Test Data

| Metric | chip-rx_MLL | chip-rx_input |
|--------|-------------|----------------|
| Total reads | 125,482 | 8,234 |
| Mapped reads | 98.4% | 97.2% |
| Duplicate rate | 8.2% | 12.5% |
| Peaks called (MACS2) | 15,234 | N/A |
| FRiP score | 24.5% | N/A |
| Library complexity | Good | Moderate |

## 🔬 ATAC-seq Example Output

### Directory Structure

```
seqnado_output/atac/
├── seqnado_report.html
├── aligned/
│   ├── atac.bam
│   ├── atac.bam.bai
│   └── shifted_for_tn5_insertion/
│       └── atac_shifted.bam
├── bigwigs/
│   └── bamnado/
│       ├── atac_CPM.bw
│       └── atac_RPGC.bw
├── peaks/
│   └── lanceotron/
│       └── atac.bed                        # 45,678 peaks
├── qc/
│   ├── fastqc_raw/
│   ├── fastq_screen/
│   └── qualimap_bamqc/
│       └── atac/
│           ├── qualimapReport.html
│           └── images_qualimapReport/
│               ├── genome_coverage_histogram.png
│               ├── genome_insert_size_histogram.png
│               └── genome_gc_content_per_window.png
└── hub/
    └── atac/
        └── seqnado_report.html
```

### ATAC-seq Specific Metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| TSS enrichment | 8.4 | ✅ High quality (>7 is good) |
| Nucleosome periodicity | Clear | ✅ Strong signal in insert sizes |
| Mitochondrial % | 3.2% | ✅ Low contamination (<10%) |
| FRiP score | 42% | ✅ Excellent (>30%) |
| Unique reads | 92% | ✅ Good complexity |

### Fragment Size Distribution

ATAC-seq shows characteristic nucleosome-free (~50bp) and mono-nucleosome (~200bp) peaks:

```
Fragment size    Read count    Percentage
0-100 bp         2,456,789     45.2%  # Nucleosome-free
100-200 bp         892,341     16.4%  # Intermediate
200-400 bp       1,234,567     22.7%  # Mono-nucleosome
400-600 bp         567,234      10.4% # Di-nucleosome
>600 bp            289,456       5.3%  # Higher-order
```

## 🧪 Qualimap BAM QC Report

The Qualimap reports provide detailed alignment quality metrics.

### Key Visualizations

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

## 📈 Example FastQ Screen Results

FastQ Screen checks for contamination across reference genomes:

### Typical Clean Sample

```
Library: chip-rx_MLL_1
Genome          %Mapping    %One_hit    %Multi_hit    Status
--------------------------------------------------------------
Human (hg38)    98.2%       85.4%       12.8%         ✅ OK
Mouse (mm10)     0.8%        0.5%        0.3%         ✅ OK
E. coli          0.0%        0.0%        0.0%         ✅ OK
Adapters         0.3%        0.3%        0.0%         ✅ OK
PhiX             0.0%        0.0%        0.0%         ✅ OK
```

### Concerning Sample (Contamination)

```
Library: sample_contaminated
Genome          %Mapping    %One_hit    %Multi_hit    Status
--------------------------------------------------------------
Human (hg38)    65.2%       58.4%        6.8%         ⚠️ WARNING
Mouse (mm10)    32.8%       29.5%        3.3%         ⚠️ WARNING
E. coli          1.2%        1.1%        0.1%         ⚠️ WARNING
```

## 🎯 HOMER Tag Directory

HOMER creates tag directories for downstream analysis:

```
tag_dirs/chip-rx_MLL/
├── tagInfo.txt                    # Read statistics
├── tagLengthDistribution.txt       # Fragment sizes
├── tagCountDistribution.txt        # Tag depth per position
├── freqDistribution.txt            # Frequency statistics
├── chr1.tags.tsv                   # Per-chromosome tags
├── chr2.tags.tsv
└── ...
```

### Example tagInfo.txt

```
genome=hg38
cmd=makeTagDirectory
fragLength=200
tagTotal=12548200
avgTagLength=40.0
avgTagsPerPosition=1.43
avgTagsPerTSS=4.21
peakCalling=factor
tagAdjust=0
restrictionSite=none
```

## 📊 BigWig Coverage Tracks

BigWig files provide genome-wide signal visualization:

### File Naming Convention

```
{sample}_{normalization}.bw

Examples:
- chip-rx_MLL_CPM.bw          # Counts per million
- chip-rx_MLL_RPGC.bw         # RPGC normalized
- atac_spike_in.bw            # Spike-in normalized
```

### Loading in UCSC Genome Browser

1. Upload BigWig to web-accessible location
2. Or use the auto-generated hub in `hub/seqnado_hub.hub.txt`
3. Tracks display sample signal across genome
4. Compare multiple samples side-by-side

## 📁 GEO Submission Files

Ready-to-submit files for GEO/SRA:

```
geo_submission/
├── metadata.xlsx                   # Sample information table
│   ├── Sample Name
│   ├── Organism
│   ├── Tissue/Cell Type
│   ├── Treatment
│   ├── Sequencing Platform
│   └── ...
├── processed_files/
│   ├── bigwigs/
│   │   ├── sample1_CPM.bw → ../../bigwigs/sample1_CPM.bw
│   │   └── sample2_CPM.bw → ../../bigwigs/sample2_CPM.bw
│   └── peaks/
│       ├── sample1_peaks.bed → ../../peaks/macs2/sample1_peaks.narrowPeak
│       └── sample2_peaks.bed → ../../peaks/macs2/sample2_peaks.narrowPeak
└── raw_fastq/
    ├── sample1_R1.fastq.gz → ../../fastqs/sample1_R1.fastq.gz
    ├── sample1_R2.fastq.gz → ../../fastqs/sample1_R2.fastq.gz
    └── ...
```

## 🎨 Genome Browser Plots (PlotNado)

Publication-ready visualizations of genomic regions:

### Example Plot Configuration

```yaml
regions:
  - chr1:1000000-1005000
  - chr2:5000000-5010000
  
samples:
  - chip-rx_MLL
  - chip-rx_input
  
tracks:
  - bigwig: chip-rx_MLL_CPM.bw
    color: "#2E86AB"
    label: "MLL ChIP"
  - bigwig: chip-rx_input_CPM.bw  
    color: "#A23B72"
    label: "Input"
  - genes: hg38_refseq
```

### Output Files

```
genome_browser_plots/
├── chr1_1000000-1005000_MLL_vs_input.pdf
├── chr1_1000000-1005000_MLL_vs_input.png
├── chr2_5000000-5010000_MLL_vs_input.pdf
└── regions.bed
```

## 💡 Tips for Exploring Outputs

### Quick Quality Check

```bash
# Check main report
firefox seqnado_output/chip/seqnado_report.html

# Count peaks called  
wc -l seqnado_output/chip/peaks/macs2/*_peaks.narrowPeak

# View alignment stats
samtools flagstat seqnado_output/chip/aligned/sample.bam

# Check bigwig file
bigWigInfo seqnado_output/chip/bigwigs/bamnado/sample_CPM.bw
```

### Finding Specific Results

```bash
# All HTML reports
find seqnado_output/ -name "*.html"

# All peak files
find seqnado_output/ -name "*peaks*"

# All coverage tracks
find seqnado_output/ -name "*.bw" -o -name "*.bigWig"

# QC images
find seqnado_output/ -name "*.png" | grep qc
```

## 📚 Understanding File Formats

### BAM Files
- **Binary alignment format** (use `samtools view` to inspect)
- Stores aligned sequencing reads
- Includes alignment quality, CIGAR strings, and flags

### BigWig Files
- **Binary coverage track format**
- Efficient genome browser visualization
- Contains normalized signal values

### BED/NarrowPeak Files
- **Tab-delimited genomic coordinates**
- BED: chr, start, end, name, score, strand
- NarrowPeak: BED6 + fold-change, pvalue, qvalue, summit

### FastQ Files
- **Raw sequencing reads** (if retained)
- Four lines per read: header, sequence, +, quality scores
- Usually gzip compressed (.gz)

---

For more information on interpreting these outputs for your specific experiment, consult the [Pipeline Overview](pipeline.md) or [FAQ](faq.md).