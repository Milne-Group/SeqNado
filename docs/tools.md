# Third-Party Tools Reference

SeqNado integrates multiple best-in-class bioinformatics tools to provide comprehensive genomics analysis pipelines. This page documents the tools used, their purposes, and key references.

## 🧬 Core Analysis Tools

### Alignment & Mapping

#### Bowtie2
**Purpose**: Short read aligner for ChIP-seq, ATAC-seq, and CUT&Tag  
**Version**: Latest via container  
**Usage**: Aligns reads to reference genomes with flexible parameters  
**Reference**: [Langmead & Salzberg, 2012](http://www.nature.com/nmeth/journal/v9/n4/full/nmeth.1923.html)  
**Citation**: Langmead B, Salzberg S. Fast gapped-read alignment with Bowtie 2. *Nature Methods*. 2012, 9:357-359.

**Key Options in SeqNado:**
- `--very-sensitive`: High accuracy mode
- `--no-mixed/--no-discordant`: Paired-end specific settings
- `--maxins`: Maximum insert size for PE data

#### HISAT2
**Purpose**: Splicing-aware aligner for RNA-seq  
**Version**: Latest via container  
**Usage**: Maps RNA-seq reads accounting for splice junctions  
**Reference**: [Kim et al., 2019](https://www.nature.com/articles/s41587-019-0201-4)  
**Citation**: Kim D, Paggi JM, Park C, Bennett C, Salzberg SL. Graph-based genome alignment and genotyping with HISAT2 and HISAT-genotype. *Nature Biotechnology*. 2019, 37:907-915.

**Key Features:**
- Graph-based alignment
- Splice junction detection
- SNP-aware alignment

### Quality Control

#### FastQC
**Purpose**: Quality assessment of raw sequencing data  
**Version**: 0.12.1+  
**Usage**: Generates per-base quality scores, GC content, adapter contamination metrics  
**Reference**: [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)  

**Metrics Reported:**
- Per base sequence quality
- Per sequence quality scores
- Per base sequence content
- Sequence duplication levels
- Overrepresented sequences
- Adapter content

#### MultiQC
**Purpose**: Aggregate QC reports across samples  
**Version**: Latest  
**Usage**: Combines FastQC, alignment stats, and tool outputs into interactive HTML  
**Reference**: [Ewels et al., 2016](https://academic.oup.com/bioinformatics/article/32/19/3047/2196507)  
**Citation**: Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*. 2016, 32(19):3047-3048.

#### Qualimap
**Purpose**: BAM file quality control  
**Version**: 2.2.1+  
**Usage**: Computes alignment metrics, coverage statistics, insert size distributions  
**Reference**: [Okonechnikov et al., 2016](https://academic.oup.com/bioinformatics/article/32/2/292/1743969)  

**Key Metrics:**
- Coverage distribution
- Insert size distribution
- Mapping quality distribution
- GC content
- Genomic origin of reads

#### FastQ Screen
**Purpose**: Contamination screening  
**Version**: 0.15+  
**Usage**: Screens reads against multiple reference genomes to detect contamination  
**Reference**: [Babraham Bioinformatics](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/)  

**Screens Against:**
- Human (hg38)
- Mouse (mm10)
- E. coli
- PhiX
- Adapters
- Custom genomes

### Read Processing

#### Cutadapt
**Purpose**: Adapter trimming and quality filtering  
**Version**: 4.0+  
**Usage**: Removes adapter sequences, trims low-quality bases  
**Reference**: [Martin, 2011](https://journal.embnet.org/index.php/embnetjournal/article/view/200)  
**Citation**: Martin M. Cutadapt removes adapter sequences from high-throughput sequencing reads. *EMBnet.journal*. 2011, 17(1):10-12.

**SeqNado Configuration:**
- Adapter auto-detection
- Quality trimming (Q20)
- Minimum length filtering
- Poly-A/T tail removal

#### Trimmomatic
**Purpose**: Alternative adapter trimming (optional)  
**Version**: 0.39+  
**Usage**: Flexible trimming with sliding window approach  
**Reference**: [Bolger et al., 2014](https://academic.oup.com/bioinformatics/article/30/15/2114/2390096)  

### SAM/BAM Processing

#### SAMtools
**Purpose**: SAM/BAM file manipulation  
**Version**: 1.17+  
**Usage**: Sorting, indexing, filtering, statistics  
**Reference**: [Li et al., 2009](https://academic.oup.com/bioinformatics/article/25/16/2078/204688)  
**Citation**: Li H, Handsaker B, Wysoker A, et al. The Sequence Alignment/Map format and SAMtools. *Bioinformatics*. 2009, 25(16):2078-2079.

**Common Operations:**
- `samtools sort`: Coordinate sorting
- `samtools index`: BAM indexing
- `samtools flagstat`: Alignment statistics
- `samtools view`: BAM filtering

#### Picard
**Purpose**: Java-based BAM processing  
**Version**: 2.27+  
**Usage**: Duplicate marking, insert size metrics  
**Reference**: [Broad Institute](http://broadinstitute.github.io/picard/)  

**Key Tools:**
- MarkDuplicates: PCR duplicate identification
- CollectInsertSizeMetrics: PE fragment analysis
- CollectAlignmentSummaryMetrics: Detailed alignment stats

## 🔬 Peak Calling & Analysis

### MACS2/MACS3
**Purpose**: Peak calling for ChIP-seq and ATAC-seq  
**Version**: MACS2 2.2.7+ / MACS3 3.0+  
**Usage**: Identifies enriched genomic regions  
**Reference**: [Zhang et al., 2008](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2008-9-9-r137) (MACS2)  
**Citation**: Zhang Y, Liu T, Meyer CA, et al. Model-based Analysis of ChIP-Seq (MACS). *Genome Biology*. 2008, 9:R137.

**Peak Types:**
- `narrowPeak`: Sharp peaks (TFs, ATAC)
- `broadPeak`: Broad domains (histone marks)

**SeqNado Parameters:**
- Auto genome size detection
- FDR cutoff (default 0.05)
- Optional input/control
- Spike-in normalization

### HOMER
**Purpose**: Comprehensive ChIP-seq analysis suite  
**Version**: Latest  
**Usage**: Peak calling, motif discovery, annotation  
**Reference**: [Heinz et al., 2010](https://www.cell.com/molecular-cell/fulltext/S1097-2765(10)00334-6)  
**Citation**: Heinz S, Benner C, et al. Simple combinations of lineage-determining transcription factors prime cis-regulatory elements required for macrophage and B cell identities. *Molecular Cell*. 2010, 38(4):576-589.

**Features:**
- Tag directory creation
- Peak calling (findPeaks)
- Motif discovery (findMotifsGenome)
- Peak annotation (annotatePeaks.pl)

### SEACR
**Purpose**: Sparse Enrichment Analysis for CUT&RUN  
**Version**: 1.3+  
**Usage**: Peak calling optimized for CUT&Tag/CUT&RUN  
**Reference**: [Meers et al., 2019](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)  
**Citation**: Meers MP, Tenenbaum D, Henikoff S. Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. *Epigenetics & Chromatin*. 2019, 12:42.

**Advantages:**
- Low background noise tolerance
- No input control required
- Optimized for sparse data

### LanceOtron
**Purpose**: Machine learning-based peak calling  
**Version**: Latest  
**Usage**: Neural network peak caller for ChIP/ATAC  
**Reference**: [Lance, 2021](https://www.biorxiv.org/content/10.1101/2021.01.07.425361v1)  

**Benefits:**
- Learns from data patterns
- Reduces false positives
- Works across assay types

## 📊 Coverage & Visualization

### deepTools
**Purpose**: Comprehensive genomics data analysis and visualization  
**Version**: 3.5+  
**Usage**: BigWig generation, heatmaps, profile plots  
**Reference**: [Ramírez et al., 2016](https://academic.oup.com/nar/article/44/W1/W160/2499308)  
**Citation**: Ramírez F, Ryan DP, et al. deepTools2: a next generation web server for deep-sequencing data analysis. *Nucleic Acids Research*. 2016, 44(W1):W160-W165.

**Key Tools:**
- `bamCoverage`: BAM to BigWig conversion
- `bamCompare`: Compute log2 ratios
- `computeMatrix`: Generate matrices for plotting
- `plotHeatmap`: Heatmap generation
- `plotProfile`: Meta-plots

**Normalization Methods:**
- CPM (Counts Per Million)
- RPKM (Reads Per Kilobase per Million)
- RPGC (Reads Per Genomic Content)
- BPM (Bins Per Million)

### bedtools
**Purpose**: Genomic interval manipulation  
**Version**: 2.30+  
**Usage**: Intersect, merge, coverage operations  
**Reference**: [Quinlan & Hall, 2010](https://academic.oup.com/bioinformatics/article/26/6/841/244688)  
**Citation**: Quinlan AR, Hall IM. BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*. 2010, 26(6):841-842.

**Common Operations:**
- `bedtools intersect`: Find overlapping features
- `bedtools merge`: Combine overlapping intervals
- `bedtools coverage`: Calculate coverage depth
- `bedtools genomecov`: Genome-wide coverage

### UCSC Tools
**Purpose**: Genome browser data format conversion  
**Version**: Latest  
**Usage**: BigWig creation and manipulation  
**Reference**: [Kent et al., 2010](https://academic.oup.com/bioinformatics/article/26/17/2204/199001)  

**Key Tools:**
- `bedGraphToBigWig`: Convert bedGraph to BigWig
- `bigWigInfo`: BigWig file information
- `bigWigMerge`: Combine BigWig files

## 🧮 RNA-seq Specific

### STAR
**Purpose**: Ultra-fast RNA-seq aligner  
**Version**: 2.7+  
**Usage**: High-speed splicing-aware alignment  
**Reference**: [Dobin et al., 2013](https://academic.oup.com/bioinformatics/article/29/1/15/272537)  
**Citation**: Dobin A, Davis CA, et al. STAR: ultrafast universal RNA-seq aligner. *Bioinformatics*. 2013, 29(1):15-21.

**Features:**
- Two-pass alignment mode
- Splice junction discovery
- Chimeric read detection
- Gene/transcript quantification

### featureCounts
**Purpose**: Read summarization for RNA-seq  
**Version**: 2.0+  
**Usage**: Gene-level quantification  
**Reference**: [Liao et al., 2014](https://academic.oup.com/bioinformatics/article/30/7/923/232889)  
**Citation**: Liao Y, Smyth GK, Shi W. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. *Bioinformatics*. 2014, 30(7):923-930.

### DESeq2
**Purpose**: Differential expression analysis  
**Version**: 1.40+  
**Usage**: Statistical testing for RNA-seq  
**Reference**: [Love et al., 2014](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8)  
**Citation**: Love MI, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. *Genome Biology*. 2014, 15:550.

## 🔐 Spike-in Normalization

### Bamnado
**Purpose**: Spike-in aware BAM processing  
**Version**: Latest  
**Usage**: Calculate spike-in normalization factors  
**Reference**: Part of SeqNado ecosystem  

**Workflow:**
1. Map to spike-in genome
2. Calculate scaling factors
3. Normalize experimental data
4. Generate comparable BigWigs

## 🎨 Visualization

### PlotNado
**Purpose**: Genome browser-style plots  
**Version**: Latest (TrackNado)  
**Usage**: Publication-ready genomic region visualizations  
**Reference**: Part of SeqNado/TrackNado ecosystem  

**Features:**
- Multi-track plots
- Gene annotation overlay
- Custom color schemes
- PDF/PNG output

### IGV (Integrative Genomics Viewer)
**Purpose**: Interactive genome browser  
**Version**: Latest  
**Usage**: View BAM, BigWig, and peak files  
**Reference**: [Robinson et al., 2011](https://www.nature.com/articles/nbt.1754)  
**Citation**: Robinson JT, Thorvaldsdóttir H, et al. Integrative genomics viewer. *Nature Biotechnology*. 2011, 29:24-26.

## 📦 Container Technology

### Apptainer/Singularity
**Purpose**: Container runtime for HPC  
**Version**: Latest  
**Usage**: Ensures reproducibility and portability  
**Reference**: [Kurtzer et al., 2017](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0177459)  

**Benefits:**
- No root privileges required
- HPC-friendly
- Bundles all dependencies
- Version locked for reproducibility

## 📝 Workflow Management

### Snakemake
**Purpose**: Workflow management system  
**Version**: 9.12+  
**Usage**: Orchestrates analysis pipelines  
**Reference**: [Mölder et al., 2021](https://f1000research.com/articles/10-33)  
**Citation**: Mölder F, Jablonski KP, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021, 10:33.

**Features:**
- Automatic parallelization
- Resource management
- Checkpoint/restart capability
- Conda/container integration

## 🔢 Statistical & Methylation Tools

### Bismark
**Purpose**: Bisulfite sequencing alignment and methylation calling  
**Version**: 0.24+  
**Usage**: WGBS and RRBS analysis  
**Reference**: [Krueger & Andrews, 2011](https://academic.oup.com/bioinformatics/article/27/11/1571/216956)  

### Preseq
**Purpose**: Library complexity estimation  
**Version**: 3.1+  
**Usage**: Predict sequencing saturation  
**Reference**: [Daley & Smith, 2013](https://www.nature.com/articles/nmeth.2375)  

## 🌐 Tool Versions & Updates

All tools are version-locked in SeqNado containers to ensure reproducibility. To check versions:

```bash
# View container manifest
apptainer exec seqnado.sif cat /opt/versions.txt

# Or check specific tool
apptainer exec seqnado.sif macs2 --version
```

## 📚 Citation Guidelines

When publishing results from SeqNado, please cite:

1. **SeqNado** itself (citation pending)
2. **Key tools used** in your specific analysis (see above)
3. **Snakemake** workflow manager
4. **Reference genomes** used

Example acknowledgment:

> "Data analysis was performed using SeqNado v1.0, which incorporates Bowtie2 for alignment, MACS2 for peak calling, and deepTools for coverage track generation. Workflows were managed with Snakemake within Apptainer containers."

## 🔗 External Resources

- [SeqNado GitHub](https://github.com/Milne-Group/SeqNado)
- [Snakemake Documentation](https://snakemake.readthedocs.io/)
- [deepTools Documentation](https://deeptools.readthedocs.io/)
- [MACS GitHub](https://github.com/macs3-project/MACS)
- [HOMER Website](http://homer.ucsd.edu/homer/)

---

*This page is automatically updated with each SeqNado release to reflect current tool versions and integrations.*