# Third-Party Tools Reference

SeqNado integrates multiple best-in-class bioinformatics tools to provide comprehensive genomics analysis pipelines. This page documents the tools used, their purposes, and key references.

## Tool Versions & Updates

All tools are version-locked in SeqNado containers to ensure reproducibility. To check versions and get tool information, use the `seqnado tools` command:

```bash
# List all tools with versions
seqnado tools

# List tools in a specific category (e.g., Download, Alignment, Analysis)
seqnado tools --category

# View detailed information about a specific tool
seqnado tools macs2

# Show tool help/options from the container
seqnado tools macs2 --options
```

See the [CLI Reference](cli.md#cli-seqnado-tools) for complete documentation of the `tools` command.

## Citation Guidelines

When publishing results from SeqNado, please cite:

1. **SeqNado** itself (citation pending - check [GitHub releases](https://github.com/Milne-Group/SeqNado/releases) for latest version)
2. **Key tools** used in your specific analysis (see [References](#references) section below)
3. **Snakemake** workflow manager: Mölder F, Jablonski KP, Letcher B, et al. Sustainable data analysis with Snakemake. *F1000Research*. 2021;10:33.
4. **Reference genomes** used (e.g., hg38, mm10)

### Example Acknowledgment

*"Data analysis was performed using SeqNado v1.0.2 (Chahrour, C and Smith, AL, 2024), which incorporates Bowtie2 (Langmead & Salzberg, 2012) for alignment, MACS2 (Zhang et al., 2008) for peak calling, and deepTools (Ramírez et al., 2016) for coverage track generation. Workflows were managed with Snakemake (Mölder et al., 2021) within Apptainer containers."*

## Tools

Tools are organized by category, matching the structure of the `seqnado tools` CLI command. For more information about any tool, use `seqnado tools <toolname>`.

<!-- AUTO-GENERATED TOOL SECTIONS - DO NOT EDIT MANUALLY -->
<!-- This section is automatically updated by docs/scripts/generate_tool_citations.py -->
<!-- To update, run: python docs/scripts/generate_tool_citations.py --update -->

### Download

#### SRA Toolkit (fasterq-dump)
**Purpose**: Fast extraction of sequences from SRA files  
**Version**: 3.0.10  
**Usage**: Download sequencing data from NCBI Sequence Read Archive  
**Reference**: National Center for Biotechnology Information. SRA Toolkit: Sequence Read Archive tools and libraries. 2023. [https://github.com/ncbi/sra-tools.](https://github.com/ncbi/sra-tools.)


### Quality Control

#### FastQ Screen
**Purpose**: Screen sequencing reads against a set of databases  
**Version**: 0.16.0  
**Usage**: Screen reads against multiple reference genomes  
**Reference**: Steven W. Wingett and Simon Andrews. FastQ Screen: A tool for multi-genome mapping and quality control. F1000Research, 7(1338):1338, 2018. [https://doi.org/10.12688/f1000research.15931.2](https://doi.org/10.12688/f1000research.15931.2)


#### FastQC
**Purpose**: Quality control for high-throughput sequencing data  
**Version**: 0.12.1  
**Usage**: Generate quality metrics for sequencing data  
**Reference**: Simon Andrews. FastQC. 2010. [https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/.)


#### Qualimap
**Purpose**: Quality assessment tool for next-generation sequencing data  
**Version**: 2.3  
**Usage**: Compute alignment metrics and coverage statistics  
**Reference**: Konstantin Okonechnikov, Ana Conesa, and Fernando García-Alcalde. Qualimap 2: advanced multi-sample quality control for high-throughput sequencing data. Bioinformatics, 32(2):292–294, 2016. [https://doi.org/10.1093/bioinformatics/btv566](https://doi.org/10.1093/bioinformatics/btv566)


### Preprocessing

#### Cutadapt
**Purpose**: Remove adapter sequences from sequencing reads  
**Version**: 5.1  
**Usage**: Remove adapter sequences and trim low-quality bases  
**Reference**: Marcel Martin. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal, 17(1):10–12, 2011. [https://doi.org/10.14806/ej.17.1.200](https://doi.org/10.14806/ej.17.1.200)


#### FLASH
**Purpose**: Fast length adjustment of short reads (FLASH) to improve genome assemblies  
**Version**: 1.2.11  
**Usage**: Merge overlapping paired-end reads  
**Reference**: Tanja Magoč and Steven L. Salzberg. FLASH: fast length adjustment of short reads to improve genome assemblies. Bioinformatics, 27(21):2957–2963, 2011. [https://doi.org/10.1093/bioinformatics/btr507](https://doi.org/10.1093/bioinformatics/btr507)


#### Trim Galore
**Purpose**: Wrapper around Cutadapt and FastQC for quality and adapter trimming  
**Version**: 0.6.10  
**Usage**: Combines Cutadapt with quality control  
**Reference**: Felix Krueger. Trim galore. 2023. [https://www.bioinformatics.babraham.ac.uk/projects/trim\_galore/.](https://www.bioinformatics.babraham.ac.uk/projects/trim\_galore/.)


### Alignment

#### Bowtie2
**Purpose**: Fast and sensitive read mapping to large genomes  
**Version**: 2.5.4  
**Usage**: Align reads to reference genomes  
**Reference**: Ben Langmead and Steven L. Salzberg. Fast gapped-read alignment with Bowtie 2. Nature Methods, 9(4):357–359, 2012. [https://doi.org/10.1038/nmeth.1923](https://doi.org/10.1038/nmeth.1923)


#### minimap2
**Purpose**: Fast pairwise sequence alignment tool  
**Version**: 2.30  
**Usage**: Fast alignment for long reads (PacBio/Nanopore)  
**Reference**: Heng Li. Minimap2: pairwise alignment for nucleotide sequences. Bioinformatics, 34(18):3094–3100, 2018. [https://doi.org/10.1093/bioinformatics/bty191](https://doi.org/10.1093/bioinformatics/bty191)


#### Picard
**Purpose**: Tools for manipulating high-throughput sequencing data and formats  
**Version**: 3.4.0  
**Usage**: Duplicate marking and BAM processing  
**Reference**: Picard toolkit. 2019. [https://broadinstitute.github.io/picard/.](https://broadinstitute.github.io/picard/.)


#### SAMtools
**Purpose**: Tools for interacting with SAM/BAM format files  
**Version**: 1.22.1  
**Usage**: Sorting, indexing, filtering SAM/BAM files  
**Reference**: Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, and Heng Li. Twelve years of SAMtools and BCFtools. GigaScience, 10(2):giab008, 2021. [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)


#### STAR
**Purpose**: Ultrafast universal RNA-seq aligner  
**Version**: 2.7.11b  
**Usage**: High-speed splicing-aware RNA-seq alignment  
**Reference**: Alexander Dobin, Carrie A. Davis, Felix Schlesinger, Jorg Drenkow, Chris Zaleski, Sonali Jha, Philippe Batut, Mark Chaisson, and Thomas R. Gingeras. STAR: ultrafast universal RNA-seq aligner. Bioinformatics, 29(1):15–21, 2013. [https://doi.org/10.1093/bioinformatics/bts635](https://doi.org/10.1093/bioinformatics/bts635)


### Analysis

#### BamNado
**Purpose**: SeqNado BAM file manipulation and analysis tool  
**Version**: 0.4.4  
**Usage**: Calculate scaling factors and spike-in normalization  
**Reference**: Smith, Alastair L. BamNado: BAM file manipulation and analysis tool. 2024. [https://pypi.org/project/bamnado/.](https://pypi.org/project/bamnado/.)


#### BCFtools
**Purpose**: Tools for BCF/VCF format variant files  
**Version**: 1.22  
**Usage**: Process VCF/BCF variant files  
**Reference**: Petr Danecek, James K Bonfield, Jennifer Liddle, John Marshall, Valeriu Ohan, Martin O Pollard, Andrew Whitwham, Thomas Keane, Shane A McCarthy, Robert M Davies, and Heng Li. Twelve years of SAMtools and BCFtools. GigaScience, 10(2):giab008, 2021. [https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)


#### BEDTools
**Purpose**: Tools for genomic arithmetic and feature analysis  
**Version**: 2.31.1  
**Usage**: Genomic interval manipulation and analysis  
**Reference**: Aaron R. Quinlan and Ira M. Hall. BEDTools: a flexible suite of utilities for comparing genomic features. Bioinformatics, 26(6):841–842, 2010. [https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)


#### Cooler
**Purpose**: Tools for high-resolution interactions and HiC contact matrices  
**Version**: 0.10.4  
**Usage**: Handle genomically-labeled arrays and Hi-C data  
**Reference**: Nezar Abdennur and Leonid A. Mirny. Cooler: scalable storage for Hi-C data and other genomically labeled arrays. Bioinformatics, 36(1):311–316, 2020. [https://doi.org/10.1093/bioinformatics/btz540](https://doi.org/10.1093/bioinformatics/btz540)


#### findPeaks (HOMER)
**Purpose**: Identify genomic peaks in ChIP-seq or ATAC-seq data (HOMER)  
**Version**: 5.1  
**Usage**: Identify genomic peaks in ChIP-seq or ATAC-seq data  
**Reference**: Sven Heinz, Christopher Benner, Nathanael Spann, Eric Bertolino, Yin C. Lin, Peter Laslo, Jason X. Cheng, Cornelis Murre, Harinder Singh, and Christopher K. Glass. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell, 38(4):576–589, 2010. [https://doi.org/10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)


#### HOMER
**Purpose**: Tools for ChIP-seq, DNase-seq, and other genomic analysis  
**Version**: 5.1  
**Usage**: Peak calling, motif discovery, annotation  
**Reference**: Sven Heinz, Christopher Benner, Nathanael Spann, Eric Bertolino, Yin C. Lin, Peter Laslo, Jason X. Cheng, Cornelis Murre, Harinder Singh, and Christopher K. Glass. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell, 38(4):576–589, 2010. [https://doi.org/10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)


#### LanceOtron
**Purpose**: Peak caller for high-resolution chromatin analysis  
**Version**: 1.2.7  
**Usage**: Machine learning-based peak calling  
**Reference**: Lance D Hentges, Martin J Sergeant, Christopher B Cole, Damien J Downes, Jim R Hughes, and Stephen Taylor. LanceOtron: a deep learning peak caller for genome sequencing experiments. Bioinformatics, pages btac525, 2022. [https://doi.org/10.1093/bioinformatics/btac525](https://doi.org/10.1093/bioinformatics/btac525)


#### MACS2
**Purpose**: Model-based analysis of ChIP-Seq data  
**Version**: 2.2.9.1  
**Usage**: Peak calling for ChIP-seq and ATAC-seq  
**Reference**: Yong Zhang, Tao Liu, Clifford A Meyer, Jérôme Eeckhoute, David S Johnson, Bradley E Bernstein, Chad Nusbaum, Richard M Myers, Myles Brown, Wei Li, and X Shirley Liu. Model-based analysis of ChIP-seq (MACS). Genome Biology, 9(9):R137, 2008. [https://doi.org/10.1186/gb-2008-9-9-r137](https://doi.org/10.1186/gb-2008-9-9-r137)


#### MAGeCK
**Purpose**: Model-based Analysis of Genome-wide CRISPR-Cas9 Knockout data  
**Version**: 0.5.9.5  
**Usage**: CRISPR screen analysis for essential gene identification  
**Reference**: Wei Li, Han Xu, Tengfei Xiao, Le Cong, Michael I Love, Feng Zhang, Rafael A Irizarry, Jun S Liu, Myles Brown, and X Shirley Liu. MAGeCK enables robust identification of essential genes from genome-scale CRISPR/Cas9 knockout screens. Genome Biology, 15(12):554, 2014. [https://doi.org/10.1186/s13059-014-0554-4](https://doi.org/10.1186/s13059-014-0554-4)


#### makeTagDirectory (HOMER)
**Purpose**: Prepare and organize ChIP-seq tag directories (HOMER)  
**Version**: 5.1  
**Usage**: Prepare and organize ChIP-seq tag directories  
**Reference**: Sven Heinz, Christopher Benner, Nathanael Spann, Eric Bertolino, Yin C. Lin, Peter Laslo, Jason X. Cheng, Cornelis Murre, Harinder Singh, and Christopher K. Glass. Simple Combinations of Lineage-Determining Transcription Factors Prime cis-Regulatory Elements Required for Macrophage and B Cell Identities. Molecular Cell, 38(4):576–589, 2010. [https://doi.org/10.1016/j.molcel.2010.05.004](https://doi.org/10.1016/j.molcel.2010.05.004)


#### MCCNado
**Purpose**: Micro-Capture-C (MCC) sequencing data processing tool  
**Version**: 0.1.6  
**Usage**: Analyze chromatin 3D interactions (Micro-Capture-C)  
**Reference**: Alastair L. Smith. MCCNado: Micro-Capture-C data processing tool. 2024. [https://pypi.org/project/mccnado.](https://pypi.org/project/mccnado.)


#### MEME Suite
**Purpose**: MEME Suite motif discovery and analysis tools  
**Version**: 5.5.9  
**Usage**: DNA motif discovery and analysis  
**Reference**: T. L. Bailey, M. Boden, F. A. Buske, M. Frith, C. E. Grant, L. Clementi, J. Ren, W. W. Li, and W. S. Noble. Meme SUITE: tools for motif discovery and searching. Nucleic Acids Research, 37(Web Server):W202–W208, 2009. [https://doi.org/10.1093/nar/gkp335](https://doi.org/10.1093/nar/gkp335)


#### MethylDackel
**Purpose**: Tools for analyzing DNA methylation from bisulfite sequencing data  
**Version**: 0.6.1  
**Usage**: Extract methylation metrics from bisulfite-seq  
**Reference**: Devon Ryan. MethylDackel. 2021. [https://github.com/dpryan79/MethylDackel.](https://github.com/dpryan79/MethylDackel.)


#### SEACR
**Purpose**: Sparse enrichment analysis for CUT&RUN  
**Version**: 1.3  
**Usage**: Peak calling optimized for CUT&Tag/CUT&RUN  
**Reference**: Michael P. Meers, Dan Tenenbaum, and Steven Henikoff. Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling. Epigenetics & Chromatin, 12(1):42, 2019. [https://doi.org/10.1186/s13072-019-0287-4](https://doi.org/10.1186/s13072-019-0287-4)


### Visualization

#### deepTools
**Purpose**: Tools for processing and visualizing deep sequencing data  
**Version**: 3.5.6  
**Usage**: BigWig generation, heatmaps, and profile plots  
**Reference**: Fidel Ramírez, Devon P Ryan, Björn Grüning, Vivek Bhardwaj, Fabian Kilpert, Andreas S Richter, Steffen Heyne, Friederike Dündar, and Thomas Manke. deepTools2: a next generation web server for deep-sequencing data analysis. Nucleic Acids Research, 44(W1):W160–W165, 2016. [https://doi.org/10.1093/nar/gkw257](https://doi.org/10.1093/nar/gkw257)


#### PlotNado
**Purpose**: SeqNado genome browser visualization and plotting tool  
**Version**: 0.1.dev101  
**Usage**: Publication-ready genomic region visualizations  
**Reference**: Alastair L. Smith. PlotNado: genome browser visualization tool. 2023. [https://github.com/alsmith151/plotnado.](https://github.com/alsmith151/plotnado.)


#### trackhub
**Purpose**: Python library to work with track hubs  
**Version**: 1.0  
**Usage**: UCSC Genome Browser track hub library  
**Reference**: Ryan K. Dale, Laura H. Matzat, and Elissa P. Lei. Trackhub: A Library for creating and remotely hosting UCSC Genome Browser track hubs. F1000Research, 10:121, 2021. [https://doi.org/10.12688/f1000research.50107.2](https://doi.org/10.12688/f1000research.50107.2)


#### TrackNado
**Purpose**: Track hub creation tool for UCSC Genome Browser  
**Version**: ≥0.3.1  
**Usage**: Create and host UCSC track hubs  
**Reference**: Alastair L. Smith. TrackNado: a python library and cli tool to rapidly generate complex ucsc genome browser track hubs. 2024. [https://pypi.org/project/tracknado.](https://pypi.org/project/tracknado.)


### Reporting

#### MultiQC
**Purpose**: Aggregate quality control reports from multiple QC tools  
**Version**: 1.31  
**Usage**: Aggregate QC reports across samples into interactive HTML  
**Reference**: Philip Ewels, Måns Magnusson, Sverker Lundin, and Max Käller. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics, 32(19):3047–3048, 2016. [https://doi.org/10.1093/bioinformatics/btw354](https://doi.org/10.1093/bioinformatics/btw354)


#### Quarto
**Purpose**: Scientific and technical publishing system for reports  
**Version**: 1.8.25  
**Usage**: Generate analysis reports with code and narrative  
**Reference**: Posit Software, PBC. Quarto: An open-source scientific and technical publishing system. 2022. [https://quarto.org/.](https://quarto.org/.)


### Quantification

#### featureCounts
**Purpose**: Count reads overlapping genomic features (genes/exons)  
**Version**: 2.1.1  
**Usage**: Gene-level read quantification for RNA-seq  
**Reference**: Yang Liao, Gordon K. Smyth, and Wei Shi. featureCounts: an efficient general purpose program for assigning sequence reads to genomic features. Bioinformatics, 30(7):923–930, 2014. [https://doi.org/10.1093/bioinformatics/btt656](https://doi.org/10.1093/bioinformatics/btt656)


#### QuantNado
**Purpose**: SeqNado multiomics quantification and dataset creation tool  
**Version**: 0.2.5  
**Usage**: Generate quantified datasets for ML analysis  
**Reference**: Catherine Chahrour and Alastair L. Smith. QuantNado: multiomics machine learning dataset generation tool. 2024. [https://pypi.org/project/QuantNado.](https://pypi.org/project/QuantNado.)


#### Salmon
**Purpose**: Quantification of gene expression from RNA-seq data  
**Version**: 1.10.3  
**Usage**: Alignment-free transcript quantification  
**Reference**: Rob Patro, Geet Duggal, Michael I. Love, Rafael A. Irizarry, and Carl Kingsford. Salmon provides fast and bias-aware quantification of transcript expression. Nature Methods, 14(4):417–419, 2017. [https://doi.org/10.1038/nmeth.4197](https://doi.org/10.1038/nmeth.4197)


#### Subread
**Purpose**: High-performance read alignment, quantification and mutation discovery  
**Version**: 2.1.1  
**Usage**: High-performance read alignment and quantification  
**Reference**: Yang Liao, Gordon K. Smyth, and Wei Shi. The Subread aligner: fast, accurate and scalable read mapping by seed-and-vote. Nucleic Acids Research, 41(10):e108, 2013. [https://doi.org/10.1093/nar/gkt214](https://doi.org/10.1093/nar/gkt214)


### Utilities

#### Apptainer
**Purpose**: Container runtime for running Singularity/Apptainer images on HPC clusters  
**Version**: User-managed (see installation docs)  
**Usage**: Run containerized tools in HPC environments  
**Reference**: Apptainer Project Contributors. Apptainer. 2023. [https://apptainer.org/.](https://apptainer.org/.)


#### bedToBigBed
**Purpose**: Convert BED format files to compressed BigBed format  
**Version**: 482  
**Usage**: Convert BED to BigBed format  
**Reference**: W. James Kent, Charles W. Sugnet, Terrence S. Furey, Krishna M. Roskin, Tom H. Pringle, Alan M. Zahler, and David Haussler. The human genome browser at UCSC. Genome Research, 12(5):996–1006, 2002. [https://doi.org/10.1101/gr.229102](https://doi.org/10.1101/gr.229102)


#### pigz
**Purpose**: Parallel gzip compression utility  
**Version**: 2.8  
**Usage**: Fast parallel compression/decompression  
**Reference**: Mark Adler. Pigz: Parallel Implementation of GZip. 2015. [https://zlib.net/pigz/.](https://zlib.net/pigz/.)


#### Snakemake
**Purpose**: Workflow management system for reproducible and scalable data analysis  
**Version**: ≥9.12.0  
**Usage**: Workflow management and execution  
**Reference**: Felix Mölder, Kim Philipp Jablonski, Brice Letcher, Michael B. Hall, Christopher H. Tomkins-Tinch, Vanessa Sochat, Jan Forster, Soohyun Lee, Sven O. Twardziok, Alexander Kanitz, Andreas Wilm, Manuel Holtgrewe, Sven Rahmann, Sven Nahnsen, and Johannes Köster. Sustainable data analysis with Snakemake. F1000Research, 10:33, 2021. [https://doi.org/10.12688/f1000research.29032.3](https://doi.org/10.12688/f1000research.29032.3)


#### UCSC Tools
**Purpose**: UCSC genome browser command-line utilities  
**Version**: 482  
**Usage**: BigWig creation and genomic format conversion  
**Reference**: W. James Kent, Charles W. Sugnet, Terrence S. Furey, Krishna M. Roskin, Tom H. Pringle, Alan M. Zahler, and David Haussler. The human genome browser at UCSC. Genome Research, 12(5):996–1006, 2002. [https://doi.org/10.1101/gr.229102](https://doi.org/10.1101/gr.229102)

