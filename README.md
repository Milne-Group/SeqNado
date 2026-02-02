<!-- Header and badges -->
# SeqNado

[![Documentation](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml/badge.svg)](https://github.com/alsmith151/SeqNado/actions/workflows/build_docs.yml)
[![Bioconda](https://anaconda.org/bioconda/seqnado/badges/version.svg)](https://anaconda.org/bioconda/seqnado)
[![Bioconda Updated](https://anaconda.org/bioconda/seqnado/badges/latest_release_date.svg)](https://anaconda.org/bioconda/seqnado)
[![PyPI Downloads](https://static.pepy.tech/badge/seqnado)](https://pepy.tech/projects/seqnado)

<p align="center">
  <img src="https://raw.githubusercontent.com/alsmith151/SeqNado/main/containers/pipeline/seqnado.png" alt="SeqNado logo" />
</p>

*A Snakemake-based bioinformatics toolkit for analyzing sequencing data from ATAC-seq, ChIP-seq, CUT&Tag, RNA-seq, SNP analysis, Methylation, CRISPR screens, and Micro-Capture-C experiments.*

## Table of Contents

- [Key Features](#key-features)
- [Supported Assays](#supported-assays)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Documentation](#documentation)
- [License](#license)

---

Modular, reproducible, and container-ready pipelines powered by Snakemake that take you from raw data to publication-ready results.

## Key Features

- **Comprehensive Assay Support**: Single framework for multiple sequencing assays
- **Customizable Workflows**: Easily modify parameters, use different tools for peak calling, bigwig generation etc.
- **User-Friendly CLI**: Intuitive command-line interface that guides you through setup and execution
- **Multiomics Support**: Analyze and integrate data from multiple sequencing assays in a single workflow
- **Snakemake-Powered**: Modular workflows with automatic parallelization and resource management
- **Container-Ready**: Fully containerized pipelines using Apptainer/Singularity for reproducibility
- **HPC-Optimized**: Seamless integration with SLURM and local execution modes
- **Advanced Analysis**: 
  - Comprehensive QC with MultiQC reports
  - Peak calling with MACS2, SEACR, HOMER, and LanceOtron
  - Consensus peakset generation and quantification across samples
  - Spike-in normalization for ChIP-seq, ATAC-seq, and RNA-seq
  - Automated differential expression with DESeq2 for RNA-seq
  - Genome browser style plots with `PlotNado`
  - UCSC genome browser hub generation
  - ML-ready dataset creation
- **Flexible Configuration**: Interactive CLI for setup, or scriptable non-interactive mode
- **Machine Learning Ready**: Tools for preparing datasets for ML applications

## Supported Assays

- **ATAC-seq** (`atac`) - Chromatin accessibility profiling with TSS enrichment and fragment analysis
- **ChIP-seq** (`chip`) - Protein-DNA interaction mapping with spike-in support
- **CUT&Tag** (`cat`) - Low-input epigenomic profiling optimized for sparse signals
- **RNA-seq** (`rna`) - Transcriptome analysis with automated DESeq2 differential expression
- **SNP Analysis** (`snp`) - Variant detection and genotyping workflows
- **Methylation** (`meth`) - Bisulfite/TAPS sequencing for DNA methylation analysis
- **CRISPR Screens** (`crispr`) - Guide-level quantification and screen statistics
- **Micro-Capture-C** (`mcc`) - Chromatin conformation capture analysis
- **Multiomics** - Run multiple assay types together in a single integrated workflow

→ [View detailed assay workflows](https://Milne-Group.github.io/SeqNado/pipeline/)

## Installation

### Via Mamba (Recommended)

Install from the Bioconda channel:

```bash
mamba create -n seqnado -c bioconda seqnado
mamba activate seqnado
```

### Via uv (Fast Alternative)

Install using [uv](https://docs.astral.sh/uv/), a fast Python package installer:

```bash
uv venv seqnado-env
source seqnado-env/bin/activate  # On macOS/Linux; use 'seqnado-env\Scripts\activate' on Windows
uv pip install seqnado
```

### Via Pip

Alternatively, install using pip:

```bash
pip install seqnado
```

### Initialize SeqNado

After installation, initialize your SeqNado environment:

```bash
seqnado init
```

**What this does:**
- Sets up genome configuration templates in `~/.config/seqnado/`
- Configures Apptainer/Singularity containers (if available)
- Installs Snakemake execution profiles for local and cluster execution

→ [Learn more about initialization](https://Milne-Group.github.io/SeqNado/initialisation/)

## Quick Start

Complete workflow from installation to results in 5 steps:

### 1. Set Up Genome References

Before processing data, configure reference genomes for alignment:

```bash
# List available genomes
seqnado genomes list atac

# Build a custom genome
seqnado genomes build rna --fasta hg38.fasta --name hg38 --outdir /path/to/genomes
```

→ [Complete genome setup guide](https://Milne-Group.github.io/SeqNado/genomes/)

### 2. Create Project Configuration

Generate a configuration file and project directory for your experiment:

```bash
seqnado config atac
```

**Output:** A dated project directory with configuration file and FASTQ folder:
```
YYYY-MM-DD_ATAC_project/
├── config_atac.yaml    # Edit this to customize analysis parameters
└── fastqs/             # Place your FASTQ files here
```

→ [Configuration options guide](https://Milne-Group.github.io/SeqNado/configuration/)

### 3. Add FASTQ Files

Symlink your raw sequencing data into the project directory:

```bash
ln -s /path/to/fastq/*.fastq.gz YYYY-MM-DD_ATAC_project/fastqs/
```

**Note:** Use symbolic links to avoid duplicating large files.

### 4. Generate Sample Metadata

Create a metadata CSV that describes your experimental design:

```bash
seqnado design atac
```

**Output:** `metadata_atac.csv` — Edit this file to specify:
- Sample names and groupings
- Experimental conditions
- Control/treatment relationships
- DESeq2 comparisons (for RNA-seq)

→ [Design file specification](https://Milne-Group.github.io/SeqNado/design/)

### 5. Run the Pipeline

Execute the analysis workflow (choose one based on your environment):

```bash
# Local machine (uses all available cores)
seqnado pipeline atac --preset le

# HPC cluster with SLURM scheduler
seqnado pipeline atac --preset ss --queue short

# Multiomics mode (processes multiple assays together)
seqnado pipeline  --preset ss# Detects all config files in current directory
```
→ [Pipeline execution details](https://Milne-Group.github.io/SeqNado/pipeline/) | [Output files explained](https://Milne-Group.github.io/SeqNado/outputs/)
### Common Pipeline Options

**Execution Presets:**
- `--preset le` - Local execution (default, recommended for workstations)
- `--preset lc` - Local execution using conda environments
- `--preset ss` - SLURM scheduler (for HPC clusters)

**Resource Management:**
- `--queue short` - Specify SLURM partition/queue name
- `--scale-resources 1.5` - Multiply memory/time requirements by 1.5×

**Debugging & Testing:**
- `-n` - Dry run to preview commands without executing
- `--unlock` - Unlock directory after interrupted runs

→ [All CLI options](https://Milne-Group.github.io/SeqNado/cli/) | [HPC cluster setup](https://Milne-Group.github.io/SeqNado/cluster_config/)

## Documentation

For comprehensive guides and API documentation, visit:

**📚 [SeqNado Documentation](https://Milne-Group.github.io/SeqNado/)**

### Key Topics

- [Installation Guide](https://Milne-Group.github.io/SeqNado/installation/)
- [Genome Setup](https://Milne-Group.github.io/SeqNado/genomes/)
- [Configuration](https://Milne-Group.github.io/SeqNado/configuration/)
- [Design Files](https://Milne-Group.github.io/SeqNado/design/)
- [Pipeline Details](https://Milne-Group.github.io/SeqNado/pipeline/)
- [Outputs](https://Milne-Group.github.io/SeqNado/outputs/)
- [CLI Reference](https://Milne-Group.github.io/SeqNado/cli/)
- [HPC Cluster Configuration](https://Milne-Group.github.io/SeqNado/cluster_config/)

## License

This project is licensed under [GPL3](LICENSE).
