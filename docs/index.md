# SeqNado

Welcome to SeqNado! SeqNado is a powerful bioinformatics tool designed to simplify the analysis and integration of high-throughput sequencing data. It supports a wide range of sequencing assays and provides advanced features for multiomics data processing.

## Key Features

- **Multiomics Data Processing**: Analyze and integrate data from multiple omics layers.
- **Configurable Workflows**: Predefined workflows for ATAC, ChIP, RNA, and more.
- **Third-Party Tool Integration**: Seamless integration with popular bioinformatics tools.
- **Advanced Data Processing**: Includes features like spike-in normalization, blacklist removal, and genome tiling.

## Supported Assays

SeqNado supports the following sequencing assays:

- **RNA-seq**: Transcriptome analysis.
- **ATAC-seq**: Chromatin accessibility profiling.
- **SNP Analysis**: Single nucleotide polymorphism detection.
- **ChIP-seq**: Protein-DNA interaction analysis.
- **CUT&Tag**: Epigenomic profiling.
- **Methylation**: DNA methylation analysis.
- **MCC**: Multiplex chromatin conformation capture.
- **CRISPR**: CRISPR screening analysis.

## Use Cases

SeqNado is ideal for:

- Researchers analyzing high-throughput sequencing data.
- Bioinformaticians integrating multiomics datasets.
- Labs requiring reproducible and scalable workflows.

## Additional Features

SeqNado includes the following advanced capabilities:

- **Spike-in Normalization**: Calculate normalization factors for spike-in controls.
- **Data Visualization**: Generate publication-ready plots and figures.
- **UCSC Hub Generation**: Create UCSC genome browser hubs for data sharing.
- **Genome Browser Plots**: Generate visualizations for genome-wide data, including UCSC genome browser hubs.
- **Quantification Methods**: Comprehensive tools for read count quantification, grouped read counts, and combined read counts.
- **Machine Learning Dataset Creation**: Prepare datasets for machine learning applications.

These features are fully customizable through configuration files, making SeqNado adaptable to a variety of research needs.

# Get Started

Follow the step-by-step guide to get up and running:

1. **[Installation](installation.md)**: Set up the SeqNado environment.
2. **[Initialisation](initialisation.md)**: Configure your local environment (CLI: [seqnado init](cli.md#cli-seqnado-init)).
3. **[Genome Setup](genomes.md)**: Manage reference genomes and indexes (CLI: [seqnado genomes](cli.md#cli-seqnado-genomes)).
4. **[Configuration](configuration.md)**: Define your experiment parameters (CLI: [seqnado config](cli.md#cli-seqnado-config)).
5. **[Design Guide](design.md)**: Create your sample metadata design (CLI: [seqnado design](cli.md#cli-seqnado-design)).
6. **[Pipeline Overview](pipeline.md)**: Run the workflow (CLI: [seqnado pipeline](cli.md#cli-seqnado-pipeline)) and explore [Outputs](outputs.md).
7. **[Output Examples](examples.md)**: See real example outputs from test data runs.

For a quick end-to-end example, see the **[Quick Start](quick_start.md)** guide.
