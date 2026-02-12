[← Back to main page](index.md)

# Configuration

After successful genome configuration of SeqNado ([Genomes](genomes.md)).

## Using `seqnado config`

The `seqnado config` command builds a workflow configuration YAML file for the selected assay. If no assay is provided, the tool operates in multiomics mode.

### Assay Types

| Assay | CLI name |
|-------|----------|
| ATAC-seq | `atac` |
| ChIP-seq | `chip` |
| CRISPR analysis | `crispr` |
| CUT&Tag | `cat` |
| MCC | `mcc` |
| Methylation | `meth` |
| RNA-seq | `rna` |
| SNP analysis | `snp` |

For all available arguments and flags, see the CLI reference: [seqnado config](cli.md#cli-seqnado-config).

### Example Usage

#### Build a Configuration for RNA-seq
```bash
seqnado config rna --output rna_config.yaml
```

#### Multiomics Mode
```bash
seqnado config --make-dirs --interactive
```

### FASTQ Files

After generating the configuration and project directory using `seqnado config`, you need to link your FASTQ files into the `fastqs` directory. This ensures that the pipeline can locate and process your input data.

#### Symlinking FASTQ Files

Use the following command to create symbolic links for your FASTQ files:

```bash
ln -s /path/to/your/fastq/files/* <project_directory>/fastqs/
```

Replace `/path/to/your/fastq/files/` with the directory containing your FASTQ files and `<project_directory>` with the path to the project directory created by `seqnado config`.

#### Example

If your FASTQ files are located in `/data/fastq/` and your project directory is `rna_project`, run:

```bash
ln -s /data/fastq/* rna_project/fastqs/
```

This will create symbolic links to all FASTQ files in the `fastqs` directory of your project.

#### Safe Naming Strategies for FASTQ Files

To ensure compatibility with the pipeline and avoid errors, use consistent and descriptive naming conventions for your FASTQ files. Below are some examples of safe naming strategies:


- **ATAC-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  ```

- **ChIP-seq**:
  ```
  sample-name-rep1_Antibody_R1.fastq.gz
  sample-name-rep1_Antibody_R2.fastq.gz
  sample-name-rep2_Input_R1.fastq.gz
  sample-name-rep2_Input_R2.fastq.gz
  ```
  - `Antibody`: Name of the antibody used for ChIP.
  - `Input`: Control sample.

- **RNA-seq**:
  ```
  sample-name-rep1_R1.fastq.gz
  sample-name-rep1_R2.fastq.gz
  sample-name-rep2_R1.fastq.gz
  sample-name-rep2_R2.fastq.gz
  ```
  - `sample-name`: Unique identifier for the sample.
  - `rep1`, `rep2`: Biological or technical replicate number.
  - `R1`, `R2`: Read pair (forward and reverse).

Using these naming conventions ensures that the pipeline can correctly parse and process your data.

### Third Party Tools

SeqNado integrates with several third-party tools to enhance its functionality.
The configuration YAML file can be edited to set the parameters for each tool.

!!! warning "Review tool parameters before running"
    The defaults are sensible for most experiments, but some parameters should be adjusted to match your specific assay and library prep. Common examples:

    - **featureCounts (subread)**: The defaults count reads at the **exon** level (`-t exon`) and group them by **gene_id** (`-g gene_id`). Check that these match your GTF annotation — some GTFs use `gene_name` instead of `gene_id`, and if the attribute doesn't match, featureCounts will assign zero counts to everything. Paired-end mode (`-p --countReadPairs`) is auto-detected from your FASTQs. Strandedness should match your library prep — `0` for unstranded, `1` for forward-stranded, `2` for reverse-stranded (e.g., dUTP-based kits). Wrong strandedness will dramatically undercount genes.
    - **MACS2**: Use `--broad` for diffuse histone marks (H3K27me3, H3K36me3) and narrow mode (default) for transcription factors. The wrong mode will either miss broad domains or call noise as peaks.
    - **SEACR**: Adjust the threshold stringency — `stringent` mode is recommended for CUT&Tag to reduce false positives from low-background data, `relaxed` if you're missing expected peaks.
    - **deeptools (bamCoverage)**: Choose the normalisation method for BigWig generation — `RPGC` (reads per genomic content) is standard for most comparisons, `CPM` if you prefer a simpler normalisation.

    See the [Pipeline Overview](pipeline.md#supported-assays) for guidance on which tools to use for each assay type.

Below is a categorized list of supported tools:

#### Alignment Tools
- **bowtie2**: For aligning sequencing reads to a reference genome.
- **star**: A fast splice-aware aligner for RNA-seq data.

#### Processing Tools
- **samtools**: Utilities for manipulating alignments in the SAM/BAM format.
- **picard**: A set of Java command-line tools for manipulating high-throughput sequencing data.

#### Trimming Tools
- **cutadapt**: Removes adapter sequences from high-throughput sequencing reads.
- **trimgalore**: A wrapper tool around Cutadapt and FastQC for quality control and adapter trimming.

#### Peak Calling Tools
- **macs**: Model-based analysis of ChIP-seq data. Use narrow mode for transcription factors, broad mode for histone marks.
- **lanceotron**: Deep-learning peak caller, good for broad or weak peaks.
- **lanceotronmcc**: Specialized for MCC interaction peak calling.
- **seacr**: A peak caller designed for sparse, low-background data (CUT&Tag).

#### Quantification Tools
- **subread (featureCounts)**: Assigns aligned reads to genomic features (genes) for RNA-seq quantification. Ensure strandedness matches your library prep.
- **salmon**: Pseudoalignment-based RNA-seq quantification — faster than alignment-based methods.

#### Analysis Tools
- **deeptools**: For the analysis and visualization of high-throughput sequencing data (BigWig generation, heatmaps, metaplots).
- **homer**: A suite of tools for motif discovery and next-gen sequencing analysis.
- **bamnado**: A tool for BAM file manipulation and BigWig generation.

#### Specialized Tools
- **methyldackel**: For processing bisulfite sequencing data and extracting methylation calls.
- **bcftools**: Utilities for variant calling and manipulating VCF/BCF files.
- **mageck**: Statistical analysis of CRISPR screen data.

For more details on configuring these tools, refer to the individual tool documentation.

For configuration command options and usage patterns, see [seqnado config](cli.md#cli-seqnado-config).

## Next Steps

Once you have configured SeqNado, proceed to design your workflows:

[Design](design.md)