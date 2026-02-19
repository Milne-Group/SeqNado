[← Back to main page](index.md)

# CLI Reference

This page provides a complete reference of all available SeqNado commands and options.

For quick examples and typical usage patterns, see the [Quick Start Guide](quick_start.md).

---

# `seqnado`

**SeqNado CLI**

Initialize your environment, build configs, create design files, and run pipelines.
Use --help on any subcommand for details.

**Usage**:

```console
$ seqnado [OPTIONS] COMMAND [ARGS]...
```

**Options**:

* `-v, --version`: Show version and exit.
* `--install-completion`: Install completion for the current shell.
* `--show-completion`: Show completion for the current shell, to copy it or customize the installation.
* `--help`: Show this message and exit.

**Commands**:

* `init`: Initialize SeqNado user environment.
* `genomes`: Manage genome configurations (list, edit, build, fastqscreen)
* `config`: Build a workflow configuration YAML
* `download`: Download FASTQ files from GEO/SRA
* `design`: Generate a SeqNado design CSV from FASTQ files
* `pipeline`: Run the data processing pipeline
* `tools`: Explore bioinformatics tools available in SeqNado

## `seqnado init` {#cli-seqnado-init}

Initialize SeqNado user environment.

- Logs the current Conda environment if active (optional).
- Runs packaged Apptainer/Singularity init (if `apptainer` on PATH).
- Ensures ~/.config/seqnado/genome_config.json exists (template or preset).

**Usage**:

```console
$ seqnado init [OPTIONS]
```

**Options**:

* `--preset / --no-preset`: Use packaged preset genomes instead of the editable template.  [default: no-preset]
* `--dry-run / --no-dry-run`: Show actions without writing files or running scripts.  [default: no-dry-run]
* `-v, --verbose`: Increase logging verbosity.
* `--help`: Show this message and exit.

## `seqnado genomes` {#cli-seqnado-genomes}

Manage genome configurations (list, edit, build, or generate fastq-screen config)

**Usage**:

```console
$ seqnado genomes [OPTIONS] COMMAND [ARGS]...
```

**Commands**:

- **list**: Show packaged and user genome presets
- **edit**: Open user genome config in $EDITOR
- **build**: Download genome and build indices via Snakemake
- **fastqscreen**: Generate FastqScreen configuration file

**Options**:

* `--help`: Show this message and exit.

### `seqnado genomes list`

Show packaged and user genome presets.

**Usage**:

```console
$ seqnado genomes list [OPTIONS] [ASSAY]
```

**Arguments**:

- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics` (default: `atac`)

**Options**:

* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

**Example**:

```bash
seqnado genomes list atac
seqnado genomes list rna
```

### `seqnado genomes edit`

Open user genome config in $EDITOR.

**Usage**:

```console
$ seqnado genomes edit [OPTIONS]
```

**Options**:

* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

**Example**:

```bash
seqnado genomes edit
```

### `seqnado genomes build`

Download genome and build indices via Snakemake.

**Usage**:

```console
$ seqnado genomes build [OPTIONS]
```

**Options**:

* `-n, --name TEXT`: Genome name(s), comma-separated for multiple (e.g., hg38 or hg38,mm39,dm6) **[required]**
* `-o, --outdir PATH`: Output directory for build  [default: genome_build]
* `-sp, --spikein TEXT`: Spike-in genome name for composite builds (e.g., mm39)
* `--preset [lc|le|ls|ld|ss|t]`: Snakemake job profile preset.  [default: le]
* `--profile PATH` / `--profiles PATH`: Path to a Snakemake profile directory (overrides --preset).
* `-c, --cores INTEGER`: Number of Snakemake cores.  [default: 4]
* `--scale-resources FLOAT`: Scale memory/time.  [default: 1.0]
* `--dry-run`: Preview the Snakemake DAG without executing.
* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

**Examples**:

```bash
# Build single genome
seqnado genomes build --name hg38 --outdir /path/to/genomes

# Build multiple genomes
seqnado genomes build --name hg38,mm39,dm6 --outdir /path/to/genomes

# Build with spike-in
seqnado genomes build --name hg38 --spikein dm6 --outdir /path/to/genomes

# Dry run
seqnado genomes build --name hg38 --outdir /path/to/genomes --dry-run
```

### `seqnado genomes fastqscreen`

Generate FastqScreen configuration file.

**Usage**:

```console
$ seqnado genomes fastqscreen [OPTIONS]
```

**Options**:

* `-s, --screen PATH`: Output path for fastqscreen config  [default: ~/.config/seqnado/fastq_screen.conf]
* `-t, --threads INTEGER`: Number of threads for Bowtie2  [default: 8]
* `--no-contaminants`: Exclude contaminant databases
* `--contaminant-path PATH`: Path to contaminant reference files
* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

**Example**:

```bash
seqnado genomes fastqscreen --screen /path/to/fastq_screen.conf --threads 16
```

## `seqnado config` {#cli-seqnado-config}

Build a workflow configuration YAML for the selected ASSAY. If no assay is provided, multiomics mode is used.

**Usage**:

```console
$ seqnado config [OPTIONS] [ASSAY]
```

**Arguments**:

- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics`; omitted → multiomics mode

**Options**:

* `--make-dirs / --no-make-dirs`: Create/don't create the output project directory or fastq subdir.  [default: make-dirs]
* `--render-options / --no-render-options`: Render all options (even if not used by the workflow).  [default: no-render-options]
* `-o, --output PATH`: Explicit path for the rendered config file.
* `-v, --verbose`: Increase logging verbosity.
* `--interactive / --no-interactive`: Interactively prompt for config values. Non-interactive mode only works for single assay configs (except MCC and multiomics).  [default: interactive]
* `--help`: Show this message and exit.

## `seqnado download` {#cli-seqnado-download}

Download FASTQ files from GEO/SRA using a metadata TSV file.

**Usage**:

```console
$ seqnado download [OPTIONS] METADATA_TSV
```

**Arguments**:

- **METADATA_TSV**: Path to TSV file with GEO/SRA metadata (required)

**Required TSV Columns:**

- `run_accession`: SRA run ID (e.g., SRR123456)
- `sample_title`: Sample name for files
- `library_name`: Library name
- `library_layout`: Either "PAIRED" or "SINGLE"

**Options**:

* `-o, --outdir PATH`: Output directory for downloaded FASTQ files.  [default: fastqs]
* `-a, --assay [rna|atac|chip|cat|snp|meth|mcc|crispr]`: Assay type for auto-generating design file after download.
* `-d, --design-output PATH`: Path for generated design file (only with --assay).  [default: metadata_{assay}.csv]
* `-c, --cores INTEGER`: Number of parallel download jobs.  [default: 4]
* `--preset [lc|le|ls|ld|ss|t]`: Snakemake job profile preset.  [default: le]
* `--profile PATH` / `--profiles PATH`: Path to a Snakemake profile directory (overrides --preset).
* `-n, --dry-run`: Show what would be run without executing downloads.
* `-v, --verbose`: Increase logging verbosity.
* `--help`: Show this message and exit.

**Example**:

```bash
# Download paired-end RNA-seq data and generate design file
seqnado download metadata.tsv -o geo_fastqs -a rna --cores 4

# Download with SLURM preset for HPC
seqnado download metadata.tsv --preset ss --cores 8
```

→ [Complete download guide](geo_download.md)

## `seqnado design` {#cli-seqnado-design}

Generate a SeqNado design CSV from FASTQ files for ASSAY. If no assay is provided, multiomics mode is used.

**Usage**:

```console
$ seqnado design [OPTIONS] [ASSAY] [FASTQ ...]
```

**Arguments**:

- **ASSAY**: `rna` | `atac` | `snp` | `chip` | `cat` | `meth` | `mcc` | `crispr` | `multiomics`; omitted → multiomics mode
- **FASTQ**: one or more FASTQ files

**Options**:

* `-o, --output PATH`: Output CSV filename (default: metadata_{assay}.csv).
* `--ip-to-control TEXT`: List of antibody,control pairings for IP assays (e.g. ChIP). Format: 'antibody1:control1,antibody2:control2'. If provided will assign a control with a specified name to that ip in the metadata. If not provided, controls will be assigned based on a best-effort matching of sample names.
* `--group-by`: Group samples by a regular expression or a column.
* `--auto-discover / --no-auto-discover`: Search common folders if none provided.  [default: auto-discover]
* `--interactive / --no-interactive`: Interactively offer to add missing columns using schema defaults.  [default: interactive]
* `--accept-all-defaults`: Non-interactive: auto-add only columns that have a schema default.
* `--deseq2-pattern TEXT`: Regex pattern to extract DESeq2 groups from sample names. First capture group will be used. Example: r'-(\w+)-rep' for 'sample-GROUP-rep1'
* `-v, --verbose`: Increase logging verbosity.
* `--help`: Show this message and exit.

## `seqnado pipeline` {#cli-seqnado-pipeline}

Run the data processing pipeline for ASSAY (Snakemake under the hood). Any additional arguments are passed to Snakemake (e.g., `seqnado pipeline rna -n` for dry-run, `--unlock`, etc.).

**Usage**:

```console
$ seqnado pipeline [OPTIONS] [ASSAY]
```

**Arguments**:

- **ASSAY**: required for single-assay; optional in multiomics mode

**Options**:

* `--configfile PATH`: Path to a SeqNado config YAML (default: config_<ASSAY>.yaml).
* `--version`: Print SeqNado version and exit.
* `--preset [lc|le|ls|ld|ss|t]`: Snakemake job profile preset.  [default: le]
* `--profile PATH` / `--profiles PATH`: Path to a Snakemake profile directory (overrides --preset).
* `--clean-symlinks / --no-clean-symlinks`: Remove symlinks created by previous runs.  [default: no-clean-symlinks]
* `-s, --scale-resources FLOAT`: Scale memory/time (env: SCALE_RESOURCES).  [default: 1.0]
* `-v, --verbose`: Increase logging verbosity.
* `-q, --queue TEXT`: Slurm queue/partition for the `ss` preset.  [default: short]
* `--print-cmd`: Print the Snakemake command before running it.
* `-c, --cores INTEGER`: Number of CPU cores for Snakemake to use.  [default: 1]
* `-n, --dry-run`: Snakemake dry-run: show what would be executed without running.
* `--unlock`: Unlock the Snakemake working directory (use after a failed/interrupted run).
* `--rerun-incomplete`: Re-run jobs left incomplete from a previous run.
* `--help`: Show this message and exit.

Any additional arguments not listed above are passed directly to Snakemake.

## `seqnado tools` {#cli-seqnado-tools}

Explore bioinformatics tools available in SeqNado. List all tools, filter by category, get version information, and view tool help/options from containers.

**Usage**:

```console
$ seqnado tools [OPTIONS] [TOOL]
```

**Arguments**:

- **TOOL**: Optional tool name to get details about a specific tool (e.g., `fastqc`, `bowtie2`, `samtools`)

**Options**:

* `-l, --list`: List all available tools in the pipeline
* `-c, --category [TEXT]`: Filter tools by category (name or number). Omit value for interactive selection
* `--options`: Show tool help/options from container (requires tool argument and apptainer)
* `--citation`: Show the BibTeX citation for a tool (requires tool argument)
* `-s, --subcommand TEXT`: Specify subcommand for tools that have multiple functions
* `-v, --verbose`: Increase logging verbosity
* `--help`: Show this message and exit.

**Tool Categories**:

Tools are organized into the following categories, matching the pipeline workflow:

1. **Download** - Data acquisition (e.g., `fasterq-dump`)
2. **Quality Control** - QC and reporting (e.g., `fastqc`, `fastq-screen`, `qualimap`)
3. **Preprocessing** - Read trimming and filtering (e.g., `trim-galore`, `cutadapt`)
4. **Alignment** - Read mapping (e.g., `bowtie2`, `STAR`, `minimap2`)
5. **Analysis** - Peak calling, variant analysis (e.g., `macs`, `homer`, `bedtools`)
6. **Visualization** - Generating tracks and plots (e.g., `deeptools`, `plotnado`)
7. **Reporting** - Report generation (e.g., `multiqc`, `quarto`)
8. **Quantification** - Expression quantification (e.g., `salmon`, `featureCounts`)
9. **Utilities** - Helper tools (e.g., `samtools`, `snakemake`, `pigz`)

### Examples

**List all available tools:**

```bash
seqnado tools --list
```

Output shows all tools organized by category with descriptions and version information.

**Interactive category selection:**

```bash
seqnado tools -c
# or
seqnado tools --category
```

Displays a numbered list of categories and prompts you to select one.

**Filter by category name:**

```bash
# Show all alignment tools
seqnado tools --category "Alignment"

# Case-insensitive matching
seqnado tools --category alignment
```

**Filter by category number:**

```bash
# Show category 4 (Alignment)
seqnado tools --category 4
```

**Get details about a specific tool:**

```bash
# Show information about FastQC
seqnado tools fastqc
```

Shows the tool's description, category, command, and version (detected from container or local installation).

**View tool help from container:**

```bash
# Get FastQC help
seqnado tools fastqc --options

# Get help for a deeptools subcommand
seqnado tools deeptools --subcommand plotHeatmap --options
```

This runs the tool in the SeqNado container and displays its help output, useful for:

- Checking available command-line options
- Understanding tool-specific parameters
- Verifying tool behavior in the container environment

**Tools with subcommands:**

Some tools like `deeptools`, `bamnado`, and `homer` expose multiple subcommands:

```bash
# Get deeptools bamCoverage help
seqnado tools deeptools --subcommand bamCoverage --options

# Get bamnado bam-coverage help
seqnado tools bamnado --subcommand bam-coverage --options

# Get HOMER makeTagDirectory help
seqnado tools homer --subcommand makeTagDirectory --options
```

### Version Detection

The `tools` command attempts to detect tool versions from:

1. **Container environment** (if Apptainer/Singularity is available)
2. **Local installation** (if tool is in PATH)

Version information helps ensure reproducibility and track tool updates across pipeline runs.

### Use Cases

- **Pipeline exploration**: Discover which tools are available for your analysis
- **Tool documentation**: Quick access to tool help without running containers manually
- **Version tracking**: Check which tool versions are used in the pipeline
- **Troubleshooting**: Verify tool availability and configuration
- **Learning**: Explore bioinformatics tools organized by analysis workflow
