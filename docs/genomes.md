[← Back to main page](index.md)

# Genome Setup

After successful initialization of SeqNado ([Initialisation](initialisation.md)), use `seqnado genomes` to manage reference genomes.

## Subcommands

| Subcommand     | Description |
|----------------|-------------|
| `list`         | Show available genome configurations |
| `edit`         | Open the user genome config in `$EDITOR` |
| `build`        | Download a genome from UCSC and build indices via Snakemake |
| `fastqscreen`  | Generate a FastqScreen configuration file |

## Building a Genome

`seqnado genomes build` downloads the FASTA, GTF, chromosome sizes, and blacklist from UCSC, then builds Bowtie2 and STAR indices.

### Single genome

```bash
seqnado genomes build --name hg38 --outdir /path/to/genomes
```

### Multiple genomes

Comma-separate names to build several genomes in one run:

```bash
seqnado genomes build --name hg38,mm39,dm6 --outdir /path/to/genomes
```

### Spike-in (composite) genome

Combine a primary genome with a spike-in. This downloads both genomes, concatenates their FASTA and GTF, and builds composite indices:

```bash
seqnado genomes build --name hg38 --spikein dm6 --outdir /path/to/genomes
```

The composite genome is named `hg38_dm6` and spike-in chromosomes are prefixed (e.g. `dm6_chr2L`).

### Dry run

Preview planned jobs without executing them:

```bash
seqnado genomes build --name hg38 --outdir /path/to/genomes --dry-run
```

### What gets built

For each genome, the workflow produces the following output structure:

```
<outdir>/<genome>/
├── sequence/
│   ├── <genome>.fa          # FASTA (downloaded from UCSC)
│   ├── <genome>.fa.fai      # samtools faidx index
│   └── <genome>.chrom.sizes # chromosome sizes
├── genes/
│   └── <genome>.ncbiRefSeq.gtf  # gene annotations
├── bt2_index/
│   └── <genome>.*.bt2      # Bowtie2 index
├── STAR_2.7.10b/            # STAR index directory
└── <genome>-blacklist.bed.gz  # ENCODE blacklist regions
```

### Genome config auto-update

On successful completion, the build automatically registers the genome in `~/.config/seqnado/genome_config.json`. This makes it immediately available for pipeline runs — no manual editing needed.

### Build options

| Option              | Default         | Description |
|---------------------|-----------------|-------------|
| `--name`, `-n`      | *(required)*    | Genome name(s), comma-separated |
| `--outdir`, `-o`    | `./genome_build`| Output directory |
| `--spikein`, `-sp`  | —               | Spike-in genome name for composite builds |
| `--preset`          | `le`            | Snakemake profile preset (see below) |
| `--cores`, `-c`     | `4`             | Number of Snakemake cores |
| `--scale-resources` | `1.0`           | Scale memory/time requests |
| `--dry-run`         | off             | Preview the Snakemake DAG without executing |
| `--verbose`, `-v`   | off             | Print the full Snakemake command |

### Presets

Presets select a Snakemake execution profile:

| Preset | Profile                | Description |
|--------|------------------------|-------------|
| `le`   | `local_environment`    | Local execution, no containers (default) |
| `ls`   | `local_singularity`    | Local execution with Apptainer/Singularity |
| `lc`   | `local_conda`          | Local execution with Conda + Apptainer |
| `ss`   | `slurm_singularity`    | SLURM cluster with Apptainer |

## Listing Genomes

Show all configured genomes and their paths:

```bash
seqnado genomes list
```

## Editing Genome Config

Open `~/.config/seqnado/genome_config.json` in your editor:

```bash
seqnado genomes edit
```

Set `$EDITOR` to your preferred editor (defaults to `nano`).

## FastqScreen Configuration

[FastQ Screen](https://www.bioinformatics.babraham.ac.uk/projects/fastq_screen/) checks for sample contamination by aligning reads against a set of reference genomes. SeqNado can auto-generate the `fastq_screen.conf` configuration file from your built genomes.

### Generating the config

```bash
seqnado genomes fastqscreen
```

This reads your genome config (`~/.config/seqnado/genome_config.json`), finds each genome's Bowtie2 index, and writes a `fastq_screen.conf` with a `DATABASE` entry per genome. Organism names (Human, Mouse, Drosophila, etc.) are inferred automatically from the genome prefix (e.g. `hg38` → Human, `mm39` → Mouse).

### Adding contaminant databases

By default, the command prompts for a path to contaminant reference files. If provided, it adds common contaminant screens (E. coli, PhiX, rRNA, adapters, etc.):

```bash
# Provide path directly
seqnado genomes fastqscreen --contaminant-path /path/to/contaminant_refs

# Skip contaminants entirely
seqnado genomes fastqscreen --no-contaminants
```

The contaminant directory should contain Bowtie2 indices organised in subdirectories:

```
contaminant_refs/
├── E_coli/Ecoli.*bt2
├── PhiX/phi_plus_SNPs.*bt2
├── rRNA/GRCm38_rRNA.*bt2
├── Vectors/Vectors.*bt2
├── Adapters/Contaminants.*bt2
└── ...
```

### Options

| Option                | Default | Description |
|-----------------------|---------|-------------|
| `--screen`, `-s`      | `~/.config/seqnado/fastq_screen.conf` | Output path |
| `--threads`, `-t`     | `8`     | Bowtie2 threads used by FastQ Screen |
| `--no-contaminants`   | off     | Skip contaminant databases |
| `--contaminant-path`  | *(prompted)* | Path to contaminant Bowtie2 indices |

---

**See Also:**

- [Configuration Guide](configuration.md) - Configure your analysis
- [CLI Reference](cli.md#cli-seqnado-genomes) - Complete genomes command options
- [Troubleshooting](troubleshooting.md#genome-setup-seqnado-genomes) - Genome setup issues
