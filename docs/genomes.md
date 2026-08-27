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

### ⚠️ Dependencies

The genome build workflow requires **samtools**, **Bowtie2**, and **STAR**. These are not installed by default in a basic SeqNado environment.

**Use `--preset` to ensure dependencies are available:**

- **`--preset ls`** (recommended) — Uses Apptainer/Singularity to provide all tools
- **`--preset ss`** — On SLURM clusters with Apptainer/Singularity
- **`--preset lc`** — Local execution with Conda + Apptainer
- **`--preset le`** (default) — Requires manual installation of dependencies

All examples below include `--preset ls`. If using a different environment, adjust accordingly.

### Single genome

```bash
seqnado genomes build --name hg38 --outdir /path/to/genomes --preset ls
```

### Multiple genomes

Comma-separate names to build several genomes in one run:

```bash
seqnado genomes build --name hg38,mm39,dm6 --outdir /path/to/genomes --preset ls
```

### Spike-in (composite) genome

Combine a primary genome with a spike-in. This downloads both genomes, concatenates their FASTA and GTF, and builds composite indices:

```bash
seqnado genomes build --name hg38 --spikein dm6 --outdir /path/to/genomes --preset ls
```

The composite genome is named `hg38_dm6` and spike-in chromosomes are prefixed (e.g. `dm6_chr2L`).

### Custom FASTA

If your genome isn't on UCSC (a custom assembly, a non-model organism, or a
custom spike-in), point the build at your own FASTA/GTF instead of
downloading. `--name` still gives the genome an identifier — it's just used to
name the output directory and register it in `genome_config.json`, not to look
anything up on UCSC:

```bash
seqnado genomes build --name mygenome --fasta /path/to/custom.fa --gtf /path/to/custom.gtf --outdir /path/to/genomes --preset ls
```

`--gtf` is only required if you need the STAR index (e.g. for RNA-seq); the
blacklist download is always best-effort and falls back to an empty file when
nothing is found for the name.

The same flags exist for a custom spike-in: `--spikein-fasta` and
`--spikein-gtf` (require `--spikein`). Primary and spike-in genomes are
resolved independently, so you can mix a UCSC download with a custom source
on either side:

```bash
# UCSC primary genome + custom spike-in FASTA/GTF
seqnado genomes build --name hg38 --spikein myspikein \
    --spikein-fasta /path/to/spikein.fa --spikein-gtf /path/to/spikein.gtf \
    --outdir /path/to/genomes --preset ls

# Custom primary genome + UCSC spike-in
seqnado genomes build --name mygenome --fasta /path/to/custom.fa --gtf /path/to/custom.gtf \
    --spikein dm6 \
    --outdir /path/to/genomes --preset ls

# Both custom
seqnado genomes build --name mygenome --fasta /path/to/custom.fa --gtf /path/to/custom.gtf \
    --spikein myspikein --spikein-fasta /path/to/spikein.fa --spikein-gtf /path/to/spikein.gtf \
    --outdir /path/to/genomes --preset ls
```

Genome/spike-in names may contain letters, digits, and hyphens, but **not
underscores** — the composite genome name is built as `{primary}_{spikein}`,
so an underscore inside either name would make that ambiguous.

### Dry run

Preview planned jobs without executing them:

```bash
seqnado genomes build --name hg38 --outdir /path/to/genomes --preset ls --dry-run
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
| `--fasta`           | —               | Custom FASTA for the primary genome, instead of downloading from UCSC (single genome only) |
| `--gtf`             | —               | Custom GTF for the primary genome, instead of downloading from UCSC (single genome only) |
| `--spikein-fasta`   | —               | Custom FASTA for the `--spikein` genome, instead of downloading from UCSC |
| `--spikein-gtf`     | —               | Custom GTF for the `--spikein` genome, instead of downloading from UCSC |
| `--preset`          | `le`            | Snakemake profile preset (see below) |
| `--profile`         | —               | Path to a Snakemake profile directory (overrides --preset) |
| `--cores`, `-c`     | `4`             | Number of cores available to Snakemake; controls max parallelism and threads per rule |
| `--scale-resources` | `1.0`           | Scale memory/time requests |
| `--dry-run`         | off             | Preview the Snakemake DAG without executing |
| `--verbose`, `-v`   | off             | Print the full Snakemake command |

### Presets

Presets select a Snakemake execution profile. For genome builds, use a preset that includes samtools, Bowtie2, and STAR:

| Preset | Tools included? | Environment | Description |
|--------|-----------------|-------------|-------------|
| `le`   | ❌ No | Local | Requires manual tool installation (not recommended for genome builds) |
| `ls`   | ✅ Yes | Local | Singularity/Apptainer containers (recommended) |
| `lc`   | ✅ Yes | Local | Conda + Singularity |
| `ld`   | ✅ Yes | Local | Docker + Conda |
| `ss`   | ✅ Yes | SLURM cluster | Singularity/Apptainer containers on HPC |
| `t`    | ✅ Yes | Local | Testing/development |

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

**Key distinction:**
- `--contaminant-path` = **input** (where your contaminant reference indices are located)
- `--screen` = **output** (where to save the generated config file)

### Generating the config

```bash
seqnado genomes fastqscreen
```

This reads your genome config (`~/.config/seqnado/genome_config.json`), finds each genome's Bowtie2 index, and writes a `fastq_screen.conf` with a `DATABASE` entry per genome. Organism names (Human, Mouse, Drosophila, etc.) are inferred automatically from the genome prefix (e.g. `hg38` → Human, `mm39` → Mouse).

### Adding contaminant databases

By default, the command prompts for a path to contaminant reference files. If provided, it adds common contaminant screens (E. coli, PhiX, rRNA, adapters, etc.).

A pre-built set of contaminant references is available. Download and extract it:

```bash
wget https://userweb.molbiol.ox.ac.uk/public/project/milne_group/seqnado/genomes/fastqscreen_reference.tar.gz
tar -xzf fastqscreen_reference.tar.gz
```

Then generate the config with contaminants:

```bash
# Quick example: contaminants in current directory
seqnado genomes fastqscreen --contaminant-path ./fastqscreen_reference

# Specify custom output path
seqnado genomes fastqscreen --contaminant-path ./fastqscreen_reference --screen /path/to/my_fastq_screen.conf

# Skip contaminants entirely
seqnado genomes fastqscreen --no-contaminants
```

The contaminant directory should contain Bowtie2 indices organised in subdirectories:

```
fastqscreen_reference/
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
| `--screen`, `-s`      | `~/.config/seqnado/fastq_screen.conf` | **Output path**: Where to write the generated `fastq_screen.conf` |
| `--contaminant-path`  | *(prompted)* | **Input path**: Directory containing contaminant Bowtie2 indices |
| `--threads`, `-t`     | `8`     | Bowtie2 threads used by FastQ Screen |
| `--no-contaminants`   | off     | Skip contaminant databases |

---

**See Also:**

- [Configuration Guide](configuration.md) - Configure your analysis
- [CLI Reference](cli.md#cli-seqnado-genomes) - Complete genomes command options
- [Troubleshooting](troubleshooting.md#genome-setup-seqnado-genomes) - Genome setup issues
