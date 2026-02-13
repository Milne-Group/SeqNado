[← Back to main page](index.md)

# Troubleshooting

This page covers common errors you may encounter at each stage of the SeqNado workflow, along with their causes and fixes.

!!! tip
    For HPC/cluster-specific issues (SLURM, resource limits, container downloads), see the [HPC Clusters](cluster_config.md#troubleshooting) guide.

---

## Installation

### Conda/Mamba cannot find the `seqnado` package

```
PackagesNotFoundError: The following packages are not available from current channels: seqnado
```

**Cause:** The bioconda channel is not configured.

**Fix:** Add the required channels and retry:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
mamba create -n seqnado seqnado
```

### Dependency conflicts during installation

```
LibMambaUnsatisfiableError: Encountered problems while solving...
```

**Cause:** An existing environment has conflicting packages.

**Fix:** Create a fresh environment rather than installing into an existing one:

```bash
mamba create -n seqnado -c bioconda seqnado
```

If the problem persists, try installing with `pip` in a clean environment:

```bash
mamba create -n seqnado python=3.12
mamba activate seqnado
pip install seqnado
```

---

## Initialisation (`seqnado init`)

### `apptainer` / `singularity` not found

```
`apptainer` not found on PATH.
```

**Cause:** Apptainer (or Singularity) is not installed or not loaded.

**Fix:** On most HPC systems, you need to load the module first:

```bash
module load apptainer
# or
module load singularity
```

If neither is available, ask your cluster admin to install Apptainer, or use a local environment preset (`le`) which does not require containers:

```bash
seqnado pipeline atac --preset le
```

### Failed to pull singularity image — remote has no library client

```
FATAL: Unable to get library client configuration:
remote has no library client
```

**Cause:** The default Apptainer remote endpoint is not configured.

**Fix:** Re-run `seqnado init`, or manually add the SylabsCloud remote:

```bash
apptainer remote add --no-login SylabsCloud cloud.sylabs.io
apptainer remote use SylabsCloud
```

Then re-run `seqnado init`.

---

## Genome Setup (`seqnado genomes`)

### No Bowtie2 index files found

```
ValueError: No Bowtie2 index files found for prefix '...'
```

**Cause:** The genome was not built yet, or the path in the genome config points to the wrong location.

**Fix:** Check your genome configuration and rebuild if needed:

```bash
seqnado genomes list           # check what's configured
seqnado genomes edit           # fix paths if wrong
seqnado genomes build atac \
  --fasta /path/to/genome.fa \
  --name hg38 \
  --outdir /path/to/output     # rebuild if missing
```

### STAR index directory does not exist

```
ValueError: The directory ... does not exist or is not a directory
```

**Cause:** The STAR index path in the genome config is incorrect or the index has not been built.

**Fix:** Build the genome with the `rna` assay type, which generates the STAR index:

```bash
seqnado genomes build rna --fasta /path/to/genome.fa --name hg38 --outdir /path/to/output
```

Then verify the path with `seqnado genomes edit`.

### Genome build runs out of memory

**Cause:** STAR index building requires significant RAM (~32 GB for the human genome).

**Fix:** Run the build on a node with enough memory. On an HPC cluster, request an interactive session with sufficient resources before building:

```bash
srun --mem=40G --cpus-per-task=8 --time=2:00:00 --pty bash
mamba activate seqnado
seqnado genomes build rna --fasta genome.fa --name hg38 --outdir /path/to/output --threads 8
```

---

## Configuration (`seqnado config`)

### YAML syntax errors

```
yaml.scanner.ScannerError: while scanning ...
```

**Cause:** The YAML config file has a formatting error. YAML is sensitive to indentation and special characters.

**Fix:** Common YAML pitfalls to check:

- Use **spaces, not tabs** for indentation
- Keep indentation consistent (2 spaces per level)
- Strings containing colons or special characters need quoting: `name: "my:sample"`
- Boolean values are `true`/`false` (lowercase)

!!! tip
    Run `seqnado config` interactively to regenerate the config file rather than editing YAML by hand. This avoids most syntax issues:
    ```bash
    seqnado config atac --interactive
    ```

### Invalid strandedness value

```
ValueError: strandedness must be 0 (unstranded), 1 (forward), or 2 (reverse).
```

**Cause:** The strandedness field in the config has an invalid value.

**Fix:** Set it to one of the accepted values:

- `0` — unstranded (most common for standard RNA-seq library preps)
- `1` — forward stranded
- `2` — reverse stranded (common for dUTP-based library preps)

If unsure, check with whoever prepared the library, or use `0` as a starting point and examine the infer_experiment results in the QC report.

### Project name contains invalid characters

```
ValueError: Name contains invalid characters. Use alphanumerics, hyphens, and underscores.
```

**Cause:** The project name contains spaces or special characters.

**Fix:** Use only letters, numbers, hyphens (`-`), and underscores (`_`):

```yaml
# Bad
name: My ChIP Experiment (2024)

# Good
name: my-chip-experiment-2024
```

---

## Design Files (`seqnado design`)

### Sample ID contains invalid characters

```
ValueError: sample_id must match ^[a-zA-Z0-9_-]+$
```

**Cause:** Sample names in your FASTQ filenames or design CSV contain spaces, dots, or other special characters.

**Fix:** Rename your FASTQ files to use only letters, numbers, hyphens, and underscores. For example:

```
# Bad
Sample 1.Rep1_R1.fastq.gz
sample.name_R1.fastq.gz

# Good
Sample1-Rep1_R1.fastq.gz
sample_name_R1.fastq.gz
```

### Duplicate sample IDs

```
ValueError: sample_id values are not unique
```

**Cause:** Two or more rows in the design CSV have the same sample ID, or multiple FASTQ file pairs resolve to the same sample name.

**Fix:** Check your design CSV for duplicate entries. For IP-based assays (ChIP-seq, CUT&Tag), duplicates are allowed only if the `sample_id` + `ip` combination is unique (e.g., same sample with different antibodies).

### Invalid number of FASTQ files for sample

```
ValueError: Invalid number of FASTQ files for sample_name: 3. Expected 1 or 2 files.
```

**Cause:** SeqNado found more than 2 FASTQ files matching the same sample name. This usually happens when file naming is inconsistent or extra files are in the `fastqs/` directory.

**Fix:** Check the `fastqs/` directory for unexpected files:

```bash
ls fastqs/ | sort
```

Each sample should have either 1 file (single-end) or exactly 2 files (paired-end, `_R1` and `_R2`). Remove or move any extra files.

### No FASTQ files found

```
FileNotFoundError: No FASTQ files found in ...
```

**Cause:** The `fastqs/` directory is empty or the files don't have a recognised extension (`.fastq.gz`, `.fq.gz`).

**Fix:** Check that your symlinks are valid and point to actual files:

```bash
ls -la fastqs/
```

If the symlinks are broken (shown in red), recreate them with the correct source path:

```bash
ln -s /correct/path/to/fastq/files/* fastqs/
```

### Antibody must be specified for ChIP/CUT&Tag assays

```
ValueError: Antibody must be specified for ChIP-seq assays.
```

**Cause:** The design file for a ChIP-seq or CUT&Tag run is missing the `ip` column, or it contains empty values.

**Fix:** Ensure your design CSV has an `ip` column specifying the antibody or "input" for each sample. See the [Design Guide](design.md) for examples.

### IP control pairing errors

```
ValueError: Multiple control samples matched ..., but no manual mapping provided
```

**Cause:** SeqNado found more than one "input" sample and cannot automatically determine which control belongs to which IP sample.

**Fix:** Explicitly map controls in your design CSV by ensuring each IP sample has a unique control, or add a `control` column to specify the pairing manually.

---

## Pipeline Execution (`seqnado pipeline`)

### `snakemake` not found

```
`snakemake` not found on PATH. Install/activate the environment that provides it.
```

**Cause:** The SeqNado conda environment is not activated.

**Fix:**

```bash
mamba activate seqnado
```

### Config file not found

```
Workflow defines configfile config_atac.yaml but it is not present or accessible.
```

**Cause:** You are running the pipeline from the wrong directory, or `seqnado config` was not run first.

**Fix:** Make sure you are in the project directory that contains your config file:

```bash
ls config_*.yaml   # check if config exists in current directory
```

If not, either `cd` into the correct project directory, or point to the config explicitly:

```bash
seqnado pipeline atac --configfile /path/to/config_atac.yaml
```

### Snakemake rule fails with a cryptic error

**Cause:** A specific step in the pipeline failed. Snakemake error output can be verbose and hard to parse.

**Fix:** Look for the actual error by scrolling up past the Snakemake traceback. The key information is usually in these lines:

1. **The rule name** — tells you which step failed (e.g., `rule align_bowtie2`)
2. **The log file path** — Snakemake prints `log: seqnado_output/.../logs/...` which contains the tool's actual error output
3. **The return code** — a non-zero exit code from the underlying tool

Read the log file for the specific error:

```bash
cat seqnado_output/atac/logs/<rule_name>/<sample>.log
```

!!! tip
    Use `--verbose` and `--print-cmd` to get more diagnostic information:
    ```bash
    seqnado pipeline atac --preset le --verbose --print-cmd
    ```

### Pipeline killed — out of memory

```
slurmstepd: error: Detected 1 oom_kill event in StepId=...
```

**Cause:** A pipeline step exceeded the allocated memory on the cluster.

**Fix:** Increase the resource scaling factor:

```bash
seqnado pipeline atac --preset ss --scale-resources 2.0
```

For persistent memory issues, see the [HPC Clusters troubleshooting guide](cluster_config.md#troubleshooting).

### Pipeline hangs or seems stuck

**Cause:** This is usually normal — some steps (alignment, peak calling) take a long time on large datasets.

**Fix:** Check if jobs are actually running:

```bash
# On SLURM clusters
squeue -u $USER

# For local execution, check CPU usage
top
```

If jobs are queued but not starting, your cluster partition may be busy. Consider switching to a less-busy queue with `--queue`.

---

## GEO/SRA Downloads (`seqnado download`)

### Missing required columns in TSV

```
Missing required columns in TSV: run_accession, sample_title
```

**Cause:** The metadata TSV downloaded from ENA is missing expected columns.

**Fix:** Ensure your TSV file contains these columns: `run_accession`, `sample_title`, `library_name`, and ideally `library_layout`. Download the file from the [ENA Browser](https://www.ebi.ac.uk/ena/browser/) using the "Download report" option with the correct column selection.

### Unknown library_layout value

```
Unknown library_layout '...' for sample_name
```

**Cause:** The `library_layout` column contains a value other than `PAIRED` or `SINGLE`.

**Fix:** Edit the TSV to correct the layout values. Valid values are `PAIRED` or `SINGLE` (case-insensitive). If unknown, check the GEO page for the experiment to determine the sequencing type.

### Downloads fail or are very slow

**Cause:** Network issues, or the SRA servers are under load.

**Fix:** Try reducing the number of parallel downloads and retrying:

```bash
seqnado download metadata.tsv --cores 2
```

If downloads consistently fail, check that you have internet access from your compute environment and that no firewall rules are blocking SRA/ENA traffic.

---

## General Tips

### Enable verbose logging

Add `--verbose` (or `-v`) to any SeqNado command to see detailed log output. This is the single most useful thing you can do when debugging:

```bash
seqnado pipeline atac --preset le --verbose
```

### Print the underlying Snakemake command

Use `--print-cmd` to see exactly what Snakemake command SeqNado is running. This is helpful when asking for help or reporting bugs:

```bash
seqnado pipeline atac --preset le --print-cmd
```

### Check Snakemake log files

Pipeline logs are stored in the output directory. Each rule writes its own log:

```
seqnado_output/<assay>/logs/<rule_name>/<sample>.log
```

These log files contain the actual output from the underlying tools (Bowtie2, MACS2, DESeq2, etc.) and are usually more informative than the Snakemake error summary.

### Resume a failed pipeline run

Snakemake automatically tracks completed steps. If a run fails partway through, simply re-run the same command — it will pick up where it left off:

```bash
seqnado pipeline atac --preset le
```

### Dry run to check your setup

Use Snakemake's dry-run mode to verify that your config, design, and files are correct without actually running anything:

```bash
seqnado pipeline atac --preset le -- --dry-run
```

This shows you what steps would be executed and can catch configuration errors early.

---

## Still stuck?

1. Check the [FAQ](faq.md) for other common questions
2. Open an issue on [GitHub](https://github.com/alsmith151/SeqNado/issues) with the error message and your `--verbose --print-cmd` output
