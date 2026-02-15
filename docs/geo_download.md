# GEO/SRA Data Download

This document describes how to use SeqNado to download data from GEO/SRA.

## Quick Start

### 1. Get metadata from GEO/ENA

First, obtain the metadata TSV file from GEO or ENA. For example:
- Go to [ENA Browser](https://www.ebi.ac.uk/ena/browser/)
- Search for your project (e.g., PRJNA1234567)
- Download the "TSV" report with run information

**Required columns:**
- `run_accession` (e.g., SRR123456)
- `sample_title` (sample name)
- `library_name` (e.g., GSM identifier)
- `library_layout` (PAIRED or SINGLE) - **Required for download**

The `library_layout` column is essential for determining whether to create one or two FASTQ files per sample.

It is advisable to alter the `sample_title` column to use seqnado naming conventions, removing any spaces. 
See [Configuration](configuration.md#safe-naming-strategies-for-fastq-files) for more details.

### 2. Download FASTQs

Run the download command:

```bash
seqnado download filereport_read_run_PRJNA1234567.tsv \
    --outdir fastqs \
    --cores 8
```

This will:
- Parse the metadata TSV
- Download all FASTQ files using prefetch/fasterq-dump
- Retry failed downloads with scaled resources on each attempt
- Compress and rename files to: `{GSM}-{sample}_R1.fastq.gz`

### 3. Generate design file (optional)

To also generate a SeqNado design file:

```bash
seqnado download filereport_read_run_PRJNA1234567.tsv \
    --outdir fastqs \
    --assay rna \
    --design-output metadata_rna.csv \
    --cores 8
```

Or run separately after download:

```bash
seqnado design rna --output metadata_rna.csv fastqs/*.fastq.gz
```

## Command Reference

```
seqnado download [OPTIONS] METADATA_TSV
```

### Required Arguments

- `METADATA_TSV`: Path to TSV file from GEO/ENA with run information

### Options

- `-o, --outdir PATH`: Output directory for FASTQ files (default: fastqs)
- `-a, --assay TEXT`: Assay type for design file generation (rna, atac, chip, etc.)
- `-d, --design-output PATH`: Output path for design CSV
- `-c, --cores INT`: Number of parallel downloads (default: 4)
- `--preset TEXT`: Snakemake profile preset (le/lsf/ss) (default: le)
- `--profile PATH` / `--profiles PATH`: Path to a Snakemake profile directory (overrides --preset)
- `-n, --dry-run`: Show what would be downloaded without downloading
- `-v, --verbose`: Increase logging verbosity

## Examples

### Download only

```bash
seqnado download filereport.tsv --outdir fastqs -c 10
```

### Download + generate RNA-seq design

```bash
seqnado download filereport.tsv \
    --outdir fastqs \
    --assay rna \
    -c 8
```

### Download ChIP-seq data + design

```bash
seqnado download filereport.tsv \
    --outdir fastqs \
    --assay chip \
    --design-output metadata_chip.csv
```

### Dry run to see what would happen

```bash
seqnado download filereport.tsv --dry-run
```

## File Naming

Downloaded files are automatically named based on library layout:

### Paired-End Data
```
{library_name}-{sample_title}_R1.fastq.gz
{library_name}-{sample_title}_R2.fastq.gz
```

### Single-End Data
```
{library_name}-{sample_title}.fastq.gz
```

For example:
- Paired: `GSM12345-WT_rep1_R1.fastq.gz`, `GSM12345-WT_rep1_R2.fastq.gz`
- Single: `GSM12345-WT_rep1.fastq.gz`

## Library Layout Detection

The download command uses the `library_layout` column from your TSV to determine how to process each sample:

- **PAIRED**: Uses `geo_download_paired` rule → creates `_R1.fastq.gz` and `_R2.fastq.gz`
- **SINGLE**: Uses `geo_download_single` rule → creates `.fastq.gz` (no R1/R2 suffix)

This approach ensures proper file structure without creating empty placeholder files.

## Troubleshooting

### Download fails

The download rule includes automatic retry logic via Snakemake's resource scaling — memory and time are doubled on each retry attempt. Check logs for details:
- Full logging in `logs/geo_download/{sample}.log`

### Missing columns in TSV

Make sure your TSV has these required columns:
- `run_accession`
- `sample_title`
- `library_name`
- `library_layout` (with values 'PAIRED' or 'SINGLE')

The `library_layout` column is typically included in TSV downloads from ENA. If it's missing from your source:
1. Check the SRA/GEO database for the layout information
2. Add it manually to your TSV file
3. Or download the full metadata from ENA which includes this column

### Memory issues

If downloads run out of memory, adjust the resources in the Snakemake profile or reduce the number of parallel downloads with `-c`.

## Integration with SeqNado Pipeline

After downloading and generating a design file, you can run the full SeqNado pipeline:

```bash
# 1. Download data
seqnado download filereport.tsv --outdir fastqs --assay rna

# 2. Generate config
seqnado config rna

# 3. Run pipeline
seqnado pipeline rna -c 20
```
