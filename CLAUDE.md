# SeqNado — CLAUDE.md

## Project Overview

SeqNado is a modular **Snakemake-based bioinformatics toolkit** for analysing high-throughput sequencing data. It provides automated, reproducible workflows from raw FASTQ files to publication-ready results across multiple assay types.

## Supported Assays

| Assay | Code | Description |
|-------|------|-------------|
| ATAC-seq | `atac` | Chromatin accessibility — TSS enrichment, peak calling, TN5 shifting |
| ChIP-seq | `chip` | Protein–DNA interactions — spike-in normalisation, IP/input pairing |
| CUT&Tag | `cat` | Low-input chromatin profiling |
| RNA-seq | `rna` | Gene expression — DESeq2 DE, featureCounts/Salmon quantification |
| SNP Analysis | `snp` | Variant calling (BCFtools / DeepVariant) |
| Methylation | `meth` | TAPS/Bisulfite — Methyldackel integration |
| CRISPR Screens | `crispr` | Guide-level quantification, MAGeCK |
| Micro-Capture-C | `mcc` | 3D chromatin conformation |
| Multiomics | — | Integrated multi-assay processing in a single workflow |

## Repository Structure

```
SeqNado/
├── seqnado/
│   ├── cli/              # Typer-based CLI (app.py, commands/, snakemake_builder.py)
│   ├── config/           # Pydantic config models (core.py, configs.py, multiomics.py, mixins.py)
│   ├── workflow/         # Snakemake workflows
│   │   ├── Snakefile             # Main workflow
│   │   ├── Snakefile_multi       # Multiomics workflow
│   │   ├── Snakefile_genome      # Genome-building workflow
│   │   ├── rules/                # 24+ rule directories by function
│   │   ├── scripts/              # Python/shell rule scripts
│   │   ├── envs/                 # Conda environment YAMLs per tool
│   │   ├── config/               # Config templates
│   │   └── helpers/              # Workflow helper modules
│   ├── inputs/           # FastqCollection classes & validation
│   ├── outputs/          # SeqnadoOutputFactory — generates output file paths
│   ├── data/             # Config/preset/genome data templates
│   ├── tools/            # Tool management & container strategies
│   ├── core.py           # Core enums: Assay, PeakCallingMethod, SpikeInMethod, etc.
│   └── utils.py          # Utilities
├── tests/
│   ├── unit/             # Fast, isolated unit tests
│   └── pipeline/         # End-to-end Snakemake tests
├── docs/                 # MkDocs documentation
├── containers/           # Container definitions (Singularity/Apptainer)
├── pyproject.toml        # Package metadata & dependencies
└── environment.yml       # Conda environment spec
```

## Key Technologies

- **Python 3.12+** — core implementation
- **Snakemake ≥9.12.0** — workflow orchestration
- **Pydantic ≤2.12.5** — config validation and data modelling
- **Typer** — CLI framework
- **Singularity/Apptainer** — container runtime (optional)
- **deepTools, Bowtie2, STAR, MACS2/3, SEACR, HOMER, LanceOtron, DESeq2, edgeR, Salmon, featureCounts, Samtools, Picard, Trim Galore** — integrated bioinformatics tools

## Configuration System

### Three layers of configuration

1. **Project config** (`seqnado config <assay>`) — YAML per project; tool selection, feature flags, output options. Validated via Pydantic models in `seqnado/config/`.
2. **Genome config** (`seqnado genomes`) — reference genomes, alignment indices, annotation files stored in `~/.config/seqnado/`.
3. **Execution profiles** (`~/.config/snakemake/profile_*/`) — 6 presets:
   - `le` local, `ls` local+Singularity, `lc` local+Conda, `ld` local+Docker
   - `ss` SLURM HPC, `t` test

### Pydantic config hierarchy

- `BaseAssayConfig` — shared features (bigwigs, heatmaps, spike-in flags)
- Assay-specific: `ATACAssayConfig`, `ChIPAssayConfig`, `RNAAssayConfig`, …
- Tool configs: `PeakCallingConfig`, `BigwigConfig`, `SpikeInConfig`, `RNAQuantificationConfig`
- Mixins in `mixins.py` for shared validation logic

## FASTQ Naming Convention

Files must follow strict naming patterns (parsed by `seqnado/inputs/`):

```
# ATAC / RNA-seq
sample-name-rep1_R1.fastq.gz
sample-name-rep1_R2.fastq.gz

# ChIP-seq / CUT&Tag (antibody vs input)
sample-name-rep1_Antibody_R1.fastq.gz
sample-name-rep1_Input_R1.fastq.gz
```

Parsed fields: `sample_id`, `replicate`, `condition`, `control_type`.

## CLI Commands

```bash
seqnado init                          # Set up genomes, containers, Snakemake profiles
seqnado config <assay>                # Generate project config YAML
seqnado design <assay> fastqs/*       # Auto-generate sample metadata from FASTQs
seqnado pipeline <assay> --preset le  # Run workflow (--preset ss for SLURM)
seqnado genomes {list,build,update,template}
seqnado download <metadata.tsv>       # Download from GEO/SRA
seqnado tools {list,install,init}
```

## Testing

Framework: **pytest**

```bash
pytest tests/unit/          # Fast unit tests
pytest tests/pipeline/      # End-to-end pipeline tests
pytest -m snakemake         # Tests that invoke Snakemake via subprocess
```

Markers: `unit`, `integration`, `pipeline`, `snakemake`, `slow`, `requires_data`, `requires_apptainer`.

## Installation

```bash
# Recommended (Bioconda)
mamba create -n seqnado -c bioconda seqnado
mamba activate seqnado

# Alternative
pip install seqnado        # or: uv pip install seqnado

# Post-install setup
seqnado init
```

## Workflow Patterns & Conventions

- **Output paths**: `{output_dir}/{assay}/{sample_id}/{step}/`
- **Rule organisation**: rules are split by function — `alignment/`, `peaks/`, `quantification/`, `normalisation/`, `visualisation/`, `multiomics/`, etc.
- **Resource scaling**: set `SCALE_RESOURCES` env variable
- **Sample grouping keys**: `normalisation_groups`, `consensus_groups`, `condition_groups`
- Input files are **symlinked** into the project directory to avoid duplication
- Each Snakemake rule declares its own conda env (`envs/`) or container image

## Best Practices

### Adding a new Snakemake rule

Every rule **must** include `log:`, `benchmark:`, and `message:` directives — this is enforced by the unit test `test_snakemake_rule_directives.py`. A minimal rule looks like:

```python
rule my_rule:
    input:  OUTPUT_DIR + "/aligned/{sample}.bam"
    output: OUTPUT_DIR + "/processed/{sample}.bed"
    params: options=CONFIG.third_party_tools.mytool.command_line_arguments
    threads: CONFIG.third_party_tools.mytool.threads
    resources:
        runtime=lambda wildcards, attempt: define_time_requested(initial_value=2, attempts=attempt, scale=SCALE_RESOURCES),
        mem=lambda wildcards, attempt: define_memory_requested(initial_value=4, attempts=attempt, scale=SCALE_RESOURCES),
    container: "oras://ghcr.io/alsmith151/seqnado_pipeline:latest"
    log:    OUTPUT_DIR + "/logs/my_rule/{sample}.log"
    benchmark: OUTPUT_DIR + "/.benchmark/my_rule/{sample}.tsv"
    message: "Running my_rule for {wildcards.sample}"
    shell: "mytool {params.options} {input} > {output} 2> {log}"
```

- Place the rule file in the appropriate `rules/<category>/` directory
- Use `temp()` on intermediate outputs to save disk space
- Use `ruleorder:` when paired/single-end variants of a rule exist
- Import helpers from `seqnado/workflow/helpers/` rather than reimplementing them
- Use `CONFIG` and `OUTPUT_DIR` globals that are set in the Snakefile — do not hardcode paths

### Adding or modifying configuration

Config models live in `seqnado/config/`. Follow these patterns:

- Subclass `BaseAssayConfig` (or the relevant mixin) — do not create standalone Pydantic models for assay configs
- Use `BeforeValidator` with `none_str_to_none` to handle YAML `"none"` strings gracefully
- Use `field_validator(mode="before")` for normalisation/coercion; `model_validator(mode="after")` for cross-field logic
- Path validation should use `PathValidatorMixin.validate_path_exists()` with `skip_path_validation` context support, so tests can bypass filesystem checks
- Raise `UserFriendlyError` (not raw `ValueError`) for errors the user will see directly in the CLI
- Every new config field needs a corresponding test in `tests/unit/config/`

### Writing tests

- **Unit tests** go in `tests/unit/`, **CLI tests** in `tests/cli/`, **end-to-end** in `tests/pipeline/`
- Mark tests with the appropriate pytest marker: `@pytest.mark.unit`, `@pytest.mark.pipeline`, `@pytest.mark.snakemake`, etc.
- Use `tmp_path` (pytest fixture) for any tests that need to create temporary files — never hardcode `/tmp` paths
- Config tests should use `skip_path_validation=True` context where real index/genome files are unavailable:

  ```python
  model = MyConfig.model_validate(data, context={"skip_path_validation": True})
  ```

- Parametrize across assay types using `@pytest.mark.parametrize` rather than copy-pasting test bodies
- Pipeline tests invoke Snakemake via subprocess using helpers in `tests/pipeline/helpers/` — do not construct Snakemake commands by hand

### CLI commands

- CLI commands are Typer functions in `seqnado/cli/commands/`; register them in `cli/app.py`
- Use `rich` for console output — do not use bare `print()`
- Validate user-provided paths and sample sheets early (in the CLI layer via `data_validation.py`), not deep inside config models, so errors surface with helpful messages
- The `snakemake_builder.py` module constructs Snakemake invocations — add new execution flags there rather than in individual commands

### General Python conventions

- Type-annotate all function signatures; use `pathlib.Path` for file paths, not raw strings
- Enums for fixed vocabularies (assay types, methods, etc.) are defined in `seqnado/core.py` — add new values there and nowhere else
- Avoid adding new top-level dependencies without updating both `pyproject.toml` and `environment.yml`
- Keep workflow scripts in `seqnado/workflow/scripts/` thin — move reusable logic into `seqnado/` Python modules that can be unit-tested independently of Snakemake

## Core Enums (`seqnado/core.py`)

| Enum | Values |
|------|--------|
| `Assay` | atac, chip, cat, rna, snp, meth, crispr, mcc |
| `PeakCallingMethod` | MACS2, MACS3, SEACR, HOMER, LanceOtron |
| `SpikeInMethod` | Orlando, WITH_INPUT, DESeq2, edgeR |
| `DataScalingTechnique` | Unscaled, CSAW, CPM, RPKM, SPIKEIN |
| `PileupMethod` | deeptools, HOMER, BAMNADO, METHYLDACKEL |


## Development and conventions

- Seqnado only requires Python packages to run. Conda is only needed if users wants to run the Snakemake workflows locally without containers. We discourage this because it can lead to dependency conflicts and reproducibility issues, but it is supported for users who prefer it.

- For a quick development setup, recommend just creating a uv virtual environment and installing the package with `uv pip install -e .`.
- Run tests with `uv run pytest` to ensure the correct environment is used.
- To run pipeline tests (anything in `tests/pipeline/`), use `--run-pipeline` when invoking pytest.
- If having issues with `uv`, use `conda activate base` first to ensure the base conda environment is active.