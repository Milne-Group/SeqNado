# SeqNado contributor guidance

This file is the shared repository guidance for coding agents. `CLAUDE.md` imports
it; keep common instructions here so the two entry points stay consistent.

## Project and source of truth

SeqNado is a Python package and Typer CLI that orchestrates Snakemake workflows
for ATAC-seq (`atac`), ChIP-seq (`chip`), CUT&Tag (`cat`), RNA-seq (`rna`), variants
(`snp`), methylation (`meth`), CRISPR screens (`crispr`), and Micro-Capture-C (`mcc`).
It also supports combined multiomics workflows.

- Use `pyproject.toml` as the source of truth for Python support, package
  dependencies, extras, entry points, and pytest settings. Python currently
  requires **3.11 or newer**. Avoid copying dependency pins into this guide.
- Check implementation and tests when README, documentation, or generated
  environment files disagree. Report relevant discrepancies without expanding
  an unrelated change into a repository-wide cleanup.
- Preserve existing user changes and keep edits focused on the requested task.

## Development setup

A Python virtual environment is enough for package development. Conda is optional
for executing workflow tools locally; container execution is preferred for
reproducibility.

```bash
uv venv
uv pip install -e . pytest
uv run seqnado --help
```

Use `uv run` for commands and tests so they use the project environment. If uv has
environment issues on a Conda-managed machine, try `conda activate base` first.
Do not run `seqnado init` merely to develop the package: it configures user-level
genomes, containers, and execution profiles.

## Repository map

- `seqnado/cli/`: Typer app, commands, input validation, and invocation helpers.
  Commands use the shared app from `app_instance.py`; `app.py` imports modules
  to trigger registration and attaches the genomes sub-application.
- `seqnado/cli/snakemake_builder.py`: construction of Snakemake invocations.
- `seqnado/config/`: Pydantic assay/tool models, validation mixins, and user input.
- `seqnado/core.py`: shared enums and fixed vocabularies.
- `seqnado/inputs/`: FASTQ/BAM/bigWig collections, sample grouping, validation.
- `seqnado/outputs/`: output factories and file definitions.
- `seqnado/analysis/`: project analysis helpers.
- `seqnado/workflow/Snakefile`: main assay workflow; `Snakefile_multi` and
  `Snakefile_genome` handle multiomics and genome building.
- `seqnado/workflow/rules/`, `scripts/`, `helpers/`, `envs/`: workflow rules,
  scripts, reusable helpers, and tool environments/execution profiles.
- `seqnado/data/`: configuration templates and packaged reference data.
- `seqnado/tools/` and `containers/`: tool metadata and container management.
- `tests/unit/`, `tests/cli/`, `tests/integration/`, `tests/pipeline/`: test suites.
- `docs/`, `mkdocs.yml`: user and API documentation.

## Implementation conventions

- Follow adjacent code, annotate function signatures, and prefer `pathlib.Path`
  for filesystem operations. Reuse shared enums from `seqnado/core.py`.
- Keep workflow scripts thin; place reusable logic in importable Python modules.
- CLI commands belong in `seqnado/cli/commands/`. Ensure new command modules are
  imported in `app.py`. Use Rich for console output and existing friendly-error
  handling for failures shown to users.
- Validate CLI paths and sample sheets early using the existing validation
  helpers. Extend `snakemake_builder.py` for execution flags.
- Extend the existing assay config hierarchy (`BaseAssayConfig` and relevant
  mixins). Follow existing handling of YAML `"none"` values, before validators
  for coercion, and after validators for cross-field checks.
- Preserve `skip_path_validation` context support when changing path validation.
  Cover new configuration fields and changed validation behavior in
  `tests/unit/config/`; update `seqnado/data/config_template.jinja` when applicable.
- Treat sample naming, IP/control pairing, group membership, and output paths
  as compatibility-sensitive. Check the collection parsers and output factories
  rather than assuming one universal filename or directory pattern.
- Update user documentation when CLI or workflow behavior changes.

## Snakemake rules

- Put rules under the appropriate functional directory in `workflow/rules/`.
- Every rule must declare `log:`, `benchmark:`, and `message:`. This is enforced
  by `tests/unit/test_snakemake_rule_directives.py`.
- Reuse `CONFIG`, `OUTPUT_DIR`, existing output factories, and workflow helpers.
  The main Snakefile defaults to `seqnado_output/{assay}`; output layouts vary
  by processing step.
- Follow neighboring rules for container/Conda declarations, configured tool
  arguments, threads, and retry-aware resource requests. Respect `SCALE_RESOURCES`.
- Use `temp()` only for disposable intermediate outputs; preserve requested
  deliverables. Follow existing paired/single-end selection and rule ordering.

## Tests and verification

Run the smallest relevant suite first and broaden for shared behavior changes:

```bash
uv run pytest tests/unit/ -q
uv run pytest tests/cli/ -q
uv run pytest tests/integration/ -q
uv run pytest tests/unit/test_snakemake_rule_directives.py -q
```

Use `tmp_path` for test files, parametrize assay variants, and apply appropriate
registered markers (see `pyproject.toml` and `tests/conftest.py`). Config models
can bypass reference-file checks in tests via:

```python
model = MyConfig.model_validate(data, context={"skip_path_validation": True})
```

Pipeline tests are opt-in and can download data, build images, and consume
substantial resources. Use the helpers in `tests/pipeline/helpers/` and explicitly
select the affected assay:

```bash
uv run pytest tests/pipeline/test_pipelines.py::test_pipeline \
  --run-pipeline --assays chip --cores 2 -vv -s
```

On macOS without Apptainer/Singularity, add `--preset ld` to use the test-only
Docker fallback. It requires a running Docker daemon and builds the devcontainer
image for `linux/amd64`, then executes Apptainer inside a privileged Linux
container. See `tests/conftest.py` for runtime detection and skip behavior.
The actual `profile_local_docker` currently configures Apptainer: do not describe
`seqnado pipeline --preset ld` as native Docker execution.

Report what was verified and any tests that could not run. Documentation-only
changes normally need content/path checks and `git diff --check`, not pipelines.

## Dependencies and formatting

- Change package requirements in `pyproject.toml`. Environment YAMLs identify
  themselves as generated; regenerate them with
  `python .github/scripts/generate_environments.py` when changing dependencies,
  and review the resulting diff. Keep `uv.lock` consistent with dependency edits.
- `.pre-commit-config.yaml` defines Black, snakefmt, and basic whitespace,
  YAML/JSON, and conflict checks. Use the existing configuration; avoid unrelated
  formatting or dependency upgrades.
- Do not hand-edit generated `seqnado/_version.py`.
