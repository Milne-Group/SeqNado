#!/usr/bin/env python3
from __future__ import annotations

import contextlib
import json
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from datetime import date
from importlib import resources
from pathlib import Path
from typing import TYPE_CHECKING, Any, List, Optional, Tuple

import click
import typer
from loguru import logger

from seqnado.outputs.multiomics import find_assay_config_paths
from seqnado.utils import (
    create_flag_filter,
    get_profile_name,
    get_preset_profiles,
    resolve_profile_path,
)

# Import from refactored modules
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    verbose_option,
    preset_option,
    dry_run_option,
    cores_option,
    _RICH_CONSOLE,
    Text,
    _resolve_working_dir,
)
from seqnado.cli.autocomplete import (
    _assay_names,
    _find_fastqs,
    _profile_autocomplete,
    assay_autocomplete,
    fastq_autocomplete,
)
from seqnado.cli.tools_helpers import (
    _display_tools_list,
    _OptionalCategoryCommand,
    _resolve_category,
    _validate_subcommand,
)
from seqnado.cli.multiomics_helpers import (
    get_multiomics_config_paths,
    should_use_multiomics_mode,
    validate_multiomics_setup,
)

# Optional: prettier tracebacks/console with rich if available
try:
    from rich.traceback import install as _rich_tb_install

    _rich_tb_install(show_locals=False)
except Exception:
    pass

from seqnado import Assay
from seqnado._version import __version__
from seqnado.cli.commands.init import init as init_command
from seqnado.cli.commands.genomes import app as genomes_app
from seqnado.cli.commands.design import design as design_command
from seqnado.cli.commands.config import config as config_command
from seqnado.cli.commands.download import download as download_command
from seqnado.cli.commands.tools import tools as tools_command
from seqnado.cli.commands.pipeline import pipeline as pipeline_command


def version_callback(value: bool):
    """Print version and exit."""
    if value:
        typer.echo(f"SeqNado version {__version__}")
        raise typer.Exit()


app = typer.Typer(
    add_completion=True,
    no_args_is_help=True,
    help="""
[bold]SeqNado CLI[/bold]

Initialize your environment, build configs, create design files, and run pipelines.
Use --help on any subcommand for details.
""",
    callback=lambda version: version_callback(version) if version else None,
)


@app.callback()
def main(
    version: bool = typer.Option(
        False,
        "--version",
        "-v",
        help="Show version and exit.",
        callback=version_callback,
        is_eager=True,
    ),
):
    """SeqNado CLI main entry point."""
    pass


# --------------------------------- init ------------------------------------- #


@app.command(
    help="""
Initialize SeqNado user environment.

- Logs the current Conda environment if active (optional).
- Runs packaged Apptainer/Singularity init (if `apptainer` on PATH).
- Ensures ~/.config/seqnado/genome_config.json exists (template or preset).
"""
)
def init(
    preset: bool = typer.Option(
        False, help="Use packaged preset genomes instead of the editable template."
    ),
    dry_run: bool = dry_run_option(),
    verbose: bool = verbose_option(),
) -> None:
    """Wrapper for init command."""
    init_command(preset=preset, dry_run=dry_run, verbose=verbose)


# -------------------------------- utils ----------------------------------- #

# Add genomes as a sub-application
app.add_typer(genomes_app, name="genomes")


# -------------------------------- config ------------------------------------ #


@app.command(
    help="Build a workflow configuration YAML for the selected ASSAY. If no assay is provided, multiomics mode is used."
)
def config(
    assay: Optional[str] = typer.Argument(
        None,
        metavar="[ASSAY]",
        autocompletion=assay_autocomplete,
        show_choices=True,
        help=", ".join(_assay_names()) + ". If omitted, multiomics mode is used.",
    ),
    make_dirs: bool = typer.Option(
        True,
        "--make-dirs/--no-make-dirs",
        help="Create/don't create the output project directory or fastq subdir.",
    ),
    render_options: bool = typer.Option(
        False, help="Render all options (even if not used by the workflow)."
    ),
    output: Optional[Path] = typer.Option(
        None, "-o", "--output", help="Explicit path for the rendered config file."
    ),
    verbose: bool = verbose_option(),
    interactive: bool = typer.Option(
        True,
        "--interactive/--no-interactive",
        help="Interactively prompt for config values. Non-interactive mode only works for single assay configs (except MCC and multiomics).",
    ),
) -> None:
    """Wrapper for config command."""
    config_command(
        assay=assay,
        make_dirs=make_dirs,
        render_options=render_options,
        output=output,
        verbose=verbose,
        interactive=interactive,
    )


# -------------------------------- tools ------------------------------------ #


@app.command(
    cls=_OptionalCategoryCommand,
    help="List and explore bioinformatics tools available in the SeqNado pipeline.",
)
def tools(
    tool_name: Optional[str] = typer.Argument(
        None,
        metavar="[TOOL]",
        help="Specific tool name to get help for (e.g., fastqc, deeptools)",
    ),
    list_tools: bool = typer.Option(
        False, "--list", "-l", help="List all available tools with descriptions."
    ),
    category: Optional[str] = typer.Option(
        None,
        "--category",
        "-c",
        help=(
            "Filter tools by category name or number. Use without a value to "
            "interactively select a category."
        ),
    ),
    show_options: bool = typer.Option(
        False,
        "--options",
        help="Show tool options/help from the container (requires tool name and apptainer).",
    ),
    show_citation: bool = typer.Option(
        False,
        "--citation",
        help="Show the BibTeX citation for a tool (requires tool name).",
    ),
    subcommand: Optional[str] = typer.Option(
        None,
        "--subcommand",
        "-s",
        help="Specify a tool subcommand for help (e.g. plotHeatmap).",
    ),
    verbose: bool = verbose_option(),
) -> None:
    """Wrapper for tools command."""
    tools_command(
        tool_name=tool_name,
        list_tools=list_tools,
        category=category,
        show_options=show_options,
        show_citation=show_citation,
        subcommand=subcommand,
        verbose=verbose,
    )


# ----------------------------------- download ------------------------------------ #


@app.command(
    help="Download FASTQ files from GEO/SRA using a metadata TSV file and optionally generate a design file."
)
def download(
    metadata_tsv: Path = typer.Argument(
        ...,
        exists=True,
        help="TSV file from GEO/ENA with columns: run_accession, sample_title, library_name, and library_layout (PAIRED or SINGLE).",
    ),
    outdir: Path = typer.Option(
        Path("fastqs"),
        "-o",
        "--outdir",
        help="Output directory for downloaded FASTQ files.",
    ),
    assay: Optional[str] = typer.Option(
        None,
        "-a",
        "--assay",
        autocompletion=assay_autocomplete,
        help="Assay type for generating design file after download. If not provided, only downloads FASTQs.",
    ),
    design_output: Optional[Path] = typer.Option(
        None,
        "-d",
        "--design-output",
        help="Output path for design CSV (default: metadata_{assay}.csv in outdir).",
    ),
    cores: int = cores_option(),
    preset: str = preset_option(),
    profile: Optional[Path] = typer.Option(
        None,
        "--profile",
        "--profiles",
        help="Path to a Snakemake profile directory (overrides --preset).",
    ),
    dry_run: bool = dry_run_option(),
    verbose: bool = verbose_option(),
) -> None:
    """Wrapper for download command."""
    download_command(
        metadata_tsv=metadata_tsv,
        outdir=outdir,
        assay=assay,
        design_output=design_output,
        cores=cores,
        preset=preset,
        profile=profile,
        dry_run=dry_run,
        verbose=verbose,
    )


# -------------------------------- design ------------------------------------ #


@app.command(
    help="Generate a SeqNado design CSV from FASTQ files for ASSAY. If no assay is provided, multiomics mode is used."
)
def design(
    assay: Optional[str] = typer.Argument(
        None,
        metavar="[ASSAY]",
        autocompletion=assay_autocomplete,
        show_choices=True,
        help="Assay type. Options: "
        + ", ".join(_assay_names())
        + ". If omitted, multiomics mode is used.",
    ),
    files: List[Path] = typer.Argument(
        None, metavar="[FASTQ ...]", autocompletion=fastq_autocomplete
    ),
    output: Optional[Path] = typer.Option(
        None,
        "-o",
        "--output",
        help="Output CSV filename (default: metadata_{assay}.csv).",
    ),
    ip_to_control: Optional[str] = typer.Option(
        None,
        "--ip-to-control",
        help="""List of antibody,control pairings for IP assays (e.g. ChIP). Format: 'antibody1:control1,antibody2:control2'
        If provided will assign a control with a specified name to that ip in the metadata. If not provided, controls will be assigned based on a best-effort matching of sample names.
        """,
    ),
    group_by: bool = typer.Option(
        False, "--group-by", help="Group samples by a regular expression or a column."
    ),
    auto_discover: bool = typer.Option(
        True,
        "--auto-discover/--no-auto-discover",
        help="Search common folders if none provided.",
    ),
    interactive: bool = typer.Option(
        True,
        "--interactive/--no-interactive",
        help="Interactively offer to add missing columns using schema defaults.",
    ),
    accept_all_defaults: bool = typer.Option(
        False,
        "--accept-all-defaults",
        help="Non-interactive: auto-add only columns that have a schema default.",
    ),
    deseq2_pattern: Optional[str] = typer.Option(
        None,
        "--deseq2-pattern",
        help="Regex pattern to extract DESeq2 groups from sample names. "
        "First capture group will be used. Example: r'-(\\w+)-rep' for 'sample-GROUP-rep1'",
    ),
    verbose: bool = verbose_option(),
) -> None:
    """Wrapper for design command."""
    design_command(
        assay=assay,
        files=files,
        output=output,
        ip_to_control=ip_to_control,
        group_by=group_by,
        auto_discover=auto_discover,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
        deseq2_pattern=deseq2_pattern,
        verbose=verbose,
    )


# -------------------------------- pipeline ---------------------------------- #


@app.command(
    help="Run the data processing pipeline for ASSAY (Snakemake under the hood). Any additional arguments are passed to Snakemake (e.g., `seqnado pipeline rna -n` for dry-run, `--unlock`, etc.).",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def pipeline(
    ctx: typer.Context,
    assay: Optional[str] = typer.Argument(
        None,
        metavar="[ASSAY]",
        autocompletion=assay_autocomplete,
        help="Assay type (required for single-assay, optional for Multiomic mode)",
    ),
    config_file: Optional[Path] = typer.Option(
        None,
        "--configfile",
        help="Path to a SeqNado config YAML (default: config_<ASSAY>.yaml).",
    ),
    show_version: bool = typer.Option(
        False, "--version", help="Print SeqNado version and exit."
    ),
    preset: str = preset_option(),
    profile: Optional[Path] = typer.Option(
        None,
        "--profile",
        "--profiles",
        help="Path to a Snakemake profile directory (overrides --preset).",
    ),
    clean_symlinks: bool = typer.Option(
        False, help="Remove symlinks created by previous runs."
    ),
    scale_resources: float = typer.Option(
        1.0, "-s", "--scale-resources", help="Scale memory/time (env: SCALE_RESOURCES)."
    ),
    verbose: bool = verbose_option(),
    queue: Optional[str] = typer.Option(
        "short", "-q", "--queue", help="Slurm queue/partition for the `ss` preset."
    ),
    print_cmd: bool = typer.Option(
        False, "--print-cmd", help="Print the Snakemake command before running it."
    ),
) -> None:
    """Wrapper for pipeline command."""
    pipeline_command(
        ctx=ctx,
        assay=assay,
        config_file=config_file,
        show_version=show_version,
        preset=preset,
        profile=profile,
        clean_symlinks=clean_symlinks,
        scale_resources=scale_resources,
        verbose=verbose,
        queue=queue,
        print_cmd=print_cmd,
    )


# -------------------------------- Entrypoint --------------------------------
if __name__ == "__main__":
    _resolve_working_dir()
    app()
