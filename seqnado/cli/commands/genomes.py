"""Manage genome configurations command."""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
from importlib import resources
from pathlib import Path
from typing import List, Optional

import typer
from loguru import logger

from seqnado.cli.snakemake_builder import SnakemakeCommandBuilder
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    require_snakemake,
    TOP_LEVEL_PASS_THROUGH,
    resolve_profile,
    verbose_option,
    preset_option,
    dry_run_option,
    cores_option,
    execute_snakemake,
)
from seqnado.cli.autocomplete import _assay_names, assay_autocomplete
from seqnado.utils import create_flag_filter


# Create genomes app
app = typer.Typer(
    help="Manage genome configurations",
    no_args_is_help=True,
)


@app.command(name="list")
def list_genomes(
    assay: str = typer.Argument(
        "atac",
        metavar="ASSAY",
        autocompletion=assay_autocomplete,
        help="Assay type. Options: " + ", ".join(_assay_names()),
    ),
    verbose: bool = verbose_option(),
) -> None:
    """Show packaged and user genome presets."""
    _configure_logging(verbose)
    
    # Import locally for snappy startup
    from seqnado.config import load_genome_configs  # local import

    try:
        cfg = load_genome_configs(assay=assay)
    except Exception as e:
        logger.error(f"Failed to load genome config: {e}")
        raise typer.Exit(code=1)

    if not cfg:
        logger.warning("No genome config found.")
        raise typer.Exit(code=0)

    for genome_name, details in cfg.items():
        typer.echo(f"[bold]{genome_name}[/bold]")
        try:
            items = details.dict()
        except Exception:
            items = dict(details)
        for k, v in items.items():
            typer.echo(f"  {k}: {v or '[not set]'}")
        typer.echo("")
    raise typer.Exit(code=0)


@app.command(name="edit")
def edit_genomes(
    verbose: bool = verbose_option(),
) -> None:
    """Open user genome config in $EDITOR."""
    _configure_logging(verbose)
    
    env = os.environ.get("SEQNADO_GENOME_CONFIG")
    cfg_path = (
        Path(env)
        if env
        else Path.home().joinpath(".config", "seqnado", "genome_config.json")
    )
    if not cfg_path.exists():
        logger.error(
            f"Genome config not found: {cfg_path} (try `seqnado init` first)"
        )
        raise typer.Exit(code=1)

    editor = (
        os.environ.get("EDITOR")
        or os.environ.get("VISUAL")
        or ("notepad" if sys.platform == "win32" else "nano")
    )
    editor_cmd = shlex.split(editor)
    if shutil.which(editor_cmd[0]) is None:
        logger.error(
            f"Editor '{editor_cmd[0]}' not found. Please set $EDITOR to your preferred editor."
        )
        raise typer.Exit(code=4)

    try:
        subprocess.check_call(editor_cmd + [str(cfg_path)])
    except subprocess.CalledProcessError as e:
        logger.error(f"Editor returned non-zero exit status: {e}")
        raise typer.Exit(code=3)
    except FileNotFoundError:
        logger.error(
            f"Editor '{editor_cmd[0]}' not found. Please set $EDITOR to your preferred editor."
        )
        raise typer.Exit(code=4)
    except Exception as e:
        logger.error(f"Failed to launch editor: {e}")
        raise typer.Exit(code=5)

    typer.echo(f"Edited genome config: [italic]{cfg_path}[/italic]")
    raise typer.Exit(code=0)


@app.command(
    name="build",
    context_settings={"allow_extra_args": True, "ignore_unknown_options": True},
)
def build_genomes(
    ctx: typer.Context,
    name: str = typer.Option(
        ...,
        "--name",
        "-g",
        help="Genome name(s), comma-separated for multiple (e.g., hg38 or hg38,mm39,dm6)",
    ),
    outdir: Path = typer.Option(
        Path.cwd() / "genome_build", "--outdir", "-o", help="Output directory for build"
    ),
    spikein: Optional[str] = typer.Option(
        None,
        "--spikein",
        "-sp",
        help="Spike-in genome name for composite builds (e.g., mm39)",
    ),
    preset: str = preset_option(),
    profile: Optional[Path] = typer.Option(
        None,
        "--profile",
        "--profiles",
        help="Path to a Snakemake profile directory (overrides --preset).",
    ),
    cores: int = cores_option(),
    scale_resources: float = typer.Option(
        1.0, "--scale-resources", help="Scale memory/time (for build subcommand)."
    ),
    dry_run: bool = dry_run_option(),
    verbose: bool = verbose_option(),
) -> None:
    """Download genome and build indices via Snakemake."""
    _configure_logging(verbose)
    
    genomes = [g.strip() for g in name.split(",") if g.strip()]
    if not genomes:
        logger.error("No valid genome names provided.")
        raise typer.Exit(code=2)

    if spikein and len(genomes) > 1:
        logger.error(
            "Spike-in (--spikein) is only supported with a single genome, not multiple."
        )
        raise typer.Exit(code=2)

    require_snakemake()

    outdir = Path(outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    pkg_root_trav = _pkg_traversable("seqnado")
    snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile_genome")

    profile_ctx, profile_path, is_custom = resolve_profile(
        preset, profile, pkg_root_trav, warn_unknown=True
    )

    try:
        with (
            resources.as_file(snake_trav) as snakefile_path,
            profile_ctx as resolved_profile_path,
        ):
            # Use resolved_profile_path from context if profile_path was None (Traversable case)
            final_profile_path = profile_path or resolved_profile_path
            
            if not Path(snakefile_path).exists():
                logger.error(f"Snakefile_genome not found: {snakefile_path}")
                raise typer.Exit(code=1)

            # Initialize builder
            builder = SnakemakeCommandBuilder(Path(snakefile_path), cores)
            
            # Prepare config values
            config_values = {
                "genome": ",".join(genomes),
                "output_dir": str(outdir),
            }
            if spikein:
                config_values["spikein"] = spikein
            
            # Add all config values in one call
            builder.add_config(**config_values)
            
            # Add profile with logging
            builder.add_profile_with_logging(final_profile_path, preset, is_custom)
            
            # Add dry-run flag if requested
            if dry_run:
                builder.add_dry_run()
            
            # Create flag filter for allowed pass-through flags
            should_pass_to_snakemake = create_flag_filter(TOP_LEVEL_PASS_THROUGH)
            
            # Pass through validated extra args from ctx
            if ctx and ctx.args:
                filtered_args = [
                    arg for arg in ctx.args if should_pass_to_snakemake(arg)
                ]
                builder.add_pass_through_args(filtered_args)

            # Build command
            cmd = builder.build()
            
            genome_label = (
                f"{genomes[0]}_{spikein}" if spikein else ",".join(genomes)
            )
            logger.info(f"Building genome(s): {genome_label}")
            logger.info(f"Output directory: {outdir}")
            
            # Execute snakemake
            cwd = str(Path(".").resolve())
            exit_code = execute_snakemake(cmd, cwd, verbose)
            raise typer.Exit(code=exit_code)

    except typer.Exit:
        raise
    except Exception as e:
        logger.exception("Failed to run genome build: %s", e)
        raise typer.Exit(code=1)


@app.command(name="fastqscreen")
def fastqscreen_config(
    output: Optional[Path] = typer.Option(
        None,
        "-s",
        "--screen",
        help="Output path for fastqscreen config (default: ~/.config/seqnado/fastq_screen.conf)",
    ),
    threads: int = typer.Option(
        8,
        "-t",
        "--threads",
        help="Number of threads for Bowtie2",
    ),
    no_contaminants: bool = typer.Option(
        False,
        "--no-contaminants",
        help="Exclude contaminant databases",
    ),
    contaminant_path: Optional[Path] = typer.Option(
        None,
        "--contaminant-path",
        help="Path to contaminant reference files",
    ),
    verbose: bool = verbose_option(),
) -> None:
    """Generate FastqScreen configuration file."""
    from seqnado.config.fastq_screen_generator import (
        generate_fastq_screen_config,
        load_genome_configs_for_fastqscreen,
    )

    _configure_logging(verbose)

    try:
        # Set default output path
        if output is None:
            output = Path.home() / ".config" / "seqnado" / "fastq_screen.conf"

        # Prompt for contaminant path if not provided and not explicitly disabled
        if not no_contaminants and contaminant_path is None:
            contaminant_input = typer.prompt(
                "Path to contaminant references (leave empty to skip contaminants)",
                default="",
                show_default=False,
            )
            if contaminant_input.strip():
                contaminant_path = Path(contaminant_input.strip())
            else:
                logger.info("Skipping contaminant databases")

        # Load genome configs
        logger.info("Loading genome configurations...")
        genome_configs = load_genome_configs_for_fastqscreen()
        logger.info(f"Found {len(genome_configs)} genome configurations")

        # Generate config
        logger.info(f"Generating FastqScreen config: {output}")

        generate_fastq_screen_config(
            genome_configs=genome_configs,
            output_path=output,
            threads=threads,
            include_contaminants=not no_contaminants
            and contaminant_path is not None,
            contaminant_base_path=contaminant_path,
        )

        logger.success(f"FastqScreen config generated successfully: {output}")

    except FileNotFoundError as e:
        logger.error(str(e))
        raise typer.Exit(code=1)
    except Exception as e:
        logger.error(f"Failed to generate FastqScreen config: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise typer.Exit(code=1)
