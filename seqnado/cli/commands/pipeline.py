"""Run the data processing pipeline command."""
from __future__ import annotations

import contextlib
import os
import subprocess
from importlib import resources
from pathlib import Path
from typing import List, Optional

import typer
from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.autocomplete import assay_autocomplete
from seqnado.cli.snakemake_builder import SnakemakeCommandBuilder
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    require_snakemake,
    TOP_LEVEL_PASS_THROUGH,
    execute_snakemake,
    resolve_profile,
    print_logo,
    verbose_option,
    preset_option,
)
from seqnado.outputs.multiomics import find_assay_config_paths
from seqnado.utils import (
    resolve_profile_path,
    extract_cores_from_options,
    create_flag_filter,
)


def _run_multiomics_pipeline(
    config_files: List[Path],
    preset: str,
    profile: Optional[Path],
    cores: int,
    queue: Optional[str],
    cleaned_opts: List[str],
    pkg_root_trav,
    verbose: bool,
    print_cmd: bool,
) -> int:
    """
    Run the multiomics pipeline with Snakefile_multi.
    
    Args:
        config_files: List of detected config files (e.g., [config_rna.yaml, config_atac.yaml])
        preset: Snakemake job profile preset
        profile: Custom Snakemake profile path (overrides preset)
        cores: Number of parallel cores
        queue: Slurm queue/partition for the `ss` preset
        cleaned_opts: Cleaned command-line options
        pkg_root_trav: Package root Traversable
        verbose: Increase logging verbosity
        print_cmd: Print the Snakemake command before running it
    
    Returns:
        Exit code from subprocess.run()
    """
    logger.info(f"Multiomic mode detected: found {len(config_files)} config files")
    logger.info(f"Assays: {', '.join([a.value for a in config_files])}")
    
    snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile_multi")
    profile_ctx, profile_path, is_custom = resolve_profile(preset, profile, pkg_root_trav)
    
    try:
        with (
            resources.as_file(snake_trav) as snakefile_path,
            profile_ctx as resolved_profile_path,
        ):
            # Use resolved_profile_path from context if profile_path was None (Traversable case)
            final_profile_path = profile_path or resolved_profile_path
                
            if not Path(snakefile_path).exists():
                logger.error(f"Snakefile_multi not found: {snakefile_path}")
                raise typer.Exit(code=1)

            # Initialize builder
            builder = SnakemakeCommandBuilder(Path(snakefile_path), cores)
            
            # Build workflow_args for nested Snakemake runs
            workflow_args: List[str] = []
            
            # Nested workflow should get profile and default-resources if applicable
            if final_profile_path:
                workflow_args.append(f"--profile {final_profile_path}")
            if queue and preset and preset.startswith("s"):
                workflow_args.append(f"--default-resources slurm_partition={queue}")
            
            # Always enable these nested-run friendly flags
            workflow_args.extend(
                [
                    "--printshellcmds",
                    "--rerun-incomplete",
                    "--show-failed-logs",
                ]
            )
            
            # Add all cleaned options to the nested workflow args
            workflow_args.extend(cleaned_opts)
            
            # Policy: some flags must also be present on the top-level snakemake invocation.
            # TOP_LEVEL_PASS_THROUGH can be adjusted if you want more/fewer flags forwarded.
            should_pass_to_top_level = create_flag_filter(TOP_LEVEL_PASS_THROUGH)
            
            # Filter cleaned_opts into top-level opts (preserving order)
            top_level_opts = [o for o in cleaned_opts if should_pass_to_top_level(o)]
            
            # Add top-level opts to builder
            if top_level_opts:
                builder.add_pass_through_args(top_level_opts)
            
            # Add workflow_args config for nested runs
            builder.add_workflow_args(workflow_args)
            
            # Add profile if present
            if profile_path:
                builder.add_profile_from_path(profile_path)
                if profile:
                    logger.info(f"Using custom Snakemake profile: {profile_path}")
                else:
                    logger.info(
                        f"Using Snakemake profile preset '{preset}' -> {profile_path}"
                    )
            
            # Run in project directory for multiomics to avoid lock conflicts
            builder.add_directory(".")
            
            # Build command
            cmd = builder.build()
            
            # Final working directory setup and execution
            cwd = str(Path(".").resolve())
            exit_code = execute_snakemake(cmd, cwd, print_cmd)
            return exit_code
            
    except typer.Exit:
        raise
    except Exception as e:
        logger.exception("Failed to run multiomics snakemake: %s", e)
        raise typer.Exit(code=1)


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
    """
    Run the data processing pipeline for ASSAY.

    This command wraps Snakemake to execute the SeqNado workflows. Any additional arguments
    provided after the known options are passed directly to Snakemake.

    Args:
        ctx: Typer context for capturing extra arguments
        assay: Assay type (required for single-assay, optional for Multiomic mode)
        config_file: Path to a SeqNado config YAML
        show_version: Print SeqNado version and exit
        preset: Snakemake job profile preset
        profile: Path to a Snakemake profile directory (overrides preset)
        clean_symlinks: Remove symlinks created by previous runs
        scale_resources: Scale memory/time (env: SCALE_RESOURCES)
        verbose: Increase logging verbosity
        queue: Slurm queue/partition for the `ss` preset
        print_cmd: Print the Snakemake command before running it

    Example:
        seqnado pipeline rna -n --unlock
    will run a dry-run of the RNA-seq pipeline and unlock the working directory if needed.
    """
    _configure_logging(verbose)

    if show_version:
        from importlib.metadata import version as _pkg_version

        logger.info(f"SeqNado version {_pkg_version('seqnado')}")
        raise typer.Exit(code=0)

    require_snakemake()

    # Detect multiomics configs early
    config_files = find_assay_config_paths(Path("."))
    use_multiomics = len(config_files) > 1 and not config_file and not assay

    # If the user accidentally put a flag into the assay position (e.g. `-n`),
    # treat it as not providing an assay and restore that token to args.
    if assay and assay.startswith("-"):
        logger.debug(
            "Argument provided in 'assay' position looks like a flag; treating as no assay."
        )
        flag_in_assay = assay
        assay = None
        # recompute multiomics condition (if config_file specified, it's single-assay)
        use_multiomics = len(config_files) > 1 and not config_file
    else:
        flag_in_assay = None

    if not assay and not use_multiomics:
        logger.error(
            "No assay provided. Use `seqnado pipeline ASSAY` or run from a directory with multiple config_*.yaml files for Multiomic mode."
        )
        raise typer.Exit(code=2)

    # Now collect extra args from Typer context (this contains only unknown/extra args).
    raw_extra_args = list(ctx.args)  # tokens Typer didn't map to declared params
    # If we restored a flag from the assay position, prepend it to preserve ordering.
    if flag_in_assay:
        raw_extra_args.insert(0, flag_in_assay)

    # Extract cores and produce cleaned options
    cleaned_opts, cores = extract_cores_from_options(raw_extra_args)

    # Sensible default cores logic for multiomics: at least one core per assay unless user requested more.
    if use_multiomics:
        cores = max(cores, len(config_files))

    os.environ["SCALE_RESOURCES"] = str(scale_resources)

    if clean_symlinks:
        # If assay is None (multiomics), we still have to pick a path; skip if not present
        if assay:
            target = Path(f"seqnado_output/{assay}/fastqs")
            logger.info(f"Cleaning symlinks in {target} ...")
            for link in target.glob("*"):
                if link.is_symlink():
                    link.unlink(missing_ok=True)
        else:
            logger.info("clean_symlinks requested but no assay specified; skipping.")

    pkg_root_trav = _pkg_traversable("seqnado")

    # Dispatch to multiomics if detected
    if use_multiomics:
        exit_code = _run_multiomics_pipeline(
            config_files=config_files,
            preset=preset,
            profile=profile,
            cores=cores,
            queue=queue,
            cleaned_opts=cleaned_opts,
            pkg_root_trav=pkg_root_trav,
            verbose=verbose,
            print_cmd=print_cmd,
        )
        raise typer.Exit(code=exit_code)

    # Single-assay mode
    snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile")
    if not config_file:
        config_file = Path(f"config_{assay}.yaml")
        if not config_file.exists():
            logger.error(
                f"No config file provided and default not found: {config_file}"
            )
            raise typer.Exit(code=1)

    profile_ctx, profile_path, is_custom = resolve_profile(preset, profile, pkg_root_trav)

    try:
        with (
            resources.as_file(snake_trav) as snakefile_path,
            profile_ctx as resolved_profile_path,
        ):
            # Use resolved_profile_path from context if profile_path was None (Traversable case)
            final_profile_path = profile_path or resolved_profile_path
                
            if not Path(snakefile_path).exists():
                logger.error(
                    f"Snakefile for assay '{assay}' not found: {snakefile_path}"
                )
                raise typer.Exit(code=1)

            # Initialize builder for single-assay mode
            builder = SnakemakeCommandBuilder(Path(snakefile_path), cores)
            
            # Add configfile
            builder.add_configfile(config_file)
            
            # Add profile with logging
            builder.add_profile_with_logging(final_profile_path, preset, is_custom)
            
            # Add queue if applicable
            if queue and preset.startswith("s"):
                builder.add_queue(queue, preset)
            
            # Add pass-through args
            builder.add_pass_through_args(cleaned_opts)
            
            # Build command
            cmd = builder.build()
            
            # Optional: print nice ASCII logo if present
            print_logo(pkg_root_trav)

            # Final working directory setup and execution
            cwd = str(Path(".").resolve())
            exit_code = execute_snakemake(cmd, cwd, print_cmd)
            raise typer.Exit(code=exit_code)
    except typer.Exit:
        raise
    except Exception as e:
        logger.exception("Failed to run snakemake: %s", e)
        raise typer.Exit(code=1)
