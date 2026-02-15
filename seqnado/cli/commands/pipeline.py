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

from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    require_snakemake,
    TOP_LEVEL_PASS_THROUGH,
)
from seqnado.outputs.multiomics import find_assay_config_paths
from seqnado.utils import (
    resolve_profile_path,
    extract_cores_from_options,
    create_flag_filter,
)


def pipeline(
    ctx: typer.Context,
    assay: Optional[str] = None,
    config_file: Optional[Path] = None,
    show_version: bool = False,
    preset: str = "le",
    profile: Optional[Path] = None,
    clean_symlinks: bool = False,
    scale_resources: float = 1.0,
    verbose: bool = False,
    queue: Optional[str] = "short",
    print_cmd: bool = False,
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
    # Note: we don't materialize extra args until we've handled `assay`.
    if assay and assay.startswith("-"):
        logger.debug(
            "Argument provided in 'assay' position looks like a flag; treating as no assay."
        )
        # Put that token back into the argv list for proper passthrough.
        # Typer removed it from ctx.args when it assigned to `assay`.
        # We'll restore it below when building raw_args.
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
        # Insert at 0 to approximate original position; append could bury it.
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

    # Choose Snakefile
    if use_multiomics:
        logger.info(f"Multiomic mode detected: found {len(config_files)} config files")
        logger.info(f"Assays: {', '.join([c.name for c in config_files])}")
        snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile_multi")
        config_file = None
    else:
        snake_trav = pkg_root_trav.joinpath("workflow").joinpath("Snakefile")
        if not config_file:
            config_file = Path(f"config_{assay}.yaml")
            if not config_file.exists():
                logger.error(
                    f"No config file provided and default not found: {config_file}"
                )
                raise typer.Exit(code=1)

    profile_ctx: contextlib.AbstractContextManager
    if profile:
        profile_path = Path(profile).expanduser()
        if not profile_path.exists():
            logger.error(f"Snakemake profile not found: {profile_path}")
            raise typer.Exit(code=1)
        profile_ctx = contextlib.nullcontext(profile_path)
    else:
        # Resolve profile using refactored helper
        profile_trav = resolve_profile_path(preset, pkg_root_trav)
        profile_ctx = (
            resources.as_file(profile_trav) if profile_trav else contextlib.nullcontext()
        )

    try:
        with (
            resources.as_file(snake_trav) as snakefile_path,
            profile_ctx as profile_path,
        ):
            if not Path(snakefile_path).exists():
                logger.error(
                    f"Snakefile for assay '{assay}' not found: {snakefile_path}"
                )
                raise typer.Exit(code=1)

            # Base command
            cmd: List[str] = [
                "snakemake",
                "--snakefile",
                str(snakefile_path),
                "--show-failed-logs",
            ]

            # Ensure cores is passed
            cmd += ["-c", str(cores)]

            # Build workflow_args for nested Snakemake runs (only used in Multiomic mode)
            workflow_args: List[str] = []
            if use_multiomics:
                # Nested workflow should get profile and default-resources if applicable
                if profile_path:
                    workflow_args.append(f"--profile {profile_path}")
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

                # Add all cleaned options to the nested workflow args so nested snakemake sees them.
                # They will be joined into a single string and injected via --config workflow_args="..."
                workflow_args.extend(cleaned_opts)

                # Policy: some flags must also be present on the top-level snakemake invocation.
                # Example: -n/--dry-run, --unlock, printshellcmds, --rerun-incomplete are meaningful at top-level.
                # TOP_LEVEL_PASS_THROUGH can be adjusted if you want more/fewer flags forwarded.
                should_pass_to_top_level = create_flag_filter(TOP_LEVEL_PASS_THROUGH)

                # Filter cleaned_opts into top-level opts (preserving order)
                top_level_opts = [
                    o for o in cleaned_opts if should_pass_to_top_level(o)
                ]

                # Append top-level opts to the top-level command so those flags take effect immediately.
                if top_level_opts:
                    cmd += top_level_opts

                # Pack workflow_args into a single config key for the multi-run Snakefile to consume.
                workflow_args_str = " ".join(workflow_args)
                cmd += ["--config", f"workflow_args={workflow_args_str}"]

            # Non-multiomics: pass through cleaned options directly to the single-assay snakemake
            if not use_multiomics:
                # For single-assay, it's safe to pass everything through; users expect flags to be honored.
                cmd += cleaned_opts

            # Add configfile for single-assay mode only
            if config_file:
                cmd += ["--configfile", str(config_file)]

            # Run in project directory for multiomics to avoid lock conflicts
            if use_multiomics:
                cmd += ["--directory", "."]

            # Add top-level profile if requested (both modes)
            if profile_path:
                cmd += ["--profile", str(profile_path)]
                if profile:
                    logger.info(f"Using custom Snakemake profile: {profile_path}")
                else:
                    logger.info(
                        f"Using Snakemake profile preset '{preset}' -> {profile_path}"
                    )

            if queue and preset.startswith("s") and not use_multiomics:
                cmd += ["--default-resources", f"slurm_partition={queue}"]

            # Optional: print nice ASCII logo if present
            logo_trav = pkg_root_trav.joinpath("data").joinpath("logo.txt")
            try:
                with resources.as_file(logo_trav) as lp:
                    try:
                        print(Path(lp).read_text())
                    except Exception:
                        pass
            except Exception:
                pass

            # Final working directory setup
            cwd = str(Path(".").resolve())
            os.chdir(cwd)
            os.environ["PWD"] = cwd

            if print_cmd:
                logger.info("Snakemake command:\n$ " + " ".join(map(str, cmd)))

            # Execute and propagate exit code
            completed = subprocess.run(cmd, cwd=cwd)
            raise typer.Exit(code=completed.returncode)
    except typer.Exit:
        raise
    except Exception as e:
        logger.exception("Failed to run snakemake: %s", e)
        raise typer.Exit(code=1)
