"""Initialize SeqNado user environment command."""
from __future__ import annotations

import os
import shutil
import subprocess
from importlib import resources
from pathlib import Path

import typer
from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    _read_json,
    _write_json,
    dry_run_option,
    verbose_option,
)
from seqnado.utils import get_preset_profiles


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
    """
    Initialize SeqNado user environment.

    - Logs the current Conda environment if active (optional).
    - Runs packaged Apptainer/Singularity init (if `apptainer` on PATH).
    - Ensures ~/.config/seqnado/genome_config.json exists (template or preset).
    """
    _configure_logging(verbose)

    conda_env = os.environ.get("CONDA_DEFAULT_ENV")
    logger.info(f"Conda environment: {conda_env or 'none'}")

    # Apptainer/Singularity bootstrap
    if shutil.which("apptainer"):
        init_script_trav = _pkg_traversable("seqnado").joinpath("init.sh")
        try:
            with resources.as_file(init_script_trav) as init_script:
                if init_script.exists():
                    logger.info(f"Configuring Apptainer/Singularity via {init_script}")
                    if not dry_run:
                        try:
                            subprocess.run(["bash", str(init_script)], check=True)
                        except subprocess.CalledProcessError as e:
                            logger.warning(
                                f"Apptainer init script failed (continuing): {e}"
                            )
                        except Exception as e:
                            logger.warning(f"Skipping Apptainer init due to error: {e}")
                    else:
                        logger.info(f"[dry-run] Would execute: bash {init_script}")
                else:
                    logger.warning(
                        "Apptainer init script not found in package; skipping."
                    )
        except Exception as e:
            logger.warning("Could not access package init script: %s", e)
    else:
        logger.info("Apptainer not found on PATH; skipping container setup.")

    # Genome config
    cfg_dir = Path.home().joinpath(".config", "seqnado")
    cfg_dir.mkdir(parents=True, exist_ok=True)
    genome_config = cfg_dir.joinpath("genome_config.json")

    data_pkg = "seqnado.data"
    template_genomes = "preset_genomes.json"
    template_config = "genomes_template.json"
    template_name = template_genomes if preset else template_config
    template_trav = _pkg_traversable(data_pkg).joinpath(template_name)
    template_config_trav = _pkg_traversable(data_pkg).joinpath(template_config)

    # move to ~/.config/snakemake.
    profile_target_dir = Path.home().joinpath(".config", "snakemake")
    profile_target_dir.mkdir(parents=True, exist_ok=True)
    profiles = get_preset_profiles()
    for profile in profiles.values():
        profile_src_trav = _pkg_traversable("seqnado.workflow.envs.profiles").joinpath(
            profile
        )
        profile_dest = profile_target_dir.joinpath(profile)
        if profile_dest.exists():
            logger.info(f"Snakemake profile already exists: {profile_dest}")
            continue
        try:
            with resources.as_file(profile_src_trav) as profile_src:
                if dry_run:
                    logger.info(
                        f"[dry-run] Would copy Snakemake profile {profile_src.name} to {profile_dest}"
                    )
                else:
                    shutil.copytree(profile_src, profile_dest)
                    logger.info(f"Copied Snakemake profile to {profile_dest}")
        except Exception as e:
            logger.error(f"Failed to copy Snakemake profile {profile}: {e}")

    if genome_config.exists():
        logger.info(f"Found genome config: {genome_config}")
        try:
            genome_data = _read_json(genome_config)
        except Exception as e:
            logger.warning(
                f"Could not parse existing genome config ({e}); leaving as-is."
            )
            genome_data = None

        try:
            with resources.as_file(template_config_trav) as template_path:
                template_data = _read_json(Path(template_path))
            if genome_data == template_data:
                logger.warning(
                    "Genome config matches the template exactly. Please update paths as needed."
                )
            else:
                logger.info("Genome config appears to be customized; leaving as-is.")
        except Exception as e:
            logger.debug("Could not read packaged template to compare: %s", e)
    else:
        try:
            with resources.as_file(template_trav) as template_path:
                if not Path(template_path).exists():
                    logger.error(
                        "Packaged genome templates missing—installation may be corrupted."
                    )
                    raise typer.Exit(code=1)
                if dry_run:
                    logger.info(
                        f"[dry-run] Would create {genome_config} from {template_path.name}"
                    )
                else:
                    template = _read_json(Path(template_path))
                    _write_json(genome_config, template)
                    logger.info(
                        "Created genome config "
                        + (
                            "from preset genomes."
                            if preset
                            else "template (please update paths)."
                        )
                    )
        except typer.Exit:
            raise
        except Exception as e:
            logger.error(f"Failed to write genome config: {e}")
            raise typer.Exit(code=1)

    logger.success("Initialization complete.")
