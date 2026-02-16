"""Build workflow configuration YAML command."""
from __future__ import annotations

import typer
from datetime import date
from importlib import resources
from pathlib import Path
from typing import Optional

from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.autocomplete import _assay_names, assay_autocomplete
from seqnado.cli.utils import _configure_logging, _pkg_traversable, validate_assay, verbose_option


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
    """
    Build a workflow configuration YAML for the selected ASSAY.
    
    If no assay is provided, multiomics mode is used.
    
    Args:
        assay: Assay type (e.g., 'rna', 'chip', etc.) or None for multiomics
        make_dirs: Create output project directory or fastq subdir
        render_options: Render all options (even if not used by the workflow)
        output: Explicit path for the rendered config file
        verbose: Increase logging verbosity
        interactive: Interactively prompt for config values. Non-interactive 
            mode only works for single assay configs (except MCC and multiomics)
    """
    _configure_logging(verbose)

    # Local imports to keep help fast
    from importlib.metadata import version as _pkg_version

    from seqnado.config import (
        build_default_workflow_config,
        build_multiomics_config,
        build_workflow_config,
        render_config,
        render_multiomics_configs,
    )
    from seqnado.inputs import Assay

    seqnado_version = _pkg_version("seqnado")

    # If no assay provided, use multiomics mode
    if assay is None or assay.lower() == "multiomics":
        logger.info("Building multiomics configuration with multiple assays")

        try:
            multiomics_config, assay_configs = build_multiomics_config(
                seqnado_version, interactive=interactive
            )
        except Exception as e:
            logger.error(f"Failed to build multiomics configuration: {e}")
            raise typer.Exit(code=1)

        # Determine output directory
        if make_dirs:
            # Use first assay's project name for directory
            first_config = next(iter(assay_configs.values()))
            dirname = f"{date.today().isoformat()}_{first_config.project.name}"
            outdir = Path(dirname)

            # Create assay-specific fastq subdirectories
            for assay_name in assay_configs.keys():
                fastq_dir = outdir / "fastqs" / assay_name
                fastq_dir.mkdir(parents=True, exist_ok=True)
                logger.info(f"Created fastq directory: {fastq_dir}")
        else:
            outdir = Path(".")

        # Render all config files
        tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")
        try:
            with resources.as_file(tpl_trav) as tpl_path:
                if not Path(tpl_path).exists():
                    logger.error(
                        "Packaged config template missing—installation may be corrupted."
                    )
                    raise typer.Exit(code=1)

                generated_files = render_multiomics_configs(
                    multiomics_config=multiomics_config,
                    assay_configs=assay_configs,
                    template=Path(tpl_path),
                    output_dir=outdir,
                )

                logger.success(
                    f"Generated {len(generated_files)} config files in {outdir}:"
                )
                for f in generated_files:
                    logger.info(f"  - {f.name}")

        except typer.Exit:
            raise
        except Exception as e:
            logger.error(f"Failed to render multiomics configs: {e}")
            raise typer.Exit(code=1)

        return

    # Regular single-assay config
    validate_assay(assay)

    assay_obj = Assay.from_clean_name(assay)

    if not interactive:
        logger.info("Running in non-interactive mode; using defaults where possible.")
        workflow_config = build_default_workflow_config(assay_obj)
        render_options = True
    else:
        workflow_config = build_workflow_config(assay_obj, seqnado_version)
        if not workflow_config:
            logger.error("Failed to build workflow configuration.")
            raise typer.Exit(code=1)

    if not make_dirs:
        config_output = output or Path(f"config_{assay_obj.clean_name}.yaml")
        config_output.parent.mkdir(parents=True, exist_ok=True)
    else:
        dirname = f"{date.today().isoformat()}_{workflow_config.project.name}"
        outdir = Path(dirname)
        (outdir / "fastqs").mkdir(parents=True, exist_ok=True)
        logger.info(f"Created output directory: {outdir / 'fastqs'}")
        config_output = output or (outdir / f"config_{assay_obj.clean_name}.yaml")

    tpl_trav = _pkg_traversable("seqnado.data").joinpath("config_template.jinja")
    try:
        with resources.as_file(tpl_trav) as tpl_path:
            if not Path(tpl_path).exists():
                logger.error(
                    "Packaged config template missing—installation may be corrupted."
                )
                raise typer.Exit(code=1)
            try:
                render_config(
                    template=Path(tpl_path),
                    workflow_config=workflow_config,
                    outfile=config_output,
                    all_options=render_options,
                )
            except Exception as e:
                logger.error(f"Failed to render config: {e}")
                raise typer.Exit(code=1)
    except typer.Exit:
        raise
    except Exception as e:
        logger.error("Failed to locate or access packaged template: %s", e)
        raise typer.Exit(code=1)

    logger.success(f"Wrote config → {config_output}")

    # Check if FastqScreen config exists at default location
    default_fastqscreen_config = Path.home() / ".config/seqnado/fastq_screen.conf"
    if not default_fastqscreen_config.exists():
        logger.warning(f"FastqScreen config not found at {default_fastqscreen_config}")
        if interactive:
            generate_fs = typer.confirm(
                "Would you like to generate a FastqScreen config now?", default=True
            )
            if generate_fs:
                contaminant_path_input = typer.prompt(
                    "Path to contaminant references (leave empty to skip contaminants)",
                    default="",
                    show_default=False,
                )

                from seqnado.config.fastq_screen_generator import (
                    generate_fastq_screen_config,
                    load_genome_configs_for_fastqscreen,
                )

                try:
                    genome_configs = load_genome_configs_for_fastqscreen()
                    generate_fastq_screen_config(
                        genome_configs=genome_configs,
                        output_path=default_fastqscreen_config,
                        threads=8,
                        include_contaminants=bool(contaminant_path_input),
                        contaminant_base_path=contaminant_path_input
                        if contaminant_path_input
                        else None,
                    )
                except Exception as e:
                    logger.warning(f"Failed to generate FastqScreen config: {e}")
                    logger.info(
                        "You can generate it later with: seqnado genomes fastqscreen"
                    )
        else:
            logger.info(
                "Run 'seqnado genomes fastqscreen' to generate FastqScreen configuration"
            )
