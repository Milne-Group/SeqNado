"""Download FASTQ files from GEO/SRA command."""
from __future__ import annotations

import shutil
import subprocess
from importlib import resources
from pathlib import Path
from typing import Optional

import typer
from loguru import logger

from seqnado.cli.app_instance import app
from seqnado.cli.autocomplete import assay_autocomplete
from seqnado.cli.snakemake_builder import SnakemakeCommandBuilder
from seqnado.cli.utils import (
    _configure_logging,
    _pkg_traversable,
    require_snakemake,
    generate_design_dataframe,
    resolve_profile,
    execute_snakemake,
    verbose_option,
    preset_option,
    dry_run_option,
    cores_option,
)


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
    """
    Download FASTQ files from GEO/SRA and optionally generate a design file.

    This command:
    1. Reads a metadata TSV file from GEO/ENA with run accession information
    2. Uses Snakemake to download FASTQs with retry logic via prefetch/fasterq-dump
    3. Optionally generates a SeqNado design file from the downloaded FASTQs

    The metadata TSV file should contain at least these columns:
    - run_accession (e.g., SRR123456)
    - sample_title (sample name)
    - library_name (e.g., GSM identifier)

    Args:
        metadata_tsv: TSV file from GEO/ENA with run accession information
        outdir: Output directory for downloaded FASTQ files
        assay: Assay type for generating design file after download
        design_output: Output path for design CSV
        cores: Number of parallel download jobs
        preset: Snakemake job profile preset for downloads
        profile: Path to a Snakemake profile directory (overrides preset)
        dry_run: Show what would be downloaded without downloading
        verbose: Increase logging verbosity

    Example:
        seqnado download filereport.tsv --outdir geo_downloads --assay rna -c 8
    """
    _configure_logging(verbose)

    # Local imports
    import pandas as pd
    import yaml

    from seqnado.inputs import Assay as AssayEnum

    require_snakemake()

    # Log assay parameter early if provided
    if assay:
        logger.info(f"Assay type: {assay}")

    # Read and parse metadata TSV
    logger.info(f"Reading metadata from {metadata_tsv}")
    try:
        samples_df = pd.read_csv(metadata_tsv, sep="\t")
        required_cols = ["run_accession", "sample_title", "library_name"]
        missing_cols = [col for col in required_cols if col not in samples_df.columns]
        if missing_cols:
            logger.error(f"Missing required columns in TSV: {', '.join(missing_cols)}")
            logger.info(f"Required columns: {', '.join(required_cols)}")
            logger.info(f"Available columns: {', '.join(samples_df.columns)}")
            raise typer.Exit(code=1)

        # Check for library_layout column
        has_layout = "library_layout" in samples_df.columns
        if not has_layout:
            logger.warning(
                "No 'library_layout' column found in TSV. "
                "Will attempt to detect layout during download, but this may be slower."
            )
            logger.info(
                "For best results, include 'library_layout' column with values 'PAIRED' or 'SINGLE'"
            )
    except Exception as e:
        logger.error(f"Failed to read metadata TSV: {e}")
        raise typer.Exit(code=1)

    # Build sample dictionaries separated by layout
    geo_samples_paired = {}
    geo_samples_single = {}
    geo_samples_unknown = {}

    for _, row in samples_df.iterrows():
        sample_name = f"{row['library_name']}-{row['sample_title']}"
        sample_info = {
            "srr": row["run_accession"],
            "gsm": row["library_name"],
            "sample": row["sample_title"],
        }

        if has_layout:
            layout = str(row["library_layout"]).upper()
            if layout == "PAIRED":
                geo_samples_paired[sample_name] = sample_info
            elif layout == "SINGLE":
                geo_samples_single[sample_name] = sample_info
            else:
                logger.warning(
                    f"Unknown library_layout '{layout}' for {sample_name}, treating as unknown"
                )
                geo_samples_unknown[sample_name] = sample_info
        else:
            # No layout info - will need to detect during download
            geo_samples_unknown[sample_name] = sample_info

    logger.info(
        f"Found {len(geo_samples_paired)} paired-end, "
        f"{len(geo_samples_single)} single-end, "
        f"{len(geo_samples_unknown)} unknown layout samples"
    )

    if geo_samples_unknown:
        logger.error(
            f"Cannot proceed with {len(geo_samples_unknown)} samples without layout information. "
            "Please add 'library_layout' column to your TSV with values 'PAIRED' or 'SINGLE'."
        )
        raise typer.Exit(code=1)

    # Create output directory
    outdir.mkdir(parents=True, exist_ok=True)

    # Create Snakemake config for downloads
    temp_config = {
        "geo_samples_paired": geo_samples_paired,
        "geo_samples_single": geo_samples_single,
        "geo_outdir": str(outdir.resolve()),
    }

    config_file = Path("seqnado_output/logs/geo_download/geo_download_config.yaml")
    config_file.parent.mkdir(parents=True, exist_ok=True)
    with open(config_file, "w") as f:
        yaml.dump(temp_config, f)

    # Get the download.smk file from package
    pkg_root_trav = _pkg_traversable("seqnado")
    download_smk_trav = (
        pkg_root_trav.joinpath("workflow")
        .joinpath("rules")
        .joinpath("geo")
        .joinpath("download.smk")
    )

    profile_ctx, profile_path, is_custom = resolve_profile(
        preset, profile, pkg_root_trav, warn_unknown=True
    )

    # Build Snakemake command
    with (
        resources.as_file(download_smk_trav) as download_smk,
        profile_ctx as resolved_profile_path,
    ):
        # Use resolved_profile_path from context if profile_path was None (Traversable case)
        final_profile_path = profile_path or resolved_profile_path
        
        # Initialize builder
        builder = SnakemakeCommandBuilder(Path(download_smk), cores)
        
        # Add configfile
        builder.add_configfile(config_file)
        
        # Add target rule
        builder.add_target("geo_download_all")
        
        # Add container support for Apptainer/Singularity
        builder.add_container_support()

        # Add dry-run if requested
        if dry_run:
            builder.add_dry_run()

        # Add profile with logging
        builder.add_profile_with_logging(final_profile_path, preset, is_custom)

        # Build command
        cmd = builder.build()
        
        logger.info("Starting GEO download with Snakemake...")
        
        # Execute Snakemake
        cwd = str(Path.cwd())
        exit_code = execute_snakemake(cmd, cwd, verbose or dry_run)

        if exit_code != 0:
            logger.error(f"Snakemake failed with exit code {exit_code}")
            raise typer.Exit(code=exit_code)

        logger.success("GEO download completed!")

        # Generate design file if assay is specified
        if assay and not dry_run:
            logger.info(f"\nGenerating design file for {assay}...")

            # Find downloaded FASTQ files
            fastq_files = sorted(outdir.glob("*.fastq.gz"))
            if not fastq_files:
                logger.error(f"No FASTQ files found in {outdir}")
                raise typer.Exit(code=1)

            logger.info(f"Found {len(fastq_files)} FASTQ files")

            df = generate_design_dataframe(
                assay=assay,
                fastq_files=fastq_files,
                ip_to_control_map={},
                interactive=True,
                accept_all_defaults=True,
                deseq2_pattern=None,
            )

            # Save design file
            if design_output is None:
                design_output = outdir / f"metadata_{assay}.csv"

            design_output.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(design_output, index=False)
            logger.success(f"Design file saved → {design_output}")
