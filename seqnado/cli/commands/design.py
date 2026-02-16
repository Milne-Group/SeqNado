"""Generate SeqNado design CSV from FASTQ files command."""
from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

import typer
from loguru import logger

from seqnado import Assay
from seqnado.cli.app_instance import app
from seqnado.cli.autocomplete import _find_fastqs, _assay_names, fastq_autocomplete, assay_autocomplete
from seqnado.cli.utils import (
    _configure_logging,
    validate_assay,
    generate_design_dataframe,
    verbose_option,
)


def _parse_ip_to_control_pairings(
    pairing_str: str,
) -> dict[str, str]:
    """
    Parse ip-to-control pairings from a string of the form:
    'antibody1:control1,antibody2:control2'
    """
    pairings = {}
    for pair in pairing_str.split(","):
        try:
            ip, control = pair.split(":")
            pairings[ip.strip()] = control.strip()
        except ValueError:
            logger.error(
                f"Invalid ip-to-control pairing format: '{pair}'. Expected 'antibody:control'."
            )
            raise typer.Exit(code=2)
    return pairings


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
    group_by: Optional[str] = typer.Option(
        None, "--group-by", help="Group samples by a regular expression or a column name."
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
    """
    Generate a SeqNado design CSV from FASTQ files for ASSAY.
    
    If no assay is provided, multiomics mode is used to detect assay directories.
    
    Args:
        assay: Assay type. Options: " + ", ".join(_assay_names()) + ". If omitted, multiomics mode is used.
        files: FASTQ files to process
        output: Output CSV filename (default: metadata_{assay}.csv)
        ip_to_control: List of antibody,control pairings for IP assays
        group_by: Group samples by a regular expression or a column name
        auto_discover: Search common folders if none provided
        interactive: Interactively offer to add missing columns using schema defaults
        accept_all_defaults: Non-interactive: auto-add only columns with schema default
        deseq2_pattern: Regex pattern to extract DESeq2 groups from sample names
        verbose: Increase logging verbosity
    """
    _configure_logging(verbose)

    # Local imports
    import pandas as pd

    from seqnado.inputs import Assay as AssayEnum
    from seqnado.inputs import FastqCollection, FastqCollectionForIP
    from seqnado.inputs.validation import DesignDataFrame

    # Parse ip-to-control pairings
    ip_to_control_map: dict[str, str] = {}
    if ip_to_control:
        ip_to_control_map = _parse_ip_to_control_pairings(ip_to_control)

    # Handle multiomics mode
    if assay is None or assay == Assay.MULTIOMICS.clean_name:
        logger.info(
            "Multiomics mode: searching for assay-specific fastq subdirectories"
        )

        # Look for fastqs/<assay>/ directories
        fastqs_base = Path("fastqs")
        if not fastqs_base.exists():
            logger.error(
                "No 'fastqs/' directory found. Run 'seqnado config' first to create the directory structure."
            )
            raise typer.Exit(code=1)

        # Find all subdirectories in fastqs/ that match known assay names
        available_assays = AssayEnum.all_assay_clean_names()
        found_assay_dirs = {}

        for assay_dir in fastqs_base.iterdir():
            if assay_dir.is_dir() and assay_dir.name in available_assays:
                # Check if there are any fastq files in this directory
                fastq_files = list(assay_dir.glob("*.fastq.gz"))
                if fastq_files:
                    found_assay_dirs[assay_dir.name] = fastq_files
                    logger.info(f"Found {len(fastq_files)} FASTQ files in {assay_dir}")

        if not found_assay_dirs:
            logger.error(
                "No FASTQ files found in any assay subdirectories under fastqs/"
            )
            logger.info("Expected structure: fastqs/{assay}/{files}.fastq.gz")
            logger.info(f"Valid assay names: {', '.join(available_assays)}")
            raise typer.Exit(code=1)

        # Generate metadata for each assay
        generated_files = []
        for assay_name, fastq_files in found_assay_dirs.items():
            logger.info(f"\n=== Generating metadata for {assay_name} ===")

            df = generate_design_dataframe(
                assay=assay_name,
                fastq_files=fastq_files,
                ip_to_control_map=ip_to_control_map,
                interactive=interactive,
                accept_all_defaults=accept_all_defaults,
                deseq2_pattern=deseq2_pattern,
            )

            # Save metadata file
            metadata_file = Path(f"metadata_{assay_name}.csv")
            metadata_file.parent.mkdir(parents=True, exist_ok=True)
            df.to_csv(metadata_file, index=False)
            generated_files.append(metadata_file)
            logger.success(f"Design file saved → {metadata_file}")

        logger.success(f"\nGenerated {len(generated_files)} metadata files:")
        for f in generated_files:
            logger.info(f"  - {f}")

        return

    # Regular single-assay mode
    validate_assay(assay)

    fastq_paths: List[Path] = []
    for p in files or []:
        if p.suffixes[-2:] == [".fastq", ".gz"]:
            fastq_paths.append(p)

    if not fastq_paths and auto_discover:
        hints = [".", "fastqs", "fastq", "data", "data/fastqs"]
        logger.info(f"No FASTQs provided; searching {hints}")
        fastq_paths = _find_fastqs(hints)

    if not fastq_paths:
        logger.error("No FASTQ files provided or found.")
        typer.echo(
            "Tip: provide paths explicitly or place *.fastq.gz files in one of: ./, fastqs/, fastq/, data/, data/fastqs/",
            err=True,
        )
        raise typer.Exit(code=1)

    logger.info(f"Found {len(fastq_paths)} FASTQ files for assay '{assay}'")

    # Set default output filename if not provided
    if output is None:
        output = Path(f"metadata_{assay}.csv")

    df = generate_design_dataframe(
        assay=assay,
        fastq_files=fastq_paths,
        ip_to_control_map=ip_to_control_map,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
        deseq2_pattern=deseq2_pattern,
    )

    if group_by:
        # Try matching against a column name first
        if group_by in df.columns:
            df["consensus_group"] = df[group_by].astype(str)
            logger.info(
                f"Grouped samples by column '{group_by}' into 'consensus_group'."
            )
        else:
            # Treat group_by as a regex pattern to extract from sample_id
            try:
                # For IP assays (ChIP, CAT), concatenate sample_id and ip
                if assay and assay.lower() in ["chip", "cat"]:
                    samples = df["sample_id"] + df["ip"]
                else:
                    samples = df["sample_id"]

                df["consensus_group"] = samples.str.extract(group_by, expand=False)
                if df["consensus_group"].isnull().all():
                    raise ValueError(
                        f"No matches found with the provided regex '{group_by}'"
                    )

                df["consensus_group"] = df["consensus_group"].fillna("unknown")
                logger.info(
                    f"Grouped samples by regex '{group_by}' into 'consensus_group'."
                )
            except Exception as e:
                logger.error(f"Failed to group by '{group_by}': {e}")
                raise typer.Exit(code=3)

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    logger.success(f"Design file saved → {output}")
