"""Generate SeqNado design CSV from FASTQ files command."""
from __future__ import annotations

import re
from pathlib import Path
from typing import List, Optional

import typer
from loguru import logger
import pandas as pd

from seqnado import Assay
from seqnado.cli.app_instance import app
from seqnado.cli.autocomplete import _find_fastqs, _assay_names, fastq_autocomplete, assay_autocomplete
from seqnado.cli.utils import (
    _configure_logging,
    validate_assay,
    generate_design_dataframe,
    verbose_option,
)


_PROTECTED_EXTRACTION_COLUMNS = {"uid", "r1", "r2", "r1_control", "r2_control"}


def _compile_metadata_pattern(expression: str, option: str) -> re.Pattern[str]:
    """Compile a metadata pattern, accepting Python and PCRE named groups."""
    # Python spells a named group ``(?P<name>...)`` while PCRE users commonly
    # write ``(?<name>...)``. Convert only valid group names, leaving lookbehind
    # expressions such as ``(?<=...)`` and ``(?<!...)`` untouched.
    python_expression = re.sub(
        r"\(\?<([A-Za-z_]\w*)>", r"(?P<\1>", expression
    )
    try:
        return re.compile(python_expression)
    except re.error as exc:
        raise ValueError(f"Invalid {option} {expression!r}: {exc}") from exc


def _parse_metadata_extractors(
    metadata_regex: list[str] | None,
    metadata_field: list[str] | None,
) -> list[tuple[str, re.Pattern[str], int | str]]:
    """Validate extraction options and return their destination columns.

    Named-group patterns contribute one extractor per named group. Explicit
    ``NAME=REGEX`` patterns must contain exactly one capture group.
    """
    extractors: list[tuple[str, re.Pattern[str], int | str]] = []
    destinations: set[str] = set()

    def add(destination: str, pattern: re.Pattern[str], group: int | str) -> None:
        if destination in _PROTECTED_EXTRACTION_COLUMNS:
            raise ValueError(
                f"'{destination}' is a protected FASTQ path/index column and cannot be extracted"
            )
        if destination in destinations:
            raise ValueError(f"Duplicate metadata destination column '{destination}'")
        destinations.add(destination)
        extractors.append((destination, pattern, group))

    for expression in metadata_regex or []:
        pattern = _compile_metadata_pattern(expression, "--metadata-regex")
        if not pattern.groupindex:
            raise ValueError(
                f"--metadata-regex {expression!r} must contain at least one named capture group"
            )
        for name in pattern.groupindex:
            add(name, pattern, name)

    for specification in metadata_field or []:
        if "=" not in specification:
            raise ValueError(
                f"Invalid --metadata-field {specification!r}; expected NAME=REGEX"
            )
        name, expression = specification.split("=", 1)
        if not name or not expression:
            raise ValueError(
                f"Invalid --metadata-field {specification!r}; expected NAME=REGEX"
            )
        pattern = _compile_metadata_pattern(
            expression, f"--metadata-field {specification!r}"
        )
        if pattern.groups != 1:
            raise ValueError(
                f"--metadata-field {name!r} must contain exactly one capture group"
            )
        add(name, pattern, 1)

    return extractors


def _apply_metadata_extraction(
    df: pd.DataFrame,
    metadata_regex: list[str] | None,
    metadata_field: list[str] | None,
) -> pd.DataFrame:
    """Add regex-derived metadata columns using each row's R1 FASTQ basename."""
    try:
        extractors = _parse_metadata_extractors(metadata_regex, metadata_field)
    except ValueError as exc:
        logger.error(str(exc))
        raise typer.Exit(code=2) from exc

    for destination, pattern, group in extractors:
        df[destination] = [
            (match.group(group) if (match := pattern.search(Path(path).name)) else None)
            for path in df["r1"]
        ]
        logger.info(f"Extracted '{destination}' from FASTQ filenames.")
    return df


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


def _populate_grouping_column(
    df: pd.DataFrame,
    assay: Optional[str],
    by: str,
    target_column: str,
) -> pd.DataFrame:
    """
    Populate `target_column` (e.g. 'consensus_group', 'condition') from `by`, which is
    either an existing column name (used as-is) or a regex extracted from sample names.
    """
    if by in df.columns:
        df[target_column] = df[by].astype(str)
        logger.info(f"Set '{target_column}' from column '{by}'.")
        return df

    try:
        # For IP assays (ChIP, CAT), match against sample_id + ip so patterns can target antibodies
        if assay and assay.lower() in ["chip", "cat"] and "ip" in df.columns:
            samples = df["sample_id"] + df["ip"]
        else:
            samples = df["sample_id"]

        extracted = samples.str.extract(by, expand=False)
        if extracted.isnull().all():
            raise ValueError(f"No matches found with the provided regex '{by}'")

        df[target_column] = extracted.fillna("unknown")
        logger.info(f"Set '{target_column}' from regex '{by}'.")
    except Exception as e:
        logger.error(f"Failed to set '{target_column}' from '{by}': {e}")
        raise typer.Exit(code=3)

    return df


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
        None,
        "--consensus-by",
        "--group-by",
        help="Populate 'consensus_group' from an existing column name or a regex extracted from sample names. "
        "e.g. '--consensus-by ip' groups ChIP-seq/CUT&Tag samples by antibody for consensus peak calling/counting.",
    ),
    condition_by: Optional[str] = typer.Option(
        None,
        "--condition-by",
        help="Populate 'condition' from an existing column name or a regex extracted from sample names. "
        "e.g. '--condition-by \"-(control|treated)-\"' extracts the condition from sample names for bigwig comparisons.",
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
    metadata_regex: List[str] = typer.Option(
        None,
        "--metadata-regex",
        help="Regex matched against each R1 FASTQ basename; named capture groups become metadata columns. Repeatable.",
    ),
    metadata_field: List[str] = typer.Option(
        None,
        "--metadata-field",
        help="Metadata extraction in NAME=REGEX form; REGEX must contain exactly one capture group. Repeatable.",
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
        group_by: Populate 'consensus_group' from a column name or regex (e.g. 'ip' for antibody-based grouping)
        condition_by: Populate 'condition' from a column name or regex
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
                metadata_extractor=lambda frame: _apply_metadata_extraction(
                    frame, metadata_regex, metadata_field
                ),
            )

            if group_by:
                df = _populate_grouping_column(df, assay_name, group_by, "consensus_group")
            if condition_by:
                df = _populate_grouping_column(df, assay_name, condition_by, "condition")

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
        metadata_extractor=lambda frame: _apply_metadata_extraction(
            frame, metadata_regex, metadata_field
        ),
    )

    if group_by:
        df = _populate_grouping_column(df, assay, group_by, "consensus_group")
    if condition_by:
        df = _populate_grouping_column(df, assay, condition_by, "condition")

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output, index=False)
    logger.success(f"Design file saved → {output}")
