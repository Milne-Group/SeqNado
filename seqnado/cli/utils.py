"""Shared CLI utilities for logging, file I/O, and package resources."""

from __future__ import annotations

import contextlib
import json
import os
import sys
import tempfile
from importlib import resources
from pathlib import Path
from typing import Optional, Tuple

import click
import typer
from loguru import logger

from seqnado.utils import get_preset_profiles, resolve_profile_path

# Rich console capture may or may not be available — keep an opt-in safe reference
try:
    from rich.console import Console  # type: ignore
    from rich.text import Text  # type: ignore

    _RICH_CONSOLE = Console(force_terminal=True)
except Exception:
    _RICH_CONSOLE = None
    Text = None  # type: ignore

# ==================== Constants ====================

TOP_LEVEL_PASS_THROUGH = (
    "-n",
    "--dry-run",
    "--printshellcmds",
    "--unlock",
    "--rerun-incomplete",
    "--show-failed-logs",
)


def _pkg_traversable(pkg: str):
    """
    Return the importlib.resources Traversable for a package.
    Use resources.as_file(...) when you need a filesystem Path for a particular file.
    """
    return resources.files(pkg)


def _read_json(path: Path) -> dict:
    """Read and parse JSON file."""
    return json.loads(path.read_text(encoding="utf-8"))


def _write_json(path: Path, data: dict) -> None:
    """
    Atomically write JSON to `path` using a temporary file + os.replace.
    Attempt to set secure permissions but ignore if not supported.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = Path(tempfile.mktemp(dir=str(path.parent)))
    try:
        tmp.write_text(json.dumps(data, indent=4), encoding="utf-8")
        os.replace(str(tmp), str(path))
        try:
            os.chmod(path, 0o600)
        except Exception:
            # best-effort — ignore permission errors
            logger.debug("Unable to set file permissions on %s", path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink(missing_ok=True)
            except Exception:
                # best-effort cleanup
                pass


def _snakemake_available() -> bool:
    """Check if snakemake is available on PATH."""
    import shutil

    return shutil.which("snakemake") is not None


def _configure_logging(verbose: bool) -> None:
    """Configure loguru logging level based on verbosity."""
    logger.remove()
    logger.add(sys.stderr, level="DEBUG" if verbose else "INFO", colorize=True)


def _style_name_with_rich(name: str, style: str = "bold cyan") -> str:
    """Return `name` styled with Rich as an ANSI string, or plain name on fallback."""
    if _RICH_CONSOLE is None or Text is None:
        return name
    with _RICH_CONSOLE.capture() as cap:
        _RICH_CONSOLE.print(Text(name, style=style), end="")
    return cap.get()

# ==================== CLI Helper Functions ====================


def resolve_profile_context(
    preset: Optional[str],
    profile: Optional[Path],
    pkg_root_trav=None,
    *,
    warn_unknown: bool = False,
) -> Tuple[contextlib.AbstractContextManager, bool]:
    """
    Resolve a Snakemake profile preset or custom path into a context manager.

    Returns a context manager that yields a profile path (or None) and a flag
    indicating whether the profile is custom (True) or preset (False).
    """
    if profile:
        profile_path = Path(profile).expanduser()
        if not profile_path.exists():
            logger.error(f"Snakemake profile not found: {profile_path}")
            raise typer.Exit(code=1)
        return contextlib.nullcontext(profile_path), True

    if not preset:
        return contextlib.nullcontext(None), False

    profile_result = resolve_profile_path(preset, pkg_root_trav)
    if profile_result:
        # If it's already a Path (user-installed profile), use nullcontext
        # Otherwise it's a Traversable (bundled profile), use resources.as_file
        if isinstance(profile_result, Path):
            return contextlib.nullcontext(profile_result), False
        else:
            return resources.as_file(profile_result), False

    if warn_unknown:
        profiles = get_preset_profiles()
        logger.warning(f"Unknown preset '{preset}'. Available: {list(profiles.keys())}")

    return contextlib.nullcontext(None), False


def require_snakemake() -> None:
    """
    Check that snakemake is available on PATH.
    Raises typer.Exit with code 127 if not found.
    """
    if not _snakemake_available():
        logger.error(
            "`snakemake` not found on PATH. Install/activate the environment that provides it."
        )
        raise typer.Exit(code=127)


def validate_assay(assay: str) -> None:
    """
    Validate an assay name against available assay types.
    
    Args:
        assay: Assay name to validate
        
    Raises:
        typer.Exit with code 2 if assay is not recognized
    """
    from seqnado.inputs import Assay as AssayEnum
    
    if assay not in AssayEnum.all_assay_clean_names():
        allowed = ", ".join(AssayEnum.all_assay_clean_names())
        logger.error(f"Unknown assay '{assay}'. Allowed: {allowed}")
        raise typer.Exit(code=2)


def generate_design_dataframe(
    assay: str,
    fastq_files: list[Path],
    ip_to_control_map: dict[str, str] | None = None,
    interactive: bool = True,
    accept_all_defaults: bool = False,
    deseq2_pattern: str | None = None,
) -> "pandas.DataFrame":
    """
    Generate a design dataframe from FASTQ files for a given assay.
    
    This encapsulates the repeated pattern of:
    - Creating FastqCollection/FastqCollectionForIP
    - Converting to dataframe and sorting
    - Extracting schema candidates
    - Applying interactive defaults
    
    Args:
        assay: Assay name (clean name, e.g. 'rna', 'chip')
        fastq_files: List of Path objects to FASTQ files
        ip_to_control_map: For IP assays, mapping of antibody to control sample names
        interactive: Whether to interactively prompt for defaults
        accept_all_defaults: Non-interactive: use only columns with schema defaults
        deseq2_pattern: Regex pattern to extract DESeq2 groups from sample names
        
    Returns:
        Pandas DataFrame with design metadata
    """
    import pandas as pd
    
    from seqnado.inputs import Assay as AssayEnum
    from seqnado.inputs import FastqCollection, FastqCollectionForIP
    from seqnado.inputs.validation import DesignDataFrame
    from seqnado.cli.data_validation import (
        _apply_interactive_defaults,
        _extract_candidate_defaults_from_schema,
    )
    
    assay_obj = AssayEnum.from_clean_name(assay)
    
    # Create design object (IP assays use FastqCollectionForIP, others use FastqCollection)
    if assay_obj in {AssayEnum.CHIP, AssayEnum.CAT}:
        design_obj = FastqCollectionForIP.from_fastq_files(
            assay=assay_obj,
            files=fastq_files,
            ip_to_control_map=ip_to_control_map or {},
        )
    else:
        design_obj = FastqCollection.from_fastq_files(
            assay=assay_obj, files=fastq_files
        )
    
    # Convert to dataframe and sort
    df = design_obj.to_dataframe().sort_values("sample_id")
    
    # Extract schema candidates and apply defaults
    schema_candidates = _extract_candidate_defaults_from_schema(DesignDataFrame, assay_obj)
    df = _apply_interactive_defaults(
        df,
        schema_candidates,
        interactive=interactive,
        accept_all_defaults=accept_all_defaults,
        deseq2_pattern=deseq2_pattern,
        assay=assay_obj,
    )
    
    return df


# ==================== Typer Option Factories ====================


def verbose_option() -> bool:
    """Factory for the --verbose/-v option used across all commands."""
    return typer.Option(
        False, "--verbose", "-v", help="Increase logging verbosity."
    )


def dry_run_option() -> bool:
    """Factory for the --dry-run option."""
    return typer.Option(
        False, "-n", "--dry-run", help="Show actions without executing them."
    )


def cores_option() -> int:
    """Factory for the -c/--cores option."""
    return typer.Option(
        4, "-c", "--cores", help="Number of cores/parallel jobs."
    )


def preset_option() -> str:
    """Factory for the --preset option with available profile choices."""
    from seqnado.cli.autocomplete import _profile_autocomplete
    
    return typer.Option(
        "le",
        "--preset",
        click_type=click.Choice(_profile_autocomplete(), case_sensitive=False),
        help="Snakemake job profile preset.",
        case_sensitive=False,
    )

def _resolve_working_dir() -> None:
    """
    Corrects a common issue where that prevents apptainer containers from seeing the current working directory.

    This is due to symlinks or bind mounts that cause the container's view of the filesystem to differ from the host's PWD.

    We resolve the current working directory to an absolute path and set PWD to that, which should be visible inside the container regardless of how it's launched.
    """
    cwd = str(Path(".").resolve())
    os.chdir(cwd)
    os.environ["PWD"] = cwd